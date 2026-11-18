!     checked into cvs archive
      subroutine ptide_noaa(nebcm,numebc,numc,damp,dg,z0,
     &     idy0,imo0,iyr0,ihr,nday,tzero)
!
!     Tidal prediction program modified based on NOS/NOAA pred.f code
!     The NOS/NOAA pred.f code is transfered as a subroutine
!
!     modified by:
!     Shejun Fan: Stevens Institute of Technology
!     Nickitas Georgas: Stevens Institute of Technology
!
      parameter(nrec=40000) ! NG: Fits a year's worth of 15min data
      real damp(nebcm,numc-1),dg(nebcm,numc-1)
     &     ,etabc(nebcm,nrec),z0(nebcm),tzero
 
      common/ptide/amp(37),epoc(37),storn(nrec),ptime(nrec)

!     damp and dg:  S2,  M2, N2,  K1,   P1,  O1
!     amp and epoc: M2,  S2, N2,  K1,   M4,  O1,   M6, MK3, S4,  MN4,
!                   NU2, S6, MU2, 2N2,  OO1, LAM2, S1, M1,  J1,  MM,
!                   SSA, SA, MSF, MF,   RHO, Q1,   T2, R2,  2Q1, P1,
!                   2SM2,M3, L2,  2MK3, K2,  M8,   MS4
      dt=0.5        ! Create series every 30 minutes BUT BE CAREFUL!!!!!:
!     dt=0.25       ! Create series every 15 minutes BUT BE CAREFUL!!!!!:
!     dt=1.         ! Create series every hour BE CAREFUL!!!!!:
      imn=0         ! Since start date in run_data does not have minutes
                    ! and no minutes are passed in from ecom3d to bcdata to 
                    ! ptide, the duration of all sequential (cold/hot/hot) runs
                    ! that use ptide_noaa constituents need to be an integer
                    ! multiple of hours (say 20 days, 3 hrs long), not longer
                    ! than a year each, and the start date has to be at hr 0. NG01062010

      do ic=1,37
         amp(ic)=0.0
         epoc(ic)=0.0
      enddo
      do iebc=1,numebc
         amp(2)=damp(iebc,1)
         amp(1)=damp(iebc,2)
         amp(3)=damp(iebc,3)
         amp(4)=damp(iebc,4)
         amp(30)=damp(iebc,5)
         amp(6)=damp(iebc,6)
         epoc(2)=dg(iebc,1)
         epoc(1)=dg(iebc,2)
         epoc(3)=dg(iebc,3)
         epoc(4)=dg(iebc,4)
         epoc(30)=dg(iebc,5)
         epoc(6)=dg(iebc,6)
         call pretide(IYR0,IMO0,IDY0,IHR,IMN,NDAY,DT,Z0(iebc),ncount)
         do itime=1,ncount
!NG04212008            etabc(iebc,itime)=storn(itime)+tzero
            etabc(iebc,itime)=storn(itime)
         enddo
      enddo
      write(*,*) '****************************************************'
      write(*,*) 'Generating elevation BC every 15min with ptide_noaa.'
      write(*,'(a,i4,"/",i2,"/",i2," ",i2,":",i2)') ' Beginning:',
     +            IYR0,IMO0,IDY0,IHR,IMN
      write(*,'(f10.4,a)') tzero,' hrs after startdate in run_data.'
      write(*,*) 'Careful: The duration of all runs (cold or hot) that'
      write(*,*) 'use ptide_noaa constituents need to be an integer'
      write(*,*) 'multiple of *hours* (say 20 days, 3 hrs long), not'
      write(*,*) 'longer than a year each, and the starting IHR in the'
      write(*,*) 'run_data has to be zero, or they will be wrong.'
      write(*,*) '****************************************************'
      do itime=1,ncount
!NG04212008      do itime=1,ncount
!NG04212008         write(90,77) ptime(itime)
         write(90,77) ptime(itime)+tzero
         write(*,77) ptime(itime)+tzero,etabc(1,itime)
         write(90,77) (ETABC(I,itime),I=1,NUMEBC)
!         write(*,77) ptime(itime)
!         write(*,77) (ETABC(I,itime),I=1,NUMEBC)
      enddo
 77   FORMAT(8E14.7)
      write(*,*) '...Generated elevation BC with ptide_noaa.'

      return
      end

      subroutine pretide(IYRS,IMMS,IDDS,IHHS,MNS,NDAY,DELT,DATUM,ncount)

!CCCCCC   this program is modified from pred.f of Chris Zervas so that it can make prediction of crossing year
!CCCCC    also call equarg.f to calculate XODE and VPU, instead of reading from data file 'yr'
!CCCCC    modified on Feb. 13, 2005 by Aijun Zhang
!   Purpose:  predict water level or current using 37 NOS standard tidal constituents.
!   Author:   Aijun Zhang       
!            CSDL/NOS/NOAA 
!  compile:   lf95 pred.f -o pred.x
!  run at unix/linux command line
!         pred.x "$BEGINDATE" "$ENDDATE" $KINDAT $DELT $filein $fileout
!         BEGINDATE="2005 01 01 12 30"
!         ENDDATE=  "2005 12 31 12 30"
!         KINDAT=1, for current prediction; =2 for water level prediction
!         DELT is time interval of output time series in hours
!         XMAJOR is principle current direction in degrees, it is ignored for water level
!               ,for this version, XMAJOR is read in from filein 
!         filein is input file name which includes tide constituents with standard format of NOS
!         fileout is output name which contains predicted water level or current time series
!                              wl
!    2003 01 02 00 00 00     0.5085
!    2003 01 02 00 06 00     0.5169
!                            sp        dir  
!   2003 01 02 00 00 00     0.7091   258.0500
!   2003 01 02 00 06 00     0.7237   258.3601
           
!    1    2    1.   0.    0              ! nsta ipredk conv tconv il2
!   tss.out                                 ! Output time series file
!   0  4  15 0  6 30 0 0.1 1998 1998 106.0   IEL,IMMS,IDDS,TIMEF,IMME,IDDE,TIMEL,DELT,IYRS,IYRE,XMAJOR
!   Harmonic Analysis of Data in  325j4b05.dat                                      
!   29-Day H.A.  Beginning  4-15-1998  at Hour 17.30  along 106 degrees             
!   12718
!       1931621828140022115186641641117973376 68733224 81743495 5868 163
!       2    0   0  804  92    0   0 36211666  4431590    0   0 24821455
!       3  3513257  6521961    0   0  5803436  6463317    0   0    0   0
!       4    0   0    0   0    0   0  3113546 15863554  8262104  1122127
!       5  213  13 39053385    0   0    0   0 26691641    0   0 38082138
!       6 28911704    0   0
!   Harmonic Analysis of Data in  325j4b05.dat                                      
!   29-Day H.A.  Beginning  4-15-1998  at Hour 17.30  along 196 degrees             
! -5820
!       1 57222694 14073452 22192976 1203 154 47261638 1292 478 26243445
!       2    0   0  658 796    0   0  4302938  6232349    0   0  2953258
!       3   563430   403046    0   0   92 316  1023593    0   0    0   0
!       4    0   0    0   0    0   0   49 617  251 639   833422   113482
!       5   34 800  398 178    0   0    0   0  3172976    0   0  3833513
!       6 12661394    0   0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!****          Purpose:
!                 1. Calculates residuals
!                 2. Calculates a predicted series
! **************************************************************
!
!               UNIT=10 is the output file written in ASCII format.
!               UNIT=11 is the input observations in CDF or ASCII format.
!
!
!
!      UNIT=         PATH NAME=                    INPUT      OUTPUT
!
!       5       Redirected std. input (< pathname)   X
!       6                                                       X
!      10                                                       X
!      11                                            X
!************************************************************************************
!      NSTA =     NUMBER OF STATIONS TO PREDICT
!      CONV =     FACTOR FOR CONVERTING PREDICTED TIME SERIES TO NEW UNITS
!**** Conversion options for time ***
!
!     TCONV = 0   NO CONVERSION OF PREDICTED TIMES
!     TCONV = TIME MERIDIAN FOR WHICH THE KAPPA PRIMES IN THE HARMONIC
!             CONSTANTS WERE DERIVED. THIS OPTION IS USED TO CONVERT
!             THE TIMES TO GREENWICH IF THE CONSTANTS WERE CALCULATED
!             FOR A LOCAL TIME MERIDIAN. A REASON FOR USING THIS OPTION
!             IS IF COMPARISONS WITH ACTUAL DATA IS REQUIRED WHICH IS
!             IN GREENWICH TIME AND THE HARMONIC CONSTANTS WERE OBTAINED
!             FROM THE PREDICTION BRANCH WHERE HARMONIC CONSTANTS ARE
!             FOR LOCAL MERIDIANS ALWAYS
!
!     TCONV changed to hours to shift predicted time series 
!     (positive = later) --- Chris Zervas (7/97)
!
!**** conversion option for using 2mn2 in the harmonic constants  
!     file versus the standard L2
!
!     ***note*** it is node (33) which is re-calculated
!
!     IL2 = 0 --- use the standard <L2> harmonic constants
!           1 --- use the <2MN2> harmonic constants
!
!
!*****************************************************************
!
!     **NOTE(1)** USE FORMAT NO. 531 FOR HARMONIC CONSTANTS
!        TEMPORARY FORMULATION FORMAT 532 IS USED
!
!******************************************************************
!
!
!     IPREDK = 0 -- USED TO CALCULATE THE DIFFERENCES BETWEEN A SET OF
!                   PREDICTED AND OBSERVED SERIES. THIS CALCULATION IS
!                   DONE USING AN INPUT OF OBSERVED VALUES IN ASCII
!                   FORMAT AND OUTPUT IS WRITTEN IN AN ASCII FILE
!
!     IPREDK = 1 -- USED TO CALCULATE THE DIFFERENCES BETWEEN A SET OF
!                   PREDICTED AND OBSERVED SERIES. THIS CALCULATION IS
!                   DONE USING AN INPUT OF OBSERVED VALUES IN CDF
!                   FORMAT AND OUTPUT IS WRITTEN IN AN ASCII FILE
!
!     **NOTE(1)** STARTING TIME FOR PREDICTIONS IS SET EQUAL TO
!                 THE STARTING TIME DEFINED BY THE VALUE OF ISTART
!     **NOTE(2)** NOS1 IS CALCULATED BY THE TIMES OF THE OBSERVED
!                 SERIES
!
!
!      IPRED = 2 -- USED TO CALCULATE A PREDICTED SERIES
!
!     **NOTE(3)** THE PREDICTIONS AND DIFFERENCES CAN NOT BE PERFORMED
!                 ACROSS 2 YEARS
!     **NOTE(4)** THE YEARS FOR EACH NSTA MUST BE THE SAME
!     **NOTE(5)** USE FORMAT NO. 532 FOR HARMONIC CONSTANTS
! ***************************** modified IEL to match KINDAT
!     IEL =    THE ELEMENT OF THE DATA SERIES TO PERFORM CALCULATIONS
!              1   MAJOR/MINOR COMPONENTS OF VECTOR VARIABLE (I.E. CURRENT)
!              2   SCALAR VARIABLE (I.E. TIDAL HEIGHT)
!              3   TEMPERATURE (CDF INPUT FIELD)
!              4   CONDUCTIVITY (CDF INPUT FIELD)
!              6   PRESSURE (CDF INPUT FIELD)
!     IYRS=YEAR OF THE FIRST DATA POINT(CAN NOT GO ACROSS YRS)
!     IMMS=MONTH OF FIRST DATA POINT
!     IDDS=DAY OF FIRST DATA POINT
!     TIMEF=TIME OF FIRST DATA POINT
!     MON=MONTH OF LAST DATA POINT
!     IDDE=DAY OF LAST DATA POINT
!     TIMEL=TIME OF LAST DATA POINT
!     DELT= desired time interval in hours, delt=0.1 for 6 minutes data
!     NOS1=NUMBER OF DATA POINTS PER HOUR
!     XMAJOR ---- AXIS FOR MAJOR/MINOR COMPONENTS (MAKE SURE HARMONIC
!              CONSTANTS WERE DERIVED ALONG THIS AXIS)
!                       IF A 0.0 IS READ --- 0 DEGRESS TRUE IS
!                       ASSUMED. ALWAYS READ IN THE CONSTITUENTS
!                       FOR THE MAJOR AXIS FIRST.
!///////////////////////////////////////////////////////////////
!
! 
! 
      PARAMETER (MXDIM=40000)
      PARAMETER (XMISS=999.)
      character*80 cdfin,cdfout,BUFFER*100
      CHARACTER*10   ALIST
      CHARACTER*80   HEAD(2)
      DIMENSION   XDATA(14,50),ALIST(37),IHEAD(34),DAT(6),SUM(6),SMN(6),
     #TIM(MXDIM),SPEED(MXDIM),DIREC(MXDIM),STORX(MXDIM)
      common/ptide/amp(37),epoc(37),storn(mxdim),ptime(mxdim)

      common/virt1/
     1            A(37),XODE(114),VPU(114),
     2            XCOS(1025),ARG(37),TABHR(24),ANG(37),SPD0(37),
     3            EPOCH(37),AMPA(37),IYR(15),NUM(15),ISTA(6),NO(6),
     4            JODAYS(12),C(37)
!
      real*8 jday,jday0,jday1,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      COMMON /XDATA/XDATA
      DIMENSION TSTART(4)
!CCCCCCCCCCCCC
      common /speeds/spd
      common /names/lable
      real*8 spd(180),fff(180),vau(180)
      character*10 lable(180)
      integer order(37)
      data order/1,3,2,5,21,4,22,28,24,30,13,25,12,8,16,11,27,15,14,
     1  35,37,36,34,33,20,18,10,9,19,17,32,26,7,29,6,23,31/

!CCCCCCCCCCCCC
      REAL NOS1
!     DATA (ALIST(I),I=1,37) /'M(2)','S(2)','N(2)','K(1)','M(4)','O(1)',
!    1'M(6)','MK(3)','S(4)','MN(4)','NU(2)','S(6)','MU(2)','2N(2)','OO(1
!    2)','LAMDA(2)','S(1)','M(1)','J(1)','MM','SSA','SA','MSF','MF','RHO
!    3(1)','Q(1)','T(2)','R(2)','2Q(1)','P(1)','2SM(2)','M(3)','L(2)','2
!    4MK(3)','K(2)','M(8)','MS(4)'/
      DATA (ALIST(I),I=1,37) /'M(2)      ','S(2)      ','N(2)      ',
     1                        'K(1)      ','M(4)      ','O(1)      ',
     2                        'M(6)      ','MK(3)     ','S(4)      ',
     3                        'MN(4)     ','NU(2)     ','S(6)      ',
     4                        'MU(2)     ','2N(2)     ','OO(1)     ',
     5                        'LAMDA(2)  ','S(1)      ','M(1)      ',
     6                        'J(1)      ','MM        ','SSA       ',
     7                        'SA        ','MSF       ','MF        ',
     8                        'RHO(1)    ','Q(1)      ','T(2)      ',
     9                        'R(2)      ','2Q(1)     ','P(1)      ',
     1                        '2SM(2)    ','M(3)      ','L(2)      ',
     2                        '2MK3(3)   ','K(2)      ','M(8)      ',
     3                        'MS(4)     '/
!      write(*,*) IYRS,IMMS,IDDS,IHHS,MNS,NDAY,DELT,DATUM,
!     &     amp,epoc
      LIN = 5
      LOUT = 6
!
!      read (lin,*)   nsta,ipredk,conv,tconv,il2
      NSTA=1
      IPREDK=2
      CONV=1.0
      TCONV=0.0
      IL2=0
!
      IF (IPREDK.GT.2.OR.IPREDK.LT.0)   then
      print*,'Error in IPREDK value'
      stop
      endif
!
!     DEVELOP COSINE TABLE
!
      H = 0.00153398078789
      DO 35 I = 1,1024
      XCOS(I) = COS(H*(I-1))
   35 CONTINUE
      XCOS(1025) = 0.0
      ms0 = 1
      CON = 1024. / 90.
!        ---
      DO 612   JOB=1,NSTA
!      read(lin,16)cdfout
16    format(a80) 
6000  continue
!      READ (LIN,*) IEL,IMMS,IDDS,TIMEF,IMME,IDDE,TIMEL,DELT,IYRS,IYRE,
!     1  XMAJOR
!sfan      CALL GETARG(1,BUFFER)
!sfan      READ(BUFFER,*)IYRS,IMMS,IDDS,IHHS,MNS
!      write(*,*)'IYRS,IMMS,IDDS,IHHS,MNS=',IYRS,IMMS,IDDS,IHHS,MNS
      TIMEF=IHHS+MNS/60.
!sfan      CALL GETARG(2,BUFFER)
!sfan      READ(BUFFER,*)IYRE,IMME,IDDE,IHHE,MNE
!sfan      write(*,*)'IYRE,IMME,IDDE,IHHE,MNE=',IYRE,IMME,IDDE,IHHE,MNE
!sfan      TIMEL=IHHE+MNE/60.
!sfan      CALL GETARG(3,BUFFER)
!sfan      READ(BUFFER,*)IEL
      iel=2
!      write(*,*)'iel=',IEL
!sfan      CALL GETARG(4,BUFFER)
!sfan      READ(BUFFER,*)DELT
!      write(*,*)'DELT',DELT
!      CALL GETARG(5,BUFFER)
!      READ(BUFFER,*)XMAJOR  !! This value will be overwritten by reading from harmonic data file
!sfan      CALL GETARG(5,cdfin)
!sfan      CALL GETARG(6,cdfout)
!sfan      call ncrght(cdfin,nct)
!sfan      open(30,file=cdfin(1:nct) )
!sfan: base time
      yearb=1940.
      monthb=1.
      dayb=1.
      hourb=0.
      jday1=JULIAN(yearb,monthb,dayb,hourb)
      yearb=IYRS
      monthb=1.
      dayb=1.
      hourb=0.
      jbase_date=JULIAN(yearb,monthb,dayb,hourb)
      dfday=jbase_date-jday1
!      print *,'dfday=',dfday
!      stop
      NOS1=1.0/DELT   !!!  data points per hour
  527 FORMAT(2I5,F10.0,2I5,F10.0)
!sfan      write(6,*)'run pred.f from ',IYRS,IMMS,IDDS,IHHS, ' to ',
!sfan     1   IYRE,IMME,IDDE,IHHE
!      write(6,*)'run pred.f from ',IYRS,IMMS,IDDS,IHHS, ' for ',
!     1   nday, ' days'
!      CALL COMPIN(MXDIM,NOS1,IYRS,IMMS,IDDS,TIMEF,IYRE,IMME,IDDE,TIMEL,
!     1DELT,NPTS,TSTART)
      yearb=IYRS
      monthb=IMMS
      dayb=IDDS
      hourb=IHHS+MNS/60.0
      jday0=JULIAN(yearb,monthb,dayb,hourb)
!sfan      yearb=IYRE
!sfan      monthb=IMME
!sfan      dayb=IDDE
!sfan      hourb=IHHE+MNE/60.0
!sfan      jday1=JULIAN(yearb,monthb,dayb,hourb)
!NG      jday1=jday0+nday
      jday1=jday0+nday+1 ! CNG04182008 one extra day for precision issues
      TSTART(4) = IYRS
      TSTART(3) = jday0-jbase_date
      TSTART(2) = IHHS
      TSTART(1) =MNS
      NPTS=INT((jday1-jday0)*24/DELT+1+0.1)
!****
      IZAJ=0
119   CIND = NOS1

111   continue
!sfan      READ (30,550)  HEAD(1),HEAD(2)
!sfan      READ (30,532)DATUM,ISTA(1),NO(1),(AMP(J),EPOC(J),J=1,7),
!sfan     1 ISTA(2),NO(2),
!sfan     2 (AMP(J),EPOC(J),J=8,14),ISTA(3),NO(3),(AMP(J),EPOC(J),J=15,21),
!sfan     3 ISTA(4),NO(4),(AMP(J),EPOC(J),J=22,28),ISTA(5),NO(5),(AMP(J),
!sfan     4 EPOC(J),J=29,35),ISTA(6),NO(6),(AMP(J),EPOC(J),J=36,37)
!sfan      write (*,550)  HEAD(1),HEAD(2)
!      write (*,*)DATUM,(AMP(J),EPOC(J),J=1,37)
!  read in principle durrent direction from harmonic constants data file
!sfan      IF( (IZAJ .EQ. 0) .and. (IEL .eq. 1) )THEN
!sfan        Ind=INDEX(HEAD(2),'along')+5
!sfan        read(HEAD(2)(ind:ind+4),'(I4)')IPCD
!sfan       XMAJOR=IPCD
!sfan       IZAJ=1
!sfan     ENDIF
!sfan  112 DO 113 L = 1,5
!sfan      IF (ISTA(L).NE.ISTA(L+1)) GO TO 451
!sfan  113 CONTINUE
!sfan      DO 114 L = 1,6
!sfan      IF (NO(L).NE.L) GO TO 450
!sfan  114 CONTINUE
!
!     CONVERT CONSTANTS IF TCONV IS NOT EQUAL TO ZERO
!     TCONV IS THE TIME MERIDIAN
!
      IF (TCONV.EQ.0.0)   GO TO 120
!
!sfan      PRINT 5902
!sfan      PRINT 550, HEAD(1),HEAD(2)
!sfan      PRINT 5902
      WRITE (LOUT,994)TCONV
  994 FORMAT (' Values of the Epochs before ',F6.2,' hour time shift')
!    1       ' THE VALUES OF THE EPOCHS BEFORE THE CONVERSION TO',
!    2       ' GREENWICH TIME'/)
!    3       ' THE VALUE OF THE LOCAL TIME MERIDIAN OF THE ORIGINAL'/
!    4       ' VALUES OF THE EPOCHS (OR THE VALUES OF KAPPA PRIME)',
!    5       ' WAS ',F10.2/)
!
!      PRINT 5902
!      PRINT 5007
!      PRINT 5902
      DO 650   J=1,37
      IF (AMP(J).EQ.0.0)   GO TO 651
!      PRINT 5000, ALIST(J),AMP(J),EPOC(J)
!
      EPOC(J) = EPOC(J) + A(J)*TCONV
  652 IF (EPOC(J).LE.360.)   GO TO 650
      EPOC(J) = EPOC(J) - 360.0
!
      GO TO 652
!
  651 continue    !!  PRINT 998, ALIST(J)
!
  650 CONTINUE

 120  if(ms0.eq.2) goto 150
!CCCCCCC  using equarg replace yrcrds.f
      length=int(jday1-jday0)+1
      call equarg (37,IYRS,1,1,365,lable(1),fff(1),vau(1))
!      call equarg (37,IYRS,IMMS,IDDS,length,lable(1),fff(1),vau(1))
      do j=1,37
        VPU(J)=VAU(ORDER(J))
        XODE(J)=FFF(ORDER(J))
!       WRITE(6,'(I2,2x,A10,1x,F12.7,2X,F5.1,2X,F6.4)')
!     1  J, ALIST(J),A(J),VPU(J),XODE(J)
      end do

!CCCCCCCCCCCCCC
  150 DO 155 J = 1,37
      C(J) = A(J) * (CON/CIND)
  155 CONTINUE
!     SET UP TABLES FOR NON-ZERO CONSTITUENTS
  160 K = 0
      DO 180  J = 1,37
  161 IF (AMP(J).EQ.0.0) GO TO 180
  170 K = K + 1
      AMPA(K) = AMP(J) * XODE(J)
      TEMX = VPU(J) - EPOC(J)
      IF (TEMX .GE. 0.) GO TO 171
      TEMX = TEMX + 360.
  171 EPOCH(K) = TEMX * CON
      SPD0(K) = C(J)
  180 CONTINUE
      NOCON = K
!
!     OPERATING TABLES NOW STORED AS AMPA(K),EPOCH(K),SPD0(K)
!
!**** CHECK LENGTH OF SERIES
!
  700 IF (NPTS.LE.MXDIM)   GO TO 191
      NPTS = MXDIM
      PRINT 100
      PRINT 1190, NPTS
!
  191 continue
!      WRITE (LOUT,1002)   NPTS
! 1002 FORMAT (/' Total number of prediction times = ',I10)
!
!**** DETERMINE FIRST HOUR OF TIME PERIOD  at 00:00 of Jan. 1, first=0.0
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC modified by zaj on May 13, 2004
      FIRST=(jday0-jbase_date)*24.* CIND
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc      
      DO 220   J=1,MXDIM
      STORX(J) = 0.
  220 CONTINUE
      KOUNT = 0
      KT = 0
!
!        ---
   43 DO 380 K = 1,NPTS
!        ---
!
!**** THE PREDICTION/TIME STEP IS STORED IN VARIABLE PRED
!
        PRED = 0.0
!
        IF (KOUNT.GT.0) GO TO 260
        KOUNT = 1
  231   IF (NOCON.EQ.0)   GO TO 375
        DO 250 J = 1,NOCON
          ARGU = SPD0(J) * FIRST + EPOCH(J)
          ARG(J) = AMOD(ARGU,4096.)
  250   CONTINUE
        GO TO 290
  260   IF (NOCON.EQ.0)   GO TO 375
        DO 280 J = 1,NOCON
          ARG(J) = ARG(J) + SPD0(J)
  270     IF (ARG(J).LT.4096.) GO TO 280
          ARG(J) = ARG(J) - 4096.
          GO TO 270
  280   CONTINUE
  290   DO 374 J = 1,NOCON
          IF (ARG(J) - 1024.)  320,320,300
  300     IF (ARG(J) - 2048.)  350,350,310
  310     IF (ARG(J) - 3072.)  360,360,330
  320     ANG(J) = ARG(J)
          GO TO 340
  330     ANG(J) = 4096. - ARG(J)
  340     NP = ANG(J) + 1.5
!     ----   ----   -----------------
          PRED = PRED + AMPA(J) * XCOS(NP)
!     ----   ----   ------------------
          GO TO 374
  350     ANG(J) = 2048. - ARG(J)
          GO TO 370
  360     ANG(J) = ARG(J) - 2048.
  370     NP = ANG(J) + 1.5
          PRED = PRED - AMPA(J) * XCOS(NP)
  374   CONTINUE
  375   IF (K.NE.NPTS) GO TO 381
        IF (IPREDK.EQ.0) GO TO 381
        IF (KT.EQ.1) GO TO 378
        FIRST = FIRST + NPTS - 1.
        KT = 1
        CHECK = PRED
        PRED = 0.0
        GO TO 231
  378   CKSUM = CHECK - PRED
!
!**** CONVERT RESULTS ACCORDING TO CONV
!
  381   PRED = (PRED + DATUM) * CONV
!
!**** STORE THE RESULTS
!
  621   STORX(K) = PRED
! ---
  380 CONTINUE
!      write(6,*)'it is ok'
!      stop
! ---
      NOS2 = 0
      IF(IEL.EQ.1) NOS2 = 1
      IF(IEL.EQ.1.AND.XMAJOR.NE.0) NOS2 = 2
!     PRINT 5902
!     PRINT 550, HEAD(1),HEAD(2)
!     PRINT 5902
!      PRINT 560, IYR1,DATUM,NOCON,CKSUM
      IF (NOS2.EQ.0) GO TO 5004
      IF(ms0.EQ.2) GO TO 4996
      IF(NOS2.EQ.1) GO TO 4995
      PRINT 5006
      GO TO 5004
 4995 PRINT 5008
      GO TO 5004
 4996 IF(NOS2.EQ.1) GO TO 5002
      PRINT 5011
      GO TO 5004
 5002 PRINT 5009
!NG04212008 5004 PRINT 5902
 5004 CONTINUE
!      IF(TCONV.EQ.0.0) PRINT 5007
!      PRINT 5902
      DO 5005 IZ =1,37
      IF(AMP(IZ).EQ.0.0) GO TO 4990
!      PRINT 5000,ALIST(IZ),AMP(IZ),EPOC(IZ)
      GO TO 5005
 4990 continue
  998 FORMAT (' ',A10)
 5005 CONTINUE
      IF (NOS2)400,400,395
  395 IF (ms0.EQ.2) GO TO 410
      ms0 = 2
!
!**** STORES THE MAJOR AXIS
!
  400 DO 1395 I = 1,NPTS
      STORN(I) = STORX(I)
 1395 CONTINUE
  624 IF (NOS2) 410,410,119
  410 CONTINUE
!
!**** FORM NEW DATA ARRAY
!
!sfan      call ncrght(cdfout,nct)
!sfan      open(20,file=cdfout(1:nct),status='unknown',form='formatted')
!      open(20,file='ptides.dat',status='unknown',form='formatted')
      idet = 0
      IPRED = 1
      jday=jday0
      
! ---
  427 CONTINUE
      call GREGORIAN(jday,yearb,monthb,dayb,hourb)
      IYEAR=INT(yearb)
      ICM=int(monthb+0.001)
      ICD=INT(dayb+0.001)
      IHRPT=INT(hourb+0.001)
      IMN=INT((hourb-IHRPT)*60+0.1)
      dayj=jday-jbase_date+1.+dfday

      ptime(ipred)=(jday-jday0)*24.0
      ncount=ipred
!
      IF(IEL.EQ.1) GO TO 418
!      write(20,571)dayj,IYEAR,ICM,ICD,IHRPT,IMN,storn(ipred)
!      write(10,572)IYEAR,ICM,ICD,IHRPT,IMN,ISEC,storn(ipred)
      go to 421
!
  418 continue 
      CALL VELDIR(STORN(IPRED),STORX(IPRED),DR,SP)
      DR = DR + XMAJOR
      IF (DR.GT.360.)   DR = DR - 360.
      u=sp*sin(dr*3.1415926/180.)
      v=sp*cos(dr*3.1415926/180.)
!     write(20,571)dayj,IYEAR,ICM,ICD,IHRPT,IMN,sp,dr,u,v
!      write(10,572)IYEAR,ICM,ICD,IHRPT,IMN,ISEC,sp,dr 
  571 format(f10.4,I5,4I3,4f10.4)
  572 format(i4,5i3.2,2f11.4)
  421 IPRED = IPRED + 1
      IF (IPRED.LT.NPTS)   GO TO 424
      IF (IPRED.EQ.NPTS)   GO TO 424
      GO TO 602
!           ---
424   CONTINUE      
      jday=jday+delt/24. 
      GO TO 427
!           ---
  602 ms0 = 1
  612 CONTINUE
! ---
!     WRITE (LOUT,1003)
 1003 FORMAT (' '///' ALL FILES CLOSED --- NORMAL TERMINATION')
!
!sfan      STOP
!sfan  450 PRINT 501
!sfan      STOP
!sfan  451 PRINT 502
!sfan      STOP
!sfan  452 PRINT 503
!sfan      STOP
!sfan  453 PRINT 504
!sfan      STOP
  500 FORMAT (5I5,F5.0)
  501 FORMAT (27H STATION CARDS OUT OF ORDER)
  502 FORMAT (31H STATION NUMBERS NOT CONSISTENT)
  503 FORMAT (28H YEAR NUMBERS NOT CONSISTENT)
  504 FORMAT (24H YEAR CARDS OUT OF ORDER)
  525 FORMAT (F12.8,2I2,F10.8)
  526 FORMAT (5I5,2F10.0,i2,f5.0)   
! 531 FORMAT (F6.3,6(/2I4,7(F5.2,F5.2)))
  532 FORMAT (F6.3,6(/2I4,7(F5.3,F4.1)))
  533 FORMAT (I4,2I2,F4.3,F4.1,F4.3,F4.1,F4.3,F4.1,F4.3,F4.1,
     1 F4.3,F4.1,F4.3,F4.1,F4.3,F4.1,F4.3,F4.1)
! 534 FORMAT (36I2)
! 537 FORMAT (9X,10F9.3)
! 538 FORMAT (3I4)
! 545 FORMAT (6F6.3)
! 546 FORMAT (10F6.3)
  550 FORMAT (A80)
! 551 FORMAT (9X,10F9.3)
! 552 FORMAT (1X//1X,15HPREDICTIONS FOR,2X,A9,2X,10HAT STATION,2X,I2)
! 553 FORMAT (1X//1X,15HPREDICTIONS FOR,2X,A9)
! 554 FORMAT (1X,I10,1X,I10,1X,18F6.3)
  560 FORMAT (' Year ',I4,'  Datum ',F6.3,'  No. of ',
     1'Constituents ',I4,'    checksum',F12.7/)
 2534 FORMAT (2I2,1X,2I2,1X,F5.2)
!   5 FORMAT (1X,2(I2,1X),I4,1X,F5.2,3X,I3,1X,F5.2,2X,F6.3,1X,F6.3)
!   6 FORMAT(24I3)
    7 FORMAT (1H1)
!   8 FORMAT (1X,46HLISTING OF PREDICTED NORTH AND EAST COMPONENTS/49H A
!    1ND THE CORRESPONDING VELOCITIES AND DIRECTIONS.)
!  14 FORMAT (4I5,F6.4)
!  15 FORMAT (8A10)
!  16 FORMAT (3I5,I7,I2,I10,I2,2I5)
  100 FORMAT (1H0)
! 101 FORMAT (4X,4HDATE,15X,5HSPEED,4X,10HCOMPONENTS)
! 213 FORMAT (1H+,50X,17HSPEED  IN KNOTS  /51X,25HDIRECTION IN DEGREES T
!    1RUE)
!8101 FORMAT(1X,16HDA MO YEAR  HOUR,3X,10HDIR (KNTS))
!8108 FORMAT(1H+,34X,1HN,6X,1HE/)
!8109 FORMAT (1H+,32X,5HMAJOR,2X,5HMINOR/)
 8110 FORMAT (1X/1X,31HERROR IN EITHER IPREDK OR ITYPE)
 8111 FORMAT (1X/1X,28HERROR IN EITHER IEL OR ITYPE)
 8112 FORMAT (1X/1X,20HEOF IN OBSERVED TAPE)
 8113 FORMAT (1X/1X,38HCalculated time of last data point is ,4F8.0)
 8114 FORMAT (1X/1X,38HCalculated time of first data point is,4F8.0)
 8115 FORMAT(1H1,5X,'NORMAL TERMINATION FOR IPREDK=2. THERE IS 1 FILE
     1IN ',a27)
!8201 FORMAT (2F5.3)
! 900 FORMAT (30X,F6.3,2X,F6.3)
! 901 FORMAT (30X,15HMEANS OF SERIES/33X,1HN,7X,1HE)
! 902 FORMAT(30X,27H(I.E. NONTIDAL  COMPONENTS))
! 903 FORMAT (30X,15HMEANS OF SERIES/33X,3HMAJ,5X,3HMIN)
 5902 FORMAT (1X,4H    )
!8102 FORMAT (4X,4HDATE,12X,10HMAJOR AXIS)
!8103 FORMAT (1X,16HDA MO YEAR  HOUR,3X,3HDIR,2X,5HSPEED)
!3901 FORMAT (27X,14HMEAN OF SERIES)
 5000 FORMAT (1X,A10,2X,F7.3,3X,F6.2)
 5006 FORMAT(1X,38HHarmonic Constants (Major Axis) ------)
 5011 FORMAT(1X,38HHarmonic Constants (Minor Axis) ------)
 5007 FORMAT (1X,11HConstituent,11X,5HKappa/16X,4HH(A),3X,5HPrime)
 5008 FORMAT (1X,26HHarmonic Constants (North))
 5009 FORMAT (1X,25HHarmonic Constants (East))
!5010 FORMAT (1X,16HPREDICTED VALUES)
!3224 FORMAT (47X,16HNONTIDAL CURRENT/51X,3HVEL,4X,3HDIR)
!3225 FORMAT (49X,F5.2,4X,F3.0)
!3223 FORMAT (1X,7HTFAC = ,F10.8)
 1190 FORMAT (' LENGTH OF PREDICTIONS SHORTENED TO',I10)
!      close(30)
!      close(20)

      return
      END

      SUBROUTINE VELDIR (AAVN,AAVE,AANGLE,AVEL)
      RADDEG= 57.295779513082
      IF(AAVN.NE.0.0) GO TO 113
      IF (AAVE) 110,111,112
  110 AANGLE = 270.0
      GO TO 117
  111 AANGLE = 000.0
      GO TO 117
  112 AANGLE = 090.0
      GO TO 117
  113 AANGLE = ATAN (AAVE/AAVN)
      AANGLE = AANGLE * RADDEG
      IF ( AAVN.GT.0.0) GO TO  114
      AANGLE = AANGLE + 180.0
      GO TO 117
  114 IF (AAVE.GT.0.0) GO TO 117
      AANGLE = AANGLE  + 360.0
  117 AVEL = SQRT (AAVE**2 + AAVN**2)
      RETURN
      END
      BLOCK DATA
      PARAMETER (MXDIM=200000)
!CCCCCCCCCCCCC
      common /speeds/spd
      common /names/lable
      common /mmss/ms
      character*10 lable(180)
      double precision spd(180)
      integer ms(180)
!CCCCCCCCCCCCC
      common/datec/idtbc(12,2),iltbc(12,2)
      common/datej/idtbj(12),iltbj(12)
      common/virt1/
     1            A(37),XODE(114),VPU(114),
     2            XCOS(1025),ARG(37),TABHR(24),ANG(37),SPD0(37),
     3            EPOCH(37),AMPA(37),IYR(15),NUM(15),ISTA(6),NO(6),
     4            JODAYS(12),C(37)
      DATA ((IDTBC(I,J),J=1,2),I=1,12)/1,31,32,59,60,90,91,120,121,151,
     1                                  152,181,182,212,213,243,244,273,
     2                                  274,304,305,334,335,365/
      DATA  ((ILTBC(I,J),J=1,2),I=1,12)/1,31,32,60,61,91,92,121,122,
     1                              152,153,182,183,213,214,244,245,
     2                                    274,275,305,306,335,336,366/
      DATA (IDTBJ(I),I=1,12)/1,32,60,91,121,152,182,213,244,274,305,335/
      DATA (ILTBJ(I),I=1,12)/1,32,61,92,122,153,183,214,245,275,306,336/
      DATA(TABHR(I), I=1,24)/  -24.,  720., 1392., 2136., 2856., 3600.,
     1    4320., 5064., 5808., 6528., 7272., 7992.,  -24.,  720., 1416.,
     2    2160., 2880., 3624., 4344., 5088., 5832., 6552., 7296., 8016./
      DATA (A(I), I=1,37)/        28.9841042,  30.0000000,  28.4397295,
     1   15.0410686,  57.9682084,  13.9430356,  86.9523127,  44.0251729,
     2   60.0000000,  57.4238337,  28.5125831,  90.0000000,  27.9682084,
     3   27.8953548,  16.1391017,  29.4556253,  15.0000000,  14.4966939,
     4   15.5854433,   0.5443747,   0.0821373,   0.0410686,   1.0158958,
     5    1.0980331,  13.4715145,  13.3986609,  29.9589333,  30.0410667,
     6   12.8542862,  14.9589314,  31.0158958,  43.4761563,  29.5284789,
     7   42.9271398,  30.0821373, 115.9364169,  58.9841042/
      DATA (JODAYS(LL),LL=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
!CCCCCCCCCCCCCCCCCCCCCCC
      data(spd(i),i=1,37)/28.9841042d0,   28.4397295d0,   30.0000000d0,
     1    13.9430356d0,   15.0410686d0,   30.0821373d0,   29.5284789d0,
     2    27.8953548d0,   30.0410667d0,   29.9589333d0,   29.4556253d0,
     3    27.9682084d0,   28.5125831d0,   15.5854433d0,   14.4966939d0,
     4    16.1391017d0,   14.9589314d0,   13.3986609d0,   12.8542862d0,
     5    13.4715145d0,   57.9682084d0,   86.9523127d0,  115.9364169d0,
     6    60.0000000d0,   90.0000000d0,   43.4761563d0,   15.0000000d0,
     7    44.0251729d0,   42.9271398d0,   57.4238337d0,   58.9841042d0,
     8    31.0158958d0,   01.0980331d0,   01.0158958d0,   00.5443747d0,
     9     0.0410686d0,    0.0821373d0/
      data(spd(i),i= 38,48)/12.9271398d0, 14.0251729d0,   14.5695476d0,
     a    15.9748272d0,   16.0569644d0,   30.5443747d0,   27.4238337d0,
     b    28.9019669d0,   29.0662415d0,   26.8794590d0,   26.9523126d0/
      data(spd(i),i= 49,87)/27.4966873d0, 31.0980331d0,   27.8039338d0,
     c    28.5947204d0,   29.1483788d0,   29.3734880d0,   30.7086493d0,
     d    43.9430356d0,   45.0410686d0,   42.3827651d0,   59.0662415d0,
     e    58.4397295d0,   57.4966873d0,   56.9523127d0,   58.5125831d0,
     f    56.8794590d0,   59.5284789d0,   71.3668693d0,   71.9112440d0,
     g    73.0092770d0,   74.0251728d0,   74.1073100d0,   72.9271398d0,
     h    71.9933813d0,   72.4649023d0,   88.9841042d0,   86.4079380d0,
     i    87.4238337d0,   87.9682084d0,   85.3920421d0,   85.8635632d0,
     j    88.5125831d0,   87.4966873d0,   89.0662415d0,   85.9364168d0,
     k    86.4807916d0,   88.0503457d0,  100.3509735d0,  100.9046318d0/
      data(spd(i),i= 88,124)/101.9112440d0,103.0092771d0,116.4079380d0,
     l   116.9523127d0,  117.9682084d0,  114.8476674d0,  115.3920422d0,
     m   117.4966873d0,  115.4648958d0,  116.4807916d0,  117.0344500d0,
     n   118.0503457d0,  129.8887360d0,  130.4331108d0,  130.9774855d0,
     o   131.9933813d0,  144.3761464d0,  144.9205211d0,  145.3920422d0,
     p   145.9364169d0,  146.4807916d0,  146.9523127d0,  160.9774855d0,
     q   174.3761464d0,  174.9205211d0,  175.4648958d0,  175.9364169d0,
     r    14.9178647d0,   15.0821353d0,   15.1232059d0,  15.5125897d0,
     s    30.6265119d0,   27.3416965d0,   42.9271397d0,   60.0821373d0,
     t    16.1391016d0,   12.8450026d0/
      data(spd(i),i=125,140)/ 26.9615963d0, 27.5059710d0, 28.6040041d0,
     u    57.5059710d0,   58.5218668d0,   85.4013258d0,   85.9457005d0,
     v    86.4900752d0,   87.5059710d0,   88.5218668d0,  115.4741794d0,
     w   116.4900752d0,  117.5059710d0,  146.4900752d0,  175.4741794d0,
     x    41.9205276d0/
      data(spd(i),i=141,150)/ 15.1232058d0, 14.8767942d0, 30.0000001d0,
     y    29.9178627d0,   30.1642746d0,   29.9178666d0,   59.9589333d0,
     z    59.9178627d0,   60.2464119d0,   59.8767999d0/
      data(spd(i),i=151,164)/ 28.9430356d0, 01.0569644d0, 00.5490165d0,
     1          00.5079479d0, 00.0410667d0, 00.1232059d0, 00.1642746d0,
     2          00.2464118d0, 00.3285841d0, 00.4106864d0, 00.4928237d0,
     3          00.9856473d0, 45.0000000d0, 75.0000000d0/
      data(spd(i),i=165,175)/ 27.8860712d0, 30.0410686d0, 43.4807981d0,
     4          44.9589314d0, 45.1232059d0, 56.3258007d0, 56.8701754d0,
     5          57.8860712d0,105.0000000d0,120.0000000d0,150.0000000d0/

      data(ms(n),n=1,  37)/3*2,2*1,8*2,7*1,4,6,8,4,6,3,1,2*3,2*4,2,5*0/
      data(ms(n),n=38, 48)/ 5*1,6*2/
      data(ms(n),n=49, 96)/7*2,3*3,7*4,8*5,12*6,4*7,7*8/
      data(ms(n),n=97, 124)/3*8,4*9,6*10,11,4*12,4*1,2*2,3,4,2*1/
      data(ms(n),n=125,140)/3*2,2*4,5*6,3*8,10,12,3/
      data(ms(n),n=141,150)/2*1,4*2,4*4/
      data(ms(n),n=151,164)/2,11*0,3,5/
      data(ms(n),n=165,175)/2*2,3*3,3*4,7,8,10/

      data(lable(m),m = 1,37)/ 'M(2)      ','N(2)      ','S(2)      ',
     1            'O(1)      ','K(1)      ','K(2)      ','L(2)      ',
     2            '2N(2)     ','R(2)      ','T(2)      ','Lambda(2) ',
     3            'Mu(2)     ','Nu(2)     ','J(1)      ','M(1)      ',
     4            'OO(1)     ','P(1)      ','Q(1)      ','2Q(1)     ',
     5            'Rho(1)    ','M(4)      ','M(6)      ','M(8)      ',
     6            'S(4)      ','S(6)      ','M(3)      ','S(1)      ',
     7            'MK(3)     ','2MK(3)    ','MN(4)     ','MS(4)     ',
     8            '2SM(2)    ','Mf        ','Msf       ','Mm        ',
     9            'Sa        ','Ssa       '/
      data lable(38)/'SK(3)     '/
      data(lable(m),m =39,77)/ 'Sigma(1)  ','MP(1)     ','Chi(1)    ',
     A            '2PO(1)    ','SO(1)     ','2NS(2)    ','2NK2S(2)  ',
     B            'MNS(2)    ','MNK2S(2)  ','2MS2K(2)  ','2KN2S(2)  ',
     C            'OP(2)     ','MKS(2)    ','M2(KS)(2) ','2SN(MK)(2)',
     D            'MSN(2)    ','2KM(SN)(2)','SKM(2)    ','NO(3)     ',
     E            'SO(3)     ','N(4)      ','3MS(4)    ','MNKS(4)   ',
     F            'SN(4)     ','KN(4)     ','MK(4)     ','SL(4)     ',
     G            'MNO(5)    ','2MO(5)    ','3MP(5)    ','MNK(5)    ',
     H            '2MP(5)    ','2MK(5)    ','MSK(5)    ','3KM(5)    ',
     I            '3NKS(6)   ','2NM(6)    ','2NMKS(6)  ','2MN(6)    '/
      data(lable(m),m=78,114)/ '2MNKS(6)  ','MSN(6)    ','MKN(6)    ',
     J            '2MS(6)    ','2MK(6)    ','NSK(6)    ','2SM(6)    ',
     K            'MSK(6)    ','2MNO(7)   ','2NMK(7)   ','2MSO(7)   ',
     L            'MSKO(7)   ','2(MN)(8)  ','3MN(8)    ','3MNKS(8)  ',
     M            '2MSN(8)   ','2MNK(8)   ','3MS(8)    ','3MK(8)    ',
     N            'MSNK(8)   ','2(MS)(8)  ','2MSK(8)   ','2M2NK(9)  ',
     O            '3MNK(9)   ','4MK(9)    ','3MSK(9)   ','4MN(10)   ',
     P            'M(10)     ','3MNS(10)  ','4MS(10)   ','2MNSK(10) ',
     Q            '3M2S(10)  ','4MSK(11)  ','4MNS(12)  ','5MS(12)   ',
     R            '3MNKS(12) ','4M2S(12)  '/
      data(lable(m),m=115,120)/'2NP(3)    ','ML(4)     ','2ML(6)    ',
     1            'MSL(6)    ','2MSL(8)   ','3ML(8)    '/
      data(lable(m),m=121,128)/'TK(1)     ','RP(1)     ','KP(1)     ',
     S            'THETA(1)  ','KJ(2)     ','OO(2)     ','MO(3)     ',
     T            'SK(4)     '/
      data(lable(m),m=129,138)/'MLN2S(2)  ','2ML2S(2)  ','MKL2S(2)  ',
     U            '2MLS(4)   ','2NMLS(6)  ','2MLNS(6)  ','3MLS(6)   ',
     V            '4MLS(8)   ','3MSL(10)  ','4MSL(12)  '/
      data(lable(m),m=139,140)/'2KO(1)    ','2OK(1)    '/
      data(lable(m),m=141,150)/'2KP(1)    ','2PK(1)    ','KP(2)     ',
     W            '2SK(2)    ','2KS(2)    ','2TS(2)    ','ST(4)     ',
     X            '3SK(4)    ','3KS(4)    ','3TS(4)    '/
      data(lable(m),m=151,164)/'SO(2)     ','SO(0)     ','.5Mf      ',
     1            '.5Msf     ','ST(0)     ','3Sa       ','4Sa       ',
     2            '6Sa       ','8Sa       ','10Sa      ','12Sa      ',
     3            '24Sa      ','S(3)      ','S(5)      '/
      data(lable(m),m=165,175)/'O(2)      ','SK(2)     ','NK(3)     ',
     4            'SP(3)     ','K(3)      ','NO(4)     ','MO(4)     ',
     5            'SO(4)     ','S(7)      ','S(8)      ','S(10)     '/

      end
       subroutine ntran(iunit,ifunc,ipar1,ipar2,ipar3,ipar4)                    
!                                                                               
        dimension ipar2(700)                                                    
        ipar3=0                                                                 
        go to (1,2,99,99,99,99,7,8,9,10,99,99,99,99,99,99,99,99,99,99,          
     x  99,99,99),ifunc                                                         
!                                                                               
        return                                                                  
1       write(iunit,err=100)(ipar2(i),i=1,ipar1)                                
        return                                                                  
!                                                                               
!****       read function                                                           
!                                                                               
2       read(iunit,end=101,err=100)(ipar2(i),i=1,ipar1)                         
        return                                                                  
!                                                                               
!       move block/record routine                                               
!                                                                               
7       isave=iabs(ipar1)                                                       
        if(ipar1.lt.0)then                                                      
        do 70 i=1,isave                                                         
        backspace iunit                                                         
70      continue                                                                
!                                                                               
        elseif(ipar1.gt.0)then                                                  
        do 71 i=1,isave                                                         
        read(iunit,err=100,end=101)idum                                         
71      continue                                                                
!                                                                               
        endif                                                                   
        return                                                                  
!                                                                               
!       routine to move files                                                   
!                                                                               
8       isave=iabs(ipar1)                                                       
        if(ipar1.lt.0)then                                                      
        do 80 i=1,isave                                                         
        go to 100                                                               
80      continue                                                                
!                                                                               
        elseif(ipar1.gt.0)then                                                  
        do 81 i=1,isave                                                         
82       read(iunit,end=81,err=100)idum                                         
        go to 82                                                                
81      continue                                                                
!                                                                               
        endif                                                                   
        return                                                                  
!                                                                               
!       routine to write an eof on unit                                         
!                                                                               
9       endfile iunit                                                           
        return                                                                  
!                                                                               
!       routine to rewind the unit                                              
!                                                                               
10        rewind iunit                                                          
        return                                                                  
!                                                                               
!                                                                               
!       catch all else                                                          
!                                                                               
99      return                                                                  
100      ipar3=-3                                                               
      print*,' !!!!!!!!!!! read err= ',ipar3
         return                                                                 
101      ipar3=-2                                                               
         return                                                                 
!                                                                               
        end                                                                     
!
!     program to calculate equilibrium arguments and node factors
!     by e.e.long      (slight revision by b.b.parker)
!       major revisions by len hickman 6/11/86 to make this a subroutine
!       of lsqha.  v is calculated for the beginning of the series.
!       u and f are adjusted to the midpoint of the series (regardless
!       of length).
*     More revisions by Geoff French 9/1/86
!
      subroutine equarg (nsped,IYR,ICM,ICD,length,label,fff,vau)
      implicit real*8(a-h,o-z)
!      implicit none
      integer nsped
      integer iyr, icm, icd, length
      character*10 label(180)
      real*8 vau(180),fff(180)
      real*8 jday,jday0,jday1,jbase_date,JULIAN,yearb,monthb,dayb,hourb
      real*8 speed
      character*10 lname
      real*8 spd(180),vuu(180),vou(180)
      common/locat/tm,gonl
      common/costx/cxx(30),oex(5)
      common/fad/ipick
      common/vee/tml,con,u,q,ui
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/boxs/aw,ai,ae,ae1,asp
      real*8 xyer
      integer iyer
      real*8 daym
      real*8 grbs,grms
      real*8 grbss, grmss
      real*8 tm
      real*8 gonl
      real*8 tml
      integer i,j
      integer isub, inum
      real*8 e, f
      real*8 daybb, daymm
      integer juday
      integer iz
      integer mt
      real*8 ae1, ai, ae, asp
      real*8 aw, cra
!      real*8 cxx(30), oex(5)
*     ******************************************************************
*     *                                                                *
*     *   nsped = number of constituents to be calculated              *
*     *                                                                *
*     *                                                                *
*     *     ICM  = month of first data point                        *
*     *     ICD    = day  of first data point                         *
*     *     IYR   = year of first data point                         *
*     *     length  = length of series to be analyzed (in days)        *
*     *                                                                *
*     ******************************************************************
      yearb=IYR
      monthb=1.
      dayb=1.
      hourb=0.
      jbase_date=JULIAN(yearb,monthb,dayb,hourb)

      iyer =IYR
      dayb = 0.0
      daym = 183.0
      grbs = 0.0
      grms = 12.0

      if (mod(iyer,4) .eq. 0) then
        daym = 184.0
        grms = 0.0
      end if

      xyer = iyer

*     compute basic astronimical constants for the time period wanted
      call astro(xyer,dayb,daym,grbs,grms)

      tm = 0.0
      gonl = 0.0
      tml = 0.0

*     look up constituent parameters by matching the name
      do 100 j = 1,nsped
      lname = label(j)
      speed = 0.0
      call name (speed,lname,isub,inum,2)
      spd(j) = speed
      call vanduf (speed,e,f,2)
      fff(j) = f
      vou(j) = e
  100 continue

!      month = ICM
!      iday = ICD
!60      call CONCTJ(juday,month,iday,iyer)
60    yearb=IYR
      monthb=ICM
      dayb=ICD
      hourb=0.
      jday=JULIAN(yearb,monthb,dayb,hourb)
      juday=jday-jbase_date+1
      do 200 j = 1,nsped
      vuu(j) = vou(j) + spd(j)*(float(juday-1)*24.00)
  200 continue
      call twopi(vuu,nsped)
      daybb = juday
      daymm = juday
      grbss = 0.0
      grmss = 0.0
      daymm = juday + length/2

      call astro (xyer,daybb,daymm,grbss,grmss)

      do 277 iz = 1,nsped
      speed = spd(iz)
      call vanduf(speed,e,f,2)
      vau(iz) = e
  277 fff(iz) = f

*     round node to 3 decimals, vo+u to 1 place
      do 222 mt = 1,nsped
      fff(mt) = anint(fff(mt)*1000.) * 0.001
      vau(mt) = anint(vau(mt)*10.) * 0.1
  222 continue

      return
      end
************************************************************************
      subroutine relate (aug,test,result)
      if (aug.ge.test) then
        result = aug - test
      else
        result = aug + test
      end if
      end
************************************************************************
      subroutine fitan( aus,auc,rta,spdx,jmap)
      implicit none
      real*8 aus, auc, rta, spdx
      integer jmap
      real*8 bxg
      if(auc.eq.0.0) go to 14
      bxg =  aus/auc
      go to 15
   14 if(aus)16,55,17
   16 rta = 270.0
      go to 88
   17 rta = 90.0
      go  to 88
   15 rta = atan(bxg)*57.2957795
      go to (88,38),jmap
   38 if(aus)32,32,34
   32 if(auc)33,55,37
   33 rta = 180.0 + rta
      go to 88
   34 if(auc)35,17,88
   35 rta = rta + 180.0
      go to 88
   37 rta = rta + 360.0
      go to 88
   55 rta = 0.0
   88 continue
      spdx = 0.0
      return
      end
************************************************************************
      subroutine twopi (aug,io)
      integer io
      real*8 aug(io)
      do 114 mo = 1,io
      zat = aug(mo)/360.0
      if(zat) 77,88,88
   77 aug(mo) = ((zat - aint(zat)) + 1.00)*360.0
      go to 114
   88 aug(mo) = (zat - aint(zat))*360.0
  114 continue
      return
      end
************************************************************************
      subroutine astro (xyer,dayb,daym,grbs,grms)
      implicit real*8(a-h,o-z)
      common/locat/tm,gonl
      common/costx/cxx(30),oex(5)
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/vee/tml,con,u,q,ui
      common/boxs/aw,ai,ae,ae1,asp
      common/boxxs/vib,vb,xib,vpb,vppb,cxsb,cxpb,cxhb,cxp1b

      pinv = 57.29578
      nyear = int(xyer)
      call orbit(xcen,xsx,xpx,xhx,xp1x,xnx,oex,t,xyer,5)
      xw = 23.4522944 - .0130125*t - .00000164*t**2 + .000000503*t**3
      xi = 5.14537628
      aw = xw*0.0174533
      ai = xi*0.0174533
      ae = 0.0548997
      ae1 = 0.01675104 - 0.0000418*t - .000000126*t**2
      asp = .46022931
      do 30 noe = 1,30
   30 cxx(noe) = 0.0
      if(dayb.gt.0.0) dayb = dayb - 1.0
      if(daym.gt.0.0) daym = daym - 1.0
      doby = 0.0
      amit = 0.0
      ami = xyer - xcen
      cplex = xcen/400.0 + 0.0001
      dicf = cplex - aint(cplex)
      if(ami.eq.0.0) go to 32
      xcet = xcen + 1.0
      cdif = xyer - xcet
      doby = cdif/4.0 + 0.0001
      amit = aint(doby)
      if(dicf.lt.0.001) amit = amit + 1.0
   32 farm = 0.25*ami
      farx = farm - aint(farm)
      cxx(1) = xsx + 129.384820*ami + 13.1763968*(dayb + amit) + 0.54901
     16532*grbs
      cxx(2) = xpx + 40.6624658*ami + 0.111404016*(dayb + amit) + 0.0046
     141834*grbs
      cxx(3) = xhx - 0.238724988*ami + 0.985647329*(dayb + amit) + 0.041
     1068639*grbs
      cxx(4) = xp1x + 0.01717836*ami + 0.000047064*(dayb + amit) + 0.000
     1001961*grbs
      cxx(5) = xpx + 40.6624658*ami + 0.111404016*(daym + amit) + 0.0046
     141834*grms
      cxx(6) = xnx - 19.3281858*ami - 0.052953934*(daym + amit) - 0.0022
     106414*grms
   40 cxx(7) = xpx + 40.6624658*ami + 0.111404016*(daym + amit) + 0.0046
     141834*grbs
      cxx(8) = xnx - 19.328185764*ami - 0.0529539336*(dayb + amit) - 0.0
     106414*grbs
      call twopi(cxx, 8)
   41 do 100 ii = 1,8
  100 cxx(ii) = float(int(cxx(ii)*100.0 + 0.5))*0.01
      ang = cxx(8)
      call table6(vib,vb,xib,vpb,vppb,xx,xx,xx,xx,xx,ang,anb,atb)
      cxx(26) = vib
      cxx(27) = vb
      cxx(28) = xib
      cxx(29) = vpb
      cxx(30) = vppb
      cxsb = cxx(1)
      cxpb = cxx(2)
      cxhb = cxx(3)
      cxp1b= cxx(4)
      ang = cxx(6)
      call table6(vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at)
      cxx(9 ) = vi
      cxx(10) = v
      cxx(11) = xi
      cxx(12) = vp
      cxx(13) = vpp
  230 do 333 ii = 9,13
  333 cxx(ii) = float(int(cxx(ii)*100.0 + 0.5))*0.01
      pgx = cxx(5) - cxx(11)
      pgx = float(int(pgx*100.0 + 0.5))*0.01
      call twopi(pgx, 1)
      xpg = pgx*0.0174533
      cxx(14) = pgx
      raxe = sin(2.0*xpg)
      raxn = (cos(0.5*at)**2/(6.0*sin(0.5*at)**2)) - cos(2.0*xpg)
      rxx = 0.0
      if(raxe.eq.0.0.or.raxn.eq.0.0) go to 232
      rax = raxe/raxn
      if(rax.gt.3450.0) go to 232
        rxx   = atan(rax )*pinv
      cxx(22) = rxx
  232 cra = sqrt(1.0 - 12.0*(sin(0.5*at)**2/cos(0.5*at)**2)*cos(2.0*xpg)
     1 + 36.0*(sin(0.5*at)**4/cos(0.5*at)**4))
      um2 = 2.0*(cxx(11) - cxx(10))
      cxx(21) = um2
      cxx(24) = cra
      ul2 = um2 - rxx
  404 ul2 = ul2 + 180.0
  405 cxx(15) = ul2
      zes = (5.0*cos(aw) - 1.0)*sin(xpg)
      zec = (7.0*cos(aw) + 1.0)*cos(xpg)
      call fitan(zes,zec,qxx,spxx,2)
      cxx(23) = qxx
      crav = 0.5*um2 + qxx + 090.0
      cxx(16) = crav
      cqa = sqrt(0.25 + 1.5*((cos(aw)/cos(0.5*aw)**2)*cos(2.0*xpg)) + 2.
     125*(cos(aw)**2/cos(0.5*aw)**4))
      cxx(25) = cqa
      do 444 iii = 14,23
  444 cxx(iii) = float(int(cxx(iii)*100.0 + 0.5 ))*0.01
      pm   = cxx(1)
      pl   = cxx(2)
      sl   = cxx(3)
      ps   = cxx(4)
      plm  = cxx(5)
      skyn = cxx(6)
      vi   = cxx(9)
      v    = cxx(10)
      xi   = cxx(11)
      vp   = cxx(12)
      vpp  = cxx(13)
      p    = cxx(14)
      aul  = cxx(15)
      aum  = cxx(16)
      cra = cxx(24)
      cqa = cxx(25)
      u = v*0.0174533
      q = p*0.0174533
      ui = vi*0.0174533
      return
      end
************************************************************************
      subroutine vanduf (speed,e,f,itype)
      implicit real*8(a-h,o-z)
*     Order of constituents is same as in NAMES common
      dimension spd(180),ms(180)
      common/locat/tm,gonl
      common/fad/ipick
      common/vee/tml,con,u,q,ui
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/boxs/aw,ai,ae,ae1,asp
      common /speeds/spd
      common /mmss/ms

      double precision speed

  500 format('0*** (v + u) not computed for constituent of speed',
     1 f12.7,'  ****',' Execution terminated')
      con = sl + tml
      do 600 j = 1,175
      ipick = j
      if(speed.eq.spd(j)) go to 611
  600 continue
      ipick = 176
  611 if(ipick.gt.164) go to 620
      if(ipick.gt.124) go to 618
      if(ipick.gt.88) go to 616
      if(ipick.gt.38) go to 614
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
     1,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38),ipick
  614 ipik = ipick - 38
      go to (39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58
     1,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80
     2,81,82,83,84,85,86,87,88),ipik
  616 ipikk = ipick - 88
      go to (89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,10
     16,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,
     2123,124),ipikk
  618 ipikkk = ipick - 124
      go to (125,126,127,128,129,130,131,132,133,134,135,136,137,138,139
     1,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,
     2156,157,158,159,160,161,162,163,164),ipikkk
  620 ipiikk = ipick - 164
      go to (165,166,167,168,169,170,171,172,173,174,175,488),ipiikk

    1 e = 2.0*(con - pm + xi - v)
  201 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))
      go to 888

    2 e =  2.0*(con + xi - v) - 3.0*pm + pl
      go to 201

    3 e = 2.0*tml
  203 f = 1.0

      go to 888
    4 e = con - v - 2.0*(pm - xi) - 90.0
      go to 218

    5 e = con - vp + 90.0
  205 f = 1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0
     1.1006)
      go to 888

    6 e = 2.0*con - vpp
  206 f = 1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) + 0
     1.0981)
      go to 888

    7 e = 2.0*con - pm - pl + aul
  207 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))*(1.0/cra)
      go to 888

    8 e = 2.0*(con + xi - v + pl) - 4.0*pm
      go to 201

    9 e = sl - ps + 180.0 + 2.0*tml
      go to 203

   10 e = 2.0*tml - (sl - ps)
      go to 203

   11 e = 2.0*(con + xi - v - sl) - pm + pl + 180.0
      go to 201

   12 e = 2.0*(con + xi - v + sl) - 4.0*pm
      go to 201

   13 e = 2.0*(con + xi - v + sl) - 3.0*pm - pl
      go to 201

   14 e = con + pm - pl - v + 90.0
  214 f = (sin(2.0*aw)*(1.0 - 1.5*sin(ai)**2))/sin(2.0*ui)
      go to 888

   15 e = con - pm + aum
      f =  ((sin(aw)*cos(0.5*aw)**2*cos(0.5*ai)**4)/(sin(ui)*cos(0.5*ui)
     1**2))*(1.0/cqa)
      go to 888
   16 e = con - v + 2.0*(pm - xi)+ 90.0
  216 f = (sin(aw)*sin(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*sin(.5*ui)**2)
      go to 888
   17 e = tml + 270.0 - sl
      go to 203
   18 e = con - v - 3.0*pm + 2.0*xi + pl - 90.0
  218 f = (sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
      go to 888
   19 e = con - v - 4.0*pm + 2.0*xi + 2.0*pl - 90.0
      go to 218
   20 e = con - v - 3.0*pm + 2.0*xi - pl + 2.0*sl - 90.0
      go to 218
   21 e = 4.0*(con - pm + xi - v)
  221 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**2
      go to 888
   22 e = 6.0*(con - pm + xi - v)
  222 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**3
      go to 888
   23 e = 8.0*(con - pm + xi - v)
  223 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**4
      go to 888
   24 e = 4.0*tml
      go to 203
   25 e = 6.0*tml
      go to 203
   26 e = 3.0*(con - pm + xi - v) + 180.0
  226 f = (cos(0.5*aw )**6*cos(0.5*ai)**6)/cos(0.5*ui)**6
      go to 888
   27 e = tml + 180.0
      go to 203
   28 e = 2.0*(con - pm + xi - v) + ( con - vp + 90.0)
  228 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**1*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   29 e = 4.0*(con - pm + xi - v) - (con - vp + 90.0)
  229 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**2*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   30 e = 4.0*(con + xi - v) + pl - 5.0*pm
      go to 221
   31 e = 2.0*(con - pm + xi - v) + 2.0*tml
      go to 201
   32 e = 4.0*tml - 2.0*(con - pm + xi - v)
      go to 201
   33 e = 2.0*(pm - xi)
  233 f = (sin(aw)**2*cos(0.5*ai)**4)/sin(ui)**2
      go to 888
   34 e = 2.0*tml - 2.0*(con - pm + xi - v)
      go to 201
   35 e = pm - pl
  235 f = ((2./3.- sin(aw)**2)*(1.- 1.5*sin(ai)**2))/(2./3.- sin(ui)**2)
      go to 888
   36 e = sl
      go to 203
   37 e = 2.0*sl
      go to 203
   38 e = 2.0*(con - v - 2.0*(pm - xi) - 90.0) - (tml+270.0 - sl)
  238 f = ((sin(aw)*cos(0.5*aw)**2*cos(0.5*ai)**4)/(sin(ui)*cos(0.5*ui)*
     1*2))**2
      go to 888
   39 e = 2.0*(con - pm + xi - v) - (tml + 270.0 - sl)
      go to 201
   40 e = con + 2.0*sl - v - pm - pl + 90.0
      go to 258
   41 e = 2.0*(tml + 270.0 - sl) - (con - v - 2.0*(pm - xi) - 90.0)
      go to 218
   42 e = 2.0*tml - (con - v - 2.0*(pm - xi) - 90.0)
      go to 218
   43 e = 2.0*tml + pm - pl
      go to 221
   44 e = 4.0*(con + xi - v) - 5.0*pm - 2.0*tml  + pl
      go to 221
   45 e = con - v - 2.0*(pm - xi) + tml - sl + 180.0
      go to 218
   46 e = 4.0*con - 2.0*(pm - xi + v + tml) - vpp
      go to 259
   47 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*(pl - tml)
      go to 221
   48 e = 6.0*(con - pm) + 4.0*(xi - v - tml) + aul
      go to 261
   49 e = 6.0*con - 5.0*pm + 4.0*(xi - v - tml) - pl + aul
      go to 261
   50 e = 2.0*(tml + pm - xi + v) - vpp
      go to 259
   51 e = 4.0*(xi - pm - v) + 2.0*(tml + vpp)
  251 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981)**2)
      go to 888
   52 e = 6.0*con - 4.0*tml - 3.0*pm + 2.0*(xi - v) - pl - vpp + aul
  252 cra = sqrt(1.0 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 3
     16.0*(sin(.5*ui)**4/cos(.5*ui)**4))
      f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*(1./cra)
      go to 888
   53 e = 6.0*con - 4.0*tml + 2.0*(xi - v - pm - vpp)
  253 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981)**2)
      go to 888
   54 e = 4.0*tml - 2.0*con + pl - pm + vpp
  254 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   55 e = 4.0*con - 2.0*(tml + vpp) + pm - pl
      go to 251
   56 e = 2.0*tml + con - v - 2.0*(pm - xi) - 90.0
      go to 218
   57 e = 2.0*tml + con - vp + 90.0
      go to 205
   58 e = 3.0*(con - v) + 4.0*xi - 5.0*pm + pl - 90.0
  258 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   59 e = 4.0*con - 2.0*(pm - xi + v) - vpp
  259 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   60 e = 2.0*(tml + con + xi - v) - 3.0*pm + pl
      go to 201
   61 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml - pl + aul
  261 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   62 e = 6.0*(con - pm + xi - v) - 2.0*tml
      go to 222
   63 e = 4.0*con - 3.0*pm + 2.0*(xi - v) - pl + aul
  263 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   64 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*pl
      go to 221
   65 e = 2.0*(con + tml) - pm - pl + aul
      go to 207
   66 e = 5.0*(con - v) - 7.0*pm + 6.0*xi + pl - 90.0
  266 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   67 e = 5.0*(con - v) + 6.0*(xi - pm) - 90.0
      go to 266
   68 e = 5.0*con + 4.0*(xi - pm - v) - vp + 90.0
      go to 229
   69 e = 3.0*con + 2.0*(xi - pm - v + tml) - vp + 90.0
      go to 228
   70 e = 5.0*con + 2.0*(xi - pm - v) - (vp + vpp) + 90.0
  270 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*(1./sqrt(
     2.8965*sin(2.*ui)**2 + .6001*sin(2.*ui)*cos(u) + 0.1006))
      go to 888
   71 e = 4.0*(con - pm + xi - v) + tml - sl + 270.0
      go to 221
   72 e = 6.0*(con - pm + xi - v) - tml + sl - 270.0
      go to 222
   73 e = 5.0*(con - pm) + 4.0*(xi - v) + pl - vp + 90.0
      go to 229
   74 e = 2.0*(con - pm + xi - v) + 4.0*tml
      go to 201
   75 e = 6.0*(con + xi - v) - 7.0*pm + pl
      go to 222
   76 e = 4.0*(con + xi - v) - 5.0*pm + 2.0*tml + pl
      go to 221
   77 e = 4.0*(con - pm + xi - v) + 2.0*tml
      go to 221
   78 e = 8.0*con - 9.0*pm + 6.0*(xi - v) - 2.0*tml + pl + aul
  278 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**4*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   79 e = 6.0*(con + xi - v) - 8.0*pm + 2.0*pl
      go to 222
   80 e = 4.0*con - 3.0*pm + 2.0*(xi - v + tml) - pl + aul
      go to 263
   81 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - pl + aul
      go to 261
   82 e = 4.0*con + 2.0*(xi - pm - v + tml) - vpp
      go to 259
   83 e = 8.0*(con - pm) + 6.0*(xi - v) - 2.0*tml + aul
      go to 278
   84 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - 2.0*tml - pl + aul
      go to 278
   85 e = 6.0*con + 4.0*(xi - pm - v) - vpp
      go to 254
   86 e = 7.0*(con - v) - 9.0*pm + 8.0*xi + pl - 90.0
  286 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   87 e = 7.0*con + 6.0*(xi - v) - 8.0*pm + 2.0*pl - vp + 90.0
  287 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**3*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   88 e = 5.0*(con - v) + 6.0*(xi - pm) + 2.0*tml - 90.0
      go to 266
   89 e = 5.0*con + 4.0*(xi - pm) - 3.0*v + 2.0*tml - vpp - 90.0
  289 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*( (sin(aw
     2)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   90 e = 6.0*con - 7.0*pm + 6.0*(xi - v) + 2.0*tml + pl
      go to 222
   91 e = 6.0*(con - pm + xi - v) + 2.0*tml
      go to 222
   92 e = 4.0*(con - pm + xi - v) + 4.0*tml
      go to 221
   93 e = 8.0*(con + xi - v) - 10.0*pm + 2.0*pl
      go to 223
   94 e = 8.0*(con + xi - v) - 9.0*pm + pl
      go to 223
   95 e = 6.0*con - 5.0*pm + 4.0*(xi - v) + 2.0*tml - pl + aul
      go to 261
   96 e = 10.0*con - 9.0*pm + 8.0*(xi - v) - 2.0*tml - pl + aul
  296 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**5*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   97 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - pl + aul
      go to 278
   98 e = 8.0*con + 6.0*(xi - pm - v) - vpp
  298 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   99 e = 6.0*con + 4.0*(xi - pm - v) + 2.0*tml - vpp
      go to 254
  100 e = 9.0*con - 10.0*pm + 8.0*(xi - v) + 2.0*pl - vp + 90.0
  300 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**4*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
  101 e = 9.0*(con - pm) + 8.0*(xi - v) + pl - vp + 90.0
      go to 300
  102 e = 9.0*con + 8.0*(xi - pm - v) - vp + 90.0
      go to 300
  103 e = 7.0*con + 6.0*(xi - pm - v) + 2.0*tml - vp + 90.0
      go to 287
  104 e = 10.0*(con + xi - v) - 11.0*pm + pl
  304 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**5
      go to 888
  105 e = 10.0*(con - pm + xi - v)
      go to 304
  106 e = 8.0*(con + xi - v) - 9.0*pm + pl + 2.0*tml
      go to 223
  107 e = 8.0*(con - pm + xi - v) + 2.0*tml
      go to 223
  108 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - pl + aul
      go to 278
  109 e = 6.0*(con - pm + xi - v) + 4.0*tml
      go to 222
  110 e = 9.0*con + 8.0*(xi - pm - v) + 2.0*tml - vp + 90.0
      go to 300
  111 e = 10.0*(con + xi - v) - 11.0*pm + pl
      go to 304
  112 e = 10.0*(con - pm + xi - v) + 2.0*tml
      go to 304
  113 e = 10.0*con - 9.0*pm + 8.0*(xi - v) + 2.0*tml - pl + aul
      go to 296
  114 e = 8.0*(con - pm + xi - v) + 4.0*tml
      go to 223
  115 e = tml - 2.0*sl + ps - 90.0
      go to 203
  116 e = con + sl - ps + 90.0
      go to 203
  117 e = con + 2.0*sl + 90.0
      go to 203
  118 e = tml + pm - sl + pl - v + 90.0
      go to 258
  119 e = 2.0*con + pm - v - vp - pl + 180.0
  319 f =  ((sin(2.*aw)*(1.0 - 1.5*sin(ai)**2))/sin(2.0*ui))*(1./sqrt(.8
     1965*sin(2.*ui)**2 + .6011*sin(2.*ui)*cos(u) + 0.1006))
      go to 888
  120 e = 2.0*(con - v) - 5.0*pm + 4.0*xi + pl + 180.0
      go to 238
  121 e = 3.0*(con - v) - 4.0*(pm - xi) - 90.0
      go to 258
  122 e = 2.0*(con + tml) - vpp
      go to 206
  123 e = con + 2.0*(pm - xi - vp) + v + 270.0
  323 f = ((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
     1)**1*(1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     2 0.1006))**2
      go to 888
  124 e = con - 2.0*v - 4.0*(pm - xi) + vp - 270.0
  324 f = ((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
     1)**2*(1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     2 0.1006))
      go to 888
  125 e = 6.0*(con - pm) + 4.0*(xi - v - tml) + 2.0*pl - vpp
  325 go to 254
  126 e = 6.0*con - 5.0*pm + 4.0*(xi - v - tml) + pl - vpp
  326 go to 254
  127 e = 6.0*con - 4.0*tml - 3.0*pm + 2.0*(xi - v - vpp) + pl
  327 go to 253
  128 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml + pl - vpp
  328 go to 254
  129 e = 4.0*con - 3.0*pm + 2.0*(xi - v) + pl - vpp
  329 go to 259
  130 e = 8.0*con - 9.0*pm + 6.0*(xi - v) + 3.0*pl - 2.0*tml - vpp
  330 go to 298
  131 e = 8.0*(con - pm) + 6.0*(xi - v) + 2.0*(pl - tml) - vpp
  331 go to 298
  132 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - 2.0*tml + pl - vpp
  332 go to 298
  133 e = 6.0*con - 5.0*pm + 4.0*(xi - v) + pl - vpp
  333 go to 254
  134 e = 4.0*con - 3.0*pm + 2.0*(xi - v + tml) + pl - vpp
  334 go to 259
  135 e = 10.0*con - 9.0*pm + 8.0*(xi - v) - 2.0*tml + pl - vpp
  335 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**4*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
  136 e = 8.0*con - 7.0*pm + 6.0*(xi - v) + pl - vpp
  336 go to 298
  137 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml + pl - vpp
  337 go to 254
  138 e = 8.0*con - 7.0*pm + 6.0*(xi - v) + 2.0*tml + pl - vpp
  338 go to 298
  139 e = 10.0*con - 9.0*pm + 8.0*(xi - v) + 2.0*tml + pl - vpp
  339 go to 335
  140 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*pl + sl - tml - 270.0
  340 go to 221
  141 e = 3.0*sl - 2.0*vp - tml - 90.0
  341 f = 1.0/(.8965*sin(2.0*ui)**2 + .6001*sin(2.0*ui)*cos(u) + .1006)
      go to 888
  142 e = tml + vp - 3.0*sl + 90.0
  342 go to 205
  143 e = 2.0*tml - vp
  343 go to 205
  144 e = 2.0*(tml - sl) + vpp
  344 go to 206
  145 e = 2.0*(tml - vpp) + 4.0*sl
  345 f = (1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981))**2
      go to 888
  146 e = 2.0*tml - 2.0*(sl - ps)
  346 go to 203
  147 e = 4.0*tml + ps - sl
  347 go to 203
  148 e = 4.0*tml - 2.0*sl + vpp
  348 go to 206
  149 e = 4.0*tml + 6.0*sl - 3.0*vpp
  349 f = (1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981))**3
      go to 888
  150 e = 4.0*tml - 3.0*(sl - ps)
  350 go to 203
  151 e = con - v - 2.0*(pm - xi) + tml + 90.0
  351 go to 218
  152 e = 2.0*(pm - xi) + v + tml - con - 90.0
  352 go to 218
  153 e = pm - xi - 90.0
  353 f = 0.31920/(sin(ui) - sin(ui)**3)
      go to 888
  154 e = pm - sl
  354 go to 235
  155 e = sl - ps
  355 go to 203
  156 e = 3.0*sl
  356 f = 1.0
      go to 888
  157 e = 4.0*sl
  357 go to 356
  158 e = 6.0*sl
  358 go to 356
  159 e = 8.0*sl
  359 go to 356
  160 e = 10.0*sl
  360 go to 356
  161 e = 12.0*sl
  361 go to 356
  162 e = 24.0*sl
  362 go to 356
  163 e = 3.0*tml + 180.0
  363 go to 356
  164 e = 5.0*tml + 180.0
  364 go to 356
  165 e = 2.0*(con - v) - 4.0*(pm - xi) - 180.0
  365 f =((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
     1**2
      go to 888
  166 e = tml + con - vp + 270.0
  366 go to 205
  167 e = 3.0*(con - pm) + 2.0*(xi - v) + pl - vp + 90.0
  367 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/cos(0.5*ui)**4 )*(1.0/sqrt(0.
     18965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
  168 e = 3.0*tml + 270.0 - sl
  368 go to 203
  169 e = 3.0*con - vpp - vp + 90.0
  369 f = 1.0/(sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981)*sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     20.1006))
      go to 888
  170 e = 4.0*(con - v) + 6.0*xi - 7.0*pm + pl - 180.0
  370 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))*((sin(aw)*c
     1os(0.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))**2
      go to 888
  171 e = 4.0*(con - v) - 6.0*(pm - xi) - 180.0
  371 go to 370
  172 e = 2.0*(tml + con - v) - 4.0*(pm - xi) - 180.0
  372 go to 365
  173 e = 7.0*tml + 180.0
  373 go to 356
  174 e = 8.0*tml
  374 go to 356
  175 e = 10.0*tml
  375 go to 356
  488 print 500, speed
!sfan      stop
  888 go to (400,390),itype
  390 e = e + float(ms(ipick))*gonl - spd(ipick)*(tm/15.0)
  400 call twopi(e,1)
      f = 1.0/f
      return
      end
************************************************************************
      subroutine name (spdd,itag,isub,inum,icode)
!     this subroutine identifies the constituent by its speed,
*     name label or constituent number,
!     and makes it available for labeling.
*     ICODE = 1, by speed
*     ICODE = 2, by label
*     ICODE = 3, by number
!     it also determine the subscript of the constituent
!        order of constituent speeds***  m(2),n(2),s(2),o(1),k(1),k(2)
!      l(2),2n(2)r(2),t(2),lambda(2),mu(2),nu(2),j(1),m(1),oo(1),p(1)
!      q(1),2q(1),rho(1),m(4),m(6),m(8),s(4),s(6),m(3),s(1),mk(3),2mk(3)
!      mn(4),ms(4),2sm(2),mf,msf,mm,sa,ssa
      real*8 spd,spdd
      common /mmss/ip
      common /speeds/spd
      common /names/lable
      character*10 lable(180)*10,itag
      dimension spd(180),ip(180)

    1 format(10x,'Constituent of speed ',f12.7,' not in list.')
    2 format(10x,'Constituent ', a10,' not in the list.')
    3 format(10x,'Constituent no.',i4,' not in list.')

      if (icode.eq.1) then
*       search by speed
        do 100 j = 1,175
        if(spdd.ne.spd(j)) go to 100
        itag = lable(j)
        isub = ip(j)
        inum = j
        return
100     continue
        print 1,spdd
      else if (icode.eq.2) then
*       search by name
        do 200 i = 1,175
        if(itag.ne.lable(i)) go to 200
        spdd = spd(i)
        isub = ip(i)
        inum = i
        return
200     continue
        print 2, itag
      else if (icode.eq.3) then
*       search by number
!        if (i.gt.0.and.i.le.175) then
!          itag = lable(k)
!          spdd = spd(k)
!          isub = ip(k)
!          return
!        end if
        write(*,*)'get wrong icode, icode must be 1 or 2'
        stop
        print 3, inum
      end if

      stop '**** Execution terminated in NAME (illegal icode) ****'

      end
************************************************************************
      subroutine orbit (xcen,xsx,xpx,xhx,xp1x,xnx,oex,t,xyer,nnn)
      implicit real*8(a-h,o-z)
      real*8 oex(nnn)
      s = 13.1763968
      p = 0.1114040
      xh = 0.9856473
      p1 = 0.0000471
      xn = -.0529539
      xcan = xyer*0.01 + 0.001
      xcen = aint(xcan)*100.0
      t = -3.0
      yr = 2.5
      gat = 1600.0
      do 10 jk = 1,30
      gp = gat/400.0 + 0.00001
      col = gp - aint(gp)
      if(col.lt.0.010) go to 11
      if(gat.eq.xcen) go to 12
      yr = yr - 1.0
      go to 9
   11 if( gat.eq.xcen) go to 12
    9 gat = gat + 100.0
   10 continue
   12 t = (gat - 1900.0)*0.01
      oex(1) = 270.437422 + 307.892*t + 0.002525*t**2 + .00000189*t**3 +
     1 yr*s
      oex(2) = 334.328019 + 109.032206*t - 0.01034444*t**2 - .0000125*t*
     1*3 + yr*p
      oex(3) = 279.696678 + 0.768925*t + .0003205*t**2 + yr*xh
      oex(4) = 281.220833 + 1.719175*t + 0.0004528*t**2 + .00000333*t**3
     1 + yr*p1
      oex(5) = 259.182533 - 134.142397*t + .00210556*t**2 + .00000222*t*
     1*3 + yr*xn
      call twopi(oex,5)
      do 100 i = 1,5
  100 oex(i) = float(int(oex(i)*100.0 + 0.5))*0.01
      xsx = oex(1)
      xpx = oex(2)
      xhx = oex(3)
      xp1x = oex(4)
      xnx = oex(5)
      return
      end
************************************************************************
      subroutine table6(vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at)
      implicit real*8(a-h,o-z)      
      common/boxs/aw,ai,ae,ae1,asp
      v = 0.0
      xi = 0.0
      vp = 0.0
      vpp = 0.0
      an = ang*0.0174533
      ax = ang
      eye = cos(ai)*cos(aw) - sin(ai)*sin(aw)*cos(an)
      c9 = acos(eye)*57.2957795
      vi = float(int(c9*100.0 + 0.5))*0.01
      cig = vi*0.0174533
      at = cig
      if(cig.eq.0.0) go to 230
      if(ax.eq.0.0.or.ax.eq.180.0) go to 230
      vxxe = sin(ai)*sin(an)
      vxxn = cos(ai)*sin(aw) + sin(ai)*cos(aw)*cos(an)
      if(vxxe.eq.0.0.or.vxxn.eq.0.0) go to 201
      vxx = vxxe/vxxn
      c10 = atan(vxx)*57.2957795
      v = float(int(c10*100.0 + 0.5))*0.01
      if(ax.gt.180.0.and.v.gt.0.0) v = -1.0*v
  201 cvx = v*0.0174533
      term = sin(ai)*(cos(aw)/sin(aw))
      exx = term*(sin(an)/cos(an)) + (cos(ai) - 1.0)*sin(an)
      if(exx.eq.0.0) go to 202
      ezz = term + cos(ai)*cos(an) + (sin(an)**2/cos(an))
      if(ezz.eq.0.0) go to 202
      exez = exx/ezz
      if(exez.gt.3450.0) go to 202
      c11 = atan(exez)*57.2957795
      xi = float(int(c11*100.0 + 0.5))*0.01
      if(ax.gt.180.0.and.xi.gt.0.0) xi = -1.0*xi
  202 cex = xi*0.0174533
      a22 = (0.5 + 0.75*ae**2)*sin(2.0*cig)
      b22 = (0.5 + 0.75*ae1**2)*sin(2.0*aw)*asp
      vpxe = a22*sin(cvx)
      vpxn = a22*cos(cvx) + b22
      if(vpxe.eq.0.0.or.vpxn.eq.0.0) go to 203
      vpx = vpxe/vpxn
      if(vpx.gt.3450.0) go to 203
      vp = atan(vpx)*57.2957795
      if(ax.gt.180.0.and.vp.gt.0.0) vp = -1.0*vp
  203 pvc = vp*0.0174533
      a47 = (0.5 + 0.75*ae**2)*sin(cig)**2
      b47 = (0.5 + 0.75*ae1**2)*asp*sin(aw)**2
      vpye = a47*sin(2.0*cvx)
      vpyn = a47*cos(2.0*cvx) + b47
      if(vpye.eq.0.0.or.vpyn.eq.0.0) go to 204
      vpy = vpye/vpyn
      if(vpy.gt.3450.0) go to 204
      vpp = atan(vpy)*57.2957795
      if(ax.gt.180.0.and.vpp.gt.0.0) vpp = -1.0*vpp
  204 pvcp = vpp*0.0174533
  230 return
      end
!-------------------------------------------------------------------------------
      SUBROUTINE ncrght(line,nc)
! -
! - Returns the last non blank character position in a string
! -
      CHARACTER*(*) line
      ilim = LEN (line)
      DO 100 i = 1,ilim
        IF(line(i:i) .ne. ' ') THEN
          nc = i
        END IF
  100 CONTINUE
      RETURN
      END
      FUNCTION JULIAN(yr,month,day,hour)
      real*8 jday,yr,month,day,hour
      real*8 Y, m ,julian
!      integer yr

      Y = yr
      if (yr.lt.100. .and. yr.gt.50.) then
        Y = yr+1900.
      elseif (yr.lt.100. .and. yr.le.50.) then
        Y = yr+2000.
      endif

      if (month.le. 2.) then
        Y = Y-1
        m = month +12.
      else
        Y = Y
        m = month
      endif

      JULIAN = aint(365.25*Y) + aint(30.6001*(m+1))
     &        + day + hour/24. + 1720981.50

!  1980 epic....
!      JULIAN = aint(365.25*(Y-1980.)) + aint(30.6001*(m+1))
!     &        + day + hour/24.
!      write(*,'(''Julian ='',f12.4,4f8.2)') julian, yr,month,day,hour
!      END FUNCTION JULIAN
       return
       END

      SUBROUTINE GREGORIAN(jday,yr,month,day,hour)
      real*8 jday,yr,month,day,hour,dayoweek,week

      a = aint(jday+.5)
      b = a+1537
      c = aint((b-122.1)/365.25)
      d = aint(365.25*c)
      e = aint((b-d)/30.6001)
      day = b-d-aint(30.6001*e)+mod(jday+0.5d00,1.0d00)
      month = e-1-12*aint(e/14)
      yr = c-4715-aint((7+month)/10)
      hour = mod(jday+0.5d00,1.0d00)*24.
      day = aint(day)

!      dayoweek = mod(aint(jday+0.5),7)
!      week = aint((jday-2444244.5)/7)
!  2000 leap year crashed.
      if(month.eq.2. .and. day .eq. 31. ) then
       yr=yr-1.
       day = 29.
      endif

      END
      SUBROUTINE CDAY(IDD,IMM,IYY,KD,KULL)  !hli (MOD to KULL)
!
! change 2-digit year to 4-digit number (hli, 04/07/2000)
!******************* A DAY,MONTH AND YEAR(2 DIGITS) ARE GIVEN IN
!******************* ARGUMENTS 1-3.  THE DAY COUNT KD IS RETURNED.
!******************* THE ORIGIN IS SUCH THAT KD=367 ON 1/1/1901.
!******************* THIS ROUTINE SHOULD NOT BE USED OUTSIDE THE
!******************* RANGE 1901 TO 1999.
!
      IMPLICIT INTEGER*4 (I-N)
      INTEGER*2 KR,LU,LP
      DIMENSION NDM(12),NDP(12)
      DATA NDP/0,31,59,90,120,151,181,212,243,273,304,334/
      DATA NDM/31,28,31,30,31,30,31,31,30,31,30,31/
      GO TO (100,200),KULL
100   IQUOT= IYY/4
      IREM = IYY-IQUOT*4
      L=1
      IF(IREM)3,1,3
1     IF(IMM-2)2,2,3
2     L=0
3     KD=IQUOT*1461+L+IREM*365+NDP(IMM)+IDD
      RETURN
!
!
200   CONTINUE                     !*      ENTRY DMY
!
!******************* GIVEN THE DAY COUNT KD , THE DAY,MONTH AND YEAR
!******************* ARE RETURNED .
!
      NDM(2)=28
      IQQ = (KD-1)/1461
      IYY = IQQ*4
      IDD = KD - IQQ*1461 -366
      IF(IDD)20,20,30
20    NDM(2)=29
      IDD =IDD +366
      GO TO 51
30    IYY = IYY +1
      IDD = IDD -365
      IF(IDD)50,50,30
50    IDD =IDD + 365
51    IMM = 0
60    IMM = IMM +1
      IDD = IDD -NDM(IMM)
      IF(IDD)70,70,60
70    IDD = IDD +NDM(IMM)
      RETURN
      END
