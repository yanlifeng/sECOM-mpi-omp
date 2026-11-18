      PROGRAM ECOM_3D
!
!*********************************************************
!     ECOMSED MODEL
!     VERSION 1.1
!     OCTOBER 1996
!*********************************************************
!
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
      INTEGER ROWIMAX(JM),ROWIMIN(JM)
! 
!----------------------------------------------------------------------|
!                  GENERAL CIRCULATION MODEL                           |
!          ORTHOGONAL CURVILINEAR COORDINATE VERSION                   |
!                         ECOM3D                                       |
!                                                                      |
!                                                                      |
!     THIS IS A VERSION OF THE THREE DIMENSIONAL, TIME DEPENDENT,      |
!   PRIMITIVE EQUATION, CIRCULATION MODEL DEVELOPED BY GEORGE MELLOR   |
!   AND ALAN BLUMBERG WITH SUBSEQUENT CONTRIBUTIONS BY LEO OEY AND     |
!   BORIS GALPERIN.                                                    |
!                                                                      |
!                                                                      |
!     FOR DETAILS OF THE GOVERNING EQUATIONS AND SOLUTION              |
!     TECHNIQUES THE INTERESTED READER IS REFERRED TO:                 |
!                                                                      |
!     BLUMBERG, A.F. AND G.L. MELLOR, "DIAGNOSTIC AND PROGNOSTIC       |
!     NUMERICAL CIRCULATION STUDIES IN THE SOUTH ATLANTIC BIGHT,"      |
!     J. GEOPHYS. RES., 88, 4579-4592, 1983.                           |
!                                                                      |
!     AND                                                              |
!                                                                      |
!     BLUMBERG, A.F. AND G.L. MELLOR, "A DESCRIPTION OF A THREE-       |
!     DIMENSIONAL COASTAL OCEAN CIRCULATION MODEL,"  IN:  THREE        |
!     DIMENSIONAL SHELF MODELS, COASTAL AND ESTUARINE SCIENCES,        |
!     5, N. HEAPS, ED., AMERICAN GEOPHYSICAL UNION, 1987.              |
!                                                                      |
!                                                                      |
!     IN SUBROUTINE PROFQ THE MODEL MAKES USE OF THE TURBULENCE        |
!     CLOSURE SUB-MODEL DESCRIBED IN:                                  |
!                                                                      |
!     MELLOR, G.L. AND T. YAMADA, "DEVELOMENT OF A TURBULENCE          |
!     CLOSURE MODEL FOR GEOPHYSICAL FLUID PROBLEMS," REV. GEOPHYS.     |
!     SPACE PHYS., 20, 851-875, 1982.                                  |
!                                                                      |
!     AND                                                              |
!                                                                      |
!     GALPERIN, B., L.H. KANTHA, S. HASSID AND A. ROSATI, "A QUASI-    |
!     EQUILIBRIUM TURBULENT ENERGY MODEL FOR GEOPHYSICAL FLOWS,"       |
!     J. ATMOS. SCI. 45, 55-62, 1988.                                  |
!                                                                      |
!        THIS PROGRAM CONTAINS A NUMBER OF COMMENT CARDS THAT          |
!     SHOULD ENABLE AN ASTUTE OCEAN SCIENTIST TO MAKE SIMULATIONS.     |
!     PLEASE DIRECT CRITICISMS AND SUGGESTIONS TO  ME.                 |
!                                                                      |
!                                             ALAN BLUMBERG            |
!                                             HYDROQUAL, INC           |
!                                               MAHWAH, NJ             |
!                                               12/28/89               |
!                                                                      |
!----------------------------------------------------------------------|
!                                                                      |
!                   LOOP LIMITS                                        |
!                                                                      |
!              T,S,etc :  J=2,JMM1                                     |
!                         I=2,IMM1                                     |
!                                                                      |
!                 U    :  J=2,JMM1                                     |
!                         I=3,IMM1                                     |
!                                                                      |
!                 V    :  J=3,JMM1                                     |
!                         I=2,IMM1                                     |
!                                                                      |
!----------------------------------------------------------------------|
!
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION COM(80)
      DIMENSION IVAR(IM),PRT(IM,KB)
!
      REAL DUMB1(QBCM),DUMB2(EBCM,KBM1),DUMB3(DBCM)
      DIMENSION DUMP1(IM,JM),DUMP2(IM,JM,KB)
!     Main-program grids have process lifetime and are shared by workers.
!     Keep them off the initial thread stack even in an OpenMP build.
      SAVE ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV
      SAVE DUMP1,DUMP2,COHTHK,NCOHTHK
!hli
      INTEGER TMP_IYR,TMP_IMO,TMP_IDA,TMP_IHOUR
!hli
!
!  LIMIT DEPOSITION TO BATHYMETRIC DEPTH 
!  PRAVI 08/12/99
!
      REAL COHTHK(IM,JM),NCOHTHK(IM,JM)
!
      EQUIVALENCE (IVAR,TPS),(PRT,A)
!
      CHARACTER*10 RESTAR,OPTION
      CHARACTER*10 OPTBFRIC,OPTHORCON,OPTUMOL,TROPT
      CHARACTER*4 IPLTFORM,OPTAVG,ITRNFORM
      CHARACTER OFORM*3,TRICOPT*5
!qa   use common block when HOT START and EXTERNAL
      COMMON /COAST1/ICNT,INDX(MAXWET),JNDX(MAXWET),RESTAR
*
****** NEW VARIABLES OF antidif.f BC  Khan 080999     
       DIMENSION TBDRY2(EBCM,KBM1), SBDRY2(EBCM,KBM1)
****** NEW VARIABLES OF antidif.f BC  Khan 080999     
       DIMENSION HRSRESTART(99) ! NG 2009, as ammended for multiple restarts in 2011
       CHARACTER*10 STARTFL
*
!NG
      Data UTCSHIFT /0.0/      ! DANGER:... ONLY IF START DATE in run_data
                               ! IS IN LST (E.G. if New York LST, UTCSHIFT=-5.0)
!NG      Data WDSMAX /50.0/       ! Max wind used for wind stress calculations
      Data WDSMAX /100./       ! Max wind used for wind stress calculations
      Data PATMOPT /1./        ! Set PATMOPT to 1 only if barometric load adjustment
                               ! of surface elevation is to be included
! 
      IURUN=1                 ! run_data
      IUGRD=3                 ! model_grid (opened & closed in setdom.f)
      IUTAS=5                 ! init_tands (opened & closed in tands.f)
      IUUAV=7                 ! synop_wind (opened & closed in bcdata.f)
      IURRT=8                 ! hrsrestart (opened & closed in ecom3d.f, optional) ! CNG06152009
!
      IUPRT=10                ! gcmprt
      IUPLT=12                ! gcmplt
      IUTSR=14                ! gcmtsr
      IUWRS=16                ! startup (opened & closed in ecom3d.f)
      IUTRN=18                ! gcm_geom (opened & closed in transport.f)
!                             ! gcm_tran (opened & closed in transport.f)
!
      IUT65=65                ! elevation potential cells NG04232008
!       IUT90=90                ! elevation boundary conditions
!       IUT91=91                ! river discharges
!       IUT92=92                ! diffuser intake/outfall
      IUT93=93                ! meteorological data
!       IUT94=94                ! temperature & salinity boundary conditions
      IUT95=95                ! synoptic wind stress comp (opened in bcdata.f)
!       IUT96=96                ! diffuser intake/outfall in loop
!       IUT103=103              ! synoptic wind and  heat flux input file
      IUT104=104              ! synoptic precipitation and  evaporation flux input file (synop_pet, optional) NG 2014
      IUT191=191              ! synoptic heat flux input file
      IUT192=192              ! synoptic heat flux (temp file)
!NG04092008 Not anymore      IUT193=193              ! synoptic wind velocity components 
!NG04092008c                               (opened in bcdata.f) used by wave model
      IUT194=194              ! 2d varying extinction coefficients
      IUT801=801              ! 2d extinction coefficients (user input)
!
      IUT195=195              ! wave parameters at open b.c.       
      IUT802=802              ! wave parameters at open b.c. (wave.inp, optional) NG 2007+
!
!********************************************************************
!
!  THESE FILES ARE OPENED IN bcdata.f
!
!  DISSOLVED TRACER TRANSPORT
!
      IUT501=501                ! tracer conc. at open b.c.
      IUT601=601                ! tracer conc. at river discharges
      IUT98=98                  ! tracer conc. at diffuser intake/outfall
      IUT99=99                  ! tracer conc. at diffuser intake/outfall loop
      IUT701=701                ! tracer load at point source
!
!  SEDIMENT TRANSPORT
!
      IUT502=502                !  coh sed conc. at open b.c.
      IUT503=503                !  non-coh sed conc. at open b.c.
      IUT602=602                !  coh sed conc. at river discharges
      IUT603=603                !  non-coh sed conc.river discharges
      IUT702=702                !  coh sed conc. at Diffuser discharges
      IUT703=703                !  non-coh sed conc Diffuser discharges
!
!  PARTICLE-BOUND TRACER TRANSPORT
!
      IUT504=504                !  coh sed conc. at open b.c.
      IUT505=505                !  non-coh sed conc. at open b.c.
      IUT604=604                !  coh sed conc. at river discharges
      IUT605=605                !  non-coh sed conc.river discharges
      IUT704=704                !  coh sed conc. at Diffuser discharges
      IUT705=705                !  non-coh sed conc Diffuser discharges
!
!********************************************************************
!
!
      OPEN (IURUN,FILE='run_data.txt') ! A.Jordi
!       OPEN (IURUN,FILE='run_data')
      OPEN (IUPRT,FILE='gcmprt')
!       OPEN (IUT90,FILE='gcm_temp90')
!       OPEN (IUT91,FILE='gcm_temp91')
!       OPEN (IUT92,FILE='gcm_temp92')
      OPEN (IUT93,FILE='gcm_temp93')
!      OPEN (IUT193,FILE='gcm_temp193')
!       OPEN (IUT94,FILE='gcm_temp94')
!       OPEN (IUT96,FILE='gcm_temp96')
!
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,11)  (COM(I),I=1,80)
 11   FORMAT(80A1)
 12   FORMAT(/1X,80A1/)
!
!********************************************************************
!
!  HYDTYPE = 'INTERNAL' :  USE ECOM HYDRODYNAMICS
!          = 'EXTERNAL' :  USE POM HYDRODYNAMICS INPUT FROM hqi_tran
!
!  NHYD = NO. OR TIMESTEPS BETWEEN EACH hqi_tran INPUT (ASSUMED CONSTANT)
!
!  WAVEDYN = 'NEGLECT ' ==>  NO EFFECT OF WAVES ON BOTTOM FRICTION
!          = 'SMBMODEL' ==>  WAVES EFFECT BOT FRIC COEFF - SMB THEORY
!          = 'DONMODEL' ==>  WAVES EFFECT BOT FRIC COEFF - Donelan(1977) THEORY
!          = 'DONONLY ' ==>  WAVES ONLY CALCULATION - Donelan(1977) THEORY
!          = 'EXTERNAL' ==>  WAVES EFFECT BOT FRIC COEFF - INPUT WAM RESULTS
!
!  TRACER = 'NEGLECT' :  NO CONSERVATIVE TRACER CALC.
!         = 'INCLUDE' :  CONSERVATIVE TRACER
!
!  SEDTRAN = 'NEGLECT' :  NO SEDIMENT TRANSPORT
!          = 'INCLUDE' :  SEDIMENT TRANSPORT
!
!  CHEMTRAN = 'NEGLECT' :  NO CHEMICAL TRANSPORT (PARTICLE BOUND)
!           = 'INCLUDE' :  CHEMICAL TRANSPORT (PARTICLE BOUND)
!                           (ONLY USE IF SEDTRAN = 'INCLUDE')
!
!  PARTICLE = 'NEGLECT' :  NO PARTICLE TRACKING
!           = 'INCLUDE' :  PARTICLE TRACKING
!
!  SEDTYPE = 'BOTH'  :  COHESIVE & NON-COHESIVE
!          = 'SAND'  :  NON-COHESIVE ONLY
!          = 'MUD '  :  COHESIVE ONLY
!
!  NOTE:  SEDTYPE APPLIES TO SEDIMENT TRANSPORT AND PARTICLE-BOUND
!         CHEMICAL (TRACER) TRANSPORT
!
      READ(IURUN,1) HYDTYPE,WAVEDYN,TRACER,SEDTRAN,CHEMTRAN,SEDTYPE,
     +                PARTICLE
 1    FORMAT(2X,A8,2X,A8,3X,A7,3X,A7,3X,A7,6X,A4,3X,A7)
!
          IF (WAVEDYN.NE.'SMBMODEL'.AND.WAVEDYN.NE.'DONMODEL'.AND.
     .        WAVEDYN.NE.'EXTERNAL'.AND.WAVEDYN.NE.'NEGLECT '.AND.
     .        WAVEDYN.NE.'DONONLY ')THEN
              WRITE(IUPRT,6112)WAVEDYN
              CALL SYSTEM ('rm gcm_temp*')
              STOP
          ENDIF
!qa
!
!  OPEN AUXILIARY INPUT FILES
!
!  FILE                UNIT      DESCRIPTION
!  water_trace.inp     IUT401    DISSOLVED TRACER
!  coh_sed.inp         IUT402    SEDIMENT TRANSPORT, COHESIVE
!  noncoh_sed.inp      IUT403    SEDIMENT TRANSPORT, NON-COHESIVE
!  coh_trace.inp       IUT404    PARTICLE-BOUND TRACER, COHESIVE
!  noncoh_trace.inp    IUT405    PARTICLE-BOUND TRACER, NON-COHESIVE
!  partrack.inp        IUT406    PARTICLE TRACKING
!
      IF (TRACER.EQ.'INCLUDE') THEN 
        IUT401=401
        OPEN (UNIT=IUT401,FILE='water_trace.inp',FORM='FORMATTED')
      ENDIF
!
      IF (SEDTRAN.EQ.'INCLUDE') THEN 
        IF (SEDTYPE.EQ.'MUD ') THEN
          IUT402=402
          OPEN (UNIT=IUT402,FILE='coh_sed.inp',FORM='FORMATTED')
!
          KSED=1
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND') THEN
          IUT403=403
          OPEN (UNIT=IUT403,FILE='noncoh_sed.inp',FORM='FORMATTED')
!
          KSED=1
        ENDIF
!
        IF (SEDTYPE.EQ.'BOTH') THEN
          IUT402=402
          IUT403=403
          OPEN (UNIT=IUT402,FILE='coh_sed.inp',FORM='FORMATTED')
          OPEN (UNIT=IUT403,FILE='noncoh_sed.inp',FORM='FORMATTED')
!
          KSED=2
        ENDIF
      ENDIF
!
      IF (CHEMTRAN.EQ.'INCLUDE') THEN 
        IF (SEDTYPE.EQ.'MUD ') THEN
          IUT404=404
          OPEN (UNIT=IUT404,FILE='coh_trace.inp',FORM='FORMATTED')
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND') THEN
          IUT405=405
          OPEN (UNIT=IUT405,FILE='noncoh_trace.inp',FORM='FORMATTED')
        ENDIF
!
        IF (SEDTYPE.EQ.'BOTH') THEN
          IUT404=404
          IUT405=405
          OPEN (UNIT=IUT404,FILE='coh_trace.inp',FORM='FORMATTED')
          OPEN (UNIT=IUT405,FILE='noncoh_trace.inp',FORM='FORMATTED')
        ENDIF
      ENDIF
!
      IF (PARTICLE.EQ.'INCLUDE') THEN 
        IUT406=406
        OPEN (UNIT=IUT406,FILE='partrack.inp',FORM='FORMATTED')
      ENDIF
!
!********************************************************************
!
!
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,2) DTI,ISPLIT,IRAMP,IYR,IMO,IDA,IHR,NHYD,SGW,WETEPS,
     &              WETMIN
      !IRAMP=86400*5
      IHOUR=IHR
 2    FORMAT(1F10.4,7I5,F5.0,2F10.0)
!

! keep year being 4-digit number (hli, 04/07/2000)
!      IF(IYR.GE.100)THEN
!         ICE=IYR/100
!         IYR=IYR-ICE*100       ! i.e. 1996 --> 96
!      ELSE
!         ICE=19   ! century
!      END IF
! end
      CALL CDAY(IDA,IMO,IYR,IJDAY,1)
      SDAY=float(IJDAY)+FLOAT(IHOUR)/24.
      CALL CDAY(1,1,IYR,J1JDAY,1)
      J1YR = IYR
      CALL CDAY(1,1,J1YR+1,J2JDAY,1)      ! QA
!end,   modified for heatflux calculation by jeff ji, 8/30/94
!
!
      DTE=DTI/FLOAT(ISPLIT)
      DTE2=2.*DTE
      DTI2=2.*DTI
      ISPI=1./FLOAT(ISPLIT)
      ISP2I=.5*ISPI
      DAYI=1./86400.
      GRAV=9.806

!NG   INITIALIZE HRSRESTART! If greater than 0, hrs after start that startup is generated
      NRESEND=0
      NENDLESS=1
      OPEN (IURRT,FILE='hrsrestart.txt',STATUS='old',ERR=19) 
      DO NENDLESS=1,99
       READ (IURRT,*,END=19,ERR=19) HRSINNG
       IF (HRSINNG.GT.0.) THEN 
        NRESEND=NRESEND+1
        HRSRESTART(NRESEND)=HRSINNG
        IF 
     .   (NRESEND.GT.1.AND.HRSRESTART(NRESEND).LE.HRSRESTART(NRESEND-1))
     .   THEN 
         WRITE (*,*) 'hrsrestart.txt entries not in ascending order!'
         WRITE (*,*) 'Please rearrange them and resubmit.'
         STOP
        ENDIF
       ENDIF
      ENDDO
      CLOSE (IURRT)
 19   NRESTART=1
      IF (NRESEND.GT.0) THEN 
       DO NENDLESS=1,NRESEND
        WRITE(*,*)
     + 'startup file will be generated ',HRSRESTART(NENDLESS),
     + ' hrs after start.'
       ENDDO
      ELSE
        WRITE(*,*)
     + 'startup file will be generated at the end of run.'
      ENDIF
! 
      WRITE(IUPRT,21) DTI,DTE,ISPLIT,IRAMP,IYR,IMO,IDA,IHOUR ! JI
 21   FORMAT(
     . ' BAROCLINIC TIME STEP IS                 ',F10.4,' SECONDS',//,
     . ' BAROTROPIC TIME STEP IS                 ',F10.4,' SECONDS',//,
     . ' INTERNAL/EXTERNAL MODE SPLITTING IS     ',I10//,
     . ' NUMBER OF RAMP TIME STEPS               ',I10//,
     & ' THE STARTING TIME IS  ', 4i5//)  ! JI
 211   FORMAT(
     & ' HOT STARTING TIME IS AROUND ', 4i5//)
!
      WRITE (IUPRT,5002)TRACER,SEDTRAN,CHEMTRAN,SEDTYPE,PARTICLE,HYDTYPE
     +  ,WAVEDYN,SGW,WETEPS,WETMIN
 5002 FORMAT (//5X,'MODEL OPTIONS: TRACER   =',3X,A7,/20X,'SEDTRAN =',
     + 3X,A7,/20X,'CHEMTRAN =',3X,A7/20X,'SEDTYPE =',6X,A4,
     + /20X,'PARTICLE =',3X,A7,/20X,'HYDTYPE =',2X,A8,
     + /20X,'WAVEDYN =',2X,A8,
     & /20X,'SGW     =',2X,F8.2,         !WAD (hli,08/10/05)
     & /20X,'WETEPS  =',2X,F8.2,         !WAD
     & /20X,'WETMIN  =',2X,F8.2//)       !WAD
!
!
!
!-----------------------------------------------------------------------
!         TYPE OF RUN -
!      BAROTROPIC: 2-D CALCULATION (BOTTOM STRESS CALCULATED IN ADVAVE)
!      PROGNOSTIC: 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
!      By Quamrul 7/12/99 QA
!      TEMP_ONLY : 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
!                : ONLY TEMPERATURE IS CALCULATED
!      SALT_ONLY : 3-D CALCULATION (BOTTOM STRESS CALCULATED IN PROFU,V)
!                : ONLY SALINITY IS CALCULATED
!      DIAGNOSTIC: 3-D CALCULATION WITH T AND S HELD FIXED
!----------------------------------------------------------------------- 
!         3-D - TYPE OF MOMENTUM ADVECTION AND BOTTOM FRICTION
!      LINEAR    : ALL MOMENTUM ADVECTION NEGLECTED
!      NON-LINEAR: COMPLETE PHYSICS 
!----------------------------------------------------------------------  
!         3-D - TYPE OF ADVECTION INTEGRATION SCHEME
!      CENTRAL   : CENTRAL FINITE DIFFERENCE SCHEME
!      UPWIND    : UPWIND FINITE DIFFERENCE SCHEME
!      SMOLAR_R  : FINITE DIFFERENCE SCHEME DUE TO SMOLARKIEWICZ USING
!                  RECURSIVE FORMULATION FOR THE ANTIDIFFUSIVE VELOCITIES
!                      (MOMENTUM AND TURBULENCE ARE UPWIND)
!                      ( to be implemented, now central)
!      SMOLAR_2  : FINITE DIFFERENCE SCHEME DUE TO SMOLARKIEWICZ USING
!                  TWO PASSES FOR CORRECTIONS OF THE NUMERICAL DIFFUSION
!                      (MOMENTUM AND TURBULENCE ARE UPWIND)
!                      ( to be implemented, now central)
!----------------------------------------------------------------------  
!
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,10) NSTEPS,IPRINT,IPRTSTART,RESTAR,TOR,ADVECT,SCHEME
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,4)  DEV,VSX,JROW,VSY,IROW,
     .               PTU,PTV,PTW,PTAM,PTS,PTT,PRHO,PTQ2,PTL,PTKM,PTKH
      READ(IURUN,11)  (COM(I),I=1,80)

      TNDAYS=FLOAT(NSTEPS)*DAYI*DTI
      EDAY=SDAY+TNDAYS
!   PUT WARNING IF TOR IS PROGNOSTIC AND SGW COEFFICIENT IS SET 0.0
      IF(SGW.EQ.0.0.AND.TOR.EQ.'PROGNOSTIC ')THEN
       WRITE(IUPRT,*)'YOU HAVE SPECIFIED SGW = 0.0 AND TOR=PROGNOSTIC'
       WRITE(IUPRT,*)'IF SGW = 0.0, TOR SHOULD BE DIAGNOSTIC'
       WRITE(IUPRT,*)'PLEASE REFER TO USERS MANUAL FOR MORE DETAILS OF
     .  SGW COEFFICIENT'
       WRITE(IUPRT,*)'PLEASE CORRECT THIS AND RESUBMIT RUN'
       STOP
      END IF
!
!  WAVEHYD = 'NEGLECT' ==>  NO EFFECT OF WAVES ON BOTTOM FRICTION
!          = 'INLCUDE' ==>  WAVES EFFECT BOTTOM FRICTION COEFF.
!
!  NWAVE = NUMBER OF TIMESTEPS BETWEEN UPDATING NEW BOTTOM FRICTION COEFF.
!
!      READ(IURUN,3)  BFRIC,Z0B,NU,ALPHA,TLAG,NWAVE
!      Modified by Quamrul QA 12/28/98
!
!      Introduce ALPHA as Advection Time Scale for T and S
!      and TLAG for Friction Time Scale for Barotropic Radiation
!      Boundary Conditions (PCLAMP)

      DO J=1,JM
       DO I=1,IM
        VARBF(I,J)=1.0
        VARHF(I,J)=1.0
        VARUF(I,J)=1.0
       END DO
      END DO

      READ(IURUN,3)  OPTBFRIC,Z0B,NU,THETA,ALPHA,TLAG,NWAVE,BCTYPE

!NG11302010      IF(WAVEDYN.NE.'NEGLECT '.AND.NWAVE.NE.1) ! NG04272010
!NG11302010     .   PAUSE 'NWAVE<>1 NOT YET SUPPORTED FOR WAVE SIMULATIONS.'//
!NG11302010     .         'PLEASE SET TO 1 IN RUN_DATA.'

         IF(OPTBFRIC(7:10).EQ.'VARI')THEN
           OPEN (70,FILE='bfric2d.inp',STATUS='OLD')
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           READ(70,'(I5,E10.4)')NVARBF,BFRIC
           WRITE(IUPRT,'(I5,E10.4)')NVARBF,BFRIC
           IF(NVARBF.GT.MAXWET)THEN
             WRITE(IUPRT,56)NVARBF,MAXWET
             CALL SYSTEM ('rm gcm_temp*')
             STOP
           ENDIF
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           DO N=1,NVARBF
            READ(70,55)II,JJ,VARBF(II,JJ)
            WRITE(IUPRT,55)II,JJ,VARBF(II,JJ)
           END DO
           CLOSE(70)
           OPTBFRIC='VARIABLE'
         ELSEIF(OPTBFRIC.EQ.'MANNINGS-N')THEN
           IF(TOR.NE.'BAROTROPIC') PAUSE 'Mannings-N needs BAROTROPIC!'
           OPEN (70,FILE='bfric2d.inp',STATUS='OLD')
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           READ(70,'(I5,E10.4)')NVARBF,BFRIC ! Bfric is minimum here
           NVARBF=MAXWET ! and NVARBF is redefined as the whole IxJ space
           WRITE(IUPRT,'(A,I8,A,E10.4)')'REDEFINED NVARBF:',NVARBF,
     +                                   ' WITH MINIMUM BFRIC:',BFRIC
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           DO N=1,NVARBF
            READ(70,55)II,JJ,VARMN(II,JJ) ! VARMN is the actual Manning's number here, not a scaling factor
            WRITE(IUPRT,55)II,JJ,VARMN(II,JJ)
           END DO
           CLOSE(70)
         ELSE
           READ(OPTBFRIC,'(F10.3)')BFRIC
           OPTBFRIC='UNIFORM'
         ENDIF
!

                IF(BCTYPE.EQ.'       ')BCTYPE='CLAMPED'

      IF(BCTYPE.NE.'IRANDB '.AND.BCTYPE.NE.'OCLAMP '.AND.
     1BCTYPE.NE.'CLAMPED'.AND.BCTYPE.NE.'PCLAMP '.AND.
     2BCTYPE.NE.'RANDB  '.AND.BCTYPE.NE.'MIXED  ') THEN
         WRITE(IUPRT,6110) BCTYPE
          GOTO 9100
      END IF
      IF(BCTYPE.EQ.'PCLAMP '.AND.TLAG.EQ.0.0) THEN
         WRITE(IUPRT,6111) BCTYPE,TLAG
         GOTO 9100
      ENDIF
!
      NWAVECNT=0
!
      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,23) HORZMIX,OPTHORCON,HPRNU

         IF(OPTHORCON(7:10).EQ.'VARI')THEN
           OPEN (70,FILE='horcon2d.inp',STATUS='OLD')
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           READ(70,'(I5,E10.3)')NVARHF,HORCON
           WRITE(IUPRT,5)NVARHF
           IF(NVARHF.GT.MAXWET)THEN
             WRITE(IUPRT,57)NVARHF,MAXWET
             CALL SYSTEM ('rm gcm_temp*')
             STOP
           ENDIF
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           DO N=1,NVARHF
            READ(70,55)II,JJ,VARHF(II,JJ)
            WRITE(IUPRT,55)II,JJ,VARHF(II,JJ)
           END DO
           CLOSE(70)
          ELSE
           READ(OPTHORCON,'(E10.4)')HORCON
         ENDIF

      READ(IURUN,11)  (COM(I),I=1,80)
      READ(IURUN,23) VERTMIX,OPTUMOL,VPRNU
         IF(OPTUMOL(7:10).EQ.'VARI')THEN
           OPEN (70,FILE='umol2d.inp',STATUS='OLD')
           READ(70,11)  (COM(I),I=1,80)
           READ(70,'(I5,E10.4)')NVARUF,UMOL
           WRITE(IUPRT,'(I5,E10.4)')NVARUF,UMOL
           IF(NVARUF.GT.MAXWET)THEN
             WRITE(IUPRT,58)NVARUF,MAXWET
             CALL SYSTEM ('rm gcm_temp*')
             STOP
           ENDIF
           READ(70,11)  (COM(I),I=1,80)
           WRITE(IUPRT,11)  (COM(I),I=1,80)
           DO N=1,NVARUF
            READ(70,55)II,JJ,VARUF(II,JJ)
            WRITE(IUPRT,55)II,JJ,VARUF(II,JJ)
           END DO
           CLOSE(70)
          ELSE
           READ(OPTUMOL,'(E10.4)')UMOL
         ENDIF

 55   FORMAT(2I5,E10.3)
 56   FORMAT(/' NUMBER OF GRID OF VARIABLE BFRIC ',I5,' IS GREATER THAN
     * MAXIMUM ALLOWED.',I5// ' PLEASE FIX THIS AND RESUBMIT'//)
 57   FORMAT(/' NUMBER OF GRID OF VARIABLE HORCON ',I5,' IS GREATER THAN
     * MAXIMUM ALLOWED.',I5// ' PLEASE FIX THIS AND RESUBMIT'//)
 58   FORMAT(/' NUMBER OF GRID OF VARIABLE UMOL ',I5,' IS GREATER THAN
     * MAXIMUM ALLOWED.',I5// ' PLEASE FIX THIS AND RESUBMIT'//)
 4    FORMAT(2X,A3,2(4X,A1,I5),11(4X,A1))
 3    FORMAT(A10,5E10.3,I10,3X,A7)
 10   FORMAT(3I10,1X,A10,1X,A10,1X,A10,1X,A10)
 23   FORMAT(2A10,2E10.3)
!
!*****************************************************************
!
!  CONSERVATIVE TRACER INPUT SECTION
!  CONDRAT = TRACER DECAY RATE (1/day)
!  CONINIT = INITIAL TRACER CONCENTRATION (ASSUMED SPATIALLY CONSTANT)
!
      IF (TRACER.EQ.'INCLUDE') THEN
!qa     Add decay rates 11/16/01

        READ (IUT401,11) (COM(I),I = 1,80)
        READ (IUT401,5009)TROPT,TRICOPT,CONDRAT,CONINIT,THETAT,VSRATE,
     .                    ASOL,ASAL
 5009   FORMAT (A10,5X,A5,6F10.0)
!
        WRITE (IUPRT,11) (COM(I),I = 1,80)
        WRITE (IUPRT,5119)TROPT,TRICOPT,CONDRAT,CONINIT,THETAT,VSRATE,
     .                    ASOL,ASAL
 5119   FORMAT ('Tracer Option=',A10,/
     .          'Initial Data Option=',A5,1X,/
     .  'Background/Freshwater Decay Rate=',F15.5/
     .  'Initial Concentration of Tracer=',F15.5/
     .  'Temperature Coefficient Theta=',F15.5/
     .  'Settling rate of Particulate Bacteria=',F15.5/
     .  'Proportionallity Constant for Insolation =',F15.5/
     .  'Switch for Salinity Dependance =',F15.5)
!
!  CONVERT DECAY RATE FROM 1/day TO 1/s
!
!qa        CONDRAT=CONDRAT/86400.
!qa   We use decay rate function of T,S, and Io
!qa   The CONDRAT is now decay rate at T=20 deg C and in Freshwater
!qa   So we will change unit of CONSRAT in advcon.f
!qa   CONDRAT is = 0.8 (Mancini, 1978) by Quamrul QA 11/12/01

! CHECK THE VALIDITY OF INPUT FILE
      IF(TROPT.NE.'PATHOGEN  '.AND.TROPT.NE.'SIMPDECAY '.AND.
     .   TROPT.NE.'CONSERV   ')THEN
       WRITE(IUPRT,5120)TROPT
       STOP
 5120 FORMAT('TRACER OPTION', 5X,A10,5X,' IS NOT A VALID OPTION'/
     . 'PLEASE CORRECT THIS AND RESUBMIT THE RUN')
      END IF
      IF(TROPT.EQ.'PATHOGEN  '.AND.THETAT.EQ.0.0)THEN
       WRITE(IUPRT,*)'FOR PATHOGEN SIMULATIONS THETAT CANNOT BE ZERO'
       WRITE(IUPRT,*)'FIX THIS INPUT AND RESUBMIT'
       STOP
      END IF

      ENDIF
!
!**********************************************************************
!
!  SEDIMENT TRANSPORT INPUT SECTION
!
      IF (SEDTRAN.EQ.'INCLUDE') THEN
!
!  COHESIVE SEDIMENT INPUT
!
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5010)NSEDBEG,NSBED
 5010     FORMAT (3I10)
!
          XNSBED=FLOAT(NSBED)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5010)NSEDBEG,NSBED
!
!  WS1 = ADEP * (C*G) ** DEPEXP
!
!  ADEP = CONSTANT, CLASS 1 SETTLING SPEED IN um/s
!
!  PDEPFORM = 'KRONE'
!           = 'PARTH'
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5019)ADEP,DEPEXP,TCRDEP,PDEPFORM
 5020     FORMAT (8F10.0)
 5019     FORMAT (3F10.0,5X,A5)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5019)ADEP,DEPEXP,TCRDEP,PDEPFORM
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5018)A0IN,RESEXP,EXPM,VARIA0N
 5018     FORMAT (3F10.0,3X,A7)
!
!
          IA0=0
          IF (VARIA0N.EQ.'INCLUDE') IA0=1
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020)A0IN,RESEXP,EXPM
!
          READ (IUT402,11) (COM(I),I = 1,80)
!
!  MODIFIED BY C.K.Z. ON 4/21/97
!
!  ADDED Z0 FOR WAVES
!
          READ (IUT402,5017)DENCOH,VARIBULK,P0(1),VARIP0,BFCOH,Z0BCOH,
     +                     Z0WAVE
 5017     FORMAT (F10.0,3X,A7,F10.0,3X,A7,3F10.0)
!
          IP0=0
          IF (VARIP0.EQ.'INCLUDE') IP0=1
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5017)DENCOH,VARIBULK,P0(1),VARIP0,BFCOH,Z0BCOH,
     +                      Z0WAVE
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) (FTIME(LL),LL=1,LAYMAX)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) (FTIME(LL),LL=1,LAYMAX)
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) (TSED0IN(LL),LL=1,LAYMAX)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) (TSED0IN(LL),LL=1,LAYMAX)
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) (TAUCR(LL),LL=1,LAYMAX)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) (TAUCR(LL),LL=1,LAYMAX)
!
!  INPUT INITIAL SUSPENDED SEDIMENT CONCENTRATIONS
!
!  ASSUMED TO BE SPATIALLY CONSTANT
!
          READ (IUT402,11) (COM(I),I = 1,80)
          READ (IUT402,5020) CSI(1)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) CSI(1)
        ENDIF
!
!  NON-COHESIVE SEDIMENT INPUT
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
!
          READ (IUT403,11) (COM(I),I = 1,80)
          READ (IUT403,5010)NSEDBEG,NSBED
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5010)NSEDBEG,NSBED
!
          READ (IUT403,11) (COM(I),I = 1,80)
          READ (IUT403,5021)WS2,DENNON,VARIBULK,SUSARM,BEDTHI
 5021     FORMAT (2F10.0,3X,A7,2F10.0)
!
!  CONVERT INPUT NON-COHESIVE BED THICKNESS FROM cm TO m
!
          BEDTHI=BEDTHI/100.
!
          XNSBED=FLOAT(NSBED)
!
!
!  INPUT INITIAL SUSPENDED SEDIMENT CONCENTRATIONS
!
!  ASSUMED TO BE SPATIALLY CONSTANT
!
          READ (IUT403,11) (COM(I),I = 1,80)
          READ (IUT403,5020) CSI(2)
!
          WRITE (IUPRT,11) (COM(I),I = 1,80)
          WRITE (IUPRT,5020) CSI(2)
!
          IF (SEDTYPE.EQ.'SAND') CSI(1)=CSI(2)
        ENDIF
!
!
!  CONVERT INITIAL CONCENTRATIONS FROM mg/l TO g/cm**3
!
        DO 5022 K=1,KSED
          CSI(K)=CSI(K)/1000000.
 5022   CONTINUE
!
!  INITIALIZE COUNTERS
!
        NSEDCT=0           !Check for 0 ini. Pravi in STCODE 
        NBLOW=0           !Check for 0 ini.
        N24CNT=0           !Check for 0 ini.
!
        NHR=NINT(3600./DTI)
        N24HR=24*NHR
      ENDIF
!
!**********************************************************************
!
!  CHEM TRANSPORT INPUT
!
      IF (CHEMTRAN.EQ.'INCLUDE') THEN
!
!  CHEMI(1) = INITIAL WATER COL. CHEM CONC., CLASS 1 (COHESIVE) (ug/l)
!  CHEMI(2) = INITIAL WATER COL. CHEM CONC., CLASS 2 (NON-COHESIVE) (ug/l)
!  NCHEMLAY = NO. OF LAYERS IN CHEM BED MODEL (EXCLUDING SURFACE LAYER)
!  CHEMTHIK = THICKNESS OF LAYERS IN CHEM BED MODEL (cm)
!  CHEMACT  = ACTIVE LAYER THICKNESS (cm)
!
!  CHEMDRAT1 =  PARTICLE-BOUND CHEMICAL DECAY RATE (1/day), CLASS 1
!  CHEMDRAT2 =  PARTICLE-BOUND CHEMICAL DECAY RATE (1/day), CLASS 2
!
!
!  COHESIVE PARTICLE-BOUND TRACER
!
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          READ (IUT404,11) (COM(I),I = 1,80)
          READ (IUT404,5030)CHEMI(1),NCHEMLAY,CHEMTHIK,CHEMACT,CHEMDRAT1
 5030     FORMAT (F10.0,I10,4F10.0)
        ENDIF
!
!  NON-COHESIVE PARTICLE-BOUND TRACER
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          READ (IUT405,11) (COM(I),I = 1,80)
          READ (IUT405,5030)CHEMI(2),NCHEMLAY,CHEMTHIK,CHEMACT,CHEMDRAT2
!
          IF (SEDTYPE.EQ.'SAND') CHEMI(1)=CHEMI(2)
        ENDIF
!
!  CONVERT DECAY RATE FROM 1/day TO 1/s
!
        CHEMDRAT1=CHEMDRAT1/86400.
        CHEMDRAT2=CHEMDRAT2/86400.
!
!  CALC. NUMBER OF LAYERS IN ACTIVE LAYER
!
        NACTLAY=NINT(CHEMACT/CHEMTHIK)
!
!  INPUT INITIAL CHEM BED CONCENTRATIONS (ug CHEM/g SOLIDS)
!
        If (RESTAR.EQ.'COLD START') Then
          OPEN (UNIT=85,FILE='bed_chemic',FORM='FORMATTED')
          DO 5029 N=1,NCHEMLAY
            DO 5029 I=1,IM
              READ (85,5028) (CBEDCHEM(N,I,J),J=1,JM)
 5028         FORMAT (10F8.0)
 5029     CONTINUE
          CLOSE (85)
        ENDIF
!
!  CONVERT FROM ug/l TO ug/cm**3
!
        CHEMI(1)=CHEMI(1)/1000.
        CHEMI(2)=CHEMI(2)/1000.
!
      ENDIF
!
!***********************************************************
!
!  FOR PARTICLE TRACKING     
!
!  NFREQ = FREQUENCY OF PARTICLE RELEASE (timesteps)
!  NPART = NUMBER OF PARTICLES PER RELEASE
!  NCONV = TOTAL NUMBER OF RELEASES BEFORE CONVERSION OF PARTICLES
!          TO CONCENTRATION (timesteps)
!
!  NOTE:  TOTAL NO. OF PARTICLES IN SYSTEM = NPART * NCONV
!
!  IRELST = BEGINNING TIMESTEP OF PARTICLE RELEASE
!  NPCLASS = NUMBER OF PARTICLE CLASSES
!
!  NSOURCE = NO. OF SOURCES OF PARTICLES
!
!  isource = i of particle input
!  jsource = j of particle input
!  ksource = k-level of particle input, for testing only
!
      IF (PARTICLE.EQ.'INCLUDE') THEN
!
        Read (IUT406,11) (COM(I),I = 1,80)
        WRITE (IUPRT,11) (COM(I),I = 1,80)
        READ (IUT406,5810)NFREQ,NPART,IRELST,IRELEND,
     +                          NSOURCE
        NCONV=(IRELEND-IRELST)/NFREQ+1       ! nkim 061098
        WRITE (IUPRT,5810)NFREQ,NPART,NCONV,IRELST,IRELEND,
     +                          NSOURCE
 5810   FORMAT (9I8)
!
        DO 5809 MM=1,NSOURCE
          READ (IUT406,5901)ISOURCE(MM),JSOURCE(MM),KSOURCE(MM)
          WRITE (IUPRT,5901)ISOURCE(MM),JSOURCE(MM),KSOURCE(MM)
 5809   CONTINUE
 5901   FORMAT (9I8)
!
        IF (NPART.GT.NPARTM) THEN
          WRITE (IUPRT,5811)NPART,NPARTM
 5811     FORMAT (/5X,'EXECUTION STOPPED BECAUSE NPART > NPARTM,'
     +     /6X,'NPART =',I6,3X,'NPARTM =',I6,
     +     /6X,'PLEASE SET NPART <= NPARTM AND RESUBMIT')
          STOP
        ENDIF
!
        IF (NCONV.GT.NCONVM) THEN
          WRITE (IUPRT,5812)NCONV,NCONVM
 5812     FORMAT (/5X,'EXECUTION STOPPED BECAUSE NCONV > NCONVM,'
     +     /6X,'NCONV =',I6,3X,'NCONVM =',I6,
     +     /6X,'PLEASE SET NCONV <= NCONVM AND RESUBMIT')
          STOP
        ENDIF
!
!
!  INTCONV = TIMESTEP AT WHICH TO INITIALLY START CONVERTING PARTICLES
!            TO CONSERVATIVE TRACER
!
        INTCONV=NCONV*NFREQ+IRELST
!
!  OPEN OUTPUT FILE FOR PARTICLES (TESTING ONLY)
!
        OPEN (UNIT=32,FILE='part_location',FORM='FORMATTED')
!
      ENDIF
!
!********************************************************************
!
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
!     READ(IURUN,7)  JHM,IAVGE
!     WRITE(IUPRT,7) JHM,IAVGE

      READ(IURUN,7101) JHM,AVGE,IPLTFORM,PLTZERO,OPTAVG,OFORM
      WRITE(IUPRT,710) JHM,AVGE,IPLTFORM,PLTZERO,OPTAVG,OFORM
      IF(OPTAVG.EQ.'    '.OR.OPTAVG.EQ.'NDTI')THEN
       IAVGE=IFIX(AVGE)
       IPLTZERO=IFIX(PLTZERO)
      ELSE IF(OPTAVG.EQ.'HOUR') THEN
       IF(AMOD(AVGE*3600.0,DTI).NE.0.0) THEN
        WRITE(IUPRT,*)'AVERAGING INTERVAL IS NOT AN INTEGER MULTIPLE
     .OF DTI: PLEASE ADJUST AVGE AND RESUBMIT'
        STOP
       ELSE
        IAVGE=IFIX(AVGE*3600/DTI)
        IPLTZERO=IFIX(PLTZERO*3600/DTI)
       END IF
      ELSE IF(OPTAVG.NE.'    '.OR.OPTAVG.NE.'HOUR'
     +.OR.OPTAVG.NE.'NDTI') THEN
       WRITE(IUPRT,*)'AVERAGING UNIT IS NOT SPECIFIED CORRECTLY'
       WRITE(IUPRT,*)'PLEASE CORRECT THIS AND RESUBMIT THE RUN'
       STOP
      END IF
! gcmplt output format (either regular binary or netCDF format)
       IF (OFORM.EQ.'CDF'.OR.OFORM.EQ.'cdf') THEN
         CDFOUT = .TRUE.
         WRITE(IUPRT,*)'netCDF output is on'
       ELSE
         CDFOUT = .FALSE.
         WRITE(IUPRT,*)'Regular binary gcmplt output is on'
       END IF
       WRITE(IUPRT,*)'CDFOUT was set to ',CDFOUT

 710  FORMAT(I10,F10.3,6X,A4,F10.2,1X,A4,2X,A3)
 7101 FORMAT(I10,F10.0,6X,A4,F10.0,1X,A4,2X,A3)


      IF(JHM.EQ.0) GO TO 200
      IF(JHM.GT.IHISTM) THEN
       WRITE(IUPRT,60) JHM,IHISTM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  60  FORMAT(//' JHM    (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' IHISTM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      DEI=1./FLOAT(IAVGE)

      IF(IPLTFORM.EQ.'    '.OR.IPLTFORM.EQ.'USER')THEN
       READ(IURUN,6)  (IHIST(I,2),I=1,JHM)
      ELSE IF(IPLTFORM.EQ.'AUTO')THEN
       DO I=1,JHM
        IHIST(I,2)=IPLTZERO+I*IAVGE
       END DO
      ELSE IF(IPLTFORM.NE.'    '.AND.IPLTFORM.NE.'USER'.AND.
     .    IPLTFORM.NE.'AUTO') THEN
      ! FIXME: Next 2 lines write in run_data (IURUN), should write in gcmprt (IUPRT) (A.Jordi)
       WRITE(IURUN,*)'YOUR SETTING FOR GCMPLT OUTPUT IS INCORRENT'
       WRITE(IURUN,*)'YOU SPECIFIED IPLTFFORM AS ',IPLTFORM
       STOP
      END IF

      WRITE(IUPRT,6) (IHIST(I,2),I=1,JHM)
 7    FORMAT(8I10)
 6    FORMAT(10I8)
!
 200  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,7)  ISKILL
      WRITE(IUPRT,7) ISKILL
      IF(ISKILL.EQ.0) THEN
       SKILLI=1.0
       GO TO 203
      ELSE
       SKILLI=1./FLOAT(ISKILL)
      ENDIF
      READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  EPTS
      WRITE(IUPRT,5) EPTS
 5    FORMAT(16I5)
      IF(EPTS.EQ.0) GO TO 201
      IF(EPTS.GT.EPTSM) THEN
       WRITE(IUPRT,61) EPTS,EPTSM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  61  FORMAT(//' EPTS  (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' EPTSM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      READ(IURUN,5)  (INXIE(I),INXJE(I),I=1,EPTS)
      WRITE(IUPRT,5) (INXIE(I),INXJE(I),I=1,EPTS)
201   READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  VPTS
      WRITE(IUPRT,5) VPTS
      IF(VPTS.EQ.0) GO TO 202
      IF(VPTS.GT.VPTSM) THEN
       WRITE(IUPRT,62) VPTS,VPTSM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  62  FORMAT(//' VPTS  (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' VPTSM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      READ(IURUN,5)  (INXIV(I),INXJV(I),I=1,VPTS)
      WRITE(IUPRT,5) (INXIV(I),INXJV(I),I=1,VPTS)
 202  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IURUN,5)  FPTS
      WRITE(IUPRT,5) FPTS
      IF(FPTS.EQ.0) GO TO 203
      IF(FPTS.GT.FPTSM) THEN
       WRITE(IUPRT,63) VPTS,VPTSM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  63  FORMAT(//' FPTS  (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' FPTSM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      READ(IURUN,49)  (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
 
      WRITE(IUPRT,49) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
 49   FORMAT(4(2I5,1X,A4,I5))
!
 203  READ(IURUN,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
!     READ(IURUN,7101) JHM,AVGE,IPLTFORM,PLTZERO,OPTAVG
      READ(IURUN,7102) JTM,AVGE,ITRNFORM,TRNZERO,OPTAVG,IWET
 7102 FORMAT(I10,F10.0,6X,A4,F10.0,1X,A4,I5)
      IF(OPTAVG.EQ.'    '.OR.OPTAVG.EQ.'NDTI')THEN   
         NPLPF=IFIX(AVGE)
         IZERO=IFIX(TRNZERO)
      ELSE IF(OPTAVG.EQ.'HOUR')THEN
         IF(AMOD(AVGE*3600.0,DTI).NE.0.0)THEN
         WRITE(IUPRT,*)'AVERAGING INTERVAL FOR GCM_TRAN IS NOT AN INTEG
     .ER MULTIPLE OF DTI: PLEASE ADJUST AVGE AND RESUBMIT'
         STOP
         ELSE
          NPLPF=IFIX(AVGE*3600/DTI)
          IZERO=IFIX(TRNZERO*3600/DTI)
         END IF
      ELSE IF(OPTAVG.NE.'   '.OR.OPTAVG.NE.'HOUR'
     .     .OR.OPTAVG.NE.'NDTI')THEN
         WRITE(IUPRT,*)'AVERAGING UNIT IS SPECIFIED INCORRECTLY'
         WRITE(IUPRT,*)'PLEASE CORRECT THIS AND RESUBMIT THE RUN'
         STOP
      END IF

      WRITE(IUPRT,7103) JTM,NPLPF,ITRNFORM,IZERO,IWET
 7103 FORMAT(2I10,6X,A4,2I10)
      IF(IWET.EQ.0)THEN
       WRITE(IUPRT,*)'WQ INFO FOR ENTIRE GRID OPTION SELECTED'
      ELSE
       WRITE(IUPRT,*)'WQ INFO FOR WET-GRID OPTION SELECTED'
      END IF
      IF(ITRNFORM.EQ.'    '.OR.ITRNFORM.EQ.'USER')
     . WRITE(IUPRT,*)'USER SPECIFIED TIME BREAKS FOR DUMP USED'
      IF(ITRNFORM.EQ.'AUTO')
     . WRITE(IUPRT,*)'ECOM WILL GENERATE THE TIME BREAKS VIA JTM,NPLPF,
     . AND IZERO'
      IF(JTM.EQ.0) GO TO 204
      IF(JTM.GT.ITRACM) THEN
       WRITE(IUPRT,64) JTM,ITRACM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  64  FORMAT(//' JTM    (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' ITRACM (=',I4,') SPECIFIED IN COMDECK'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      FLTWT=1./FLOAT(NPLPF)
      IF(ITRNFORM.EQ.'    '.OR.ITRNFORM.EQ.'USER')THEN
      READ(IURUN,6) (ITRAC(I,2),I=1,JTM)
      ELSE
       DO I=1,JTM
        ITRAC(I,2)=IZERO+I*NPLPF
       END DO
      END IF
      WRITE(IUPRT,6) (ITRAC(I,2),I=1,JTM)
 204  CONTINUE
!
      IF(TOR.NE.'BAROTROPIC' .AND. TOR.NE.'PROGNOSTIC' .AND.
     .   TOR.NE.'TEMP_ONLY '.AND.TOR.NE.'SALT_ONLY '.AND.
     .   TOR.NE.'DIAGNOSTIC') THEN
      WRITE(IUPRT,24) TOR
 24   FORMAT(//' TYPE OF RUN (TOR=',A10,') IS SPECIFIED INCORRECTLY',/
     .         ' PLEASE FIX AND RESUBMIT'//)
      CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
!
      IF (ADVECT.NE.'LINEAR    ' .AND. ADVECT.NE.'NON-LINEAR') THEN
      WRITE(IUPRT,25) ADVECT
  25  FORMAT(//' TYPE OF MOMENTUM ADVECTION (ADVECT=',A10,') IS',
     .'SPECIFIED INCORRECLTY',/' PLEASE FIX AND RESUBMIT'//)
      CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
!
      IF(SCHEME.NE.'CENTRAL   ' .AND. SCHEME.NE.'UPWIND    ' .AND.
     .   SCHEME.NE.'SMOLAR_R  ' .AND. SCHEME.NE.'SMOLAR_2  ') THEN
      WRITE(IUPRT,34) SCHEME
 34   FORMAT(//' T&S ADVECTION SCHEME (SCHEME=',A10,') IS SPECIFIED',
     .   ' INCORRECTLY',/' PLEASE FIX AND RESUBMIT'//)
      CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
!
      IF (HORZMIX.NE.'CLOSURE   ' .AND. HORZMIX.NE.'CONSTANT  ') THEN
      WRITE(IUPRT,26) HORZMIX
  26  FORMAT(//' TYPE OF HORIZONTAL MIXING (HORZMIX=',A10,') IS',
     .         ' SPECIFIED INCORRECTLY',/
     .         ' PLEASE FIX AND RESUBMIT'//)
      CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
!
      IF (VERTMIX.NE.'CLOSURE   ' .AND. VERTMIX.NE.'CONSTANT  '.AND.
     .  VERTMIX.NE.'EMPIRICAL ') THEN
      WRITE(IUPRT,27) VERTMIX
  27  FORMAT(//' TYPE OF VERTICAL MIXING (VERTMIX=',A10,') IS',
     .         ' SPECIFIED INCORRECTLY',/
     .         ' PLEASE FIX AND RESUBMIT'//)
      CALL SYSTEM ('rm gcm_temp*')
      STOP
      END IF
!
!qa      IF (ITRAC(1,2).NE.0 .AND. NPLPF.NE.0) THEN
!qa       IF (RESTAR.NE.'HOT START ' .AND. TOR.NE.'BAROTROPIC' .AND.
!qa     .     ITRAC(1,2).EQ.NPLPF) THEN
!qa       WRITE(IUPRT,28) ITRAC(1,2),NPLPF,RESTAR,TOR
!qa  28   FORMAT(//' THE FIRST TRANSPORT COMPUTATION INTERVAL (ITRAC)',I10,/
!qa     .          ' CANNOT BE EQUAL TO THE NUMBER OF POINTS IN NPLPF',I10,/
!qa     .          ' IN A ',A10,' 3-DIMENSIONAL (',A10,') RUN',/
!qa     .          ' PLEASE FIX AND RESUBMIT'//)
!qa      CALL SYSTEM ('rm gcm_temp*')
!qa       STOP
!qa       END IF
!qa      END IF
!
      IF(TOR.EQ.'BAROTROPIC') THEN
      WRITE(IUPRT,14) TOR
      ELSE
      WRITE(IUPRT,13) TOR
      END IF
 13   FORMAT(/' THIS IS A THREE DIMENSIONAL MODEL RUN',2X,A10/)
 14   FORMAT(/' THIS IS A TWO DIMENSIONAL MODEL RUN',2X,A10/)
!
      WRITE(IUPRT,22) ADVECT
 22   FORMAT(/' THIS SIMULATION HAS ',A10,' MOMENTUM ADVECTION '/)
!
      WRITE(IUPRT,222) SCHEME
 222  FORMAT(/' THIS SIMULATION USES ',A10,' DIFFERENCING FOR T&S'/)
!
      IF(HORZMIX.EQ.'CLOSURE   ') THEN
      WRITE(IUPRT,29) HORZMIX,HORCON,HPRNU
      ELSE
      WRITE(IUPRT,31) HORZMIX,HORCON,HPRNU
      END IF
 29   FORMAT(/' THIS SIMULATION HAS ',A10,' HORIZONTAL MIXING ',
     . ' HORCON = ',1PE10.3,' HPRNU = ',1PE10.3/)
 31   FORMAT(/' THIS SIMULATION HAS ',A10,' HORIZONTAL MIXING ',
     . ' CONSTANT = ',1PE10.3,'m**2/s  HPRNU = ',1PE10.3/)
!
      IF(VERTMIX.EQ.'CLOSURE   ') THEN
      WRITE(IUPRT,32) VERTMIX,UMOL,VPRNU
      ELSE
      WRITE(IUPRT,33) VERTMIX,UMOL,VPRNU
      END IF
 32   FORMAT(/' THIS SIMULATION HAS ',A10,' VERTICAL MIXING ',
     . ' UMOL = ',1PE10.3,' VPRNU = ',1PE10.3/)
 33   FORMAT(/' THIS SIMULATION HAS ',A10,' VERTICAL MIXING ',
     . ' CONSTANT = ',1PE10.3,'m**2/s  VPRNU = ',1PE10.3/)
!
      CALL ZEROES(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
!

!NG01052010  SET SOME PARAMETERS FOR PRINTING/SAVING ACCORDINGLY
!NG01052010  (MOVED THEM FROM PUTCDF TO COMDECK)
!NG  save wave output
      IF(WAVEDYN.NE.'NEGLECT')THEN
       PTWAVE = 'Y'
!NG       WRITE(IUPRT,'A') 'Enabling WAVE printing and saving'
      ELSE
       PTWAVE = 'N'
      END IF
!NG
      IF(RESTAR.EQ.'COLD START') THEN
       CALL SETDOM
!hli
       INTX_HOT = INTX ! NG09212011 FOR MULTIPLE RESTARTS IF COLD START
       ISTART=INTX+1
       IEND=INTX+NSTEPS
!hli
       CALL BCDATA(RESTAR)
       CALL TANDS
       CALL DENS
       CALL ARCHIVE !NG03182010 FIX COLD-START IC AVERAGING (STILL NO PRINTING, IF NO AVERAGING)
      ELSE
!QA
!QA    IF RUN TYPE "HOT START" AND "EXTERNAL" and IWET is 1 then  "WET GRID"
!QA
       IF(HYDTYPE.EQ.'EXTERNAL'.AND.IWET.EQ.1)
     . CALL SETDOM
!
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
!
!  NO HOT START CAPABILITY FOR PARTICLE TRACKING (C.K.Z. 4/22/96)
!
!ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
!
      CALL NETCDF_READ_RESTART(ADVUA,ADVVA,ADVUU,ADVVV) ! A.Jordi
!sfan
!sfan     hot start time step
!
      INTX_HOT = INTX
      write(*,*)'HOT READ, INTX', INTX_HOT
      IF(VERTMIX.EQ.'CONSTANT  ') UMOL=0.0 ! NG03182010 OTHERWISE HOTSTART ISSUE, BECAUSE UMOL->KH
!
!********************************************************************
!
!  FOR PARTICLE TRACKING (corner_loc INPUT IN setdom FOR COLD START)
!
      IF (PARTICLE.EQ.'INCLUDE') THEN
        IRELST=IRELST+INTX
!
!  OPEN CORNER LOCATIONS FILE
!
!  CORNER LOCATION CONVENTION:  XCOR(I,J) = x(i-1/2, j-1/2) 
!     (LOWER LEFT-HAND CORNER)  YCOR(I,J) = y(i-1/2, j-1/2)
!
        OPEN (UNIT=33,FILE='corner_loc',FORM='FORMATTED')
!
        DO 2300 N=1,1000000
          READ (33,*,END=2310)I,J,XCOR(I,J),YCOR(I,J)
 2300   CONTINUE
!
!  CALC. H1 AND H2 AT ELEMENT INTERFACES
!
 2310   DO 2320 J=1,JM
          DO 2320 I=1,IM
            H1P(I,J)=SQRT((XCOR(I+1,J)-XCOR(I,J))**2.+
     +                        (YCOR(I+1,J)-YCOR(I,J))**2.)
!
            H2P(I,J)=SQRT((XCOR(I,J+1)-XCOR(I,J))**2.+
     +                        (YCOR(I,J+1)-YCOR(I,J))**2.)
 2320   CONTINUE
        CLOSE(33)
      ENDIF
!
!*********************************************************************
!

!hli
       ISTART=INTX+1
       IEND=INTX+NSTEPS
!hli
!      SDAY=SDAY+FLOAT(INTX)*DAYI*DTI     ! NKIM (11/06/01)
!      EDAY=SDAY+TNDAYS                  ! NKIM (11/06/01)
       EDAY=SDAY+TNDAYS+FLOAT(INTX)*DAYI*DTI
! RESET IYR,IMO,IDA FOR PTIDE IF HOT START !hli (04/04/01)
!NGCAREFUL:NEEDS STARTING HOUR FOR THIS RUN       HOUR=AMOD(SDAY,1.0)*24.0  
       HOUR=AMOD(AMOD(SDAY,1.0)+FLOAT(INTX)*DAYI*DTI,1.0)*24 ! CNG
!NG04212008       IHOUR1=IFIX(HOUR)
       IHOUR1=NINT(HOUR) ! To the closest hour for rounding errors.
!      KD=JINT(SDAY)
!svv       KD=JINT(SDAY+FLOAT(INTX)*DAYI*DTI)
       KD=SDAY+FLOAT(INTX)*DAYI*DTI
       CALL CDAY(IDA1,IMO1,IYR1,KD,2)
!hli
       TMP_IYR=IYR
       TMP_IMO=IMO
       TMP_IDA=IDA
       TMP_IHOUR=IHOUR
!
       IYR=IYR1
       IMO=IMO1
       IDA=IDA1
       IHOUR=IHOUR1
!hli
       WRITE(IUPRT,211) IYR,IMO,IDA,IHOUR
       CALL BCDATA(RESTAR)
!hli
       IYR=TMP_IYR
       IMO=TMP_IMO
       IDA=TMP_IDA
       IHOUR=TMP_IHOUR
!hli
      ENDIF

      call RDCORLOC ! NG04242009

!NG   HARDWIRE
!NG      DO I=73,108
!NG       DO J=3,5
!NG        H1(I,J)=H1(I,6)  
!NG        H2(I,J)=H2(I,6)  
!NG        ANG(I,J)=ANG(I,6)  
!NG        COR(I,J)=2.*7.292E-5*SIN(COR(I,6)*2.*3.141593/360.)*FSM(I,J)
!NG       ENDDO
!NG      ENDDO
!NG      DO I=73,108
!NG       DO J=3,5
!NG        ART(I,J)=H1(I,J)*H2(I,J)
!NG        ARU(I,J)=.25E0*(H1(I,J)+H1(I-1,J))*(H2(I,J)+H2(I-1,J))
!NG        ARV(I,J)=.25E0*(H1(I,J)+H1(I,J-1))*(H2(I,J)+H2(I,J-1))
!NG       ENDDO
!NG      ENDDO
!NG   HARDWIRE

!
!  INITIALIZE SEDIMENT TRANSPORT VARIABLES 
!
        IF (SEDTRAN.EQ.'INCLUDE') CALL SEDIC(RESTAR)
!
!**********************************************************************
!
      IF(OPTBFRIC.EQ.'MANNINGS-N')THEN
       CALL UPDATE_CBCM ! CNG2014 Manning's N
      ELSE
       CALL UPDATE_CBC ! CNG04262011 MOVED TO SUBROUTINE TOO !!!
      ENDIF
!NG04262011      Do 6030 J = 1, JM
!NG04262011        Do 6020 I = 1, IM
!NG04272010 W&D BASED ON H1          If (H(I,J).GT.0.0) Then
!NG04272010 W&D AND A MINIMUM DEPTH OF 10cm IS ALLOWED FOR CBC CALCULATION.
!NG04262011          If (FSM(I,J).NE.0.0) Then
!
!  VARIABLE BOTTOM FRICTION FOR SEDIMENT TRANSPORT
!
!NG04262011            IF (SEDTRAN.EQ.'INCLUDE') THEN
!
!  COHESIVE ELEMENTS
!
!NG04262011              IF (IBMSK(I,J).EQ.0) THEN
!NG04262011                Z0 = Z0BCOH
!NG04262011                CBCMIN = BFCOH
!NG04262011                IF (TOR.EQ.'BAROTROPIC') THEN
!NG04262011                  CBC(I,J)=BFCOH*FSM(I,J)
!NG04262011                ELSE
!qa                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((ZZ(KBM1)-Z(KB))*
!NG04262011                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
!NG04262011     +                      AMAX1(DT(I,J),1.E-1,WETMIN)/Z0)**2)*FSM(I,J) ! NG W&D
!NG04262011                ENDIF
!
!  NON-COHESIVE ELEMENTS
!
!NG04262011              ELSE
!NG04262011                Z0 = Z0B
!NG04262011                CBCMIN = BFRIC
!NG04262011                IF (TOR.EQ.'BAROTROPIC') THEN
!NG04262011                  CBC(I,J)=BFRIC*FSM(I,J)
!NG04262011                ELSE
!qa                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((ZZ(KBM1)-Z(KB))*
!NG04262011                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
!NG04262011     +                      AMAX1(DT(I,J),1.E-1,WETMIN)/Z0)**2)*FSM(I,J) ! NG W&D
!NG04262011                 ENDIF
!NG04262011               ENDIF
!NG04262011            ELSE
!NG04262011              Z0 = Z0B
!NG04262011              CBCMIN = BFRIC
!NG04262011              IF (TOR.EQ.'BAROTROPIC') THEN
!NG04262011                CBC(I,J)=BFRIC*FSM(I,J)
!NG04262011              ELSE
!qa                CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((ZZ(KBM1)-Z(KB))*
!NG04262011                CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
!NG04262011     +                      AMAX1(DT(I,J),1.E-1,WETMIN)/Z0)**2)*FSM(I,J) ! NG W&D
!NG04262011              ENDIF
!NG04262011            ENDIF
!NG04262011          END IF
!NG04262011          CBC(I,J) = AMIN1(CBC(I,J)*VARBF(I,J),1.) ! FACTOR IN 2D VARIABLE BFRIC, AND DO NOT ALLOW CBC TO EXCEED UNITY.
!NG04232008         IF(VARBF(I,J).GT.0.0)WRITE(IUPRT,*)'VARBF',VARBF(I,J)
!NG04262011         IF(VARBF(I,J).NE.1.)WRITE(IUPRT,*)'VARBF',VARBF(I,J)
!NG04262011 6020   CONTINUE
!NG04262011 6030 CONTINUE

!
!  SET INITIAL TRACER CONCENTRATIONS
!
      IF (TRACER.EQ.'INCLUDE'.AND.RESTAR.EQ.'COLD START') THEN
        IF(TRICOPT.NE.'DATA ')THEN
        DO 6034 K=1,KBM1
          DO 6034 J=2,JMM1
            DO 6034 I=2,IMM1
!NG              IF (H(I,J).GT.0.0) THEN
              IF (H1(I,J).GT.0.0) THEN ! BASE IT ON H1
                CONC1(I,J,K)=CONINIT
                CONC1B(I,J,K)=CONINIT
              ENDIF
 6034   CONTINUE
        ELSE
        OPEN(402,FILE='init_tracer',STATUS='OLD')
        DO N=1,IM*JM
         READ(402,'(2I5,100F5.0)',END=6033)I,J,(CONC1(I,J,K),K=1,KB)
        END DO
 6033   WRITE(IUPRT,*)'NUM OF CONS. TRACER INITIAL CONDITION',N-1
        DO K=1,KBM1
         DO J=2,JMM1
          DO I=2,IMM1
           CONC1B(I,J,K)=CONC1(I,J,K)
          END DO
         END DO
        END DO
        END IF
      ENDIF
!
!  SET INITIAL SEDIMENT CONCENTRATIONS
!
      IF (SEDTRAN.EQ.'INCLUDE'.AND.RESTAR.EQ.'COLD START') THEN
        DO 6035 K=1,KBM1
          DO 6035 J=2,JMM1
            DO 6035 I=2,IMM1
!NG              IF (H(I,J).GT.0.0) THEN
              IF (H1(I,J).GT.0.0) THEN ! BASE IT ON H1
                IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
                  CSED1(I,J,K)=CSI(1)
                  CSED1B(I,J,K)=CSI(1)
                ENDIF
!
                IF (SEDTYPE.EQ.'SAND') THEN
                  CSED2(I,J,K)=CSI(1)
                  CSED2B(I,J,K)=CSI(1)
                ELSE
                  CSED2(I,J,K)=CSI(2)
                  CSED2B(I,J,K)=CSI(2)
                ENDIF
              ENDIF
 6035   CONTINUE
      ENDIF
!
!***********************************************************
!
!  CHEM TRANSPORT 
!
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.RESTAR.EQ.'COLD START') THEN
!
!  SET:  1. INITIAL CHEM WATER COLUMN CONCENTRATIONS (ug/cm**3)
!        2. INITIAL CHEM MASS IN BED LAYERS (ug CHEM/cm**2)
!        3. INITIAL SEDIMENT MASS IN BED LAYERS (g solids/cm**2)
!
        DO 36 K=1,KBM1
          DO 36 J=2,JMM1
            DO 36 I=2,IMM1
!NG              IF (H(I,J).GT.0.0) THEN
              IF (H1(I,J).GT.0.0) THEN ! BASE IT ON H1
                IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
                  CHEM1(I,J,K)=CHEMI(1)
                  CHEM1B(I,J,K)=CHEMI(1)
                ENDIF
!
                IF (SEDTYPE.EQ.'SAND') THEN
                  CHEM2(I,J,K)=CHEMI(1)
                  CHEM2B(I,J,K)=CHEMI(1)
                ELSE
                  CHEM2(I,J,K)=CHEMI(2)
                  CHEM2B(I,J,K)=CHEMI(2)
                ENDIF
              ENDIF
 36     CONTINUE
!
        DO 37 N=1,NCHEMLAY
          DO 37 J=1,JM
            DO 37 I=1,IM
              IF (FSM(I,J).GT.0.0) THEN
                SEDMASS(N,I,J)=CBED(I,J)*CHEMTHIK
!
!  TOP LAYER IS INITIALLY VERY THIN
!
!pl 02/24/00              IF (N.EQ.1) SEDMASS(N,I,J)=0.001*SEDMASS(N,I,J)
                IF (N.EQ.1) SEDMASS(N,I,J)=0.0
                CHEMMASS(N,I,J)=SEDMASS(N,I,J)*CBEDCHEM(N,I,J)
              ELSE
                CBEDCHEM(N,I,J)=0.0
!pl add next two lines
                SEDMASS(N,I,J)=0.0
                CHEMMASS(N,I,J)=0.0
              ENDIF
 37     CONTINUE
      ENDIF
!
!***********************************************************
!
!  CALC. FETCH & MEAN DEPTH FOR WIND WAVE MODEL
!
      IF (WAVEDYN.EQ.'SMBMODEL') CALL FHCALC
!
!  INPUT INITIAL FLOW FIELD FROM gcm_tran
!
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
        OPEN (IUTRN,FORM='unformatted',FILE='gcm_tran')
        NHYDCNT=0
        CALL TRANINP(0)
      ENDIF
!
!*******************************************************************
!
!
      IF(RESTAR.EQ.'COLD START') CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
!
      IRELEND=IRELEND+INTX   !hli
!
      ISTART=INTX+1
!hli  TNDAYS=FLOAT(NSTEPS)*DAYI*DTI
      IEND=INTX+NSTEPS
!QA   SEDIMENT RUNS 1 AVERAGING PERIOD LESS THAN HYDRODYNAMICS Quamrul 10/1/98
!QA           IF (HYDTYPE.EQ.'EXTERNAL') THEN
!QA            IEND=IEND - NPLPF-1
!QA           ENDIF
      TPRT=FLOAT(IPRINT)*DAYI*DTI
      TAVG=FLOAT(IAVGE)*DAYI*DTI
      TSKILL=FLOAT(ISKILL)*DAYI*DTI
      TRACE=FLOAT(NPLPF)*DAYI*DTI
      DO 16 I=1,JHM
      IHIST(I,2)=INTX+IHIST(I,2)
 16   IHIST(I,1)=IHIST(I,2)-IAVGE+1
      DO 18 I=1,JTM
      ITRAC(I,2)=INTX+ITRAC(I,2)
 18   ITRAC(I,1)=ITRAC(I,2)-NPLPF+1
! 
      WRITE(IUPRT,15) ISTART,IEND 
 15   FORMAT(//30H MODEL STARTING UP...ISTART = ,I6,8H IEND = ,I6/)
!svv      KD=JINT(EDAY)
      KD=EDAY
      WRITE(IUPRT,*)'SDAY AND EDAY',SDAY,EDAY
      CALL CDAY(IDAE,IMOE,IYRE,KD,2)
!      CALL CDAY(IDAE,IMOE,IYRE,KD,IIDINT(2.D0))
      WRITE(IUPRT,210)IYR,IMO,IDA,IYRE,IMOE,IDAE
 210  FORMAT('SIMULATION FROM',I5,2I3,4X,'TO',3X,I5,2I3) !hli
      WRITE(IUPRT,20) TNDAYS
 20   FORMAT(//32H NUMBER OF DAYS IN SIMULATION = ,F6.2/) 
      WRITE(IUPRT,30) TPRT,IPRINT,IPRTSTART,TAVG,IAVGE,TSKILL,ISKILL
 30   FORMAT(//' TPRT =   ',F10.3,'  IPRINT =  ',I10,
     .        '  IPRTSTART =  ',I10,//
     .         ' TAVG =   ',F10.3,'  IAVGE =   ',I10,//
     .         ' TSKILL = ',F10.3,'  ISKILL =  ',I10,//)
      WRITE(IUPRT,35)
      WRITE(IUPRT,40) (IHIST(I,1),IHIST(I,2),I=1,JHM)
 35   FORMAT(//' HISTORY TAKEN AT TIMESTEPS   START ----- STOP ')
 40   FORMAT(27X,I8,3X,I8)
      WRITE(IUPRT,47)
      WRITE(IUPRT,48) (ITRAC(I,1),ITRAC(I,2),I=1,JTM)
47    FORMAT(//' QUALITY PARAMETERS FOR RCA  INTEGRATED OVER TIMESTEPS
     .      START      STOP')
48    FORMAT(57X,I8,2X,I8)
!
      WRITE(IUPRT,41)BFRIC,Z0B,NU,THETA,ALPHA,TLAG,NWAVE,BCTYPE,OPTBFRIC
 41   FORMAT(//' BFRIC            =   ',F10.4,' nondimensional'/
     .         ' Z0B              =   ',F10.4,' m'/
     .         ' NU               =   ',F10.4,' nondimensional'/
     .         ' THETA            =   ',F10.4,' nondimensional'/
     .         ' ALPHA            =   ',F10.4,' nondimensional'/
     .         ' TLAG             =   ',F10.4,
     .         ' FRICTION TIME SCALE in PCLAMP BC (Hours)'/
     .         ' WAVE MODEL ACTIVATED EVERY  =   ',I10,  ' Time Steps'/
     .         ' BOUNDARY TYPE  =   ',3x,A7/
     .         ' BOTTOM FRICTION IS  =   ',A10/)
!
      TLAG=TLAG*3600+1.0e-10    ! CONVERT TO SECONDS
      IF(TLAG.GT.0)TLAG=1./TLAG
      ALPHA=ALPHA*3600  ! CONVERT TO SECONDS
      IF(ALPHA.GT.0)ALPHA=1./ALPHA
!qa
      WRITE(IUPRT,70)
  70  FORMAT(/1X,' K',6X,'Z',10X,'ZZ',8X,'DZ',/)
      DO 90 K=1,KBM1
      WRITE(IUPRT,80) K,Z(K),ZZ(K),DZ(K)
  90  CONTINUE
      WRITE(IUPRT,80) K,Z(KB)
  80  FORMAT(I3,3F10.3)
!
      THOUR=FLOAT(INTX)*DTI/3600.
      CALL FIRST
! 
      AREA=0.0
      EMI=0.0
      APEI=0.0
      DO 280 J=1,JM
      DO 280 I=1,IM
  280 AREA=AREA+FSM(I,J)*ART(I,J)
      DO 285 K=1,KBM1
      DO 285 J=1,JM
      DO 285 I=1,IM
      TRHO=(RMEAN(I,J,K)+1.)*1000.
!NG W&D      VOL=H(I,J)*ART(I,J)*DZ(K)
      VOL=D(I,J)*ART(I,J)*DZ(K)
      EMI=EMI+TRHO*VOL*FSM(I,J)
!NG W&D  285 APEI=APEI+GRAV*TRHO*ZZ(K)*H(I,J)*VOL*FSM(I,J)
  285 APEI=APEI+GRAV*TRHO*ZZ(K)*D(I,J)*VOL*FSM(I,J)
      WRITE(IUPRT,85) AREA,EMI,APEI
   85 FORMAT(//
     . ' SURFACE AREA                          ',1PE14.7,' m**2  ',//,
     . ' INITIAL MASS                          ',1PE14.7,' Kg    ',//,
     . ' INITIAL AVAILABLE POTENTIAL ENERGY    ',1PE14.7,' joules',//)
*
*      INITIALIZE TBDRY2 AND SBDRY2   Khan 080999
*******   NEW VARIABLES FOR USE IN antidif.f   Khan 080999
*
         DO N = 1, NUMEBC
           IE = IETA(N)
           JE = JETA(N)
           DO K = 1, KBM1
!             TBDRY2(N,K) = TB(IE,JE,K)
!             SBDRY2(N,K) = SB(IE,JE,K)
             TBDRY2(N,K) = T(IE,JE,K)
             SBDRY2(N,K) = S(IE,JE,K)
           ENDDO
         ENDDO
*
*******   NEW VARIABLES FOR USE IN antidif.f   Khan 080999
*
!
!-----------------------------------------------------------------------
!
!                 BEGIN NUMERICAL INTEGRATION
!
!-----------------------------------------------------------------------
!
                   DO 9000 INTX=ISTART,IEND
! 
!           RAMP=TANH(FLOAT(INTX)/FLOAT(IRAMP+1))
        IF(INTX.GT.IRAMP) then
           RAMP=1.0
        else
!           print *, 'WHY THE HELL ARE WE HERE?'
           RAMP = FLOAT(INTX)/FLOAT(IRAMP) ! stop divide by zero bug! dbk
        endif
      TIME=FLOAT(INTX)*DAYI*DTI
      THOUR=TIME*24.
!
!  WAVE EFFECTS ON BOTTOM FRICTION (HYDRODYNAMICS ONLY)
!
      CALL BCOND(7,DTI2,0)
!
!NKIM11232010 NWAVECNT IS NOW IN BOTTAU (and comdeck)
!NKIM11232010      IF (WAVEDYN.EQ.'SMBMODEL'.OR.
!NKIM11232010     .WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'EXTERNAL') THEN
!NKIM11232010CNG02222007     .WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'EXTERNAL') THEN
!NKIM11232010        NWAVECNT=NWAVECNT+1
!NKIM11232010        IF (NWAVECNT.EQ.NWAVE) THEN
!NKIM11232010          CALL BCOND(11,DTI2,0)
!NKIM11232010          NWAVECNT=0
      IF(OPTBFRIC.EQ.'MANNINGS-N')THEN
       CALL UPDATE_CBCM ! CNG2014 Manning's N
      ELSE
       CALL UPDATE_CBC ! NG 2014; move it here, from BOTTAU, and just call it before BOTTAU
      ENDIF
      CALL BOTTAU   ! compute shear stress using currents+wave
!NKIM11232010        ENDIF
!NKIM11232010      ELSE
!NKIM11232010          CALL BOTTAU   ! compute shear stress using currents only
!NKIM11232010      ENDIF
!
!NG04152008
      IF (WAVEDYN.EQ.'DONONLY') THEN
       CALL ARCHIVE
       GOTO 9910
      ENDIF

      IF (HYDTYPE.EQ.'EXTERNAL') THEN
!
!     Like Hydrodynamics skip first time step QA 9/25/98
!
!
        NHYDCNT=NHYDCNT+1
        IF(INTX.EQ.1) GO TO 8200  !hli
!  
!  INPUT NEW FLOW FIELD
!
        IF (NHYDCNT.GT.NHYD) THEN
!qa          write(*,*)'INTX,NHYDCNT,NHYD=',INTX,NHYDCNT,NHYD
          CALL TRANINP(1)
          NHYDCNT=1
        ENDIF
!
!qa        CALL BCOND(7,DTI2,0)
!
!  SKIP HYDRODYNAMICS
!
        GOTO 9900
      ENDIF
!
!
      IF(TOR.NE.'BAROTROPIC') THEN
!
      CALL BAROPG(DRHOX,DRHOY,TRNU,TRNV,SGW)
!
      DO 50 J=1,JM
      DO 50 I=1,IM
      TRNU(I,J)=TRNU(I,J)+ADVUU(I,J)-ADVUA(I,J)
 50   TRNV(I,J)=TRNV(I,J)+ADVVV(I,J)-ADVVA(I,J)
!
      DO 120 J=1,JM
      DO 120 I=1,IM
 120  EGF(I,J)=EL(I,J)*ISPI
      DO 400 J=2,JM
      DO 400 I=2,IM
      UTF(I,J)=UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
      VTF(I,J)=VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
400   CONTINUE
!
      ENDIF

      IF(HORZMIX.EQ.'CLOSURE   ') CALL SMAG !CNG MOVE UP HERE, BEFORE DUMWAD,DVMWAD IS RESET FOR NEXT EXTERNAL

      DO 130 J=1,JM
      DO 130 I=1,IM
      DUMWAD(I,J)=DUM(I,J) ! NG09292011
      DVMWAD(I,J)=DVM(I,J) ! NG09292011
130   CONTINUE

!
!NG      IF(HORZMIX.EQ.'CLOSURE   ') CALL SMAG
!
!-----------------------------------------------------------------------
!         BEGIN EXTERNAL MODE
!-----------------------------------------------------------------------
!
      IF(TOR.EQ.'BAROTROPIC'.AND.JTM.EQ.0.AND.ISPLIT.GT.0) THEN
        CALL BARO_SUBCYCLE(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,
     &                     VAMAX,ELMINX,IAMAX,JAMAX)
        IF(VAMAX.GT.100.) GO TO 9100
        GO TO 8001
      ENDIF
                   DO 8000 IEXT=1,ISPLIT
!
!pl      CALL EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,BFRIC,DTI2)
      CALL EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2)
      IF(TOR.EQ.'BAROTROPIC') GO TO 440
      IF(IEXT.LT.(ISPLIT-2)) GO TO 440
      IF(IEXT.EQ.(ISPLIT-2)) THEN
      DO 402 J=1,JM
      DO 402 I=1,IM
 402  ETF(I,J)=.25*NU*ELF(I,J)
      GO TO 440
      ENDIF
      IF(IEXT.EQ.(ISPLIT-1)) THEN
      DO 404 J=1,JM
      DO 404 I=1,IM
 404  ETF(I,J)=ETF(I,J)+.5*(1.-.5*NU)*ELF(I,J)
      GO TO 440
      ENDIF
      IF(IEXT.EQ.(ISPLIT-0)) THEN
      DO 406 J=1,JM
      DO 406 I=1,IM
 406  ETF(I,J)=(ETF(I,J)+.5*ELF(I,J))*FSM(I,J)
      ENDIF
 440  CONTINUE
!-----------------------------------------------------------------------
!         APPLY FILTER TO REMOVE TIME SPLIT AND RESET TIME SEQUENCE
!-----------------------------------------------------------------------
!
!     Without an intervening transport call these operations are all
!     pointwise. Read the filtered value at its original precision and
!     rotate immediately. Retain the original path when tracers exist.
      IF (TOR.NE.'BAROTROPIC'.OR.JTM.EQ.0) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
        DO J=1,JM
          DO I=1,IM
            UAB(I,J)=UA(I,J)+.5*NU*(UAB(I,J)-2.*UA(I,J)+UAF(I,J))
            VAB(I,J)=VA(I,J)+.5*NU*(VAB(I,J)-2.*VA(I,J)+VAF(I,J))
            ELB(I,J)=EL(I,J)+.5*NU*(ELB(I,J)-2.*EL(I,J)+ELF(I,J))
            EL(I,J)=ELF(I,J)
            D(I,J)=H(I,J)+EL(I,J)
            UA(I,J)=UAF(I,J)
            VA(I,J)=VAF(I,J)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ELSE
      DO 150 J=1,JM 
      DO 150 I=1,IM 
      UA(I,J)=UA(I,J)+.5*NU*(UAB(I,J)-2.*UA(I,J)+UAF(I,J)) 
      VA(I,J)=VA(I,J)+.5*NU*(VAB(I,J)-2.*VA(I,J)+VAF(I,J)) 
 150  EL(I,J)=EL(I,J)+.5*NU*(ELB(I,J)-2.*EL(I,J)+ELF(I,J)) 
      IF(TOR.EQ.'BAROTROPIC') THEN
      DO 155 JTRAC=1,JTM
      IF(INTX.GE.ITRAC(JTRAC,1).AND.INTX.LE.ITRAC(JTRAC,2))
     +   CALL TRANSPORT 
 155  CONTINUE
      END IF
      DO 160 J=1,JM 
      DO 160 I=1,IM 
!NG      DUMWAD(I,J)=DUM(I,J) ! NG05192010
!NG      DVMWAD(I,J)=DVM(I,J) ! NG05192010
      ELB(I,J)=EL(I,J)
      EL(I,J)=ELF(I,J)
      D(I,J)=H(I,J)+EL(I,J) 
      UAB(I,J)=UA(I,J)
      UA(I,J)=UAF(I,J)
      VAB(I,J)=VA(I,J)
 160  VA(I,J)=VAF(I,J)
      ENDIF
!
!
!     INTRODUCE WET/DRY SCHEME BASED ON KKY's DISSERTATION (hli, 08/03/05).
!
      CALL WANDD_KKY(EL,D)
!
!  TEST FOR MODEL BLOWUP. IF SO, PRINT AND STOP
!
     
!$OMP PARALLEL DEFAULT(SHARED)
      CALL SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
!$OMP END PARALLEL
      IF(VAMAX.GT.100.) GO TO 9100
!end,   check elevation minimum, jeff ji, 1/17/96

      IF(TOR.EQ.'BAROTROPIC') GO TO 8000
      IF(IEXT.EQ.ISPLIT) GO TO 8000
      DO 445 J=1,JM
      DO 445 I=1,IM
 445  EGF(I,J)=EGF(I,J)+EL(I,J)*ISPI
      DO 450 J=2,JM
      DO 450 I=2,IM
      UTF(I,J)=UTF(I,J)+UA(I,J)*(D(I,J)+D(I-1,J))*ISP2I
 450  VTF(I,J)=VTF(I,J)+VA(I,J)*(D(I,J)+D(I,J-1))*ISP2I
!
 8000              CONTINUE
 8001 CONTINUE
!
!-----------------------------------------------------------------------
!         END EXTERNAL (2-D) MODE CALCULATION
!      AND CONTINUE WITH INTERNAL (3-D) MODE CALCULATION
!-----------------------------------------------------------------------
! 
      IF(INTX.EQ.1) GO TO 8200    !hli
!qa       IF(INTX.EQ.1)THEN
!qa        DO 305 JTRAC=1,JTM
!qa        IF(INTX.GE.ITRAC(JTRAC,1).AND.INTX.LE.ITRAC(JTRAC,2)) 
!qa     .    CALL TRANSPORT 
!qa 305    CONTINUE
!qa       GO TO 8200 
!qa       ENDIF
!
      IF (TOR.EQ.'BAROTROPIC') THEN
        CALL ARCHIVE
        GO TO 9900
      END IF
!
!-----------------------------------------------------------------------
!         ADJUST U(Z) AND V(Z) SUCH THAT
!      VERTICAL AVERAGE OF (U,V) = (UA,VA)
!-----------------------------------------------------------------------
!
      DO 299 J=1,JM
      DO 299 I=1,IM
 299  TPS(I,J)=0.0
      DO 300 K=1,KBM1
      DO 300 J=2,JM
      DO 300 I=2,IM
 300  TPS(I,J)=TPS(I,J)+U(I,J,K)*DZ(K)
      DO 302 K=1,KBM1
      DO 302 J=2,JM
      DO 302 I=2,IM
 302  U(I,J,K)=(U(I,J,K)-TPS(I,J))+
     .     (UTB(I,J)+UTF(I,J))/(DT(I,J)+DT(I-1,J))
      DO 303 J=1,JM
      DO 303 I=1,IM
 303  TPS(I,J)=0.0
      DO 304 K=1,KBM1
      DO 304 J=2,JM
      DO 304 I=2,IM
 304  TPS(I,J)=TPS(I,J)+V(I,J,K)*DZ(K)
      DO 306 K=1,KBM1
      DO 306 J=2,JM
      DO 306 I=2,IM
 306  V(I,J,K)=(V(I,J,K)-TPS(I,J))+
     .     (VTB(I,J)+VTF(I,J))/(DT(I,J)+DT(I,J-1))
!
!-------- CORRECT FOR RIVER BOUNDARY CONDITIONS ------------------------
      DO 320 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 320 K=1,KBM1
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            IF(VQDIST(N,K).EQ.0.0) U(IC,JC,K)=0.0
         ELSE
            IF(VQDIST(N,K).EQ.0.0) U(ID,JD,K)=0.0
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            IF(VQDIST(N,K).EQ.0.0) V(IC,JC,K)=0.0
         ELSE
            IF(VQDIST(N,K).EQ.0.0) V(ID,JD,K)=0.0
         ENDIF
      ENDIF
  320 CONTINUE
!
!-----------------------------------------------------------------------
!         CALCULATE HORIZONTAL MASS FLUXES, (H2*U*D) AND (H1*V*D)
!-----------------------------------------------------------------------
      DO 310 K=1,KBM1
      DO 311 J=2,JMM1
      DO 311 I=2,IM
      XMFL3D(I,J,K)=0.25*(H2(I-1,J)+H2(I,J))*(DT(I-1,J)+DT(I,J))*
     .   U(I,J,K)
  311 CONTINUE
      DO 312 J=2,JM
      DO 312 I=2,IMM1
      YMFL3D(I,J,K)=0.25*(H1(I,J-1)+H1(I,J))*(DT(I,J-1)+DT(I,J))*
     .   V(I,J,K)
  312 CONTINUE
  310 CONTINUE
!
!-----------------------------------------------------------------------
!         VERTVL INPUT = U,V,DT(=H+ET),ETF,ETB; OUTPUT = W
!----------------------------------------------------------------------- 
      CALL VERTVL(DTI2)
      CALL BCOND(5,DTI2,0)
      DO 307 JTRAC=1,JTM
      IF(INTX.GE.ITRAC(JTRAC,1).AND.INTX.LE.ITRAC(JTRAC,2))
     +   CALL TRANSPORT 
 307  CONTINUE
**
**      MOVED TO BE CONSITENT WITH  ecom3d  Khan 080999
        CALL ARCHIVE
**      MOVED TO BE CONSITENT WITH  ecom3d  Khan 080999
**
!
!-----------------------------------------------------------------------
!         COMPUTE Q2F AND Q2LF USING UF AND VF AS TEMPORARY VARIABLES
!-----------------------------------------------------------------------
!
      IF(VERTMIX.EQ.'CLOSURE   ') THEN
      CALL ADVQ(Q2B,Q2,DTI2,UF)
      CALL ADVQ(Q2LB,Q2L,DTI2,VF)
      CALL PROFQ(DTI2)
      CALL BCOND(6,DTI2,0)
      DO 325 K=1,KB 
      DO 325 J=1,JM 
      DO 325 I=1,IM 
      Q2 (I,J,K)=Q2 (I,J,K)+.5*NU*(UF(I,J,K)+Q2B(I,J,K)-2.*Q2(I,J,K)) 
      Q2B(I,J,K)=Q2(I,J,K)
      Q2(I,J,K)=UF(I,J,K) 
 325  CONTINUE
      DO 335 K=1,KB 
      DO 335 J=1,JM 
      DO 335 I=1,IM 
      Q2L(I,J,K)=Q2L(I,J,K)+.5*NU*(VF(I,J,K)+Q2LB(I,J,K)-
     .  2.*Q2L(I,J,K))
      Q2LB(I,J,K)=Q2L(I,J,K)
      Q2L(I,J,K)=VF(I,J,K)
 335  CONTINUE
      END IF
      IF(VERTMIX.EQ.'EMPIRICAL ')CALL PROFE
!
!-----------------------------------------------------------------------
!     COMPUTE UPDATED BOTTOM FRICTION (CBC) EVERY TIMESTEP
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!         COMPUTE TF AND SF USING UF AND VF AS TEMPORARY VARIABLES
!-----------------------------------------------------------------------
!
*
*******    CHANGES FOR antidif.f BC        Khan 080999
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY ')
!C     .   CALL ADVT(TB,T,TMEAN,DTI2,UF,TDIF,TDIS,TBDRY) !ORIGINAL
     .    CALL ADVT(TB,T,TMEAN,DTI2,UF,TDIF,TDIS,TBDRY2)
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY ')
!C     .   CALL ADVT(SB,S,SMEAN,DTI2,VF,SDIF,SDIS,SBDRY) !ORIGINAL
     .    CALL ADVT(SB,S,SMEAN,DTI2,VF,SDIF,SDIS,SBDRY2) 
*******    CHANGES FOR antidif.f BC        Khan 080999
*
!qa
!qa   Modified by Quamrul QA for SWRAD penetration on 3/10/99
!qa
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY ')
     .    CALL PROFT(UF,WTSURF,DTI2,SWRAD)
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY ')
     .    CALL PROFT(VF,WSSURF,DTI2,ZEROS)
      CALL BCOND(4,DTI2,0)
!qa
*******   NEW VARIABLES FOR USE IN antidif.f   Khan 080999
*
         DO N = 1, NUMEBC
           IE = IETA(N)
           JE = JETA(N)
           DO K = 1, KBM1
             TBDRY2(N,K) = UF(IE,JE,K)
             SBDRY2(N,K) = VF(IE,JE,K)
           ENDDO
         ENDDO 
*
*******   NEW VARIABLES FOR USE IN antidif.f   Khan 080999

      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'TEMP_ONLY ')THEN
      DO 345 K=1,KB 
      DO 345 J=1,JM 
      DO 345 I=1,IM 
      T(I,J,K)=T(I,J,K)+.5*NU*(UF(I,J,K)+TB(I,J,K)-2.*T(I,J,K))
      TB(I,J,K)=T(I,J,K)
      T(I,J,K)=UF(I,J,K)
 345  CONTINUE
      ENDIF
!qa
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY ')THEN
      DO 355 K=1,KB 
      DO 355 J=1,JM 
      DO 355 I=1,IM 
      S(I,J,K)=S(I,J,K)+.5*NU*(VF(I,J,K)+SB(I,J,K)-2.*S(I,J,K))
      SB(I,J,K)=S(I,J,K)
      S(I,J,K)=VF(I,J,K)
 355  CONTINUE
!PL      WRITE(4001,*)'INSIDE ECOM',INTX
!PL      WRITE(4001,'(8F10.3)')(S(IETA(N),JETA(N),1),N=1,NUMEBC)
!PL      WRITE(4001,'(8F10.3)')(VF(IETA(N),JETA(N),1),N=1,NUMEBC)

      ENDIF
! 
      IF(TOR.EQ.'PROGNOSTIC'.OR.TOR.EQ.'SALT_ONLY '.OR.
     .   TOR.EQ.'TEMP_ONLY ')THEN
       CALL DENS
      END IF
!
!  FOR CONSERVATIVE TRACER
!
 9900 IF (TRACER.EQ.'INCLUDE') THEN

!qa   We Introduce a decay coefficient which is function of T,S and Io
!qa   Based on Mancini 1978.
!qa   K = kd + Ki + Ks
!qa   Kd = decay rate as a function of T,S
!qa      = (0.8 + 0.006*(% of sea water))*1.07**(t-20)
!qa   Ki = decay rate due to irradiance Io
!qa      = alpha*Io/(Ke*Z) [ 1-exp (-Ke*Z)]
!qa   Ks = Vs/Z
!qa
!qa   We define decay rate DKRATE(I,J,K) = DKD(I,J,K) + DKI(I,J,K) + DKS(I,J,K)
!qa
      DO 141 K=1,KBM1
      DO 141 J=2,JMM1
      DO 141 I=2,IMM1
        IF(TROPT.EQ.'PATHOGEN  ')THEN
         DKD(I,J,K) = (CONDRAT + 0.006*ASAL*(S(I,J,K)/35.0*100.))*
     .                                   THETAT**(T(I,J,K)-20.)
         SOLRAD = SWRAD(I,J)*(-4.186E6) * 2.066115 /24.0 ! Langley/Hr
         DKI(I,J,K) = (ASOL*SOLRAD/(EXTC*(-ZZ(K))*DT(I,J)))*
     .                       (1.0-EXP(EXTC*ZZ(K)*DT(I,J)))

         DKS(I,J,K) = VSRATE/DT(I,J)
         DKRATE(I,J,K) = DKD(I,J,K) + DKI(I,J,K) + DKS(I,J,K)
         DKRATE(I,J,K) = DKRATE(I,J,K)/86400.  ! rate to 1/sec
        ELSE IF(TROPT.EQ.'SIMPDECAY ')THEN
         DKRATE(I,J,K)=CONDRAT/86400.
        END IF
141   CONTINUE
!

        CALL ADVCON(CONC1B,CONC1,CMEAN1,DTI2,UF,CDIF1,CDIS1,CBDRY1,
     +                               CPSTR)
!
        CALL PROFT(UF,WCSURF,DTI2,ZEROS)
        CALL BCOND(8,DTI2,0)
        DO 346 K=1,KB
          DO 346 J=1,JM 
            DO 346 I=1,IM 
              CONC1(I,J,K)=CONC1(I,J,K)+.5*NU*(UF(I,J,K)+CONC1B(I,J,K)-
     +              2.*CONC1(I,J,K))
              CONC1B(I,J,K)=CONC1(I,J,K)
              CONC1(I,J,K)=UF(I,J,K)
 346    CONTINUE
      ENDIF
!
!********************************************************************
!
!  FOR PARTICLE TRACKING 
!
        IF (PARTICLE.EQ.'INCLUDE') THEN
          IF (INTX.GE.IRELST) THEN
            IF (INTX.EQ.IRELST) IDUM=-1
            CALL PARTRAK(IDUM,IRELEND)
          ENDIF
        ENDIF
!
!********************************************************************
!
!  SEDIMENT TRANSPORT
!
        IF (SEDTRAN.EQ.'INCLUDE'.AND.INTX.GT.NSEDBEG) THEN
!
!qa       WHEN INTX=1 HYDRODYNAMIC SKIPS ie THIS SEGMENT SKIPS    
          IF(INTX.EQ.2)NSEDCT=1
          IF(INTX.EQ.2)N24CNT=1
          NSEDCT=NSEDCT+1
          N24CNT=N24CNT+1

!  COHESIVE CLASS
!
          IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
            DO 6300 N=1,NUMQBCSE
              DUMB1(N)=CDIS(1,N)
 6300       CONTINUE
!
            DO 6310 K=1,KBM1
              DO 6310 N=1,NUMEBCSE
                DUMB2(N,K)=CBDRY(1,N,K)
 6310      CONTINUE
!
!  MODIFIED BY C.K.Z. ON 5/12/97
!
!  FOR DIFFUSER INPUTS
!
            DO 6315 N=1,NUMDBCSE
                DUMB3(N)=CSDIFF(1,N)
 6315       CONTINUE
!
! yang,08/31/99
!           CALL ADVSED(CSED1B,CSED1,DTI2,UF,DUMB1,DUMB2,DUMB3,WSET1,
!     +                                          WCT1BOT)
            CALL ADVSED(CSED1B,CSED1,DTI2,UF,DUMB1,DUMB2,DUMB3)
!
!           CALL PROFSED(UF,DTI2)
            CALL PROFSED(UF,DTI2,WSET1,WCT1BOT)  !! yang, 08/31/99
!
            CALL BCOND(9,DTI2,1)
!
            DO 6320 K=1,KBM1 
              DO 6320 J=1,JM
                DO 6320 I=1,IM
! yang start (dec/22/1999): set minimum of conc to prevent from creat mass.
!                 CSED1(I,J,K)=CSED1(I,J,K)+.5*NU*(UF(I,J,K)+
!    +                         CSED1B(I,J,K)-2.*CSED1(I,J,K))
 
                  IF (UF(I,J,K).LE. 0.0 )THEN
                    WCT1BOT(I,J)=0.0
                    CSED1(I,J,K)=0.0
                  ELSE
                    CSED1(I,J,K)=CSED1(I,J,K)+.5*NU*(UF(I,J,K)+
     +                           CSED1B(I,J,K)-2.*CSED1(I,J,K))
                  ENDIF
! yang end
 
                  CSED1B(I,J,K)=CSED1(I,J,K)
                  CSED1(I,J,K)=UF(I,J,K)
 6320       CONTINUE
          ENDIF
!
8900  Format (f8.3,I7,10F10.5)
9007  Format (2X,'        (02,12)          (03,11)
     .          (03,12)          (03,13)')


!
!  NON-COHESIVE CLASS
!
          IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
!
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
!
            DO 6330 N=1,NUMQBCSE
              DUMB1(N)=CDIS(KK,N)
 6330       CONTINUE
!
            DO 6340 K=1,KBM1
              DO 6340 N=1,NUMEBCSE
                DUMB2(N,K)=CBDRY(KK,N,K)
 6340       CONTINUE
!
!  MODIFIED BY C.K.Z. ON 5/12/97
!
!  FOR DIFFUSER INPUTS
!
            DO 6345 N=1,NUMDBCSE
                DUMB3(N)=CSDIFF(KK,N)
 6345       CONTINUE
!
! yang. 08/31/99
!          CALL ADVSED(CSED2B,CSED2,DTI2,UF,DUMB1,DUMB2,DUMB3,
!    +                                             WSET2,WCT2BOT)
            CALL ADVSED(CSED2B,CSED2,DTI2,UF,DUMB1,DUMB2,DUMB3)
!
!
!           CALL PROFSED(UF,DTI2)
            CALL PROFSED(UF,DTI2,WSET2,WCT2BOT)   !! yang, 08/31/99
!
            CALL BCOND(9,DTI2,KK)
!
            DO 6350 K=1,KBM1
              DO 6350 J=1,JM
                DO 6350 I=1,IM
! yang start (dec/22/1999): set minimum of conc to prevent from creat mass.
!                 CSED2(I,J,K)=CSED2(I,J,K)+.5*NU*(UF(I,J,K)+
!    +                         CSED2B(I,J,K)-2.*CSED2(I,J,K))
 
                  IF (UF(I,J,K).LE. 0.0 )THEN
                    WCT2BOT(I,J)=0.0
                    CSED2(I,J,K)=0.0
                  ELSE
                    CSED2(I,J,K)=CSED2(I,J,K)+.5*NU*(UF(I,J,K)+
     +                           CSED2B(I,J,K)-2.*CSED2(I,J,K))
                  ENDIF
! yang end
 
                  CSED2B(I,J,K)=CSED2(I,J,K)
                  CSED2(I,J,K)=UF(I,J,K)

 6350       CONTINUE
          ENDIF
!
!  CHECK IF MODEL IS BLOWING UP
!
!  STOP PROGRAM IF CONC. > 100,000 mg/l (0.10 g/cm**3)
!
          NBLOW=NBLOW+1
          IF (NBLOW.EQ.100) THEN
            NBLOW=0
            CMAX=0.10
            DO 6360 J=3,JM-2
              DO 6360 I=2,IM-1
                IF (H(I,J).GT.0.01) THEN
                  TOTCON=CSED1(I,J,1)+CSED2(I,J,1)
                  IF (ABS(TOTCON).GT.CMAX) THEN
                  WRITE (IUPRT,6370)INTX,I,J,CSED1(I,J,1),CSED2(I,J,1)
 6370             FORMAT (/5X,'******* PROGRAM EXECUTION STOPPED AT',
     +   I8,' TIMESTEP',/8X,'CONCENTRATION > 100,000 mg/l AT ',2I3,/8X,
     +   'C1 AND C2 =',2E12.4)
                  WRITE (6,6380)TIME
 6380             FORMAT (/5X,'$$$$ TIME (days) =',F10.4)
!
                  CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
!
                  STOP
                  ENDIF
                ENDIF
 6360       CONTINUE
          ENDIF
!
!***********************************************************
!
!  CHEM TRANSPORT 
!
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
!
!  COHESIVE CLASS
!
            IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
              DO 401 N=1,NUMQBCCH
!pl 03/03/00                DUMB1(N)=CHEMDIS(1,N)
                DUMB1(N)=PDIS(1,N)
 401          CONTINUE
!
              DO 410 K=1,KBM1
                DO 410 N=1,NUMEBCCH
!pl 03/03/00                  DUMB2(N,K)=CHEMBDRY(1,N,K)
                  DUMB2(N,K)=PBDRY(1,N,K)
 410          CONTINUE
!
!  FOR DIFFUSER INPUTS   Quamrul QA 5/25/99
!
            DO 415 N=1,NUMDBCCH
                DUMB3(N)=PDIFF(1,N)
 415        CONTINUE
!
!             CALL ADVCHEM(CHEM1B,CHEM1,DTI2,UF,DUMB1,DUMB2,DUMB3,WSET1,
!    +                     CHEMBOT1,CHEMDRAT1)
              CALL ADVCHEM(CHEM1B,CHEM1,DTI2,UF,DUMB1,DUMB2,DUMB3,
     +                     CHEMDRAT1)             !! yang, 11/29/99
!
!             CALL PROFSED(UF,DTI2)
              CALL PROFSED(UF,DTI2,WSET1,CHEMBOT1)   !! yang, 08/31/99
!
              CALL BCOND(10,DTI2,1)
!
              DO 420 K=1,KBM1
                DO 420 J=1,JM
                  DO 420 I=1,IM
! yang start (dec/22/1999): set minimum of conc to prevent from creat mass.
!                   CHEM1(I,J,K)=CHEM1(I,J,K)+.5*NU*(UF(I,J,K)+
!    +                CHEM1B(I,J,K)-2.*CHEM1(I,J,K))
 
                    IF (UF(I,J,K).LE. 0.0 )THEN
                      CHEMBOT1(I,J)=0.0
                      CHEM1(I,J,K)=0.0
                    ELSE
                      CHEM1(I,J,K)=CHEM1(I,J,K)+.5*NU*(UF(I,J,K)+
     +                  CHEM1B(I,J,K)-2.*CHEM1(I,J,K))
                    ENDIF
! yang end
                    CHEM1B(I,J,K)=CHEM1(I,J,K)
                    CHEM1(I,J,K)=UF(I,J,K)
 420          CONTINUE
            ENDIF
!
!  NON-COHESIVE CLASS
!
            IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
              IF (SEDTYPE.EQ.'SAND') THEN
                KK=1
              ELSE
                KK=2
              ENDIF
!
              DO 430 N=1,NUMQBCCH
!pl 03/03/00       DUMB1(N)=CHEMDIS(KK,N)
                DUMB1(N)=PDIS(KK,N)
 430          CONTINUE
!
              DO 441 K=1,KBM1
                DO 441 N=1,NUMEBCCH
!pl 03/03/00                  DUMB2(N,K)=CHEMBDRY(KK,N,K)
                  DUMB2(N,K)=PBDRY(KK,N,K)
 441          CONTINUE
!
!  FOR DIFFUSER INPUTS   Quamrul QA 5/25/99
!
            DO 439 N=1,NUMDBCCH
                DUMB3(N)=PDIFF(KK,N)
 439        CONTINUE
!
!            CALL ADVCHEM(CHEM2B,CHEM2,DTI2,UF,DUMB1,DUMB2,DUMB3,WSET2,
!    +                                          CHEMBOT2,CHEMDRAT2)
              CALL ADVCHEM(CHEM2B,CHEM2,DTI2,UF,DUMB1,DUMB2,DUMB3,
     +                                      CHEMDRAT2)   !! yang, 11/29/99
!
!             CALL PROFSED(UF,DTI2)
              CALL PROFSED(UF,DTI2,WSET2,CHEMBOT2)   !! yang, 08/31/99
!
              CALL BCOND(10,DTI2,KK)
!
              DO 451 K=1,KBM1
                DO 451 J=1,JM
                  DO 451 I=1,IM
! yang start (dec/22/1999): set minimum of conc to prevent from creat mass.
!                   CHEM2(I,J,K)=CHEM2(I,J,K)+.5*NU*(UF(I,J,K)+
!    +                CHEM2B(I,J,K)-2.*CHEM2(I,J,K))
 
                    IF (UF(I,J,K).LE. 0.0 )THEN
                      CHEMBOT2(I,J)=0.0
                      CHEM2(I,J,K)=0.0
                    ELSE
                      CHEM2(I,J,K)=CHEM2(I,J,K)+.5*NU*(UF(I,J,K)+
     +                  CHEM2B(I,J,K)-2.*CHEM2(I,J,K))
                    ENDIF
! yang end
                    CHEM2B(I,J,K)=CHEM2(I,J,K)
                    CHEM2(I,J,K)=UF(I,J,K)
 451          CONTINUE
            ENDIF
          ENDIF
!
!  CALCULATE FLUXES AT SEDIMENT-WATER INTERFACE
!
!qa          if(TIME.GE.18..AND.TIME.LE.20.)
!qa     .write(2001,*)'INTX,NSEDCT,NSBED=',INTX,NSEDCT,NSBED
        IF (NSEDCT.EQ.NSBED) THEN
!
!  LIMIT DEPOSITION TO BATHYMETRIC DEPTH 
!  PRAVI 08/12/99
!
         DO 459 J=1,JM
         DO 459 I=1,IM
          COHTHK(I,J) =0.0
          NCOHTHK(I,J)=0.0
          IF (IBMSK(I,J).EQ.0) THEN
           DO 458 LL=1,LAYMAX
            COHTHK(I,J)=COHTHK(I,J)+TSED(LL,I,J)
 458       CONTINUE
           COHTHK(I,J)=COHTHK(I,J)/CBED(I,J)
          ELSE
           IF (IBMSK(I,J).EQ.1) THEN
            NCOHTHK(I,J)=NCOHTHK(I,J)+BEDTH(1,I,J)
            NCOHTHK(I,J)=100.*(NCOHTHK(I,J)/((CBED(I,J)/2.65)*
     +                   H1(I,J)*H2(I,J)))
           ENDIF
          ENDIF
          SEDTHK(I,J)=(COHTHK(I,J)+NCOHTHK(I,J))/100.
 459     CONTINUE
!
!  CALC. BOTTOM SHEAR STRESSES
!
!qa          CALL BOTTAU
!
!  COHESIVE BED MODEL
!
          IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') CALL SEDFLX
!
!  NON-COHESIVE BED MODEL
!
          IF (IBED.EQ.1) CALL SUSLOD
!
!  CALCULATE CHEM FLUX AT SEDIMENT-WATER INTERFACE
!
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
            CALL CHEMFLX
!
!  CHEM BED MODEL
!
            CALL CHEMBED(CHEMDRAT1,CHEMDRAT2)
          ENDIF
!
          NSEDCT=0
        ENDIF
!
!***********************************************************
!
!  REORDER COHESIVE SEDIMENT BED LAYERS, IF NECESSARY
!
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (N24CNT.EQ.N24HR) CALL REORDR
        ENDIF
      ENDIF
!
!  MODIFY ETA AND VELOCITIES FOR NEXT TIMESTEP (HYDRO INFO READ FROM hqi_tran)
!
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
        DO 9920 J=2,JMM1
          DO 9920 I=2,IMM1
            ETB(I,J)=ET(I,J)
            ET(I,J)=ETF(I,J)
!
            ETF(I,J)=ETB(I,J)+2.*DTI*DETA(I,J)
!
            DT(I,J)=H(I,J)+ET(I,J)
 9920   CONTINUE
!
!  CALC NEW VELOCITIES
!
        DO 9930 K=1,KBM1
          DO 9930 J=2,JMM1
            DO 9930 I=2,IM
              U(I,J,K)=0.0
              DBAR=0.25*(DT(I,J)+DT(I-1,J))*(H2(I-1,J)+H2(I,J))
            IF (DBAR.GT.0.0) U(I,J,K)=XMFL3D(I,J,K)/DBAR
 9930   CONTINUE
!
        DO 9940 K=1,KBM1
          DO 9940 J=2,JM
            DO 9940 I=2,IMM1
              V(I,J,K)=0.0
              DBAR=0.25*(DT(I,J)+DT(I,J-1))*(H1(I,J-1)+H1(I,J))
              IF (DBAR.GT.0.0) V(I,J,K)=YMFL3D(I,J,K)/DBAR
 9940   CONTINUE
!
!  SKIP HYDRODYNAMICS IF READING FROM gcm_tran
!
        GOTO 9910
      ENDIF
!
      IF (TOR.EQ.'BAROTROPIC') GO TO 9910
!
!********************************************************************
!
!-----------------------------------------------------------------------
!         COMPUTE UF AND VF
!-----------------------------------------------------------------------
!
      CALL ADVU(DRHOX,ADVUU,DTI2)
      CALL ADVV(DRHOY,ADVVV,DTI2)
      CALL PROFU(DTI2)
      CALL PROFV(DTI2)
      CALL BCOND(3,DTI2,0)
!
!JKL
      IF(BCTYPE.EQ.'MIXED  ') CALL MIXED(2)
!JKL
!
      DO 369 J=1,JM 
      DO 369 I=1,IM 
 369  TPS(I,J)=0.0
      DO 370 K=1,KBM1
      DO 370 J=1,JM 
      DO 370 I=1,IM 
 370  TPS(I,J)=TPS(I,J)+(UF(I,J,K)+UB(I,J,K)-2.*U(I,J,K))*DZ(K)
      DO 372 K=1,KBM1
      DO 372 J=1,JM 
      DO 372 I=1,IM 
 372  U(I,J,K)=U(I,J,K)+.5*NU*(UF(I,J,K)+UB(I,J,K)-2.*U(I,J,K)
     .        -TPS(I,J))
      DO 373 J=1,JM 
      DO 373 I=1,IM 
 373  TPS(I,J)=0.0
      DO 374 K=1,KBM1
      DO 374 J=1,JM 
      DO 374 I=1,IM 
 374  TPS(I,J)=TPS(I,J)+(VF(I,J,K)+VB(I,J,K)-2.*V(I,J,K))*DZ(K)
      DO 376 K=1,KBM1
      DO 376 J=1,JM 
      DO 376 I=1,IM 
 376  V(I,J,K)=V(I,J,K)+.5*NU*(VF(I,J,K)+VB(I,J,K)-2.*V(I,J,K)
     .        -TPS(I,J))
      DO 377 K=1,KB 
      DO 377 J=1,JM 
      DO 377 I=1,IM 
      UB(I,J,K)=U(I,J,K)
      U(I,J,K)=UF(I,J,K)
      VB(I,J,K)=V(I,J,K)  
 377  V(I,J,K)=VF(I,J,K)
!NG      write (*,*) time
!NG      write (*,*) U(143,6,1),U(144,6,1),U(145,6,1)
!NG      write (*,*) U(143,5,1),U(144,5,1),U(145,5,1)
!NG      write (*,*) U(143,4,1),U(144,4,1),U(145,4,1)
!NG      write (*,*) V(143,6,1),V(144,6,1),V(145,6,1)
!NG      write (*,*) V(143,5,1),V(144,5,1),V(145,5,1)
!NG      write (*,*) V(143,4,1),V(144,4,1),V(145,4,1)
!
!     INTRODUCE WET/DRY SCHEME BASED ON KKY's DIISERTATION (hli,08/03/05).
!
!wad
      DO I=1,IMM1
         DO J=2,JMM1
!NG02222011AFTERKIMUTF?            IF(UA(I+1,J).EQ.0.0) THEN
            IF(UA(I+1,J).EQ.0.0) THEN ! ONLY WHEN AVERAGE OF EXTERNAL IS ZERO; FAT CHANCE THAT AVERAGE COMES TO EXACTLY ZERO.
               DO K=1,KBM1
                  U(I+1,J,K)=0.0
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!    
      DO I=2,IMM1
         DO J=1,JMM1
!NG02222011AFTERKIMVTF?            IF(VA(I,J+1).EQ.0.0) THEN
            IF(VA(I,J+1).EQ.0.0) THEN ! ONLY WHEN AVERAGE OF EXTERNAL IS ZERO; FAT CHANCE THAT AVERAGE COMES TO EXACTLY ZERO.
               DO K=1,KBM1
                  V(I,J+1,K)=0.0
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
      CALL WREAL(DTI2)
!
 8200 CONTINUE
!

      DO 380 J=1,JM 
      DO 380 I=1,IM 
      EGB(I,J)=EGF(I,J)
      ETB(I,J)=ET(I,J)
      ET(I,J)=ETF(I,J)
      DT(I,J)=H(I,J)+ET(I,J)
      UTB(I,J)=UTF(I,J)
 380  VTB(I,J)=VTF(I,J)
!NG04232008 INITIALIZE QEPCNUDF FOR AVERAGING OF ELEVATION POTENTIAL Qs
! FOR EPC POTENTIAL
      DO N=1,NUMEPC
       QEPCNUDG(N)=QEPCNUDF(N)
       QEPCNUDF(N)=0.0
      ENDDO

!NG      write (*,*) ET(143,6),ET(144,6),ET(145,6),ET(146,6)
!NG      write (*,*) ET(143,5),ET(144,5),ET(145,5),ET(146,5)
!NG      write (*,*) ET(143,4),ET(144,4),ET(145,4),ET(146,4)
!NG      write (*,*) ET(143,3),ET(144,3),ET(145,3),ET(146,3)
!NG      write (*,*) ET(143,2),ET(144,2),ET(145,2),ET(146,2)
!wad
!wad      DO I=2,IMM1
!wad         DO J=2,JMM1
!wad            IF(DT(I,J).LE.WETMIN) THEN
!wad               DO K=1,KBM1
!wad                  U(I,J,K)=UTF(I,J)
!wad                  V(I,J,K)=VTF(I,J)
!wad               ENDDO
!wad            ENDIF
!wad         ENDDO
!wad      ENDDO
!wad
!wad
!wad      CALL WANDD_KKY(ET,DT,WETMIN,WETEPS)
!wad
!
9910  continue
!
! copied from ECOMSI, jeff ji , 1/31/94  ------------
!        II = 23
!        JJ = 55
!        II = 31
!        JJ = 10
!        SPD1 = SQRT(.25*(U(II,JJ,1)+U(II+1,JJ,1))**2+.25*(V(II,JJ,1)+
!     *      V(II,JJ+1,1))**2) * 100.
!        ANG1 = ATAN2((V(II,JJ,1)+V(II,JJ+1,1)),(U(II,JJ,1)+
!     *      U(II+1,JJ,1)+1.E-10)) * 57.296
!        SPD2 = SQRT(.25*(U(II,JJ,KBM1)+U(II+1,JJ,KBM1))**2+.25*(
!     *      V(II,JJ,KBM1)+V(II,JJ+1,KBM1))**2) * 100.
!        ANG2 = ATAN2((V(II,JJ,KBM1)+V(II,JJ+1,KBM1)),(U(II,JJ,KBM1)+
!     *      U(II+1,JJ,KBM1)+1.E-10)) * 57.296
!        IMOD = 1
!        If (MOD(INTX,50*IMOD).EQ.0) Write (IUPRT,9006)
!        If (MOD(INTX,IMOD).EQ.0) Write (IUPRT,8900) time,INTX, 100. *
!     *      ELF(II,JJ), SPD1, ANG1, SPD2, ANG2, S(II,JJ,1),
!     *      S(II,JJ,KBM1), T(II,JJ,1), T(II,JJ,KBM1)
!C
!8900  Format (f8.3,I7,5F5.0,4F9.5)
!9006  Format (2X,'  DAYS     TS ELE Usurf ANG Ubot  ANG',
!     &'  SSURF    Sbot',
!     *'     Tsurf    Tbot ')
! copied from ECOMSI, jeff ji , 1/31/94 ------------------
!
**      MOVED ABOVE TO BE CONSITENT WITH  ecom3d  Khan 080999
!qa      CALL ARCHIVE        
**      MOVED ABOVE TO BE CONSITENT WITH  ecom3d  Khan 080999
!
      IF(MOD(INTX,IPRINT).EQ.0 .AND. INTX.GE.IPRTSTART)  THEN
       CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
      ENDIF

!NG Generate restart file at HRSRESTART HOURS
!NG!!!!!!!! REMEMBER: WRITE STATEMENTS NEED INTX, NOT IEND FOR THIS !!!!
!NG
      IF (NRESEND.GT.0) THEN
       IRESTART=INTX_HOT+HRSRESTART(NRESTART)*3600/DTI
      ELSE
       IRESTART=-1
      ENDIF
!NG    WRITE(*,*) 'INTX=', INTX, 'IRESTART=', IRESTART
      IF (INTX.EQ.IRESTART) THEN
       IF (NRESEND.LE.1) THEN
        OPEN (IUWRS,FORM='unformatted',FILE='startup')
       ELSE
        STARTFL="startup.00"
        if (NRESTART.GT.9) THEN
         write (STARTFL(9:10),'(i2)') NRESTART
        else
         write (STARTFL(10:10),'(i1)') NRESTART
        endif
        OPEN (IUWRS,FORM='unformatted',FILE=STARTFL)
       ENDIF
       IF (INTX.NE.IEND) NRESTART=NRESTART+1
!
      IF (TRACER.EQ.'INCLUDE') THEN
       IF (SEDTRAN.EQ.'NEGLECT') THEN
!hli
        IF(WAVEDYN.EQ.'NEGLECT') THEN
         WRITE(IUWRS) 
     .     INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     CONC1,CONC1B
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
        ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
         IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
          WRITE(IUWRS) 
     .      INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      CONC1,CONC1B,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
         ENDIF
        ENDIF
!hli
       ELSE
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +     SEDMASS,SEDEP
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +      SEDMASS,SEDEP,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ELSE
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,TAU
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,TAU,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ENDIF
       ENDIF
      ELSE
       IF (SEDTRAN.EQ.'NEGLECT') THEN
!hli
        IF(WAVEDYN.EQ.'NEGLECT') THEN
         WRITE(IUWRS) 
     .     INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
        ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
         IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
          WRITE(IUWRS) 
     .      INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
         ENDIF
        ENDIF
!hli
       ELSE
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +     SEDMASS,SEDEP
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +      SEDMASS,SEDEP,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ELSE
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,TAU
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      INTX,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,TAU,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ENDIF
       ENDIF
      ENDIF
!
!********************************************************************
!
      CLOSE (IUWRS)
! 
      WRITE(*,*) 'RESTART FILE FINISHED, INTX=,IRESTART=',INTX,IRESTART
      endif
!
 9000              CONTINUE
!
!NG   GENERATE RESTART HERE ONLY IF REQUESTED AT END-OF-FILE
      IF (IEND.LE.IRESTART.OR.NRESEND.EQ.0) THEN
       IF (NRESEND.LE.1) THEN
        OPEN (IUWRS,FORM='unformatted',FILE='startup')
       ELSE
        STARTFL="startup.00"
        if (NRESTART.GT.9) THEN
         write (STARTFL(9:10),'(i2)') NRESTART
        else
         write (STARTFL(10:10),'(i1)') NRESTART
        endif
        OPEN (IUWRS,FORM='unformatted',FILE=STARTFL)
       ENDIF
!NG
!
      IF (TRACER.EQ.'INCLUDE') THEN
       IF (SEDTRAN.EQ.'NEGLECT') THEN
!hli
        IF(WAVEDYN.EQ.'NEGLECT') THEN
         WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     CONC1,CONC1B
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
        ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
         IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
          WRITE(IUWRS) 
     .      IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      CONC1,CONC1B,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
         ENDIF
        ENDIF
!hli
       ELSE
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +     SEDMASS,SEDEP
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +      SEDMASS,SEDEP,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ELSE
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,TAU
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      CONC1,CONC1B,N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,TAU,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ENDIF
       ENDIF
      ELSE
       IF (SEDTRAN.EQ.'NEGLECT') THEN
!hli
        IF(WAVEDYN.EQ.'NEGLECT') THEN
         WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
        ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
         IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
          WRITE(IUWRS) 
     .      IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
         ENDIF
        ENDIF
!hli
       ELSE
         IF (CHEMTRAN.EQ.'INCLUDE') THEN
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +     SEDMASS,SEDEP
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,CHEM1B,CHEM2B,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,CHEM1,CHEM2,CBEDCHEM,CHEMMASS,
     +      SEDMASS,SEDEP,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ELSE
!hli
          IF(WAVEDYN.EQ.'NEGLECT') THEN
           WRITE(IUWRS) 
     .     IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,TXBOT,TYBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     TDIF,SDIF,
     +     N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +     TAUMAX,TAUCUR,TSED,
     +     EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +     FR,NCNT,NHRCNT,NDTCNT,TAU
     +     ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +     ,DUMWAD,DVMWAD ! NGW&D
          ELSE
!NG02222007          IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
           IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
            WRITE(IUWRS) 
     .      IEND,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,  
     .      ART,ARU,ARV,DUM,DVM,FSM,FSMW,COR,CURV42D,TXBOT,TYBOT,
     .      UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .      VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,AAM2D,AAM,KM,KH,KQ,Q2,Q2B,
     .      Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .      TDIF,SDIF,
     +      N24CNT,LAYER,CSED1,CSED2,CSED1B,CSED2B,
     +      TAUMAX,TAUCUR,TSED,
     +      EBTOT,EBMAX,EBCUR,PSED1,PSED2,BEDTH,FRAC0,ACTLAY,CARMOR,
     +      FR,NCNT,NHRCNT,NDTCNT,TAU,
     .      XMOM,YMOM,CO,CS,SN,SIG,WN,NR ! NG FOR WAVES
     +      ,VARWIF,AREAICE,TXICE,TYICE ! NG FOR ICE
     +      ,DUMWAD,DVMWAD ! NGW&D
           ENDIF
          ENDIF
!hli
         ENDIF
       ENDIF
      ENDIF
!
!********************************************************************
!
      CLOSE (IUWRS)
!NG
      ENDIF
!NG
! 
 9100 CONTINUE
!       CLOSE (IUT90)
!       CLOSE (IUT91)
!       CLOSE (IUT92)
      CLOSE (IUT93)
!NG04092008      CLOSE (IUT193)
!       CLOSE (IUT94)
!
      IF (TRACER.EQ.'INCLUDE') THEN
!         CLOSE (IUT96)
!PL050800        CLOSE (IUT97)
        CLOSE (IUT98)
        CLOSE (IUT99)
      ENDIF
!
!PL 050800      IF (SEDTRAN.EQ.'INCLUDE') THEN
!PL 050800         CLOSE (IUT101)
!PL 050800         CLOSE (IUT102)
!PL 050800       ENDIF
!
!********************************************************************
!
      CALL SYSTEM ('rm gcm_temp*')
!
!      IF(VAMAX.LT.10.) THEN
!      IF(HYDTYPE.EQ.'EXTERNAL'.OR.(VAMAX.LT.10.and.elminx.gt.0.1)) THEN ! Ji, 1/17/96
      IF(HYDTYPE.EQ.'EXTERNAL'.OR.(VAMAX.LT.10)) THEN  ! Remove elminx for 
                                                       ! WAD (hli, 08/22/05)
      WRITE(IUPRT,602) TIME
602   FORMAT(/2X,'JOB SUCCESSFULLY COMPLETED; TIME = ',1P1E10.2,' DAYS',
     . //'PLEASE MAKE SURE YOUR RESULTS ARE PHYSICALLY CORRECT')
!
      ELSE
      CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
      WRITE(IUPRT,'(''*********************************************'')')
      WRITE(IUPRT,'(''************* ABNORMAL JOB END **************'')')
      WRITE(IUPRT,'(''************** RUN TERMINATED ***************'')')
      WRITE(IUPRT,'(''** PLEASE CHECK YOUR INPUTS AND PARAMETERS **'')')
      WRITE(IUPRT,'(''*********************************************'')')
      if(vamax.gt.10.0) 
     .WRITE(IUPRT,7500)IAMAX,JAMAX,VAMAX,TIME,FSM(IAMAX,JAMAX)   ! Ji, 1/17/96
 7500 FORMAT(1X,'I=',I3,3X,'J=',I3,3X,'MAX SPEED=',E12.3,'TIME=',
     .F10.4,'FSM= ',f5.0)
      if(elminx.lt.0.1) 
     .WRITE(IUPRT,7501)IAMAX,JAMAX,elmin,TIME,FSM(IAMAX,JAMAX)   ! Ji, 1/17/96
 7501 FORMAT(1X,'I=',I3,3X,'J=',I3,3X,'MIN ELEV=',E12.4,'TIME=',
     .F10.4,'FSM= ',f5.0)
      ENDIF
!--------------------------------------------------------------------
 6110 FORMAT(/' BCTYPE = ',A7,', IS INCORRECTLY SPECIFIED.'//
     *    ' PLEASE FIX AND RESUBMIT'//)
 6111 FORMAT(/' BCTYPE = ',A7,' AND TLAG = ',F10.2,
     *    ' ARE INCORRECTLY SPECIFIED.'//
     *    ' PLEASE FIX AND RESUBMIT'//)
 6112 FORMAT(/' WAVEDYN = ',A8,' IS INCORRECTLY SPECIFIED.'//
     *    ' PLEASE FIX AND RESUBMIT'//)


      STOP
      END
!
! Persistent team for the barotropic, no-transport path only.
! Every worker takes the same early-abort branch after the SINGLE
! that publishes the extrema. IEXT keeps its original exit value.
      SUBROUTINE BARO_SUBCYCLE(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,
     &                         VAMAX,ELMINX,IAMAX,JAMAX)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM),TRNU(IM,JM),TRNV(IM,JM)
      REAL ROWVMAX(JM),ROWEMIN(JM)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,ISEQ)
      DO ISEQ=1,ISPLIT
! No worker reads IEXT before the first flux-stage barrier inside
! EXTRNL_TEAM. The preceding substep publishes its own check first.
!$OMP MASTER
        IEXT=ISEQ
!$OMP END MASTER
        CALL EXTRNL_TEAM(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,.TRUE.)
        CALL WANDD_TEAM(EL,D)
        CALL SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
        IF(VAMAX.GT.100.) EXIT
      ENDDO
!$OMP END PARALLEL
      IF(.NOT.(VAMAX.GT.100.)) IEXT=ISPLIT+1
      END
