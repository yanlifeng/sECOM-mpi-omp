      SUBROUTINE WAVEDON
!******************************PROGRAM WAVE***************************  
! 
! PURPOSE:     TO PREDICT SURFACE GRAVITY WAVES FROM THE WIND.      
! 
! VARIABLES  :
!     THERE ARE SIX MAJOR ARRAYS USED TO DESCRIBE THE WAVE FIELD.
!     THEY ARE DIMENSIONED (IM,JM). THE SIX ARRAYS ARE:
!
!          XMOM - THE X-COMPONENT OF MOMENTUM (units of m*s)
!          YMOM - THE Y-COMPONENT OF MOMENTUM
!          CO   - THE MAGNITUDE OF THE MEAN PHASE SPEED (m/s)
!          SIG  - THE WAVE HEIGHT STANDARD DEVIATION (m)
!          CS   - THE COSINE OF THE MEAN WAVE DIRECTION
!          SN   - THE SINE OF THE MEAN WAVE DIRECTION
!
!     THE THREE OTHER ARRAYS ARE:
!
!          WSX2D and WSY2D   - THE X AND Y COMPONENTS OF THE WIND (m/s)
!                              positive wsx2d is a westerly wind and
!                              positive wsy2d is a southerly wind.
!
! HISTORY :  
!     WRITTEN BY J.R. BENNETT, 1982, GLERL, ANN ARBOR,MI., BY  
!     MODIFYING A PROGRAM WRITTEN BY M.A. DONELAN, NATIONAL WATER  
!     RESEARCH INSTITUTE,BURLINGTON,ONTARIO. 
! 
!     *********************************************************************
!     MODIFIED TO INCORPORATE IN HYDROQUAL'S ECOM MODEL IN CURVILINEAR GRID
!     By QUAMRUL AHSAN, QA, October, 1998
!
!     HydroQual, Inc
!     1 Lethbridge Plaza
!     Mahwah, NJ 07430
!     e-mail qahsan@hydroqual.com
!     *********************************************************************
!
      INCLUDE 'comdeck'
!
!   PHYSICAL CONSTANTS
!
      SAVE 
      REAL FSMWAD(IM,JM) !NG12202010 FOR WAD FOR WAVES
      DATA VK/0.41/ ! NG12172009 Von Karman to 0.41 from 0.4
!hli      DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.0001/
!      DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.00005/
!sfan change gamms to 0.028*2 
!NG09192006      DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.056,0.00005/
!NG09192006      Changed gamma back to 0.028, as ~3% (Lin, et.al. 2002 Continental Self Research pg 2677)
!NG09192006      and changed formulation below for consistency
!NG09192006      Also, a factor of 1/2 was brought in the bottom friction equation for consistency with linear
!NG09192006      wave theory. Thus, the opening parameter CBF is also adjusted.
!NG05152008      DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.0002/
      DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.028,0.0002/ ! 2.8% is the standard
      !DATA RHOAIR,RHOH2O,GAMMA,CBF/1.2233,1000.,0.05,0.0002/ ! 5% of the total momentum to water goes to waves... upper Donelan limit

!
!   COMPUTE CONSTANTS
!
!      write(*,*) 'start wavedon, nr= ', nr
!NG09192006      Changed gamma back to 0.028, as ~3% (Lin, et.al. 2002 Continental Self Research pg 2677)
!NG09192006      and changed formulation below for consistency
!NG09192006      FAC2=0.5*GAMMA*RHOAIR/(RHOH2O*G)
!NG04232008 HOWEVER THERE IS A 0.5 IN THE ORIGINAL CODE, WITH A 0.028 FAC2=GAMMA*RHOAIR/(RHOH2O*G)
!NG09192006      At sea level: g=9.780327[1+0.0053024sin^2(lat)-0.0000058sin^2(2lat)]
      G=GRAV ! NG04272010 Feed from ecom3d for consistency
      FAC2=GAMMA*RHOAIR/(RHOH2O*G) ! 2.8% of the wave momentum is in
      PI= 3.141593
      TPI=6.283185
      GTPI=G/TPI
!--------------------------------------------------------------------------
!  SET INITIAL CONDITIONS   THIS SHOULD GO TO first.f
!
      IF(NR.EQ.0)THEN
        DO J=1,JM
          DO I=1,IM
            XMOM(I,J)=0.0
            YMOM(I,J)=0.0
            CO(I,J)=0.0
            CS(I,J)=1.0
            SN(I,J)=0.0
            SIG(I,J)=0.0
            WN(I,J)=1.
          ENDDO 
        ENDDO 
!hli
        DO J=1,JM
           DO I=1,IM
              FSMW(I,J)=FSM(I,J)
           ENDDO
        ENDDO
        DO N=1,NUMQBC
           IC=IQC(N)
           JC=JQC(N)
           FSMW(IC,JC)=0.0
        ENDDO
!     
!sjf  let FSMW on the open boundary = 0 to void calculation on the 
!     open boundary. Wave parameters on the open boundary are forced
!     

        DO N=1,NUMEBC
           FSMW(IETA(N),JETA(N))=0.0
        ENDDO
        DO N=1,NUMECR ! SPECIAL BOUNDARY CORNER
           FSMW(ICHARD(N),JCHARD(N))=0.0
        ENDDO
!        WRITE(IUNIT,*) 'print FSMW'
!NG        call print(fsm,fsm,im,jm,var,1.0,10,'SRC')
      ENDIF
!hli
!
!   THIS PROGRAM USES TWO TIME STEPS. THE LARGE TIME STEP (DTI)
!   IS SPECIFIED BY THE USER AND CONTROLS THE FREQUENCY AT WHICH DATA
!   CAN BE PRINTED OR STORED FOR LATER ANALYSIS BY WOUTP, THE USER SUP-
!   PLIED OUTPUT ROUTINE. THE SMALL TIME STEP (DTT) IS CALCULATED BY
!   THE PROGRAM FROM THE FORMULA FOR NUMERICAL STABILITY OF THE UPWIND
!   FINITE DIFFERENCE SCHEME. THE NUMBER OF SMALL TIME STEPS PER LARGE
!   TIME STEP IS NDTT; IT IS STORED IN COMMON BLOCK /TPARM/.
!   INSIDE EACH SMALL TIME STEP THERE ARE TWO LOOPS OVER THE SPATIAL
!   COORDINATES I AND J. IN THE FIRST LOOP THE WAVE ADVECTION TERMS ARE
!   CALCULATED. IN THE SECOND LOOP THE WIND INPUT OF MOMENTUM IS ADDED
!   AND THE DIAGNOSTIC VARIABLES ARE CALCULATED.
!
!  CALCULATE SMALL TIME STEP BASED ON CONSTANT WIND FROM MIDDLE OF TIME STEP
!

       IF(CREADWBC) THEN ! Use read-in boundary condition
        WMAX=0.
!NG09192006simplified       NMAX=0
!NG09192006simplified       DMAX=99999.

!NG05312012 MOVE THIS BEFORE THE TIME STEP ASSIGNMENT TO ACCOMODATE
!NG         CALCULATION OF "EQUIVALENT-WIND" OF OBC SWELL FOR TIME STEP
!NG12212010 MOVE NUMEBC ASSIGNMENT OUT OF THE 60 LOOP. FSMW ENSURES THESE ARE NOT CHANGED.
        DO N=1,NUMEBC
!NG09192006     WHBFLX: WAVE HEIGHT STANDARD DEVIATION (Ho/4) IN METERS
         SIG(IETA(N),JETA(N))=WHBFLX(N)
!NG09192006     WDBFLX: MOMENTUM FLUX DIRECTION, W.R.T. KSI1-DIRECTION AS
!NG09192006             DONELAN AND SCHWAB'S CONVENTION. Wave is moving towards the
!NG09192006             direction not from, counter clockwise is positive.
!NG01312008        CS(IETA(N),JETA(N))=COSD(WDBFLX(N)*PI/180.)
!NG01312008        SN(IETA(N),JETA(N))=SIND(WDBFLX(N)*PI/180.)
         CS(IETA(N),JETA(N))=COS(WDBFLX(N)*PI/180.) !NG01312008 standard FORTRAN does not have COSD
         SN(IETA(N),JETA(N))=SIN(WDBFLX(N)*PI/180.) !NG01312008 standard FORTRAN does not have SIND
!NG09192006 Note that open boundary period WPBFLX or its phase speed, CO,
!NG09192006 does not enter the calculations presently, so it is not required.
         WVPD(IETA(N),JETA(N))=WPBFLX(N) !NG06012009 Store it anyway for display of open BC
         WHNG=4*WHBFLX(N)
        WMAX=AMAX1(WMAX,0.0065*WHNG*WHNG*WHNG-0.23*WHNG*WHNG+4.*WHNG+2.)! CNG12212010 BASED ON AN UPPER BOUND OF BEAUFORT
         WMAX=AMAX1(WMAX,GTPI*WPBFLX(N))! CNG12212010 BASED ON LIMIT OF PHASE SPEED
        ENDDO
!NG05312012 MOVE THIS BEFORE THE TIME STEP ASSIGNMENT TO ACCOMODATE

       ELSE ! assign from internal solution (open boundary condition).
        WMAX=0.
        DO N=1,NUMEBC
!NG10292013 FOR NO-GRADIENT BOUNDARY CONDITION, TAKE BOUNDARY FROM PRIOR INNER SOLUTION:
         WHBFLX(N)=SIG(ICON(N),JCON(N))
         WPBFLX(N)=WVPD(ICON(N),JCON(N))
         WDBFLX(N)=90. - WVDR(ICON(N),JCON(N))
     +                 - ANG(ICON(N),JCON(N))*180./PI ! KSI-1 DIR
!NG10292013 AND THEN DO THE SAME AS ABOVE
         SIG(IETA(N),JETA(N))=WHBFLX(N)
         CS(IETA(N),JETA(N))=COS(WDBFLX(N)*PI/180.) !NG01312008 standard FORTRAN does not have COSD
         SN(IETA(N),JETA(N))=SIN(WDBFLX(N)*PI/180.) !NG01312008 standard FORTRAN does not have SIND
         WVPD(IETA(N),JETA(N))=WPBFLX(N) !NG06012009 Store it anyway for display of open BC
         WHNG=4*WHBFLX(N)
        WMAX=AMAX1(WMAX,0.0065*WHNG*WHNG*WHNG-0.23*WHNG*WHNG+4.*WHNG+2.)! CNG12212010 BASED ON AN UPPER BOUND OF BEAUFORT
         WMAX=AMAX1(WMAX,GTPI*WPBFLX(N))! CNG12212010 BASED ON LIMIT OF PHASE SPEED
        ENDDO
       ENDIF

       DS=99999.
!NG       if (nr.eq.0) then
!      call print(wsx2d,fsm,im,jm,var,1.0,10,'SRC')
!      call print(wsy2d,fsm,im,jm,var,1.0,10,'SRC')
!NG      endif
!      write(1007) wsx2d,wsy2d
       DO J=2,JMM1
         DO I=2,IMM1
           IF(FSMW(I,J).NE.0.) THEN 
             IF(DT(I,J).LE.0) THEN
              WRITE(*,*)'WAVE ERROR',I,J,'DT<=0',DT(I,J)
              WRITE(IUPRT,*)'WAVE ERROR',I,J,'DT<=0',DT(I,J)
              CALL PRINTS(DRHOX,DRHOY,TRNU,TRNV)
              STOP
             ENDIF
             DS  =AMIN1(DS,SQRT(H1(I,J)**2 + H2(I,J)**2))
             WMAX=AMAX1(WMAX,SQRT(WSX2D(I,J)**2 + WSY2D(I,J)**2)
     +                 ,ABS(CO(I,J)))
!             write(*,*) 'wmax,win,co(i,j),wsx2d(i,j),wsy2d(i,j)',
!     &            wmax,win,co(i,j),wsx2d(i,j),wsy2d(i,j),i,j
!NG09192006simplified             DS  = AMIN1(DMAX,DS)
!NG09192006simplified             DMAX = DS
           ENDIF
         ENDDO       
       ENDDO       
!NG    Calculate number of timesteps: sqrt(2)*maxwindspeed*timeinterval/smallesthypotenuse + 2
       NDTT=IFIX((1.414*WMAX*DTI*NWAVE)/DS)+2
!NG09192006simplified       NDTT=MAX0(NMAX,NDTT)
       DTT=DTI*FLOAT(NWAVE)/NDTT
!NG       WRITE (*,*) NDTT,DTT
       DTTGAM=FAC2*DTT
!
!  BEGIN LOOP OVER SMALL TIME STEPS
!
!      write(1005,*) 'time=',time
!      write(1005,'(8f10.5)') (WHBFLX(N),n=1,numebc)
!      write(1005,'(8f10.5)') (WDBFLX(N),n=1,numebc)
!       write(*,*) 'win,ds',win,ds
!       write(*,*) 'wmax,dtt,ndtt,nwave,dti = ',wmax,dtt,ndtt,nwave,dti

!
!Based on sfan's (SIT) input, force sig and cs, sn on open boundaries
!(hli,03/19/04)
!

      DO J=1,JMM1
      DO I=1,IMM1
       FSMWAD(I,J)=DUMWAD(I+1,J)+DUMWAD(I,J)+DVMWAD(I,J+1)+DVMWAD(I,J) ! NG12202010 DEFINE WAVE WAD MASK HERE 0=DRY
      ENDDO
      ENDDO

      DO 60 NDT=1,NDTT
!
!   FIRST LOOP OVER I,J: CALCULATE THE WAVE ADVECTION TERMS
!
!      call smooth(sig,fsmw,im,jm,1)
!      call smooth(cs,fsmw,im,jm,1)
!      call smooth(sn,fsmw,im,jm,1)

      DO 70 I=2,IMM1
      DO 70 J=2,JMM1
       IF(FSMW(I,J)*FSMWAD(I,J).EQ.0) GOTO 70 ! NG12202010 CALCULATION ONLY FOR INTERNAL, NON-DRY, NON-BOUNDARY ELEMENTS
       DTTDS=DTT*0.25/(H1(I,J)*H2(I,J))
!
!   X-COMPONENT    UPWIND SCHEME
!
       IF(CS(I,J).GE.0.0)
!NG03022007 For Thin Dams     .  XFLXD=H2(I,J)*((SIG(I,J)*CS(I,J))**2)
!NG03022007 For Thin Dams     .  -H2(I-1,J)*((SIG(I-1,J)*CS(I-1,J))**2)
     .  XFLXD=          H2(I,J)*((SIG(I,J)*CS(I,J))**2)
     .                  +DUM(I,J)*DUMWAD(I,J)*
     .                  (-H2(I-1,J)*((SIG(I-1,J)*CS(I-1,J))**2))
       IF(CS(I,J).LT.0.0)
!NG03022007 For Thin Dams     .  XFLXD=H2(I+1,J)*((SIG(I+1,J)*CS(I+1,J))**2)
!NG03022007 For Thin Dams     .  -H2(I,J)*((SIG(I,J)*CS(I,J))**2)
     .  XFLXD=DUM(I+1,J)*DUMWAD(I+1,J)*
     .                   (H2(I+1,J)*((SIG(I+1,J)*CS(I+1,J))**2))
     .                   -H2(I,J)*((SIG(I,J)*CS(I,J))**2)
       IF(SN(I,J).LT.0.0)
!NG03022007 For Thin Dams     .  YFLXD=H1(I,J+1)*(CS(I,J+1)*SN(I,J+1)*SIG(I,J+1)**2)
!NG03022007 For Thin Dams     .  -H1(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
     .  YFLXD=DVM(I,J+1)*DVMWAD(I,J+1)*
     .                   (H1(I,J+1)*(CS(I,J+1)*SN(I,J+1)*SIG(I,J+1)**2))
     .                   -H1(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
       IF(SN(I,J).GE.0.0)
!NG03022007 For Thin Dams     .  YFLXD=H1(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
!NG03022007 For Thin Dams     .  -H1(I,J-1)*(CS(I,J-1)*SN(I,J-1)*SIG(I,J-1)**2)
     .  YFLXD=          H1(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
     .                  +DVM(I,J)*DVMWAD(I,J)*
     .                  (-H1(I,J-1)*(CS(I,J-1)*SN(I,J-1)*SIG(I,J-1)**2))

       SIP=SIG(I+1,J)
       SIM=SIG(I-1,J)
       H2F=H2 (I+1,J)
       H2B=H2 (I-1,J)
!NG03022007 For Thin Dams      IF(FSMW(I+1,J).EQ.0.) THEN
!NG03022007 For Thin Dams         SIP=2.*SIG(I,J)-SIG(I-1,J)
!NG03022007 For Thin Dams         H2F=2.*H2 (I,J)-H2 (I-1,J)
       IF(FSMW(I+1,J)*DUM(I+1,J)*DUMWAD(I+1,J).EQ.0.) THEN
!NG       IF(FSM(I+1,J)*DUM(I+1,J)*DUMWAD(I+1,J).EQ.0.) THEN
        SIP=2.*SIG(I,J)-DUM(I,J)*DUMWAD(I,J)*SIG(I-1,J)
        H2F=2.*H2 (I,J)-DUM(I,J)*DUMWAD(I,J)*H2 (I-1,J)
       ENDIF
!NG03022007 For Thin Dams      IF(FSMW(I-1,J).EQ.0.) THEN
!NG03022007 For Thin Dams         SIM=2.*SIG(I,J)-SIG(I+1,J)
!NG03022007 For Thin Dams         H2B=2.*H2 (I,J)-H2 (I+1,J)
       IF(FSMW(I-1,J)*DUM(I,J)*DUMWAD(I,J).EQ.0.) THEN
!NG       IF(FSM(I-1,J)*DUM(I,J)*DUMWAD(I,J).EQ.0.) THEN
        SIM=2.*SIG(I,J)-DUM(I+1,J)*DUMWAD(I+1,J)*SIG(I+1,J)
        H2B=2.*H2 (I,J)-DUM(I+1,J)*DUMWAD(I+1,J)*H2 (I+1,J)
       ENDIF
       XMOM(I,J)=XMOM(I,J)-DTTDS*(XFLXD+YFLXD+0.25*(H2F*SIP**2
     .          -H2B*SIM**2))
!
!   Y-COMPONENT    UPWIND SCHEME
!
       IF(CS(I,J).LT.0.0)
!NG03022007 For Thin Dams     .  XFLXD=H2(I+1,J)*(CS(I+1,J)*SN(I+1,J)*SIG(I+1,J)**2)
!NG03022007 For Thin Dams     .  -H2(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
     .  XFLXD=DUM(I+1,J)*DUMWAD(I+1,J)*
     .                   (H2(I+1,J)*(CS(I+1,J)*SN(I+1,J)*SIG(I+1,J)**2))
     .                   -H2(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
       IF(CS(I,J).GE.0.0)
!NG03022007 For Thin Dams     .  XFLXD=H2(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
!NG03022007 For Thin Dams     .  -H2(I-1,J)*(CS(I-1,J)*SN(I-1,J)*SIG(I-1,J)**2)
     .  XFLXD=           H2(I,J)*(CS(I,J)*SN(I,J)*SIG(I,J)**2)
     .       +DUM(I,J)*DUMWAD(I,J)*
     .                 (-H2(I-1,J)*(CS(I-1,J)*SN(I-1,J)*SIG(I-1,J)**2))
       IF(SN(I,J).GE.0.0)
!NG03022007 For Thin Dams     .  YFLXD=H1(I,J)*((SN(I,J)*SIG(I,J))**2)
!NG03022007 For Thin Dams     .  -H1(I,J-1)*((SN(I,J-1)*SIG(I,J-1))**2)
     .  YFLXD=           H1(I,J)*((SN(I,J)*SIG(I,J))**2)
     .       +DVM(I,J)*DVMWAD(I,J)*
     .                 (-H1(I,J-1)*((SN(I,J-1)*SIG(I,J-1))**2))
       IF(SN(I,J).LT.0.0)
!NG03022007 For Thin Dams     .  YFLXD=H1(I,J+1)*((SN(I,J+1)*SIG(I,J+1))**2)
!NG03022007 For Thin Dams     .  -H1(I,J)*((SN(I,J)*SIG(I,J))**2)
     .  YFLXD=DVM(I,J+1)*DVMWAD(I,J+1)*
     .                   (H1(I,J+1)*((SN(I,J+1)*SIG(I,J+1))**2))
     .                   -H1(I,J)*((SN(I,J)*SIG(I,J))**2)

       SJP=SIG(I,J+1)
       SJM=SIG(I,J-1)
       H1F=H1 (I,J+1)
       H1B=H1 (I,J-1)
!NG03022007 For Thin Dams      IF(FSMW(I,J+1).EQ.0.) THEN
       IF(FSMW(I,J+1)*DVM(I,J+1)*DVMWAD(I,J+1).EQ.0.) THEN
!NG       IF(FSM(I,J+1)*DVM(I,J+1)*DVMWAD(I,J+1).EQ.0.) THEN
        SJP=2.*SIG(I,J)-DVM(I,J)*DVMWAD(I,J)*SIG(I,J-1)
        H1F=2.*H1 (I,J)-DVM(I,J)*DVMWAD(I,J)*H1 (I,J-1)
       ENDIF
!NG03022007 For Thin Dams      IF(FSMW(I,J-1).EQ.0.) THEN
       IF(FSMW(I,J-1)*DVM(I,J)*DVMWAD(I,J).EQ.0.) THEN
!NG       IF(FSM(I,J-1)*DVM(I,J)*DVMWAD(I,J).EQ.0.) THEN
        SJM=2.*SIG(I,J)-DVM(I,J+1)*DVMWAD(I,J+1)*SIG(I,J+1)
        H1B=2.*H1 (I,J)-DVM(I,J+1)*DVMWAD(I,J+1)*H1 (I,J+1)
       ENDIF
       YMOM(I,J)=YMOM(I,J)-DTTDS*(XFLXD+YFLXD+0.25*(H1F*SJP**2
     .          -H1B*SJM**2))
!
   70 CONTINUE
!      call smooth(xmon,fsmw,im,jm,1)
!      call smooth(ymon,fsmw,im,jm,1)
!
!   SECOND LOOP OVER I,J: CALCULATE THE WIND INPUT OF MOMENTUM
!   AND COMPUTE THE DIAGNOSTIC VARIABLES (CS,SN,S,C)
!
      DO 90 J=2,JMM1
      DO 90 I=2,IMM1
       IF(FSMW(I,J)*FSMWAD(I,J).EQ.0) GOTO 90 ! NG12202010 CALCULATION ONLY FOR INTERNAL, NON-DRY, NON-BOUNDARY ELEMENTS
!
!   CALCULATE THE VERTICAL MOMENTUM FLUX IN THE 2 DIRECTIONS,
!   THE WIND AND THE WAVE FIELD
!NG09192006 NOTE WIND IS ROTATED TO KSI-1 AND KSI-2 DIRECTIONS FROM E-W AND N-S
!NG09192006 TO FORM THE REQUIRED DOT PRODUCT OF WAVES.WIND, AS IN THE ORIGINAL
!NG09192006 SCHWAB-BENNETT CODE.
!
!     wsx2d and wsy2d are E-W and N-S components. positive u is a westerly
!     wind and positive v is a southerly wind.
!NG04232008 UPDATED FORMULATION TO BRING IT CLOSER TO THE SCHWAB ORIGINAL GLERL-51
!NG04232008 ALSO TRIED THE GLERL-51 FORMULATION (NOT CODE), BUT IT CRASHES ON HIGH WINDS
       UWIND=(WSX2D(I,J)*COS(ANG(I,J))+
     +        WSY2D(I,J)*SIN(ANG(I,J)))
       VWIND=(-WSX2D(I,J)*SIN(ANG(I,J))+
     +        WSY2D(I,J)*COS(ANG(I,J)))
       WSPDSQ=UWIND*UWIND+VWIND*VWIND
       WNDSPD=SQRT(WSPDSQ)+1.E-8 ! NG12172009 increased from 1.E-10
       CSWIND=UWIND/WNDSPD
       SNWIND=VWIND/WNDSPD
       CTHE1=CSWIND*CS(I,J)+SNWIND*SN(I,J)
       B=CO(I,J)*.83/WNDSPD
       SG=AMAX1(0.005,SIG(I,J)*ABS(CTHE1)) ! NG The magnitude of covariance enters the drag coeff...
                                           ! so that if wind is normal to waves, minimal drag.
       Z0WAVE=AMIN1(SG,49.999)/5. ! NG12172009 Original surface roughness Z0 in GLERL
!NG12172009 TRIED, BUT GAMMA SEEMS TO NEED ADJUSTMENT TO 6-7%... Z0WAVE=1.38E-4*4.*SIG(I,J)*(WNDSPD/CO(I,J))**2.66 +1.E-5 ! NG12172009 Z0 from Donelan 1990, a function of inverse wave age
       DRAG1=(VK/ALOG(10./Z0WAVE))**2
!NG03232007 NG's TAKE OF SCHWAB'S PUBLISHED METHOD (GLERL-51) BASED ON FORMULA PAGE 3. CRASHES.
!NG03232007      AA=SQRT(1+B*B-2*B*CTHE1)
!NG03232007      WI1=DRAG1*AA/WNDSPD
!NG03232007      WI2=-B*WI1
       AA=1.-B*CTHE1
       WI1=DRAG1*AA*ABS(AA)
       AA=CTHE1-B ! AS PER ORIGINAL
       WI2=DRAG1*AA*ABS(AA)
!NG09192006      VFLXA=WI1*CSWIND+WI2*CS(I,J)
!NG09192006      VFLYA=WI1*SNWIND+WI2*SN(I,J)
!NG09192006 INCLUDE * WSPSSQ HERE, NOT IN DTFAC
!NG09192006 ADD SKIN FRICTION CONTRIBUTION
       DRAGS = 0.7E-03              ! SKIN DRAG COEFFICIENT
       VFLXA = DRAGS*WNDSPD*UWIND +
     +        (WI1*CSWIND+WI2*CS(I,J)) * WSPDSQ
       VFLYA = DRAGS*WNDSPD*VWIND +
     +        (WI1*SNWIND+WI2*SN(I,J)) * WSPDSQ
!
!   COMPUTE BOTTOM FRICTION LOSS FOR SHALLOW WATER USING STANDARD DEVIATION
!    OF WAVEHEIGHT FOR AMPLITUDE AND WAVENUMBER OF PEAK ENERGY FREQUENCY
!
       IF(DT(I,J)*WN(I,J).LT.10.) THEN
!NG09192006       UBOT=SIG(I,J)*WN(I,J)*CO(I,J)/SINH(WN(I,J)*DT(I,J))
!NG09192006 A 1/2 FACTOR WAS BROUGHT IN THE EQUATION FOR CONSISTENCY WITH LINEAR
!NG09192006 WAVE THEORY. THUS, THE OPENING PARAMETER CBF IS ALSO ADJUSTED.
!NG09192006 10.0 is just an arbitrary limit for when 1/sinh->0
        UBOT=SIG(I,J)*WN(I,J)*CO(I,J)/2./(SINH(WN(I,J)*DT(I,J))+1.E-10)
       ELSE
        UBOT=0.
       ENDIF
!NG09192006      DTFAC=DTTGAM*WSPDSQ
       DTFAC=DTTGAM
       XMOM(I,J)=XMOM(I,J)+DTFAC*VFLXA-DTT*CBF*UBOT*UBOT*CS(I,J)
       YMOM(I,J)=YMOM(I,J)+DTFAC*VFLYA-DTT*CBF*UBOT*UBOT*SN(I,J)
   90 CONTINUE
!
!     using laplacian smoother every time step
!      call smooth(xmon,fsmw,im,jm,1)
!      call smooth(ymon,fsmw,im,jm,1)
      
      DO 91 J=2,JMM1
      DO 91 I=2,IMM1
       IF(FSMW(I,J)*FSMWAD(I,J).EQ.0) GOTO 91 ! NG12202010 CALCULATION ONLY FOR INTERNAL, NON-DRY, NON-BOUNDARY ELEMENTS
!     
!   COMPUTE THE NEW COSINES AND SINES; AND COMPUTE THE PHASE
!   SPEED (C) AND WAVE HEIGHT STANDARD DEVIATION (S).

!NG09192006 This uses Hasselmann et al. (1973) JONSWAP spectrum:
!NG09192006 Joint North Sea Wave Observation Project

       CM=SQRT(XMOM(I,J)**2+YMOM(I,J)**2)+1.E-10
       CS(I,J)=XMOM(I,J)/CM
       SN(I,J)=YMOM(I,J)/CM
!NG12172009 BIG ERROR: UWIND and VWIND were not defined as (I,J) matrices? Added...
       UWIND=(WSX2D(I,J)*COS(ANG(I,J))+
     +        WSY2D(I,J)*SIN(ANG(I,J)))
       VWIND=(-WSX2D(I,J)*SIN(ANG(I,J))+
     +        WSY2D(I,J)*COS(ANG(I,J)))
!NG09212006 UU is the wind in the direction of the wave momentum flux
!NG09212006 (simply rotated). Then, if wind opposing the wave, do 150.
       UU=UWIND*CS(I,J)+VWIND*SN(I,J)
       IF (UU.LE.0.0) GO TO 150
       FM=0.01788735*(((UU*UU)/(CM*CM*CM+1.E-10))**0.14285714)
       IF((UU*FM).GT.0.760545) GOTO 160
  150  FM=(CM*14343.09)**(-.3333333333)
  160  CONTINUE
!NG09212006 THROW THIS AFTER IF STATEMENT:SIG(I,J)=SQRT(CM * GTPI/FM)
!NG09192006      WN(I,J)=WNUM(FM,DT(I,J))
!NG12022008      WN(I,J)=WNUM(FM,DT(I,J),G) ! NG09192006 INCLUDE G AS INPUT
       WN(I,J)=WNUMYOU(FM,DT(I,J),G) ! NG12022008 USED BETTER SOLUTION
       IF(DT(I,J)*WN(I,J).LT.3.14) THEN
!NG09192006 LINEAR WAVE THEORY. STRAIGHT FROM THE DISPERSION EQUATION.
!NG04272010 w^2=gktanh(kh) and c=w/k.
!NG09192006 3.14 is just the limit when tanh->1
        CO(I,J)=SQRT(G*TANH(WN(I,J)*DT(I,J))/(WN(I,J)+1.E-10))
       ELSE
        CO(I,J)=GTPI/(FM+1.E-10) ! NG09192006 classic limit w^2=gk, g/(2pif)=g/w=w/k
       ENDIF
!NG09192006      SIG(I,J)=SQRT(CM * GTPI/FM)
!NG12202010      SIG(I,J)=SQRT(CM * CO(I,J)) ! THROWN HERE, AFTER CO IS CALCULATED
       SIG(I,J)=AMAX1(0.,SQRT(CM * CO(I,J))) ! DO NOT ALLOW NEGATIVE HEIGHTS

!
   91 CONTINUE
!      call smooth(sig,fsmw,im,jm,1)
!      call smooth(cs,fsmw,im,jm,1)
!      call smooth(sn,fsmw,im,jm,1)
!
!   END LOOP OVER SMALL TIME STEPS
!
!
!Based on sfan's (SIT) input, force sig and cs, sn on open boundaries
!(hli,03/19/04)
!
!      DO N=1,NUMEBC
!        SIG(IETA(N),JETA(N))=WHBFLX(N)
!        CS(IETA(N),JETA(N))=COSD(WDBFLX(N))
!        SN(IETA(N),JETA(N))=SIND(WDBFLX(N))
!        FM=1/(WPBFLX(N)+1.E-10)
!        WN(IETA(N),JETA(N))=WNUM(FM,DT(IETA(N),JETA(N)))
!        IF(DT(IETA(N),JETA(N))*WN(IETA(N),JETA(N)).LT.3.14) THEN
!           CO(IETA(N),JETA(N))=SQRT(G*TANH(WN(IETA(N),JETA(N))
!     &          *DT(IETA(N),JETA(N)))/WN(IETA(N),JETA(N)))
!        ELSE
!           CO(IETA(N),JETA(N))=GTPI/FM
!        ENDIF
!        XMOM(IETA(N),JETA(N))=0.0
!        YMOM(IETA(N),JETA(N))=0.0
!      ENDDO
! END
   60 CONTINUE
!
!  CALCULATE WAVE PARAMETERS  (Significant Wave Height, Wave Period
!                              and Wave Direction)
! 
      DO 21 J=1,JM  
      DO 21 I=1,IM  
!NG03092007PreserveOpenBC      IF(FSMW(I,J).EQ.0.) GO TO 21
!NG       IF(FSM(I,J)*FSMWAD(I,J).EQ.0) THEN !NG12212010
       IF(FSM(I,J).EQ.0) THEN !NG12212010
        WVHT(I,J)=-99999.0 ! mask out exterior points and currently dry points (but not OBC, unless dry!)
        GOTO 21
       ENDIF
! 2 below worked but took out WAD and Open BC points in output at least...
!       WVHT(I,J)=-999.0
!       IF(FSMW(I,J)*FSMWAD(I,J).EQ.0) GOTO 21
       WVHT(I,J)=4.*SIG(I,J) 
   21 CONTINUE  
!PL
      DO J=1,JM
      DO I=1,IM
!NG       IF(FSMW(I,J)*FSMWAD(I,J).NE.0) THEN ! NG12212010 if interior non-boundary, or wet (now), then
       IF(FSMW(I,J).NE.0) THEN ! NG12212010 if interior non-boundary, or wet (now), then
!NG10202008          HBREAK=0.78*DT(I,J) ! Outdated McCowan, 1894 expression
        HBREAK=0.8261*DT(I,J) ! NG10202008 Use Longuet-Higgins and Fenton, 1974
        !WVHT(I,J)=AMIN1(WVHT(I,J),HBREAK)
        WVHT(I,J)=AMIN1(ABS(WVHT(I,J)),HBREAK) ! NG12212010 FOR PREVIOUSLY DRY MASKED-OUT POINTS THAT ARE NOW WET USE ABS
        SIG(I,J)=WVHT(I,J)/4. ! NG09212006 update SIG too! Be mindful of the 21 loop!
       END IF
      END DO
      END DO
! 
!  SAVE WAVE PERIODS  
! 
      DO 22 J=1,JM  
      DO 22 I=1,IM  
      IF(FSMW(I,J).EQ.0) GOTO 22 ! NG12202010 CALCULATION ONLY FOR INTERNAL, NON-BOUNDARY ELEMENTS
!NG      IF(FSMW(I,J).EQ.0) GOTO 22 ! NG12202010 CALCULATION ONLY FOR INTERNAL, NON-DRY, NON-BOUNDARY ELEMENTS
!NG04272010 THIS WAS DEEP WATER APPROXIMATION:      WVPD(I,J)=TPI*CO(I,J)/GRAV ! NG09192006
!NG04272010 INSTEAD USE: c=w/k => w=ck => 2pi/T=ck => T=2pi/(ck)
!NG12202010      WVPD(I,J)=TPI/(CO(I,J)*WN(I,J)+1.E-10)
       WVPD(I,J)=AMAX1(0.,TPI/(CO(I,J)*WN(I,J)+1.E-10)) ! DO NOT ALLOW NEGATIVE PERIODS
       IF (FSMWAD(I,J).EQ.0) WVPD(I,J)=0. ! ZERO OUT PERIOD OF DRY POINTS FOR OUTPUT
22    CONTINUE  
! 
!   SAVE WAVE DIRECTIONS  
! 
      DO 23 J=1,JM  
      DO 23 I=1,IM  
!NG03092007forOpenBoundaryPlotting      IF(FSMW(I,J).EQ.0.) GO TO 23 
!NG03092007forOpenBoundaryPlotting      UU=XMOM(I,J) 
!NG03092007forOpenBoundaryPlotting      VV=YMOM(I,J) 
!NG      IF(FSM(I,J)*FSMWAD(I,J).EQ.0) GOTO 23 ! NG12202010 CALCULATION ONLY FOR ALL WET ELEMENTS, INCLUDING OBC
      IF(FSM(I,J).EQ.0) GOTO 23 ! NG12202010 CALCULATION FOR ALL WET AND DRY ELEMENTS, INCLUDING OBC
       UU=CS(I,J) 
       VV=SN(I,J) 
       WVDR(I,J) = 0.
       IF(UU .EQ. 0. .AND. VV .EQ. 0.) GOTO 23

!     Donelan and Schwab's convention of wave direction w.r.t. zeta1 axis
!     wave is moving towards the direction not from, counter clock is positive 

       WVDR(I,J)=ATAN2(VV,UU)*180./PI  ! w.r.t. zeta1-direction
!NG       IF(WVDR(I,J) .LT. 0.) WVDR(I,J) = WVDR(I,J) + 360.
!NG       IF(WVDR(I,J) .GT. 360.) WVDR(I,J) = WVDR(I,J) - 360. 

!     Direction with respect to east (x-direction)

       WVDR(I,J) = WVDR(I,J) + ANG(I,J)*180./PI

!     IN stress.f and in SMB Model WAVE DIR IS IN w.r.t North Clock Wise
       WVDR(I,J)= 90.-WVDR(I,J)
       WVDR(I,J) = AMOD(WVDR(I,J),360.)
       IF(WVDR(I,J).LE.0.) WVDR(I,J) = WVDR(I,J) + 360.
!     final wave direction is w.r.t North clock wise and towards (oceanographic)
23    CONTINUE  
!
!NG03072008      NR=NR+1
      NR=1 ! stop counter

!     HARDWIRE for Checking Purposes************
!      WRITE(1001)TIME
!      WRITE(1001)WVHT
!      WRITE(1001)WVPD
!      WRITE(1001)WVDR
!      write(1001) time,wvht(720,79),wvht(574,89),wvht(320,111)
!     &     ,wvht(144,26),wvpd(720,79),wvpd(574,89),wvpd(320,111)
!     &     ,wvpd(144,26),wvdr(720,79),wvdr(574,89)
!     &     ,wvdr(320,111),wvdr(144,26)
!
  540 FORMAT(1X, 'DEPTH RELATIVE TO MEAN WATER LEVEL')
!      write(*,*) 'end wavedon, nr= ', nr
!      write(1008) wsx2d,wsy2d

      RETURN
      END

      FUNCTION WNUM(F,DD,G) ! NG09192006 Include G as input
!NG09192006      FUNCTION WNUM(F,DD)
!NG I HAVE CHECKED THIS VERSUS CETN-1-17 6/85 Coastal Engineering
!NG Technical Note; Direct Methods for Calculating Wavelength.
!NG Provides 0.1% accuracy.
!
!  APPROXIMATE SOLUTION OF WAVE DISPERSION EQUATION
!
!  F = FREQUENCY (HZ)
!  DD = DEPTH (M)
!  WNUM = WAVENUMBER (1/M)
!
!  REFERENCE: HUNT, J.N. 1979. DIRECT SOLUTION OF WAVE DISPERSION EQUATION.
!   JWPCOD, ASCE, 105(WW4):457-459.
!NG ACTUALLY, AS REFINED BY CHEN AND THOMPSON, 1985. Iterative and
!NG Pade solution for the water wave dispersion relation CERC-85-4.
!
      DATA TPI/6.283185/
!NG09192006      DATA G/9.806/
      Y=DD*(TPI*F)**2/G
      X=Y*(Y+1./(1.+Y*(0.6522+Y*(0.4622+Y*Y*(0.0864+Y*0.0675)))))
      WNUM=SQRT(X)/(DD+1.E-10)
      RETURN
      END

      FUNCTION WNUMYOU(F,DD,G) ! NG12022009 
!NG
!NG  APPROXIMATE SOLUTION OF WAVE DISPERSION EQUATION
!NG
!NG  F = FREQUENCY (HZ)
!NG  DD = DEPTH (M)
!NG  WNUMYOU = WAVENUMBER (1/M)
!NG
!NG Reference: Zai-Jin You, 2008
!NG "A close approximation of wave dispersion relation for
!NG direct calculation of wavelength in any coastal water depth."
!NG Equations 3 and 2 (or 3 and 4). 
!NG Applied Ocean Research
!NG Provides 0.01% (0.1%) maximum relative error.
!NG
!NG Note that the YOU (2) equation can be analytically rewritten BUT
!NG the resulting simpler expression becomes numerically indefinite
!NG because it involves ~inf/~inf expressions like cosh*sinh.
!NG
      DATA TPI/6.283185/
      HK0=DD*TPI*F*TPI*F/G
      X0=SQRT(HK0)*(1.+HK0/6.+HK0*HK0/30.) ! YOU (3)
      Y=X0*( (HK0+X0/COSH(X0)*X0/COSH(X0)) /
     +       (X0*TANH(X0)+X0/COSH(X0)*X0/COSH(X0)) ) ! YOU (2) NEWTON-RAPHSON
!NG   Y=HK0/TANH(X0) ! YOU (4), SOMEWHAT FASTER ~ HUNT IN ACCURACY AND SPEED
      WNUMYOU=Y/(DD+1.E-10)
      RETURN
      END

