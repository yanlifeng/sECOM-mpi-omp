      SUBROUTINE SEDIC(RESTAR)
!
!  SETS INITIAL CONDITIONS AND PARAMETERS FOR SEDIMENT TRANSPORT
!
!  REVISION DATE :  SEP 29, 1998  
!
!***********************************************************************
!
      INCLUDE 'comdeck'
!
      CHARACTER*10 RESTAR
      REAL PINIT(IM,JM)
      REAL D50VAR(IM,JM)
!
      WS1MIN=6.6E-06  ! MINIMUM SETTLING VELOCITY FOR COHESIVE SEDIMENT IN m/s
!
!  VARIABLE INITIALIZATION
!
!  BED MASK
!
      IF (RESTAR.EQ.'COLD START') THEN
        DO 4 J=1,JM
          DO 4 I=1,IM
            IBMSK(I,J)=0
            DO 3 K=1,KSMAX
              FALAY(K,I,J)=0.0
 3          CONTINUE
 4      CONTINUE
      ENDIF
!
!  INPUT BED LOAD MASK
!  MASK FOR FINE SEDIMENT AREAS
!
!  IBED = 0 :  NO NON-COHESIVE SUSPENDED LOAD
!  IBED = 1 :  NON-COHESIVE SUSPENDED LOAD
!
      IBED=0
!
        OPEN (UNIT=85,FILE='bed_mask',FORM='FORMATTED')
        DO 5 I=1,IM
          READ (85,6) (IBMSK(I,J),J=1,JM)
!
!  CHECK FOR NON-COHESIVE SUSPENDED LOAD
!
          DO 95 J=1,JM
            IF (IBMSK(I,J).EQ.1) IBED=1
 95       CONTINUE
!
!  SEARCH FOR INCORRECT BED MASK FOR SAND TRANSPORT ONLY
!
          IF (SEDTYPE.EQ.'SAND') THEN
            DO 96 J=1,JM
              IF (IBMSK(I,J).EQ.0) THEN
                WRITE (IUPRT,97)I,J
 97             FORMAT (//5X,'**** PROGRAM EXECUTION STOPPED ****',/
     +     /5X,'INCORRECT BED MASK FOR SAND TRANSPORT ONLY',
     +     /5X,'BED MASK SET TO 0 AT i,j =',2I5,
     +     /5X,'MUST SPECIFY BED MASK VALUE AS 1 (SUSPENDED TRANSPORT)',
     +     /5X,'OR BED MASK VALUE AS -1 (HARD BOTTOM)',
     +     /5X,'PLEASE CORRECT bed_mask FILE AND RESUBMIT'/)
                STOP
              ENDIF
 96         CONTINUE
          ENDIF
 5      CONTINUE
 6      FORMAT (40I3)
        CLOSE (85)
!
!  NO SUSPENDED LOAD TRANSPORT IN COHESIVE ONLY RUN
!
        IF (SEDTYPE.EQ.'MUD ') IBED=0
!
!  VARIABLE BED BULK DENSITY (g/cm**3)
!
!  FOR COMPUTING BED ELEVATION CHANGES AND FOR WASTOX/WASP USE
!
      IF (VARIBULK.EQ.'INCLUDE') THEN
        OPEN (UNIT=85,FILE='bed_bulkden',FORM='FORMATTED')
        DO 9019 I=1,IM
          READ (85,9018) (CBED(I,J),J=1,JM)
 9019   CONTINUE
 9018   FORMAT (20F6.0)
        CLOSE(85)
      ELSE
        DO 9020 J=1,JM
          DO 9020 I=1,IM
!
!  COHESIVE ELEMENTS
!
            IF (IBMSK(I,J).EQ.0) THEN
              CBED(I,J)=DENCOH
            ELSE
              CBED(I,J)=DENNON
            ENDIF
 9020   CONTINUE
      ENDIF
!
!  INPUT SPATIALLY VARIABLE P0 
!
      IF (IP0.EQ.1) THEN
        OPEN (UNIT=86,FILE='p0_init',FORM='FORMATTED')
        DO 7 I=1,IM
          READ (86,8) (PINIT(I,J),J=1,JM)
 7      CONTINUE
 8      FORMAT (20F6.2)
        CLOSE (86)
      ENDIF
!
!  INPUT SPATIALLY VARIABLE A0 AND REXP
!
!  A0 INPUT FOR mg/cm**2 EROSION POTENTIAL
!
      IF (IA0.EQ.1) THEN
        OPEN (UNIT=53,FILE='a0_init',FORM='FORMATTED')
        DO 9000 I=1,IM
          READ (53,9001) (A0(1,I,J),J=1,JM)
!
!  CONVERT FROM mg/cm**2 TO g/cm**2
!
          DO 9000 J=1,JM
            A0(1,I,J)=A0(1,I,J)/1000.
 9000   CONTINUE
 9001   FORMAT (8E10.3)
        DO 9002 J=1,JM
          DO 9002 I=1,IM
            DO 9002 LL=2,LAYMAX
              A0(LL,I,J)=A0(1,I,J)
 9002   CONTINUE
        CLOSE (53)
!
        OPEN (UNIT=53,FILE='exp_init',FORM='FORMATTED')
        DO 9010 I=1,IM
          READ (53,9001) (REXP(I,J),J=1,JM)
 9010   CONTINUE
        CLOSE (53)
      ENDIF
!
      DO 10 J=1,JM
        DO 10 I=1,IM
          IF (RESTAR.EQ.'COLD START') THEN  
          DO 15 K=1,KB
            TAU(I,J,K)=0.0
 15       CONTINUE
          ENDIF                            
          DO 20 LL=1,LAYMAX
!
!  SET CRITICAL SHEAR STRESS AS SPATIALLY CONSTANT OVER EACH LAYER
!
            TCRIT(LL,I,J)=TAUCR(LL)
!
            IF (IBMSK(I,J).EQ.0) THEN
              TSED0(LL,I,J)=TSED0IN(LL)*CBED(I,J)
            ELSE
              TSED0(LL,I,J)=0.0
            ENDIF
!
!  BED LOAD THICKNESS INPUT IN m, CONVERT TO VOLUME (m**3)
!  NOTE: V(SED) = (1. - BPOR) * V(TOTAL)
!        WHERE V(TOTAL) = AREA * BEDTHI
!
           IF (RESTAR.EQ.'COLD START') THEN
            N24CNT=0
            IF (IBMSK(I,J).EQ.0) THEN
              TSED(LL,I,J)=TSED0(LL,I,J)
              BEDTH(1,I,J)=0.0
              BEDTH(2,I,J)=0.0
            ELSE
              IF (IBMSK(I,J).EQ.1) THEN
                TSED(LL,I,J)=0.0
                BEDTH(1,I,J)=(CBED(I,J)/2.65)*BEDTHI*H1(I,J)*H2(I,J)
                BEDTH(2,I,J)=(CBED(I,J)/2.65)*BEDTHI*H1(I,J)*H2(I,J)
              ELSE
                TSED(LL,I,J)=0.0
                BEDTH(1,I,J)=0.0
                BEDTH(2,I,J)=0.0
              ENDIF
            ENDIF
!
!
            TAUMAX(LL,I,J)=0.0
            TAUCUR(LL,I,J)=0.0
            EBMAX(LL,I,J)=0.0
            EBCUR(LL,I,J)=0.0
!
            IF (TSED(LL,I,J).GT.0.0) THEN
              LAYER(LL,I,J)=1
            ELSE
              LAYER(LL,I,J)=0
            ENDIF
            IF (KSED.EQ.1) THEN
              IF (LAYER(LL,I,J).EQ.1) THEN
                PSED1(LL,I,J)=1.0
                PSED2(LL,I,J)=0.0
              ELSE
                PSED1(LL,I,J)=0.0
                PSED2(LL,I,J)=0.0
              ENDIF
            ELSE
              IF (LAYER(LL,I,J).EQ.1) THEN
                IF (IP0.EQ.0) THEN
                  PSED1(LL,I,J)=P0(1)
                  PSED2(LL,I,J)=1.-P0(1)
                ELSE
                  PSED1(LL,I,J)=PINIT(I,J)
                  PSED2(LL,I,J)=1.-PINIT(I,J)
                ENDIF
              ELSE
                PSED1(LL,I,J)=0.0
                PSED2(LL,I,J)=0.0
              ENDIF
            ENDIF
!
            IF (IBMSK(I,J).EQ.1.OR.IBMSK(I,J).LT.0) THEN
              PSED1(LL,I,J)=0.0
              PSED2(LL,I,J)=0.0
            ENDIF
           ENDIF
 20       CONTINUE
!
          WCT1BOT(I,J)=0.0
          WCT2BOT(I,J)=0.0
 10   CONTINUE
!
      DO 26 K=1,KB
        DO 26 J=1,JM
          DO 26 I=1,IM
            WSET1(I,J,K)=0.0
            WSET2(I,J,K)=0.0
 26   CONTINUE
!
!      TCRDP2=1.
       TCRDP2=(WS2/10000.)**2.
!
      WRITE (IUPRT,27)TCRDP2
 27   FORMAT (/5X,'CLASS 2 CRITICAL SHEAR STRESS FOR DEP (dyne/cm**2)',
     +            F10.4)
!
!  WATER-COLUMN SETTLING SPEED FOR COARSE SEDIMENT 
!  AND NON-COHESIVE SUSPENDED SEDIMENTS
!  CONVERT FROM microns/s TO m/s
!
      WMIN=1000000.
      WMAX=-1000000.
!
      DO 30 K=1,KBM1
        DO 30 J=1,JM
          DO 30 I=1,IM
            WSET2(I,J,K)=WS2/1000000.
!
!  MODIFIED BY C.K.Z. ON 4/4/97
!
!  SET MAX. SETTLING SPEED FOR CLASS 2 EQUAL TO 0.25*CFL(VERT)
!
            IF (FSM(I,J).GT.0.0.AND.K.GE.2) THEN
!NG04272010 W&D              CFLZ1=WSET2(I,J,K)*DTI/(H(I,J)*DZ(K-1))
!NG04272010 W&D              CFLZ2=WSET2(I,J,K)*DTI/(H(I,J)*DZ(K))
              CFLZ1=WSET2(I,J,K)*DTI/(D(I,J)*DZ(K-1))
              CFLZ2=WSET2(I,J,K)*DTI/(D(I,J)*DZ(K))
!
              WST1=WSET2(I,J,K)
              WST2=WSET2(I,J,K)
!
              IF (CFLZ1.GT.0.1) THEN
!NG04272010 W&D                WST1=0.25*H(I,J)*DZ(K-1)/DTI
                WST1=0.25*D(I,J)*DZ(K-1)/DTI
              ENDIF
!
              IF (CFLZ2.GT.0.1) THEN
!NG04272010 W&D                WST2=0.25*H(I,J)*DZ(K)/DTI
                WST2=0.25*D(I,J)*DZ(K)/DTI
              ENDIF
!
              WSET2(I,J,K)=AMIN1(WST1,WST2)
!
              IF (WSET2(I,J,K).LT.WMIN) THEN
                IMIN=I
                JMIN=J
                KMIN=K
                WMIN=WSET2(I,J,K)
              ENDIF
!
              IF (WSET2(I,J,K).GT.WMAX) THEN
                IMAX=I
                JMAX=J
                KMAX=K
                WMAX=WSET2(I,J,K)
              ENDIF
            ENDIF
 30   CONTINUE
!
!  SET WSET2 EQUAL TO MINIMUM SETTLING SPEED
!
      DO 301 K=1,KBM1
        DO 301 J=1,JM
          DO 301 I=1,IM
            WSET2(I,J,K)=WMIN
 301  CONTINUE
        WMIN=1000000.*WMIN
        WMAX=1000000.*WMAX
        WRITE (IUPRT,31)IMIN,JMIN,KMIN,WMIN,H(IMIN,JMIN)
 31     FORMAT (/5X,'**** MIN WSET2 & H: ',3I5,F10.1,F10.2)
        WRITE (IUPRT,32)IMAX,JMAX,KMAX,WMAX,H(IMAX,JMAX)
 32     FORMAT (/5X,'**** MAX WSET2 & H: ',3I5,F10.1,F10.2)
!
!  WS2 FOR COMP. 2 (COARSE), INPUT IN microns/s
!  FOR SEDIMENT-WATER INTERFACE FLUX
!
!  FOR SPEED IMPROVEMENT
!  MULT. WS2 BY XNSBED
!
      WS2=XNSBED*WS2*DTI/10000.
!
!  CALC. SETTLING SPEEDS FOR FLOCCULATION
!
!  SEDZL SETTLING SPEEDS:  BASED ON REANALYSIS OF LICK DATA
! 
!pl 10/19/98      DO 40 II=1,62500
      DO 40 II=1,6250000
        CT=FLOAT(II)/100.
!
!  SETTLING SPEED IN cm/s
!
        WS1(II)=ADEP*(CT**DEPEXP)/10000.
 40   CONTINUE
!
!  CALC. PROBABILITY OF DEPOSITION
!
!  INCREMENT EVERY 0.01 dyne/cm**2
! 
      DO 45 II=1,20000 
        TAUBOT=FLOAT(II-1)*0.01
!
!  KRONE FORMULATION
!
        IF (PDEPFORM.EQ.'KRONE') THEN
          IF (TAUBOT.LT.TCRDEP) THEN
            PDEP(II)=1.-TAUBOT/TCRDEP
          ELSE
            PDEP(II)=0.0
          ENDIF
        ELSE
!
!  PARTHENIADES FORMULATION
!
!  ASSUME:  SIGY = 0.49
!           (TAUB - 1)50 = 0.84
!
          IF (TAUBOT.LE.TCRDEP) THEN
            PDEP(II)=1.
          ELSE
!
!  CALC. Y PARAMETER
!
          YPARM=2.04*ALOG10(0.25*((TAUBOT/TCRDEP)-1.)*EXP(1.27*TCRDEP))
!
            IF (YPARM.LT.0.0) THEN
              YPARM=-YPARM
              IFLAG=1
            ELSE
              IFLAG=0
            ENDIF
!
!  CALC. Z(Y)
!
            ZY=EXP(-YPARM*YPARM/2.)/SQRT(6.283185)
!
!  APPROXIMATE EVALUATION OF CEQ INTEGRAL
!
            TY=1./(1.+0.3327*YPARM)
!
            CEQ=1.-ZY*(0.4362*TY-0.1202*(TY**2.)+0.9373*(TY**3.))
!
            IF (IFLAG.EQ.1) CEQ=1.-CEQ
!
!  CALC. PROBABILITY OF DEPOSITION
!
            PDEP(II)=1.-CEQ
          ENDIF
        ENDIF
 45   CONTINUE
!
!  CALC. TAUD2
!
      TAU1=0.0
      DO 60 I=1,2000
        TAU1=TAU1+.1
        DO 70 LL=1,LAYMAX 
          RAT1=(TAU1-TAUCR(LL))/TAUCR(LL)
          IF (RAT1.LT.0.0) RAT1=0.0
          TAUD2(LL,I)=RAT1**RESEXP
 70     CONTINUE
 60   CONTINUE
!
!********************************************************************
!
!  COHESIVE SEDIMENT CONSTANTS
!
!  FOR SPEED IMPROVEMENT
!  MULT. S1 BY XNSBED
!
      S1=XNSBED*DTI/3600.
!
!  CONVERT A0IN FROM mg/cm**2 TO g/cm**2
!
      A0IN=A0IN/1000.
!
      DO 90 LL=1,LAYMAX
        FTIME2(LL)=FTIME(LL)**EXPM
!
        DO 90 J=1,JM
          DO 90 I=1,IM
            IF (IBMSK(I,J).EQ.0) THEN
              IF (IA0.EQ.0) THEN
                A0(LL,I,J)=A0IN/FTIME2(LL)
              ELSE
                A0(LL,I,J)=A0(LL,I,J)/FTIME2(LL)
              ENDIF
              ASED(LL,I,J)=S1*A0(LL,I,J)
            ENDIF
 90   CONTINUE
!
!
!  SUSPENDED LOAD CONSTANTS
!
!  REF:  VAN RIJN, 1984, J. HYDR. ENGR., p. 1431-1456
!        VAN RIJN, 1984, J. HYDR. ENGR., p. 1613-1641
!
!  MKS SYSTEM USED
!  ASSUMED THAT SEDIMENT DENSITY IS 2.65 g/cm**3 AND KINEMATIC
!  VISCOSITY IS 0.011 cm**2/s
!
      IF (IBED.GE.1) THEN
!
!************************************************
!
!  INPUT VARIABLE D50
!
        OPEN (UNIT=85,FILE='bed_d50',FORM='FORMATTED')
        REWIND(85)
        DO 104 I=1,IM
          READ (85,103) (D50VAR(I,J),J=1,JM)
 104    CONTINUE
 103    FORMAT (20F6.0)
        CLOSE (85)
!
!  CONVERT FROM MICRONS TO METERS
!
        DO 102 J=1,JM
          DO 102 I=1,IM
            D50VAR(I,J)=D50VAR(I,J)/1000000.
 102    CONTINUE
!
!  COMPUTE EFFECTIVE D50 OF CLASS 2 FROM SETTLING SPEED (IN cm)
!
        D50C2=0.05
        ACON1=133636.
        ACON2=9.091*WSET2(1,1,1)*100.
!
!  USE NEWTON'S METHOD TO SOLVE FOR D50, 100 ITERATIONS
!
        DO 105 N=1,20 
          FD50=SQRT(1.+ACON1*(D50C2**3.))-1.-ACON2*D50C2
          FPD50=((1.5*ACON1*D50C2*D50C2)/SQRT(1.+ACON1*(D50C2**3.)))
     +             -ACON2
          D50C2=D50C2-FD50/FPD50
 105    CONTINUE
        D50C2=D50C2*10000.
        WRITE (IUPRT,106)D50C2
 106    FORMAT (/5X,'NON-COHESIVE EFFECTIVE D50 (microns) =',F8.0
     +          /5X,'(BASED UPON SETTLING SPEED)'/)
        D50C2=D50C2/1000000.
!
!  VARIABLE D50
!
        DO 117 J=1,JM
          DO 117 I=1,IM
            IF (IBMSK(I,J).EQ.1.AND.FSM(I,J).GT.0.0) THEN
!
!
!  CALC. DSTAR (EQ. 1)
!
              DSTAR=23731.*D50VAR(I,J)
!
!  CALC. THETA CRITICAL (FIG. 1)
!
              IF (DSTAR.LE.4.) THEN
                THCR=0.24/DSTAR
                GOTO 110
              ENDIF
!
              IF (DSTAR.LE.10.) THEN
                THCR=0.14/(DSTAR**0.64)
                GOTO 110
              ENDIF
!
              IF (DSTAR.LE.20) THEN
                THCR=0.04/(DSTAR**0.10)
                GOTO 110
              ENDIF
!
              IF (DSTAR.LE.150) THEN 
                THCR=0.013*(DSTAR**0.29)
              ELSE
                THCR=0.055
              ENDIF
!
!  CALC. USTAR CRITICAL (FIG. 1)
!
 110          UCRVAR(I,J)=SQRT(16.17*D50VAR(I,J)*THCR)
            ENDIF
 117    CONTINUE
!
        CSED(8)=DSTAR
!
!  CLASS 2
!
!
!  CALC. DSTAR (EQ. 1)
!
        DSTAR2=23731.*D50C2
!
!  CALC. THETA CRITICAL (FIG. 1)
!
        IF (DSTAR2.LE.4.) THEN
          THCR=0.24/DSTAR2
          GOTO 111
        ENDIF
!
        IF (DSTAR2.LE.10.) THEN
          THCR=0.14/(DSTAR2**0.64)
          GOTO 111
        ENDIF
!
        IF (DSTAR2.LE.20) THEN
          THCR=0.04/(DSTAR2**0.10)
          GOTO 111
        ENDIF
!
        IF (DSTAR2.LE.150) THEN 
          THCR=0.013*(DSTAR2**0.29)
        ELSE
          THCR=0.055
        ENDIF
!
!  CALC. USTAR CRITICAL (FIG. 1)
!
 111    UCR2=SQRT(16.17*D50C2*THCR)
!
!  SUSPENDED LOAD CONSTANTS
!
!  NOTE:  CSED(1) AND CSED(2) ARE USED IN BEDLOAD CALC.
!
        CSED(1)=SQRT(9.8)
!
        CSED(3)=0.015*D50C2/(DSTAR2**0.3)
!
!  INPUT NON-COHESIVE MASS FRACTIONS
!
       IF (RESTAR.EQ.'COLD START') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (UNIT=85,FILE='bed_frac.mud',FORM='FORMATTED')
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (UNIT=86,FILE='bed_frac.sand',FORM='FORMATTED')
        ENDIF
!
        DO 120 K=1,2
          DO 120 I=1,IM
            DO 125 J=1,JM
              FPBED(K,I,J)=0.0
 125        CONTINUE
!
            IF((SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH').AND.K.EQ.1) THEN
              READ (85,130) (FPBED(K,I,J),J=1,JM)
            ENDIF
!
            IF((SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH').AND.K.EQ.2) THEN
              READ (86,130) (FPBED(K,I,J),J=1,JM)
            ENDIF
!
              DO 140 J=1,JM
                IF (IBMSK(I,J).EQ.1.AND.FSM(I,J).GT.0.0) THEN
                  FALAY(K,I,J)=0.0
!
                  ACTLAY(2,I,J)=0.0
!
                  CARMOR(I,J)=(2.*D50VAR(I,J)*(CBED(I,J)/2.65)*
     +              H1(I,J)*H2(I,J))/(10000.*UCRVAR(I,J)*UCRVAR(I,J))
                ENDIF
 140          CONTINUE
 120    CONTINUE
 130    FORMAT (20F6.3)
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') CLOSE (85)
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') CLOSE (86)
       ENDIF
      ENDIF
      RETURN
      END
