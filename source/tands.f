      SUBROUTINE TANDS
!     VERSION(06/28/91)
!
      INCLUDE 'comdeck'
      SAVE
!
      DIMENSION ZM(KBM1),TI(KB),SI(KB)
      DIMENSION TA(KSL),SA(KSL),TS(IM,JM,KSL),SS(IM,JM,KSL)
!
!  CHECK THAT NO DEPTH IS DEEPER THAN BOTTOM STANDARD LEVEL
!
      DO 9000 J=2,JMM1
        DO 9000 I=2,IMM1
          IF (H(I,J).GT.0.0) THEN
            IF (H(I,J).GT.-DPTHSL(KSL)) THEN
              WRITE (IUPRT,9001)I,J,H(I,J),DPTHSL(KSL)
 9001         FORMAT (/5X,'****** PROGRAM EXECUTION STOPPED *****',
     +  /7X,'SPECIFIED DEPTH IN model_grid AT (i,j) =',2I4,
     +  /7X,'IS ',F10.2,3X,'WHICH IS DEEPER THAN THE LOWEST STANDARD',
     +  /7X,'LEVEL, WHICH IS ',F10.2,
     +  /7X,'PLEASE CORRECT AND RESUBMIT')
!
             STOP
           ENDIF
         ENDIF
 9000 CONTINUE
!
!-------- SPREAD STANDARD LEVEL TEMPERATURE AND SALINITY ---------------
      IF(OPTTSI(1:5).EQ.'FIXED') THEN
       DO 120 K=1,KSL
       DO 120 J=1,JM
       DO 120 I=1,IM
       TS(I,J,K)=TSI(K)*FSM(I,J)
       SS(I,J,K)=SSI(K)*FSM(I,J)
  120  CONTINUE
!
      ELSEIF(OPTTSI(1:5).EQ.'SIGMA') THEN
       OPEN (IUTAS,FILE='init_tands',FORM='formatted',status='OLD')
       DO 121 N=1,(IM*JM)
       READ(IUTAS,132,END=122)
     .            I,J,(T(I,J,K),K=1,KBM1),(S(I,J,K),K=1,KBM1)
  121  CONTINUE
  122  CLOSE(IUTAS)
       GOTO 360
!
      ELSE
!-------- INPUT TEMPERATURE AND SALINITY DATA FROM FILE ----------------
       OPEN (IUTAS,FILE='init_tands',FORM='formatted',status='OLD')
       DO 130 N=1,(IM*JM)
       READ(IUTAS,132,END=1310)
     .            I,J,(TS(I,J,K),K=1,KSL),(SS(I,J,K),K=1,KSL)
  130  CONTINUE
  131  CONTINUE
       CLOSE(IUTAS)
  132  FORMAT (2I5,100F5.0)
      ENDIF
!
      DO 220 K=1,KSL
      TA(K)=0.0
      SA(K)=0.0
!
      COUNT=0.0
      DO 210 J=1,JM
      DO 210 I=1,IM
!
            IF (K.EQ.1 )THEN
             IF(FSM(I,J).GT.0.0)THEN
              COUNT = COUNT + FSM(I,J)*ART(I,J)
              TA(K) = TA(K) + ART(I,J)*TS(I,J,K) * FSM(I,J)
              SA(K) = SA(K) + ART(I,J)*SS(I,J,K) * FSM(I,J)
             ENDIF
            ELSE
             IF(H(I,J).GT.-DPTHSL(K-1)) THEN 
              COUNT = COUNT + FSM(I,J)*ART(I,J)
              TA(K) = TA(K) + ART(I,J)*TS(I,J,K) * FSM(I,J)
              SA(K) = SA(K) + ART(I,J)*SS(I,J,K) * FSM(I,J)
             ENDIF
            ENDIF
!
  210 CONTINUE
      IF (COUNT.LT.1.0) GO TO 220
      TA(K)=TA(K)/COUNT
      SA(K)=SA(K)/COUNT
!
  220 CONTINUE
      DO 250 J=1,JM
      DO 250 I=1,IM
      IF (FSM(I,J).EQ.0.0) GO TO 250
      DO 230 K=1,KBM1
      ZM(K)=ZZ(K)*H(I,J)
  230 CONTINUE
      CALL SINTER (DPTHSL,TA,ZM,TI,KSL,KBM1)
      CALL SINTER (DPTHSL,SA,ZM,SI,KSL,KBM1)
!
      DO 240 K=1,KBM1
      T(I,J,K)=TI(K)
      S(I,J,K)=SI(K)
!
  240 CONTINUE
  250 CONTINUE
!
      CALL DENS
!
      DO 280 K=1,KBM1
      DO 280 J=1,JM
      DO 280 I=1,IM
      TMEAN(I,J,K)=T(I,J,K)*FSM(I,J)
      SMEAN(I,J,K)=S(I,J,K)*FSM(I,J)
      RMEAN(I,J,K)=RHO(I,J,K)*FSM(I,J)
  280 CONTINUE
!
      DO 350 J=1,JM
      DO 350 I=1,IM
      IF (FSM(I,J).EQ.0.0) GO TO 350
      DO 320 K=1,KSL
      TA(K)=TS(I,J,K)
      SA(K)=SS(I,J,K)
!
  320 CONTINUE
      DO 330 K=1,KBM1
      ZM(K)=ZZ(K)*H(I,J)
  330 CONTINUE
      CALL SINTER (DPTHSL,TA,ZM,TI,KSL,KBM1)
      CALL SINTER (DPTHSL,SA,ZM,SI,KSL,KBM1)
!
      DO 340 K=1,KBM1
      T(I,J,K)=TI(K)*FSM(I,J)
      S(I,J,K)=SI(K)*FSM(I,J)
!
  340 CONTINUE
  350 CONTINUE
  360 CALL DENS
!
! Introduce SGW's scheme to get density (NKIM, 12/l7/02)
!
      DO K=1,KBM1
         DO J=1,JM
            DO I=1,IM
               RMEAN(I,J,K)=RHO(I,J,K)*FSM(I,J)
            ENDDO
         ENDDO
      ENDDO
!end
      DO 430 K=1,KBM1
      DO 430 J=1,JM
      DO 430 I=1,IM
      TB(I,J,K)=T(I,J,K)
      SB(I,J,K)=S(I,J,K)
!
  430 CONTINUE
!
      DO 440 K=1,KB
      DO 440 J=1,JM
      DO 440 I=1,IM
      A(I,J,K)=0.0
      C(I,J,K)=0.0
      VH(I,J,K)=0.0
      VHP(I,J,K)=0.0
      PROD(I,J,K)=0.0
      DTEF(I,J,K)=0.0
  440 CONTINUE
!
      DO 450 K=1,KSL
      DO 450 J=1,JM
      DO 450 I=1,IM
      TS(I,J,K)=0.0
      SS(I,J,K)=0.0
!
  450 CONTINUE
!
!  SET INITIAL SETTLING SPEEDS FOR SEDIMENT TRANSPORT
!
!  SETTLING SPEED IN m/s
!
      DO 201 K=1,KBM1
        DO 201 J=2,JMM1
          DO 201 I=2,IMM1
            IF (D(I,J).GT.0.0) THEN
!
!  REMOVE DEPENDENCE ON TAU INITIALLY
!
              IC=NINT(1000000.*CSED1(I,J,K)*100.)
              IF (IC.LT.1) IC=1
!jah          IF (IC.GT.6250000) IC=6250000
              IF (IC.GT.625000) IC=625000
! 
!   NO DEPOSITION WHEN BED THICKNESS > BATHYMETRIC DEPTH
!   PRAVI 08/12/99 
!
!                IF (CSED1(I,J,K).GE.0.0) THEN
              IF (CSED1(I,J,K).GT.0.0 
     +           .AND. SEDTHK(I,J).LT.H(I,J)) THEN
                WSET1(I,J,K)=WS1(IC)/100.
              ELSE
                WSET1(I,J,K)=0.0
              ENDIF
            ENDIF
 201  CONTINUE
!
!
      RETURN
1310    WRITE (iuprt,*) 'check init_tands file and resubmit the run'
        stop
      END
