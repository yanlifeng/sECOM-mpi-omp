      SUBROUTINE CHEMBED(CHEMDRAT1,CHEMDRAT2)
!
!  CHEM BED MODEL
!
!  NOTE:  MUST USE OPPOSITE SIGN OF FLUXES FOR DEPOSITION TO BED
!
!  REVISION DATE:  OCTOBER 9, 1996
!
!************************************************************
!
      INCLUDE 'comdeck'
!
      REAL*8 CHEMM,SEDM,SUMTEM,TLAY1,TOTLAY,DM,FRAC,SEDTEMP
!
      CDEC1=EXP(-CHEMDRAT1*DTI*XNSBED)
      CDEC2=EXP(-CHEMDRAT2*DTI*XNSBED)
!
      DO 10 J=2,JMM1
        DO 10 I=2,IMM1
          IF (FSM(I,J).GT.0.0) THEN
!
            CHEMM=-XNSBED*(CHEMBOT1(I,J)+CHEMBOT2(I,J))
            SEDM=-XNSBED*(WCT1BOT(I,J)+WCT2BOT(I,J))
!
            IF (SEDM.EQ.0.0) GOTO 10
!
!  DEPOSITION 
!
            IF (SEDM.GT.0.0) THEN
!
!  ADD MASS OF CHEM TO LAYER 1 (ug CHEM/cm**2)
!
              CHEMMASS(1,I,J)=CHEMMASS(1,I,J)+CHEMM
!
!  ADD MASS OF SOLIDS TO LAYER 1 (g SOLIDS/cm**2)
!
              SEDMASS(1,I,J)=SEDMASS(1,I,J)+SEDM
!
!  CALC. THICKNESS OF LAYER 1 (cm)
!
              TLAY1=SEDMASS(1,I,J)/CBED(I,J)
!
!  MIX ACTIVE LAYER
!
              IF (NACTLAY.GT.1) THEN
!
!  TOTAL CHEM MASS IN ACTIVE LAYER (ug CHEM/cm**2)
!
                SUMTEM=0.0
                DO 20 N=1,NACTLAY
                  SUMTEM=SUMTEM+CHEMMASS(N,I,J)
 20             CONTINUE
!
!  TOTAL THICKNESS OF ACTIVE LAYER (cm)
!
                TOTLAY=TLAY1+FLOAT(NACTLAY-1)*CHEMTHIK
!
!  THICKNESS-WEIGHTED DISTRIBUTION OF MASS CHEM IN ACTIVE LAYER
!
                CHEMMASS(1,I,J)=TLAY1*SUMTEM/TOTLAY
!
                DO 30 N=2,NACTLAY
                  CHEMMASS(N,I,J)=CHEMTHIK*SUMTEM/TOTLAY
 30             CONTINUE
              ENDIF
!
!  CHECK IF LAYER 1 THICKNESS IS GREATER THAN BED LAYER THICKNESS
!
!  REORDER ALL LAYERS
!
              IF (SEDMASS(1,I,J).GT.(CBED(I,J)*CHEMTHIK)) THEN
!
                IF (NCHEMLAY.GT.2) THEN
                  CHEMMASS(NCHEMLAY,I,J)=CHEMMASS(NCHEMLAY,I,J)+
     +                                   CHEMMASS(NCHEMLAY-1,I,J)
                  SEDMASS(NCHEMLAY,I,J)=SEDMASS(NCHEMLAY,I,J)+
     +                                   SEDMASS(NCHEMLAY-1,I,J)
!
                  IF (NCHEMLAY.GT.3) THEN
                    DO 50 N=NCHEMLAY-1,3,-1
                      CHEMMASS(N,I,J)=CHEMMASS(N-1,I,J)
                      SEDMASS(N,I,J)=SEDMASS(N-1,I,J)
 50                 CONTINUE
                  ENDIF
!
!  CALC. FRACTION OF SEDIMENT LEFT IN LAYER 1
!
                  DM=SEDMASS(1,I,J)-(CBED(I,J)*CHEMTHIK)
                  FRAC=DM/SEDMASS(1,I,J)
!
!  RESET MASSES IN LAYERS 1 & 2
!
                  SEDMASS(1,I,J)=DM
                  SEDMASS(2,I,J)=CBED(I,J)*CHEMTHIK
!
                  CHEMMASS(2,I,J)=(1.-FRAC)*CHEMMASS(1,I,J)
                  CHEMMASS(1,I,J)=FRAC*CHEMMASS(1,I,J) 
                ELSE
!
!  TWO LAYER MODEL
!
                  DM=SEDMASS(1,I,J)-(CBED(I,J)*CHEMTHIK)
                  FRAC=DM/SEDMASS(1,I,J)
!
                  SEDMASS(1,I,J)=DM
                  SEDMASS(2,I,J)=SEDMASS(2,I,J)+CBED(I,J)*CHEMTHIK
!
                  CHEMMASS(2,I,J)=CHEMMASS(2,I,J)+(1.-FRAC)*
     +                   CHEMMASS(1,I,J)
                  CHEMMASS(1,I,J)=FRAC*CHEMMASS(1,I,J) 
                ENDIF
              ENDIF
!
!  EROSION
!
            ELSE
!
!  REMOVE MASS OF SOLIDS FROM LAYER 1 (g SOLIDS/cm**2)
!
              SEDTEMP=SEDMASS(1,I,J)+SEDM
!
!  LAYER 1 ERODED
! 
              IF (SEDTEMP.LE.0.0) THEN
!
!  REORDER LAYERS
!
!  FIND DEEPEST LAYER WITH MASS
!
                DO 60 N=2,NCHEMLAY
                  IF (SEDMASS(N,I,J).LE.0.0) THEN
                    NMAX=N-1
                    GOTO 70
                  ENDIF
 60             CONTINUE
                NMAX=NCHEMLAY
!
 70             SEDMASS(1,I,J)=SEDMASS(2,I,J)+SEDTEMP
                CHEMMASS(1,I,J)=CHEMMASS(2,I,J)+SEDTEMP*CBEDCHEM(2,I,J)
!
                IF (NMAX.GT.2) THEN
                  DO 80 N=2,NMAX-1
                    SEDMASS(N,I,J)=SEDMASS(N+1,I,J)
                    CHEMMASS(N,I,J)=CHEMMASS(N+1,I,J)
 80               CONTINUE
                ENDIF
!
                SEDMASS(NMAX,I,J)=0.0
                CHEMMASS(NMAX,I,J)=0.0
!
!  LAYER 1 NOT ERODED
!
              ELSE
                SEDMASS(1,I,J)=SEDTEMP
                CHEMMASS(1,I,J)=CHEMMASS(1,I,J)+CHEMM
              ENDIF
            ENDIF
!
!  RE-CALC. CHEM CONC. (ug CHEM/g SOLIDS)
!
!  INCLUDE FIRST-ORDER DECAY RATE
!
            IF (IBMSK(I,J).EQ.0) THEN
              CDEC=PSED1(LAYMAX,I,J)*CDEC1+PSED2(LAYMAX,I,J)*CDEC2
            ELSE
              CDEC=FPBED(1,I,J)*CDEC1+FPBED(2,I,J)*CDEC2
            ENDIF
!
            DO 90 N=1,NCHEMLAY
              IF (SEDMASS(N,I,J).GT.0.0) THEN
!hli(04/02/01)
                if(CHEMMASS(N,I,J).le.0.0) then
                   CHEMMASS(N,I,J)=0.0
                endif
!hli
                CBEDCHEM(N,I,J)=CDEC*CHEMMASS(N,I,J)/SEDMASS(N,I,J)
              ELSE
                CBEDCHEM(N,I,J)=0.0
              ENDIF
 90         CONTINUE
          ENDIF
 10   CONTINUE
!
      RETURN
      END
