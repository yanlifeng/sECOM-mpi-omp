      SUBROUTINE REORDR
!
!  CHANGES AGE OF LAYERS AT END OF EACH DAY AND RESETS TAUMAX
!
!  REVISION DATE :  JUNE 8, 1995
!
!****************************************************************
!
      INCLUDE 'comdeck'
!
!  24 HOUR REORDERING
!
      LMAXM=LAYMAX-1
      CRAT=FTIME2(LMAXM)/FTIME2(LAYMAX)
      DO 10 J=2,JM-1
        DO 10 I=2,IM-1
!
          IF (IBMSK(I,J).LT.0.OR.IBMSK(I,J).EQ.1.OR.
     +                              FSM(I,J).EQ.0.0) GOTO 10
!
!  ADD MATERIAL TO 7 DAY OLD LAYER
!
          IF (LAYER(LMAXM,I,J).EQ.1) THEN
            IF (KSED.GT.1) THEN
              IF (LAYER(LAYMAX,I,J).EQ.1) THEN        
                PSED1(LAYMAX,I,J)=(PSED1(LMAXM,I,J)*TSED(LMAXM,I,J)
     +          +PSED1(LAYMAX,I,J)*TSED(LAYMAX,I,J))/
     +          (TSED(LAYMAX,I,J)+TSED(LMAXM,I,J))
                PSED2(LAYMAX,I,J)=1.-PSED1(LAYMAX,I,J)
              ELSE
                PSED1(LAYMAX,I,J)=PSED1(LMAXM,I,J)
                PSED2(LAYMAX,I,J)=PSED2(LMAXM,I,J)
              ENDIF
            ELSE
              PSED1(LAYMAX,I,J)=1.
              PSED2(LAYMAX,I,J)=0.
            ENDIF
            LAYER(LAYMAX,I,J)=1
            TSED(LAYMAX,I,J)=TSED(LAYMAX,I,J)+TSED(LMAXM,I,J)
!
!PL ADD THE FOLLOWING LINES PRAVI 09/08/99
!
            EBMAX(LAYMAX,I,J)=CRAT*EBMAX(LMAXM,I,J)
            EBTOT(LAYMAX,I,J)=CRAT*EBTOT(LMAXM,I,J)
            EBCUR(LAYMAX,I,J)=CRAT*EBCUR(LMAXM,I,J)
            TAUMAX(LAYMAX,I,J)=TAUMAX(LMAXM,I,J)
            TAUCUR(LAYMAX,I,J)=TAUCUR(LMAXM,I,J)
!PL 09/08/99            EBTOT(LAYMAX,I,J)=EBTOT(LAYMAX,I,J)-TSED(LMAXM,I,J)
!PL 09/08/99            IF (EBTOT(LAYMAX,I,J).LT.0.0) EBTOT(LAYMAX,I,J)=0.0
!
            PSED1(LMAXM,I,J)=0.0
            PSED2(LMAXM,I,J)=0.0
            TSED(LMAXM,I,J)=0.0
            LAYER(LMAXM,I,J)=0
            TAUMAX(LMAXM,I,J)=0.0
            EBTOT(LMAXM,I,J)=0.0
            TAUCUR(LMAXM,I,J)=0.0
            EBMAX(LMAXM,I,J)=0.0
            EBCUR(LMAXM,I,J)=0.0
          ENDIF
!
!  ADD MATERIAL TO YOUNGER LAYERS 
!
         DO 16 LL=LAYMAX-1,2,-1
          IF (LAYER(LL-1,I,J).EQ.1) THEN
            IF (KSED.GT.1) THEN
              IF (LAYER(LL,I,J).EQ.1) THEN        
                PSED1(LL,I,J)=(PSED1(LL-1,I,J)
     +          *TSED(LL-1,I,J)+PSED1(LL,I,J)
     +          *TSED(LL,I,J))/
     +          (TSED(LL,I,J)+TSED(LL-1,I,J))
                PSED2(LL,I,J)=1.-PSED1(LL,I,J)
              ELSE
                PSED1(LL,I,J)=PSED1(LL-1,I,J)
                PSED2(LL,I,J)=PSED2(LL-1,I,J)
              ENDIF
            ELSE
              PSED1(LL,I,J)=1.
              PSED2(LL,I,J)=0.
            ENDIF
            LAYER(LL,I,J)=1
            TSED(LL,I,J)=TSED(LL,I,J)+TSED(LL-1,I,J)
!
            CRAT=FTIME2(LL-1)/FTIME2(LL)
!
            EBMAX(LL,I,J)=CRAT*EBMAX(LL-1,I,J)
            EBTOT(LL,I,J)=CRAT*EBTOT(LL-1,I,J)
            EBCUR(LL,I,J)=CRAT*EBCUR(LL-1,I,J)
            TAUMAX(LL,I,J)=TAUMAX(LL-1,I,J)
            TAUCUR(LL,I,J)=TAUCUR(LL-1,I,J)
!
            PSED1(LL-1,I,J)=0.0
            PSED2(LL-1,I,J)=0.0
            TSED(LL-1,I,J)=0.0
            LAYER(LL-1,I,J)=0
            TAUMAX(LL-1,I,J)=0.0
            EBTOT(LL-1,I,J)=0.0
            TAUCUR(LL-1,I,J)=0.0
            EBMAX(LL-1,I,J)=0.0
            EBCUR(LL-1,I,J)=0.0
          ENDIF
 16      CONTINUE
!
          PSED1(1,I,J)=0.0
          PSED2(1,I,J)=0.0
          TSED(1,I,J)=0.0
          LAYER(1,I,J)=0
          TAUMAX(1,I,J)=0.0
          EBTOT(1,I,J)=0.0
          TAUCUR(1,I,J)=0.0
          EBMAX(1,I,J)=0.0
          EBCUR(1,I,J)=0.0
 10   CONTINUE
!
      N24CNT=0
!
      RETURN
      END
