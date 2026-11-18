      RECURSIVE SUBROUTINE HALO43(AU43,AV43,E43,T43,N43)
      INCLUDE 'comdeck'
      INTEGER E43,T43,N43,G43,NG43,LO43,HI43
      REAL AU43(IM,JM),AV43(IM,JM)
      REAL, SAVE :: UROW43(IM,2),VROW43(IM,2)
      NG43=N43/2
      G43=1+T43/NG43
      LO43=1
      HI43=JM/2
      IF(G43.EQ.2) THEN
       LO43=HI43+1
       HI43=JM
      ENDIF
      IF(MOD(T43,NG43).EQ.0) THEN
       IF(HI43.LT.JM) CALL UFACE43(HI43+1,UROW43(1,G43))
       IF(LO43.GT.1) CALL VFACE43(LO43-1,VROW43(1,G43))
      ENDIF
      CALL WORK43(AU43,AV43,LO43,HI43,UROW43(1,G43),VROW43(1,G43),T43,NG43,E43)
!$OMP SINGLE
      CALL EDGE43(AU43,AV43)
!$OMP END SINGLE
      END
      RECURSIVE SUBROUTINE WORK43(ADVUA,ADVVA,LO43,HI43,YU43,YV43,T43,NG43,E43)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM)
      INTEGER LO43,HI43,T43,NG43,E43,L43,H43
      REAL YU43(IM),YV43(IM)
      CALL RANGE43(LO43,HI43,MOD(T43,NG43),NG43,L43,H43)
      DO J=L43,H43
       DO I=1,IM
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
      ADVUA(I,J)=0.0
      TPS(I,J)=0.0
       ENDDO
       IF(J.GE.2) THEN
        IF(ADVECT.EQ.'NON-LINEAR') THEN
         IF(J.LE.JMM1) THEN
          DO I=2,IMM1
      FLUXUA(I,J)=.125*((D(I+1,J)+D(I,J))*UA(I+1,J)
     .                 +(D(I,J)+D(I-1,J))*UA(I,J))
     .                  *(UA(I+1,J)+UA(I,J))
          ENDDO
         ENDIF
         DO I=2,IMM1
      FLUXVA(I,J)=.125*((D(I,J)+D(I,J-1))*VA(I,J)
     .                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     .                 *(UA(I,J)+UA(I,J-1))
         ENDDO
        ENDIF
        IF(J.LE.JMM1) THEN
         DO I=2,IMM1
      FLUXUA(I,J)=FLUXUA(I,J)
     .         -D(I,J)*2.*AAM2D(I,J)*(UAB(I+1,J)-UAB(I,J))/H1(I,J)
         ENDDO
        ENDIF
        DO I=2,IMM1
      TPS(I,J)=.25*(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     .            *(AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1))
     .                *((UAB(I,J)-UAB(I,J-1))
     .                /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .                 +(VAB(I,J)-VAB(I-1,J))
     .               /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)) )
!-------- NOTE TPS CARRIED OVER TO V EQUATION !!! ----------------------
      FLUXUA(I,J)=FLUXUA(I,J)*H2(I,J)
      FLUXVA(I,J)=(FLUXVA(I,J)-TPS(I,J)*DUM(I,J)*DUM(I,J-1)
     .            *DUMWAD(I,J)*DUMWAD(I,J-1)*DVMWAD(I,J)*DVMWAD(I-1,J) ! Not in original NYHOPS W&D
     .                                 *DVM(I,J)*DVM(I-1,J))! THIN DAMS NG08192008
     .            *.25*(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1))
        ENDDO
       ENDIF
      ENDDO
!$OMP FLUSH
      CALL LOCALSYNC43(T43,NG43,3*(E43-1)+1)
!$OMP FLUSH
      CALL RANGE43(MAX(2,LO43),MIN(JMM1,HI43),MOD(T43,NG43),NG43,L43,H43)
      DO J=L43,H43
       IF(J.EQ.HI43.AND.HI43.LT.JM) THEN
        DO I=3,IMM1
      ADVUA(I,J)=FLUXUA(I,J)-FLUXUA(I-1,J)
     .           +YU43(I)-FLUXVA(I,J)
        ENDDO
       ELSE
        DO I=3,IMM1
      ADVUA(I,J)=FLUXUA(I,J)-FLUXUA(I-1,J)
     .           +FLUXVA(I,J+1)-FLUXVA(I,J)
        ENDDO
       ENDIF
      ENDDO
!$OMP FLUSH
      CALL LOCALSYNC43(T43,NG43,3*(E43-1)+2)
!$OMP FLUSH
      CALL RANGE43(LO43,HI43,MOD(T43,NG43),NG43,L43,H43)
      DO J=L43,H43
       DO I=1,IM
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
      ADVVA(I,J)=0.0
       ENDDO
       IF(J.GE.2.AND.J.LE.JMM1) THEN
        IF(ADVECT.EQ.'NON-LINEAR') THEN
         IF(J.GE.3) THEN
          DO I=2,IM
      FLUXUA(I,J)=.125*((D(I,J)+D(I-1,J))*UA(I,J)
     .         +(D(I,J-1)+D(I-1,J-1))*UA(I,J-1))*
     .                        (VA(I-1,J)+VA(I,J))
          ENDDO
         ENDIF
         DO I=2,IMM1
      FLUXVA(I,J)=.125*((D(I,J+1)+D(I,J))
     .       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     .      *(VA(I,J+1)+VA(I,J))
         ENDDO
        ENDIF
        DO I=2,IMM1
      FLUXVA(I,J)=FLUXVA(I,J)
     .        -D(I,J)*2.*AAM2D(I,J)*(VAB(I,J+1)-VAB(I,J))/H2(I,J)
        ENDDO
        DO I=2,IM
      FLUXVA(I,J)=FLUXVA(I,J)*H1(I,J)
      FLUXUA(I,J)=(FLUXUA(I,J)-TPS(I,J)*DVM(I,J)*DVM(I-1,J)
     .            *DVMWAD(I,J)*DVMWAD(I-1,J)*DUMWAD(I,J)*DUMWAD(I,J-1) ! Not in original NYHOPS W&D
     .                                 *DUM(I,J)*DUM(I,J-1))! THIN DAMS NG08192008
     .             *.25*(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
        ENDDO
       ENDIF
      ENDDO
!$OMP FLUSH
      CALL LOCALSYNC43(T43,NG43,3*(E43-1)+3)
!$OMP FLUSH
      CALL RANGE43(MAX(3,LO43),MIN(JMM1,HI43),MOD(T43,NG43),NG43,L43,H43)
      DO J=L43,H43
       IF(J.EQ.LO43.AND.LO43.GT.1) THEN
        DO I=2,IMM1
      ADVVA(I,J)=FLUXUA(I+1,J)-FLUXUA(I,J)
     .          +FLUXVA(I,J)-YV43(I)
        ENDDO
       ELSE
        DO I=2,IMM1
      ADVVA(I,J)=FLUXUA(I+1,J)-FLUXUA(I,J)
     .          +FLUXVA(I,J)-FLUXVA(I,J-1)
        ENDDO
       ENDIF
      ENDDO
!$OMP BARRIER
      END
      RECURSIVE SUBROUTINE UFACE43(J,YY43)
      INCLUDE 'comdeck'
      REAL YY43(IM),TT43
      YY43=0.0
       IF(J.GE.2) THEN
        IF(ADVECT.EQ.'NON-LINEAR') THEN
         DO I=2,IMM1
      YY43(I)=.125*((D(I,J)+D(I,J-1))*VA(I,J)
     .                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     .                 *(UA(I,J)+UA(I,J-1))
         ENDDO
        ENDIF
        DO I=2,IMM1
      TT43=.25*(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     .            *(AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1))
     .                *((UAB(I,J)-UAB(I,J-1))
     .                /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .                 +(VAB(I,J)-VAB(I-1,J))
     .               /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)) )
!-------- NOTE TPS CARRIED OVER TO V EQUATION !!! ----------------------
      YY43(I)=(YY43(I)-TT43*DUM(I,J)*DUM(I,J-1)
     .            *DUMWAD(I,J)*DUMWAD(I,J-1)*DVMWAD(I,J)*DVMWAD(I-1,J) ! Not in original NYHOPS W&D
     .                                 *DVM(I,J)*DVM(I-1,J))! THIN DAMS NG08192008
     .            *.25*(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1))
        ENDDO
       ENDIF
      END
      RECURSIVE SUBROUTINE VFACE43(J,YY43)
      INCLUDE 'comdeck'
      REAL YY43(IM),TT43
      YY43=0.0
       IF(J.GE.2.AND.J.LE.JMM1) THEN
        IF(ADVECT.EQ.'NON-LINEAR') THEN
         DO I=2,IMM1
      YY43(I)=.125*((D(I,J+1)+D(I,J))
     .       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     .      *(VA(I,J+1)+VA(I,J))
         ENDDO
        ENDIF
        DO I=2,IMM1
      YY43(I)=YY43(I)
     .        -D(I,J)*2.*AAM2D(I,J)*(VAB(I,J+1)-VAB(I,J))/H2(I,J)
        ENDDO
        DO I=2,IM
      YY43(I)=YY43(I)*H1(I,J)
        ENDDO
       ENDIF
      END
      RECURSIVE SUBROUTINE EDGE43(ADVUA,ADVVA)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM)
      IF (OPTEBC(15:20).EQ.'LINEAR') THEN ! IF LINEAR OPTION ENABLED AT THE END OF OPTEBC
       DO 309 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        IC=ICON(N)
        JC=JCON(N)
        IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE,JE)=0.0
         ADVVA(IE,JE)=0.0   ! LINEARIZE THIS AND BCOND IF NEEDED
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE+1,JE)=0.0
         ADVVA(IE+1,JE)=0.0 ! LINEARIZE THIS AND BCOND IF NEEDED
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE)=0.0
         ADVUA(IE,JE)=0.0   ! LINEARIZE THIS AND BCOND IF NEEDED
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE+1)=0.0
         ADVUA(IE,JE+1)=0.0 ! LINEARIZE THIS AND BCOND IF NEEDED
        ENDIF
 309   CONTINUE
      ELSE
       DO 310 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        IC=ICON(N)
        JC=JCON(N)
        IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE,JE)=0.0
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE+1,JE)=0.0
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE)=0.0
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE+1)=0.0
        ENDIF
 310   CONTINUE
      ENDIF
      END

      RECURSIVE SUBROUTINE RANGE43(LO,HI,ME,N,L,H)
      INTEGER LO,HI,ME,N,L,H,LEN43,K,EXTRA
      LEN43=MAX(0,HI-LO+1)
      K=LEN43/N
      EXTRA=MOD(LEN43,N)
      L=LO+ME*K+MIN(ME,EXTRA)
      H=L+K-1
      IF(ME.LT.EXTRA) H=H+1
      END
