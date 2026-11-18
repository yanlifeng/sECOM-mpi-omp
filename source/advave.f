      SUBROUTINE ADVAVE(ADVUA,ADVVA)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM)
      CALL ADVAVE_METRICS(ADVUA,ADVVA,.TRUE.)
      END
      SUBROUTINE ADVAVE_METRICS(ADVUA,ADVVA,REFRESH)
!     VERSION(09/27/90)
!
      INCLUDE 'comdeck'
      SAVE
!
!-----------------------------------------------------------------------
!      THIS SUBROUTINE CALCULATES ADVECTION & DIFFUSION
!-----------------------------------------------------------------------
!
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM)
      LOGICAL REFRESH
      REAL H1SUM(IM,JM),H2SUM(IM,JM),AAMSUM(IM,JM)
!
! Keep the existing U/V scratch arrays, but finish each row's fluxes
! while it is hot. Only divergences read other rows, and U divergence
! must finish before the V stage reuses the shared FLUX arrays.
!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J)
      DO J=1,JM
       DO I=1,IM
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
      ADVUA(I,J)=0.0
      TPS(I,J)=0.0
       ENDDO
       IF(J.GE.2) THEN
! Cache only the input geometry and viscosity sums. The raw operands
! are constant within the caller's ten (in general ISPLIT) substeps;
! SMAG can change AAM2D on the next outer step. Never cache D or masks.
        IF(REFRESH) THEN
         DO I=2,IMM1
          H1SUM(I,J)=H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)
          H2SUM(I,J)=H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1)
          AAMSUM(I,J)=AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1)
         ENDDO
         IF(J.LE.JMM1) H2SUM(IM,J)=H2(IM,J)+H2(IM-1,J)+H2(IM,J-1)+H2(IM-1,J-1)
        ENDIF
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
     .            *AAMSUM(I,J)
     .                *((UAB(I,J)-UAB(I,J-1))
     .                /H2SUM(I,J)
     .                 +(VAB(I,J)-VAB(I-1,J))
     .               /H1SUM(I,J) )
!-------- NOTE TPS CARRIED OVER TO V EQUATION !!! ----------------------
      FLUXUA(I,J)=FLUXUA(I,J)*H2(I,J)
      FLUXVA(I,J)=(FLUXVA(I,J)-TPS(I,J)*DUM(I,J)*DUM(I,J-1)
     .            *DUMWAD(I,J)*DUMWAD(I,J-1)*DVMWAD(I,J)*DVMWAD(I-1,J) ! Not in original NYHOPS W&D
     .                                 *DVM(I,J)*DVM(I-1,J))! THIN DAMS NG08192008
     .            *.25*H1SUM(I,J)
        ENDDO
       ENDIF
      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J)
      DO 480 J=2,JMM1
      DO 480 I=3,IMM1
 480  ADVUA(I,J)=FLUXUA(I,J)-FLUXUA(I-1,J)
     .           +FLUXVA(I,J+1)-FLUXVA(I,J)
!$OMP END DO
!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J)
      DO J=1,JM
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
     .             *.25*H2SUM(I,J)
        ENDDO
       ENDIF
      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J)
      DO 880 J=3,JMM1
      DO 880 I=2,IMM1
 880  ADVVA(I,J)=FLUXUA(I+1,J)-FLUXUA(I,J)
     .          +FLUXVA(I,J)-FLUXVA(I,J-1)
!$OMP END DO
!
! Keep potentially overlapping boundary entries in original order.
!$OMP SINGLE
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
!$OMP END SINGLE
!
      RETURN
      END
