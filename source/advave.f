      SUBROUTINE ADVAVE(ADVUA,ADVVA)
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
!
      CALL ADV_SLAB(ADVUA,ADVVA)
!$OMP BARRIER
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
!
! Each thread owns a contiguous set of divergence rows. Process small
! tiles with one privately recomputed halo row on either side. The
! input state is read-only here, and only the owning thread writes the
! shared diagnostic/work arrays. No neighbor barriers inside the tile.
      RECURSIVE SUBROUTINE ADV_SLAB(ADVUA,ADVVA)
!$    USE OMP_LIB, ONLY: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM)
      INTEGER, PARAMETER :: NROW=8
      REAL XF(IM,0:NROW+1),YF(IM,0:NROW+1)
      REAL XG(IM,0:NROW+1),YG(IM,0:NROW+1),ST(IM)
      INTEGER NT,ME,ILO,IHI,JB,JE,J,I,JLOC
      NT=1
      ME=0
!$    NT=OMP_GET_NUM_THREADS()
!$    ME=OMP_GET_THREAD_NUM()
      ILO=1+(JM*ME)/NT
      IHI=(JM*(ME+1))/NT
      DO JB=ILO,IHI,NROW
       JE=MIN(IHI,JB+NROW-1)
       DO J=MAX(1,JB-1),MIN(JM,JE+1)
        JLOC=J-JB+1
        DO I=1,IM
         XF(I,JLOC)=0.
         YF(I,JLOC)=0.
         XG(I,JLOC)=0.
         YG(I,JLOC)=0.
         ST(I)=0.
        ENDDO
        IF(J.GE.2) THEN
         IF(ADVECT.EQ.'NON-LINEAR') THEN
          IF(J.LE.JMM1) THEN
           DO I=2,IMM1
      XF(I,JLOC)=.125*((D(I+1,J)+D(I,J))*UA(I+1,J)
     .                 +(D(I,J)+D(I-1,J))*UA(I,J))
     .                  *(UA(I+1,J)+UA(I,J))
      YG(I,JLOC)=.125*((D(I,J+1)+D(I,J))
     .       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     .      *(VA(I,J+1)+VA(I,J))
           ENDDO
          ENDIF
          DO I=2,IMM1
      YF(I,JLOC)=.125*((D(I,J)+D(I,J-1))*VA(I,J)
     .                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     .                 *(UA(I,J)+UA(I,J-1))
          ENDDO
          IF(J.GE.3.AND.J.LE.JMM1) THEN
           DO I=2,IM
      XG(I,JLOC)=.125*((D(I,J)+D(I-1,J))*UA(I,J)
     .         +(D(I,J-1)+D(I-1,J-1))*UA(I,J-1))*
     .                        (VA(I-1,J)+VA(I,J))
           ENDDO
          ENDIF
         ENDIF
         IF(J.LE.JMM1) THEN
          DO I=2,IMM1
      XF(I,JLOC)=XF(I,JLOC)
     .         -D(I,J)*2.*AAM2D(I,J)*(UAB(I+1,J)-UAB(I,J))/H1(I,J)
      YG(I,JLOC)=YG(I,JLOC)
     .        -D(I,J)*2.*AAM2D(I,J)*(VAB(I,J+1)-VAB(I,J))/H2(I,J)
          ENDDO
         ENDIF
         DO I=2,IMM1
      ST(I)=.25*(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     .            *(AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1))
     .                *((UAB(I,J)-UAB(I,J-1))
     .                /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .                 +(VAB(I,J)-VAB(I-1,J))
     .               /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)) )
      XF(I,JLOC)=XF(I,JLOC)*H2(I,J)
      YF(I,JLOC)=(YF(I,JLOC)-ST(I)*DUM(I,J)*DUM(I,J-1)
     .            *DUMWAD(I,J)*DUMWAD(I,J-1)*DVMWAD(I,J)*DVMWAD(I-1,J) ! Not in original NYHOPS W&D
     .                                 *DVM(I,J)*DVM(I-1,J))! THIN DAMS NG08192008
     .            *.25*(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1))
         ENDDO
         IF(J.LE.JMM1) THEN
          DO I=2,IM
      YG(I,JLOC)=YG(I,JLOC)*H1(I,J)
      XG(I,JLOC)=(XG(I,JLOC)-ST(I)*DVM(I,J)*DVM(I-1,J)
     .            *DVMWAD(I,J)*DVMWAD(I-1,J)*DUMWAD(I,J)*DUMWAD(I,J-1) ! Not in original NYHOPS W&D
     .                                 *DUM(I,J)*DUM(I,J-1))! THIN DAMS NG08192008
     .             *.25*(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
          ENDDO
         ENDIF
        ENDIF
        IF(J.GE.JB.AND.J.LE.JE) THEN
         DO I=1,IM
          TPS(I,J)=ST(I)
          FLUXUA(I,J)=XG(I,JLOC)
          FLUXVA(I,J)=YG(I,JLOC)
         ENDDO
        ENDIF
       ENDDO
       DO J=JB,JE
        JLOC=J-JB+1
        DO I=1,IM
         ADVUA(I,J)=0.
         ADVVA(I,J)=0.
        ENDDO
        IF(J.GE.2.AND.J.LE.JMM1) THEN
         DO I=3,IMM1
          ADVUA(I,J)=XF(I,JLOC)-XF(I-1,JLOC)+YF(I,JLOC+1)-YF(I,JLOC)
         ENDDO
        ENDIF
        IF(J.GE.3.AND.J.LE.JMM1) THEN
         DO I=2,IMM1
          ADVVA(I,J)=XG(I+1,JLOC)-XG(I,JLOC)+YG(I,JLOC)-YG(I,JLOC-1)
         ENDDO
        ENDIF
       ENDDO
      ENDDO
      END
