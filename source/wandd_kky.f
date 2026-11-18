      SUBROUTINE WANDD_KKY(ELEV,TD)
      INCLUDE 'comdeck'
      REAL ELEV(IM,JM),TD(IM,JM)
!$OMP PARALLEL DEFAULT(SHARED)
      CALL WANDD_TEAM(ELEV,TD)
!$OMP END PARALLEL
      END
      RECURSIVE SUBROUTINE WANDD_TEAM(ELEV,TD)
!
!     INTRODUCE WET/DRY SCHEME BASED ON KKY's DIISERTATION (hli, 08/03/05).
!NGMAR2009   SIMPLIFIED AND ADDED AN EXTRA CONDITION (NO FLOW OUT OF "DRY" CELL).
!
!
      INCLUDE 'comdeck'
!
      DIMENSION ELEV(IM,JM),TD(IM,JM)
!
!     Owner (i,j) handles both tests of its east and north faces. Its
!     ET/ETF scatter may also touch row j+1. Alternate row colors with
!     a barrier so there are no concurrent writes (even equal writes).
!     The sole caller supplies EL,D, neither aliases the written fields.
      DO ICOLOR=1,2
!$OMP DO SCHEDULE(STATIC)
      DO J = ICOLOR, JMM1, 2
         DO I = 1, IMM1
               IF(TD(I,J).LE.WETMIN.AND.FSM(I,J).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I+1,J)).LE.WETMIN).OR.
     &            ((ELEV(I+1,J)-ELEV(I,J)).LE.WETEPS).OR.
     *             (UA(I+1,J).GT.0).OR.DUMWAD(I+1,J).EQ.0.0) THEN
                     UA(I+1,J)=0.0
                     UAF(I+1,J)=0.0
                     DUMWAD(I+1,J)=0.0
                     ET(I,J)=ETB(I,J) ! NG04012011
                     ETF(I,J)=ETB(I,J) ! NG04012011
                  ENDIF
                  IF((0.5*(TD(I,J)+TD(I,J+1)).LE.WETMIN).OR.
     &            ((ELEV(I,J+1)-ELEV(I,J)).LE.WETEPS).OR.
     *             (VA(I,J+1).GT.0).OR.DVMWAD(I,J+1).EQ.0.0) THEN
                     VA(I,J+1)=0.0
                     VAF(I,J+1)=0.0
                     DVMWAD(I,J+1)=0.0
                     ET(I,J)=ETB(I,J) ! NG04012011
                     ETF(I,J)=ETB(I,J) ! NG04012011
                  ENDIF
!                  write(*,*)'W&D1 ',time,I,J,TD(I,J),UA(I+1,J),VA(I,J+1)
               ENDIF
!
               IF(TD(I+1,J).LE.WETMIN.AND.DUM(I+1,J).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I+1,J)).LE.WETMIN).OR.
     &            ((ELEV(I,J)-ELEV(I+1,J)).LE.WETEPS).OR.
     *             (UA(I+1,J).LT.0).OR.DUMWAD(I+1,J).EQ.0.0) THEN
                     UA(I+1,J)=0.0
                     UAF(I+1,J)=0.0
                     DUMWAD(I+1,J)=0.0
                     ET(I+1,J)=ETB(I+1,J) ! NG04012011
                     ETF(I+1,J)=ETB(I+1,J) ! NG04012011
!                  write(*,*)'W&D2 ',time,I+1,J,TD(I+1,J),UA(I+1,J)
                  ENDIF
               ENDIF
!
               IF(TD(I,J+1).LE.WETMIN.AND.DVM(I,J+1).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I,J+1)).LE.WETMIN).OR.
     &            ((ELEV(I,J)-ELEV(I,J+1)).LE.WETEPS).OR.
     *             (VA(I,J+1).LT.0).OR.DVMWAD(I,J+1).EQ.0.0) THEN
                     VA(I,J+1)=0.0
                     VAF(I,J+1)=0.0
                     DVMWAD(I,J+1)=0.0
                     ET(I,J+1)=ETB(I,J+1) ! NG04012011
                     ETF(I,J+1)=ETB(I,J+1) ! NG04012011
!                  write(*,*)'W&D3 ',time,I,J+1,TD(I,J+1),VA(I,J+1)
                  ENDIF
               ENDIF
         ENDDO 
      ENDDO 

!$OMP END DO
      ENDDO
      RETURN
      END
!
! Only the persistent, no-transport BARO_SUBCYCLE calls this entry.
! EL/D are not ET/ETF/ETB. The public WANDD_KKY keeps its old semantics.
      RECURSIVE SUBROUTINE WANDD_SEAM37(SEAM37)
      USE ISO_FORTRAN_ENV, ONLY: INT32
!$    USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
      INCLUDE 'comdeck'
      INTEGER(INT32) SEAM37((IMM1+31)/32+16,JM)
      INTEGER NT37,MY37,NR37,NX37,JLO37,JHI37,IW37,I37
      NT37=1
      MY37=0
!$    NT37=OMP_GET_NUM_THREADS()
!$    MY37=OMP_GET_THREAD_NUM()
      IF(NT37.LE.1.OR.NT37.GT.JMM1) THEN
! Serial and teams larger than the number of face-owner rows use
! the original two colors (including arbitrary empty workshares).
        CALL WANDD_TEAM(EL,D)
        RETURN
      ENDIF
      NR37=JMM1/NT37
      NX37=MOD(JMM1,NT37)
      JLO37=1+MY37*NR37+MIN(MY37,NX37)
      JHI37=JLO37+NR37-1
      IF(MY37.LT.NX37) JHI37=JHI37+1
      IF(MY37.EQ.NT37-1) THEN
! Last worker also owns ET/ETF(:,JM). There is no receiver north.
        CALL WANDD_ROWS37(EL,D,JLO37,JHI37,.FALSE.,SEAM37(1,MY37+1))
      ELSE
        DO IW37=1,(IMM1+31)/32
          SEAM37(IW37,MY37+1)=0_INT32
        ENDDO
        CALL WANDD_ROWS37(EL,D,JLO37,JHI37-1,.FALSE.,SEAM37(1,MY37+1))
        CALL WANDD_ROWS37(EL,D,JHI37,JHI37,.TRUE.,SEAM37(1,MY37+1))
      ENDIF
! Publish all speed/face-mask updates and all deferred North requests.
!$OMP BARRIER
      IF(MY37.GT.0) THEN
        DO IW37=1,(IMM1+31)/32
          IF(SEAM37(IW37,MY37).NE.0_INT32) THEN
            DO I37=32*(IW37-1)+1,MIN(32*IW37,IMM1)
              IF(BTEST(SEAM37(IW37,MY37),MOD(I37-1,32))) THEN
                ET(I37,JLO37)=ETB(I37,JLO37)
                ETF(I37,JLO37)=ETB(I37,JLO37)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
! No exit barrier: immediately following CHECK reads only face masks,
! UAF/VAF, FSM, ELF. Its reduction barrier precedes any buffer reuse.
      END
      RECURSIVE SUBROUTINE WANDD_ROWS37(ELEV,TD,JLO37,JHI37,DEFER37,MASK37)
      USE ISO_FORTRAN_ENV, ONLY: INT32
      INCLUDE 'comdeck'
      REAL ELEV(IM,JM),TD(IM,JM)
      LOGICAL DEFER37
      INTEGER(INT32) MASK37(*)
      INTEGER JLO37,JHI37,IW37
      DO J=JLO37,JHI37
         DO I = 1, IMM1
               IF(TD(I,J).LE.WETMIN.AND.FSM(I,J).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I+1,J)).LE.WETMIN).OR.
     &            ((ELEV(I+1,J)-ELEV(I,J)).LE.WETEPS).OR.
     *             (UA(I+1,J).GT.0).OR.DUMWAD(I+1,J).EQ.0.0) THEN
                     UA(I+1,J)=0.0
                     UAF(I+1,J)=0.0
                     DUMWAD(I+1,J)=0.0
                     ET(I,J)=ETB(I,J) ! NG04012011
                     ETF(I,J)=ETB(I,J) ! NG04012011
                  ENDIF
                  IF((0.5*(TD(I,J)+TD(I,J+1)).LE.WETMIN).OR.
     &            ((ELEV(I,J+1)-ELEV(I,J)).LE.WETEPS).OR.
     *             (VA(I,J+1).GT.0).OR.DVMWAD(I,J+1).EQ.0.0) THEN
                     VA(I,J+1)=0.0
                     VAF(I,J+1)=0.0
                     DVMWAD(I,J+1)=0.0
                     ET(I,J)=ETB(I,J) ! NG04012011
                     ETF(I,J)=ETB(I,J) ! NG04012011
                  ENDIF
!                  write(*,*)'W&D1 ',time,I,J,TD(I,J),UA(I+1,J),VA(I,J+1)
               ENDIF
!
               IF(TD(I+1,J).LE.WETMIN.AND.DUM(I+1,J).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I+1,J)).LE.WETMIN).OR.
     &            ((ELEV(I,J)-ELEV(I+1,J)).LE.WETEPS).OR.
     *             (UA(I+1,J).LT.0).OR.DUMWAD(I+1,J).EQ.0.0) THEN
                     UA(I+1,J)=0.0
                     UAF(I+1,J)=0.0
                     DUMWAD(I+1,J)=0.0
                     ET(I+1,J)=ETB(I+1,J) ! NG04012011
                     ETF(I+1,J)=ETB(I+1,J) ! NG04012011
!                  write(*,*)'W&D2 ',time,I+1,J,TD(I+1,J),UA(I+1,J)
                  ENDIF
               ENDIF
!
               IF(TD(I,J+1).LE.WETMIN.AND.DVM(I,J+1).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I,J+1)).LE.WETMIN).OR.
     &            ((ELEV(I,J)-ELEV(I,J+1)).LE.WETEPS).OR.
     *             (VA(I,J+1).LT.0).OR.DVMWAD(I,J+1).EQ.0.0) THEN
                     VA(I,J+1)=0.0
                     VAF(I,J+1)=0.0
                     DVMWAD(I,J+1)=0.0
                     IF(DEFER37) THEN
                        IW37=(I-1)/32+1
                        MASK37(IW37)=IBSET(MASK37(IW37),MOD(I-1,32))
                     ELSE
                     ET(I,J+1)=ETB(I,J+1) ! NG04012011
                     ETF(I,J+1)=ETB(I,J+1) ! NG04012011
                     ENDIF
!                  write(*,*)'W&D3 ',time,I,J+1,TD(I,J+1),VA(I,J+1)
                  ENDIF
               ENDIF
         ENDDO 
      ENDDO
      END
