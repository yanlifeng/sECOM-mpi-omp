      SUBROUTINE WANDD_KKY(ELEV,TD)
      INCLUDE 'comdeck'
      REAL ELEV(IM,JM),TD(IM,JM)
!$OMP PARALLEL DEFAULT(SHARED)
      CALL WANDD_TEAM(ELEV,TD)
!$OMP END PARALLEL
      END
      RECURSIVE SUBROUTINE WANDD_TEAM(ELEV,TD)
!$    USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
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
!     ET/ETF scatter may also touch row j+1. Divide rows into contiguous
!     slabs; hold the final row of each slab back until the others finish.
!     The sole caller supplies EL,D, neither aliases the written fields.
      NTEAM=1
      ME=0
!$    NTEAM=OMP_GET_NUM_THREADS()
!$    ME=OMP_GET_THREAD_NUM()
      IF(NTEAM.EQ.1) THEN
        CALL WANDD_ROWS(ELEV,TD,1,JMM1)
      ELSEIF(JMM1/NTEAM.GE.2) THEN
        NROW=JMM1/NTEAM
        NREM=MOD(JMM1,NTEAM)
        JLO=1+ME*NROW+MIN(ME,NREM)
        JHI=JLO+NROW-1
        IF(ME.LT.NREM) JHI=JHI+1
        CALL WANDD_ROWS(ELEV,TD,JLO,JHI-1)
!$OMP BARRIER
!       Seam rows differ by at least two: their write regions are disjoint.
        CALL WANDD_ROWS(ELEV,TD,JHI,JHI)
!$OMP BARRIER
      ELSE
!     Too few rows for disjoint seams. The original two-color schedule
!     is race free for arbitrary team sizes.
      DO ICOLOR=1,2
!$OMP DO SCHEDULE(STATIC)
      DO J = ICOLOR, JMM1, 2
        CALL WANDD_ROWS(ELEV,TD,J,J)
      ENDDO
!$OMP END DO
      ENDDO
      ENDIF
      END
      RECURSIVE SUBROUTINE WANDD_ROWS(ELEV,TD,JFIRST,JLAST)
      INCLUDE 'comdeck'
      DIMENSION ELEV(IM,JM),TD(IM,JM)
      DO J=JFIRST,JLAST
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

      RETURN
      END
