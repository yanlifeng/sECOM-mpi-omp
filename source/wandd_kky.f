      SUBROUTINE WANDD_KKY(ELEV,TD)
      INCLUDE 'comdeck'
      REAL ELEV(IM,JM),TD(IM,JM)
      INTEGER(SELECTED_INT_KIND(2)) IWDFLAG(IM,JM)
      SAVE IWDFLAG
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J)
!$OMP DO SCHEDULE(STATIC)
      DO J=1,JM
       DO I=1,IM
        IWDFLAG(I,J)=0
        IF(TD(I,J).LE.WETMIN) IWDFLAG(I,J)=
     &     MERGE(1,0,FSM(I,J).NE.0.)+MERGE(2,0,DUM(I,J).NE.0.)+MERGE(4,0,DVM(I,J).NE.0.)
       ENDDO
      ENDDO
!$OMP END DO
      CALL WANDD_FLAGS_TEAM(ELEV,TD,IWDFLAG)
!$OMP END PARALLEL
      END
      RECURSIVE SUBROUTINE WANDD_FLAGS_TEAM(ELEV,TD,IWDFLAG)
!
!     INTRODUCE WET/DRY SCHEME BASED ON KKY's DIISERTATION (hli, 08/03/05).
!NGMAR2009   SIMPLIFIED AND ADDED AN EXTRA CONDITION (NO FLOW OUT OF "DRY" CELL).
!
!
      INCLUDE 'comdeck'
!
      DIMENSION ELEV(IM,JM),TD(IM,JM)
      INTEGER(SELECTED_INT_KIND(2)) IWDFLAG(IM,JM)
!
!     Owner (i,j) handles both tests of its east and north faces. Its
!     ET/ETF scatter may also touch row j+1. Alternate row colors with
!     a barrier so there are no concurrent writes (even equal writes).
!     The sole caller supplies EL,D, neither aliases the written fields.
      DO ICOLOR=1,2
!$OMP DO SCHEDULE(STATIC)
      DO J = ICOLOR, JMM1, 2
         DO IB0=1,IMM1,16
          ILAST=MIN(IB0+15,IMM1)
! A block may be skipped only if none of the three source predicates
! can hold (the east and north face read their neighbor's flag).
          IF(.NOT.ANY(IWDFLAG(IB0:ILAST+1,J).NE.0).AND.
     &       .NOT.ANY(IWDFLAG(IB0:ILAST,J+1).NE.0)) CYCLE
          DO I=IB0,ILAST
               IF(BTEST(IWDFLAG(I,J),0)) THEN
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
               IF(BTEST(IWDFLAG(I+1,J),1)) THEN
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
               IF(BTEST(IWDFLAG(I,J+1),2)) THEN
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
      ENDDO 

!$OMP END DO
      ENDDO
      RETURN
      END
