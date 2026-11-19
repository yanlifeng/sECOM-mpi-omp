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
      INTEGER WAD_SEND(IM),JSTART,JSTOP
      LOGICAL MULTI_WAD,NORTH_GHOST
! One outgoing row per process. This SAVE is intentional: all OpenMP
! callers contribute to and publish the same buffer after the DO join.
      SAVE WAD_SEND
!
!     Owner (i,j) handles both tests of its east and north faces. Its
!     ET/ETF scatter may also touch row j+1. Alternate row colors with
!     a barrier so there are no concurrent writes (even equal writes).
!     The sole caller supplies EL,D, neither aliases the written fields.
      MULTI_WAD=SMPI_SIZE.GT.1
      IF(MULTI_WAD) THEN
        CALL SMPI_HALO(ELEV,IM,JM)
        CALL SMPI_HALO(TD,IM,JM)
        CALL SMPI_HALO(VA,IM,JM)
        CALL SMPI_HALO(DVMWAD,IM,JM)
      ENDIF
      DO ICOLOR=1,2
      JSTART=ICOLOR
      JSTOP=JMM1
      IF(MULTI_WAD) THEN
! Use global parity, including when the local strip has just one row.
        JSTART=SMPI_COLOR(ICOLOR,JMM1,2)
        JSTOP=SMPI_LAST(JMM1)
!$OMP SINGLE
        WAD_SEND=0
!$OMP END SINGLE
      ENDIF
!$OMP DO SCHEDULE(STATIC) PRIVATE(NORTH_GHOST)
      DO J = JSTART, JSTOP, 2
         NORTH_GHOST=MULTI_WAD.AND.J.EQ.SMPI_JHI.AND.J.LT.JM
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
! Bit 0: deliver this exact face closure, even if its old bits were +0.
                     IF(NORTH_GHOST) WAD_SEND(I)=IBSET(WAD_SEND(I),0)
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
                     IF(NORTH_GHOST) THEN
! Bits 0 and 1: zero the face and copy the target OWNER'S ETB into
! ET and ETF. Do not read or send a possibly stale ghost ETB.
                       WAD_SEND(I)=3
                     ELSE
                       ET(I,J+1)=ETB(I,J+1) ! NG04012011
                       ETF(I,J+1)=ETB(I,J+1) ! NG04012011
                     ENDIF
!                  write(*,*)'W&D3 ',time,I,J+1,TD(I,J+1),VA(I,J+1)
                  ENDIF
               ENDIF
         ENDDO 
      ENDDO 

!$OMP END DO
      IF(MULTI_WAD) THEN
! Only marked columns change on the remote owner, before the next
! color. Stale unmarked ghost values must never overwrite owner data.
        CALL SMPI_WAD_PUBLISH(WAD_SEND,ET,ETF,VA,VAF,DVMWAD,ETB,IM,JM)
        CALL SMPI_HALO(VA,IM,JM)
        CALL SMPI_HALO(DVMWAD,IM,JM)
        IF(ICOLOR.EQ.2) CALL SMPI_HALO(VAF,IM,JM)
      ENDIF
      ENDDO
      RETURN
      END
