! Differential MPI kernel fixture. No solver, input files or NetCDF calls.
      PROGRAM TEST_WAD_OWNED
      USE MPI, ONLY: MPI_ALLREDUCE,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD
      USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS
      INCLUDE 'comdeck'
      INTEGER EXPECT(IM,JM,8),NFAIL,NALL,NTEST,NTEXP,IERR,M,IGLOBAL
      CHARACTER*16 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) NTEST
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) NTEXP
      CALL SMPI_INIT(IM,JM)
      NFAIL=0
      IF(SMPI_SIZE.NE.NTEST) NFAIL=NFAIL+1
!$OMP PARALLEL
!$OMP SINGLE
      IF(OMP_GET_NUM_THREADS().NE.NTEXP) NFAIL=NFAIL+1
!$OMP END SINGLE
!$OMP END PARALLEL
      DO M=1,8
! Each rank independently computes the frozen v46 oracle **only in QA**.
        CALL WAD_FIXTURE(M)
        CALL V46_WANDD_KKY(EL,D)
        DO J=1,JM
          DO I=1,IM
            EXPECT(I,J,:)=(/TRANSFER(UA(I,J),0),TRANSFER(UAF(I,J),0),
     &          TRANSFER(DUMWAD(I,J),0),TRANSFER(VA(I,J),0),
     &          TRANSFER(VAF(I,J),0),TRANSFER(DVMWAD(I,J),0),
     &          TRANSFER(ET(I,J),0),TRANSFER(ETF(I,J),0)/)
          ENDDO
        ENDDO
        CALL WAD_FIXTURE(M)
! All non-owned mutable rows are explicitly stale. Static masks are
! replicated; ETB must be read on the receiving owner, not its neighbor.
        DO J=1,JM
          IF(J.GE.SMPI_JLO.AND.J.LE.SMPI_JHI) CYCLE
          D(:,J)=-271.
          EL(:,J)=682.
          UA(:,J)=487.
          UAF(:,J)=-425.
          DUMWAD(:,J)=650.
          VA(:,J)=-517.
          VAF(:,J)=-719.
          DVMWAD(:,J)=385.
          ET(:,J)=-233.
          ETF(:,J)=128.
          ETB(:,J)=-1024.
        ENDDO
        CALL WANDD_KKY(EL,D)
        DO J=SMPI_JLO,SMPI_JHI
          DO I=1,IM
            IF(ANY(EXPECT(I,J,:).NE.(/TRANSFER(UA(I,J),0),TRANSFER(UAF(I,J),0),
     &          TRANSFER(DUMWAD(I,J),0),TRANSFER(VA(I,J),0),
     &          TRANSFER(VAF(I,J),0),TRANSFER(DVMWAD(I,J),0),
     &          TRANSFER(ET(I,J),0),TRANSFER(ETF(I,J),0)/))) THEN
              IF(NFAIL.LT.3) PRINT *, 'BAD OWNER',SMPI_RANK,M,I,J
              NFAIL=NFAIL+1
            ENDIF
          ENDDO
        ENDDO
! EXTRNL/CHECK callers are promised current immediate neighbors for
! all three north-face data; both directions matter after the last color.
        DO J=MAX(1,SMPI_JLO-1),MIN(JM,SMPI_JHI+1)
          IF(J.GE.SMPI_JLO.AND.J.LE.SMPI_JHI) CYCLE
          DO I=1,IM
            IF(ANY(EXPECT(I,J,4:6).NE.(/TRANSFER(VA(I,J),0),
     &          TRANSFER(VAF(I,J),0),TRANSFER(DVMWAD(I,J),0)/))) THEN
              IF(NFAIL.LT.3) PRINT *, 'BAD HALO',SMPI_RANK,M,I,J
              NFAIL=NFAIL+1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CALL MPI_ALLREDUCE(NFAIL,NALL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERR)
      IF(SMPI_RANK.EQ.0) PRINT *, 'WAD_OWNED',SMPI_SIZE,NTEXP,NALL
      CALL SMPI_FINALIZE()
      IF(NALL.NE.0) STOP 91
      END

      SUBROUTINE WAD_FIXTURE(M)
      INCLUDE 'comdeck'
      INTEGER M
      WETMIN=.30
      WETEPS=.03
      DO J=1,JM
        DO I=1,IM
          FSM(I,J)=1.
          DUM(I,J)=1.
          DVM(I,J)=1.
          DUMWAD(I,J)=1.
          DVMWAD(I,J)=1.
          D(I,J)=1.
! Include wet/wet (no flags), dry/wet (face only) and wet/dry
! (face plus ET/ETF), in both parities and at every rank seam.
          IF(MOD(I,4).EQ.1.AND.MOD(J+M,2).EQ.0) D(I,J)=WETMIN
          IF(MOD(I,4).EQ.2.AND.MOD(J+M,2).EQ.1) D(I,J)=WETMIN-.03
          IF(MOD(I,4).EQ.3) D(I,J)=WETMIN-.05
          EL(I,J)=0.
          UA(I,J)=.031*I-.047*J
          UAF(I,J)=UA(I,J)*.71
          VA(I,J)=.059*I-.071*J
          VAF(I,J)=VA(I,J)*.83
          ET(I,J)=-.89*I-J
          ETF(I,J)=1.03*I+J
          ETB(I,J)=.125*I+J
          IF(MOD(I+J,4).EQ.0) ETB(I,J)=-0.0
          IF(M.EQ.3) THEN
            IF(MOD(I+2*J,3).EQ.0) FSM(I,J)=0.
            IF(MOD(I+J,4).EQ.0) DVM(I,J)=0.
            IF(MOD(I+4*J,5).EQ.0) DUM(I,J)=0.
          ELSEIF(M.EQ.4) THEN
            IF(MOD(I+J,3).EQ.0) DVMWAD(I,J)=0.
            EL(I,J)=.13*J
            IF(MOD(I,2).EQ.0) VAF(I,J)=-0.0
          ELSEIF(M.EQ.5) THEN
            D(I,J)=2.
          ELSEIF(M.EQ.6) THEN
            D(I,J)=WETMIN
          ELSEIF(M.EQ.7) THEN
            D(I,J)=.31
            IF(MOD(I*7+J*13,5).EQ.0) D(I,J)=.11
            IF(MOD(I+J*11,5).EQ.0) DVMWAD(I,J)=0.
            EL(I,J)=.08*MOD(I-J,3)
          ELSEIF(M.EQ.8) THEN
! 3-way split: rank 0 ends at J=3. Only its north-facing closure
! resets the dry cell (5,4): ensure that the other three faces stay
! open and cannot accidentally conceal a missing ET/ETF publication.
            D(I,J)=1.
            EL(I,J)=1.
            IF(I.EQ.5.AND.J.EQ.4) THEN
              D(I,J)=.1
              EL(I,J)=0.
              VA(I,J)=-.25
              UA(I,J)=.25
            ENDIF
            IF(I.EQ.6.AND.J.EQ.4) UA(I,J)=-.25
            IF(I.EQ.5.AND.J.EQ.5) VA(I,J)=-.25
          ENDIF
        ENDDO
      ENDDO
      END
