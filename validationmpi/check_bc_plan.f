      PROGRAM TEST_MPI_BOUNDARY_PLAN
! No solver, files, or model steps. argv: exact expected rank count.
      INCLUDE 'comdeck'
      INTEGER NRCBC,CASEBC,JSPLBC,NERRBC,ITBC
      LOGICAL WANTBC
      CHARACTER*16 ARGBC
      CALL GET_COMMAND_ARGUMENT(1,ARGBC)
      READ(ARGBC,*) NRCBC
      CALL SMPI_INIT(IM,JM)
      NERRBC=0
      IF(SMPI_SIZE.NE.NRCBC.OR.SMPI_LOCALBC) NERRBC=1
      JSPLBC=JM/SMPI_SIZE+MIN(1,MOD(JM,SMPI_SIZE))
      IF(SMPI_SIZE.EQ.1) JSPLBC=JM/2
      DO CASEBC=1,34
        FSM=1.
        TOR='BAROTROPIC'
        BCTYPE='CLAMPED'
        OPTEBC=' '
        OPTEBC(15:20)='LINEAR'
        NUMQBC=0
        NUMDBC=0
        NUMEPC=0
        NUMEBC=1
        NUMECR=1
        IETA(1)=10
        JETA(1)=JSPLBC
        ICON(1)=10
        JCON(1)=JSPLBC
        ICHARD(1)=10
        JCHARD(1)=5
        NCHARD(1)=1
        WANTBC=.FALSE.
        SELECT CASE(CASEBC)
        CASE(1:4)
          NCHARD(1)=CASEBC
          WANTBC=.TRUE.
        CASE(5)
          TOR='PROGNOSTIC'
        CASE(6)
          BCTYPE='MIXED'
        CASE(7)
          OPTEBC(15:20)=' '
        CASE(8)
          NUMQBC=1
        CASE(9)
          NUMDBC=1
        CASE(10)
          NUMEPC=1
        CASE(11)
! N-S 310 copy from the next rank.
          FSM(10,JSPLBC-1)=-0.
        CASE(12)
! N-S 310 copy from the preceding rank.
          JETA(1)=JSPLBC+1
          JCON(1)=JSPLBC+1
          FSM(10,JSPLBC+2)=0.
        CASE(13:14)
! Earlier E/W case wins; the N-S copy would cross the seam if executed.
          FSM(10,JSPLBC-1)=0.
          IF(CASEBC.EQ.13) FSM(11,JSPLBC)=-0.
          IF(CASEBC.EQ.14) FSM(9,JSPLBC)=0.
          WANTBC=.TRUE.
        CASE(15:16)
! The bottom corner's j and j+1 share an owner, but j+2 does not.
          NCHARD(1)=CASEBC-14
          JCHARD(1)=JSPLBC-1
        CASE(17:18)
          NCHARD(1)=CASEBC-14
          JCHARD(1)=JSPLBC+1
        CASE(19)
          IF(SMPI_RANK.EQ.SMPI_SIZE-1) NUMQBC=1
        CASE(20)
          NUMEBC=EBCM+1
        CASE(21)
          NUMECR=-1
        CASE(22)
          IETA(1)=IM
        CASE(23)
          JETA(1)=1
        CASE(24)
          JCON(1)=JM+1
        CASE(25)
          NCHARD(1)=5
        CASE(26)
          ICHARD(1)=IM-1
        CASE(27)
          NCHARD(1)=2
          ICHARD(1)=1
        CASE(28)
          JCHARD(1)=JM-1
        CASE(29)
          NCHARD(1)=4
          JCHARD(1)=1
        CASE(30)
! sNaN integer bit injection must not trap, even under -ffpe-trap.
          FSM(11,JSPLBC)=TRANSFER(INT(Z'7F800001'),0.)
        CASE(31)
          FSM(11,JSPLBC)=TRANSFER(1,0.)
        CASE(32)
! No possible copies with both indices different.
          ICON(1)=12
          JCON(1)=JSPLBC+1
          FSM(10,JSPLBC-1)=0.
          WANTBC=.TRUE.
        CASE(33)
          NUMEBC=0
          NUMECR=0
          WANTBC=.TRUE.
        CASE(34)
! Ordered duplicates are still local, not a reason to reject.
          NUMEBC=2
          IETA(2)=IETA(1)
          JETA(2)=JETA(1)
          ICON(2)=ICON(1)
          JCON(2)=JCON(1)
          FSM(11,JSPLBC)=0.
          WANTBC=.TRUE.
        END SELECT
        IF(SMPI_SIZE.EQ.1) WANTBC=.FALSE.
        CALL SECOM_MPI_BC_PLAN()
        IF(SMPI_LOCALBC.NEQV.WANTBC) THEN
          WRITE(*,*) 'FAIL_MPI_BC_PLAN rank/case/got/want=',
     &               SMPI_RANK,CASEBC,SMPI_LOCALBC,WANTBC
          NERRBC=NERRBC+1
        END IF
      END DO
      CALL SMPI_INDEX_RAW(NERRBC)
      ITBC=SMPI_RANK
      CALL SMPI_FINALIZE()
      IF(NERRBC.NE.0.OR.SMPI_LOCALBC) STOP 21
      IF(ITBC.EQ.0) WRITE(*,*) 'PASS_MPI_BC_PLAN cases/ranks=',34,NRCBC
      END
