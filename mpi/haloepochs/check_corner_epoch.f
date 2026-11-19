      PROGRAM CHECK_CORNER_EPOCH
      INCLUDE 'comdeck'
      INTEGER NSLAB,CASECE,ERRCE,XRANK,JCUT
      LOGICAL WBC,WU,WV
      CHARACTER*16 TARG
      CALL GET_COMMAND_ARGUMENT(1,TARG)
      READ(TARG,*) NSLAB
      CALL SMPI_INIT(IM,JM)
      IF(SMPI_SIZE.NE.NSLAB) STOP 81
      TOR='BAROTROPIC'
      BCTYPE='CLAMPED'
      OPTEBC=' '
      OPTEBC(15:20)='LINEAR'
      NUMQBC=0
      NUMDBC=0
      NUMEPC=0
      NUMEBC=0
      FSM=1.
      ERRCE=0
      JCUT=JM/NSLAB+MIN(1,MOD(JM,NSLAB))
      IF(NSLAB.EQ.1) JCUT=JM/2
      DO CASECE=1,6
        NUMECR=0
        WBC=.TRUE.
        WU=.FALSE.
        WV=.FALSE.
        IF(CASECE.NE.1) THEN
          NUMECR=1
          ICHARD(1)=5
          JCHARD(1)=JCUT
          NCHARD(1)=3
          SELECT CASE(CASECE)
          CASE(2)
! Top-at-last-owned-J writes the UA needed at next rank's J-1.
            WU=.TRUE.
          CASE(3)
! First row on second rank has a remote top-corner source at J-1.
            JCHARD(1)=JCUT+1
            WBC=.FALSE.
            WU=.TRUE.
            WV=.TRUE.
          CASE(4)
! Bottom at first row of second rank writes VA on J+1 and UA on J;
! neither is an exported VA-low/UA-high edge (slabs have >=3 rows).
            NCHARD(1)=1
            JCHARD(1)=JCUT+1
          CASE(5)
! Only rank zero reports an unsafe UA publication. All ranks agree.
            IF(SMPI_RANK.NE.0) NUMECR=0
            WU=.TRUE.
          CASE(6)
! Same-owner top in the middle, no touched halo edge.
            JCHARD(1)=JCUT-1
          END SELECT
        END IF
        IF(NSLAB.EQ.1) THEN
          WBC=.FALSE.
          WU=.TRUE.
          WV=.TRUE.
        END IF
        CALL SECOM_MPI_BC_PLAN()
        IF((SMPI_LOCALBC.NEQV.WBC).OR.(SMPI_KEEP_UA.NEQV.WU).OR.
     &     (SMPI_KEEP_VA.NEQV.WV)) THEN
          WRITE(*,*) 'FAIL_CORNER_EPOCH rank/case/got/want=',SMPI_RANK,
     &      CASECE,SMPI_LOCALBC,SMPI_KEEP_UA,SMPI_KEEP_VA,WBC,WU,WV
          ERRCE=ERRCE+1
        END IF
      END DO
      CALL SMPI_INDEX_RAW(ERRCE)
      XRANK=SMPI_RANK
      CALL SMPI_FINALIZE()
      IF(ERRCE.NE.0) STOP 82
      IF(XRANK.EQ.0) WRITE(*,*) 'PASS_CORNER_EPOCH cases/ranks=',6,NSLAB
      END
