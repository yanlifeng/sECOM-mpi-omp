      SUBROUTINE SECOM_MPI_BC_PLAN()
! Serial startup, called by every rank after FSM and the boundary lists
! have their final values. They are geometry, not WAD's evolving masks.
! Re-plan if that geometry ever changes. Any rank can veto the shortcut.
      INCLUDE 'comdeck'
      INTEGER IBADBC,INEBC,IEBC,JEBC,ICBC,JCBC,JSBC,ITSBC
      INTEGER IMASKBC(4),IMAGBC(4),IKBC,ICORB,JCORB,ITYPBC
      SMPI_LOCALBC=.FALSE.
      IF(SMPI_SIZE.LE.1) RETURN
      IBADBC=1
      IF(TOR.NE.'BAROTROPIC'.OR.BCTYPE.NE.'CLAMPED') GOTO 90
      IF(OPTEBC(15:20).NE.'LINEAR') GOTO 90
      IF(NUMQBC.NE.0.OR.NUMDBC.NE.0.OR.NUMEPC.NE.0) GOTO 90
      IF(NUMEBC.LT.0.OR.NUMEBC.GT.EBCM) GOTO 90
      IF(NUMECR.LT.0.OR.NUMECR.GT.EBCM) GOTO 90
      IF(STORAGE_SIZE(0.0).NE.32.OR.BIT_SIZE(0).NE.32) GOTO 90
      IF(TRANSFER(1.0,0).NE.INT(Z'3F800000')) GOTO 90
      DO INEBC=1,NUMEBC
        IEBC=IETA(INEBC)
        JEBC=JETA(INEBC)
        ICBC=ICON(INEBC)
        JCBC=JCON(INEBC)
        IF(IEBC.LT.2.OR.IEBC.GT.IM-1) GOTO 90
        IF(JEBC.LT.2.OR.JEBC.GT.JM-1) GOTO 90
        IF(ICBC.LT.1.OR.ICBC.GT.IM) GOTO 90
        IF(JCBC.LT.1.OR.JCBC.GT.JM) GOTO 90
! Use integer loads only: do not raise an early FP exception on an sNaN
! and do not assume how a compiler/runtime handles subnormal comparisons.
        IMASKBC(1)=TRANSFER(FSM(IEBC+1,JEBC),0)
        IMASKBC(2)=TRANSFER(FSM(IEBC-1,JEBC),0)
        IMASKBC(3)=TRANSFER(FSM(IEBC,JEBC+1),0)
        IMASKBC(4)=TRANSFER(FSM(IEBC,JEBC-1),0)
        DO IKBC=1,4
          IMAGBC(IKBC)=IBCLR(IMASKBC(IKBC),31)
          IF(IMAGBC(IKBC).GE.INT(Z'7F800000')) GOTO 90
          IF(IMAGBC(IKBC).GT.0.AND.
     &       IMAGBC(IKBC).LT.INT(Z'00800000')) GOTO 90
        END DO
! Mirror the ELSE IF priority of EXTRNL_TEAM 310, including corners
! where more than one mask is zero. Preserve original list order at use.
        JSBC=JEBC
        IF(IMAGBC(1).EQ.0.AND.JEBC.EQ.JCBC) THEN
          JSBC=JEBC
        ELSE IF(IMAGBC(2).EQ.0.AND.JEBC.EQ.JCBC) THEN
          JSBC=JEBC
        ELSE IF(IMAGBC(3).EQ.0.AND.IEBC.EQ.ICBC) THEN
          JSBC=JEBC-1
        ELSE IF(IMAGBC(4).EQ.0.AND.IEBC.EQ.ICBC) THEN
          JSBC=JEBC+1
        END IF
        IF(SMPI_OWNER_RAW(JEBC).NE.SMPI_OWNER_RAW(JSBC)) GOTO 90
      END DO
      DO INEBC=1,NUMECR
        ICORB=ICHARD(INEBC)
        JCORB=JCHARD(INEBC)
        ITYPBC=NCHARD(INEBC)
        IF(ITYPBC.LT.1.OR.ITYPBC.GT.4) GOTO 90
        IF(ITYPBC.EQ.1.OR.ITYPBC.EQ.4) THEN
          IF(ICORB.LT.1.OR.ICORB.GT.IM-2) GOTO 90
        ELSE
          IF(ICORB.LT.2.OR.ICORB.GT.IM) GOTO 90
        END IF
! Bottom entries touch j,j+1,j+2 (two target rows). Top entries use
! j-1,j (one target row). Both statements must execute on the same rank.
        IF(ITYPBC.LE.2) THEN
          IF(JCORB.LT.1.OR.JCORB.GT.JM-2) GOTO 90
          ITSBC=SMPI_OWNER_RAW(JCORB)
          IF(ITSBC.NE.SMPI_OWNER_RAW(JCORB+1)) GOTO 90
          IF(ITSBC.NE.SMPI_OWNER_RAW(JCORB+2)) GOTO 90
        ELSE
          IF(JCORB.LT.2.OR.JCORB.GT.JM) GOTO 90
          IF(SMPI_OWNER_RAW(JCORB).NE.SMPI_OWNER_RAW(JCORB-1)) GOTO 90
        END IF
      END DO
      IBADBC=0
   90 CALL SMPI_INDEX_RAW(IBADBC)
      SMPI_LOCALBC=IBADBC.EQ.0
      END
