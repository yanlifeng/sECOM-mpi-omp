      PROGRAM HEALTH_CONTRACT
! Direct contract: fabricated summaries isolate integer preflight and
! actual winner-row rescan from the pre-existing dense SIMD stage.
      INCLUDE 'comdeck'
      REAL RVH(JM),REH(JM),WHV,WHE
      INTEGER IHCASE,IFH,IHR,JHR,NRH,JSMP
      LOGICAL OLDH,WHOLD
      CHARACTER*16 AH,MODEH
      CALL GET_COMMAND_ARGUMENT(1,AH)
      READ(AH,*) NRH
      CALL GET_COMMAND_ARGUMENT(2,MODEH)
      CALL SMPI_INIT(IM,JM)
      IFH=0
      IF(NRH.NE.SMPI_SIZE) STOP 88
      DO IHCASE=1,18
        RVH=0.
        REH=1.
        ELF=1.
        FSM=1.
        UAF=0.
        VAF=0.
        DUMWAD=1.
        DVMWAD=1.
        VAMAX=17.
        ELMINX=TRANSFER(IBSET(0,31),0.)
        IAMAX=-97
        JAMAX=-98
        WHV=0.
        WHE=1.
        IHR=1
        JHR=1
        WHOLD=.FALSE.
        SELECT CASE(IHCASE)
        CASE(1:2)
! Opposite zero signs on different ranks; spatial priority precedes sign.
          REH=2.
          ELF=2.
          REH(3)=0.
          REH(JM-1)=0.
          ELF(45,3)=MERGE(-0.,+0.,IHCASE.EQ.1)
          ELF(2,JM-1)=MERGE(+0.,-0.,IHCASE.EQ.1)
          WHE=ELF(45,3)
          IHR=45
          JHR=3
        CASE(3:6)
! Exceptional MASK or ELF at prefix of a losing rank's candidate row:
! integer veto, not a speculative floating-point exception.
          IF(SMPI_RANK.EQ.SMPI_SIZE-1) THEN
            IF(IHCASE.EQ.3) FSM(1,SMPI_JLO)=TRANSFER(INT(Z'7F800001'),0.)
            IF(IHCASE.EQ.4) FSM(1,SMPI_JLO)=TRANSFER(1,0.)
            IF(IHCASE.EQ.5) ELF(1,SMPI_JLO)=TRANSFER(INT(Z'7F800001'),0.)
            IF(IHCASE.EQ.6) ELF(1,SMPI_JLO)=TRANSFER(IBSET(1,31),0.)
          END IF
          WHOLD=.TRUE.
        CASE(7)
          RVH(3)=NEAREST(100.,1.)
          WHOLD=.TRUE.
        CASE(8)
          RVH(3)=100.
          UAF(15:16,3)=100.
          WHV=100.
        CASE(9)
          REH=1.E10
          FSM=0.
          WHOLD=.TRUE.
        CASE(10)
          IF(SMPI_RANK.EQ.SMPI_SIZE-1) ELF(:,SMPI_JLO)=2.
          WHOLD=.TRUE.
        CASE(11)
          RVH(JMM1)=TRANSFER(IBSET(0,31),0.)
          WHV=TRANSFER(IBSET(0,31),0.)
        CASE(12)
          RVH(5)=6.
          UAF(8:9,5)=6.
! A wrong local-per-rank speed rescan would trap; only original winner J=5 is legal.
          IF(SMPI_RANK.EQ.SMPI_SIZE-1) UAF(IM,JMM1)=TRANSFER(INT(Z'7F800001'),0.)
          WHV=6.
        CASE(13:14)
          REH=4.
          ELF=4.
          REH(3)=-1.
          REH(JM-1)=NEAREST(-1.,-1.)
          IF(IHCASE.EQ.14) THEN
            REH(3)=NEAREST(1.,1.)
            REH(JM-1)=1.
          END IF
          ELF(37,3)=REH(3)
          ELF(51,JM-1)=REH(JM-1)
          WHE=REH(JM-1)
          IHR=51
          JHR=JM-1
        CASE(15)
          REH(1)=TRANSFER(INT(Z'7F800001'),0.)
          WHOLD=.TRUE.
        CASE(16)
          REH(1)=TRANSFER(IBSET(1,31),0.)
          WHOLD=.TRUE.
        CASE(17)
          RVH=-1.E10
          WHOLD=.TRUE.
        CASE(18)
! An sNaN beyond the already matched first I must not be newly examined.
          ELF(2,SMPI_JLO)=TRANSFER(INT(Z'7F800001'),0.)
        END SELECT
        CALL SMPI_TRY_HEALTH(VAMAX,ELMINX,IAMAX,JAMAX,RVH,REH,OLDH)
        IF(OLDH.NEQV.WHOLD) IFH=IFH+1
        IF(WHOLD) THEN
          IF(VAMAX.NE.17..OR.TRANSFER(ELMINX,0).NE.IBSET(0,31)
     &      .OR.IAMAX.NE.-97.OR.JAMAX.NE.-98) IFH=IFH+1
        ELSE
          IF(TRANSFER(VAMAX,0).NE.TRANSFER(WHV,0).OR.
     &      TRANSFER(ELMINX,0).NE.TRANSFER(WHE,0).OR.IAMAX.NE.IHR.OR.JAMAX.NE.JHR) IFH=IFH+1
        END IF
      END DO
      CALL SMPI_INDEX_RAW(IFH)
      IF(IFH.NE.0) STOP 86
      IF(MODEH.EQ.'trap-new'.OR.MODEH.EQ.'trap-old') THEN
! This test MUST trap during the real selected speed row rescan in both
! routines. ROWVMAX deliberately bypasses the producer's own FP trap.
        RVH=0.
        REH=1.
        ELF=1.
        FSM=1.
        UAF=0.
        VAF=0.
        IF(SMPI_OWNS(JMM1)) UAF(IM,JMM1)=TRANSFER(INT(Z'7F800001'),0.)
        IF(SMPI_RANK.EQ.0) THEN
          WRITE(*,*) 'EXPECT_WINNER_SIGFPE ',TRIM(MODEH)
          FLUSH(6)
        END IF
        IF(MODEH.EQ.'trap-new') THEN
          CALL SMPI_TRY_HEALTH(VAMAX,ELMINX,IAMAX,JAMAX,RVH,REH,OLDH)
        ELSE
          CALL SMPI_VERDICT2(VAMAX,ELMINX,IAMAX,JAMAX,RVH,REH,OLDH)
        END IF
! Peers must not exit early while the selected rank is raising SIGFPE.
        JSMP=0
        CALL SMPI_INDEX_RAW(JSMP)
        STOP 87
      END IF
      JSMP=SMPI_RANK
      CALL SMPI_FINALIZE()
      IF(JSMP.EQ.0) WRITE(*,*) 'PASS_HEALTH_CONTRACT cases/ranks=',18,NRH
      END
