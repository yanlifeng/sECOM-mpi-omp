      PROGRAM TEST_EXTREMA_MPI
! MPI/OMP QA only; not a solver or timing sample.
! argv: expected MPI ranks, expected OpenMP threads/rank.
      USE MPI_F08, ONLY: MPI_COMM_WORLD,MPI_INTEGER,MPI_SUM,
     &                   MPI_ALLREDUCE,MPI_SUCCESS
      USE OMP_LIB, ONLY: OMP_GET_NUM_THREADS
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM),WANTV,WANTZ
      INTEGER NRTEST,NTTEST,NTUSED,JSEAM,JLATE,IE_EXPECT,JE_EXPECT
      INTEGER IERRTEST,IFAILTEST,ITOTALERR,ICASE
      CHARACTER*16 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) NRTEST
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) NTTEST
      CALL SMPI_INIT(IM,JM)
      IFAILTEST=0
      JSEAM=JM/SMPI_SIZE+MIN(1,MOD(JM,SMPI_SIZE))
      IF(SMPI_SIZE.EQ.1) JSEAM=JM/2
      JLATE=JM/2+2
      DO ICASE=1,9
! FSM is a replicated static field. All dynamic non-owned rows, including
! incoming ghost slots, begin poisoned. Only the true owner seeds a pulse.
        FSM=1.
        IF(ICASE.EQ.5) FSM=0.
        ROWVMAX=9000.+SMPI_RANK
        ROWEMIN=-9000.-SMPI_RANK
        DUMWAD=3.
        DVMWAD=4.
        UAF=400.+SMPI_RANK
        VAF=800.+SMPI_RANK
        ELF=-1234.-SMPI_RANK
        DO J=SMPI_JLO,SMPI_JHI
          DUMWAD(:,J)=1.
          DVMWAD(:,J)=1.
          UAF(:,J)=0.
          VAF(:,J)=0.
          ELF(:,J)=-1.
          IF(ICASE.EQ.1.OR.ICASE.EQ.2) THEN
            IF(J.EQ.5) UAF(7:8,J)=128.
            IF(J.EQ.5.OR.J.EQ.6) VAF(7,J)=.046875
            IF(ICASE.EQ.1.AND.J.EQ.JM-2) UAF(3:4,J)=128.
            IF(ICASE.EQ.2.AND.J.EQ.5) UAF(120:121,J)=128.
          ELSEIF(ICASE.EQ.3.OR.ICASE.EQ.4) THEN
            IF(J.EQ.5) UAF(7:8,J)=100.
            IF(J.EQ.5.OR.J.EQ.6) THEN
              IF(ICASE.EQ.3) VAF(7,J)=.015625
              IF(ICASE.EQ.4) VAF(7,J)=.03125
            ENDIF
          ELSEIF(ICASE.EQ.6.OR.ICASE.EQ.7) THEN
            ELF(:,J)=1.
            IF(J.EQ.JLATE.OR.J.EQ.JM-1) THEN
              ELF(13,J)=-0.
              ELF(14,J)=+0.
              IF(ICASE.EQ.7) THEN
                ELF(13,J)=+0.
                ELF(14,J)=-0.
              ENDIF
            ENDIF
          ELSEIF(ICASE.EQ.8) THEN
! Only the +J ghost's (VAF=100,DVMWAD=2) makes the seam's speed
! 120; the neighboring row is exactly 100, not an abort.
            IF(J.EQ.JSEAM) VAF(7,J)=40.
            IF(J.EQ.JSEAM+1) THEN
              VAF(7,J)=100.
              DVMWAD(7,J)=2.
            ENDIF
          ELSEIF(ICASE.EQ.9) THEN
! Equal minima on first and last owner: first global J outranks lower I.
            IF(J.EQ.3) ELF(23,J)=-2.
            IF(J.EQ.JM-2) ELF(2,J)=-2.
          ENDIF
        ENDDO
        ELMINX=123.
        IAMAX=-99
        JAMAX=-99
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
        NTUSED=OMP_GET_NUM_THREADS()
!$OMP END MASTER
        CALL SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
!$OMP END PARALLEL
        WANTV=0.
        WANTZ=123.
        SELECT CASE(ICASE)
        CASE(1)
          WANTV=128.
          IE_EXPECT=3
          JE_EXPECT=JM-2
        CASE(2)
          WANTV=128.
          IE_EXPECT=120
          JE_EXPECT=5
        CASE(3)
          WANTV=100.
          WANTZ=-1.
          IE_EXPECT=1
          JE_EXPECT=1
        CASE(4)
          WANTV=NEAREST(100.,1.)
          IE_EXPECT=7
          JE_EXPECT=5
        CASE(5)
          WANTZ=1.E10
          IE_EXPECT=IMM1
          JE_EXPECT=JMM1
        CASE(6,7)
          WANTZ=-0.
          IF(ICASE.EQ.7) WANTZ=+0.
          IE_EXPECT=13
          JE_EXPECT=JLATE
        CASE(8)
          WANTV=120.
          IE_EXPECT=7
          JE_EXPECT=JSEAM
        CASE(9)
          WANTZ=-2.
          IE_EXPECT=23
          JE_EXPECT=3
        END SELECT
        IF(TRANSFER(VAMAX,0).NE.TRANSFER(WANTV,0).OR.
     &     TRANSFER(ELMINX,0).NE.TRANSFER(WANTZ,0).OR.
     &     IAMAX.NE.IE_EXPECT.OR.JAMAX.NE.JE_EXPECT.OR.
     &     NTUSED.NE.NTTEST.OR.SMPI_SIZE.NE.NRTEST) THEN
          WRITE(*,*) 'FAIL rank/case/vi/zi/i/j/ranks/threads:',
     &     SMPI_RANK,ICASE,TRANSFER(VAMAX,0),TRANSFER(ELMINX,0),
     &     IAMAX,JAMAX,SMPI_SIZE,NTUSED
          IFAILTEST=IFAILTEST+1
        ENDIF
      ENDDO
! Independent assertion reduction; never silently accept rank zero only.
      CALL MPI_ALLREDUCE(IFAILTEST,ITOTALERR,1,MPI_INTEGER,MPI_SUM,
     &                   MPI_COMM_WORLD,IERRTEST)
      IF(IERRTEST.NE.MPI_SUCCESS) ITOTALERR=-1
      IF(SMPI_RANK.EQ.0.AND.ITOTALERR.EQ.0)
     & WRITE(*,*) 'PASS_MPI_EXTREMA cases/ranks/threads=',9,NRTEST,NTTEST
      CALL SMPI_FINALIZE()
      IF(ITOTALERR.NE.0) STOP 42
      END
