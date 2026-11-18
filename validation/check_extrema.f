      PROGRAM TEST_EXTREMA
! Small deterministic regression for the production safety diagnostic.
! Calling inside a team is required: the routine uses orphaned OMP DO.
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
      DO M=1,5
        DUMWAD=1.
        DVMWAD=1.
        FSM=1.
        UAF=0.
        VAF=0.
        ELF=-1.
        ELMINX=123.
        IAMAX=-99
        JAMAX=-99
        IF(M.LE.2) THEN
! sqrt(128**2+0.046875**2) and sqrt(128**2) round to the same
! binary32 speed; selecting coordinates by squared speed is wrong.
          UAF(7:8,5)=128.
          VAF(7,5:6)=.046875
          IF(M.EQ.1) UAF(3:4,JM-2)=128.
          IF(M.EQ.2) UAF(120:121,5)=128.
        ELSEIF(M.LE.4) THEN
          UAF(7:8,5)=100.
          IF(M.EQ.3) VAF(7,5:6)=.015625
          IF(M.EQ.4) VAF(7,5:6)=.03125
        ELSE
          FSM=0.
        ENDIF
!$OMP PARALLEL DEFAULT(SHARED)
        CALL SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
!$OMP END PARALLEL
        IF(M.EQ.1.AND.(VAMAX.NE.128..OR.IAMAX.NE.3.OR.
     &   JAMAX.NE.JM-2.OR.ELMINX.NE.123.)) STOP 1
        IF(M.EQ.2.AND.(VAMAX.NE.128..OR.IAMAX.NE.120.OR.
     &   JAMAX.NE.5.OR.ELMINX.NE.123.)) STOP 2
        IF(M.EQ.3.AND.(VAMAX.NE.100..OR.IAMAX.NE.1.OR.
     &   JAMAX.NE.1.OR.ELMINX.NE.-1.)) STOP 3
        IF(M.EQ.4.AND.(VAMAX.NE.NEAREST(100.,1.).OR.IAMAX.NE.7.OR.
     &   JAMAX.NE.5.OR.ELMINX.NE.123.)) STOP 4
        IF(M.EQ.5.AND.(VAMAX.NE.0..OR.IAMAX.NE.IMM1.OR.
     &   JAMAX.NE.JMM1.OR.ELMINX.NE.1.E10)) STOP 5
      ENDDO
      PRINT *, 'PASS: rounded speed ties, both abort thresholds, and all-land mask'
      END
