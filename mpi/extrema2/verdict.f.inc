
      RECURSIVE SUBROUTINE SMPI_VERDICT2(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLDTAIL2)
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM),VR2,ER2
      INTEGER REC2(4),ALL2(4,SMPI_SIZE),PICK2(3),IV2,IE2,MAG2
      LOGICAL OLDTAIL2
! Record layout: speed bits, speed row, water bits, water row. Negative
! water row is an all-rank request to execute the unmodified old tail.
! SMPI_INIT enforces default REAL32/INT32 for MPI; reject non-IEEE layout.
      OLDTAIL2=.TRUE.
      REC2=0
      VR2=-1.E10
      ER2=1.E10
      IV2=0
      IE2=0
      IF(TRANSFER(1.,0).NE.INT(Z'3F800000')) GOTO 62
      DO J=SMPI_JLO,SMPI_JHI
! Never introduce an early water comparison on sNaN, or assume DAZ
! semantics. These checks operate on bits of existing row summaries
! only. Comparisons on finite normals and signed zeros cannot trap.
        MAG2=IBCLR(TRANSFER(ROWEMIN(J),0),31)
        IF(MAG2.GE.INT(Z'7F800000').OR.
     &    (MAG2.GT.0.AND.MAG2.LT.INT(Z'00800000'))) GOTO 62
        IF(ROWEMIN(J).LT.ER2) THEN
          ER2=ROWEMIN(J)
          IE2=J
        END IF
        IF(J.EQ.JM) CYCLE
        MAG2=IBCLR(TRANSFER(ROWVMAX(J),0),31)
        IF(MAG2.GE.INT(Z'7F800000').OR.
     &    (MAG2.GT.0.AND.MAG2.LT.INT(Z'00800000'))) GOTO 62
        IF(ROWVMAX(J).GE.VR2.AND.ROWVMAX(J).GE.0.) THEN
          VR2=ROWVMAX(J)
          IV2=J
        END IF
      END DO
      REC2=[TRANSFER(VR2,0),IV2,TRANSFER(ER2,0),IE2]
      GOTO 63
   62 REC2(4)=-1
   63 CALL SMPI_CHECKREC_RAW(REC2,ALL2)
      IF(ANY(ALL2(4,:).LT.0)) RETURN
      OLDTAIL2=.FALSE.
      VAMAX=-1.E10
      JVMAX=0
      DO IR2=1,SMPI_SIZE
        IF(ALL2(2,IR2).EQ.0) CYCLE
        VR2=TRANSFER(ALL2(1,IR2),0.)
! Ranks in ascending order own strictly ascending contiguous J intervals.
        IF(VR2.GE.VAMAX) THEN
          VAMAX=VR2
          JVMAX=ALL2(2,IR2)
        END IF
      END DO
      PICK2=0
      IF(JVMAX.GT.0) THEN
        J=JVMAX
        IF(SMPI_OWNS(J)) THEN
        DO I=IMM1,1,-1
      VMAXTMP=SQRT((0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2)
          IF(VMAXTMP.EQ.VAMAX) THEN
            PICK2(1)=I
            EXIT
          END IF
        END DO
        END IF
      END IF
      JWATER=0
      IF(.NOT.(VAMAX.GT.100.)) THEN
        ELMINX=1.E10
        DO IR2=1,SMPI_SIZE
          IF(ALL2(4,IR2).EQ.0) CYCLE
          ER2=TRANSFER(ALL2(3,IR2),0.)
          IF(ER2.LT.ELMINX) THEN
            JWATER=ALL2(4,IR2)
            ELMINX=ER2
          END IF
        END DO
        IF(JWATER.GT.0) THEN
          J=JWATER
          IF(SMPI_OWNS(J)) THEN
          DO I=1,IM
            IF(FSM(I,J).NE.0..AND.ELF(I,J).EQ.ELMINX) THEN
              PICK2(2)=I
              PICK2(3)=TRANSFER(ELF(I,J),0)
              EXIT
            END IF
          END DO
          END IF
        END IF
      END IF
      CALL SMPI_CHECKPICKS_RAW(PICK2)
      IF(PICK2(1).GT.0) THEN
        IAMAX=PICK2(1)
        JAMAX=JVMAX
      END IF
      IF(PICK2(2).GT.0) THEN
        IAMAX=PICK2(2)
        JAMAX=JWATER
        ELMINX=TRANSFER(PICK2(3),0.)
      END IF
      END
