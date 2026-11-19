      RECURSIVE SUBROUTINE SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
! Independent callers may supply poisoned/stale ghosts: always refresh.
      CALL SECOM_CHECK_READY_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,
     &                            ROWVMAX,ROWEMIN,.FALSE.)
      END
      RECURSIVE SUBROUTINE SECOM_CHECK_READY_TEAM(VAMAX,ELMINX,
     &             IAMAX,JAMAX,ROWVMAX,ROWEMIN,HALO_READY)
      LOGICAL HALO_READY
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
      LOGICAL OLDTAIL2
! Wet/dry owns and commits all face changes first. Only these two
! dynamic J+1 inputs cross a row owner for the speed calculation.
      IF(SMPI_SIZE.GT.1.AND..NOT.HALO_READY) THEN
        CALL SMPI_HALO(DVMWAD,IM,JM)
        CALL SMPI_HALO(VAF,IM,JM)
      ENDIF
! First reduce only values, making the per-cell loop vectorizable.
! Recover indices only on the two winning rows, preserving the legacy
! tie order and speed-before-minimum abort behavior (including all-land).
!$OMP DO SCHEDULE(STATIC)
      DO J=SMPI_FIRST(1),SMPI_LAST(JMM1)
        VROWMAX=-1.E10
!$OMP SIMD REDUCTION(MAX:VROWMAX) PRIVATE(VMAXTMP)
        DO I=1,IMM1
! SQRT is monotone on a finite nonnegative sum of squares. Reduce the
! original rounded squared norm, then round the root once per row.
! The final scan below still compares rounded roots: distinct squared
! norms may tie after SQRT, even in different rows or at the abort limit.
      VMAXTMP=(0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2
          IF(VMAXTMP.GT.VROWMAX) VROWMAX=VMAXTMP
        ENDDO
        ROWVMAX(J)=VROWMAX
        IF(VROWMAX.GE.0.) ROWVMAX(J)=SQRT(VROWMAX)
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO J=SMPI_FIRST(1),SMPI_LAST(JM)
! SMPI_ROWS transfers JM entries of each vector; JM has no speed cell.
        IF(SMPI_SIZE.GT.1.AND.J.EQ.JM) ROWVMAX(J)=-1.E10
        ZROWMIN=1.E10
!$OMP SIMD REDUCTION(MIN:ZROWMIN)
        DO I=1,IM
          IF(FSM(I,J).NE.0..AND.ELF(I,J).LT.ZROWMIN)
     &      ZROWMIN=ELF(I,J)
        ENDDO
        ROWEMIN(J)=ZROWMIN
      ENDDO
!$OMP END DO
      IF(SMPI_SIZE.GT.1) THEN
!$OMP SINGLE
        CALL SMPI_TRY_HEALTH(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLDTAIL2)
        IF(OLDTAIL2) THEN
        CALL SMPI_VERDICT2(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLDTAIL2)
        END IF
!$OMP END SINGLE COPYPRIVATE(OLDTAIL2)
        IF(.NOT.OLDTAIL2) RETURN
      END IF
      IF(SMPI_SIZE.GT.1) CALL SMPI_ROWS(ROWVMAX,ROWEMIN,JM)
!$OMP SINGLE
      VAMAX=-1.E10
      JVMAX=0
      DO J=1,JMM1
        IF(ROWVMAX(J).GE.VAMAX.AND.ROWVMAX(J).GE.0.) THEN
          JVMAX=J
          VAMAX=ROWVMAX(J)
        ENDIF
      ENDDO
      IF(JVMAX.NE.0) THEN
        J=JVMAX
        IPICK=0
        IF(SMPI_OWNS(J)) THEN
        DO I=IMM1,1,-1
      VMAXTMP=SQRT((0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2)
          IF(VMAXTMP.EQ.VAMAX) THEN
            IPICK=I
            EXIT
          ENDIF
        ENDDO
        ENDIF
! A remote row's dense data need not exist locally. Only the owner
! reports its descending-I winner. Keep old coordinates if no match.
        IF(SMPI_SIZE.GT.1) CALL SMPI_INDEX_RAW(IPICK)
        IF(IPICK.GT.0) THEN
          IAMAX=IPICK
          JAMAX=JVMAX
        ENDIF
      ENDIF
      IF(.NOT.(VAMAX.GT.100.)) THEN
      ELMINX=1.E10
      JWATER=0
      DO J=1,JM
        IF(ROWEMIN(J).LT.ELMINX) THEN
          JWATER=J
          ELMINX=ROWEMIN(J)
        ENDIF
      ENDDO
      IF(JWATER.NE.0) THEN
        J=JWATER
        IPICK=0
        IF(SMPI_OWNS(J)) THEN
        DO I=1,IM
          IF(FSM(I,J).NE.0..AND.ELF(I,J).EQ.ELMINX) THEN
            IPICK=I
            ELMINX=ELF(I,J)
            EXIT
          ENDIF
        ENDDO
        ENDIF
        IF(SMPI_SIZE.GT.1) CALL SMPI_INDEX_RAW(IPICK)
        IF(IPICK.GT.0) THEN
          IAMAX=IPICK
          JAMAX=JWATER
! Preserve the *matched cell*, including its original sign of zero.
          IF(SMPI_SIZE.GT.1) CALL SMPI_RBCAST_ROW_RAW(ELMINX,JWATER)
        ENDIF
      ENDIF
      ENDIF
!$OMP END SINGLE
      END

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

      RECURSIVE SUBROUTINE SMPI_TRY_HEALTH(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,OLD1)
      USE ISO_FORTRAN_ENV, ONLY: INT64
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
      LOGICAL OLD1
      INTEGER(INT64) KEYS1(3),BASE1,TOTAL1,SCORE1,TOP1,LOC1
      PARAMETER(BASE1=2147483648_INT64)
      INTEGER MV1,MW1,MF1,MZ1,BV1,BW1,BF1,BZ1,JW1,JVEL1
      INTEGER, VOLATILE :: IPSPEED1
! BASE1 is 2^31. Positive high water score uses at most 32 bits;
! low bits carry reverse global cell priority and the actual cell sign.
      OLD1=.TRUE.
      KEYS1=0_INT64
      TOTAL1=INT(IM,INT64)*INT(JM,INT64)
      IF(TOTAL1.GT.1073741824_INT64.OR.JM.GE.1073741824) GOTO 71
      IF(TRANSFER(1.,0).NE.INT(Z'3F800000')) GOTO 71
      JW1=0
      BW1=TRANSFER(1.E10,0)
      TOP1=BASE1-INT(BW1,INT64)
      DO J=SMPI_JLO,SMPI_JHI
        BZ1=TRANSFER(ROWEMIN(J),0)
        MZ1=IBCLR(BZ1,31)
        IF(MZ1.GE.INT(Z'7F800000').OR.
     &    (MZ1.GT.0.AND.MZ1.LT.INT(Z'00800000'))) GOTO 71
        IF(BTEST(BZ1,31)) THEN
          SCORE1=BASE1+INT(MZ1,INT64)
        ELSE
          SCORE1=BASE1-INT(MZ1,INT64)
        END IF
! Higher sortable water score means smaller numerical value. Signed
! zeros share a score, so the legacy first J then first I decides.
        IF(SCORE1.GT.TOP1) THEN
          TOP1=SCORE1
          JW1=J
          BW1=BZ1
        END IF
        IF(J.EQ.JM) CYCLE
        BV1=TRANSFER(ROWVMAX(J),0)
        MV1=IBCLR(BV1,31)
        IF(MV1.GE.INT(Z'7F800000').OR.
     &    (MV1.GT.0.AND.MV1.LT.INT(Z'00800000'))) GOTO 71
! Skip the old negative speed sentinel, but retain -0 and its bits.
        IF(BTEST(BV1,31).AND.MV1.NE.0) CYCLE
        IF(MV1.GT.TRANSFER(100.,0)) GOTO 71
        SCORE1=INT(MV1,INT64)*BASE1+2_INT64*J+IBITS(BV1,31,1)
        KEYS1(1)=MAX(KEYS1(1),SCORE1)
      END DO
      IF(JW1.NE.0) THEN
        MW1=IBCLR(BW1,31)
        DO I=1,IM
! Inspect only integer bits. No speculative FP operation on a losing
! rank, and no FP comparison even for an sNaN in a masked cell.
          BF1=TRANSFER(FSM(I,JW1),0)
          MF1=IBCLR(BF1,31)
          BZ1=TRANSFER(ELF(I,JW1),0)
          MZ1=IBCLR(BZ1,31)
          IF(MF1.GE.INT(Z'7F800000').OR.MZ1.GE.INT(Z'7F800000')) GOTO 71
          IF((MF1.GT.0.AND.MF1.LT.INT(Z'00800000')).OR.
     &       (MZ1.GT.0.AND.MZ1.LT.INT(Z'00800000'))) GOTO 71
          IF(MF1.EQ.0) CYCLE
          IF(BZ1.NE.BW1.AND..NOT.(MW1.EQ.0.AND.MZ1.EQ.0)) CYCLE
          LOC1=INT(JW1-1,INT64)*IM+I
          KEYS1(2)=TOP1*BASE1+2_INT64*(TOTAL1-LOC1)+IBITS(BZ1,31,1)
          EXIT
        END DO
        IF(KEYS1(2).EQ.0_INT64) GOTO 71
      END IF
      GOTO 72
   71 KEYS1(3)=1_INT64
   72 CALL SMPI_CHECKKEYS_RAW(KEYS1)
! Exact older tail handles abnormal summaries, abort, no valid water,
! or a missing speed winner. Absolutely no public output has been set.
      IF(KEYS1(3).NE.0.OR.KEYS1(2).EQ.0.OR.KEYS1(1).EQ.0) RETURN
      OLD1=.FALSE.
      MV1=INT(KEYS1(1)/BASE1)
      BV1=MV1
      IF(BTEST(KEYS1(1),0)) BV1=IBSET(BV1,31)
      VAMAX=TRANSFER(BV1,0.)
      JVEL1=INT(MOD(KEYS1(1),BASE1)/2_INT64)
      IF(SMPI_OWNS(JVEL1)) THEN
        J=JVEL1
! This is the ORIGINAL global winning speed row only. The volatile
! assignment makes the successful search observable to the compiler,
! so it cannot discard any original SQRT merely because the guaranteed
! water winner overwrites the externally returned coordinates.
        IPSPEED1=0
        DO I=IMM1,1,-1
      VMAXTMP=SQRT((0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2)
          IF(VMAXTMP.EQ.VAMAX) THEN
            IPSPEED1=I
            EXIT
          END IF
        END DO
      END IF
      TOP1=KEYS1(2)/BASE1
      LOC1=TOTAL1-MOD(KEYS1(2),BASE1)/2_INT64
      MV1=INT(ABS(TOP1-BASE1))
      IF(BTEST(KEYS1(2),0)) MV1=IBSET(MV1,31)
      ELMINX=TRANSFER(MV1,0.)
      IAMAX=INT(MOD(LOC1-1_INT64,INT(IM,INT64))+1_INT64)
      JAMAX=INT((LOC1-1_INT64)/IM+1_INT64)
      END
