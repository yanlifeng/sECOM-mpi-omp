
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
