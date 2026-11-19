"""Extend the unmodified 9-case native MPI regression; leave the oracle alone."""
import argparse
from pathlib import Path
from prepare import once


def extend(orig):
    s = once(orig, '      CHARACTER*16 ARG', '      LOGICAL VETO2\n      CHARACTER*16 ARG')
    s = once(s, 'DO ICASE=1,9', 'DO ICASE=1,16')
    s = once(s, 'IF(ICASE.EQ.5) FSM=0.', 'IF(ICASE.EQ.5.OR.ICASE.EQ.14) FSM=0.')
    s = once(s, '          ELSEIF(ICASE.EQ.9) THEN', '''          ELSEIF(ICASE.EQ.10.OR.ICASE.EQ.11) THEN
            ELF(:,J)=1.
            IF(J.EQ.3) ELF(45,J)=MERGE(-0.,+0.,ICASE.EQ.10)
            IF(J.EQ.JM-2) ELF(2,J)=MERGE(+0.,-0.,ICASE.EQ.10)
          ELSEIF(ICASE.EQ.12) THEN
            IF(J.EQ.JM-2) UAF(20:21,J)=100.
            IF(J.EQ.3) ELF(45,J)=-2.
          ELSEIF(ICASE.EQ.13.OR.ICASE.EQ.16) THEN
            IF(J.EQ.JM-2) UAF(20:21,J)=NEAREST(100.,1.)
            IF(J.EQ.3) ELF(45,J)=-99.
          ELSEIF(ICASE.EQ.14) THEN
            IF(J.EQ.3) UAF(199:200,J)=3.
            IF(J.EQ.JM-2) UAF(16:17,J)=3.
          ELSEIF(ICASE.EQ.15) THEN
            ELF(:,J)=1.E10
          ELSEIF(ICASE.EQ.9) THEN''')
    s = once(s, '        ELMINX=123.', '''        ELMINX=123.
        IF(ICASE.EQ.16) ELMINX=TRANSFER(IBSET(0,31),0.)''')
    s = once(s, '        CASE(9)', '''        CASE(10,11)
          WANTZ=MERGE(-0.,+0.,ICASE.EQ.10)
          IE_EXPECT=45
          JE_EXPECT=3
        CASE(12)
          WANTV=100.
          WANTZ=-2.
          IE_EXPECT=45
          JE_EXPECT=3
        CASE(13,16)
          WANTV=NEAREST(100.,1.)
          IF(ICASE.EQ.16) WANTZ=TRANSFER(IBSET(0,31),0.)
          IE_EXPECT=20
          JE_EXPECT=JM-2
        CASE(14)
          WANTV=3.
          WANTZ=1.E10
          IE_EXPECT=16
          JE_EXPECT=JM-2
        CASE(15)
          WANTZ=1.E10
          IE_EXPECT=IMM1
          JE_EXPECT=JMM1
        CASE(9)''')
    s = once(s, '! Independent assertion reduction;', '''! Deliberately inject exceptional ROW SUMMARIES, not solver inputs. The
! optimized tail must request the old path without touching ANY outputs
! or raising invalid, even if the unusual bits belong to another rank.
      DO ICASE=1,4
        ROWVMAX=101.
        ROWEMIN=1.
        IF(SMPI_RANK.EQ.SMPI_SIZE-1) THEN
          J=SMPI_JLO
          IF(ICASE.EQ.1) ROWEMIN(J)=TRANSFER(INT(Z'7F800001'),0.)
          IF(ICASE.EQ.2) ROWEMIN(J)=TRANSFER(IBSET(1,31),0.)
          IF(ICASE.EQ.3) ROWVMAX(J)=TRANSFER(INT(Z'7F800000'),0.)
          IF(ICASE.EQ.4) ROWVMAX(J)=TRANSFER(INT(Z'7F800001'),0.)
        END IF
        VAMAX=5.
        ELMINX=-11.
        IAMAX=-21
        JAMAX=-22
        CALL SMPI_VERDICT2(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,VETO2)
        IF(.NOT.VETO2.OR.VAMAX.NE.5..OR.ELMINX.NE.-11..OR.
     &     IAMAX.NE.-21.OR.JAMAX.NE.-22) IFAILTEST=IFAILTEST+1
      END DO
! Independent assertion reduction;''')
    return once(s, "'PASS_MPI_EXTREMA cases/ranks/threads=',9,NRTEST,NTTEST",
                "'PASS_EXTREMA2 cases/veto/ranks/threads=',16,4,NRTEST,NTTEST")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('original', type=Path)
    ap.add_argument('out', type=Path)
    args = ap.parse_args()
    data = extend(args.original.read_text())
    with args.out.open('x') as f:
        f.write(data)
