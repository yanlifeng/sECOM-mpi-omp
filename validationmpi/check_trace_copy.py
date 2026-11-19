"""Compile actual generated early/late captures; catch ghost repair invisible in the frame.

Local-only data-motion test. A deterministic substitute for the MPI allgather
provides the other rank's owned bits and overwrites every nonowned row. This
tests what is passed to the collective, not MPI/network correctness.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
from unittest.mock import patch

HERE = Path(__file__).resolve().parents[1]
FIELDS = ('DUMWAD DVMWAD EL ELB ELF ET ETF UA UAB UAF VA VAB VAF D '
          'CURV42D CBC TXBOT TYBOT FLUXUA FLUXVA TPS TRNU TRNV ADVUA ADVVA').split()
INC = """      INTEGER IM,JM,INTX,IEXT,SMPI_RANK,SMPI_JLO,SMPI_JHI,NCALL
      PARAMETER(IM=5,JM=7)
      COMMON /QCONTROL/ INTX,IEXT,SMPI_RANK,SMPI_JLO,SMPI_JHI,NCALL
      REAL QSCI(IM,JM,25)
      COMMON /QPLANES/ QSCI
""" + ''.join(f'      REAL {f}(IM,JM)\n      EQUIVALENCE({f}(1,1),QSCI(1,1,{k}))\n'
              for k, f in enumerate(FIELDS[:21], 1))
DRIVER = """      PROGRAM CHECK_CAPTURE
      INCLUDE 'comdeck'
      INTEGER QBITS,BAD,IA,IB
      CHARACTER*6 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) SMPI_RANK
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) IA
      CALL GET_COMMAND_ARGUMENT(3,ARG)
      READ(ARG,*) IB
      SMPI_JLO=1
      SMPI_JHI=4
      IF(SMPI_RANK.EQ.1) THEN
        SMPI_JLO=5
        SMPI_JHI=7
      ENDIF
      NCALL=0
      DO K=1,25
      DO J=1,JM
      DO I=1,IM
        QSCI(I,J,K)=TRANSFER(QBITS(K,I,J,J.LT.SMPI_JLO.OR.J.GT.SMPI_JHI),0.)
      ENDDO
      ENDDO
      ENDDO
      DO L=1,2
        INTX=2
        IEXT=3
        IF(L.EQ.2) THEN
          INTX=IA
          IEXT=IB
        ENDIF
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP SINGLE
        CALL CAPTURE_MPI_SUBSTEP(QSCI(1,1,24),QSCI(1,1,25),QSCI(1,1,22),QSCI(1,1,23),3.25,TRANSFER(IBSET(0,31),0.),3,5)
!$OMP END SINGLE
!$OMP END PARALLEL
        IF(NCALL.NE.25*(L-1)) STOP 91
      ENDDO
      BAD=0
      DO K=1,25
      DO J=1,JM
      DO I=1,IM
        IF(TRANSFER(QSCI(I,J,K),0).NE.QBITS(K,I,J,J.LT.SMPI_JLO.OR.J.GT.SMPI_JHI)) BAD=BAD+1
      ENDDO
      ENDDO
      ENDDO
      PRINT '(A,I0)', 'MUTATED ',BAD
      IF(BAD.NE.0) STOP 93
      END
      SUBROUTINE SMPI_ALL_RAW(G,NX,NY)
      INCLUDE 'comdeck'
      REAL G(NX,NY)
      INTEGER QBITS
      NCALL=NCALL+1
      IF(NX.NE.IM.OR.NY.NE.JM.OR.NCALL.GT.25) STOP 94
      DO J=1,NY
      DO I=1,NX
        IF(J.GE.SMPI_JLO.AND.J.LE.SMPI_JHI) THEN
          IF(TRANSFER(G(I,J),0).NE.QBITS(NCALL,I,J,.FALSE.)) STOP 95
        ELSE
          G(I,J)=TRANSFER(QBITS(NCALL,I,J,.FALSE.),0.)
        ENDIF
      ENDDO
      ENDDO
      END
      INTEGER FUNCTION QBITS(K,I,J,GHOST)
      LOGICAL GHOST
      INTEGER PAY
      PAY=K*64+J*8+I
      IF(GHOST) THEN
        QBITS=IOR(IBSET(INT(Z'7FC00000'),31),PAY)
      ELSE
        SELECT CASE(MOD(K+I+J,5))
        CASE(0)
          QBITS=IBSET(0,31)
        CASE(1)
          QBITS=IOR(INT(Z'7FC00000'),PAY)
        CASE(2)
          QBITS=IOR(INT(Z'7F800000'),PAY)
        CASE(3)
          QBITS=PAY
        CASE DEFAULT
          QBITS=IOR(INT(Z'3F800000'),PAY)
        END SELECT
      ENDIF
      END
"""


def oracle(outer, inner):
    data = bytearray(struct.pack('<4i2f2i', outer, inner, 5, 7, 3.25, -0.0, 3, 5))
    for k in range(1, 26):
        for j in range(1, 8):
            for i in range(1, 6):
                p = k*64 + j*8 + i
                x = (0x80000000, 0x7fc00000 | p, 0x7f800000 | p, p, 0x3f800000 | p)[(k+i+j) % 5]
                data.extend(struct.pack('<I', x))
    return data


def test(tmp, src, fc):
    path = HERE/'mpi/build_trace.py'
    spec = importlib.util.spec_from_file_location('trace_builder', path)
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    inc = tmp/'mini'
    inc.mkdir()
    (inc/'comdeck').write_text(INC)
    (inc/'driver.f').write_text(DRIVER)
    flags = [fc, '-O3', '-fno-fast-math', '-ffp-contract=off', '-ffpe-trap=invalid,zero,overflow',
             '-std=legacy', '-ffixed-line-length-132', '-fopenmp', '-fcheck=bounds']
    report = dict(kind='local_mock_collective_not_solver', builder_sha=hashlib.sha256(path.read_bytes()).hexdigest(),
                  compiler=subprocess.check_output([fc, '--version'], text=True).splitlines()[0],
                  results=[])
    class Rendered(Exception):
        pass
    for mode, oi, ii in [('early', 153, 5), ('late', 3600, 4)]:
        d = tmp/mode
        with patch.object(sys, 'argv', ['trace', str(src), str(d)]+(['--late'] if mode=='late' else [])), \
             patch.object(builder.subprocess, 'check_output', side_effect=Rendered):
            try:
                builder.main()
            except Rendered:
                pass
        cap = d/'capture.f'
        subprocess.run([*flags, '-I'+str(src), '-fsyntax-only', str(cap)], cwd=d, check=True)
        txt = cap.read_text()
        anchor = '      CALL SMPI_ALL_RAW(BUFFER,IM,JM)'
        assert txt.count(anchor) == 1
        mutant = d/'mutant.f'
        mutant.write_text(txt.replace(anchor, '      CALL SMPI_ALL_RAW(CAPFIELD,IM,JM)\n      BUFFER=CAPFIELD'))
        expected = oracle(oi,ii)
        for kind, ff in [('copy',cap),('inplace_negative',mutant)]:
            exe = d/(kind+'.x')
            subprocess.run([*flags, '-I'+str(inc), str(ff), str(inc/'driver.f'), '-o', str(exe)], check=True)
            for rank in (0,1):
                for nth in (1,3):
                    run = d/f'{kind}-{rank}-{nth}'
                    run.mkdir()
                    got = subprocess.run([str(exe),str(rank),str(oi),str(ii)], cwd=run,
                                         env={**os.environ,'OMP_NUM_THREADS':str(nth),'OMP_DYNAMIC':'FALSE'},
                                         capture_output=True, text=True)
                    nbad = (0 if kind=='copy' else (3,4)[rank]*5*25)
                    if got.returncode != (0 if kind=='copy' else 93) or 'MUTATED '+str(nbad) not in got.stdout:
                        raise AssertionError((got.args, got.returncode, got.stdout, got.stderr))
                    files = sorted(run.glob('transient-*.bin'))
                    assert [p.name for p in files] == ([f'transient-{oi}-{ii:02}.bin'] if rank==0 else [])
                    assert not files or files[0].read_bytes() == expected
                    report['results'].append(dict(mode=mode,variant=kind,rank=rank,omp=nth,
                        mutated_words=nbad, frame_sha=hashlib.sha256(expected).hexdigest() if files else None))
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path, default=HERE/'source')
    p.add_argument('--fc', default=shutil.which('gfortran-8') or 'gfortran')
    p.add_argument('--out', type=Path)
    a=p.parse_args()
    with tempfile.TemporaryDirectory(prefix='trace-copy-') as t:
        result=test(Path(t),a.source.resolve(),a.fc)
    if a.out:
        with a.out.open('x') as f:
            json.dump(result, f, indent=2)
    print('PASS', len(result['results']), 'cases; inplace mutated 375/500 words but identical frame')
