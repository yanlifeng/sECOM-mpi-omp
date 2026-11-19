"""Build a separate diagnostic MPI executable; never edit the timed objects."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess

FIELDS = ('DUMWAD DVMWAD EL ELB ELF ET ETF UA UAB UAF VA VAB VAF D '
          'CURV42D CBC TXBOT TYBOT FLUXUA FLUXVA TPS TRNU TRNV ADVUA ADVVA').split()
EARLY = [(1, 1), (1, 10), (17, 3), (153, 5), (300, 10)]
LATE = [(601, 2), (1199, 10), (1200, 10), (1201, 1), (2400, 7), (3600, 4), (4800, 10)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('source', type=Path)
    ap.add_argument('dest', type=Path)
    ap.add_argument('--late', action='store_true')
    args = ap.parse_args()
    src, dest = args.source.resolve(), args.dest.resolve()
    dest.mkdir()
    code = (src / 'ecom3d.f').read_text()
    anchor = '        IF(VAMAX.GT.100.) EXIT'
    assert code.count(anchor) == 1
    capture = """!$OMP SINGLE
        CALL CAPTURE_MPI_SUBSTEP(ADVUA,ADVVA,TRNU,TRNV,VAMAX,ELMINX,IAMAX,JAMAX)
!$OMP END SINGLE
"""
    (dest / 'main.f').write_text(code.replace(anchor, capture + anchor))
    frames = LATE if args.late else EARLY
    checks = [f'(INTX.EQ.{a}.AND.IEXT.EQ.{b})' for a, b in frames]
    condition = '      IF(.NOT.(' + checks[0] + '\n' + ''.join('     & .OR.' + s + '\n' for s in checks[1:]) + '     & )) RETURN\n'
    text = """      SUBROUTINE CAPTURE_MPI_SUBSTEP(ADVUA,ADVVA,TRNU,TRNV,VAMAX,ELMINX,IAMAX,JAMAX)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM),TRNU(IM,JM),TRNV(IM,JM)
      CHARACTER*32 TRNAME
      INTEGER TRUNIT
""" + condition
    text += """      IF(SMPI_RANK.NE.0) RETURN
      WRITE(TRNAME,'("transient-",I3.3,"-",I2.2,".bin")') INTX,IEXT
      OPEN(NEWUNIT=TRUNIT,FILE=TRNAME,FORM='UNFORMATTED',ACCESS='STREAM',STATUS='NEW')
      WRITE(TRUNIT) INTX,IEXT,IM,JM,VAMAX,ELMINX,IAMAX,JAMAX
"""
    text = text.replace('      IF(SMPI_RANK.NE.0) RETURN',
                        '      TRUNIT=-1\n      IF(SMPI_RANK.EQ.0) THEN')
    text += '      ENDIF\n'
    text += ''.join('      CALL CAPTURE_MPI_PLANE(TRUNIT,%s)\n' % f for f in FIELDS)
    text += """      IF(SMPI_RANK.EQ.0) CLOSE(TRUNIT)
      END
      SUBROUTINE CAPTURE_MPI_PLANE(TRUNIT,CAPFIELD)
      INCLUDE 'comdeck'
      INTEGER TRUNIT
      REAL CAPFIELD(IM,JM),BUFFER(IM,JM)
      SAVE BUFFER
! Gather a COPY of owned rows: a diagnostic must not refresh the
! scientific ghost rows and silently repair a missing halo later.
      DO J=SMPI_JLO,SMPI_JHI
        BUFFER(:,J)=CAPFIELD(:,J)
      ENDDO
      CALL SMPI_ALL_RAW(BUFFER,IM,JM)
      IF(SMPI_RANK.EQ.0) WRITE(TRUNIT) BUFFER
      END
"""
    if args.late:
        text = text.replace('I3.3', 'I4.4')
    (dest / 'capture.f').write_text(text)
    nf = '/home/user_home/ylf/kt5-licom-secom-ylf/deps/bin/nf-config'
    flags = shlex.split('-O3 -fno-fast-math -ffp-contract=off -std=legacy -ffixed-line-length-132 '
                        '-g -fbacktrace -ffpe-trap=invalid,zero,overflow -march=native -mprefer-vector-width=256 -fopenmp')
    fc = '/usr/lib64/openmpi/bin/mpifort'
    incs = shlex.split(subprocess.check_output([nf, '--fflags'], universal_newlines=True))
    libs = shlex.split(subprocess.check_output([nf, '--flibs'], universal_newlines=True))
    commands = []
    for unit in ('main', 'capture'):
        commands.append([fc, *flags, *incs, '-I' + str(src), '-c', str(dest / (unit + '.f')), '-o', str(dest / (unit + '.o'))])
    objects = sorted((src / 'build/mpi').glob('*.o'))
    assert len(objects) > 90 and any(p.name == 'secom_mpi.o' for p in objects)
    commands.append([fc, *flags, '-o', str(dest / 's--pom.x'), str(dest / 'capture.o'),
                     *[str(dest / 'main.o' if p.name == 'ecom3d.o' else p) for p in objects], *libs])
    with (dest / 'build.log').open('x') as log:
        for cmd in commands:
            log.write(' '.join(map(shlex.quote, cmd)) + '\n')
            log.flush()
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    proof = dict(kind='instrumented_not_for_timing', snapshot_mode='owned_copy_nondestructive',
                 commands=commands, frames=frames, fields=FIELDS,
                 binary_sha=hashlib.sha256((dest / 's--pom.x').read_bytes()).hexdigest())
    (dest / 'build.json').write_text(json.dumps(proof, indent=2))


if __name__ == '__main__':
    main()
