"""Read-only r014 science audit and actual non-MPI serial/OMP regressions."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess as sp
import tempfile

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];SCI=ROOT/'source'
FREEZE='46f5da3de7d5904ed20e902e1c1ba501f2ca40fa090e11c07a8d7754a6b9d1a4'
BASE='c371fa3';V46='1b58f33639538c2d214f14a8f9a575ceffed5267'
FLAGS='-O3 -std=legacy -ffixed-line-length-132 -ffpe-trap=invalid,zero,overflow -fno-fast-math -ffp-contract=off -fcheck=bounds'.split()
def sha(b):return hashlib.sha256(b).hexdigest()
def mapping():
    m={p.name:sha(p.read_bytes()) for p in SCI.iterdir() if p.is_file()};assert len(m)==107
    h=sha(json.dumps(m,sort_keys=True,separators=(',',':')).encode());assert h==FREEZE,h
    return m,h
def git(ref,path):return sp.check_output(['git','-C',str(ROOT),'show',f'{ref}:source/{path}'])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fc',default='gfortran');p.add_argument('--result',type=Path);a=p.parse_args()
    pre,h=mapping();oldpub=git(BASE,'wandd_kky.f').decode();now=(SCI/'wandd_kky.f').read_text().split('      RECURSIVE SUBROUTINE SECOM_WAD_COLOR_TEAM',1)[0]
    assert [s.rstrip() for s in oldpub.rstrip().splitlines()]==[s.rstrip() for s in now.rstrip().splitlines()]
    for n in ('comdeck','Makefile.portable','archive.f'):assert git(BASE,n)==(SCI/n).read_bytes()
    old46=git(V46,'wandd_kky.f').decode().replace('WANDD_KKY','V46_WANDD_KKY').replace('WANDD_TEAM','V46_WANDD_TEAM')
    fixture=(ROOT/'validationmpi/wad_owned_driver.f').read_text();init=fixture[fixture.index('      SUBROUTINE WAD_FIXTURE(M)'):]
    prog='''      PROGRAM DEFAULTS
      INCLUDE 'comdeck'
      INTEGER NT,SAW
      INTEGER,ALLOCATABLE::EXPECTED(:,:,:)
      CHARACTER*8 ARG
      LOGICAL SMPI_ARC_REQUESTED
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) NT
      IF(SMPI_OUTER.OR.SMPI_ARC_REQUESTED().OR.SMPI_LOCALBC) STOP 51
      IF(.NOT.SMPI_KEEP_UA.OR..NOT.SMPI_KEEP_VA) STOP 52
      IF(SMPI_SIZE.NE.1.OR.SMPI_RANK.NE.0.OR.SMPI_JLO.NE.1.OR.SMPI_JHI.NE.HUGE(0)) STOP 53
      CALL SMPI_INIT(IM,JM)
      CALL SECOM_MPI_BC_PLAN()
      IF(SMPI_OUTER.OR.SMPI_LOCALBC.OR.SMPI_JLO.NE.1.OR.SMPI_JHI.NE.JM) STOP 54
      ALLOCATE(EXPECTED(IM,JM,8))
      SAW=0
!$OMP PARALLEL REDUCTION(+:SAW)
      SAW=SAW+1
      IF(SMPI_OUTER.OR.SMPI_LOCALBC.OR..NOT.SMPI_KEEP_UA.OR..NOT.SMPI_KEEP_VA) STOP 55
!$OMP END PARALLEL
      IF(SAW.NE.NT) STOP 56
      DO M=1,8
        CALL WAD_FIXTURE(M)
        CALL V46_WANDD_KKY(EL,D)
        DO J=1,JM
          DO I=1,IM
            EXPECTED(I,J,:)=[TRANSFER(UA(I,J),0),TRANSFER(UAF(I,J),0),TRANSFER(DUMWAD(I,J),0),
     &       TRANSFER(VA(I,J),0),TRANSFER(VAF(I,J),0),TRANSFER(DVMWAD(I,J),0),TRANSFER(ET(I,J),0),TRANSFER(ETF(I,J),0)]
          ENDDO
        ENDDO
        CALL WAD_FIXTURE(M)
        CALL WANDD_KKY(EL,D)
        DO J=1,JM
          DO I=1,IM
            IF(ANY(EXPECTED(I,J,:).NE.[TRANSFER(UA(I,J),0),TRANSFER(UAF(I,J),0),TRANSFER(DUMWAD(I,J),0),
     &       TRANSFER(VA(I,J),0),TRANSFER(VAF(I,J),0),TRANSFER(DVMWAD(I,J),0),TRANSFER(ET(I,J),0),TRANSFER(ETF(I,J),0)])) STOP 57
          ENDDO
        ENDDO
      ENDDO
      SMPI_KEEP_UA=.FALSE.
      SMPI_KEEP_VA=.FALSE.
      SMPI_LOCALBC=.TRUE.
      CALL SMPI_FINALIZE()
      IF(.NOT.SMPI_KEEP_UA.OR..NOT.SMPI_KEEP_VA.OR.SMPI_LOCALBC.OR.SMPI_OUTER) STOP 58
      IF(SMPI_SIZE.NE.1.OR.SMPI_RANK.NE.0.OR.SMPI_JLO.NE.1.OR.SMPI_JHI.NE.HUGE(0)) STOP 59
      PRINT *, 'PASS WAD/defaults',IM,JM,SAW
      END
      SUBROUTINE PUTCDF(T)
      STOP 60
      END
      SUBROUTINE TSCDF(T,I)
      STOP 60
      END
'''+init
    rows={}
    with tempfile.TemporaryDirectory(prefix='smpi-r014-defaults-') as t:
        d=Path(t);(d/'driver.f').write_text(prog);(d/'v46.f').write_text(old46)
        for mode in ['serial','omp']:
            f=FLAGS+(['-fopenmp'] if mode=='omp' else [])
            sp.run([a.fc,*f,'-I'+str(SCI),'-J.','-cpp','-c',str(SCI/'secom_mpi.f90'),'-o','rt.o'],cwd=d,check=True)
            sp.run([a.fc,*f,'-I'+str(SCI),str(SCI/'archive.f'),str(SCI/'mpi_boundary.f'),str(SCI/'wandd_kky.f'),'v46.f','driver.f','rt.o','-o','defaults'],cwd=d,check=True)
            sp.run([a.fc,*f,'-I'+str(SCI),str(SCI/'baro_check.f'),str(ROOT/'validation/check_extrema.f'),'rt.o','-o','extrema'],cwd=d,check=True)
            for nt in ([1] if mode=='serial' else [1,3]):
                en=dict(os.environ,OMP_NUM_THREADS=str(nt),OMP_DYNAMIC='FALSE');r={}
                for name in ['defaults','extrema']:
                    q=sp.run([str(d/name),str(nt)],cwd=d,env=en,stdout=sp.PIPE,stderr=sp.STDOUT,timeout=35)
                    assert q.returncode==0 and b'PASS' in q.stdout,(name,nt,q.stdout);r[name]=q.stdout.decode().strip()
                rows[f'{mode}{nt}']=r;print(mode,nt,r,flush=True)
    post,h2=mapping();assert pre==post
    changed=[n for n,v in pre.items() if sha(git(BASE,n))!=v]
    out=dict(source_pre=h,source_post=h2,nsource=len(pre),changed_since_r004b=sorted(changed),public_wad_equal_ignoring_trailing_whitespace=True,
             shape=[372,479],compiler=sp.check_output([a.fc,'--version']).decode().splitlines()[0],common_flags=FLAGS,
             original_extrema_fixture=sha((ROOT/'validation/check_extrema.f').read_bytes()),wad_fixture=sha(init.encode()),cases=rows)
    if a.result:
        with a.result.open('x') as o:json.dump(out,o,indent=2)
if __name__=='__main__':main()
