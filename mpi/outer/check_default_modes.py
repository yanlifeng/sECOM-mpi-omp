"""Link the real non-MPI runtime and archive, and check unconfigured defaults.

Compare focused raw history and skill records from v46 and the candidate in
both serial and OMP x3. Does not run the solver or require MPI/NetCDF.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess as sp
import tempfile
from check_archive_local import HERE,FLAGS,BASE,fixture,once,sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fixture',type=Path,default=HERE/'fixture-v39.f')
    p.add_argument('--source',type=Path,default=HERE.parents[1]/'source');p.add_argument('--fc',default='gfortran');p.add_argument('--result',type=Path);a=p.parse_args()
    source=a.source.resolve();base=sp.check_output(['git','-C',str(source.parent),'show',BASE+':source/archive.f'])
    assert sha(base)=='e38b83224cfb10b38ec1eca494b6fff291966371110e1330ba60f0f1a5d4041b'
    new=(source/'archive.f').read_bytes();assert sha(new)==sha((HERE/'archive.f').read_bytes())
    src=fixture(a.fixture.read_bytes());rows={};gold={}
    with tempfile.TemporaryDirectory(prefix='smpi-defaults-') as t:
        d=Path(t);(d/'comdeck').write_text(once((source/'comdeck').read_text(),'IM=372,JM=479','IM=21,JM=19'))
        shutil.copyfile(source/'mpideck',d/'mpideck')
        for nth in (1,3):
            flags=[v for v in FLAGS if v!='-fopenmp']+(['-fopenmp'] if nth>1 else [])
            sp.run([a.fc,'-cpp',*flags,'-J.','-c',str(source/'secom_mpi.f90'),'-o','rt.o'],cwd=d,check=True)
            for role,b in [('v46',base),('default',new)]:
                s=once(src,'      CALL QA_ARC_SET(ILAYOUT)','! No policy assignment anywhere; actual initializer is under test.')
                if role=='default':
                    s=once(s,'      CALL SMPI_INIT(IM,JM)','      CALL QA_DEFAULT(0,'+str(nth)+')\n      CALL SMPI_INIT(IM,JM)\n      CALL QA_DEFAULT(JM,'+str(nth)+')')
                    s+='''      SUBROUTINE QA_DEFAULT(JEND,NTARGET)
      INCLUDE 'mpideck'
      INTEGER JEND,NTARGET,NWORK,NBAD
      LOGICAL SMPI_ARC_REQUESTED
      IF(SMPI_OUTER.OR.SMPI_ARC_REQUESTED()) STOP 71
      IF(SMPI_SIZE.NE.1.OR.SMPI_RANK.NE.0.OR.SMPI_JLO.NE.1) STOP 72
      IF(JEND.NE.0.AND.SMPI_JHI.NE.JEND) STOP 73
      IF(JEND.EQ.0.AND.SMPI_JHI.NE.HUGE(0)) STOP 74
      NBAD=0
      NWORK=0
!$OMP PARALLEL REDUCTION(+:NBAD,NWORK)
      NWORK=NWORK+1
      IF(SMPI_ARC_REQUESTED().OR.SMPI_OUTER) NBAD=NBAD+1
!$OMP END PARALLEL
      IF(NBAD.NE.0.OR.NWORK.NE.NTARGET) STOP 75
      END
'''
                (d/'driver.f').write_text(s);(d/'archive.f').write_bytes(b)
                sp.run([a.fc,*flags,'-I.','archive.f','driver.f','rt.o','-o','exe'],cwd=d,check=True)
                for mode in [2,10,11]:
                    ident=f'{role}-threads{nth}-mode{mode}';run=d/ident;run.mkdir()
                    q=sp.run([str(d/'exe'),str(mode),'0'],cwd=run,env=dict(os.environ,OMP_NUM_THREADS=str(nth),OMP_DYNAMIC='FALSE'),stdout=sp.PIPE,stderr=sp.STDOUT,timeout=25)
                    assert q.returncode==0,(ident,q.stdout)
                    hashes={f.name:sha(f.read_bytes()) for f in run.iterdir() if f.is_file()}
                    assert len(hashes)==(2 if mode==10 else 1),hashes
                    if nth==1 and role=='v46':gold[mode]=hashes
                    assert hashes==gold[mode],ident
                    rows[ident]=hashes
                print(role,'threads',nth,'PASS',flush=True)
    out=dict(compiler=sp.check_output([a.fc,'--version']).decode().splitlines()[0],source_sha=sha(new),runtime_sha=sha((source/'secom_mpi.f90').read_bytes()),cases=rows)
    if a.result:
        with a.result.open('x') as f:json.dump(out,f,indent=2)
if __name__=='__main__':main()
