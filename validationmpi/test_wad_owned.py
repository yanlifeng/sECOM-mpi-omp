"""Stage an isolated, miniature WANDD kernel for MPI; never starts a solver.

  python3 validationmpi/test_wad_owned.py stage /new/directory
  python3 validationmpi/test_wad_owned.py build /existing/stage --mpifort mpifort
  # With the campaign's approved MPI rank binding, run the executable
  #   wad-qa <expected-MPI-size> <expected-OpenMP-team-size>
  # Test 2/3 ranks (even/odd seams), also 7 ranks (single-row domains).
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARENT='1b58f33639538c2d214f14a8f9a575ceffed5267'

def stage(d,src):
    d.mkdir(parents=False,exist_ok=False)
    for name in ('wandd_kky.f','mpideck','secom_mpi.f90'):
        shutil.copyfile(src/name,d/name)
    com=(src/'comdeck').read_text()
    assert com.count('IM=372,JM=479')==1
    (d/'comdeck').write_text(com.replace('IM=372,JM=479','IM=12,JM=7'))
    shutil.copyfile(HERE/'wad_owned_driver.f',d/'wad_owned_driver.f')
    old=subprocess.check_output(['git','-C',str(src.parent),'show',PARENT+':source/wandd_kky.f']).decode()
    assert old.count('WANDD_KKY')==1 and old.count('WANDD_TEAM')==2
    (d/'v46_wandd.f').write_text(old.replace('WANDD_KKY','V46_WANDD_KKY').replace('WANDD_TEAM','V46_WANDD_TEAM'))
    files={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(d.iterdir())}
    (d/'fixture.json').write_text(json.dumps(dict(reference=PARENT,files=files),indent=2)+'\n')

def build(d,fc,drop_north_reset=False):
    manifest=json.loads((d/'fixture.json').read_text())['files']
    assert all(hashlib.sha256((d/n).read_bytes()).hexdigest()==h for n,h in manifest.items())
    flags=[fc,'-O3','-fopenmp','-ffp-contract=off','-fno-fast-math','-ffpe-trap=invalid,zero,overflow','-fcheck=bounds']
    subprocess.run(flags+['-cpp','-DSECOM_USE_MPI','-I.','-J.','-c','secom_mpi.f90','-o','runtime.o'],cwd=d,check=True)
    wad='wandd_kky.f'
    if drop_north_reset:
        # Deliberately lose ET/ETF updates at the seam, but not face closures.
        f=(d/wad).read_text()
        assert f.count('WAD_SEND(I)=3')==1
        wad='wad_missing_reset.f'
        with (d/wad).open('x') as out:out.write(f.replace('WAD_SEND(I)=3','WAD_SEND(I)=1'))
    subprocess.run(flags+['-I.','-std=legacy','-ffixed-line-length-132','wad_owned_driver.f',
                          'v46_wandd.f',wad,'runtime.o','-o','wad-mutant' if drop_north_reset else 'wad-qa'],cwd=d,check=True)

def run_p3(d,ranks,threads,rank_exec,negative,first_cpu):
    # p3 campaign budget. Do not rely on mpirun to infer the binding.
    assert (ranks,threads) in {(2,1),(3,2),(7,1)}
    cpus=list(range(first_cpu,first_cpu+ranks*threads));assert min(cpus)>=4 and max(cpus)<=35
    top=Path('/sys/devices/system/cpu')
    assert len({tuple((top/('cpu%d'%c)/'topology'/k).read_text().strip()
              for k in ('physical_package_id','core_id')) for c in cpus})==len(cpus)
    if rank_exec:
        me=int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK'])
        assert int(os.environ['OMPI_COMM_WORLD_SIZE'])==ranks
        owned=cpus[me*threads:(me+1)*threads];os.sched_setaffinity(0,owned)
        assert set(os.sched_getaffinity(0))==set(owned)
        env=dict(os.environ,OMP_NUM_THREADS=str(threads),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',
                 OMP_PLACES=','.join('{%d}'%i for i in owned))
        print('p3 rank',me,'CPU',owned,flush=True)
        exe=str(d/('wad-mutant' if negative else 'wad-qa'))
        os.execve(exe,[exe,str(ranks),str(threads)],env)
    os.sched_setaffinity(0,cpus)
    cmd=['/usr/lib64/openmpi/bin/mpirun','--host','localhost:%d'%ranks,'-np',str(ranks),
         '--nooversubscribe','--bind-to','none','--mca','plm','isolated','--mca','pml','ob1',
         '--mca','btl','self,vader',sys.executable,str(Path(__file__).resolve()),'p3-rank',str(d),
         '--ranks',str(ranks),'--threads',str(threads),'--first-cpu',str(first_cpu)]+(['--negative'] if negative else [])
    print('launcher',cmd,flush=True)
    if not negative:subprocess.run(cmd,check=True,timeout=45)
    else:
        q=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=45)
        print(q.stdout,end='')
        assert q.returncode!=0 and 'BAD OWNER' in q.stdout and 'WAD_OWNED' in q.stdout
        print('PASS: discarded north ET/ETF reset was rejected; mpirun exit',q.returncode)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['stage','build','p3-run','p3-rank']);p.add_argument('directory',type=Path)
    p.add_argument('--source',type=Path,default=HERE.parent/'source');p.add_argument('--mpifort',default='mpifort')
    p.add_argument('--ranks',type=int);p.add_argument('--threads',type=int)
    p.add_argument('--negative',action='store_true')
    p.add_argument('--first-cpu',type=int,default=4)
    a=p.parse_args()
    if a.action=='stage':stage(a.directory.resolve(),a.source.resolve())
    elif a.action=='build':build(a.directory.resolve(),a.mpifort,a.negative)
    else:run_p3(a.directory.resolve(),a.ranks,a.threads,a.action=='p3-rank',a.negative,a.first_cpu)
