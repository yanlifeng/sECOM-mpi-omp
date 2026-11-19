"""Launcher for staged positive and negative p3 *units*, never a solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess as sp
import sys

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dir',type=Path);p.add_argument('--cpus',required=True);p.add_argument('--worker',action='store_true')
    p.add_argument('--ranks',type=int);p.add_argument('--threads',type=int);p.add_argument('--exe',type=Path);a=p.parse_args()
    d=a.dir.resolve();cs=list(map(int,a.cpus.split(',')));assert len(cs)==len(set(cs)) and 4<=min(cs)<=max(cs)<=35
    if a.worker:
        rank=int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK']);assert int(os.environ['OMPI_COMM_WORLD_SIZE'])==a.ranks
        subset=cs[rank*a.threads:(rank+1)*a.threads];assert len(subset)==a.threads and set(subset)<=os.sched_getaffinity(0)
        os.sched_setaffinity(0,subset);assert set(os.sched_getaffinity(0))==set(subset)
        print('rank',rank,'CPU',subset,flush=True)
        en=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='CLOSE',OMP_PLACES=','.join('{%d}'%c for c in subset))
        exe=str(a.exe.resolve());os.execve(exe,[exe,str(a.ranks),str(a.threads)],en)
    assert len(cs)>=7 and set(cs)<=os.sched_getaffinity(0)
    top=Path('/sys/devices/system/cpu');phys={tuple((top/f'cpu{c}'/'topology'/k).read_text().strip() for k in ['physical_package_id','core_id']) for c in cs}
    assert len(phys)==len(cs);os.sched_setaffinity(0,cs)
    root=Path(__file__).resolve().parent;py=sys.executable;report={};fc=shutil.which('mpifort');mp=shutil.which('mpirun');assert fc and mp
    for typ in ['checkready','both']:
        for mut in (['none'] if typ=='checkready' else ['none','public-stale','drop-reset','drop-vaf']):
            sp.run([py,'-B',str(root/'test_epochs.py'),'build',str(d/typ),'--fc',fc,'--mutant',mut],check=True)
            for ranks,nt in ([(2,1),(3,2),(7,1)] if typ=='both' and mut=='none' else [(3,2)]):
                exe=d/typ/('qa-'+mut);desc=f'{typ}-{mut}-{ranks}x{nt}';assert ranks*nt<=len(cs)
                cmd=[mp,'--nooversubscribe','--bind-to','none','--host',f'localhost:{ranks}','--mca','plm','isolated','--mca','pml','ob1','--mca','btl','self,vader','-np',str(ranks),
                    py,'-B',str(Path(__file__).resolve()),str(d),'--worker','--cpus',a.cpus,'--ranks',str(ranks),'--threads',str(nt),'--exe',str(exe)]
                q=sp.run(cmd,cwd=d,stdout=sp.PIPE,stderr=sp.STDOUT,timeout=45);s=q.stdout.decode(errors='replace')
                with (d/(desc+'.log')).open('x') as f:f.write(s)
                expect={'none':'WAD_OWNED','public-stale':'BAD PUBLIC','drop-reset':'BAD OWNER','drop-vaf':'BAD HALO'}[mut]
                if (q.returncode!=0)!=(mut!='none') or expect not in s:raise RuntimeError((desc,q.returncode,s))
                report[desc]={'rc':q.returncode,'marker':expect,'stdout_sha256':hashlib.sha256(q.stdout).hexdigest()}
                print(desc,'PASS',q.returncode,expect,flush=True)
    result=dict(cpus=cs,compiler=sp.check_output([fc,'--version']).decode().splitlines()[0],mpi=sp.check_output([mp,'--version']).decode().splitlines()[0],
        manifests={typ:json.loads((d/typ/'inputs.json').read_text()) for typ in ['checkready','both']},tests=report)
    with (d/'report.json').open('x') as f:json.dump(result,f,indent=2)
