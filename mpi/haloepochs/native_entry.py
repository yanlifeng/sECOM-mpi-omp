"""Run staged checked MPI kernels on explicitly assigned p3 cores. No solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess as sp
import sys
from test_entry import FIELDS,FLAGS

def shax(b):return hashlib.sha256(b).hexdigest()
def oracle_reference(allrows,n,rank):
    nx,ny=372,19;na,extra=divmod(ny,n);lo=rank*na+min(rank,extra);length=na+(rank<extra)
    globalplane=4*nx*ny;block=24+len(FIELDS)*globalplane;assert len(allrows)==8*block
    out=[]
    for k in range(8):
        start=block*k;out.append(allrows[start:start+24]);start+=24
        for f in FIELDS:
            base=start+4*nx*lo;out.append(allrows[base:base+4*nx*length]);start+=globalplane
    return b''.join(out)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dir',type=Path);p.add_argument('--cpus',required=True);p.add_argument('--wad5',action='store_true');p.add_argument('--corner',action='store_true');p.add_argument('--raw',choices=['pair','entry','corner']);p.add_argument('--worker',action='store_true');p.add_argument('--ranks',type=int);p.add_argument('--threads',type=int);p.add_argument('--science-ranks',type=int);p.add_argument('--mode',type=int);p.add_argument('--exe',type=Path);a=p.parse_args()
    d=a.dir.resolve();cs=list(map(int,a.cpus.split(',')));assert len(cs)==len(set(cs)) and 4<=min(cs)<=max(cs)<=35
    if a.worker:
        rank=int(os.environ['OMPI_COMM_WORLD_LOCAL_RANK']);assert int(os.environ['OMPI_COMM_WORLD_SIZE'])==a.ranks
        sel=cs[rank*a.threads:(rank+1)*a.threads];assert len(sel)==a.threads and set(sel)<=os.sched_getaffinity(0)
        os.sched_setaffinity(0,sel);assert set(os.sched_getaffinity(0))==set(sel)
        en=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='CLOSE',OMP_PLACES=','.join('{%d}'%c for c in sel))
        print('rank',rank,'CPU',sel,flush=True);os.chdir('rank'+str(rank));exe=str(a.exe.resolve());args=[str(a.mode),str(a.science_ranks)]
        if a.raw=='entry':args=[str(a.science_ranks),str(a.ranks),str(a.threads)]
        if a.raw=='pair':args=[str(a.ranks),str(a.threads),str(a.science_ranks),'reverse' if a.mode==2 else 'none']
        if a.raw=='corner':args=[str(a.ranks)]
        os.execve(exe,[exe]+args,en)
    assert len(cs)>=4 and set(cs)<=os.sched_getaffinity(0)
    top=Path('/sys/devices/system/cpu');assert len({tuple((top/f'cpu{c}'/'topology'/k).read_text().strip() for k in ('core_id','physical_package_id')) for c in cs})==len(cs);os.sched_setaffinity(0,cs)
    fc=shutil.which('mpifort');mp=shutil.which('mpirun');assert fc and mp;py=sys.executable;here=Path(__file__).resolve().parent
    roles=['old','cert','bad-skip-ua'] if a.corner else ['old','five','bad-el-south','bad-vaf-north'] if a.wad5 else ['old','entry','wadinput','bad-first-el','bad-next-vab','bad-wad-mask'];report={}
    for role in roles:sp.run([py,'-B',str(here/('test_corner.py' if a.corner else 'test_entry.py')),'build',str(d/role),'--fc',fc],check=True)
    def run(name,role,n,nt,mode,phys,raw=None):
        z=d/name;z.mkdir()
        for r in range(n):(z/f'rank{r}').mkdir()
        cmd=[mp,'--nooversubscribe','--bind-to','none','--host',f'localhost:{n}','--mca','plm','isolated','--mca','pml','ob1','--mca','btl','self,vader','-np',str(n),py,'-B',str(Path(__file__).resolve()),str(d),'--worker','--cpus',a.cpus,'--ranks',str(n),'--threads',str(nt),'--mode',str(mode),'--science-ranks',str(phys),'--exe',str(d/role/'unit')]
        if raw:cmd=cmd[:-1]+[str(d/role/('qa-'+raw)),'--raw',raw]
        q=sp.run(cmd,cwd=z,stdout=sp.PIPE,stderr=sp.STDOUT,timeout=45);(z/'placement.log').write_bytes(q.stdout)
        if raw:
            if raw=='pair' and mode==2:
                assert q.returncode==41 and b'DETECTED_WAD_PAIR_WRONG_DIRECTION' in q.stdout,(name,q.returncode,q.stdout)
            else:assert q.returncode==0 and b'PASS' in q.stdout,(name,q.returncode,q.stdout)
            report[name]={'rc':q.returncode,'log':q.stdout.decode()};print('PASS',name,q.returncode,flush=True);return
        if q.returncode:raise RuntimeError((name,q.returncode,q.stdout))
        return [(z/f'rank{r}'/'owner.raw').read_bytes() for r in range(n)]
    for n,nt in [(3,1),(2,2)]:
        for mode in (([1,3,4,5,6] if n==3 else [1,4]) if a.corner else ([1,2,3] if n==3 else [1])):
            whole=run(f'oracle-{n}-{mode}','old',1,1,mode,n)[0];gold=[oracle_reference(whole,n,r) for r in range(n)]
            for role in (roles if mode==1 and n==3 else [r for r in roles if not r.startswith('bad-')]):
                name=f'{role}-{n}x{nt}-mode{mode}';out=run(name,role,n,nt,mode,n);mismatches=[]
                for r,(x,y) in enumerate(zip(out,gold)):
                    assert len(x)==len(y);i=next((k for k,(u,v) in enumerate(zip(x,y)) if u!=v),None) if x!=y else None
                    if i is not None:
                        width=(19//n+(r<19%n))*372*4;step,inside=divmod(i,24+len(FIELDS)*width)
                        field='scalar' if inside<24 else FIELDS[(inside-24)//width];mismatches.append((r,step+1,field,i))
                if bool(mismatches)!=role.startswith('bad-'):raise RuntimeError((name,'unexpected mismatch/match',mismatches))
                report[name]={'files':[shax(x) for x in out],'oracle_files':[shax(x) for x in gold],'detected':mismatches};print('PASS',name,mismatches,flush=True)
    if a.wad5:
        for raw,filename in [('pair','check_wad_pair.f90'),('entry','check_entry_epoch.f90')]:
            z=d/'five';sp.run([fc,*FLAGS,'-I.','-J.',str(here/filename),'rt.o','-o','qa-'+raw],cwd=z,check=True)
            for n,nt,jm,mode in [(4,1,4,1),(3,2,8,1),(1,2,7,1)]+([(2,1,5,2)] if raw=='pair' else []):
                run(f'raw-{raw}-{n}x{nt}-jm{jm}-mode{mode}','five',n,nt,mode,jm,raw)
    if a.corner:
        z=d/'cert';sp.run([fc,*FLAGS,'-I.','-J.',str(here/'check_corner_epoch.f'),'mpi_boundary.f','rt.o','-o','qa-corner'],cwd=z,check=True)
        for n in [1,2,3]:run(f'corner-meta-{n}','cert',n,1,1,19,'corner')
    result=dict(cpus=cs,compiler=sp.check_output([fc,'--version']).decode().splitlines()[0],mpi=sp.check_output([mp,'--version']).decode().splitlines()[0],flags=FLAGS,
        inputs={r:json.loads((d/r/'inputs.json').read_text()) for r in roles},cases=report)
    with (d/'report.json').open('x') as f:json.dump(result,f,indent=2)
