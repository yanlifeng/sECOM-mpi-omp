"""Make normally unobserved domain-energy and other skill outputs observable."""
from pathlib import Path
import subprocess,os,shutil,hashlib,json,sys
root=Path(__file__).resolve().parent
old=Path('/home/user_home/ylf/secom-parallel-check-20261118.jNqsOL')
src=root/sys.argv[1]
std=old/'code-delivery/source/build/serial/s--pom.x'
template=old/'code-delivery/case20m/input'
out=root/('qa-skill-'+sys.argv[1]); out.mkdir()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
records=[]
for name,binary,n in [('serial',std,1),('omp',src/'build/omp/s--pom.x',32)]:
    d=out/name;shutil.copytree(template,d)
    p=d/'run_data.txt'; s=p.read_bytes()
    for a,b in [(b'      4800     21600',b'       300     21600'),(b'1200.0000',b' 100.0000'),(b'    999999\n',b'       100\n')]:
        assert s.count(a)==1,(a,p)
        s=s.replace(a,b)
    p.write_bytes(s)
    inputs={p.name:sha(p) for p in d.iterdir() if p.is_file()}
    cpus=list(range(40,40+n))
    env=dict(os.environ,LD_PRELOAD=str(old/'h5_preinit.so'),LD_LIBRARY_PATH='/home/user_home/ylf/kt5-licom-secom-ylf/deps/lib:'+os.getenv('LD_LIBRARY_PATH',''),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',OMP_PLACES=','.join('{%d}'%c for c in cpus),OMP_NUM_THREADS=str(n))
    with (out/(name+'.log')).open('x') as f:
        subprocess.run(['taskset','-c',','.join(map(str,cpus)),str(binary)],cwd=d,stdout=f,stderr=subprocess.STDOUT,env=env,check=True,timeout=180)
    files={p.name:sha(p) for p in d.iterdir() if p.name not in inputs}
    assert (d/'gcmtsr').stat().st_size>0 and len(files)==7 and b'JOB SUCCESSFULLY COMPLETED' in (d/'gcmprt').read_bytes()
    records.append(dict(binary=sha(binary),inputs=inputs,outputs=files,bytes_skill=(d/'gcmtsr').stat().st_size))
    assert records[0]['inputs']==inputs and records[0]['outputs']==files
(out/'results.json').write_text(json.dumps(records,indent=2))
