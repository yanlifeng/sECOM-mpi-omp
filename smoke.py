"""Standalone same-machine port check. Usage: python3 smoke.py /new/empty/path 4"""
import hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def main():
    kit=Path(__file__).resolve().parent
    target=Path(sys.argv[1]).resolve()
    nt=int(sys.argv[2]) if len(sys.argv)>2 else 4
    if nt<1:raise ValueError('threads must be positive')
    ref=None; output={}
    paths={m:kit/'source/build'/m/'s--pom.x' for m in ('serial','omp')}
    if not all(p.is_file() for p in paths.values()):raise RuntimeError('build both MODE=serial and MODE=omp first')
    target.mkdir(parents=True,exist_ok=False)
    for mode,exe in paths.items():
        d=target/mode
        shutil.copytree(kit/'case20m/input',d)
        p=d/'run_data.txt';b=p.read_bytes()
        for orig,repl in [(b'      4800     21600',b'       300     21600'),(b'1200.0000',b' 100.0000')]:
            if b.count(orig)!=1:raise RuntimeError('unexpected run_data template')
            b=b.replace(orig,repl)
        p.write_bytes(b)
        inputs={p.name:digest(p) for p in d.iterdir() if p.is_file()}
        if mode=='serial': refinputs=inputs
        if inputs!=refinputs:raise RuntimeError('input mismatch')
        env=dict(os.environ,OMP_NUM_THREADS=str(nt if mode=='omp' else 1),OMP_DYNAMIC='FALSE')
        start=time.monotonic()
        with (target/(mode+'.log')).open('x') as log:
            subprocess.run([str(exe)],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        elapsed=time.monotonic()-start
        hashes={p.name:digest(p) for p in d.iterdir() if p.is_file() and p.name not in inputs}
        if set(hashes)!={'gcmprt','gcmplt.cdf','heat.dat','startup','fort.1002','fort.1006'}:
            raise RuntimeError('unexpected output set, inspect logs')
        if ref is None:ref=hashes
        output[mode]=dict(elapsed=elapsed,outputs=hashes,executable_sha=digest(exe),matches_serial=hashes==ref)
    print(json.dumps(output,indent=2))
    return 0 if output['omp']['matches_serial'] else 1
if __name__=='__main__':sys.exit(main())
