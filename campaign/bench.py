"""p3 regression/performance; read-only v1 and closed-form named run lists."""
import hashlib,json,os,resource,shutil,subprocess,sys,time
from pathlib import Path

root=Path(__file__).resolve().parent
old=Path('/home/user_home/ylf/secom-parallel-check-20261118.jNqsOL')
template=old/'code-delivery/case20m/input'
bindir=old/'code-delivery/source/build'
bins={'s':bindir/'serial/s--pom.x','a':bindir/'omp/s--pom.x','b':root/'source/build/omp/s--pom.x'}
bins.update({'h':root/'mix-history','t':root/'mix-transport'})
prev=Path('/home/user_home/ylf/secom-omp2-20261118.ZtxsGz')
bins.update({'b':prev/'source/build/omp/s--pom.x',
             'c':root/'source/build/omp/s--pom.x','p':root/'profile-v2','q':root/'profile-v3'})
bins.update({'d':root/'src4/build/omp/s--pom.x','r':root/'profile-v4'})
bins['e']=root/'src5/build/omp/s--pom.x'
bins['f']=root/'src6/build/omp/s--pom.x'
outnames={'fort.1002','fort.1006','startup','gcmplt.cdf','gcmprt','heat.dat'}
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''): h.update(block)
    return h.hexdigest()
def snapshot(d):
    return {p.name:sha(p) for p in sorted(d.iterdir()) if p.is_file()}
def proc_times():
    return {a[0]:[int(x) for x in a[1:9]] for a in
            (s.split() for s in Path('/proc/stat').read_text().splitlines()) if a[0].startswith('cpu')}
phase=sys.argv[1]
# Validation at one thread; counterbalanced old/new at both widths.
phases={'short':['b1','a4','b4','b8','a8','b4','a4','a8','b8'],
        'extra':['a8','t8','h8','b8','b16','a16','a16','b16','b8','h8','t8','a8'],
        'long':['b8','a8','s1']}
phases.update({'base32':['b32'],'profile':['p32'],'pilot':['c1','c8','c32'],
               'diag3':['q32'],'full32':['c32','b32','s1','c32']})
if len(sys.argv)>2: phases[phase]=sys.argv[2:]
steps=4800 if phase.startswith(('long','full')) else 300
shortref=json.loads((old/'h5warm/serial.json').read_text())
fullinputs=json.loads((prev/'frozen.json').read_text())['input_sha256']
longref=json.loads((prev/'long-02-s1.json').read_text())
for k,entry in enumerate(phases[phase]):
    variant,nt=entry[0],int(entry[1:])
    cpus=list(range(40,40+nt))
    assert set(cpus)<=os.sched_getaffinity(0)
    tag='%s-%02d-%s'%(phase,k,entry)
    d=root/tag
    shutil.copytree(template,d)
    p=d/'run_data.txt'
    if steps==300:
        buff=p.read_bytes()
        for a,b in [(b'      4800     21600',b'       300     21600'),(b'1200.0000',b' 100.0000')]:
            assert buff.count(a)==1
            buff=buff.replace(a,b)
        p.write_bytes(buff)
    inputs=snapshot(d)
    assert inputs==(shortref['input_hashes'] if steps==300 else fullinputs)
    if variant in ('s','a'):
        ref_bin=json.loads((old/('h5warm/%s.json'%('serial' if variant=='s' else 'omp4'))).read_text())['bin_sha']
        assert sha(bins[variant])==ref_bin
    env=dict(os.environ,LD_PRELOAD=str(old/'h5_preinit.so'),
             LD_LIBRARY_PATH='/home/user_home/ylf/kt5-licom-secom-ylf/deps/lib:'+os.getenv('LD_LIBRARY_PATH',''),
             OMP_NUM_THREADS=str(nt),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='TRUE',
             OMP_PLACES=','.join('{%d}'%i for i in cpus),OPENBLAS_NUM_THREADS='1')
    seen=set(); max_threads=0; t0=proc_times(); load0=os.getloadavg()
    start=time.monotonic()
    with (root/(tag+'.log')).open('x') as log:
        pr=subprocess.Popen([str(bins[variant])],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,
                            preexec_fn=lambda: os.sched_setaffinity(0,cpus))
        while pr.poll() is None:
            threads=list(Path('/proc/%d/task'%pr.pid).glob('[0-9]*'))
            max_threads=max(max_threads,len(threads))
            if len(seen)<nt:
                for t in threads:
                    try: mask=os.sched_getaffinity(int(t.name))
                    except ProcessLookupError: continue
                    if len(mask)==1: seen.update(mask)
            if time.monotonic()-start>1200:
                pr.kill();pr.wait();raise RuntimeError('timeout '+tag)
            time.sleep(0.08)
    elapsed=time.monotonic()-start
    assert pr.returncode==0,(tag,pr.returncode)
    assert seen==set(cpus) and max_threads==nt,(tag,seen,max_threads)
    t1=proc_times()
    outputs={p.name:sha(p) for p in d.iterdir() if p.name not in inputs}
    assert set(outputs)==outnames
    assert all(sha(d/name)==sig for name,sig in inputs.items())
    result=dict(tag=tag,variant=variant,threads=nt,elapsed=elapsed,steps=steps,
                input_sha=inputs,binary_sha=sha(bins[variant]),outputs=outputs,
                single_cpu_bindings_seen=sorted(seen),max_threads=max_threads,
                load0=load0,load1=os.getloadavg(),cpu_tick_delta={i:[y-x for x,y in zip(t0['cpu%d'%i],t1['cpu%d'%i])] for i in cpus})
    if steps==300:
        result['matches_p3_serial']=outputs=={n:v['sha256'] for n,v in shortref['outputs'].items()}
    else:
        result['matches_p3_serial']=outputs==longref['outputs']
    with (root/(tag+'.json')).open('x') as f:json.dump(result,f,indent=2)
    print(tag,round(elapsed,3),result['matches_p3_serial'],flush=True)
    assert result['matches_p3_serial'],tag
