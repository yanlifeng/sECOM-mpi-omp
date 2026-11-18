"""Independent re-read of the actual (completed) full-run directories and topology."""
from pathlib import Path
import hashlib,json,subprocess
root=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
def top():
    s=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET,NODE'],universal_newlines=True)
    return {int(a[0]):[int(x) for x in a[1:]] for a in (l.split(',') for l in s.splitlines() if l and not l.startswith('#'))}
topology=top()
tags=['fullconfirm4-01-s1','full6-00-f32','full6-02-f32','full6-03-f32','full6same-00-f1','full6same-01-f32']
base=None;runs=[]
for tag in tags:
    rec=json.loads((root/(tag+'.json')).read_text())
    d=root/tag
    names={'startup','gcmplt.cdf','gcmprt','fort.1002','fort.1006','heat.dat'}
    actual_files={f.name:f for f in d.iterdir() if f.is_file()}
    assert set(actual_files)-names==set(rec['input_sha'])
    actual={n:sha(f) for n,f in actual_files.items()}
    assert all(actual[n]==v for n,v in dict(rec['input_sha'],**rec['outputs']).items())
    if base is None:base=actual
    assert base==actual
    assert b'JOB SUCCESSFULLY COMPLETED' in (d/'gcmprt').read_bytes()
    assert (d/'startup').stat().st_size>0 and (d/'gcmplt.cdf').stat().st_size>0
    places=rec['single_cpu_bindings_seen']
    assert len({tuple(topology[c][:2]) for c in places})==rec['threads']==rec['max_threads']
    assert len({topology[c][1] for c in places})==1
    runs.append(dict(tag=tag,steps=rec['steps'],wall_s=rec['elapsed'],physical_cores=len(places),socket=topology[places[0]][1],binary_sha=rec['binary_sha'],files_checked=len(actual)))
files=list((root/'src6').glob('*.f'))+[root/'src6'/n for n in ('comdeck','netcdf.inc','makefile','Makefile.portable') if (root/'src6'/n).exists()]
testbin=sha(root/'src6/build/omp/s--pom.x')
assert all(r['binary_sha']==testbin for r in runs if 'f' in r['tag'].split('-')[-1])
out=dict(runs=runs,output_sha={n:base[n] for n in sorted(names)},input_sha={n:v for n,v in base.items() if n not in names},source_sha={f.name:sha(f) for f in files},cpu_map={c:topology[c] for c in range(40,72)},preload_sha=sha(Path('/home/user_home/ylf/secom-parallel-check-20261118.jNqsOL/h5_preinit.so')))
print(json.dumps(out,indent=2,sort_keys=True))
