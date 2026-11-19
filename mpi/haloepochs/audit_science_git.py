"""Compare an actual Git-object source tree with separately hashed working files.

Read-only. Before final commit use --allow-uncommitted to show the delta from
the current Git object. Afterwards omit this flag: any mismatch is a failure.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess as sp

ROOT=Path(__file__).resolve().parents[2];EXPECTED='46f5da3de7d5904ed20e902e1c1ba501f2ca40fa090e11c07a8d7754a6b9d1a4'
def sha(x):return hashlib.sha256(x).hexdigest()
def canon(x):return sha(json.dumps(x,sort_keys=True,separators=(',',':')).encode())
def git(*args,stdin=None):return sp.check_output(['git','-C',str(ROOT),*args],input=stdin)

def git_sources(ref):
    ref=git('rev-parse','--verify',ref+'^{commit}').decode().strip()
    records={}
    for row in git('ls-tree','-rz',ref,'--','source').split(b'\0'):
        if not row:continue
        info,path=row.split(b'\t');mode,kind,oid=info.split();assert kind==b'blob' and mode in {b'100644',b'100755'}
        assert path.startswith(b'source/');name=path[7:].decode();assert '/' not in name and name not in records
        records[name]=git('cat-file','blob',oid.decode())
    return ref,records

def proc(data,name):
    t=data.decode('latin1');m=re.search(r'(?im)^\s*(?:RECURSIVE\s+)?SUBROUTINE\s+'+name+r'\s*\(',t)
    assert m,name;t=t[m.start():];end=re.search(r'(?im)^\s+END[ \t]*$',t);assert end,name
    return '\n'.join(line.rstrip() for line in t[:end.end()].strip().splitlines()).encode('latin1')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--commit',default='HEAD');p.add_argument('--allow-uncommitted',action='store_true');p.add_argument('--result',type=Path);a=p.parse_args()
    ref,store=git_sources(a.commit);pr,base=git_sources('c371fa3');src=ROOT/'source'
    work={}
    for f in src.iterdir():
        assert f.is_file() and not f.is_symlink(),f;work[f.name]=f.read_bytes()
    assert set(base)==set(store)==set(work) and len(work)==107
    wm={k:sha(v) for k,v in work.items()};gm={k:sha(v) for k,v in store.items()};old={k:sha(v) for k,v in base.items()}
    assert canon(wm)==EXPECTED,(canon(wm),EXPECTED)
    public={}
    for f,n in [('wandd_kky.f','WANDD_KKY'),('wandd_kky.f','WANDD_TEAM'),('advave.f','ADVAVE'),('extrnl.f','EXTRNL'),('extrnl.f','SECOM_ADV34_CORNERS')]:
        b,c=proc(base[f],n),proc(work[f],n);assert b==c,n;public[n]=sha(c)
    diff=sorted(k for k in wm if gm[k]!=wm[k]);delta={k:{'r004b':old[k],'frozen':wm[k]} for k in sorted(wm) if old[k]!=wm[k]}
    report=dict(worktree_canonical=canon(wm),git_commit=ref,git_canonical=canon(gm),git_worktree_mismatches=diff,matched_commit=not diff,
                r004b_commit=pr,r004b_canonical=canon(old),delta_from_r004b=delta,normalized_public_procedure_hashes=public,source_count=len(work))
    print(json.dumps({k:report[k] for k in ['source_count','git_commit','worktree_canonical','git_canonical','matched_commit','git_worktree_mismatches']},indent=2))
    if a.result:
        with a.result.open('x') as f:json.dump(report,f,indent=2)
    if diff and not a.allow_uncommitted:raise SystemExit('Git object does not contain the frozen 107-file source mapping')
if __name__=='__main__':main()
