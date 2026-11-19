"""Audit every full-plane MPI actual in mpi_outer.f against compiler-declared types.

Missing explicit shape in legacy fixed-form Fortran can silently mean scalar REAL.
Also run a negative control reintroducing the erroneous local PRECTEMP actual.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess as sp

def symbols(fc,source,src):
    q=sp.run([fc,'-fsyntax-only','-std=legacy','-ffixed-line-length-132','-fdump-fortran-original',f'-I{source}','-x','f77','-'],input=src.encode(),stdout=sp.PIPE,stderr=sp.PIPE)
    if q.returncode:raise RuntimeError(q.stderr.decode())
    dump=q.stdout.decode()
    syms=dict(re.findall(r"^  symtree: '(\w+)'[^\n]*\n((?:    [^\n]*\n)+)",dump,re.M))
    return syms

def statements(src):
    current='';first=None
    for n,line in enumerate(src.splitlines(),1):
        if not line.strip() or line[:1].lower() in ('c','*','!'):continue
        body=line[6:].split('!')[0].lower().strip()
        if not body:continue
        if line[5:6] not in ('',' ','0'):
            current+=' '+body
        else:
            if current:yield first,current
            first=n;current=body
    if current:yield first,current

def audit(fc,source,src):
    syms=symbols(fc,source,src)
    def param(name):return int(re.search(r'value: (\d+)',syms[name]).group(1))
    im,jm,kb=(param(i) for i in ['im','jm','kb'])
    accepted=[];reject=[];loops=[]
    for lineno,s in statements(src):
        compact=re.sub(r'\s+','',s)
        if re.match(r'do\s+\w+\s*=',s):
            loops.append(compact)
        if compact=='enddo' and loops:loops.pop()
        if not re.search(r'\bcall\s+smpi_all_raw\b',s):continue
        # Unknown syntax is a failure, not an unobserved call. Be deliberately
        # conservative with future guards, arithmetic designators, or kinds.
        m=re.fullmatch(r'callsmpi_all_raw\((\w+)(\(1,1,k\))?,im,jm\)',compact)
        if not m:reject.append((lineno,s,'unrecognized call/shape'));continue
        name=m[1];rank=3 if m[2] else 2;block=syms[name]
        kind=re.search(r'type spec\s*:\s*\(([^)]+)\)',block).group(1)
        a=re.search(r'Array spec:\((\d+)\s+\[0\]\s+AS_EXPLICIT\s+([\d\s]+)\)',block)
        bounds=list(map(int,a[2].split())) if a else []
        expected=[1,im,1,jm]+([1,kb] if rank==3 else [])
        if kind!='REAL 4' or not a or int(a[1])!=rank or bounds!=expected:
            reject.append((lineno,name,kind,bounds));continue
        if rank==3 and 'dok=1,kb' not in loops:
            reject.append((lineno,name,'K not inside DO K=1,KB'));continue
        if not ('IN-COMMON' in block or 'DUMMY' in block):
            reject.append((lineno,name,'not a common/dummy'));continue
        accepted.append(dict(line=lineno,name=name,dimensions=bounds[1::2],element_type=kind,placement='dummy' if 'DUMMY' in block else 'common'))
    return accepted,reject

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[2]/'source')
    p.add_argument('--fc',default='gfortran');p.add_argument('--result',type=Path);a=p.parse_args()
    s=(a.source/'mpi_outer.f').read_text();valid,bad=audit(a.fc,a.source,s)
    assert not bad and valid,(bad,valid)
    needle='      CALL SMPI_ALL_RAW(CBC,IM,JM)';assert s.count(needle)==1
    _,negative=audit(a.fc,a.source,s.replace(needle,needle+'\n      CALL SMPI_ALL_RAW(PRECTEMP,IM,JM)'))
    assert len(negative)==1 and negative[0][1]=='prectemp',negative
    _,badkind=audit(a.fc,a.source,s.replace(needle,'      CALL SMPI_ALL_RAW(ARCTHIK,IM,JM)'))
    assert len(badkind)==1 and badkind[0][1:3]==('arcthik','REAL 8'),badkind
    _,badk=audit(a.fc,a.source,s.replace('      DO K=1,KB','      DO K=0,KB'))
    assert len(badk)==7 and all(x[2]=='K not inside DO K=1,KB' for x in badk),badk
    out=dict(compiler=sp.check_output([a.fc,'--version']).decode().splitlines()[0],source_sha=hashlib.sha256(s.encode()).hexdigest(),
             accepted=valid,rejected_mutation=negative,rejected_kind=badkind,rejected_k_bounds=badk)
    print(f'PASS {len(valid)} calls; 2D={sum(len(x["dimensions"])==2 for x in valid)}; scalar PRECTEMP, REAL8 and invalid K bounds rejected')
    if a.result:
        with a.result.open('x') as f:json.dump(out,f,indent=2)
if __name__=='__main__':main()
