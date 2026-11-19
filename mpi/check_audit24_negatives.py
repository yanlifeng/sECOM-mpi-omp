"""Targeted negative controls using real r014 records; substitutions are in memory only."""
import argparse
import json
from pathlib import Path
from unittest.mock import patch
import audit_acceptance24 as a
import audit_evidence as ref
import audit_trace_sources as ts
import check_regress24 as regress


def check(bundle,side):
    result={};prefix='r014-p32-t1-'; tags=[prefix+'full%02d'%j for j in range(1,6)]
    science='8398560238451a7114df8ccd47917751a8b56a6b2aea515ed5a2f4eb3f40cf65'
    realread=ref.read; proof=realread(bundle/'reference/day19-serial-proof.json')
    def must_reject(name,fn,reason):
        try:fn()
        except ValueError as ex:
            if reason not in str(ex):raise
            result[name]=str(ex)
        else:raise ValueError('negative control was accepted: '+name)
    def replace_receipt(path,mutate):
        def read(p):
            v=realread(p)
            if p==path:mutate(v)
            return v
        return patch.object(ref,'read',side_effect=read)
    def formal0():return a.run(bundle,'r014final',tags[0],'science',science,proof)
    r0=bundle/'runs'/tags[0]/'result.json'
    def badmask(v):
        o=v['rank_binding_samples']['17'];o['masks'][str(o['pid'])]=[36]
    with replace_receipt(r0,badmask):must_reject('formal_rank_on_cpu36',formal0,'observed OS affinity rank 17')
    def badsha(v):v['checks']['0']['outputs']['startup']='0'*64
    with replace_receipt(r0,badsha):must_reject('formal_bad_restart_sha',formal0,'input/output receipt rank 0')
    # Replace actual clock pairs for the second AB/BA pair; every binding, SHA,
    # elapsed value and frozen plan remains correct. Reject specifically order.
    times={};it={}
    for t in tags:
        r=realread(bundle/'runs'/t/'result.json'); times[t]=r['elapsed'];it[t]=(r['started'],r['started']+r['elapsed'])
    c2=realread(bundle/'runs/r014-control02/result.json');old=it[tags[1]][0]
    it[tags[1]]=(c2['started'],c2['started']+times[tags[1]])
    orig=a.run
    def swap(root,revision,tag,*args):
        e,s,k=orig(root,revision,tag,*args)
        return (e,old if tag=='r014-control02' else s,k)
    with patch.object(a,'run',side_effect=swap):
        must_reject('swap_second_registered_pair',lambda:a.paired_r014(bundle,tags,times,it,proof,a.R014_PLAN),
                    'actual crossover order differs')
    one=bundle/'runs'/(prefix+'serial300')
    def wrongio(v):v['root_only_large_outputs']=True
    def both(p):
        v=realread(p)
        if p.parent==one and p.name in ('result.json','started.json'):wrongio(v)
        return v
    with patch.object(ref,'read',side_effect=both):
        must_reject('misstated_one_rank_regression_mode',lambda:regress.audit(bundle,'r014final',prefix,proof,science),
                    'shape/CPU bounds')
    # Bypass the independent transfer-hash layer in this simulation only;
    # mutate read_bytes, whereas ref.sha hashes the unmodified file stream.
    realbytes=Path.read_bytes
    def readbytes(p):
        raw=realbytes(p)
        if p==side/'early/capture.f':
            raw=raw.replace(b'CALL SMPI_ALL_RAW(BUFFER,IM,JM)',b'CALL SMPI_ALL_RAW(CAPFIELD,IM,JM)')
        return raw
    with patch.object(Path,'read_bytes',readbytes):
        must_reject('early_capture_changed_to_inplace',lambda:ts.audit(bundle,side),'regenerated capture.f differs')
    return dict(passed=True,modified_evidence=False,expected_rejections=result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('evidence',type=Path);p.add_argument('sidecar',type=Path)
    p.add_argument('--output',type=Path);u=p.parse_args()
    if u.output and u.evidence.resolve() in u.output.resolve().parents:p.error('report outside immutable evidence only')
    s=json.dumps(check(u.evidence.resolve(),u.sidecar.resolve()),indent=2)+'\n'
    if u.output:
        with u.output.open('x') as f:f.write(s)
    print(s,end='')
