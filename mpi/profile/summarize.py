"""Aggregate rank-local diagnostic files; no cross-rank clock subtraction."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import statistics


def stats(values):
    return dict(min=min(values),median=statistics.median(values),max=max(values),sum_over_ranks=sum(values))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run',type=Path)
    ap.add_argument('--ranks',type=int,required=True)
    ap.add_argument('--threads',type=int,required=True)
    ap.add_argument('--steps',type=int,required=True)
    ap.add_argument('--split',type=int,required=True)
    ap.add_argument('--expect-owned-epochs',action='store_true',
                    help='Require each rank to have measured the unpacked-entry, healthy-key and WAD5 macros')
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    files=sorted(a.run.glob('rank*/smpi-profile-*.log'))
    if len(files)!=a.ranks:raise ValueError('missing/extra profile files')
    samples={}
    for f in files:
        lines=f.read_text().splitlines()
        fields=lines[0].split()
        rank,n,threads,it,done=map(int,fields[1:])
        if (fields[0]!='SECOM_DIAG_V1' or rank in samples or
            (n,threads,it,done)!=(a.ranks,a.threads,a.steps,1) or f.parent.name!=f'rank{rank:03d}' or
            f.name!=f'smpi-profile-{rank:06d}.log'):
            raise ValueError('invalid or partial profile: '+str(f))
        sample={}
        for l in lines[1:]:
            kind,life,site,op,n,s=l.split()
            key=(kind,life,site,op)
            if key in sample or not math.isfinite(float(s)) or float(s)<0 or int(n)<1:
                raise ValueError('invalid metric: '+l)
            sample[key]=(int(n),float(s))
        if len([k for k in sample if k[0]=='life'])!=5:raise ValueError('missing lifecycle phases')
        for site in ['bcond7','cbc','bottau','baro','extrnl','wetdry','check']:
            expect=a.steps*(a.split if site in ('extrnl','wetdry','check') else 1)
            if sample.get(('site','iteration',site,'inclusive'),(0,))[0]!=expect:
                raise ValueError(f'{f}: unmatched {site} count; expected {expect}')
        if a.expect_owned_epochs:
            subs=a.steps*a.split
            for site,name,want in [('extrnl','entry-sendrecv',a.steps*(9+5*(a.split-1))),
                                   ('check','checkkeys-max',subs),
                                   ('wetdry','wad-pair-sendrecv',4*subs),
                                   ('wetdry','wad-sendrecv',subs),
                                   ('wetdry','oneway-sendrecv',subs)]:
                actual=sample.get(('mpi','iteration',site,name),(0,))[0]
                if actual!=want:raise ValueError(f'{f}: {site}/{name}: expected {want}, got {actual}')
        mpi_sum=sum(s for (k,life,site,op),(n,s) in sample.items() if k=='mpi' and life=='iteration')
        loop=sample['life','iteration','all','wall'][1]
        if loop<=0 or mpi_sum>loop*(1+1e-6):raise ValueError('invalid serialized MPI/loop denominator')
        sample['derived','iteration','all','mpi-seconds']=(1,mpi_sum)
        sample['derived','iteration','all','mpi-fraction-of-local-loop']=(1,mpi_sum/loop)
        samples[rank]=sample
    if sorted(samples)!=list(range(a.ranks)):raise ValueError('rank ids not contiguous')
    keys=sorted({k for s in samples.values() for k in s})
    rows=[]
    for key in keys:
        vals=[s.get(key,(0,0.)) for s in samples.values()]
        unit='fraction' if key[-1].startswith('mpi-fraction') else 'seconds'
        rows.append(dict(zip(['kind','life','site','event'],key),count=stats([n for n,_ in vals]),
                         unit=unit,value=stats([s for _,s in vals])))
    report=dict(kind='instrumented_diagnostic_only',ranks=a.ranks,threads=a.threads,steps=a.steps,split=a.split,
                checked_owned_epoch_counts=a.expect_owned_epochs,
                note='Rank maxima/medians are not full-job throughput. Phase/subphase/MPI timings are nested. '
                     'Sum over ranks is rank-seconds, not seconds elapsed; a sum of fractions is not a fraction of wall.',
                metrics=rows)
    with a.out.open('x') as f: json.dump(report,f,indent=2,allow_nan=False)
    print('profiled ranks:',a.ranks,'valid outer steps per rank:',a.steps,'metrics:',len(rows))


if __name__=='__main__':
    main()
