"""Offline validation of one predeclared five-sample 24x candidate. Does not run MPI.

Supply SHA-256 of the science, early and late binaries and a frozen source-tree
SHA from *outside* the result bundle. Source digest: canon({basename: file-SHA})
of all regular files immediately inside the frozen science source directory.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import statistics
import struct

import audit_evidence as ref


# Pinned outside the bundle before its transfer. Never select the signed-off
# sequence or control ELF based on strings found in the candidate evidence.
R014_PLAN='172063efa3a978d4d99d6901dae353236351761e5d923ce2143d8749149c3b5a'


def scoped(prefix):
    """Use original inventory/trace validation with this invocation's names.

    The original implementations use module-level run-name constants; isolate
    them in a private module instance, leaving the old verifier unchanged.
    No path or reference hashes are read from the candidate to select the oracle.
    """
    spec = importlib.util.spec_from_file_location('_accept24_original', Path(ref.__file__))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.PREFIX = prefix
    m.REP = prefix + 'full03'
    return m


def digestarg(s):
    if not re.fullmatch(r'[a-f0-9]{64}', s):
        raise argparse.ArgumentTypeError('explicit lowercase SHA-256 required')
    return s


def namearg(s):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,100}', s) or '..' in s:
        raise argparse.ArgumentTypeError('expected a safe single path component')
    return s


def mode_receipt(path, rank):
    stdout=path.read_text().splitlines()
    lines = [s for s in stdout if 'MPI DOMAIN' in s]
    ref.need(len(lines) == 1, str(path) + ': expected one MPI DOMAIN receipt')
    m = re.fullmatch(r'\s*MPI DOMAIN\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+LOCALBC\s+T\s+OWNED_OUTER\s+T\s*', lines[0])
    q,r=divmod(ref.HEIGHT,ref.RANKS)
    lo=1+rank*q+min(rank,r)
    hi=lo+q+int(rank<r)-1
    ref.need(m and tuple(map(int,m.groups()))==(rank,ref.RANKS,lo,hi),str(path)+': wrong domain or deactivated mode')
    corners=[s for s in stdout if 'MPI CORNER PLAN' in s]
    if not corners:
        return None
    ref.need(rank==0 and len(corners)==1,str(path)+': expected exactly one corner receipt on rank 0 only')
    c=re.fullmatch(r'\s*MPI CORNER PLAN LOCALBC=([TF]) KEEP_UA=([TF]) KEEP_VA=([TF])\s*',corners[0])
    ref.need(c is not None,str(path)+': malformed corner plan')
    return dict(zip(('LOCALBC','KEEP_UA','KEEP_VA'),c.groups()))


def run(root, revision, tag, kind, expected_sha, proof):
    d=root/'runs'/tag
    a,b=ref.read(d/'started.json'),ref.read(d/'result.json')
    def check(x,msg): ref.need(x,tag+': '+msg)
    check(all(k in b and b[k]==v for k,v in a.items()),'started/result disagree')
    size='short' if kind=='early' else 'long'
    check(a['ranks']==ref.RANKS and a['threads']==1 and a['steps']==(300 if size=='short' else 4800),'shape')
    check(a['root_only_large_outputs'] is True and a['cpus']==ref.CPU,'I/O, cores')
    check(len(a['input_sha'])==24 and ref.canon(a['input_sha'])==ref.INPUT_SHA[size],'inputs != original')
    sub='source/build/mpi' if kind=='science' else 'trace-'+kind
    binary=ref.REMOTE+revision+'/'+sub+'/s--pom.x'
    cmd=a['command']
    def opt(k):
        check(cmd.count(k)==1 and cmd.index(k)+1<len(cmd),'missing/duplicate '+k)
        return cmd[cmd.index(k)+1]
    check(cmd[0]=='/usr/lib64/openmpi/bin/mpirun' and '--nooversubscribe' in cmd and
          opt('-np')=='32' and opt('--ranks')=='32' and opt('--threads')=='1' and
          opt('--binary')==binary and opt('--output')==ref.REMOTE+'runs/'+tag and
          a['binary_sha']==expected_sha,'launch path, count, or digest')
    topo=a['topology']
    check(set(topo)==set(map(str,ref.CPU)) and
          len({(v['physical_package_id'],v['core_id']) for v in topo.values()})==ref.RANKS,'physical cores')
    check(b['passed'] is True and b['error'] is None and b['returncode']==0,'solver return')
    rankset=set(map(str,range(ref.RANKS)))
    check(set(b['rank_binding_samples'])==set(b['checks'])==rankset,'incomplete receipts')
    pids=set(); corner=None
    for i,c in enumerate(ref.CPU):
        rd=d/('rank%03d'%i)
        plan=mode_receipt(rd/'stdout.log',i)
        if i==0:corner=plan
        launch=ref.read(rd/'launch.json')
        obs=b['rank_binding_samples'][str(i)]
        masks=obs['masks']
        check((launch['rank'],launch['size'],launch['threads'],launch['cpus'],launch['binary'],launch['binary_sha'])==
              (i,ref.RANKS,1,[c],binary,expected_sha),'launcher rank '+str(i))
        check(launch['pid']==obs['pid'] and obs['pid'] not in pids and str(obs['pid']) in masks and
              obs['application_omp_threads']==1 and obs['os_threads_including_mpi']==len(masks) and
              all(mask==[c] for mask in masks.values()),'observed OS affinity rank '+str(i))
        pids.add(obs['pid'])
        q=b['checks'][str(i)]; names=ref.FILES if i==0 else ('gcmprt','heat.dat')
        check(q['outputs']=={k:proof['output_sha'][size][k] for k in names} and
              sorted(q['expected_output_files'])==sorted(names) and q['unchanged_inputs'] is True and
              q['matches_same_day_serial'] is True,'input/output receipt rank '+str(i))
    check((d/'trace-audit.json').is_file()==(kind!='science'),'diagnostic/timing separation')
    t,start=b['elapsed'],a['started']
    check(all(type(x) in (float,int) and math.isfinite(x) and x>0 for x in (t,start)),'clocks')
    return t,start,corner


def paired_r014(root, tags, times, intervals, proof, plan_sha):
    """Validate the predeclared crossover against actual receipts, not printed medians."""
    plan=root/'build/formal-series-plan.json'
    ref.need(plan_sha==R014_PLAN and ref.sha(plan)==plan_sha,'formal precommit file digest != independent pin')
    ctrl=['r014-control%02d'%k for k in range(1,6)]
    ref.need({p.name for p in (root/'runs').iterdir() if re.fullmatch(r'r014-control\d+',p.name)}
             ==set(ctrl),'missing or extra crossover control sample')
    expected=[x for k,(c,t) in enumerate(zip(ctrl,tags)) for x in ((c,t) if k%2==0 else (t,c))]
    p=ref.read(plan)
    ref.need(p['sequence']==[['r004b' if t in ctrl else 'r014final',t] for t in expected]
             and (p['steps'],p['ranks'],p['threads'])==(4800,32,1)
             and p['target_seconds']==ref.FROZEN/24
             and p['sample_rule']=='all five candidate full samples retained; interleaved AB/BA; no pilot in formal median',
             'formal plan content != adjudicated ten-trial design')
    control={}
    elf=ref.ELF['science'][1]  # independently frozen r004b, not its claimed file name
    ce=root/'build/control-executable'
    with ce.open('rb') as f:ref.need(f.read(4)==b'\x7fELF' and ref.sha(ce)==elf,'wrong actual parent control ELF')
    for tag in ctrl:
        elapsed,start,_=run(root,'r004b',tag,'science',elf,proof)
        control[tag]=elapsed;intervals[tag]=(start,start+elapsed)
    seq=sorted(expected,key=lambda k:intervals[k][0])
    ref.need(seq==expected,'actual crossover order differs from AB/BA precommit')
    for a,b in zip(seq,seq[1:]):
        ref.need(intervals[a][1]<=intervals[b][0],'overlapping crossover trials: '+a+', '+b)
    vals=list(control.values()); med=statistics.median(vals); newm=statistics.median(times.values())
    return dict(precommit_sha256=plan_sha,control_binary_sha256=elf,observed_order=seq,
                control_seconds=control,control_median_seconds=med,control_range_seconds=[min(vals),max(vals)],
                control_to_candidate_median_speedup=med/newm,candidate_median_wall_reduction=1-newm/med,
                candidate_pair_wins=sum(times[t]<control[c] for t,c in zip(tags,ctrl)))


def audit(args):
    root=args.evidence.resolve()
    prefix=args.prefix
    ref.need(prefix.endswith('-'),'--prefix must include trailing dash')
    lib=scoped(prefix)
    members=lib.inventory(root)
    pp=root/'reference/day19-serial-proof.json'
    ref.need(ref.sha(pp)==ref.PROOF_SHA,'serial proof digest')
    proof=ref.read(pp)
    ref.need(proof['both_verified'] is True and proof['date']=='20261119' and
             all(set(proof['output_sha'][s])==set(ref.FILES) for s in ('short','long')),'serial proof schema')
    expected={'science':args.science_sha,'early':args.early_sha,'late':args.late_sha}
    ref.need(len(set(expected.values()))==3,'diagnostic and timed ELF identities must be distinct')
    tags=[prefix+'full%02d'%i for i in range(1,6)]
    ref.need({p.name for p in (root/'runs').iterdir() if re.fullmatch(re.escape(prefix)+r'full\d+',p.name)}
             ==set(tags), 'missing/extra same-prefix full timing sample')
    times,intervals,traces,snapshot_modes,corners={},{},{},{},{}
    for tag in tags:
        t,s,corners[tag]=run(root,args.revision,tag,'science',expected['science'],proof)
        times[tag]=t;intervals[tag]=(s,s+t)
    for label in ('early','late'):
        tag=prefix+label
        t,s,corners[tag]=run(root,args.revision,tag,label,expected[label],proof)
        intervals[tag]=(s,s+t)
        traces[label]=lib.trace(root,label)
        b=ref.read(root/'build'/('trace-'+label+'.json'))
        snapshot_modes[label]=b.get('snapshot_mode','legacy_in_place')
        ref.need(b['kind']=='instrumented_not_for_timing' and b['binary_sha']==expected[label] and
                 b['frames']==[list(p) for p in ref.POSITIONS[label]] and b['fields']==ref.FIELDS,
                 'trace build mismatch: '+label)
    paired=regressions=None
    if args.revision=='r014final':
        paired=paired_r014(root,tags,times,intervals,proof,args.formal_plan_sha)
        from check_regress24 import audit as audit_regressions
        regressions=audit_regressions(root,args.revision,prefix,proof,expected['science'])
        for rec in regressions.values():
            intervals[rec['tag']]=(rec['started'],rec['started']+rec['elapsed_seconds'])
    seq=sorted(intervals,key=intervals.__getitem__)
    ref.need([k for k in seq if k in times]==tags,'full sample timestamps out of order')
    for a,b in zip(seq,seq[1:]):ref.need(intervals[a][1]<=intervals[b][0],'overlapping runs: '+a+', '+b)
    out={k:ref.sha(root/'outputs'/k) for k in ref.FILES}
    ref.need(out==proof['output_sha']['long'],'actual six files != frozen full serial')
    src={}
    for p in (root/'source').iterdir():
        ref.need(p.is_file() and not p.is_symlink(),'unexpected source entry: '+str(p))
        src[p.name]=ref.sha(p)
    sd=ref.canon(src)
    ref.need(src and sd==args.source_tree_sha,'actual source file SHA mapping != independently supplied tree digest')
    files={'science':root/'build/mpi-executable','early':args.early_elf or root/'build/trace-early-executable',
           'late':args.late_elf or root/'build/trace-late-executable'}
    actual,missing={},[]
    for k,p in files.items():
        if not p.is_file():
            missing.append(k)
            continue
        ref.need(not p.is_symlink(),'symlink is not an ELF witness: '+str(p))
        with p.open('rb') as f:ref.need(f.read(4)==b'\x7fELF','missing ELF header: '+str(p))
        actual[k]=ref.sha(p)
        ref.need(actual[k]==expected[k],'actual ELF digest mismatch: '+k)
    vals=list(times.values()); m=statistics.median(vals); slow=max(vals)
    threshold=ref.FROZEN/24
    nondestructive=all(v=='owned_copy_nondestructive' for v in snapshot_modes.values())
    corner_required=(args.require_corner_elision or args.revision.startswith('r014') or any(v is not None for v in corners.values()))
    corner_ok=all(v==dict(LOCALBC='T',KEEP_UA='F',KEEP_VA='F') for v in corners.values())
    return dict(accepted=not missing and nondestructive and (not corner_required or corner_ok) and m<=threshold,
        numerical_pass=True,all_three_actual_elf_pass=not missing,corner_elision_required=corner_required,
        corner_elision_pass=corner_ok if corner_required else None,root_corner_plans=corners,
        snapshot_modes=snapshot_modes,nondestructive_snapshot_pass=nondestructive,
        missing_actual_elf=missing,revision=args.revision,prefix=prefix,manifest_members=members,
        manifest_sha256=ref.sha(root/'MANIFEST.json'),source_files=len(src),source_tree_sha256=sd,
        supplied_elf_sha256=expected,actual_elf_sha256=actual,serial_proof_sha256=ref.PROOF_SHA,
        original_input_mapping_sha256=ref.INPUT_SHA,original_frame_mapping_sha256=ref.FRAME_SHA,
        output_sha256=out,traces=traces,mode_and_binding_receipts=ref.RANKS*(len(tags)+len(traces)),
        control_binding_receipts=ref.RANKS*5 if paired else 0,paired_control=paired,compatibility_regressions=regressions,
        timing_samples_seconds=times,median_seconds=m,range_seconds=[min(vals),slow],slowest_seconds=slow,
        frozen_original_seconds=ref.FROZEN,same_day_original_seconds=ref.SAME_DAY,
        target_seconds=threshold,frozen_median_24x_pass=m<=threshold,
        frozen_slowest_24x_pass=slow<=threshold,same_day_median_24x_pass=m<=ref.SAME_DAY/24,
        same_day_slowest_24x_pass=slow<=ref.SAME_DAY/24,
        frozen_median_speedup=ref.FROZEN/m,frozen_slowest_speedup=ref.FROZEN/slow,
        same_day_median_speedup=ref.SAME_DAY/m,same_day_slowest_speedup=ref.SAME_DAY/slow)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('evidence',type=Path)
    p.add_argument('--revision',required=True,type=namearg)
    p.add_argument('--prefix',required=True,type=namearg,help='exact candidate run prefix, including trailing dash')
    p.add_argument('--formal-plan-sha',type=digestarg,default=R014_PLAN,help='out-of-band SHA-256 of predeclared r014final crossover')
    p.add_argument('--require-corner-elision',action='store_true',
                   help='require T/F/F corner-plan receipts in all seven runs (automatic for r014 or when any receipt exists)')
    for name in ['science','early','late','source-tree']:p.add_argument('--'+name+'-sha',required=True,type=digestarg)
    for name in ['early','late']:p.add_argument('--'+name+'-elf',type=Path,help='actual trace ELF when not in bundle')
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.output and (a.evidence.resolve()==a.output.resolve() or a.evidence.resolve() in a.output.resolve().parents):
        p.error('report must be a fresh file outside immutable evidence')
    try:
        report=audit(a)
    except (ValueError,IndexError,KeyError,TypeError,OSError,struct.error) as e:
        report=dict(accepted=False,error=type(e).__name__+': '+str(e))
    result=json.dumps(report,indent=2,allow_nan=False)+'\n'
    if a.output:
        with a.output.open('x') as f:f.write(result)
    print(result,end='')
    return 0 if report['accepted'] else 1


if __name__=='__main__':
    raise SystemExit(main())
