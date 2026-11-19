"""Local isolated EXTRNL/ADVAVE test, never links or starts the solver."""
import argparse
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPO = ROOT.parent / 'secom-omp'
PROBE = ROOT/'v34-local/probe.f'
sys.path.insert(0, str(HERE.parent))
from prepare import COMMIT, generate, once

def location_args(parser):
    parser.add_argument('--repo',type=Path,help='Git repo containing pinned v44 commit (unless SECOM45_QA_SOURCE is set)')
    parser.add_argument('--probe',type=Path,help='path to the archived v34 probe.f fixture')

def use_locations(args):
    global REPO,PROBE
    if args.repo:REPO=args.repo.resolve()
    if args.probe:PROBE=args.probe.resolve()

def probe_text():
    data=PROBE.read_bytes()
    check(hashlib.sha256(data).hexdigest()=='8e68d6a27ab832e89d6a712bfa6668eee785d72937da781f7d20b9ef6fa811d0',
          'wrong v34 probe.f')
    return data.decode()

def base(name):
    source=os.environ.get('SECOM45_QA_SOURCE')
    raw=(Path(source)/name).read_bytes() if source else subprocess.check_output(
        ['git','-C',str(REPO),'show',COMMIT+':source/'+name])
    pins={'extrnl.f':'8099efb76e98f44b69af1b5a830d87497592ab835cc6bb68c88168119545b7d3',
          'advave.f':'4d2a44c8b944d691fd5e4cd5e0846ce4a6df73c3eb0db52dbf841c38f8d749aa',
          'comdeck':'5bc82dedffb8cf4f2781e39b9e0bd078b8b3f3a63bbe19295d62b808c5fae564'}
    check(hashlib.sha256(raw).hexdigest()==pins[name],'wrong '+name)
    return raw.decode()

def run(cmd, cwd, **kw):
    p=subprocess.run(cmd,cwd=cwd,capture_output=True,text=True,timeout=60,**kw)
    return p.returncode,p.stdout+p.stderr

def check(test, why):
    if not test: raise AssertionError(why)

def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--fc',default='gfortran')
    location_args(arg)
    opts=arg.parse_args()
    use_locations(opts)
    flags=[opts.fc,'-O3','-fopenmp','-ffp-contract=off','-fno-fast-math','-std=legacy',
           '-ffixed-line-length-132','-ffpe-trap=invalid,zero,overflow','-fcheck=bounds']
    ext=base('extrnl.f'); cand=generate(ext); witness=generate(ext,True)
    probe=probe_text()
    probe=once(probe,'      DO IQ=0,23','      DO IQ=0,23\n        INTX=IQ+1')
    probe=once(probe,'        DO KQSTEP=1,2', '''        COR(1,1)=TRANSFER(INT(Z'80000000'),COR(1,1))
        IF(IQ.EQ.1) COR(1,1)=TRANSFER(INT(Z'7fc00123'),COR(1,1))
        IF(IQ.EQ.2) COR(1,1)=TRANSFER(INT(Z'00000001'),COR(1,1))
        IF(IQ.EQ.3) COR(1,1)=TRANSFER(INT(Z'01800000'),COR(1,1))
        DO KQSTEP=1,2''')
    # The original fixture already has atomic publication/observation checks.
    # Delay the SINGLE writer, publish immediately after BCOND, require at mask.
    delayed=once(witness,'      CALL BCOND(12,DTI2,0)', '''      IF(FAST45) CALL DELAY34
      CALL BCOND(12,DTI2,0)
      IF(FAST45) CALL PUBLISH34''')
    anchor='      ENDIF\n!$OMP DO SCHEDULE(STATIC)\n      DO J=1,JM\n        DO I=1,IM\n          ELF(I,J)=ELF(I,J)*FSM(I,J)'
    delayed=once(delayed,anchor,anchor.replace('      ENDIF\n','      ENDIF\n      IF(FAST45) CALL REQUIRE34\n'))
    mutant=once(delayed,'!$OMP END DO\n      ELSE\n!$OMP BARRIER','!$OMP END DO NOWAIT\n      ELSE\n!$OMP BARRIER')
    # Exact source predicate. No floating-point operation on the fixture inputs.
    guard=cand[cand.index('      LOGICAL FUNCTION SECOM45_COR_OK'):]
    with tempfile.TemporaryDirectory(prefix='secom45-safety-') as tmp:
        d=Path(tmp)
        (d/'guard.f').write_text(guard)
        def build(src, output, extra=()):
            rc,txt=run([*flags,*extra,'-I.',*src,'-o',output],d)
            check(rc==0,txt)
        build([str(HERE/'guard.f'),'guard.f'],'guard')
        env=dict(os.environ,SECOM_V45_OVERLAP='1',OMP_DYNAMIC='FALSE')
        for setting in ('1','0','10',''):
            r,log=run([str(d/'guard'),'yes' if setting=='1' else 'no'],d,
                      env=dict(env,SECOM_V45_OVERLAP=setting))
            check(r==0,(setting,log))
        env.pop('SECOM_V45_OVERLAP')
        r,log=run([str(d/'guard'),'yes'],d,env=env);check(r==0,log)
        # With a different default real ABI, even healthy inputs must fail closed.
        build([str(HERE/'guard.f'),'guard.f'],'wide',['-fdefault-real-8'])
        r,log=run([str(d/'wide'),'no'],d,env=env);check(r==0,log)
        print('PASS: binary32 guard (sNaN included, traps enabled), opt-out, absent env, real-8 fallback',flush=True)
        (d/'comdeck').write_text(once(base('comdeck'),'IM=372,JM=479','IM=61,JM=67'))
        (d/'advave.f').write_text(base('advave.f'))
        (d/'probe.f').write_text(probe)
        digest=None
        for role,s in [('baseline',ext),('candidate',cand),('witness',delayed),('broken',mutant)]:
            (d/'unit.f').write_text(s)
            build(['probe.f','advave.f','unit.f'],role)
            for n in ([3] if role=='broken' else [1,3]):
                for setting in (['1','0'] if role=='witness' else ['1']):
                    sub=d/f'{role}-{n}-{setting}';sub.mkdir()
                    rc,log=run([str(d/role),str(n)],sub,env=dict(env,OMP_NUM_THREADS=str(n),SECOM_V45_OVERLAP=setting))
                    if role=='broken':
                        check(rc==94,('missing barrier mutant survived',rc,log))
                        print('PASS: removing early END DO join produces STOP 94',flush=True);continue
                    check(rc==0,(role,n,rc,log))
                    if role=='witness':
                        seq=[tuple(map(int,x)) for x in re.findall(r'V45_STEP\s+(\d+)\s+(\d+)\s+([01])',log)]
                        expected=[(iq+1,k,int(setting=='1' and iq in (0,3) and k==2)) for iq in range(24) for k in (1,2)]
                        check(seq==expected,('wrong eligibility/path',seq))
                    f=sub/'fields.bin'
                    check(f.stat().st_size==61*67*4*24*2*20,'partial state')
                    actual=hashlib.sha256(f.read_bytes()).hexdigest()
                    if digest is None:digest=actual
                    check(actual==digest,('state mismatch',role,n,actual,digest))
                    print(role,n,setting,actual,flush=True)
        print('PASS: 20 bit-exact fields, 24 configurations x 2 substeps, 1 and 3 threads; OMP publication checked')
if __name__=='__main__':main()
