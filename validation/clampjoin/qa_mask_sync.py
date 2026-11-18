"""Standalone v30 QA, not a timing harness. Default: read immutable local Git v27.

    python3 validation/clampjoin/qa_mask_sync.py
    python3 qa_mask_sync.py --v27-source /path/to/v27 --threads 1 3

Tests the actual extracted mask-through-ADVAVE synchronization window, with
unaltered v27 ADVAVE/real special-clamp arithmetic and QA-only scheduler
perturbations. No NetCDF or SSH. Full solver regression is separate.
"""
import argparse
import hashlib
import os
from pathlib import Path
import resource
import signal
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
BASE = 'ed2defca7b14a40ede9b5869e67ae2cd5f144f8d72a7f47c1d3f3a9bba50a16f'


def once(s, old, new):
    if s.count(old) != 1:
        raise ValueError('source mismatch: ' + repr(old[:80]))
    return s.replace(old, new)


def window(src):
    mask = '          ELF(I,J)=ELF(I,J)*FSM(I,J)'
    marker = '      CALL BCOND(12,DTI2,0)\n!$OMP END SINGLE\n'
    if src.count(marker) != 1 or src.count('      CALL ADVAVE(ADVUA,ADVVA)') != 1:
        raise ValueError('unknown mask window')
    s = src.split(marker)[1].split('      CALL ADVAVE(ADVUA,ADVVA)')[0]
    s = once(s,'      DO J=1,JM\n','      DO J=1,JM\n        IF(J.EQ.JMM1) CALL STALL30\n')
    s = once(s,mask+'\n        ENDDO',mask+'\n        ENDDO\n        IF(J.EQ.JMM1) CALL MASKED30')
    return ("""      RECURSIVE SUBROUTINE WINDOW30(ADVUA,ADVVA,ROTATE)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM)
      LOGICAL ROTATE
""" + s+'      CALL ADVAVE(ADVUA,ADVVA)\n      END\n')


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--v27-source', type=Path)
    cli.add_argument('--threads', nargs='+',type=int,default=[1,3])
    cli.add_argument('--fc',default=os.environ.get('FC','gfortran'))
    a = cli.parse_args()
    if any(n < 1 or n > 46 for n in a.threads) or max(a.threads) < 2:
        cli.error('provide at least one thread count >=2; the compact QA grid supports 1..46')
    def read(f):
        if a.v27_source:
            return (a.v27_source/f).read_bytes().decode('ascii')
        return subprocess.check_output(['git','-C',str(HERE.parents[1]),
                                        'show','secom-omp-v27:source/'+f]).decode('ascii')
    old = read('extrnl.f')
    if hashlib.sha256(old.encode()).hexdigest()!=BASE:
        raise ValueError('expected original v27 EXTRNL')
    candidate = subprocess.check_output(['git','-C',str(HERE.parents[1]),'show',
                                         'b3284b5:source/extrnl.f']).decode('ascii')
    flags = [a.fc,'-O2','-fopenmp','-std=legacy','-ffixed-line-length-132',
             '-ffp-contract=off','-fno-fast-math','-ffpe-trap=invalid,zero,overflow','-fcheck=bounds']
    with tempfile.TemporaryDirectory(prefix='qa-v30-') as t:
        d=Path(t)
        # Verify patch against original; also compile the entire production
        # EXTRNL with untouched full dimensions, independently of the probes.
        (d/'source').mkdir()
        pf=d/'source'/'extrnl.f'
        pf.write_text(old)
        patch=subprocess.check_output(['git','-C',str(HERE.parents[1]),'diff','--no-ext-diff',
                                       'secom-omp-v27','b3284b5','--','source/extrnl.f'])
        subprocess.run(['git','apply','--no-index','-'],input=patch,cwd=d,check=True)
        if pf.read_bytes() != candidate.encode('ascii'):
            raise ValueError('delivered source and patch differ')
        (d/'comdeck').write_text(read('comdeck'))
        subprocess.run(flags+['-I.','-c','source/extrnl.f'],cwd=d,check=True)
        # Compact synthetic grid; all tested branches and real loop bodies stay intact.
        (d/'comdeck').write_text(once(read('comdeck'),'IM=372,JM=479','IM=37,JM=47'))
        adv = once(read('advave.f'),'      DO J=1,JM\n       DO I=1,IM\n      FLUXUA(I,J)=0.0\n'
                   '      FLUXVA(I,J)=0.0\n      ADVUA(I,J)=0.0',
                   '      DO J=1,JM\n       IF(J.EQ.1) CALL ENTER_ADV30\n       DO I=1,IM\n'
                   '      FLUXUA(I,J)=0.0\n      FLUXVA(I,J)=0.0\n      ADVUA(I,J)=0.0')
        (d/'advave.f').write_text(adv)
        for f,anchor in [('pclamp.f','        DO 100 N = 1, NUMEBC'),('oclamp.f','      CLAM=1.0')]:
            s=once(read(f),anchor,'      CALL BEFORE_CLAMP30\n'+anchor)
            if f=='pclamp.f':
                s=once(s,' 100  CONTINUE','      CALL CHECK_IEEE30(ELF(IE,JE))\n 100  CONTINUE')
            (d/f).write_text(s)
        case = {'reference':window(old), 'candidate':window(candidate)}
        case['badgate'] = once(case['candidate'],'!$OMP BARRIER\n','')
        for role, code in case.items():
            (d/'window.f').write_text(code)
            subprocess.run(flags+['-I.',str(HERE/'window-driver.f'),'window.f','advave.f',
                                  'pclamp.f','oclamp.f','-o',role],cwd=d,check=True)
        env=dict(os.environ,OMP_DYNAMIC='FALSE')
        def run(role,mode,rotate,nt,fast=False,exception='normal'):
            # Deliberate FP fault is part of negative QA; never generate a core dump.
            def nocore():
                resource.setrlimit(resource.RLIMIT_CORE,(0,0))
            return subprocess.run([str(d/role),mode,'rotate' if rotate else 'fallback',
                                   str(nt),str(int(fast)),exception],
                                  env=env,cwd=d,timeout=60,stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT,preexec_fn=nocore)
        total=nstates=nfast=0
        for nt in a.threads:
            for mode,rotate in [('CLAMPED',True),('CLAMPED',False),('OCLAMP',True),
                                ('IRANDB',True),('MIXED',True),('PCLAMP',True),('RANDB',True)]:
                gold=None
                for role in ('reference','candidate'):
                    fast=role=='candidate' and mode=='CLAMPED' and rotate
                    ret=run(role,mode,rotate,nt,fast)
                    if ret.returncode:
                        raise RuntimeError((role,mode,rotate,nt,ret.returncode,ret.stdout.decode(errors='replace')))
                    state=(d/'states.bin').read_bytes()
                    if len(state)!=10*6*37*47*4:
                        raise RuntimeError('missing state')
                    if gold is not None and state!=gold:
                        raise RuntimeError(('different state',mode,rotate,nt))
                    gold=state
                nfast+=10*int(nt>1 and mode=='CLAMPED' and rotate)
                total+=len(gold); nstates+=10
        nt=max(a.threads)
        # The same fixture must reject a deliberately broken fallback barrier.
        # Non-deterministic SINGLE election permits a late worker to win once;
        # retry only this detuned mutant, never a candidate or solver benchmark.
        for _ in range(3):
            bad=run('badgate','PCLAMP',True,nt)
            if bad.returncode==66: break
        else:
            raise RuntimeError(('negative control missed the race',bad.returncode,bad.stdout))
        for setting in ('invalid','trap'):
            rc=[]
            for role in ('reference','candidate'):
                ex=run(role,'PCLAMP',True,nt,exception=setting)
                rc.append(ex.returncode)
            # 67 = NaN with sticky invalid, 68 = NaN without observable
            # sticky invalid (possible under translation). Require BOTH
            # builds to agree on the signal or on the NaN AND flag; never
            # call a nontrapping result a demonstrated SIGFPE.
            if rc[0]!=rc[1] or rc[0] not in ([67,68] if setting=='invalid' else [-signal.SIGFPE,67,68]):
                raise RuntimeError(('PCLAMP invalid exception mismatch',setting,rc))
            print('PCLAMP',setting,'v27 and v30 return codes:',rc)
        print('PASS: full unit compiled; patch verified;',nstates,'substeps with six byte-equal fields,',
              total,'bytes;',nfast,'forced fast-path overlaps. Fallback race mutant rejected.')


if __name__=='__main__':
    main()
