"""v46 actual kernel + true persistent team; runs no solver. Linux --trap is strict."""
import argparse
import hashlib
import os
import re
import signal
import subprocess
import sys
import tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from prepare import COMMIT,generate,once
def sha(data):return hashlib.sha256(data).hexdigest()
def need(c,m):
    if not c:raise AssertionError(m)
def main():
    p=argparse.ArgumentParser(description=__doc__)
    src=p.add_mutually_exclusive_group(required=True)
    src.add_argument('--repo',type=Path);src.add_argument('--base',type=Path)
    p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--science-file',type=Path,required=True)
    p.add_argument('--fc',default='gfortran');p.add_argument('--trap',action='store_true')
    a=p.parse_args()
    def raw(n):
        return (a.base/n).read_bytes() if a.base else subprocess.check_output(
            ['git','-C',str(a.repo),'show',COMMIT+':source/'+n])
    old=raw('extrnl.f').decode();cand=generate(old);wit=generate(old,True)
    adv=raw('advave.f');com=raw('comdeck');probe=a.probe.read_bytes()
    need(cand.encode()==a.science_file.read_bytes() and sha(cand.encode())=='baf5718f127fa3da4ff808398c3a2a8250f22a4b907a5d59efab0beb129bfe47','wrong science')
    need([sha(x) for x in [adv,com,probe]]==['4d2a44c8b944d691fd5e4cd5e0846ce4a6df73c3eb0db52dbf841c38f8d749aa',
        '5bc82dedffb8cf4f2781e39b9e0bd078b8b3f3a63bbe19295d62b808c5fae564',
        '8e68d6a27ab832e89d6a712bfa6668eee785d72937da781f7d20b9ef6fa811d0'],'wrong frozen fixtures')
    driver=once(probe.decode(),'      DO IQ=0,23','      DO IQ=0,23\n        INTX=IQ+1')
    driver=once(driver,'''        DO KQSTEP=1,2
          IEXT=KQSTEP
          QFLAG=0
!$OMP PARALLEL DEFAULT(SHARED)''', '''        COR(1,1)=TRANSFER(IBSET(0,31),COR(1,1))
        IF(IQ.EQ.1) COR(1,1)=TRANSFER(INT(Z'7FC00101'),COR(1,1))
        IF(IQ.EQ.2) COR(1,1)=TRANSFER(INT(Z'00000001'),COR(1,1))
        IF(IQ.EQ.3) COR(1,1)=TRANSFER(INT(Z'7F800000'),COR(1,1))
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(KQSTEP)
        DO KQSTEP=1,2
!$OMP MASTER
          IEXT=KQSTEP
          QFLAG=0
!$OMP END MASTER''')
    driver=once(driver,'!$OMP END PARALLEL\n          WRITE(91)','!$OMP SINGLE\n          WRITE(91)')
    driver=once(driver,'        ENDDO\n      ENDDO\n      CLOSE(91)', '''!$OMP END SINGLE
        ENDDO
!$OMP END PARALLEL
      ENDDO
      CLOSE(91)''')
    flags=[a.fc,'-O3','-fopenmp','-std=legacy','-ffixed-line-length-132','-fno-fast-math',
           '-ffp-contract=off','-ffpe-trap=invalid,zero,overflow','-fcheck=bounds']
    print('compiler',subprocess.check_output([a.fc,'--version']).decode().splitlines()[0],
          'affinity',sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else 'unavailable',
          'flags',flags,'science',sha(cand.encode()),flush=True)
    env=dict(os.environ,OMP_DYNAMIC='FALSE')
    with tempfile.TemporaryDirectory(prefix='micro46-') as tmp:
        d=Path(tmp);(d/'comdeck').write_text(once(com.decode(),'IM=372,JM=479','IM=61,JM=67'))
        (d/'advave.f').write_bytes(adv)
        def build(s,drv=driver):
            (d/'ext.f').write_text(s);(d/'drv.f').write_text(drv)
            subprocess.run(flags+['ext.f','advave.f','drv.f','-I.','-o','qa'],cwd=d,check=True)
        def trial(tag,nt,setting):
            dest=d/tag;dest.mkdir()
            q=subprocess.run([str(d/'qa'),str(nt)],cwd=dest,env=dict(env,OMP_NUM_THREADS=str(nt),SECOM_V46_OVERLAP=setting),
                stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
            return q.returncode,q.stdout,dest
        gold=None
        for role,s in [('v44',old),('science',cand),('witness',wit)]:
            build(s)
            for n in (1,3):
                for setting in (('1','0') if role=='witness' else ('1',)):
                    rc,log,dest=trial(role+str(n)+setting,n,setting);need(rc==0,(rc,log))
                    blob=(dest/'fields.bin').read_bytes();need(len(blob)==15694080,'partial output')
                    h=sha(blob)
                    if gold is None:gold=h
                    need(h==gold,('unequal fields',role,n,h,gold))
                    if role=='witness':
                        seq=[tuple(map(int,z)) for z in re.findall(r'V46_STEP\s+(\d+)\s+(\d+)\s+([01])',log)]
                        need(seq==[(i+1,k,int(i<4 and k==2 and setting=='1')) for i in range(24) for k in (1,2)],seq)
                    print(role,n,setting,h,flush=True)
        delayed=once(cand,'      CALL BCOND(12,DTI2,0)','''      IF(FAST46) CALL DELAY34
      CALL BCOND(12,DTI2,0)
      IF(FAST46) CALL PUBLISH34''')
        pat='      ENDIF\n!$OMP DO SCHEDULE(STATIC)\n      DO J=1,JM\n        DO I=1,IM\n          ELF(I,J)=ELF(I,J)*FSM(I,J)'
        delayed=once(delayed,pat,pat.replace('      ENDIF\n','      ENDIF\n      IF(FAST46) CALL REQUIRE34\n'))
        mutant=once(delayed,'!$OMP END DO\n      ELSE\n!$OMP BARRIER','!$OMP END DO NOWAIT\n      ELSE\n!$OMP BARRIER')
        for role,s,expected in [('delay_join',delayed,0),('missing_join',mutant,94)]:
            build(s);rc,log,_=trial(role,3,'1');need(rc==expected,(rc,log));print(role,'rc',rc,flush=True)
        if a.trap:
            bad=once(driver,"IBSET(0,31)","INT(Z'7F800001')")
            pos='!$OMP DO SCHEDULE(STATIC)\n      DO 95'
            mark="      WRITE(*,*) 'AFTER_ADV_BEFORE_OLD_DO95'\n      FLUSH(6)\n"
            mut=once(cand,'      CALL BCOND(12,DTI2,0)','      IF(COR(1,1).GT.0) CALL PUBLISH34\n      CALL BCOND(12,DTI2,0)')
            for role,s in [('old',old),('new',cand),('premature',mut)]:
                build(once(s,pos,mark+pos),bad)
                rc,log,_=trial('trap_'+role,1,'1')
                need(rc==-signal.SIGFPE and ('AFTER_ADV_BEFORE_OLD_DO95' in log)==(role!='premature'),(rc,log))
                print('sNaN',role,'rc',rc,'phase',('early' if role=='premature' else 'original DO95'),flush=True)
if __name__=='__main__':main()
