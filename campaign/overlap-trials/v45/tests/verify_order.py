"""Negative controls for v44 deferred-edge publication and original sNaN trap order."""
import os
import hashlib
import argparse
from pathlib import Path
import signal
import tempfile
from verify import base, check, generate, once, run, location_args, use_locations, probe_text

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fc',nargs='?',default='gfortran')
    location_args(parser)
    args=parser.parse_args()
    fc=args.fc
    use_locations(args)
    flags=[fc,'-O3','-fopenmp','-ffp-contract=off','-fno-fast-math','-std=legacy',
           '-ffixed-line-length-132','-ffpe-trap=invalid,zero,overflow','-fcheck=bounds']
    ext=base('extrnl.f')
    cand=generate(ext)
    adv=base('advave.f')
    probe=probe_text()
    print('compiler:',run([fc,'--version'],Path.cwd())[1].splitlines()[0],flush=True)
    print('flags:',*flags[1:],flush=True)
    print('CPU affinity:',sorted(os.sched_getaffinity(0)) if hasattr(os,'sched_getaffinity') else 'not available',flush=True)
    def digest(s):return hashlib.sha256(s.encode()).hexdigest()
    print('base/science/adv/driver SHA256:',*[digest(s) for s in (ext,cand,adv,probe)],flush=True)
    if os.environ.get('SECOM45_SCI_FILE'):
        check(Path(os.environ['SECOM45_SCI_FILE']).read_text()==cand,'remote science source mismatch')
    # Use the actual shared-team kernel but only its first mode. No production input.
    probe=once(probe,'      DO IQ=0,23','      DO IQ=0,0')
    with tempfile.TemporaryDirectory(prefix='secom45-order-') as tmp:
        d=Path(tmp)
        (d/'comdeck').write_text(once(base('comdeck'),'IM=372,JM=479','IM=61,JM=67'))
        def trial(name,ext_src,adv_src,driver,n):
            (d/'unit.f').write_text(ext_src);(d/'adv.f').write_text(adv_src)
            (d/'probe.f').write_text(driver)
            rc,log=run([*flags,'-I.','unit.f','adv.f','probe.f','-o','qa'],d);check(rc==0,log)
            print(name,'instrumented source/binary SHA256:',digest(ext_src),hashlib.sha256((d/'qa').read_bytes()).hexdigest(),flush=True)
            r=d/name;r.mkdir()
            return run([str(d/'qa'),str(n)],r,env=dict(os.environ,OMP_NUM_THREADS=str(n),
                            OMP_DYNAMIC='FALSE',SECOM_V45_OVERLAP='1'))
        bcp='      IF(IDX.EQ.12) ELF(2,JM-1)=ELF(3,JM-1)+.031'
        bad=once(probe,bcp,"""      IF(IDX.EQ.12) THEN
        WRITE(*,*) 'ENTERED_BC12'
        FLUSH(6)
      ENDIF
"""+bcp)
        bad=once(bad,'        DO KQSTEP=1,2',"""        COR(1,1)=TRANSFER(INT(Z'7F800001'),COR(1,1))
        IF(TRANSFER(COR(1,1),0).NE.INT(Z'7F800001')) STOP 81
        DO KQSTEP=1,2""")
        # Gfortran's bounds and FPE traps stay on. Only one lane to observe an
        # unambiguous order, independent of asynchronous multi-lane diagnostics.
        pos='!$OMP DO SCHEDULE(STATIC)\n      DO 95'
        original_trap="""      WRITE(*,*) 'ENTERED_OLD_DO95_AFTER_ADV'
      FLUSH(6)
"""
        corrupt=once(cand,'      IF(MODE45.AND.IEXT.EQ.1) OK45=',"""      IF(MODE45.AND.IEXT.EQ.1.AND.COR(1,1).GT.0) CALL PUBLISH34
      IF(MODE45.AND.IEXT.EQ.1) OK45=""")
        for role,x in [('v44',ext),('v45',cand),('float_guard_mutant',corrupt)]:
            rc,log=trial(role,once(x,pos,original_trap+pos),adv,bad,1)
            check(rc in (-signal.SIGFPE,128+signal.SIGFPE) and 'ENTERED_BC12' in log,
                  (role,rc,log))
            check(('ENTERED_OLD_DO95_AFTER_ADV' in log)==(role!='float_guard_mutant'),(role,log))
            print(role,'rc',rc,'SIGFPE at', 'pre-ADV guard' if role=='float_guard_mutant' else 'original DO95',flush=True)
        # Delay lane 0 BEFORE its final V divergence. Publish after that lane
        # has completed; verify from the ordered edge SINGLE, not from another DO.
        tail='!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J)\n      DO 880'
        late=once(adv,tail,"""      IF(EDGE_LATE.AND.IEXT.EQ.2) CALL DELAY_LANE45()
"""+tail)
        late=once(late,'      IF(EDGE_LATE) RETURN',"""      IF(EDGE_LATE.AND.IEXT.EQ.2) CALL PUBLISH_LANE45()
      IF(EDGE_LATE) RETURN""")
        observed=once(cand,'      IF(EDGE_LATE) CALL SECOM_ADV34_EDGES',"""      IF(FAST45) CALL REQUIRE34()
      IF(EDGE_LATE) CALL SECOM_ADV34_EDGES""")
        ty='!$OMP END DO\n!$OMP SINGLE\n! The existing stress join'
        unsafe=once(observed,ty,ty.replace('!$OMP END DO','!$OMP END DO NOWAIT'))
        driver=probe+"""      SUBROUTINE DELAY_LANE45()
      USE OMP_LIB
      IF(OMP_GET_THREAD_NUM().EQ.0) CALL DELAY34()
      END
      SUBROUTINE PUBLISH_LANE45()
      USE OMP_LIB
      IF(OMP_GET_THREAD_NUM().EQ.0) CALL PUBLISH34()
      END
"""
        for role,x,want in [('v45_ty_join',observed,0),('ty_nowait_mutant',unsafe,94)]:
            rc,log=trial(role,x,late,driver,3);check(rc==want,(role,rc,log))
            print(role,'exit',rc,flush=True)
if __name__=='__main__':main()
