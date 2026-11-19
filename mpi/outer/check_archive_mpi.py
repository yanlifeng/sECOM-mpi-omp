"""Validate real owned ARCHIVE+MPI with stale non-owner rows and event captures.

MPI mode is for authorized OpenMPI allocations only; --smoke needs no MPI.
Input: frozen v39 fixture + frozen v46 archive. No solver or NetCDF required.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess as sp
import sys
import tempfile
from check_archive_local import BASE,FLAGS,HERE,fixture,once,sha

ACC2='ARCET ARCPATM ARCWSX2D ARCWSY2D ARCAIRT ARCCLD ARCRHUM ARCSWOBS ARCQPREC ARCQEVAP ARCWVHT ARCWVPD ARCWVDR'.split()
ACC3='ARCU ARCV ARCAM ARCUX ARCVX'.split()
SRC2='CBC VARWIF AREAICE QPREC2D QEVAP2D DUMWAD DVMWAD AAM2D'.split()
ACTIVE='EL UA VA UAF VAF D WSX2D WSY2D PATM CLDFRC2D AIRTM2D RELHU2D SWOBS2D WVPD WVHT WVDR TXSURF TYSURF'.split()+SRC2

def driver(probe):
    s=fixture(probe)
    s=once(s,'      IF(ACT) N=2','! Actual SMPI_INIT supplied rank, size and inclusive global slab.')
    s=once(s,'      CLOSE(23)','      CLOSE(23)\n      CALL SMPI_FINALIZE()')
    s=once(s,'        INXJV(1)=JM/2','        INXJV(1)=JM/3+1')
    s=once(s,'      JHM=3','      JHM=4')
    s=once(s,'1+(IH-1)*IAVGE','1+(IH-2)*IAVGE');s=once(s,'=IH*IAVGE','=(IH-1)*IAVGE')
    s=once(s,'      DAYI=1./86400.','      IHIST(1,1)=0\n      IHIST(1,2)=0\n      DAYI=1./86400.')
    # Mode 2: only the actual owner has the signed zero (other ranks see benign +0).
    s=once(s,'      INTX=2','''      IF(ILAYOUT.GT.0.AND.SMPI_JHI.LT.JM.AND.MODE.EQ.2)
     &  WSET1(3,JM,1)=0.
      INTX=2''')
    s=once(s,'      TIME=0.','      CALL QA_ARC_FIELDS(0,0)\n      TIME=0.')
    s=once(s,'        CALL SECOM_ARCHIVE_BODY(PROOF)','        CALL QA_ARC_FIELDS(INTX,ILAYOUT)\n        CALL SECOM_ARCHIVE_BODY(PROOF)')
    # Both reference and trial write only authoritative rows at non-event times.
    # This reads without publishing: the test does not repair stale state for the code.
    a=s.index('      WRITE(23) INTX,ARCU,');b=s.index('      END\n',a)
    s=s[:a]+'      WRITE(23) INTX\n'+''.join('      WRITE(23) '+x+'(:,SMPI_JLO:SMPI_JHI)\n' for x in ACC2+['ARCTHIK'])+''.join('      WRITE(23) '+x+'(:,SMPI_JLO:SMPI_JHI,:)\n' for x in ACC3+'ARCWS1 ARCWS2 ARCSED1 ARCSED2 ARCTAU ARCHEM1 ARCHEM2'.split())+s[b:]
    # At PUTCDF the real consumer sees full planes: unconditionally capture all 28.
    s=once(s,'      WRITE(23) TMIDDLE\n      CALL SNAP39',"      INCLUDE 'comdeck'\n      WRITE(23) TMIDDLE\n"+''.join('      WRITE(23) '+x+'\n' for x in ACC2+SRC2)+''.join('      WRITE(23) '+x+'(:,:,1)\n' for x in ACC3)+'      WRITE(23) AAM(:,:,2),KM(:,:,2)')
    s=once(s,'      WRITE(23) TMIDDLE,ESAVE(1),UZSAVE(1,1),VZSAVE(1,1),CCFLUX(1:2,1),TXSAVE(1)',
           '      WRITE(23) TMIDDLE,CCFLUX(1:2,1)\n'+''.join('      WRITE(23) '+x+'(1)\n' for x in 'ESAVE TXSAVE TYSAVE PATMSAVE WSX2DSAVE WSY2DSAVE CLDSAVE AIRTSAVE RHUMSAVE SWOBSSAVE QPRECSAVE QEVAPSAVE WVPDSAVE WVHTSAVE WVDRSAVE THSAVE'.split())+''.join('      WRITE(23) '+x+'(1,1)\n' for x in 'UZSAVE VZSAVE AMZSAVE WS1SAVE WS2SAVE C1SAVE C2SAVE TAUSAVE P1SAVE P2SAVE'.split()))
    return s+'''      SUBROUTINE QA_ARC_FIELDS(IT,LAY)
      INCLUDE 'comdeck'
      INTEGER IT,LAY
      DO J=1,JM
       DO I=1,IM
        X=I*.013+J*.17+IT*.7
        IF(LAY.GT.0.AND.(J.LT.SMPI_JLO.OR.J.GT.SMPI_JHI))
     &    X=-1000.-SMPI_RANK-IT
'''+''.join(f'        {x}(I,J)=X*{n+1}.0/11.0\n' for n,x in enumerate(ACTIVE))+'''        TAU(I,J,1)=X*31.
        AAM(I,J,2)=X*32.
        KM(I,J,2)=X*33.
       ENDDO
      ENDDO
      END
'''

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE.parents[1]/'source')
    p.add_argument('--fixture',type=Path,default=HERE/'fixture-v39.f');p.add_argument('--reference',type=Path)
    p.add_argument('--fc',default='mpifort');p.add_argument('--mpiexec',default='mpirun')
    p.add_argument('--ranks',type=int,default=3);p.add_argument('--threads',type=int,default=1)
    p.add_argument('--smoke',action='store_true');p.add_argument('--result',type=Path);p.add_argument('--cpus')
    p.add_argument('--skip-negative',action='store_true');a=p.parse_args()
    n=1 if a.smoke else a.ranks
    assert 1<=n<=19 and (a.smoke or n>=2)
    cpu=list(map(int,a.cpus.split(','))) if a.cpus else []
    if not a.smoke and (len(cpu)!=n*a.threads or len(set(cpu))!=len(cpu)):
        p.error('--cpus must specify one distinct, allocated physical CPU for each rank/thread')
    old=(a.reference.read_bytes() if a.reference else sp.check_output(['git','-C',str(a.source.parent),'show',BASE+':source/archive.f']))
    assert sha(old)=='e38b83224cfb10b38ec1eca494b6fff291966371110e1330ba60f0f1a5d4041b'
    new=(HERE/'archive.f').read_bytes(); raw=driver(a.fixture.read_bytes());rows={}
    def capture(d):return {str(f.relative_to(d)):sha(f.read_bytes()) for f in sorted(d.glob('rank*/*')) if f.is_file()}
    with tempfile.TemporaryDirectory(prefix='smpi-archive-poison-') as t:
        d=Path(t);src=a.source.resolve()
        (d/'comdeck').write_text(once((src/'comdeck').read_text(),'IM=372,JM=479','IM=21,JM=19'))
        shutil.copyfile(src/'mpideck',d/'mpideck');(d/'driver.f').write_text(raw)
        sp.run([a.fc,'-cpp',*(['-DSECOM_USE_MPI'] if not a.smoke else []),'-fopenmp','-J.','-c',str(src/'secom_mpi.f90'),'-o','rt.o'],cwd=d,check=True)
        # Launcher obeys inherited CPU allocation; it cannot select extra CPUs.
        (d/'worker.py').write_text('''import os,sys
from pathlib import Path
rank=int(os.environ.get('OMPI_COMM_WORLD_RANK','0'));cpus=sys.argv[4]
if cpus:
 group=list(map(int,cpus.split(','))); nt=int(sys.argv[5]); sel=set(group[rank*nt:(rank+1)*nt]); assert sel<=os.sched_getaffinity(0);os.sched_setaffinity(0,sel)
os.chdir('rank'+str(rank));os.execv(sys.argv[1],[sys.argv[1],sys.argv[2],sys.argv[3]])
''')
        roles=[('reference',old,0),('candidate',new,1)]
        if not a.smoke and not a.skip_negative:
            roles.append(('missing_CBC',once(new.decode(),'      CALL SMPI_ALL_RAW(CBC,IM,JM)','! QA mutation: omit global output CBC').encode(),1))
        for role,data,on in roles:
            (d/'archive.f').write_bytes(data)
            sp.run([a.fc,*FLAGS,'-I.','archive.f','driver.f','rt.o','-o','exe'],cwd=d,check=True)
            for mode in ([11] if role=='missing_CBC' else [2,10,11]):
                z=d/(role+str(mode));z.mkdir()
                for k in range(n):(z/('rank'+str(k))).mkdir()
                cmd=[sys.executable,str(d/'worker.py'),str(d/'exe'),str(mode),str(on),a.cpus or '',str(a.threads)]
                if not a.smoke:cmd=[a.mpiexec,'--bind-to','none','-np',str(n)]+cmd
                q=sp.run(cmd,cwd=z,env=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='FALSE',OMP_PROC_BIND='FALSE'),stdout=sp.PIPE,stderr=sp.STDOUT,timeout=45)
                assert q.returncode==0,(role,mode,q.returncode,q.stdout)
                files=capture(z); assert len(files)==(2*n if mode==10 else n),(role,mode,files)
                rows[role+'/'+str(mode)]=files
                if role!='reference':
                    mismatches=sum(files[k]!=v for k,v in rows[f'reference/{mode}'].items())
                    assert (mismatches==n if role=='missing_CBC' else mismatches==0),(role,mode,mismatches,files)
                print(role,mode,'PASS',flush=True)
    report=dict(kind='serial-build smoke only' if a.smoke else 'real MPI slab poison + raw branch captures',ranks=n,threads=a.threads,cpus=cpu,
          compiler=sp.check_output([a.fc,'--version']).decode().splitlines()[0],mpi=None if a.smoke else sp.check_output([a.mpiexec,'--version']).decode().splitlines()[0],flags=FLAGS,
          source_sha=sha(new),driver_sha=sha(raw.encode()),runtime_sha=sha((a.source/'secom_mpi.f90').read_bytes()),cases=rows)
    if a.result:
        with a.result.open('x') as f:json.dump(report,f,indent=2)
if __name__=='__main__':main()
