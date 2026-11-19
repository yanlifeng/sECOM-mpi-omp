"""Exercise real archive code/branches locally. No MPI transport or solver.

This compares the baseline to candidate single-rank and to rank zero with a
synthetic full owner range and no-op communication. It is not a multi-rank QA.
The independently pinned v39 fixture supplies raw captures, not NetCDF output.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
BASE='1b58f33639538c2d214f14a8f9a575ceffed5267'
FLAGS='-O3 -fopenmp -fno-fast-math -ffp-contract=off -std=legacy -ffixed-line-length-132 -ffpe-trap=invalid,zero,overflow -fcheck=bounds'.split()

def sha(b):return hashlib.sha256(b).hexdigest()
def once(s,a,b):
    assert s.count(a)==1,a
    return s.replace(a,b)

def fixture(probe):
    assert sha(probe)=='662bd6e05c16448c37dfc36941f40e7aba5ef66e0320d1d0b4e678fb159f363f'
    src=once(probe.decode(),"      TOR='BAROTROPIC'",'''      CALL SMPI_INIT(IM,JM)
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) ILAYOUT
      CALL QA_ARC_SET(ILAYOUT)
      INTX=0
      TOR='BAROTROPIC'
      HYDTYPE='INTERNAL'
      BCTYPE='CLAMPED'
      TRACER='NEGLECT'
      WAVEDYN='NEGLECT'
      IUTSR=77
      CONSTSR=.TRUE.
      FSM=1.
      ART=95.
      GRAV=9.81
      D=14.
      UA=.753
      VA=-.1999''')
    src=once(src,'      SKILLI=1.0/FLOAT(ISKILL)', '''      IF(MODE.EQ.10.OR.MODE.EQ.11) ISKILL=3
      IF(MODE.EQ.11) THEN
        EPTS=1
        INXIE(1)=5
        INXJE(1)=JM-1
        VPTS=1
        INXIV(1)=7
        INXJV(1)=JM/2
        FPTS=2
        ISFLX(1:2)=3
        JSFLX(1:2)=4
        DIRFLX(1)='IDIR'
        DIRFLX(2)='JDIR'
        NFLXE(1:2)=7
      ENDIF
      SKILLI=1.0/FLOAT(ISKILL)''')
    src=once(src,'      DO INTX=1,9','''! Exercise the unchanged cold-start INTX=0 path first.
      TIME=0.
      CALL ARCHIVE
      DO INTX=1,9''')
    src=once(src,'!@@CERT@@','''      INTX=2
      PROOF=SECOM_HIST_CERT0()
      IF(PROOF.NEQV.(MODE.EQ.0.OR.MODE.EQ.5.OR.MODE.EQ.9.OR.MODE.EQ.10.OR.MODE.EQ.11)) STOP 92
      INTX=0''')
    src=once(src,'!@@CALL@@','        CALL SECOM_ARCHIVE_BODY(PROOF)')
    src=once(src,'      STOP 94','''      INCLUDE 'comdeck'
      WRITE(23) TMIDDLE,ESAVE(1),UZSAVE(1,1),VZSAVE(1,1),CCFLUX(1:2,1),TXSAVE(1)
''')
    # This routine deliberately has no comdeck include: its policy common
    # interoperates with the future runtime declaration in mpideck.
    src+='''      SUBROUTINE QA_ARC_SET(ISIM)
      INTEGER ISIM
      LOGICAL ACT
      COMMON /SMPI_OUTER_FLAG/ ACT
      INTEGER R,N,LO,HI
      COMMON /SECOM_MPISLAB/ R,N,LO,HI
      ACT=ISIM.NE.0
      IF(ACT) N=2
      END
'''
    return src

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE.parents[1]/'source')
    p.add_argument('--fixture',type=Path,default=HERE/'fixture-v39.f')
    p.add_argument('--fc',default='gfortran');p.add_argument('--result',type=Path)
    a=p.parse_args()
    old=subprocess.check_output(['git','-C',str(a.source.parent),'show',BASE+':source/archive.f'])
    probe=a.fixture.read_bytes();new=(HERE/'archive.f').read_bytes()
    assert sha(old)=='e38b83224cfb10b38ec1eca494b6fff291966371110e1330ba60f0f1a5d4041b'
    src=fixture(probe)
    rows={};gold={}
    with tempfile.TemporaryDirectory(prefix='smpi-archive-local-') as t:
        d=Path(t)
        (d/'comdeck').write_text(once((a.source/'comdeck').read_text(),'IM=372,JM=479','IM=21,JM=19'))
        shutil.copyfile(a.source/'mpideck',d/'mpideck')
        (d/'driver.f').write_text(src)
        subprocess.run([a.fc,'-cpp','-fopenmp','-J.','-c',str((a.source/'secom_mpi.f90').resolve()),'-o','rt.o'],cwd=d,check=True)
        for role,n,data in [('v46',0,old),('off',0,new),('full_owner_no_transport',1,new)]:
            (d/'archive.f').write_bytes(data)
            subprocess.run([a.fc,*FLAGS,'-I.','archive.f','driver.f','rt.o','-o','exe'],cwd=d,check=True)
            for mode in [0,1,2,5,6,8,9,10,11]:
                sub=d/(role+str(mode));sub.mkdir()
                q=subprocess.run([str(d/'exe'),str(mode),str(n)],cwd=sub,env=dict(os.environ,OMP_NUM_THREADS='3',
                         OMP_DYNAMIC='FALSE'),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
                assert q.returncode==(95 if mode==9 else 0),(role,mode,q.returncode,q.stdout)
                captures={f.name:sha(f.read_bytes()) for f in sub.iterdir() if f.is_file()}
                assert captures or mode==9,(role,mode)
                r=dict(rc=q.returncode,outputs=captures)
                if role=='v46':gold[mode]=r
                assert r==gold[mode],(role,mode,r,gold[mode])
                rows[f'{role}/{mode}']=r
            print(role, 'matched',len(gold),'modes',flush=True)
    report=dict(kind='local reference vs policy-off/full-owner with NO MPI transport',compiler=subprocess.check_output(
         [a.fc,'--version']).decode().splitlines()[0],flags=FLAGS,source_sha=sha(new),cases=rows)
    if a.result:
        with a.result.open('x') as f:json.dump(report,f,indent=2)
if __name__=='__main__':main()
