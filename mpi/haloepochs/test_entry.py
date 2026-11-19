"""Stage and build full tiny EXTRNL/WAD/CHECK epoch fixtures; never a solver."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess as sp
from entry import HERE,entry,wadinput,once,source

FLAGS='-O3 -fopenmp -fno-fast-math -ffp-contract=off -ffpe-trap=invalid,zero,overflow -fcheck=bounds -std=legacy -ffixed-line-length-132'.split()
FIELDS='D UA VA UAB VAB EL ELB DUMWAD DVMWAD UAF VAF ELF ET ETF CURV42D TXBOT TYBOT TXICE TYICE FLUXUA FLUXVA TPS AX AY'.split()

def stage(d,src,role):
    d.mkdir();files={f:source(src,f) for f in ['comdeck','mpideck','secom_mpi.f90','mpi_boundary.f','extrnl.f','advave.f','wandd_kky.f','baro_check.f']}
    if role in ['entry','bad-next-vab','bad-first-el']:
        files.update(entry(src))
        if role!='entry':
            # Intentionally discard the north VAB or south EL version just
            # produced. An actually-read packet component must catch this.
            n='     &                        DVMWAD,IM,JM,IEXT.EQ.1)\n'
            f='VAB' if role=='bad-next-vab' else 'EL';cond='IEXT.GT.1.AND.SMPI_JHI.LT.JM' if f=='VAB' else 'IEXT.EQ.1.AND.SMPI_JLO.GT.1';j='SMPI_JHI+1' if f=='VAB' else 'SMPI_JLO-1'
            files['extrnl.f']=once(files['extrnl.f'],n,n+'!$OMP SINGLE\n        IF('+cond+') '+f+'(:,'+j+')=91.756\n!$OMP END SINGLE\n')
    if role in ('wadinput','bad-wad-mask'):
        files.update({f:t for f,t in wadinput(src).items() if f!='ecom3d.f'})
        if role=='bad-wad-mask':
            f=files['wandd_kky.f'];p=f.index('! WAD_INPUT_EPOCH:');f=f[:p]+'''!$OMP SINGLE
        IF(SMPI_JHI.LT.JM) DVMWAD(:,SMPI_JHI+1)=0.
!$OMP END SINGLE
'''+f[p:];files['wandd_kky.f']=f
    files['comdeck']=once(files['comdeck'],'IM=372,JM=479','IM=372,JM=19')
    # Use the two actual zero-discharge LINEAR branch bodies from BCOND.
    b=source(src,'bcond.f');b=b[b.index(' 10   CONTINUE'):b.index('!-------- INTERNAL VELOCITY BOUNDARY CONDITIONS')]
    assert b.count('IF(IDX.EQ.12) RETURN')==1 and b.count('IF(IDX.EQ.13) RETURN')==1
    files['two_bconds.f']="      SUBROUTINE BCOND(IDX,DTI2,IZERO)\n      INCLUDE 'comdeck'\n      IF(IDX.EQ.12) GOTO 10\n      IF(IDX.EQ.13) GOTO 20\n      STOP 73\n"+b+'      END\n'
    rs=files['wandd_kky.f'];s=(HERE/'entry_driver.f').read_text()
    routine=('SECOM_WAD_COLOR_TEAM' if 'SUBROUTINE SECOM_WAD_COLOR_TEAM' in rs else 'SECOM_WAD_INREADY_TEAM' if 'SUBROUTINE SECOM_WAD_INREADY_TEAM' in rs else 'WANDD_TEAM')
    s=once(s,'!@@WAD@@','''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL '''+routine+'''(EL,D)
        ELSE
          CALL WANDD_TEAM(EL,D)
        ENDIF''')
    s=once(s,'!@@CAPTURE@@','\n'.join(f'        WRITE(22) {f}(:,SMPI_JLO:SMPI_JHI)' for f in FIELDS))
    files['driver.f']=s
    dig={}
    for name,t in files.items():
        b=t.encode('latin1');(d/name).write_bytes(b);dig[name]=hashlib.sha256(b).hexdigest()
    (d/'inputs.json').write_text(json.dumps(dict(role=role,files=dig),indent=2)+'\n')

def build(d,fc,smoke):
    sig=json.loads((d/'inputs.json').read_text())['files'];assert all(hashlib.sha256((d/f).read_bytes()).hexdigest()==v for f,v in sig.items())
    sp.run([fc,*FLAGS,'-cpp',*([] if smoke else ['-DSECOM_USE_MPI']),'-J.','-I.','-c','secom_mpi.f90','-o','rt.o'],cwd=d,check=True)
    sp.run([fc,*FLAGS,'-I.','-J.','extrnl.f','advave.f','wandd_kky.f','baro_check.f','two_bconds.f','driver.f','rt.o','-o','unit'],cwd=d,check=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['stage','build']);p.add_argument('directory',type=Path);p.add_argument('--source',default=HERE.parents[1]/'source',type=Path);p.add_argument('--role',choices=['old','entry','wadinput','bad-next-vab','bad-first-el','bad-wad-mask'],default='old');p.add_argument('--fc',default='mpifort');p.add_argument('--smoke',action='store_true');a=p.parse_args()
    if a.action=='stage':stage(a.directory,a.source,a.role)
    else:build(a.directory,a.fc,a.smoke)
