"""Small real eight-step MPI corner publication proof, with a false-UA-skip control."""
import argparse
import hashlib
import json
import subprocess as sp
from pathlib import Path
from entry import HERE,once
from test_entry import stage as init,FLAGS
from corner import build as overlay

def stage(d,src,role):
    init(d,src,'old');s=(d/'driver.f').read_text()
    if role!='old':
        for f,t in overlay(src).items():(d/f).write_bytes(t.encode('latin1'))
        if role=='bad-skip-ua':
            f=d/'extrnl.f';s0=f.read_text(encoding='latin1');f.write_bytes(once(s0,'          IF(SMPI_KEEP_UA) CALL SMPI_HALO(UA,IM,JM)','! QA mutation: discard even the required outgoing UA update').encode('latin1'))
    s=once(s,'      FSM(4,2)=0.','''      FSM(4,2)=0.
      IF(LAYOUT.EQ.4) JCHARD(1)=JE-1
      IF(LAYOUT.EQ.5) JCHARD(1)=JE+1
      IF(LAYOUT.EQ.6) NUMECR=0
! Trigger an observable corner correction at the seam.
      IF(NUMECR.GT.0) THEN
        UA(7,JCHARD(1))=3.17
        UAF(7,JCHARD(1))=2.23
        DUM(7,JCHARD(1))=1.
        DVM(7,MIN(JM,JCHARD(1)+1))=1.
      ENDIF
      COR=.13
      CALL SECOM_MPI_BC_PLAN()
      SMPI_OUTER=SMPI_LOCALBC.AND.LAYOUT.NE.3
''')
    if role!='old':
        s=once(s,'      COR=.13','      COR=.13')
        s=once(s,'      SMPI_OUTER=SMPI_LOCALBC.AND.LAYOUT.NE.3','''      IF(SMPI_SIZE.GT.1) THEN
        IF(LAYOUT.EQ.5) THEN
          IF(SMPI_LOCALBC.OR..NOT.SMPI_KEEP_UA.OR..NOT.SMPI_KEEP_VA) STOP 85
        ELSE
          IF(.NOT.SMPI_LOCALBC.OR.SMPI_KEEP_VA) STOP 86
          IF(SMPI_KEEP_UA.NEQV.(LAYOUT.EQ.1.OR.LAYOUT.EQ.3)) STOP 87
        ENDIF
      ENDIF
      SMPI_OUTER=SMPI_LOCALBC.AND.LAYOUT.NE.3''')
    if 'SUBROUTINE SECOM_WAD_FIVE_TEAM' in (d/'wandd_kky.f').read_text():
        s=once(s,'        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN\n          CALL SECOM_WAD_COLOR_TEAM(EL,D)',
            '        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC.AND.SMPI_OUTER) THEN\n          CALL SECOM_WAD_FIVE_TEAM(EL,D)\n        ELSEIF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN\n          CALL SECOM_WAD_COLOR_TEAM(EL,D)')
    s=once(s,'        CALL SECOM_CHECK_TEAM(VMAX,EMIN,IIT,JIT,RV,RE)','        CALL SECOM_CHECK_READY_TEAM(VMAX,EMIN,IIT,JIT,RV,RE,SMPI_LOCALBC)')
    (d/'driver.f').write_bytes(s.encode('latin1'));p=d/'inputs.json';v=json.loads(p.read_text());v['role']=role;v['files']={f:hashlib.sha256((d/f).read_bytes()).hexdigest() for f in v['files']};p.write_text(json.dumps(v,indent=2)+'\n')

def build(d,fc,smoke):
    files=json.loads((d/'inputs.json').read_text())['files'];assert all(hashlib.sha256((d/f).read_bytes()).hexdigest()==h for f,h in files.items())
    sp.run([fc,*FLAGS,'-cpp',*([] if smoke else ['-DSECOM_USE_MPI']),'-I.','-J.','-c','secom_mpi.f90','-o','rt.o'],cwd=d,check=True)
    sp.run([fc,*FLAGS,'-I.','-J.','extrnl.f','advave.f','wandd_kky.f','baro_check.f','two_bconds.f','mpi_boundary.f','driver.f','rt.o','-o','unit'],cwd=d,check=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['stage','build']);p.add_argument('directory',type=Path);p.add_argument('--source',type=Path,default=HERE.parents[1]/'source');p.add_argument('--role',choices=['old','cert','bad-skip-ua'],default='old');p.add_argument('--fc',default='mpifort');p.add_argument('--smoke',action='store_true');a=p.parse_args()
    if a.action=='stage':stage(a.directory,a.source,a.role)
    else:build(a.directory,a.fc,a.smoke)
