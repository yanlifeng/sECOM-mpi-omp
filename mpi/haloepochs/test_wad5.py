"""Extend the real 372-column entry physics fixture with directed WAD probes."""
import argparse
import hashlib
import json
from pathlib import Path
from entry import HERE,once
from test_entry import stage as baseline,build
from wad5 import build as overlay

def stage(d,src,role):
    baseline(d,src,'old');s=(d/'driver.f').read_text();w=(d/'wandd_kky.f').read_text();rt=(d/'secom_mpi.f90').read_text()
    if role!='old':
        q=overlay(src);w=q['wandd_kky.f'];rt=q['secom_mpi.f90']
        mark='      RECURSIVE SUBROUTINE SECOM_WAD_FIVE_TEAM(ELEV,TD)';left,right=w.split(mark)
        if role=='bad-el-south':
            a='        CALL SMPI_WAD_PAIR(ELEV,TD,IM,JM,.TRUE.)';right=once(right,a,a+'\n!$OMP SINGLE\n        IF(SMPI_JLO.GT.1) ELEV(:,SMPI_JLO-1)=71.43\n!$OMP END SINGLE')
        if role=='bad-vaf-north':
            right=once(right,'        CALL SMPI_FLUX_HALO_N(VAF,IM,JM)','! wrong: omitted final VAF publication')
        w=left+mark+right
        # Actual source trace counts: three pair_dir calls + one south for
        # ELEV/D + directional single VAF. The flags Sendrecv is separate.
        rt=once(rt,'    call MPI_Sendrecv(wpair_send,','    call qa_wadcount()\n    call MPI_Sendrecv(wpair_send,')
    rt=once(rt,'  if(src==MPI_PROC_NULL)then\n    ! No out-of-range actual argument','  call qa_wadcount()\n  if(src==MPI_PROC_NULL)then\n    ! No out-of-range actual argument')
    rt+='''
subroutine qa_wadcount()
 implicit none
 integer nc
 common /qa_wad_calls/ nc
 nc=nc+1
end subroutine
'''
    s=once(s,'      LOGICAL READY','      INTEGER NWC\n      COMMON /QA_WAD_CALLS/ NWC\n      LOGICAL READY')
    old='''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL SECOM_WAD_COLOR_TEAM(EL,D)
        ELSE
          CALL WANDD_TEAM(EL,D)
        ENDIF'''
    before='''!$OMP SINGLE
        NWC=0
        DO J=1,JM
          IF(J.GE.SMPI_JLO.AND.J.LE.SMPI_JHI) CYCLE
          EL(:,J)=-62.51
          D(:,J)=67.86
          VA(:,J)=-83.33
          VAF(:,J)=511.25
          DVMWAD(:,J)=.37
        ENDDO
!$OMP END SINGLE
'''
    fast='''SMPI_SIZE.GT.1.AND.SMPI_LOCALBC.AND.SMPI_OUTER.AND.
     &     NUMEPC.EQ.0.AND.NUMQBC.EQ.0.AND.NUMDBC.EQ.0'''
    call=old if role=='old' else once(old,'        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN','        IF('+fast+') THEN\n          CALL SECOM_WAD_FIVE_TEAM(EL,D)\n        ELSEIF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN')
    guard=''
    if role!='old':
        guard='''
!$OMP SINGLE
        IF('''+fast+''') THEN
! North-only messages MUST NOT repair south ghosts by accident, even
! for a single-row owner. The following real CHECK and next EXTRNL
! have to agree with serial with those stale values left in place.
          IF(SMPI_JLO.GT.1) THEN
            J=SMPI_JLO-1
            IF(ANY(VA(:,J).NE.-83.33)) STOP 81
            IF(ANY(VAF(:,J).NE.511.25)) STOP 82
            IF(ANY(DVMWAD(:,J).NE..37)) STOP 83
          ENDIF
          IF(NWC.NE.'''+('4' if role=='bad-vaf-north' else '5')+''') STOP 84
        ENDIF
!$OMP END SINGLE'''
    s=once(s,old,before+call+guard)
    s=once(s,'        CALL SECOM_CHECK_TEAM(VMAX,EMIN,IIT,JIT,RV,RE)',
             '        CALL SECOM_CHECK_READY_TEAM(VMAX,EMIN,IIT,JIT,RV,RE,SMPI_LOCALBC)')
    for n,x in [('driver.f',s),('secom_mpi.f90',rt),('wandd_kky.f',w)]:(d/n).write_bytes(x.encode('latin1'))
    p=d/'inputs.json';data=json.loads(p.read_text());data['role']=role
    data['files']={n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in data['files']};p.write_text(json.dumps(data,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['stage','build']);p.add_argument('directory',type=Path);p.add_argument('--source',type=Path,default=HERE.parents[1]/'source');p.add_argument('--role',choices=['old','five','bad-el-south','bad-vaf-north'],default='old');p.add_argument('--fc',default='mpifort');p.add_argument('--smoke',action='store_true');a=p.parse_args()
    if a.action=='stage':stage(a.directory,a.source,a.role)
    else:build(a.directory,a.fc,a.smoke)
