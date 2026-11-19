"""Veto unnecessary post-corner UA/VA refresh once at geometry initialization."""
import argparse
import hashlib
from pathlib import Path
import re
from entry import HERE,source,active,once

def build(src):
    b=source(src,'mpi_boundary.f');e=source(src,'extrnl.f');rt=source(src,'secom_mpi.f90');m=source(src,'mpideck')
    if 'SMPI_KEEP_UA' in m:raise ValueError('already has corner metadata')
    # Must see the actual full consumer/producer sequence, not a same-name API.
    if 'CALL SMPI_ENTRY_EPOCH(' not in e:raise ValueError('epoch entry parent required')
    x=e.split('      SUBROUTINE SECOM_ADV34_CORNERS()')[1].split('      RETURN')[0]
    lhs=re.findall(r'(?mi)^\s*(ua|va)\(ichard\(n\)(?:\+1)?,jchard\(n\)(\+1)?\)\s*=',x)
    assert [(n.upper(),k) for n,k in lhs]==[('VA','+1'),('UA',''),('VA','+1'),('UA',''),('VA',''),('UA',''),('VA',''),('UA','')],lhs
    assert 'IF(SMPI_LOCALBC.AND..NOT.SMPI_OWNS(JCHARD(N))) CYCLE' in active(x)
    early=e[:e.index('      IF(EDGE_LATE) CALL SECOM_ADV34_CORNERS()')];ad=source(src,'advave.f');bc=source(src,'bcond.f');bc12=bc[bc.index(' 10   CONTINUE'):bc.index('      IF(IDX.EQ.12) RETURN')]
    for s in [early,ad,bc12]:assert not re.search(r'(?m)^(?:\d+\s*)?(UA|VA)\([^=]*?\)\s*=',active(s)), 'unexpected pre-corner velocity writer'
    mom=active(e[e.index('      DO 425 J='):e.index('      CALL BCOND(13,')]);rr={v:set(re.findall(r'\b'+v+r'\(([^()]+)\)',mom)) for v in ('UA','VA')}
    assert rr['UA']=={'I+1,J','I,J','I+1,J-1','I,J-1'} and rr['VA']=={'I,J+1','I,J','I-1,J+1','I-1,J'},rr
    m=once(m,'      SAVE /SECOM_MPI_BC_FLAG/\n','''      SAVE /SECOM_MPI_BC_FLAG/
      LOGICAL SMPI_KEEP_UA,SMPI_KEEP_VA
      COMMON /SECOM_MPI_CORNER_KEEP/ SMPI_KEEP_UA,SMPI_KEEP_VA
      SAVE /SECOM_MPI_CORNER_KEEP/
''')
    b=once(b,'      SMPI_LOCALBC=.FALSE.\n','''      INTEGER IKEEP_U,IKEEP_V,JTRG_V
      SMPI_KEEP_UA=.TRUE.
      SMPI_KEEP_VA=.TRUE.
      IKEEP_U=0
      IKEEP_V=0
      SMPI_LOCALBC=.FALSE.
''')
    b=once(b,'      END DO\n      IBADBC=0','''! Middle momentum reads UA only to its south (export rank's last row)
! and VA only to its north (export rank's first row). Check geometric
! owners, not our own bounds, so any rank's list can veto the elision.
        JTRG_V=JCORB
        IF(ITYPBC.LE.2) JTRG_V=JCORB+1
        IF(JCORB.LT.JM) THEN
          IF(SMPI_OWNER_RAW(JCORB).NE.SMPI_OWNER_RAW(JCORB+1))
     &      IKEEP_U=1
        ENDIF
        IF(JTRG_V.GT.1) THEN
          IF(SMPI_OWNER_RAW(JTRG_V).NE.SMPI_OWNER_RAW(JTRG_V-1))
     &      IKEEP_V=1
        ENDIF
      END DO
      IBADBC=0''')
    b=once(b,'      SMPI_LOCALBC=IBADBC.EQ.0\n','''      SMPI_LOCALBC=IBADBC.EQ.0
      IF(SMPI_LOCALBC) THEN
        CALL SMPI_INDEX_RAW(IKEEP_U)
        CALL SMPI_INDEX_RAW(IKEEP_V)
        SMPI_KEEP_UA=IKEEP_U.NE.0
        SMPI_KEEP_VA=IKEEP_V.NE.0
      ENDIF
      IF(SMPI_RANK.EQ.0) WRITE(*,'(A,L1,A,L1,A,L1)')
     &  'MPI CORNER PLAN LOCALBC=',SMPI_LOCALBC,
     &  ' KEEP_UA=',SMPI_KEEP_UA,' KEEP_VA=',SMPI_KEEP_VA
''')
    e=once(e,'''        CALL SMPI_HALO(CURV42D,IM,JM)
        CALL SMPI_HALO(UA,IM,JM)
        CALL SMPI_HALO(VA,IM,JM)''','''        CALL SMPI_HALO(CURV42D,IM,JM)
! Geometry was checked once. Keep both publications on public, static
! fallback, or uncertified calls; do not change the curve or ELF epochs.
        IF(SMPI_OUTER.AND.EDGE_LATE.AND.SMPI_SIZE.GT.1) THEN
          IF(SMPI_KEEP_UA) CALL SMPI_HALO(UA,IM,JM)
          IF(SMPI_KEEP_VA) CALL SMPI_HALO(VA,IM,JM)
        ELSE
          CALL SMPI_HALO(UA,IM,JM)
          CALL SMPI_HALO(VA,IM,JM)
        ENDIF''')
    rt=once(rt,'  save /secom_mpi_bc_flag/\n','''  save /secom_mpi_bc_flag/
  logical :: smpi_keep_ua,smpi_keep_va
  common /secom_mpi_corner_keep/ smpi_keep_ua,smpi_keep_va
  save /secom_mpi_corner_keep/
''')
    rt=once(rt,'  data localbc /.false./','''  logical::keep_u,keep_v
  common /secom_mpi_corner_keep/ keep_u,keep_v
  data keep_u,keep_v /.true.,.true./
  data localbc /.false./''')
    assert rt.count('  smpi_localbc=.false.\n')==2
    rt=rt.replace('  smpi_localbc=.false.\n','  smpi_localbc=.false.\n  smpi_keep_ua=.true.;smpi_keep_va=.true.\n')
    return {'mpideck':m,'mpi_boundary.f':b,'extrnl.f':e,'secom_mpi.f90':rt}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();files=build(a.source);a.out.mkdir()
    for n,s in files.items():
        b=s.encode('latin1');(a.out/n).write_bytes(b);print(hashlib.sha256(b).hexdigest(),n)
