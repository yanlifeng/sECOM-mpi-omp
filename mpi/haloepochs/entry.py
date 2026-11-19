"""Generate strictly gated EXTRNL entry epochs OR independent WAD-input reuse."""
import argparse
from pathlib import Path
import hashlib
import re
from prepare import once

HERE=Path(__file__).resolve().parent
PLANES='D UA VA UAB VAB EL ELB DUMWAD DVMWAD'.split()

def source(src,name):return (src/name).read_bytes().decode('latin1')
def active(s):
    return '\n'.join(x.split('!')[0].strip().upper() for x in s.splitlines() if x and x[0] not in '!cC*' and x.split('!')[0].strip())

def entry(src):
    e=source(src,'extrnl.f'); rt=source(src,'secom_mpi.f90');w=source(src,'wandd_kky.f');adv=source(src,'advave.f');mn=source(src,'ecom3d.f');bc=source(src,'mpi_boundary.f')
    if 'SMPI_ENTRY_EPOCH' in e or 'WAD_INPUT_EPOCH' in w:raise ValueError('not stackable with previous entry/WAD-input reduction')
    # Fail closed on touched cross-epoch invariants: neither the active WAD
    # nor intervening checker may start mutating the reused input versions.
    for b in w.split('RECURSIVE SUBROUTINE')[1:]:
        p=active(b)
        first=p.split('DO ICOLOR=1,2')[0]
        for f in ['ELEV','TD','VA','DVMWAD']:
            assert 'CALL SMPI_HALO('+f+',IM,JM)' in first,(f,'WAD entry')
        for f in ['VA','DVMWAD','VAF']:
            assert 'CALL SMPI_HALO('+f+',IM,JM)' in p.rsplit('CALL SMPI_WAD_PUBLISH',1)[1],(f,'WAD return')
        assert not re.search(r'(?m)^(ELEV|TD|EL|D)\([^=]*\)\s*=',p)
    assert 'CALL EXTRNL_TEAM(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,.TRUE.)' in mn
    chain=active(mn.split('      DO ISEQ=1,ISPLIT')[1].split('      ENDDO')[0]);assert 'WANDD_TEAM(EL,D)' in chain and 'CHECK_' in chain
    assert 'IF(NUMQBC.NE.0.OR.NUMDBC.NE.0.OR.NUMEPC.NE.0) GOTO 90' in active(bc)
    assert 'CALL SECOM_ADV34_CORNERS()' in e and 'IF(SMPI_LOCALBC.AND..NOT.SMPI_OWNS(JCHARD(N))) CYCLE' in e
    joined=active(e+'\n'+adv);ref=set(re.findall(r'\bDVMWAD\(([^()]+)\)',joined))
    assert ref=={'I,J','I-1,J'},('DVMWAD nonowned access or unfamiliar expression',ref)
    assert all(c in e for c in ['CALL SMPI_HALO(CURV42D,IM,JM)','CALL SMPI_HALO(UA,IM,JM)','CALL SMPI_HALO(VA,IM,JM)','CALL SMPI_HALO(ELF,IM,JM)'])
    old='      IF(SMPI_LOCALBC) THEN\n'+''.join(f'      CALL SMPI_HALO({f},IM,JM)\n' for f in PLANES)+'      ELSE\n'
    rep='''      IF(SMPI_LOCALBC) THEN
      IF(SMPI_SIZE.GT.1.AND.SMPI_OUTER.AND.ROTATE.AND.
     &   NUMEPC.EQ.0.AND.NUMQBC.EQ.0.AND.NUMDBC.EQ.0) THEN
! The caller's MASTER sets global IEXT without an exit barrier. Do
! not evaluate the FIRST flag before that store is published to all.
!$OMP BARRIER
        CALL SMPI_ENTRY_EPOCH(D,UA,VA,UAB,VAB,EL,ELB,DUMWAD,
     &                        DVMWAD,IM,JM,IEXT.EQ.1)
      ELSE
'''+''.join(f'      CALL SMPI_HALO({f},IM,JM)\n' for f in PLANES)+'      ENDIF\n      ELSE\n'
    rt=once(rt,'  if(nx==0)return\n','  if(nx==0)return\n  call smpi_entry_epoch_release()\n')
    return {'extrnl.f':once(e,old,rep),'secom_mpi.f90':rt+'\n'+(HERE/'entry_epoch.f90.inc').read_text()}

def wadinput(src):
    e=source(src,'extrnl.f');w=source(src,'wandd_kky.f');m=source(src,'ecom3d.f')
    if 'SMPI_ENTRY_EPOCH' in e or 'WAD_INPUT_EPOCH' in w:raise ValueError('refuse stacked/duplicate WAD-input reduction')
    entryblock=e.split('      IF(SMPI_LOCALBC) THEN',1)[1].split('      ELSE',1)[0]
    assert 'CALL SMPI_HALO(DVMWAD,IM,JM)' in entryblock or 'CALL SMPI_HALO9(D,UA,VA,UAB,VAB,EL,ELB,DUMWAD,DVMWAD,IM,JM)' in entryblock
    assert not re.search(r'(?m)^DVMWAD\([^=]*\)\s*=',active(e))
    cut='      RECURSIVE SUBROUTINE SECOM_WAD_COLOR_TEAM(ELEV,TD)'
    if cut in w:
        left,right=w.split(cut);assert m.count('CALL SECOM_WAD_COLOR_TEAM(EL,D)')==1
        right=once(right,'        CALL SMPI_HALO(VA,IM,JM)\n        CALL SMPI_HALO(DVMWAD,IM,JM)\n! A process',
                         '        CALL SMPI_HALO(VA,IM,JM)\n! WAD_INPUT_EPOCH: DVMWAD is fresh from the immediately prior EXTRNL entry.\n! A process')
        return {'wandd_kky.f':left+cut+right}
    cut='      RECURSIVE SUBROUTINE WANDD_TEAM(ELEV,TD)';right=w.split(cut)[1]
    right=once(right,'        CALL SMPI_HALO(DVMWAD,IM,JM)\n      ENDIF\n      DO ICOLOR',
                     '! WAD_INPUT_EPOCH: EXTRNL published this mask; no intervening writes.\n      ENDIF\n      DO ICOLOR')
    priv='      RECURSIVE SUBROUTINE SECOM_WAD_INREADY_TEAM(ELEV,TD)'
    m=once(m,'        CALL WANDD_TEAM(EL,D)\n','''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL SECOM_WAD_INREADY_TEAM(EL,D)
        ELSE
          CALL WANDD_TEAM(EL,D)
        ENDIF
''')
    return {'wandd_kky.f':w+'\n'+priv+right,'ecom3d.f':m}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--variant',choices=['entry','wadinput'],required=True);a=p.parse_args()
    files=(entry if a.variant=='entry' else wadinput)(a.source);a.out.mkdir()
    for f,t in files.items():
        b=t.encode('latin1');(a.out/f).write_bytes(b);print(hashlib.sha256(b).hexdigest(),f)
