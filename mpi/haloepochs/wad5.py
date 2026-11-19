"""Five directional WAD Sendrecvs only inside the validated ENTRY→WAD→READY chain."""
import argparse
import hashlib
from pathlib import Path
import re
from entry import HERE,source,active,once

COND='''SMPI_SIZE.GT.1.AND.SMPI_LOCALBC.AND.SMPI_OUTER.AND.
     &     NUMEPC.EQ.0.AND.NUMQBC.EQ.0.AND.NUMDBC.EQ.0'''

def build(src):
    ext=source(src,'extrnl.f');wad=source(src,'wandd_kky.f');check=source(src,'baro_check.f');main=source(src,'ecom3d.f');rt=source(src,'secom_mpi.f90')
    if 'CALL SMPI_ENTRY_EPOCH(D,UA,VA,UAB,VAB,EL,ELB,DUMWAD,' not in ext or 'SECOM_CHECK_READY_TEAM' not in main or 'subroutine smpi_flux_one_raw(' not in rt.lower():
        raise ValueError('require complete validated ENTRY + fluxdir + checkready parent')
    old='''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL SECOM_WAD_COLOR_TEAM(EL,D)
'''
    if 'WAD_INPUT_EPOCH' in wad or 'SECOM_WAD_FIVE_TEAM' in wad:raise ValueError('incompatible WAD input/version contract')
    assert 'CALL EXTRNL_TEAM(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,.TRUE.)\n'+old in main
    assert 'ROWVMAX,ROWEMIN,SMPI_LOCALBC)' in main
    jref=set(re.findall(r'\b(?:VAF|DVMWAD)\(([^()]+)\)',active(check)))
    assert jref=={'I,J','I,J+1'},('unreviewed CHECK incoming-neighbor read',jref)
    public,private=wad.split('      RECURSIVE SUBROUTINE SECOM_WAD_COLOR_TEAM(ELEV,TD)')
    # Fail closed on any new vertical source in the WAD loop.
    q=active(private.split('DO ICOLOR=1,2',1)[1].split('CALL SMPI_WAD_PUBLISH',1)[0]);rd=set(re.findall(r'\b(?:VA|VAF|DVMWAD)\(([^()]+)\)',q));assert rd=={'I,J+1'},rd
    private=once(private,'''        CALL SMPI_HALO(ELEV,IM,JM)
        CALL SMPI_HALO(TD,IM,JM)
        CALL SMPI_HALO(VA,IM,JM)
        CALL SMPI_HALO(DVMWAD,IM,JM)''','''! ELEV=EL and TD=D: south sides are kept CURRENT for the next entry
! even though this WAD reads only j+1. Velocities and v-face masks
! for this WAD have no south-halo consumer under the private contract.
        CALL SMPI_WAD_PAIR(ELEV,TD,IM,JM,.TRUE.)
        CALL SMPI_WAD_PAIR(VA,DVMWAD,IM,JM,.FALSE.)''')
    private=once(private,'''        CALL SMPI_HALO(VA,IM,JM)
        CALL SMPI_HALO(DVMWAD,IM,JM)
        CALL SMPI_HALO(VAF,IM,JM)''','''! CHECK sees current north DVMWAD,VAF. Next ENTRY retains north VA.
! No public caller observes the deliberately unrefreshed south ghosts.
        CALL SMPI_WAD_PAIR(VA,DVMWAD,IM,JM,.FALSE.)
        CALL SMPI_FLUX_HALO_N(VAF,IM,JM)''')
    main=once(main,old,'        IF('+COND+') THEN\n          CALL SECOM_WAD_FIVE_TEAM(EL,D)\n        ELSEIF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN\n          CALL SECOM_WAD_COLOR_TEAM(EL,D)\n')
    rt=once(rt,'  if(nx==0)return\n  call smpi_entry_epoch_release()',
               '  if(nx==0)return\n  call smpi_wad_pair_release()\n  call smpi_entry_epoch_release()')
    return {'wandd_kky.f':wad+'\n      RECURSIVE SUBROUTINE SECOM_WAD_FIVE_TEAM(ELEV,TD)'+private,
            'ecom3d.f':main,'secom_mpi.f90':rt+'\n'+(HERE/'wad_pair.f90.inc').read_text()}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();files=build(a.source);a.out.mkdir()
    for n,s in files.items():
        b=s.encode('latin1');(a.out/n).write_bytes(b);print(hashlib.sha256(b).hexdigest(),n)
