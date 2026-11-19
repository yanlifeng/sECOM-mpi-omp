"""Stage tiny real 8-substep EXTRNL/WAD/CHECK kernels, signed-NaN dead-N ghosts."""
from pathlib import Path
import hashlib
import json
from entry import once
from test_entry import stage as base_stage
from momentum_south import build as south


def stage(dst, src, role):
    base_stage(dst, src, 'old')
    e = (dst/'extrnl.f').read_text() if role == 'old' else south(src)
    d = (dst/'driver.f').read_text()
    # The real WAD5 + ready CHECK path and both public fallbacks.
    old = '''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL SECOM_WAD_COLOR_TEAM(EL,D)'''
    d = once(d, old, '''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC.AND.SMPI_OUTER) THEN
          CALL SECOM_WAD_FIVE_TEAM(EL,D)
        ELSE IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL SECOM_WAD_COLOR_TEAM(EL,D)''')
    d = once(d, 'CALL SECOM_CHECK_TEAM(VMAX,EMIN,IIT,JIT,RV,RE)',
             'CALL SECOM_CHECK_READY_TEAM(VMAX,EMIN,IIT,JIT,RV,RE,SMPI_LOCALBC)')
    d = once(d, 'SMPI_OUTER=LAYOUT.EQ.1', 'SMPI_OUTER=LAYOUT.EQ.1.OR.LAYOUT.EQ.4')
    d = once(d, 'CALL EXTRNL_TEAM(AX,AY,RX,RY,.03,.06,.TRUE.)',
             'CALL EXTRNL_TEAM(AX,AY,RX,RY,.03,.06,LAYOUT.NE.4)')
    marker=' 310  CONTINUE\n      IF(EDGE_LATE) CALL SECOM_ADV34_CORNERS()\n!$OMP END SINGLE\n'
    p='''      IF(SMPI_LOCALBC.AND.SMPI_OUTER.AND.ROTATE) THEN
!$OMP SINGLE
        IF(SMPI_JHI.LT.JM) THEN
          CURV42D(:,SMPI_JHI+1)=TRANSFER(INT(Z'7F800011'),0.)
          ELF(:,SMPI_JHI+1)=TRANSFER(INT(Z'7F800017'),0.)
        ENDIF
        IF(SMPI_JLO.GT.1) THEN
          CURV42D(:,SMPI_JLO-1)=-19.6875
          ELF(:,SMPI_JLO-1)=2.375
        ENDIF
!$OMP END SINGLE
      ENDIF
'''
    e=once(e,marker,marker+p)
    if role=='south':
        # Make accidental bidirectional repair visible without touching sNaNs as FP.
        marker=' 5000 CONTINUE'
        e=once(e,marker,'''      IF(SMPI_LOCALBC.AND.SMPI_OUTER.AND.ROTATE.AND.SMPI_JHI.LT.JM) THEN
!$OMP SINGLE
        DO I=1,IM
          IF(TRANSFER(CURV42D(I,SMPI_JHI+1),0).NE.INT(Z'7F800011')) STOP 88
          IF(TRANSFER(ELF(I,SMPI_JHI+1),0).NE.INT(Z'7F800017')) STOP 89
        ENDDO
!$OMP END SINGLE
      ENDIF
''' + marker)
    if role.startswith('wrong-'):
        field={'wrong-curve':'CURV42D','wrong-elf':'ELF'}[role]
        anchor='CALL SMPI_FLUX_HALO_S('+field+',IM,JM)'
        e=once(e,anchor,anchor.replace('_S(','_N('))
    for n,t in [('extrnl.f',e),('driver.f',d)]:
        (dst/n).write_bytes(t.encode('latin1'))
    p=dst/'inputs.json'
    cfg=json.loads(p.read_text());cfg['role']=role
    cfg['files']={n:hashlib.sha256((dst/n).read_bytes()).hexdigest() for n in cfg['files']}
    p.write_text(json.dumps(cfg,sort_keys=True,indent=2))
