"""v46: generate v44 EXTRNL with the v45 overlap, without re-reading COR for admission."""
import argparse,difflib,hashlib,json,subprocess
from pathlib import Path

COMMIT='87490b3a06f85cdafdadce7730ff8547a2deb37a'
PIN='8099efb76e98f44b69af1b5a830d87497592ab835cc6bb68c88168119545b7d3'
def sha(s):return hashlib.sha256(s.encode()).hexdigest()
def once(s,a,b):
    if s.count(a)!=1:raise ValueError('v44 anchor changed: '+a)
    return s.replace(a,b)

def generate(extrnl_text,witness=False):
    s=extrnl_text
    if sha(s)!=PIN:raise ValueError('wrong v44 scientific source')
    d='      LOGICAL ROTATE,EDGE_LATE'
    s=once(s,d,d+'''
      LOGICAL MODE46,FAST46,ON46,SECOM46_ENABLED
      SAVE ON46
      DATA ON46 /.FALSE./''')
    a='!$OMP END DO\n!$OMP SINGLE\n      DO 406'
    witness_code=('''!$OMP MASTER
      WRITE(*,*) 'V46_STEP',INTX,IEXT,MERGE(1,0,FAST46)
!$OMP END MASTER
''' if witness else '')
    s=once(s,a,'''!$OMP END DO
      MODE46=ROTATE.AND.TOR.EQ.'BAROTROPIC'.AND.BCTYPE.EQ.'CLAMPED'
      FAST46=.FALSE.
      IF(MODE46.AND.IEXT.GT.1) FAST46=ON46
'''+witness_code+'''!$OMP SINGLE
      DO 406''')
    a='      CALL BCOND(12,DTI2,0)\n!$OMP END SINGLE'
    s=once(s,a,'''      CALL BCOND(12,DTI2,0)
      IF(MODE46.AND.IEXT.EQ.1) ON46=SECOM46_ENABLED()
!$OMP END SINGLE NOWAIT
      IF(FAST46) THEN
!$OMP DO SCHEDULE(STATIC)
        DO J=1,JM
          DO I=1,IM
            CURV42D(I,J)=.25*COR(I,J)
          ENDDO
        ENDDO
!$OMP END DO
      ELSE
!$OMP BARRIER
      ENDIF''')
    s=once(s,'!$OMP DO SCHEDULE(STATIC)\n      DO 95','      IF(.NOT.FAST46) THEN\n!$OMP DO SCHEDULE(STATIC)\n      DO 95')
    s=once(s,"!$OMP END DO\n      IF(TOR.NE.'BAROTROPIC')","!$OMP END DO\n      ENDIF\n      IF(TOR.NE.'BAROTROPIC')")
    return s+'''
! ON46 is read only at IEXT>1, when unchanged DO95 was actually completed
! by the same persistent team on substep 1. No extra COR loads or FP work.
      LOGICAL FUNCTION SECOM46_ENABLED()
      CHARACTER*1 OPT
      INTEGER L,STAT
      CALL GET_ENVIRONMENT_VARIABLE('SECOM_V46_OVERLAP',OPT,LENGTH=L,STATUS=STAT)
      SECOM46_ENABLED=STAT.EQ.1.OR.(STAT.EQ.0.AND.L.EQ.1.AND.OPT.EQ.'1')
      END
'''

def main():
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--repo',type=Path);g.add_argument('--source',type=Path)
    ap.add_argument('--out',required=True,type=Path);ap.add_argument('--witness',action='store_true')
    a=ap.parse_args()
    if a.source:
        old=(a.source/'extrnl.f' if a.source.is_dir() else a.source).read_bytes().decode()
    else:old=subprocess.check_output(['git','-C',str(a.repo),'show',COMMIT+':source/extrnl.f']).decode()
    new=generate(old,a.witness)
    patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),'v44/extrnl.f','v46/extrnl.f'))
    info=dict(parent=COMMIT,parent_sha=PIN,source_sha=sha(new),science_sha=sha(generate(old)),witness=a.witness,
              generator_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),patch_sha=sha(patch))
    a.out.mkdir(exist_ok=False)
    for n,s in [('extrnl.f',new),('v46.patch',patch),('manifest.json',json.dumps(info,indent=2)+'\n')]:
        with (a.out/n).open('x') as f:f.write(s)
    print(json.dumps(info))
if __name__=='__main__':main()
