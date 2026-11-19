"""Generate just extrnl.f for v45. No build, solver, checkout or submission."""
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
    if sha(s)!=PIN:raise ValueError('expected exact v44 extrnl.f')
    decl='      LOGICAL ROTATE,EDGE_LATE'
    s=once(s,decl,decl+'''
      LOGICAL FAST45,MODE45,OK45,SECOM45_COR_OK
      SAVE OK45
      DATA OK45 /.FALSE./''')
    # After the flux/ELF barriers IEXT is published to every lane. Recursive
    # automatic locals FAST45/MODE45 are private to each caller; OK45 is saved.
    a='!$OMP END DO\n!$OMP SINGLE\n      DO 406'
    proof=('''!$OMP MASTER
      WRITE(*,*) 'V45_STEP',INTX,IEXT,MERGE(1,0,FAST45)
!$OMP END MASTER
''' if witness else '')
    s=once(s,a,'''!$OMP END DO
      MODE45=ROTATE.AND.TOR.EQ.'BAROTROPIC'.AND.BCTYPE.EQ.'CLAMPED'
      FAST45=.FALSE.
      IF(MODE45.AND.IEXT.GT.1) FAST45=OK45
'''+proof+'''!$OMP SINGLE
      DO 406''')
    a='      CALL BCOND(12,DTI2,0)\n!$OMP END SINGLE'
    s=once(s,a,'''      CALL BCOND(12,DTI2,0)
      IF(MODE45.AND.IEXT.EQ.1) OK45=SECOM45_COR_OK(COR,IM,JM)
!$OMP END SINGLE NOWAIT
      IF(FAST45) THEN
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
    s=once(s,'!$OMP DO SCHEDULE(STATIC)\n      DO 95','      IF(.NOT.FAST45) THEN\n!$OMP DO SCHEDULE(STATIC)\n      DO 95')
    s=once(s,"!$OMP END DO\n      IF(TOR.NE.'BAROTROPIC')","!$OMP END DO\n      ENDIF\n      IF(TOR.NE.'BAROTROPIC')")
    return s+'''
! No floating operation touches COR in this guard (including signaling NaNs).
! Check once per BARO_SUBCYCLE at IEXT=1, not across restarts/outer steps.
      LOGICAL FUNCTION SECOM45_COR_OK(A,N,M)
      USE ISO_FORTRAN_ENV, ONLY: INT32,REAL32
      INTEGER N,M,I,J,E,L,STAT
      INTEGER(INT32) W
      REAL A(N,M)
      CHARACTER*1 OPT
      SECOM45_COR_OK=.FALSE.
      IF(KIND(1.0).NE.REAL32.OR.STORAGE_SIZE(1.0).NE.32) RETURN
      IF(RADIX(1.0).NE.2.OR.DIGITS(1.0).NE.24.OR.MAXEXPONENT(1.0).NE.128.OR.MINEXPONENT(1.0).NE.-125) RETURN
      IF(TRANSFER(1.0,0_INT32).NE.INT(Z'3F800000',INT32)) RETURN
      CALL GET_ENVIRONMENT_VARIABLE('SECOM_V45_OVERLAP',OPT,LENGTH=L,STATUS=STAT)
      IF(STAT.NE.1) THEN
        IF(STAT.NE.0.OR.L.NE.1.OR.OPT.NE.'1') RETURN
      ENDIF
      DO J=1,M
        DO I=1,N
          W=TRANSFER(A(I,J),0_INT32)
          E=IBITS(W,23,8)
          IF(IBCLR(W,31).NE.0_INT32.AND.(E.LT.3.OR.E.EQ.255)) RETURN
        ENDDO
      ENDDO
      SECOM45_COR_OK=.TRUE.
      END
'''

def main():
    ap=argparse.ArgumentParser(description=__doc__);g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--repo',type=Path);g.add_argument('--source',type=Path)
    ap.add_argument('--out',required=True,type=Path);ap.add_argument('--witness',action='store_true')
    a=ap.parse_args()
    if a.source:
        f=a.source/'extrnl.f' if a.source.is_dir() else a.source
        raw=f.read_bytes()
    else:raw=subprocess.check_output(['git','-C',str(a.repo),'show',COMMIT+':source/extrnl.f'])
    old=raw.decode();new=generate(old,a.witness)
    diff=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),'v44/extrnl.f','v45/extrnl.f'))
    info=dict(parent=COMMIT,parent_sha=PIN,source_sha=sha(new),science_sha=sha(generate(old)),witness=a.witness,
              generator_sha=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),patch_sha=sha(diff))
    a.out.mkdir(exist_ok=False)
    for n,s in [('extrnl.f',new),('v45.patch',diff),('manifest.json',json.dumps(info,indent=2)+'\n')]:
        with (a.out/n).open('x') as f:f.write(s)
    print(json.dumps(info))

if __name__=='__main__':main()
