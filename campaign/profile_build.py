"""Build a separate wall-time-instrumented main; never alter official timing binaries."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parent
kind=sys.argv[1]
assert kind in ('v2','v3','v4')
prev=Path('/home/user_home/ylf/secom-omp2-20261118.ZtxsGz')
src=(prev if kind=='v2' else root)/'source'
if kind=='v4': src=root/'src4'
unit=(src/'ecom3d.f').read_text()
def before(anchor,lines):
    global unit
    if not anchor.startswith('\n'): anchor='\n'+anchor
    assert unit.count(anchor)==1,anchor
    unit=unit.replace(anchor,'\n'+lines+anchor)
unit=unit.replace("      INCLUDE 'comdeck'", """      INCLUDE 'comdeck'
      REAL*8 mpclock(7),mpt0,omp_get_wtime
      EXTERNAL omp_get_wtime
      DATA mpclock /7*0.d0/""",1)
start="      mpt0=omp_get_wtime()"
def end(i):
    return f"      mpclock({i})=mpclock({i})+(omp_get_wtime()-mpt0)"
before("      CALL BCOND(7,DTI2,0)\n!",start)
before("!NG04152008\n      IF (WAVEDYN",end(1))
before("\n      CALL EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2)",start)
before("      IF(TOR.EQ.'BAROTROPIC') GO TO 440",end(2))
before("!         APPLY FILTER TO REMOVE TIME SPLIT",start)
before("      CALL WANDD_KKY(EL,D)",end(3)+'\n'+start)
before("!  TEST FOR MODEL BLOWUP",end(4)+'\n'+start)
before("!end,   check elevation minimum",end(5))
unit=unit.replace("        CALL ARCHIVE\n        GO TO 9900",start+"\n        CALL ARCHIVE\n"+end(6)+"\n        GO TO 9900",1)
before("                   DO 9000 INTX=ISTART,IEND","      mpclock(7)=omp_get_wtime()")
before(" 9100 CONTINUE\n","      mpclock(7)=omp_get_wtime()-mpclock(7)\n      WRITE(*,*) 'DIAG_SECONDS',mpclock")
out=root/('diag-'+kind+'.f')
out.write_text(unit)
deps=Path('/home/user_home/ylf/kt5-licom-secom-ylf/deps')
flags=['-O3','-fno-fast-math','-ffp-contract=off','-std=legacy','-ffixed-line-length-132','-g','-fbacktrace','-ffpe-trap=invalid,zero,overflow']
obj=root/('diag-'+kind+'.o')
subprocess.run(['gfortran',*flags,*(['-fopenmp'] if kind!='v2' else []),'-I'+str(src),'-I'+str(deps/'include'),'-c',str(out),'-o',str(obj)],check=True)
names=next(s for s in (src/'Makefile.portable').read_text().splitlines() if s.startswith('OBJECT_NAMES =')).split('=',1)[1].split()
libs=subprocess.check_output([str(deps/'bin/nf-config'),'--flibs'],universal_newlines=True).split()
subprocess.run(['gfortran',*flags,'-fopenmp','-o',str(root/('profile-'+kind)),*[str(obj if n=='ecom3d.o' else src/'build/omp'/n) for n in names],*libs],check=True)
