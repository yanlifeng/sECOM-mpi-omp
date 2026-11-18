"""Probe the actual wet/dry and safety-check bodies against the unmodified serial sources."""
from pathlib import Path
import os,sys,subprocess,hashlib,json
root=Path(__file__).resolve().parent
old=Path('/home/user_home/ylf/secom-parallel-check-20261118.jNqsOL/code-delivery/source')
new=root/sys.argv[1]
test=root/('qa-kernels-'+sys.argv[1]);test.mkdir()
def build(src,omp):
    d=test/('omp' if omp else 'serial'); d.mkdir()
    prog=(src/'ecom3d.f').read_text()
    a=prog.index('      VAMAX=-1.E10')
    b=prog.index('!end,   check elevation minimum',a)
    # No reimplementation of either loop: substitute only the outer driver's
    # abort label with a return path; report the exact production indices.
    (d/'checks.f').write_text("""      SUBROUTINE CHECKTEST(VAMAX,ELMINX,IAMAX,JAMAX)
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
      INTEGER ROWIMAX(JM),ROWIMIN(JM)
"""+prog[a:b]+"\n 9100 CONTINUE\n      END\n")
    (d/'driver.f').write_text("""      PROGRAM QATEST
      INCLUDE 'comdeck'
      WETMIN=.30
      WETEPS=.03
      OPEN(99,FILE='probe.bin',ACCESS='STREAM',FORM='UNFORMATTED',
     & STATUS='NEW')
      DO M=1,6
        DO J=1,JM
          DO I=1,IM
            D(I,J)=WETMIN+(MOD(7*I+3*J+M,7)-3)*.01
            IF(M.EQ.2) D(I,J)=WETMIN
            IF(M.GE.4) D(I,J)=2.
            EL(I,J)=(MOD(I+J,5)-2)*.023
            ELF(I,J)=EL(I,J)
            ETB(I,J)=I*.01-J*.003
            ET(I,J)=-7.
            ETF(I,J)=-9.
            FSM(I,J)=1.
            IF(MOD(I+2*J,19).EQ.0) FSM(I,J)=0.
            DUM(I,J)=1.
            DVM(I,J)=1.
            DUMWAD(I,J)=1.
            DVMWAD(I,J)=1.
            IF(MOD(I+J,7).EQ.0) DUMWAD(I,J)=0.
            IF(MOD(2*I+J,7).EQ.0) DVMWAD(I,J)=0.
            UA(I,J)=(MOD(I+J,3)-1)*.082
            VA(I,J)=(MOD(3*I+J,5)-2)*.057
            UAF(I,J)=UA(I,J)*.8
            VAF(I,J)=VA(I,J)*.9
            IF(M.GE.4) THEN
              UA(I,J)=0.
              VA(I,J)=0.
              UAF(I,J)=0.
              VAF(I,J)=0.
              ELF(I,J)=-1.
              IF(M.EQ.5) THEN
                UAF(I,J)=101.
                DUMWAD(I,J)=1.
              ENDIF
              IF(M.EQ.6) FSM(I,J)=0.
            ENDIF
          ENDDO
        ENDDO
        CALL WANDD_KKY(EL,D)
        VAMAX=-99.
        ELMINX=123.
        IAMAX=-99
        JAMAX=-99
        CALL CHECKTEST(VAMAX,ELMINX,IAMAX,JAMAX)
        IF(M.EQ.4.AND.(IAMAX.NE.1.OR.JAMAX.NE.1)) STOP 21
        IF(M.EQ.5.AND.(VAMAX.NE.101..OR.IAMAX.NE.IMM1.OR.
     &   JAMAX.NE.JMM1.OR.ELMINX.NE.123.)) STOP 22
        IF(M.EQ.6.AND.(IAMAX.NE.IMM1.OR.JAMAX.NE.JMM1.OR.
     &   ELMINX.NE.1.E10)) STOP 23
        WRITE(99) UA,VA,UAF,VAF,DUMWAD,DVMWAD,ET,ETF,EL,D
        WRITE(99) VAMAX,ELMINX,IAMAX,JAMAX
      ENDDO
      CLOSE(99)
      END
""")
    subprocess.run(['gfortran','-O3','-fno-fast-math','-ffp-contract=off','-fcheck=bounds','-ffpe-trap=invalid,zero,overflow','-std=legacy','-ffixed-line-length-132',*(['-fopenmp'] if omp else []),'-I'+str(src),'driver.f','checks.f',str(src/'wandd_kky.f'),'-o','qa'],cwd=d,check=True)
    return d/'qa'
outs=[]
for src,omp,widths in ((old,False,[1]),(new,True,[1,3,32])):
    binary=build(src,omp)
    for n in widths:
        d=binary.parent/str(n); d.mkdir()
        subprocess.run(['taskset','-c',f'40-{39+n}',str(binary)],cwd=d,check=True,env=dict(os.environ,OMP_DYNAMIC='FALSE',OMP_NUM_THREADS=str(n),OMP_PROC_BIND='TRUE'))
        path=d/'probe.bin'
        h=hashlib.sha256()
        with path.open('rb') as f:
            for block in iter(lambda:f.read(1048576),b''):h.update(block)
        outs.append(dict(src=str(src),width=n,sha256=h.hexdigest(),size=path.stat().st_size))
        assert outs[0]['sha256']==h.hexdigest()
(test/'results.json').write_text(json.dumps(outs,indent=2))
