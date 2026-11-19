"""Make one isolated halo-epoch overlay. No in-place writes or MPI execution."""
import argparse
import hashlib
from pathlib import Path

def once(src,a,b):
    if src.count(a)!=1:raise ValueError('expected exactly one source anchor: '+repr(a[:85]))
    return src.replace(a,b)

def checkready(src):
    name='SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)'
    orig='      RECURSIVE SUBROUTINE '+name+'\n'
    new=orig+'''      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
! Independent callers may supply poisoned/stale ghosts: always refresh.
      CALL SECOM_CHECK_READY_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,
     &                            ROWVMAX,ROWEMIN,.FALSE.)
      END
      RECURSIVE SUBROUTINE SECOM_CHECK_READY_TEAM(VAMAX,ELMINX,
     &             IAMAX,JAMAX,ROWVMAX,ROWEMIN,HALO_READY)
      LOGICAL HALO_READY
'''
    bar=once((src/'baro_check.f').read_text(),orig,new)
    bar=once(bar,'      IF(SMPI_SIZE.GT.1) THEN\n        CALL SMPI_HALO(DVMWAD,IM,JM)',
             '      IF(SMPI_SIZE.GT.1.AND..NOT.HALO_READY) THEN\n        CALL SMPI_HALO(DVMWAD,IM,JM)')
    m=(src/'ecom3d.f').read_bytes().decode('latin1')
    m=once(m,'        CALL WANDD_TEAM(EL,D)\n        CALL '+name,
           '''        CALL WANDD_TEAM(EL,D)
! WANDD completed its last DVMWAD and VAF halos, incl. END SINGLE
! barriers, immediately above. Nothing modifies either version here.
        CALL SECOM_CHECK_READY_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,
     &                              ROWVMAX,ROWEMIN,SMPI_LOCALBC)''')
    return {'baro_check.f':bar,'ecom3d.f':m}

def wadcolor(src):
    s=(src/'wandd_kky.f').read_text()
    s=once(s,'        CALL SMPI_HALO(DVMWAD,IM,JM)\n      ENDIF\n      DO ICOLOR=1,2',
        '''        CALL SMPI_HALO(DVMWAD,IM,JM)
! A process has just one north seam: J=SMPI_JHI executes in ONE
! color only. Keep its explicit flags until both colors are joined.
!$OMP SINGLE
        WAD_SEND=0
!$OMP END SINGLE
      ENDIF
      DO ICOLOR=1,2''')
    s=once(s,'        JSTOP=SMPI_LAST(JMM1)\n!$OMP SINGLE\n        WAD_SEND=0\n!$OMP END SINGLE',
             '        JSTOP=SMPI_LAST(JMM1)')
    s=once(s,'''!$OMP END DO
      IF(MULTI_WAD) THEN
! Only marked columns change on the remote owner, before the next
! color. Stale unmarked ghost values must never overwrite owner data.
        CALL SMPI_WAD_PUBLISH(WAD_SEND,ET,ETF,VA,VAF,DVMWAD,ETB,IM,JM)
        CALL SMPI_HALO(VA,IM,JM)
        CALL SMPI_HALO(DVMWAD,IM,JM)
        IF(ICOLOR.EQ.2) CALL SMPI_HALO(VAF,IM,JM)
      ENDIF
      ENDDO''','''!$OMP END DO
      ENDDO
      IF(MULTI_WAD) THEN
! Both vertical tests for face (I,J+1) belong to the same (I,J)
! iteration. The opposite color never reads/writes that face. ET/ETF
! may share a target but are only assigned the same target-owner ETB
! and are NOT read by another color (production actuals are EL,D).
! Retain both color DO joins, then commit marked cross-owner writes.
        CALL SMPI_WAD_PUBLISH(WAD_SEND,ET,ETF,VA,VAF,DVMWAD,ETB,IM,JM)
        CALL SMPI_HALO(VA,IM,JM)
        CALL SMPI_HALO(DVMWAD,IM,JM)
        CALL SMPI_HALO(VAF,IM,JM)
      ENDIF''')
    # Keep public WANDD_TEAM/KY semantics and any undocumented aliasing
    # contract untouched. Only the known-private BARO call passes EL,D.
    old=(src/'wandd_kky.f').read_text()
    cut='      RECURSIVE SUBROUTINE WANDD_TEAM(ELEV,TD)\n'
    if s.count(cut)!=1:raise ValueError('expected one WANDD_TEAM definition')
    body=s.split(cut,1)[1]
    private='      RECURSIVE SUBROUTINE SECOM_WAD_COLOR_TEAM(ELEV,TD)\n'+body
    m=(src/'ecom3d.f').read_bytes().decode('latin1')
    m=once(m,'        CALL WANDD_TEAM(EL,D)\n','''        IF(SMPI_SIZE.GT.1.AND.SMPI_LOCALBC) THEN
          CALL SECOM_WAD_COLOR_TEAM(EL,D)
        ELSE
          CALL WANDD_TEAM(EL,D)
        ENDIF
''')
    return {'wandd_kky.f':old+'\n'+private,'ecom3d.f':m}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--variant',choices=['checkready','wadcolor','both'],required=True);a=p.parse_args()
    files=checkready(a.source) if a.variant!='wadcolor' else wadcolor(a.source)
    if a.variant=='both':
        w=wadcolor(a.source);orig='        CALL WANDD_TEAM(EL,D)\n';full=(a.source/'ecom3d.f').read_bytes().decode('latin1')
        before,after=full.split(orig);middle=w['ecom3d.f'][len(before):len(w['ecom3d.f'])-len(after)]
        files['ecom3d.f']=once(files['ecom3d.f'],orig,middle);files['wandd_kky.f']=w['wandd_kky.f']
    a.out.mkdir()
    for name,text in files.items():
        b=text.encode('latin1');(a.out/name).write_bytes(b);print(hashlib.sha256(b).hexdigest(),name)
if __name__=='__main__':main()
