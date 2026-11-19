"""Isolated WAD→CHECK halo-epoch unit; stage locally, build/run in an allocation.

Pass the executable two arguments: expected ranks, expected OMP threads. Each
rank checks eight field bit patterns, all final seam ghosts, both CHECK APIs,
the signed extrema + indices, and the actual RAW halo/scatter call counts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess as sp
from prepare import checkready,wadcolor,once

HERE=Path(__file__).resolve().parent;BASE='1b58f33639538c2d214f14a8f9a575ceffed5267'

def stage(d,src,variant):
    root=HERE.parents[1];source=src.resolve()
    wad=(wadcolor(src)['wandd_kky.f'] if variant=='both' else (src/'wandd_kky.f').read_text())
    check=checkready(src)['baro_check.f']
    rt=(src/'secom_mpi.f90').read_text()
    rt=once(rt,'  call dims(im,jm)\n  ! Use scalar scratch',
            '  call qa_count_halo()\n  call dims(im,jm)\n  ! Use scalar scratch')
    rt=once(rt,'  call dims(im,jm)\n  call MPI_Sendrecv(flags,im,MPI_INTEGER',
            '  call qa_count_publish()\n  call dims(im,jm)\n  call MPI_Sendrecv(flags,im,MPI_INTEGER')
    rt+='''
subroutine qa_count_halo()
  implicit none
  integer nh,np
  common /qa_epoch_counts/ nh,np
  nh=nh+1
end subroutine
subroutine qa_count_publish()
  implicit none
  integer nh,np
  common /qa_epoch_counts/ nh,np
  np=np+1
end subroutine
'''
    ref=''
    for f,symbols in [('wandd_kky.f',['WANDD_KKY','WANDD_TEAM']),('baro_check.f',['SECOM_CHECK_TEAM'])]:
        t=sp.check_output(['git','-C',str(root),'show',BASE+':source/'+f]).decode()
        for n in symbols:t=t.replace(n,'V46_'+n)
        ref+=t+'\n'
    drive=(root/'validationmpi/wad_owned_driver.f').read_text()
    drive=once(drive,'      CHARACTER*16 ARG','''      REAL ROWVMAX(JM),ROWEMIN(JM),WANTV,WANTE
      INTEGER IWANT,JWANT,NH,NP
      COMMON /QA_EPOCH_COUNTS/ NH,NP
      CHARACTER*16 ARG''')
    drive=once(drive,'        CALL V46_WANDD_KKY(EL,D)','''        ELF=EL
        CALL V46_WANDD_KKY(EL,D)
        IWANT=-43
        JWANT=17
        WANTE=78.
!$OMP PARALLEL DEFAULT(SHARED)
        CALL V46_SECOM_CHECK_TEAM(WANTV,WANTE,IWANT,JWANT,ROWVMAX,ROWEMIN)
!$OMP END PARALLEL''')
    # Stale ELF rows may never provide the distributed minimum.
    drive=once(drive,'          D(:,J)=-271.','          ELF(:,J)=-1812.\n          D(:,J)=-271.')
    drive=once(drive,'        CALL WANDD_KKY(EL,D)','''        DO J=SMPI_JLO,SMPI_JHI
          ELF(:,J)=EL(:,J)
        ENDDO
        NH=0
        NP=0
        IAMAX=-43
        JAMAX=17
        ELMINX=78.
!$OMP PARALLEL DEFAULT(SHARED)
        CALL '''+('SECOM_WAD_COLOR_TEAM' if variant=='both' else 'WANDD_TEAM')+'''(EL,D)
        CALL SECOM_CHECK_READY_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN,.TRUE.)
!$OMP END PARALLEL
        IF(TRANSFER(VAMAX,0).NE.TRANSFER(WANTV,0).OR.
     &    TRANSFER(ELMINX,0).NE.TRANSFER(WANTE,0).OR.
     &    IAMAX.NE.IWANT.OR.JAMAX.NE.JWANT) THEN
          PRINT *, 'BAD READY',M,SMPI_RANK,VAMAX,WANTV,IAMAX,IWANT,JAMAX,JWANT
          NFAIL=NFAIL+1
        ENDIF
        IF(SMPI_SIZE.GT.1.AND.(NH.NE.'''+('7.OR.NP.NE.1' if variant=='both' else '9.OR.NP.NE.2')+''')) THEN
          PRINT *, 'BAD COUNTS PRIVATE',M,SMPI_RANK,NH,NP
          NFAIL=NFAIL+1
        ENDIF''')
    drive=once(drive,'        ENDDO\n      ENDDO\n      CALL MPI_ALLREDUCE','''        ENDDO
! Deliberately invalidate the ghosts before the separate public API.
! The public checker must communicate even if the previous private one did not.
        DO J=1,JM
          IF(J.GE.SMPI_JLO.AND.J.LE.SMPI_JHI) CYCLE
          DVMWAD(:,J)=3.
          VAF(:,J)=499.
        ENDDO
        NH=0
        NP=0
        IAMAX=-43
        JAMAX=17
        ELMINX=78.
!$OMP PARALLEL DEFAULT(SHARED)
        CALL SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
!$OMP END PARALLEL
        IF(TRANSFER(VAMAX,0).NE.TRANSFER(WANTV,0).OR.
     &    TRANSFER(ELMINX,0).NE.TRANSFER(WANTE,0).OR.
     &    IAMAX.NE.IWANT.OR.JAMAX.NE.JWANT) THEN
          PRINT *, 'BAD PUBLIC',M,SMPI_RANK
          NFAIL=NFAIL+1
        ENDIF
        IF(SMPI_SIZE.GT.1.AND.(NH.NE.2.OR.NP.NE.0)) THEN
          PRINT *, 'BAD COUNTS PUBLIC',M,SMPI_RANK,NH,NP
          NFAIL=NFAIL+1
        ENDIF
      ENDDO
      CALL MPI_ALLREDUCE''')
    # Exercise speed-before-minimum abort, equal maxima and signed-zero minima
    # with a trained non-aliased VA/VAF fixture in addition to dry/face masks.
    drive=once(drive,'          ELSEIF(M.EQ.6) THEN','''            IF(I.EQ.2.OR.I.EQ.3) UAF(I,J)=144.
            IF(I.EQ.6.AND.J.EQ.4) EL(I,J)=-0.
          ELSEIF(M.EQ.6) THEN''')
    d.mkdir(); shas={}
    fs={'wandd_kky.f':wad,'baro_check.f':check,'runtime.f90':rt,'v46.f':ref,'driver.f':drive,
        'comdeck':once((src/'comdeck').read_text(),'IM=372,JM=479','IM=12,JM=7'),
        'mpideck':(src/'mpideck').read_text()}
    for f,s in fs.items():
        b=s.encode();(d/f).write_bytes(b);shas[f]=hashlib.sha256(b).hexdigest()
    (d/'inputs.json').write_text(json.dumps(dict(variant=variant,v46=BASE,files=shas),indent=2)+'\n')

def build(d,fc,mutant,smoke):
    manifest=json.loads((d/'inputs.json').read_text());files=manifest['files']
    assert all(hashlib.sha256((d/n).read_bytes()).hexdigest()==v for n,v in files.items())
    flags=['-O3','-fopenmp','-ffp-contract=off','-fno-fast-math','-ffpe-trap=invalid,zero,overflow','-fcheck=bounds']
    if smoke:
        # Exactly one process; check actual compiler code but make no MPI claim.
        m='''module mpi
          integer,parameter::MPI_INTEGER=0,MPI_SUM=0,MPI_COMM_WORLD=0
          contains
          subroutine mpi_allreduce(a,b,n,k,op,c,i)
          integer a,b,n,k,op,c,i
          b=a;i=0
          end subroutine
          end module
'''
        (d/'serial_mpi_stub.f90').write_text(m);sp.run([fc,'-c','serial_mpi_stub.f90'],cwd=d,check=True)
    sp.run([fc,*flags,'-cpp',*([] if smoke else ['-DSECOM_USE_MPI']),'-I.','-J.','-c','runtime.f90','-o','runtime.o'],cwd=d,check=True)
    wad='wandd_kky.f';check='baro_check.f'
    if mutant=='public-stale':
        f=(d/check).read_text();check='check-public-stale.f';(d/check).write_text(once(f,'ROWVMAX,ROWEMIN,.FALSE.)','ROWVMAX,ROWEMIN,.TRUE.)'))
    elif mutant in ('drop-reset','drop-vaf'):
        assert manifest['variant']=='both'
        prefix='      RECURSIVE SUBROUTINE SECOM_WAD_COLOR_TEAM';left,right=(d/wad).read_text().split(prefix)
        if mutant=='drop-reset':right=once(right,'WAD_SEND(I)=3','WAD_SEND(I)=1')
        else:right=once(right,'        CALL SMPI_HALO(VAF,IM,JM)','! omitted on purpose')
        wad='wad-'+mutant+'.f';(d/wad).write_text(left+prefix+right)
    sp.run([fc,*flags,'-std=legacy','-ffixed-line-length-132','-I.',wad,check,'v46.f','driver.f','runtime.o',
            *(['serial_mpi_stub.o'] if smoke else []),'-o','qa-'+mutant],cwd=d,check=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['stage','build']);p.add_argument('directory',type=Path)
    p.add_argument('--source',default=HERE.parents[1]/'source',type=Path);p.add_argument('--variant',choices=['both','checkready'],default='both')
    p.add_argument('--fc',default='mpifort');p.add_argument('--mutant',choices=['none','public-stale','drop-reset','drop-vaf'],default='none');p.add_argument('--smoke',action='store_true');a=p.parse_args()
    if a.action=='stage':stage(a.directory,a.source,a.variant)
    else:build(a.directory,a.fc,a.mutant,a.smoke)
