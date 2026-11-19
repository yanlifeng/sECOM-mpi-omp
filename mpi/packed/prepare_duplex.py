"""Experimental two-direction concurrent Halo9; generate a fresh overlay only."""
import argparse
import hashlib
from pathlib import Path
from prepare import build, replace_once


def duplex(src):
    files = build(src)
    rt = files['secom_mpi.f90']
    rt = replace_once(rt, '  real, allocatable, save :: send9(:,:),recv9(:,:)\n',
                     '  real, allocatable, save :: send9(:,:),recv9(:,:),send9_hi(:,:),recv9_lo(:,:)\n')
    rt = replace_once(rt, '    allocate(send9(im,9),recv9(im,9),stat=ierr)',
                     '    allocate(send9(im,9),recv9(im,9),send9_hi(im,9),recv9_lo(im,9),stat=ierr)')
    rt = replace_once(rt, '  if(allocated(send9))deallocate(send9,recv9)',
                     '  if(allocated(send9))deallocate(send9,recv9,send9_hi,recv9_lo)')
    start = rt.index('subroutine smpi_halo9_raw(')
    stop = rt.index('\nsubroutine smpi_halo9(', start)
    rt = rt[:start] + '''subroutine smpi_halo9_raw(a,b,c,d,e,f,g,h,v,im,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm),b(im,jm),c(im,jm),d(im,jm),e(im,jm),f(im,jm),g(im,jm),h(im,jm),v(im,jm)
#ifdef SECOM_USE_MPI
  integer::ierr,n
  type(MPI_Request)::requests(4)
  if(smpi_size<=1)return
  call dims(im,jm)
  if(.not.allocated(send9))then
    call smpi_halo_raw(a,im,jm)
    call smpi_halo_raw(b,im,jm)
    call smpi_halo_raw(c,im,jm)
    call smpi_halo_raw(d,im,jm)
    call smpi_halo_raw(e,im,jm)
    call smpi_halo_raw(f,im,jm)
    call smpi_halo_raw(g,im,jm)
    call smpi_halo_raw(h,im,jm)
    call smpi_halo_raw(v,im,jm)
    return
  endif
  n=9*im
  ! Both send buffers and both receive buffers are disjoint and live
  ! until Waitall. No caller may observe or modify these module buffers.
  if(south/=MPI_PROC_NULL)call pack(send9,smpi_jlo)
  if(north/=MPI_PROC_NULL)call pack(send9_hi,smpi_jhi)
  call MPI_Irecv(recv9,n,MPI_REAL,north,19,comm,requests(1),ierr)
  call checked(ierr,'HALO9 Irecv north')
  call MPI_Irecv(recv9_lo,n,MPI_REAL,south,20,comm,requests(2),ierr)
  call checked(ierr,'HALO9 Irecv south')
  call MPI_Isend(send9,merge(n,0,south/=MPI_PROC_NULL),MPI_REAL,south,19,comm,requests(3),ierr)
  call checked(ierr,'HALO9 Isend south')
  call MPI_Isend(send9_hi,merge(n,0,north/=MPI_PROC_NULL),MPI_REAL,north,20,comm,requests(4),ierr)
  call checked(ierr,'HALO9 Isend north')
  call MPI_Waitall(4,requests,MPI_STATUSES_IGNORE,ierr)
  call checked(ierr,'HALO9 Waitall')
  if(north/=MPI_PROC_NULL)call unpack(recv9,smpi_jhi+1)
  if(south/=MPI_PROC_NULL)call unpack(recv9_lo,smpi_jlo-1)
contains
  subroutine pack(buf,j)
    real,intent(out)::buf(im,9)
    integer,intent(in)::j
    buf(:,1)=a(:,j);buf(:,2)=b(:,j);buf(:,3)=c(:,j)
    buf(:,4)=d(:,j);buf(:,5)=e(:,j);buf(:,6)=f(:,j)
    buf(:,7)=g(:,j);buf(:,8)=h(:,j);buf(:,9)=v(:,j)
  end subroutine
  subroutine unpack(buf,j)
    real,intent(in)::buf(im,9)
    integer,intent(in)::j
    a(:,j)=buf(:,1);b(:,j)=buf(:,2);c(:,j)=buf(:,3)
    d(:,j)=buf(:,4);e(:,j)=buf(:,5);f(:,j)=buf(:,6)
    g(:,j)=buf(:,7);h(:,j)=buf(:,8);v(:,j)=buf(:,9)
  end subroutine
#endif
end subroutine
''' + rt[stop:]
    files['secom_mpi.f90'] = rt
    return files


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    files = duplex(a.source)
    a.out.mkdir()
    for name, text in files.items():
        b = text.encode('latin1')
        (a.out / name).write_bytes(b)
        print(hashlib.sha256(b).hexdigest(), name)
