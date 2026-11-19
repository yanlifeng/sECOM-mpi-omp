
! Append after smpi_private_rt; the parent supplies smpi_flux_one_raw.
module smpi_wpair_workspace
  implicit none
#ifdef SECOM_USE_MPI
  real, allocatable, save :: wpair_send(:,:),wpair_recv(:,:)
#endif
end module

subroutine smpi_wad_pair_release()
  use smpi_wpair_workspace
  implicit none
#ifdef SECOM_USE_MPI
  if(allocated(wpair_send))deallocate(wpair_send,wpair_recv)
#endif
end subroutine

subroutine smpi_wad_pair(a,b,im,jm,both)
  implicit none
  integer,intent(in)::im,jm
  logical,intent(in)::both
  real,intent(inout)::a(im,jm),b(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_wad_pair_raw(a,b,im,jm,both)
!$omp end single
#endif
end subroutine

subroutine smpi_wad_pair_raw(a,b,im,jm,both)
  use smpi_private_rt
  use smpi_wpair_workspace
  implicit none
  integer,intent(in)::im,jm
  logical,intent(in)::both
  real,intent(inout)::a(im,jm),b(im,jm)
#ifdef SECOM_USE_MPI
  integer::rc
  if(smpi_size<=1)return
  call dims(im,jm)
  if(im>huge(im)/2)then
    call smpi_flux_one_raw(a,im,jm,.true.)
    call smpi_flux_one_raw(b,im,jm,.true.)
    if(both)then
      call smpi_flux_one_raw(a,im,jm,.false.)
      call smpi_flux_one_raw(b,im,jm,.false.)
    endif
    return
  endif
  if(allocated(wpair_send))then
    if(size(wpair_send,1)/=im)call smpi_wad_pair_release()
  endif
  if(.not.allocated(wpair_send))then
    allocate(wpair_send(im,2),wpair_recv(im,2),stat=rc)
    if(rc/=0)call bad('cannot allocate WAD pair buffers')
  endif
  call pair_dir(south,north,smpi_jlo,smpi_jhi,1,45)
  if(both)call pair_dir(north,south,smpi_jhi,smpi_jlo,-1,46)
contains
  subroutine pair_dir(dst,src,owned,edge,delta,tag)
    integer,intent(in)::dst,src,owned,edge,delta,tag
    integer::ns,nr,ghost,ierr
    ns=0;nr=merge(2*im,0,src/=MPI_PROC_NULL)
    if(dst/=MPI_PROC_NULL)then
      ns=2*im
      wpair_send(:,1)=a(:,owned)
      wpair_send(:,2)=b(:,owned)
    endif
    call MPI_Sendrecv(wpair_send,ns,MPI_REAL,dst,tag,wpair_recv,nr,MPI_REAL,src,tag,comm,MPI_STATUS_IGNORE,ierr)
    call checked(ierr,'paired WAD ghost')
    if(src==MPI_PROC_NULL)return
    ghost=edge+delta
    a(:,ghost)=wpair_recv(:,1)
    b(:,ghost)=wpair_recv(:,2)
  end subroutine
#endif
end subroutine
