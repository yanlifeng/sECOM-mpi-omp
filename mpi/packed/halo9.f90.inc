
! All nine actual planes must be distinct and contiguous; buffers belong
! to the same serialized process/communicator as the old SMPI_HALO_RAW.
subroutine smpi_halo9_raw(a,b,c,d,e,f,g,h,v,im,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm),b(im,jm),c(im,jm),d(im,jm),e(im,jm),f(im,jm),g(im,jm),h(im,jm),v(im,jm)
#ifdef SECOM_USE_MPI
  if(smpi_size<=1)return
  call dims(im,jm)
  ! If 9*IM cannot be represented as a default MPI count, all ranks take
  ! the identical nine-call fallback. A failed allocation aborts at INIT.
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
  call exchange(south,north,19,smpi_jlo,smpi_jhi,1)
  call exchange(north,south,20,smpi_jhi,smpi_jlo,-1)
contains
  subroutine exchange(dest,src,tag,owned,edge,offset)
    integer,intent(in)::dest,src,tag,owned,edge,offset
    integer::ierr,n,nsend,ghost
    n=9*im
    ! The two-direction transfers finish in order and reuse two IM*9
    ! buffers. PROC_NULL receives do not touch any out-of-grid subscript.
    nsend=0
    if(dest/=MPI_PROC_NULL)then
      nsend=n
      send9(:,1)=a(:,owned);send9(:,2)=b(:,owned);send9(:,3)=c(:,owned)
      send9(:,4)=d(:,owned);send9(:,5)=e(:,owned);send9(:,6)=f(:,owned)
      send9(:,7)=g(:,owned);send9(:,8)=h(:,owned);send9(:,9)=v(:,owned)
    endif
    call MPI_Sendrecv(send9,nsend,MPI_REAL,dest,tag,recv9,n,MPI_REAL,src,tag,comm,MPI_STATUS_IGNORE,ierr)
    call checked(ierr,'HALO9 Sendrecv')
    if(src==MPI_PROC_NULL)return
    ghost=edge+offset
    a(:,ghost)=recv9(:,1);b(:,ghost)=recv9(:,2);c(:,ghost)=recv9(:,3)
    d(:,ghost)=recv9(:,4);e(:,ghost)=recv9(:,5);f(:,ghost)=recv9(:,6)
    g(:,ghost)=recv9(:,7);h(:,ghost)=recv9(:,8);v(:,ghost)=recv9(:,9)
  end subroutine
#endif
end subroutine

subroutine smpi_halo9(a,b,c,d,e,f,g,h,v,im,jm)
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm),b(im,jm),c(im,jm),d(im,jm),e(im,jm),f(im,jm),g(im,jm),h(im,jm),v(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_halo9_raw(a,b,c,d,e,f,g,h,v,im,jm)
!$omp end single
#endif
end subroutine
