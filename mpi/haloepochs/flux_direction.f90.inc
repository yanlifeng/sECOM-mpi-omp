
! Private helpers for three separate producer epochs of the same scratch.
! recv_north means consume J+1: send first owned row south. For J-1 the
! source, destination, owned send row and ghost receive row all reverse.
subroutine smpi_flux_one_raw(a,im,jm,recv_north)
  use smpi_private_rt
  implicit none
  integer,intent(in)::im,jm
  logical,intent(in)::recv_north
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
  integer::dest,src,sendj,ghost,tag,ierr
  real::unused
  if(smpi_size<=1)return
  call dims(im,jm)
  if(recv_north)then
    src=north; dest=south; sendj=smpi_jlo; tag=43
    if(src/=MPI_PROC_NULL)ghost=smpi_jhi+1
  else
    src=south; dest=north; sendj=smpi_jhi; tag=44
    if(src/=MPI_PROC_NULL)ghost=smpi_jlo-1
  end if
  if(src==MPI_PROC_NULL)then
    ! No out-of-range actual argument, including 1-row exterior slabs.
    call MPI_Sendrecv(a(1,sendj),im,MPI_REAL,dest,tag,unused,0,MPI_REAL,src,tag,comm,MPI_STATUS_IGNORE,ierr)
  else
    call MPI_Sendrecv(a(1,sendj),im,MPI_REAL,dest,tag,a(1,ghost),im,MPI_REAL,src,tag,comm,MPI_STATUS_IGNORE,ierr)
  end if
  call checked(ierr,'directed FLUXVA ghost')
#endif
end subroutine

subroutine smpi_flux_halo_n(a,im,jm)
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_flux_one_raw(a,im,jm,.true.)
!$omp end single
#endif
end subroutine

subroutine smpi_flux_halo_s(a,im,jm)
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_flux_one_raw(a,im,jm,.false.)
!$omp end single
#endif
end subroutine
