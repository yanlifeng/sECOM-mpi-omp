
subroutine smpi_checkkeys_raw(keys)
  use smpi_private_rt
  use iso_fortran_env, only:int64
  implicit none
  integer(int64),intent(inout)::keys(3)
#ifdef SECOM_USE_MPI
  integer::ierr
  call MPI_Allreduce(MPI_IN_PLACE,keys,3,MPI_INTEGER8,MPI_MAX,comm,ierr)
  call checked(ierr,'CHECK healthy keys')
#endif
end subroutine
