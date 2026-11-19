
! Raw collectives: both called only by a SINGLE-selected thread per rank.
subroutine smpi_checkrec_raw(local,allrec)
  use smpi_private_rt
  implicit none
  integer,intent(in)::local(4)
  integer,intent(out)::allrec(4,*)
#ifdef SECOM_USE_MPI
  integer::ierr
  call MPI_Allgather(local,4,MPI_INTEGER,allrec,4,MPI_INTEGER,comm,ierr)
  call checked(ierr,'CHECK rank records')
#else
  allrec(:,1)=local
#endif
end subroutine

subroutine smpi_checkpicks_raw(picks)
  use smpi_private_rt
  implicit none
  integer,intent(inout)::picks(3)
#ifdef SECOM_USE_MPI
  integer::ierr
  ! At most one owner sets each slot; other ranks supply all-zero bits.
  call MPI_Allreduce(MPI_IN_PLACE,picks,3,MPI_INTEGER,MPI_BOR,comm,ierr)
  call checked(ierr,'CHECK selected cells')
#endif
end subroutine
