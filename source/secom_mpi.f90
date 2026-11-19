! Full-sized legacy planes, decomposed by contiguous global J rows.
! Compile with -cpp; only MODE=mpi defines SECOM_USE_MPI. No physics here.
module smpi_private_rt
#ifdef SECOM_USE_MPI
  use mpi_f08
#endif
  implicit none
  integer :: smpi_rank,smpi_size,smpi_jlo,smpi_jhi
  common /secom_mpislab/ smpi_rank,smpi_size,smpi_jlo,smpi_jhi
  save /secom_mpislab/
  logical :: smpi_localbc
  common /secom_mpi_bc_flag/ smpi_localbc
  save /secom_mpi_bc_flag/
  integer, save :: nx=0, ny=0
#ifdef SECOM_USE_MPI
  type(MPI_Comm), save :: comm
  logical, save :: owns_mpi=.false.
  integer, allocatable, save :: rows(:),offs(:),cells(:),starts(:),wad_recv(:)
  integer, save :: south,north
#endif
contains
  subroutine bad(reason)
    use iso_fortran_env, only: error_unit
    character(*),intent(in)::reason
#ifdef SECOM_USE_MPI
    integer::ierr
    logical::ready,done
#endif
    write(error_unit,'(a,i0,2a)') 'SMPI rank ',smpi_rank,': ',trim(reason)
#ifdef SECOM_USE_MPI
    call MPI_Initialized(ready,ierr)
    call MPI_Finalized(done,ierr)
    if(ready.and..not.done)call MPI_Abort(MPI_COMM_WORLD,93,ierr)
#endif
    error stop 93
  end subroutine
#ifdef SECOM_USE_MPI
  subroutine checked(rc,where)
    integer,intent(in)::rc
    character(*),intent(in)::where
    if(rc/=MPI_SUCCESS)call bad(where)
  end subroutine
  integer function owner(j)
    integer,intent(in)::j
    integer::base,nextra,nlong
    if(j<1.or.j>ny)call bad('row owner index outside global grid')
    base=ny/smpi_size; nextra=mod(ny,smpi_size); nlong=nextra*(base+1)
    if(j<=nlong)then
      owner=(j-1)/(base+1)
    else
      owner=nextra+(j-1-nlong)/base
    endif
  end function
#endif
  subroutine dims(im,jm)
    integer,intent(in)::im,jm
    if(nx/=im.or.ny/=jm)call bad('helper grid differs from initialization')
  end subroutine
end module

block data smpi_default_bounds
  implicit none
  integer::r,n,lo,hi
  common /secom_mpislab/ r,n,lo,hi
  logical::localbc
  common /secom_mpi_bc_flag/ localbc
  data r,n,lo,hi /0,1,1,2147483647/
  data localbc /.false./
end block data

subroutine smpi_init(im,jm)
  use smpi_private_rt
!$ use omp_lib, only: omp_get_level
  implicit none
  integer,intent(in)::im,jm
#ifdef SECOM_USE_MPI
  integer::provided,ierr,k,local(2),small(2),large(2)
  logical::ready,done,onmain
#endif
!$ if(omp_get_level()/=0)call bad('INIT must be called outside any parallel region')
  if(nx/=0)call bad('duplicate INIT')
  smpi_localbc=.false.
  if(im<1.or.jm<1)call bad('invalid dimensions')
  if(jm>huge(jm)/im)call bad('32-bit MPI count overflow')
#ifdef SECOM_USE_MPI
  call MPI_Finalized(done,ierr);call checked(ierr,'MPI_Finalized')
  if(done)call bad('MPI cannot be reinitialized after finalize')
  call MPI_Initialized(ready,ierr);call checked(ierr,'MPI_Initialized')
  if(.not.ready)then
    call MPI_Init_thread(MPI_THREAD_SERIALIZED,provided,ierr);call checked(ierr,'MPI_Init_thread')
    owns_mpi=.true.
  else
    call MPI_Query_thread(provided,ierr);call checked(ierr,'MPI_Query_thread')
  endif
  if(provided<MPI_THREAD_SERIALIZED)call bad('MPI_THREAD_SERIALIZED not provided')
  call MPI_Is_thread_main(onmain,ierr);call checked(ierr,'MPI_Is_thread_main')
  if(.not.onmain)call bad('INIT not on initial MPI thread')
  call MPI_Comm_dup(MPI_COMM_WORLD,comm,ierr);call checked(ierr,'MPI_Comm_dup')
  call MPI_Comm_set_errhandler(comm,MPI_ERRORS_RETURN,ierr);call checked(ierr,'MPI_Comm_set_errhandler')
  call MPI_Comm_size(comm,smpi_size,ierr);call checked(ierr,'MPI_Comm_size')
  call MPI_Comm_rank(comm,smpi_rank,ierr);call checked(ierr,'MPI_Comm_rank')
  local=[im,jm]
  call MPI_Allreduce(local,small,2,MPI_INTEGER,MPI_MIN,comm,ierr);call checked(ierr,'grid minimum')
  call MPI_Allreduce(local,large,2,MPI_INTEGER,MPI_MAX,comm,ierr);call checked(ierr,'grid maximum')
  if(any(small/=large).or.smpi_size>jm)call bad('inconsistent grid or empty slabs')
  if(storage_size(1.)/=32.or.storage_size(1)/=32)call bad('expected default real32 and integer32')
  allocate(rows(0:smpi_size-1),offs(0:smpi_size-1),cells(0:smpi_size-1),starts(0:smpi_size-1),wad_recv(im))
  do k=0,smpi_size-1
    rows(k)=jm/smpi_size+merge(1,0,k<mod(jm,smpi_size))
    offs(k)=k*(jm/smpi_size)+min(k,mod(jm,smpi_size))
  enddo
  cells=im*rows;starts=im*offs
  smpi_jlo=offs(smpi_rank)+1;smpi_jhi=offs(smpi_rank)+rows(smpi_rank)
  south=merge(smpi_rank-1,MPI_PROC_NULL,smpi_rank>0)
  north=merge(smpi_rank+1,MPI_PROC_NULL,smpi_rank+1<smpi_size)
#else
  smpi_rank=0;smpi_size=1;smpi_jlo=1;smpi_jhi=jm
#endif
  nx=im;ny=jm
end subroutine

subroutine smpi_finalize()
  use smpi_private_rt
!$ use omp_lib, only: omp_get_level
  implicit none
#ifdef SECOM_USE_MPI
  integer::ierr
  logical::onmain
#endif
!$ if(omp_get_level()/=0)call bad('FINALIZE must be called outside any parallel region')
  if(nx==0)return
  smpi_localbc=.false.
#ifdef SECOM_USE_MPI
  call MPI_Is_thread_main(onmain,ierr);call checked(ierr,'MPI_Is_thread_main')
  if(.not.onmain)call bad('FINALIZE not on initial MPI thread')
  call MPI_Comm_free(comm,ierr);call checked(ierr,'MPI_Comm_free')
  deallocate(rows,offs,cells,starts,wad_recv)
  if(owns_mpi)then
    call MPI_Finalize(ierr);call checked(ierr,'MPI_Finalize')
  endif
  owns_mpi=.false.
#endif
  smpi_rank=0;smpi_size=1;smpi_jlo=1;smpi_jhi=huge(smpi_jhi);nx=0;ny=0
end subroutine
integer function smpi_owner_raw(j)
  use smpi_private_rt
  implicit none
  integer,intent(in)::j
  if(j<1.or.j>ny)call bad('OWNER_RAW requires initialized global row')
#ifdef SECOM_USE_MPI
  smpi_owner_raw=owner(j)
#else
  smpi_owner_raw=0
#endif
end function
integer function smpi_first(lo)
  use smpi_private_rt,only:smpi_jlo
  integer,intent(in)::lo
  smpi_first=max(lo,smpi_jlo)
end function
integer function smpi_last(hi)
  use smpi_private_rt,only:smpi_jhi
  integer,intent(in)::hi
  smpi_last=min(hi,smpi_jhi)
end function
logical function smpi_owns(j)
  use smpi_private_rt,only:smpi_jlo,smpi_jhi
  integer,intent(in)::j
  smpi_owns=j>=smpi_jlo.and.j<=smpi_jhi
end function
integer function smpi_color(lo,hi,stride)
  use smpi_private_rt,only:smpi_jlo,smpi_jhi,bad
  implicit none
  integer,intent(in)::lo,hi,stride
  integer::j
  if(stride<1)call bad('COLOR needs positive stride')
  j=max(lo,smpi_jlo);j=j+modulo(lo-j,stride)
  smpi_color=j
  if(j>min(hi,smpi_jhi))smpi_color=hi+1
end function

subroutine smpi_halo_raw(a,im,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
  integer::ierr
  real::dummy
  if(smpi_size<=1)return
  call dims(im,jm)
  ! Use scalar scratch with receive count zero at exterior boundaries,
  ! never alias a send buffer or form an invalid subscript for PROC_NULL.
  if(north/=MPI_PROC_NULL)then
    call MPI_Sendrecv(a(1,smpi_jlo),im,MPI_REAL,south,1,a(1,smpi_jhi+1),im,MPI_REAL,north,1,comm,MPI_STATUS_IGNORE,ierr)
  else
    call MPI_Sendrecv(a(1,smpi_jlo),im,MPI_REAL,south,1,dummy,0,MPI_REAL,north,1,comm,MPI_STATUS_IGNORE,ierr)
  endif
  call checked(ierr,'HALO north ghost')
  if(south/=MPI_PROC_NULL)then
    call MPI_Sendrecv(a(1,smpi_jhi),im,MPI_REAL,north,2,a(1,smpi_jlo-1),im,MPI_REAL,south,2,comm,MPI_STATUS_IGNORE,ierr)
  else
    call MPI_Sendrecv(a(1,smpi_jhi),im,MPI_REAL,north,2,dummy,0,MPI_REAL,south,2,comm,MPI_STATUS_IGNORE,ierr)
  endif
  call checked(ierr,'HALO south ghost')
#endif
end subroutine
subroutine smpi_all_raw(a,im,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
  integer::ierr
  if(smpi_size<=1)return
  call dims(im,jm)
  call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_REAL,a,cells,starts,MPI_REAL,comm,ierr);call checked(ierr,'ALL')
#endif
end subroutine
subroutine smpi_rows_raw(v,e,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::jm
  real,intent(inout)::v(jm),e(jm)
#ifdef SECOM_USE_MPI
  integer::ierr
  if(smpi_size<=1)return
  call dims(nx,jm)
  call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_REAL,v,rows,offs,MPI_REAL,comm,ierr);call checked(ierr,'ROWS vmax')
  call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_REAL,e,rows,offs,MPI_REAL,comm,ierr);call checked(ierr,'ROWS emin')
#endif
end subroutine
subroutine smpi_irows_raw(v,e,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::jm
  integer,intent(inout)::v(jm),e(jm)
#ifdef SECOM_USE_MPI
  integer::ierr
  if(smpi_size<=1)return
  call dims(nx,jm)
  call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_INTEGER,v,rows,offs,MPI_INTEGER,comm,ierr);call checked(ierr,'IROWS imax')
  call MPI_Allgatherv(MPI_IN_PLACE,0,MPI_INTEGER,e,rows,offs,MPI_INTEGER,comm,ierr);call checked(ierr,'IROWS imin')
#endif
end subroutine
subroutine smpi_index_raw(i)
  use smpi_private_rt
  implicit none
  integer,intent(inout)::i
#ifdef SECOM_USE_MPI
  integer::ierr
  if(smpi_size<=1)return
  call MPI_Allreduce(MPI_IN_PLACE,i,1,MPI_INTEGER,MPI_MAX,comm,ierr);call checked(ierr,'INDEX')
#endif
end subroutine
subroutine smpi_rbcast_row_raw(x,j)
  use smpi_private_rt
  implicit none
  integer,intent(in)::j
  real,intent(inout)::x
#ifdef SECOM_USE_MPI
  integer::ierr
  if(smpi_size<=1)return
  call MPI_Bcast(x,4,MPI_BYTE,owner(j),comm,ierr);call checked(ierr,'RBCAST_ROW')
#endif
end subroutine
subroutine smpi_wad_publish_raw(flags,et,etf,va,vaf,dvmwad,etb,im,jm)
  use smpi_private_rt
  implicit none
  integer,intent(in)::im,jm,flags(im)
  real,intent(inout)::et(im,jm),etf(im,jm),va(im,jm),vaf(im,jm),dvmwad(im,jm)
  real,intent(in)::etb(im,jm)
#ifdef SECOM_USE_MPI
  integer::ierr,i,j
  if(smpi_size<=1)return
  call dims(im,jm)
  call MPI_Sendrecv(flags,im,MPI_INTEGER,north,3,wad_recv,im,MPI_INTEGER,south,3,comm,MPI_STATUS_IGNORE,ierr)
  call checked(ierr,'WAD seam publish')
  if(south==MPI_PROC_NULL)return
  j=smpi_jlo
  do i=1,im
    if(btest(wad_recv(i),0))then
      va(i,j)=0.;vaf(i,j)=0.;dvmwad(i,j)=0.
    endif
    if(btest(wad_recv(i),1))then
      et(i,j)=etb(i,j);etf(i,j)=etb(i,j)
    endif
  enddo
#endif
end subroutine
subroutine smpi_halo(a,im,jm)
  use smpi_private_rt, only:smpi_size
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_halo_raw(a,im,jm)
!$omp end single
#endif
end subroutine
subroutine smpi_all(a,im,jm)
  use smpi_private_rt, only:smpi_size
  implicit none
  integer,intent(in)::im,jm
  real,intent(inout)::a(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_all_raw(a,im,jm)
!$omp end single
#endif
end subroutine
subroutine smpi_rows(v,e,jm)
  use smpi_private_rt, only:smpi_size
  implicit none
  integer,intent(in)::jm
  real,intent(inout)::v(jm),e(jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_rows_raw(v,e,jm)
!$omp end single
#endif
end subroutine
subroutine smpi_irows(v,e,jm)
  use smpi_private_rt, only:smpi_size
  implicit none
  integer,intent(in)::jm
  integer,intent(inout)::v(jm),e(jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_irows_raw(v,e,jm)
!$omp end single
#endif
end subroutine
subroutine smpi_index(i)
  use smpi_private_rt, only:smpi_size
  implicit none
  integer,intent(inout)::i
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_index_raw(i)
!$omp end single
#endif
end subroutine
subroutine smpi_rbcast_row(x,j)
  use smpi_private_rt, only:smpi_size
  implicit none
  real,intent(inout)::x
  integer,intent(in)::j
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_rbcast_row_raw(x,j)
!$omp end single
#endif
end subroutine
subroutine smpi_wad_publish(flags,et,etf,va,vaf,dvmwad,etb,im,jm)
  use smpi_private_rt, only:smpi_size
  implicit none
  integer,intent(in)::im,jm,flags(im)
  real,intent(inout)::et(im,jm),etf(im,jm),va(im,jm),vaf(im,jm),dvmwad(im,jm)
  real,intent(in)::etb(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_wad_publish_raw(flags,et,etf,va,vaf,dvmwad,etb,im,jm)
!$omp end single
#endif
end subroutine
