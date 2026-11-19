
! Append after the parent runtime (in particular after smpi_private_rt).
! Only the rank's serialized MPI caller may use or reset this workspace.
module smpi_entry_workspace
  implicit none
#ifdef SECOM_USE_MPI
  real, allocatable, save :: epoch_send(:,:),epoch_recv(:,:)
#endif
end module

subroutine smpi_entry_epoch_release()
  use smpi_entry_workspace
  implicit none
#ifdef SECOM_USE_MPI
  if(allocated(epoch_send))deallocate(epoch_send,epoch_recv)
#endif
end subroutine

subroutine smpi_entry_epoch(d,ua,va,uab,vab,el,elb,dumwad,dvmwad,im,jm,first)
  implicit none
  integer,intent(in)::im,jm
  logical,intent(in)::first
  real,intent(inout)::d(im,jm),ua(im,jm),va(im,jm),uab(im,jm),vab(im,jm)
  real,intent(inout)::el(im,jm),elb(im,jm),dumwad(im,jm),dvmwad(im,jm)
#ifdef SECOM_USE_MPI
!$omp barrier
!$omp single
  call smpi_entry_epoch_raw(d,ua,va,uab,vab,el,elb,dumwad,dvmwad,im,jm,first)
!$omp end single
#endif
end subroutine

subroutine smpi_entry_epoch_raw(d,ua,va,uab,vab,el,elb,dumwad,dvmwad,im,jm,first)
  use smpi_private_rt
  use smpi_entry_workspace
  implicit none
  integer,intent(in)::im,jm
  logical,intent(in)::first
  real,intent(inout)::d(im,jm),ua(im,jm),va(im,jm),uab(im,jm),vab(im,jm)
  real,intent(inout)::el(im,jm),elb(im,jm),dumwad(im,jm),dvmwad(im,jm)
#ifdef SECOM_USE_MPI
  integer::lo,hi,ji,ns,nr,rc,nf
  if(smpi_size<=1)return
  call dims(im,jm)
  lo=smpi_jlo;hi=smpi_jhi
  ! Preserve exact *directional* ghost sets if a packed count would overflow.
  ! Full-bidirectional SMPI_HALO is not an admissible fallback here.
  if(im>huge(im)/6)then
    if(first)then
      call one(d,south,north,lo,hi,1,22)
      call one(va,south,north,lo,hi,1,22)
    endif
    call one(vab,south,north,lo,hi,1,22)
    if(first)call one(d,north,south,hi,lo,-1,23)
    call one(ua,north,south,hi,lo,-1,23)
    call one(uab,north,south,hi,lo,-1,23)
    if(first)call one(el,north,south,hi,lo,-1,23)
    call one(elb,north,south,hi,lo,-1,23)
    call one(dumwad,north,south,hi,lo,-1,23)
    return
  endif
  if(allocated(epoch_send))then
    if(size(epoch_send,1)/=im)call smpi_entry_epoch_release()
  endif
  if(.not.allocated(epoch_send))then
    allocate(epoch_send(im,6),epoch_recv(im,6),stat=rc)
    if(rc/=0)call bad('cannot allocate directional entry buffers')
  endif

  ! Pass 1, each low owned row travels south. Receive northern low row.
  nf=merge(3,1,first)
  ns=0;nr=merge(im*nf,0,north/=MPI_PROC_NULL)
  if(south/=MPI_PROC_NULL)then
    ns=im*nf
    epoch_send(:,1)=vab(:,lo)
    if(first)then
      epoch_send(:,2)=d(:,lo)
      epoch_send(:,3)=va(:,lo)
    endif
  endif
  call MPI_Sendrecv(epoch_send,ns,MPI_REAL,south,22,epoch_recv,nr,MPI_REAL,north,22,comm,MPI_STATUS_IGNORE,rc)
  call checked(rc,'ENTRY_EPOCH northern ghost')
  if(north/=MPI_PROC_NULL)then
    ji=hi+1
    vab(:,ji)=epoch_recv(:,1)
    if(first)then
      d(:,ji)=epoch_recv(:,2)
      va(:,ji)=epoch_recv(:,3)
    endif
  endif

  ! Pass 2, each high owned row travels north. Receive southern high row.
  nf=merge(6,4,first)
  ns=0;nr=merge(im*nf,0,south/=MPI_PROC_NULL)
  if(north/=MPI_PROC_NULL)then
    ns=im*nf
    epoch_send(:,1)=ua(:,hi)
    epoch_send(:,2)=uab(:,hi)
    epoch_send(:,3)=elb(:,hi)
    epoch_send(:,4)=dumwad(:,hi)
    if(first)then
      epoch_send(:,5)=d(:,hi)
      epoch_send(:,6)=el(:,hi)
    endif
  endif
  call MPI_Sendrecv(epoch_send,ns,MPI_REAL,north,23,epoch_recv,nr,MPI_REAL,south,23,comm,MPI_STATUS_IGNORE,rc)
  call checked(rc,'ENTRY_EPOCH southern ghost')
  if(south/=MPI_PROC_NULL)then
    ji=lo-1
    ua(:,ji)=epoch_recv(:,1)
    uab(:,ji)=epoch_recv(:,2)
    elb(:,ji)=epoch_recv(:,3)
    dumwad(:,ji)=epoch_recv(:,4)
    if(first)then
      d(:,ji)=epoch_recv(:,5)
      el(:,ji)=epoch_recv(:,6)
    endif
  endif
contains
  subroutine one(a,dst,src,sendrow,edge,offset,tag)
    integer,intent(in)::dst,src,sendrow,edge,offset,tag
    real,intent(inout)::a(im,jm)
    real::dummy
    integer::ierr
    if(src/=MPI_PROC_NULL)then
      call MPI_Sendrecv(a(1,sendrow),im,MPI_REAL,dst,tag,a(1,edge+offset),im,MPI_REAL,src,tag,comm,MPI_STATUS_IGNORE,ierr)
    else
      call MPI_Sendrecv(a(1,sendrow),im,MPI_REAL,dst,tag,dummy,0,MPI_REAL,src,tag,comm,MPI_STATUS_IGNORE,ierr)
    endif
    call checked(ierr,'ENTRY_EPOCH unpacked fallback')
  end subroutine
#endif
end subroutine
