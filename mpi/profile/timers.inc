
! Dedicated diagnostic ELF: no profiling in an ordinary science build.
module secom_diag_store
  use omp_lib
  use iso_fortran_env, only: int64,real64
  implicit none
  integer,parameter :: nl=5, ns=14, nm=8
  character(18),parameter :: life_name(0:nl)=[character(18):: &
    'outside','mpi-init','initialization','iteration','post-loop','mpi-finalize']
  character(18),parameter :: site_name(0:ns)=[character(18):: &
    'other','bcond7','cbc','bottau','archive','baro','restart-window','putcdf','prints','tscdf', &
    'extrnl','wetdry','check','advave','smag']
  character(18),parameter :: op_name(nm)=[character(18):: &
    'halo-sendrecv','plane-allgatherv','row-allgatherv','irow-allgatherv', &
    'index-allreduce','row-bcast','wad-sendrecv','halo9-sendrecv']
  integer,save::life=0,depth=0,stack(8)=0,saved_rank=-1,nr=0,nt=0,subsite=0
!$omp threadprivate(subsite)
  logical,save :: passed_enddo=.false.
  integer(int64),save::iter=0,scount(0:nl,ns)=0,mcount(0:nl,0:ns,nm)=0
  real(real64),save::lt(nl)=0,st(0:nl,ns)=0,mt(0:nl,0:ns,nm)=0,mark_l=0,mark_s(ns)=0,mark_m=0
contains
  subroutine record(u,kind,s,x,o,n,t)
    integer,intent(in)::u
    character(*),intent(in)::kind,s,x,o
    integer(int64),intent(in)::n
    real(real64),intent(in)::t
    write(u,'(4(a,1x),i0,1x,es24.16e3)') kind,trim(s),trim(x),trim(o),n,t
  end subroutine
end module
subroutine sd_id(rank,size)
  use secom_diag_store
  integer,intent(in)::rank,size
  saved_rank=rank;nr=size
end subroutine
subroutine sd_stage(id)
  use secom_diag_store
  integer,intent(in)::id
  real(real64)::now
  now=omp_get_wtime()
  if(life>0)lt(life)=lt(life)+(now-mark_l)
  if(depth/=0)error stop 'profile: stage changed with open region'
  if(id==4)passed_enddo=.true.
  life=id;mark_l=now
end subroutine
subroutine sd_step()
  use secom_diag_store
  iter=iter+1
end subroutine
subroutine sd_site(id,enter)
  use secom_diag_store
  integer,intent(in)::id,enter
  real(real64)::now
  if(id>=10.and.id<=12)then
    ! Every caller labels its own lane. Only the master times the stage;
    ! no new team join. At 32x1 this is exactly the rank's stage time.
    subsite=merge(id,0,enter==1)
    if(omp_get_thread_num()/=0)return
    nt=max(nt,omp_get_num_threads())
  else
    ! Serial caller owns the stack. Stack top is immutable during teams.
    if(enter==1)then
      depth=depth+1
      if(depth>size(stack))error stop 'profile: scope depth'
      stack(depth)=id
    else
      if(depth<1)error stop 'profile: empty scope'
      if(stack(depth)/=id)error stop 'profile: scope imbalance'
      depth=depth-1
    endif
  endif
  now=omp_get_wtime()
  if(enter==1)then
    mark_s(id)=now
  else
    st(life,id)=st(life,id)+(now-mark_s(id))
    scount(life,id)=scount(life,id)+1
  endif
end subroutine
subroutine sd_mp(id,enter)
  use secom_diag_store
  integer,intent(in)::id,enter
  integer::site
  real(real64)::now
  ! RAW callers are already serialized by SMPI's existing contract: the
  ! start and end are on the very same thread; no extra MPI/OMP barriers.
  now=omp_get_wtime()
  if(enter==1)then
    mark_m=now
  else
    site=subsite
    if(site==0.and.depth>0)site=stack(depth)
    mt(life,site,id)=mt(life,site,id)+(now-mark_m)
    mcount(life,site,id)=mcount(life,site,id)+1
  endif
end subroutine
subroutine sd_dump()
  use secom_diag_store
  integer::u,i,j,k
  character(80)::name
  call sd_stage(0)
  write(name,'(a,i6.6,a)') 'smpi-profile-',saved_rank,'.log'
  open(newunit=u,file=trim(name),status='new',action='write')
  write(u,'(a,5(1x,i0))') 'SECOM_DIAG_V1',saved_rank,nr,nt,iter,merge(1,0,passed_enddo)
  do i=1,nl
    call record(u,'life',life_name(i),'all','wall',1_int64,lt(i))
    do j=1,ns
      if(scount(i,j)>0)call record(u,'site',life_name(i),site_name(j),'inclusive',scount(i,j),st(i,j))
    enddo
    do j=0,ns
      do k=1,nm
        if(mcount(i,j,k)>0)call record(u,'mpi',life_name(i),site_name(j),op_name(k),mcount(i,j,k),mt(i,j,k))
      enddo
    enddo
  enddo
  close(u)
end subroutine
