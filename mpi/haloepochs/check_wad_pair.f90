! Standalone diagnostic: argv expected_ranks expected_threads JM {none,reverse}.
! "reverse" intentionally toggles BOTH on generation 2 (north-only -> both);
! only an exact omitted-south-ghost check, not a generic failure, is accepted.
program check_wad_pair
  use omp_lib, only: omp_get_num_threads
  implicit none
  include 'mpideck'
  integer,parameter::nx=7
  integer::nr,nt,ny,i,j,e,k,f,off,got,want,p,bad,layout,os,extra
  real,allocatable::a(:,:),b(:,:)
  logical::both,actual,mutant,owned,updated
  character(24)::arg
  call get_command_argument(1,arg);read(arg,*)nr
  call get_command_argument(2,arg);read(arg,*)nt
  call get_command_argument(3,arg);read(arg,*)ny
  call get_command_argument(4,arg)
  if(arg/='none'.and.arg/='reverse')stop 90
  mutant=arg=='reverse'
  call smpi_init(nx,ny)
  allocate(a(nx,ny),b(nx,ny))
  bad=0;extra=0;os=0
  do e=1,4
    both=mod(e,2)==1
    actual=both
    if(e==2.and.mutant)actual=.not.actual
    do j=1,ny
      do i=1,nx
        owned=j>=smpi_jlo.and.j<=smpi_jhi
        a(i,j)=transfer(raw(1,i,j,e,owned),0.)
        b(i,j)=transfer(raw(2,i,j,e,owned),0.)
      enddo
    enddo
!$omp parallel default(shared)
!$omp master
    os=omp_get_num_threads()
!$omp end master
    call smpi_wad_pair(a,b,nx,ny,actual)
!$omp end parallel
    do j=1,ny
      owned=j>=smpi_jlo.and.j<=smpi_jhi
      updated=owned.or.j==smpi_jhi+1.or.(both.and.j==smpi_jlo-1)
      do i=1,nx
        do f=1,2
          if(f==1)then
            got=transfer(a(i,j),0)
          else
            got=transfer(b(i,j),0)
          endif
          want=raw(f,i,j,e,updated)
          if(got/=want)then
            if(mutant.and.e==2.and.j==smpi_jlo-1.and.got==raw(f,i,j,e,.true.))then
              extra=extra+1
            else
              bad=bad+1
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  if(smpi_size/=nr.or.os/=nt)bad=bad+1
  call smpi_index_raw(bad)
  call smpi_index_raw(extra)
  k=smpi_rank
  call smpi_finalize()
  if(bad/=0)stop 45
  if(mutant)then
    if(extra/=2*nx)stop 46
    if(k==0)write(*,*)'DETECTED_WAD_PAIR_WRONG_DIRECTION extra_fields_cells=',extra
    stop 41
  endif
  if(k==0)write(*,*)'PASS_WAD_PAIR epochs/ranks/threads/JM=',4,nr,nt,ny
contains
  integer function raw(plane,ii,jj,epoch,is_owner_value)
    integer,intent(in)::plane,ii,jj,epoch
    logical,intent(in)::is_owner_value
    integer::salt,sign0
    salt=plane*16384+epoch*1024+jj*nx+ii
    if(.not.is_owner_value)then
      ! Rank-private poison includes field and generation. Far rows and
      ! the omitted south ghost must retain it byte for byte.
      raw=ior(int(z'7FC40000'),salt+smpi_rank*128)
    else if(ii==1)then
      raw=int(z'7FC10000')+salt
    else if(ii==2)then
      sign0=mod(plane+jj+epoch,2)
      raw=0
      if(sign0==1)raw=ibset(0,31)
    else
      raw=transfer(real(salt),0)
    endif
  end function
end program
