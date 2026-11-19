! Standalone directional-ghost test. Launch under *approved* rank placement
! with argv: JM size threads. Checks the touched AND untouched cells bitwise.
program check_entry_epoch
  use omp_lib
  use mpi_f08,only:MPI_COMM_WORLD,MPI_Abort
  implicit none
  include 'mpideck'
  integer,parameter::im=372
  real,allocatable::a(:,:,:)
  integer::jm,nr,nt,i,j,k,t,old,want,lo,hi,who,rc
  logical::first,touched
  character(20)::arg
  call get_command_argument(1,arg);read(arg,*)jm
  call get_command_argument(2,arg);read(arg,*)nr
  call get_command_argument(3,arg);read(arg,*)nt
  call smpi_init(im,jm)
  if(smpi_size/=nr)call MPI_Abort(MPI_COMM_WORLD,70,rc)
  lo=smpi_jlo;hi=smpi_jhi
  allocate(a(im,jm,9))
  do t=1,12
    first=mod(t,3)==1
!$omp parallel default(shared) private(i,j,k,old,want,touched,who,rc)
    if(omp_get_num_threads()/=nt)call MPI_Abort(MPI_COMM_WORLD,71,rc)
!$omp do collapse(2) schedule(static)
    do k=1,9
      do j=1,jm
        do i=1,im
          old=payload(i,j,k,t,smpi_rank,j<lo.or.j>hi)
          a(i,j,k)=transfer(old,0.)
        enddo
      enddo
    enddo
!$omp end do nowait
    if(mod(t,2)==0)then
!$omp barrier
!$omp single
      call smpi_entry_epoch_raw(a(1,1,1),a(1,1,2),a(1,1,3),a(1,1,4),a(1,1,5),a(1,1,6),a(1,1,7),a(1,1,8),a(1,1,9),im,jm,first)
!$omp end single
    else
      call smpi_entry_epoch(a(1,1,1),a(1,1,2),a(1,1,3),a(1,1,4),a(1,1,5),a(1,1,6),a(1,1,7),a(1,1,8),a(1,1,9),im,jm,first)
    endif
!$omp do collapse(2) schedule(static)
    do k=1,9
      do j=1,jm
        who=smpi_rank
        touched=j>=lo.and.j<=hi
        if(j==lo-1)then
          touched=any(k==[2,4,7,8]).or.(first.and.any(k==[1,6]))
          if(touched)who=smpi_rank-1
        else if(j==hi+1)then
          touched=k==5.or.(first.and.any(k==[1,3]))
          if(touched)who=smpi_rank+1
        endif
        do i=1,im
          want=payload(i,j,k,t,who,.not.touched)
          if(transfer(a(i,j,k),0)/=want)then
            print *, 'BAD ENTRY EPOCH',smpi_rank,t,k,i,j
            call MPI_Abort(MPI_COMM_WORLD,72,rc)
          endif
        enddo
      enddo
    enddo
!$omp end do
!$omp end parallel
    if(t==6)call smpi_entry_epoch_release()
  enddo
  print *, 'ENTRY_EPOCH_PASS rank,size,lo,hi,im,jm,nt',smpi_rank,nr,lo,hi,im,jm,nt
  call smpi_entry_epoch_release()
  call smpi_finalize()
contains
  integer function payload(ii,jj,kk,tt,rank,poison)
    integer,intent(in)::ii,jj,kk,tt,rank
    logical,intent(in)::poison
    integer::key
    key=ii+379*(jj+17*(kk+11*(tt+13*rank)))
    if(poison)then
      payload=ior(ibset(int(z'7fa00000'),31),iand(key,int(z'000fffff')))
    else
      select case(mod(ii+kk+tt,64))
      case(0);payload=0
      case(1);payload=ibset(0,31)
      case(2);payload=1
      case(3);payload=ior(int(z'7f800001'),iand(key,int(z'001ffffe')))
      case(4);payload=ior(int(z'7fc00000'),iand(key,int(z'001fffff')))
      case default;payload=ior(int(z'3f000000'),iand(key,int(z'007fffff')))
      end select
    endif
  end function
end program
