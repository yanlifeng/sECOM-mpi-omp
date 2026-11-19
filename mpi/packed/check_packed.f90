program check_packed
  use omp_lib
  use mpi_f08, only: MPI_Abort,MPI_COMM_WORLD,MPI_Query_thread,MPI_THREAD_SERIALIZED
  use iso_c_binding, only:c_int
  implicit none
  include 'mpideck'
  interface
    integer(c_int) function place(rank,lane,nt) bind(c)
      import c_int
      integer(c_int),value::rank,lane,nt
    end function
  end interface
  integer, parameter :: im=372
  integer :: jm,nr,nt,i,j,k,t,want,level,rc,cp(8),lane,owner,lo,hi
  real, allocatable :: old(:,:,:),new(:,:,:)
  character(16) :: arg
  call get_command_argument(1,arg);read(arg,*)jm
  call get_command_argument(2,arg);read(arg,*)nr
  call get_command_argument(3,arg);read(arg,*)nt
  call smpi_init(im,jm)
  call MPI_Query_thread(level,rc)
  if(rc/=0.or.level<MPI_THREAD_SERIALIZED.or.smpi_size/=nr.or.nt>8)call fail(70)
  allocate(old(im,jm,9),new(im,jm,9))
  lo=smpi_jlo;hi=smpi_jhi;cp=-1
!$omp parallel default(shared) private(i,j,k,t,lane,owner,want)
  lane=omp_get_thread_num()
  if(omp_get_num_threads()/=nt.or.lane>=8)call fail(71)
  cp(lane+1)=place(smpi_rank,lane,nt)
  if(cp(lane+1)<0)call fail(72)
  do t=1,12
!$omp do collapse(2) schedule(static)
    do k=1,9
      do j=1,jm
        do i=1,im
          want=poison(i,j,k,t,smpi_rank)
          if(j>=lo.and.j<=hi)want=bits(i,j,k,t,smpi_rank)
          old(i,j,k)=transfer(want,0.)
          new(i,j,k)=transfer(want,0.)
        enddo
      enddo
    enddo
!$omp end do nowait
    ! The missing producer join is intentional: entry joins in the public
    ! wrappers must publish the new arrays from every OpenMP lane.
    if(mod(t,2)==1)then
      call group_public()
    endif
    do k=1,9
      call smpi_halo(old(1,1,k),im,jm)
    enddo
    if(mod(t,2)==0)then
!$omp single
      call group_raw()
!$omp end single
    endif
!$omp do collapse(2) schedule(static)
    do k=1,9
      do j=1,jm
        owner=smpi_owner_raw(j)
        do i=1,im
          want=poison(i,j,k,t,smpi_rank)
          if(j>=max(1,lo-1).and.j<=min(jm,hi+1))want=bits(i,j,k,t,owner)
          if(transfer(old(i,j,k),0)/=want.or.transfer(new(i,j,k),0)/=want)then
            write(*,*) 'packed mismatch: rank,lane,t,field,j,i',smpi_rank,lane,t,k,j,i
            call fail(73)
          endif
        enddo
      enddo
    enddo
!$omp end do
  enddo
!$omp end parallel
  write(*,*) 'HALO9_PASS rank,nr,jlo,jhi,im,jm,omp,provided,cpu:',smpi_rank,nr,lo,hi,im,jm,nt,level,cp(1:nt)
  call smpi_finalize()
contains
  subroutine group_public()
    call smpi_halo9(new(1,1,1),new(1,1,2),new(1,1,3),new(1,1,4),new(1,1,5),new(1,1,6), &
                    new(1,1,7),new(1,1,8),new(1,1,9),im,jm)
  end subroutine
  subroutine group_raw()
    call smpi_halo9_raw(new(1,1,1),new(1,1,2),new(1,1,3),new(1,1,4),new(1,1,5),new(1,1,6), &
                        new(1,1,7),new(1,1,8),new(1,1,9),im,jm)
  end subroutine
  integer function bits(ii,jj,kk,it,peer)
    integer,intent(in)::ii,jj,kk,it,peer
    integer::salt
    salt=ii+401*(jj+17*(kk+11*(it+13*peer)))
    select case(mod(ii+jj+kk+it,32))
    case(0);bits=0
    case(1);bits=ibset(0,31)
    case(2);bits=1
    case(3);bits=ibset(1,31)
    case(4);bits=ior(int(z'7FC00000'),iand(salt,int(z'001FFFFF')))
    case(5);bits=ior(int(z'7F800001'),iand(salt,int(z'001FFFFE')))
    case default;bits=ior(int(z'3F000000'),iand(salt,int(z'007FFFFF')))
    end select
  end function
  integer function poison(ii,jj,kk,it,peer)
    integer,intent(in)::ii,jj,kk,it,peer
    poison=ior(ibset(int(z'7FA00000'),31),iand(ii+379*(jj+17*(kk+11*(it+13*peer))),int(z'000FFFFF')))
  end function
  subroutine fail(code)
    integer,intent(in)::code
    integer::ierr
    call MPI_Abort(MPI_COMM_WORLD,code,ierr)
    error stop 1
  end subroutine
end program
