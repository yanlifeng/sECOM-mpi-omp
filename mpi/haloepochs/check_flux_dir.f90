program check_flux_dir
  use iso_fortran_env, only:int32
  use omp_lib, only:omp_get_num_threads
  implicit none
  include 'mpideck'
  integer,parameter::im=7
  integer::ny,nr,nt,witness,j,i,epoch,bad,pbad,poison,expected,raw,lo,hi,k
  logical::north,actual,wrong
  real,allocatable::a(:,:)
  real::div,ref
  character(24)::arg
  call get_command_argument(1,arg);read(arg,*)nr
  call get_command_argument(2,arg);read(arg,*)nt
  call get_command_argument(3,arg);read(arg,*)ny
  call get_command_argument(4,arg);wrong=arg=='reverse'
  call smpi_init(im,ny)
  allocate(a(im,ny))
  bad=0; pbad=0; witness=0
  do epoch=1,3
    north=epoch/=3
    actual=north
    if(wrong.and.epoch==2)actual=.not.actual
    poison=int(z'7FC01000')+epoch*32+smpi_rank
    a=transfer(poison,0.)
    do j=smpi_jlo,smpi_jhi
      do i=1,im
        a(i,j)=transfer(bits(i,j,epoch),0.)
      end do
    end do
!$omp parallel default(shared)
!$omp master
    witness=omp_get_num_threads()
!$omp end master
    if(actual)then
      call smpi_flux_halo_n(a,im,ny)
    else
      call smpi_flux_halo_s(a,im,ny)
    endif
!$omp end parallel
    do j=1,ny
      do i=1,im
        expected=poison
        if((j>=smpi_jlo.and.j<=smpi_jhi).or.(north.and.j==smpi_jhi+1).or. &
           (.not.north.and.j==smpi_jlo-1))expected=bits(i,j,epoch)
        raw=transfer(a(i,j),0)
        if(raw/=expected)pbad=pbad+1
      enddo
    enddo
    ! Same J differences as physical 410/480/880; independently recompute
    ! global oracle from indices. i=1,2 are byte-only NaN/signed-zero cases.
    lo=max(merge(2,3,north),smpi_jlo)
    hi=min(ny-1,smpi_jhi)
    do j=lo,hi
      k=merge(j+1,j-1,north)
      do i=3,im-1
        ! Mutant deliberately leaves qNaN poison. Count it by bits without
        ! introducing an invalid trap or allowing NaN comparisons to pass.
        if(transfer(a(i,k),0)==poison)then
          pbad=pbad+1
        else
          div=a(i,k)-a(i,j)
          ref=transfer(bits(i,k,epoch),0.)-transfer(bits(i,j,epoch),0.)
          if(.not.north)then
            div=a(i,j)-a(i,k)
            ref=transfer(bits(i,j,epoch),0.)-transfer(bits(i,k,epoch),0.)
          endif
          if(transfer(div,0)/=transfer(ref,0))pbad=pbad+1
        endif
      enddo
    enddo
  enddo
  if(smpi_size/=nr.or.witness/=nt)bad=1
  call smpi_index_raw(bad);call smpi_index_raw(pbad)
  k=smpi_rank
  call smpi_finalize()
  if(bad/=0)stop 40
  if(pbad/=0)then
    if(k==0.and.wrong)print *,'DETECTED_WRONG_DIRECTION',pbad
    stop 41
  endif
  if(wrong)stop 42
  if(k==0)print *,'PASS_FLUX_DIR generations/ranks/threads/JM=',3,nr,nt,ny
contains
  integer function bits(ii,jj,e)
    integer,intent(in)::ii,jj,e
    real::v
    if(ii==1)then
      bits=int(z'7FC10000')+e*256+jj
    else if(ii==2)then
      bits=0
      if(mod(e+jj,2)==0)bits=ibset(0,31)
    else
      v=real(e*10000+jj*jj*7+ii)
      bits=transfer(v,0)
    endif
  end function
end program
