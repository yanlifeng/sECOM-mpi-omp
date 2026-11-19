! Local profile-accounting fixture; not an MPI program or an elapsed benchmark.
program check_timers
  use omp_lib
  implicit none
  integer::i,j
  call sd_stage(1)
  call sd_id(0,1)
  call sd_stage(2)
  call sd_site(4,1);call sd_site(7,1);call sd_site(7,0);call sd_site(4,0)
  call sd_stage(3)
  do i=1,3
    call sd_step()
    do j=1,3
      call sd_site(j,1);call sd_site(j,0)
    enddo
    call sd_site(5,1)
!$omp parallel default(shared) private(j)
    do j=1,2
      call sd_site(10,1)
!$omp single
      call sd_mp(1,1);call sd_mp(1,0)
!$omp end single
      call sd_site(10,0)
      call sd_site(11,1)
!$omp single
      call sd_mp(7,1);call sd_mp(7,0)
!$omp end single
      call sd_site(11,0)
      call sd_site(12,1)
!$omp single
      call sd_mp(5,1);call sd_mp(5,0)
!$omp end single
      call sd_site(12,0)
    enddo
!$omp end parallel
    call sd_site(5,0)
  enddo
  call sd_stage(4)
  call sd_site(6,1);call sd_site(6,0)
  call sd_stage(5)
  call sd_dump()
end program
