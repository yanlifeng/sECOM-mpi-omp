program check_runtime
  use omp_lib
  use mpi_f08,only:MPI_Query_thread,MPI_THREAD_SERIALIZED
  use iso_c_binding
  implicit none
  include 'mpideck'
  interface
    integer(c_int) function place(rank,lane,nt) bind(c)
      import c_int
      integer(c_int),value::rank,lane,nt
    end function
    integer(c_int) function getpid() bind(c)
      import c_int
    end function
  end interface
  integer,parameter::im=5,maxj=13
  integer::jm,i,j,rank,n,l,h,idx,nt,want_nt,provided,rc,lane,cp(8),nr,trial
  integer::iv(maxj),ie(maxj),flags(im)
  real::a(im,maxj),v(maxj),e(maxj),x
  real::et(im,maxj),etf(im,maxj),va(im,maxj),vaf(im,maxj),mask(im,maxj),etb(im,maxj)
  character(16)::arg
  call get_command_argument(1,arg);read(arg,*) jm
  call get_command_argument(2,arg);read(arg,*) nr
  call get_command_argument(3,arg);read(arg,*) want_nt
  if(smpi_size/=1.or.smpi_rank/=0.or.smpi_first(3)/=3.or.smpi_last(13)/=13)stop 60
  if(smpi_localbc)stop 75
  smpi_localbc=.true.
  call smpi_init(im,jm)
  if(smpi_localbc)stop 76
  call MPI_Query_thread(provided,rc)
  if(rc/=0.or.provided<MPI_THREAD_SERIALIZED.or.smpi_size/=nr.or.jm>maxj)stop 61
  rank=smpi_rank;n=smpi_size;l=smpi_jlo;h=smpi_jhi
  if(l/=1+rank*(jm/n)+min(rank,mod(jm,n)))stop 62
  if(h/=l-1+(jm/n)+merge(1,0,rank<mod(jm,n)))stop 63
  do i=1,4
    do j=1,jm
      if((j>=smpi_first(2).and.j<=smpi_last(jm)).neqv.(j>=max(2,l).and.j<=h))stop 64
      if(smpi_owns(j).neqv.(j>=l.and.j<=h))stop 65
      if((smpi_owner_raw(j)==rank).neqv.smpi_owns(j))stop 77
      if((j>=smpi_color(2,jm,i).and.j<=h.and.mod(j-2,i)==0).neqv. &
          (j>=max(2,l).and.j<=h.and.mod(j-2,i)==0))stop 66
    enddo
  enddo
  a=-999;v=-999;e=-999;iv=-999;ie=-999
  et=11.;etf=12.;va=13.;vaf=14.;mask=15.;etb=16.
  flags=[0,1,2,3,0]
  cp=-1;nt=0
!$omp parallel default(shared) private(j,i,lane,trial)
  lane=omp_get_thread_num()
  if(lane>=8)stop 67
  cp(1+lane)=place(rank,lane,want_nt)
!$omp master
  nt=omp_get_num_threads()
!$omp end master
  do trial=1,3
!$omp do schedule(static)
    do j=l,h
      do i=1,im
        a(i,j)=trial*1000+10*j+i
      enddo
      v(j)=j+0.5;e(j)=-j-0.25;iv(j)=j;ie(j)=-j
    enddo
!$omp end do nowait
    call smpi_halo(a,im,jm)
!$omp single
    do j=max(1,l-1),min(jm,h+1)
      if(any(a(:,j)/=[(real(trial*1000+10*j+i),i=1,im)]))stop 68
    enddo
!$omp end single
    call smpi_all(a,im,jm)
    call smpi_rows(v,e,jm)
    call smpi_irows(iv,ie,jm)
  enddo
!$omp single
  x=0.;idx=0
  if(smpi_owns(jm))then
    x=transfer(ibset(0,31),1.)
    idx=123
  endif
  call smpi_rbcast_row_raw(x,jm)
  call smpi_index_raw(idx)
!$omp end single
  call smpi_wad_publish(flags,et,etf,va,vaf,mask,etb,im,jm)
!$omp end parallel
  if(nt/=want_nt.or.any(cp(1:nt)<0).or.transfer(x,0)/=ibset(0,31).or.idx/=123)stop 69
  do j=1,jm
    if(any(a(:,j)/=[(real(3000+10*j+i),i=1,im)]).or.v(j)/=j+0.5.or.e(j)/=-j-0.25.or.iv(j)/=j.or.ie(j)/=-j)stop 70
    do i=1,im
      if(j==l.and.rank>0.and.btest(flags(i),0))then
        if(va(i,j)/=0..or.vaf(i,j)/=0..or.mask(i,j)/=0.)stop 71
      else
        if(va(i,j)/=13..or.vaf(i,j)/=14..or.mask(i,j)/=15.)stop 72
      endif
      if(j==l.and.rank>0.and.btest(flags(i),1))then
        if(et(i,j)/=16..or.etf(i,j)/=16.)stop 73
      else
        if(et(i,j)/=11..or.etf(i,j)/=12.)stop 74
      endif
    enddo
  enddo
  write(*,*) 'SMPI_PASS rank,size,pid,lo,hi,omp,threadlevel,cpus:',rank,n,getpid(),l,h,nt,provided,cp(1:nt)
  smpi_localbc=.true.
  call smpi_finalize()
  if(smpi_localbc)stop 78
end program
