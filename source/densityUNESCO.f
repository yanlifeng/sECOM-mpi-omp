!NG      SUBROUTINE densityUNESCO
      SUBROUTINE dens
!
!*************************************************************************
!     ECOMSED MODEL 
!     VERSION 1.3
!     FEBRUARY 2002
!
*               Copyright (C) 2003, Stevens Institute of Technology
*                                                                        *
*     Author: Dov Kruger
!*************************************************************************
      INCLUDE 'comdeck'
!      SAVE
!
!         COMPUTE DENSITY-1.0 and store into the array RHO
!
      REAL*8 RHOF(IM,JM,KB),TF(IM,JM,KB),SF(IM,JM,KB)
      EQUIVALENCE (A,RHOF),(VH,TF),(UF,SF)

      real*8 k00, k11, k12, k13, k14, k15
      real*8 k21, k22, k23, k24, k25
      real*8 k31, k32, k33

!NG   !NGafterPOMfix old:k00 = 999.842594 - 1000
      k00 = -0.157406 ! avoid 32-bit rounding error
      k11 = 6.793952E-2
!NG   !NGafterPOM old:k12 = -9.09529E-3
      k12 = -9.095290E-3
      k13 = 1.001685E-4
      k14 = -1.120083E-6
      k15 = 6.536332E-9
      
      k21 = 0.824493
      k22 = -4.0899E-3
      k23 = 7.6438E-5
      k24 = -8.2467E-7
      k25 = 5.3875E-9
      
      k31 = -5.72466E-3
      k32 = 1.0227E-4
      k33 = -1.6546E-6
      
      k41 = 4.8314E-4

      do k=1,kbm1
         do j=1,jm
            do i=1,im
               s(i,j,k)=amax1(0.,s(i,j,k))
!NG             if(s(i,j,k).lt.0.0)then
!NG              write(*,*)'negative salinity',
!NG     $    time,i,j,k,s(i,j,k)
!NG              s(i,j,k)=0.0
!NG             end if
               s2 = s(i,j,k) * s(i,j,k)
               rho(i,j,k) = k00 +
     $              ((((k15*t(i,j,k) + k14)*t(i,j,k) + k13)*t(i,j,k) +
     $              k12) * t(i,j,k) + k11) * t(i,j,k) +
     $              ((((k25*t(i,j,k) + k24)*t(i,j,k) + k23)*t(i,j,k) +
     $              k22) * t(i,j,k) + k21) * s(i,j,k) +
     $              ((k33*t(i,j,k) + k32)*t(i,j,k) + k31)
     $              * sqrt(s2 * s(i,j,k)) +
     $              k41 * s2 
              rho(i,j,k)= rho(i,j,k) * 1.0E-03 * fsm(i,j)
            enddo
         enddo
      enddo
      RETURN
      END
