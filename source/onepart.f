      Subroutine onepart(xstart,ystart,zstart,
     1    idum,deltat,xend,yend,zend,iout)     
      include 'comdeck'

      real KHj0k0,KHj1k0,KHk0,KHj0k1,KHj1k1,KHk1,KHp
      common/uvwpmax/upmax,vpmax,wpmax,deltaxmax,deltaymax,deltazmax
     .              ,udpmax,vdpmax,wdpmax

!
      HScNu=HPRNU
      VScNu=VPRNU

!     absolute coordinates in tranformed grid
!     for cell i:  x=[i,i+1) left close '[', right open ')' 
!     for cell j:  y=[j,j+1) left close '[', right open ')'
!     for cell k:  z=[z(k),z(k+1)] 

!     find cell indieces
      ic= xstart
      jc= ystart

      do k=1,kbm1
      if(zstart.le.z(k).and.zstart.ge.z(k+1)) goto 1
      enddo
      WRITE(*,*) 'trouble !',zstart,INTX
 1    kc=k

!     coordinates of cell center
      xc= ic+0.5
      yc= jc+0.5
      zc= zz(kc)


!---> interpolate u  
!     find i-range for u 
      i0u   = ic
      i1u   = ic+1
      xlocal= xstart-xc
      xdist = abs(xlocal)

!     find j-range for u
      j0u   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1u=j0u+1
      if(ylocal.le.0.) j1u=j0u-1

!     free-slip on land boundary
      if(fsm(ic,j1u).eq.0.) j1u=j0u

!     find k-range for u
      k0u   = kc
      zlocal= zstart-zc
      zdist = abs(zlocal)
      if(zlocal.gt.0.) k1u=k0u-1
      if(zlocal.le.0.) k1u=k0u+1

!     free-slip on land boundary
      if(k1u.lt.1) k1u=k0u
      if(k1u.gt.kbm1)k1u=k0u

!     interpolate u in k0 plane
      ui0k0=u(i0u,j0u,k0u)+ydist*(u(i0u,j1u,k0u)-u(i0u,j0u,k0u))
      ui1k0=u(i1u,j0u,k0u)+ydist*(u(i1u,j1u,k0u)-u(i1u,j0u,k0u))
      uk0  =ui0k0*(0.5-xlocal)+ui1k0*(0.5+xlocal)

!     interpolate u in k1 plane
      ui0k1=u(i0u,j0u,k1u)+ydist*(u(i0u,j1u,k1u)-u(i0u,j0u,k1u))
      ui1k1=u(i1u,j0u,k1u)+ydist*(u(i1u,j1u,k1u)-u(i1u,j0u,k1u))
      uk1  =ui0k1*(0.5-xlocal)+ui1k1*(0.5+xlocal)

!     interpolate u in z direction
      up=uk0 + (zdist/abs(0.5*dz(k1u)+0.5*dz(k0u))) * (uk1-uk0)

!---> interpolate v

!     find i-range for v
      i0v   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1v=i0v+1
      if(xlocal.le.0.) i1v=i0v-1

!     free-slip on land
      if(fsm(i1v,jc).eq.0.) i1v=i0v

!     find j-range for v
      j0v   = jc
      j1v   = j0v+1
      ylocal= ystart-yc
      ydist = abs(ylocal)

!     find k-range for v
      k0v   = kc
      zlocal= zstart-zc
      zdist = abs(zlocal)
      if(zlocal.gt.0.) k1v=k0v-1
      if(zlocal.le.0.) k1v=k0v+1

!     free-slip on land boundary
      if(k1v.lt.1) k1v=k0v
      if(k1v.gt.kbm1)k1v=k0v

!     interpolate v in k0 plane
      vj0k0=v(i0v,j0v,k0v)+xdist*(v(i1v,j0v,k0v)-v(i0v,j0v,k0v))
      vj1k0=v(i0v,j1v,k0v)+xdist*(v(i1v,j1v,k0v)-v(i0v,j1v,k0v))
      vk0  =vj0k0*(0.5-ylocal)+vj1k0*(0.5+ylocal)

!     interpolate v in k1 plane
      vj0k1=v(i0v,j0v,k1v)+xdist*(v(i1v,j0v,k1v)-v(i0v,j0v,k1v))
      vj1k1=v(i0v,j1v,k1v)+xdist*(v(i1v,j1v,k1v)-v(i0v,j1v,k1v))
      vk1  =vj0k1*(0.5-ylocal)+vj1k1*(0.5+ylocal)

!     interpolate v in z direction
      vp=vk0 + (zdist/abs(0.5*dz(k1v)+0.5*dz(k0v))) * (vk1-vk0)


!---> interpolate w

!     find i-range for w
      i0w   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1w=i0w+1
      if(xlocal.le.0.) i1w=i0w-1

!     free-slip on land
      if(fsm(i1w,jc).eq.0.)i1w=i0w

!     find j-range for w
      j0w   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1w=j0w+1
      if(ylocal.le.0.) j1w=j0w-1

!     free-slip on land
      if(fsm(ic,j1w).eq.0.)j1w=j0w

!xyz  this is an important change    
      if(fsm(i1w,j1w).eq.0.)then
      i1w=i0w
      j1w=j0w
      endif

!     find k-range for w
      k0w   = kc
      k1w   = kc+1
      zlocal= zstart-zc
      zdist = abs(zlocal)

!     interpolate w in k0 plane
      wj0k0= w(i0w,j0w,k0w) + xdist * (w(i1w,j0w,k0w)-w(i0w,j0w,k0w))
      wj1k0= w(i0w,j1w,k0w) + xdist * (w(i1w,j1w,k0w)-w(i0w,j1w,k0w))
      wk0  = wj0k0          + ydist * (wj1k0-wj0k0)

!     interpolate w in k1 plane
      wj0k1= w(i0w,j0w,k1w) + xdist * (w(i1w,j0w,k1w)-w(i0w,j0w,k1w))
      wj1k1= w(i0w,j1w,k1w) + xdist * (w(i1w,j1w,k1w)-w(i0w,j1w,k1w))
      wk1  = wj0k1          + ydist * (wj1k1-wj0k1)

!     interpolate w in z direction
      wp=wk0*(0.5+zlocal/dz(kc)) + wk1*(0.5-zlocal/dz(kc))



!---> interpolate h1 and h2

!     find i-range for h1 and h2
      i0h   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1h=i0h+1
      if(xlocal.le.0.) i1h=i0h-1

!     use value at one element
      if(fsm(i1h,jc).eq.0.)i1h=i0h

!     find j-range for h1 and h2
      j0h   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1h=j0h+1
      if(ylocal.le.0.) j1h=j0h-1

!     use value at one element
      if(fsm(ic,j1h).eq.0.)j1h=j0h

!xyz  this is an important change    
      if(fsm(i1h,j1h).eq.0.)then
      i1h=i0h
      j1h=j0h
      endif

!     interpolate h1
      h1j0   = h1(i0h,j0h) + xdist * (h1(i1h,j0h)-h1(i0h,j0h))
      h1j1   = h1(i0h,j1h) + xdist * (h1(i1h,j1h)-h1(i0h,j1h))
      h1pxyz = h1j0        + ydist * (h1j1-h1j0)

!     interpolate h2
      h2j0   = h2(i0h,j0h) + xdist * (h2(i1h,j0h)-h2(i0h,j0h))
      h2j1   = h2(i0h,j1h) + xdist * (h2(i1h,j1h)-h2(i0h,j1h))
      h2pxyz = h2j0        + ydist * (h2j1-h2j0)

!     interpolate D
      Dj0    = D(i0h,j0h) + xdist * (D(i1h,j0h)-D(i0h,j0h))
      Dj1    = D(i0h,j1h) + xdist * (D(i1h,j1h)-D(i0h,j1h))
      Dpxyz  = Dj0        + ydist * (Dj1-Dj0)
      if(Dpxyz.le.0) then
      WRITE(*,*)Dpxyz,xdist,ydist
      WRITE(*,*)i0h,i1h,j0h,j1h
      WRITE(*,*)D(i0h,j0h),D(i1h,j0h),D(i0h,j1h),D(i1h,j1h)
      endif



!     find k-range for AAM, KH
      k0h   = kc
      zlocal= zstart-zc
      zdist = abs(zlocal)
      if(zlocal.gt.0.) k1h=k0h-1
      if(zlocal.le.0.) k1h=k0h+1

!     free-slip on land boundary
      if(k1h.lt.1) k1h=k0h
      if(k1h.gt.kbm1)k1h=k0h

!     interpolate AAM in k0 plane
      AAMj0k0=AAM(i0h,j0h,k0h)+xdist*(AAM(i1h,j0h,k0h)-AAM(i0h,j0h,k0h))
      AAMj1k0=AAM(i0h,j1h,k0h)+xdist*(AAM(i1h,j1h,k0h)-AAM(i0h,j1h,k0h))
      AAMk0  =AAMj0k0         +ydist*(AAMj1k0-AAMj0k0)

!     interpolate AAM in k1 plane
      AAMj0k1=AAM(i0h,j0h,k1h)+xdist*(AAM(i1h,j0h,k1h)-AAM(i0h,j0h,k1h))
      AAMj1k1=AAM(i0h,j1h,k1h)+xdist*(AAM(i1h,j1h,k1h)-AAM(i0h,j1h,k1h))
      AAMk1  =AAMj0k1         +ydist*(AAMj1k1-AAMj0k1)

!     interpolate AAM in z direction
      AAMp=AAMk0 + (zdist/abs(0.5*dz(k1h)+0.5*dz(k0h))) * (AAMk1-AAMk0)


!     interpolate KH in k0 plane
      KHj0k0=KH(i0h,j0h,k0h)+xdist*(KH(i1h,j0h,k0h)-KH(i0h,j0h,k0h))
      KHj1k0=KH(i0h,j1h,k0h)+xdist*(KH(i1h,j1h,k0h)-KH(i0h,j1h,k0h))
      KHk0  =KHj0k0         +ydist*(KHj1k0-KHj0k0)

!     interpolate KH in k1 plane
      KHj0k1=KH(i0h,j0h,k1h)+xdist*(KH(i1h,j0h,k1h)-KH(i0h,j0h,k1h))
      KHj1k1=KH(i0h,j1h,k1h)+xdist*(KH(i1h,j1h,k1h)-KH(i0h,j1h,k1h))
      KHk1  =KHj0k1         +ydist*(KHj1k1-KHj0k1)

!     interpolate KH in z dirction
      KHp=KHk0 + (zdist/abs(0.5*dz(k1h)+0.5*dz(k0h))) * (KHk1-KHk0)



!---> interpolate UD and VD

!     interpolate UD at k0 plane
      UDj0k0= h2(i1h,j0h)/h1(i1h,j0h)*AAM(i1h,j0h,k0h)*D(i1h,j0h)-
     .        h2(i0h,j0h)/h1(i0h,j0h)*AAM(i0h,j0h,k0h)*D(i0h,j0h)
      UDj1k0= h2(i1h,j1h)/h1(i1h,j1h)*AAM(i1h,j1h,k0h)*D(i1h,j1h)-
     .        h2(i0h,j1h)/h1(i0h,j1h)*AAM(i0h,j1h,k0h)*D(i0h,j1h)
      UDk0  = UDj0k0        + ydist * (UDj1k0-UDj0k0)

!     interpolate UD at k1 plane
      UDj0k1= h2(i1h,j0h)/h1(i1h,j0h)*AAM(i1h,j0h,k1h)*D(i1h,j0h)-
     .        h2(i0h,j0h)/h1(i0h,j0h)*AAM(i0h,j0h,k1h)*D(i0h,j0h)
      UDj1k1= h2(i1h,j1h)/h1(i1h,j1h)*AAM(i1h,j1h,k1h)*D(i1h,j1h)-
     .        h2(i0h,j1h)/h1(i0h,j1h)*AAM(i0h,j1h,k1h)*D(i0h,j1h)
      UDk1  = UDj0k1        + ydist * (UDj1k1-UDj0k1)

!     interpolate UD in z-direction
      UDp   =UDk0 +(zdist/abs(0.5*dz(k1h)+0.5*dz(k0h)))*(UDk1-UDk0)
      UDp   =UDp/HScNu / (h1pxyz*h2pxyz*Dpxyz)
      if(xlocal.le.0.0) UDp= -UDp


!     interpolate VD at k0 plane
      VDi0k0= h1(i0h,j1h)/h2(i0h,j1h)*AAM(i0h,j1h,k0h)*D(i0h,j1h)-
     .        h1(i0h,j0h)/h2(i0h,j0h)*AAM(i0h,j0h,k0h)*D(i0h,j0h)
      VDi1k0= h1(i1h,j1h)/h2(i1h,j1h)*AAM(i1h,j1h,k0h)*D(i1h,j1h)-
     .        h1(i1h,j0h)/h2(i1h,j0h)*AAM(i1h,j0h,k0h)*D(i1h,j0h)
      VDk0  = VDi0k0        + xdist * (VDi1k0-VDi0k0)

!     interpolate VD at k1 plane
      VDi0k1= h1(i0h,j1h)/h2(i0h,j1h)*AAM(i0h,j1h,k1h)*D(i0h,j1h)-
     .        h1(i0h,j0h)/h2(i0h,j0h)*AAM(i0h,j0h,k1h)*D(i0h,j0h)
      VDi1k1= h1(i1h,j1h)/h2(i1h,j1h)*AAM(i1h,j1h,k1h)*D(i1h,j1h)-
     .        h1(i1h,j0h)/h2(i1h,j0h)*AAM(i1h,j0h,k1h)*D(i1h,j0h)
      VDk1  = VDi0k1        + xdist * (VDi1k1-VDi0k1)

!     interpolate VD in z-direction
      VDp   =VDk0 + (zdist/abs(0.5*dz(k1h)+0.5*dz(k0h)))*(VDk1-VDk0)
      VDp   =VDp/HScNu / (h1pxyz*h2pxyz*Dpxyz)
      if(ylocal.le.0.0) VDp= -VDp

!     interpolate WD at k0 plane
      WDp   =(KHk1-KHK0)/abs(0.5*dz(k1h)+0.5*dz(k0h))
      WDp   =WDp/VScNu / (Dpxyz*Dpxyz) 
      if(zlocal.le.0.0) WDp= -WDp

      if(abs(up).gt.abs(upmax))upmax=up
      if(abs(vp).gt.abs(vpmax))vpmax=vp
      if(abs(wp).gt.abs(wpmax))wpmax=wp
      if(abs(udp).gt.abs(udpmax))udpmax=udp
      if(abs(vdp).gt.abs(vdpmax))vdpmax=vdp
      if(abs(wdp).gt.abs(wdpmax))wdpmax=wdp

         
!---> advection + random walk + pseudo velocity
      deltax=deltat * (Up/h1pxyz)
     .   +gasdev(idum)*sqrt(2.*deltat*(AAMp/HScNu)/(h1pxyz*h1pxyz))
     .   +deltat*UDp 
      deltay=deltat*Vp/h2pxyz
     .   +gasdev(idum)*sqrt(2.*deltat*(AAMp/HScNu)/(h2pxyz*h2pxyz))
     .   +deltat*VDp 
      deltaz=deltat*Wp/ Dpxyz
     .   +gasdev(idum)*sqrt(2.*deltat*( KHp/VScNu)/( Dpxyz* Dpxyz))
     .   +deltat*WDp 

      if(abs(deltax).gt.abs(deltaxmax))deltaxmax=deltax
      if(abs(deltay).gt.abs(deltaymax))deltaymax=deltay
      if(abs(deltaz).gt.abs(deltazmax))deltazmax=deltaz


!---> update particle location
      xend=xstart + deltax
      yend=ystart + deltay
      zend=zstart + deltaz

!     indices for new cell
      icnew = xend
      jcnew = yend
      do k=1,kbm1
      if(zend.le.z(k).and.zend.ge.z(k+1)) goto 2
      enddo
 2    kcnew=k



      if(abs(ic-icnew).gt.1.or.abs(jc-jcnew).gt.1) then
      WRITE(*,*)INTX,ic,icnew,jc,jcnew
      WRITE(*,*)deltax,deltay,deltaz
      endif

!     check if it is on open boundary
      do loop=1,NUMEBC
         if(icnew.eq.ieta(loop).and.jcnew.eq.jeta(loop)) then
            iout=1
            WRITE(*,*)'i am out from open boundary'
            return
         endif
      enddo

!---> reflection from landboundary, surface and bottom
      ireflect=1
      if(ireflect.eq.0) goto 100

!     center coordinates for new cell 
      xcnew = icnew+0.5
      ycnew = jcnew+0.5
      zcnew = zz(kcnew)

!---> new position is on land
      if(fsm(icnew,jcnew).eq.0.0) then

!---> reflection from +x direction
      if(deltax.gt.0.) then
      do loop=ic+1,icnew
      if(fsm(loop,jc).eq.0.0) goto 10
      enddo
      goto 11
!hli10    icwall =loop
!hli      xwall  =icwall
!hli      xtowall=xwall-xstart
!hli      xbounce=deltax-xtowall
!hli      if(xbounce.ge.0.) xend=xwall-xbounce
!hli      endif
!hli
!hli  keep half distance between previous location and wall (06/04/98)
!hli
10    icwall =loop
      xwall  =icwall
      xtowall=xwall-xstart
      xend=xstart+xtowall/2.0
      endif
11    continue

!---> reflection from -x direction
      if(deltax.lt.0.) then
      do loop=ic-1,icnew,-1
      if(fsm(loop,jc).eq.0.0) goto 20
      enddo
      goto 21
!hli20    icwall =loop
!hli      xwall  =icwall+1
!hli      xtowall=xstart-xwall
!hli      xbounce=abs(deltax)-xtowall
!hli      if(xbounce.ge.0.) xend=xwall+xbounce
!hli      endif
!hli
!hli   keep half distance between previous location and wall (06/04/98)
!hli
20    icwall =loop
      xwall  =icwall
      xtowall=xstart-xwall
      xend=xstart-xtowall/2.0
      endif
21    continue


!---> reflection from +y direction
      if(deltay.gt.0.) then
      do loop=jc+1,jcnew
      if(fsm(ic,loop).eq.0.0) goto 30
      enddo
      goto 31
!hli30    jcwall =loop
!hli      ywall  =jcwall
!hli      ytowall=ywall-ystart
!hli      ybounce=deltay-ytowall
!hli      if(ybounce.gt.0.) yend=ywall-ybounce
!hli      endif
!hli
!hli  keep half distance between previous location and wall (06/04/98)
30    jcwall =loop
      ywall  =jcwall
      ytowall=ywall-ystart
      yend=ystart+ytowall/2.0
      endif
!hli
31    continue

!---> reflection from -y direction
      if(deltay.lt.0.) then
      do loop=jc-1,jcnew,-1
      if(fsm(ic,loop).eq.0.0) goto 40
      enddo
      goto 41
!hli40    jcwall =loop
!hli      ywall  =jcwall+1
!hli      ytowall=ystart-ywall
!hli      ybounce=abs(deltay)-ytowall
!hli      if(ybounce.ge.0.) yend=ywall+ybounce
!hli      endif
!hli
!hli  keep half distance between previous location and wall (06/04/98)
!hli
40    jcwall =loop
      ywall  =jcwall
      ytowall=ystart-ywall
      yend=ystart-ytowall/2.0
      endif
!hli
41    continue

      endif
!     end of horizontal reflection 


!---> reflection from surface
      if(zend.gt.0)then
!hli      zwall  =0.
!hli      ztowall=zwall-zstart
!hli      zbounce=deltaz-ztowall
!hli      zend   =zwall-zbounce
!hli      endif
!hli
!hli  keep half distance between previous location and surface (06/04/98)
      zwall  =0.
      ztowall=zwall-zstart
      zend   =zstart+ztowall/2.0
      endif
!hli

!---> reflection from bottom
      if(zend.lt.-1.0)then
!hli      zwall  =-1.0
!hli      ztowall=zstart-zwall
!hli      zbounce=abs(deltaz)-ztowall
!hli      zend   =zwall+zbounce
!hli      endif
!hli
!hli  keep half distance between previous location and bottom (06/04/98)
      zwall  =-1.0
      ztowall=zstart-zwall
      zend   =zstart-ztowall/2.0
      endif
!hli
   
!---> 'big incremental distance' due to advection and/or
!     diffusion might bounce the particle out of domain. 
!     Force particle to stay if still on land.

!     also force it to stay if it on edge of the outmost cell
      istay=1
      if(istay.eq.1)then
      icnewest = xend
      jcnewest = yend
      if(fsm(icnewest,jcnewest).eq.0.0) then
      xend=xstart
      yend=ystart
      zend=zstart
      endif
      if(zend.le.-1.0.or.zend.gt.0.0) then
      xend=xstart
      yend=ystart
      zend=zstart
      endif
      endif

100   return
      end
