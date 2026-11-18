      subroutine PARTRAK(idum,IRELEND)
!
      include 'comdeck'
!
      common/uvwpmax/upmax,vpmax,wpmax,deltaxmax,deltaymax,deltazmax
     .              ,udpmax,vdpmax,wdpmax
!
      if (INTX.lt.IRELST) RETURN
      intip=INTX-IRELST 
!
!---> calculate how many groups of particles in the system
!
      ngradeloop=intip/NFREQ + 1
      if (ngradeloop.ge.NCONV) ngradeloop=NCONV
 4321 FORMAT(I8,F10.4,3I5)
!     IF (INTX.GT.IRELEND)GOTO 31
!
!     if (mod(intip,NFREQ).eq.0.and.ngradeloop.gt.1) then     !nkim 061098
!
!   ANOTHER RELEASE OF GROUP OF PARTICLES
        if (mod(intip,NFREQ).ne.0.or.INTX.gt.irelend)go to 31
        if(ngradeloop.eq.1) go to 29     ! nkim 061098
!---> advance all groups to a higher grade before new group is introduced
        do 10 NN=ngradeloop,2,-1
          do 10 MM=1,NPART
            do 10 LL=1,nsource
              XP(LL,MM,NN)=   XP(LL,MM,NN-1)
              YP(LL,MM,NN)=   YP(LL,MM,NN-1)
              ZP(LL,MM,NN)=   ZP(LL,MM,NN-1)
              inout(LL,MM,NN)=inout(LL,MM,NN-1)
10      CONTINUE
!     endif                       ! nkim 061098
!
  29   continue
!---> introduce a new group every NFREQ timestep
        do 30 MM=1,NPART
          do 30 LL=1,nsource
            XP(LL,MM,1)= isource(LL)+0.5
            YP(LL,MM,1)= jsource(LL)+0.5
            ZP(LL,MM,1)= ZZ(ksource(LL))
            inout(LL,MM,1)  = 0
 30     CONTINUE
!
!---> advance all particles in all possible groups at the current time
!
 31   continue
      do 40 LL=1,nsource
        do 40 NN=1,ngradeloop
          do 40 MM=1,NPART
!
            if(inout(LL,MM,NN).eq.1) goto 40
!
            xstart=xp(LL,MM,NN)
            ystart=yp(LL,MM,NN)
            zstart=zp(LL,MM,NN)
            iout=inout(LL,MM,NN)
!
            call onepart(xstart,ystart,zstart,idum,dti,xend,yend,
     +                            zend,iout)
!
            xp(LL,MM,NN)=xend
            yp(LL,MM,NN)=yend
            zp(LL,MM,NN)=zend
            inout(LL,MM,NN)=iout
 40   CONTINUE
!
      return
      end
!
!**************************************************************************
!
      function gasdev(idum)
      data iset/0/
      if( iset.eq.0 )then
 1    continue
       v1=2.0*ran2(idum)-1.
       v2=2.0*ran2(idum)-1.
       r=v1*v1+v2*v2
       if(r.ge.1.0)goto 1
       fac=sqrt(-2.0*log(r)/r)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
      else
       gasdev=gset
       iset=0
      endif
!
      return
      end
!
!**************************************************************************
!
      function ran2(idum)
      parameter(m=714025, ia=1366, ic=150889, rm=1.0/m)
      dimension ir(97)
      data iff/0/
      if(idum.lt.0.or.iff.eq.0)then
         iff=1
         idum=mod(ia-idum,m)
         do 11 j=1,97
         idum=mod(ia*idum+ic,m)
         ir(j)=idum
 11      continue
         idum=mod(ia*idum+ic,m)
         iy=idum
      endif
         j=1+(97*iy)/m
         if(j.gt.97.or.j.lt.1)then
         print *, 'stop in ran2 ',j
         stop
         endif
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end
