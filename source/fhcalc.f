      SUBROUTINE FHCALC
!
!  CALC. FETCH AND AVERAGE DEPTH FOR WIND WAVE MODEL
!
!  REVISION DATE:  OCTOBER 11, 1995
!
!************************************************************
!
      INCLUDE 'comdeck'
!
      REAL UTEMP(IM,JM),VTEMP(IM,JM),FTEMP(36,IM,JM)
!
!  STORE ORIGINAL VEL. VALUES (IN CASE OF HOT START)
!
      DO 5 J=1,JM
        DO 5 I=1,IM
          UTEMP(I,J)=U(I,J,1)
          VTEMP(I,J)=V(I,J,1)
 5    CONTINUE
!
!  OPEN CORNER LOCATIONS FILE
!
!  CORNER LOCATION CONVENTION:  XWCOR(I,J) = x(i-1/2, j-1/2) 
!     (LOWER LEFT-HAND CORNER)  YWCOR(I,J) = y(i-1/2, j-1/2)
!
        OPEN (UNIT=33,FILE='corner_loc',FORM='FORMATTED',STATUS='OLD') !hli
!
        DO 300 N=1,1000000
          READ (33,*,END=310)I,J,XWCOR(I,J),YWCOR(I,J)
 300    CONTINUE
!
!  CALC. FETCH AND AVERAGE DEPTH FOR EACH OF 36 DIRECTIONS (10 deg SEPARATION)
!
 310  DO 10 N=1,36
!
        WANG=2.*3.141593*FLOAT(N)/36.
!
!  SET VELOCITIES FOR THIS DIRECTION
!
        DO 20 J=2,JMM1
          DO 20 I=2,IMM1
            DANG=2.5*3.141593-ANG(I,J)
!
            PHIANG=DANG-WANG
!
            U(I,J,1)=COS(PHIANG)
            V(I,J,1)=SIN(PHIANG)
 20     CONTINUE
!
        DO 21 I=2,IM
          U(I,JM,1)=U(I,JMM1,1)
          V(I,JM,1)=V(I,JMM1,1)
 21     CONTINUE
!
        DO 22 J=2,JM
          U(IM,J,1)=U(IMM1,J,1)
          V(IM,J,1)=V(IMM1,J,1)
 22     CONTINUE
 
!
!  CALC. FETCH FOR EACH (i,j) USING THIS WIND DIRECTION
!
        DO 30 J=2,JMM1
          DO 30 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
!
!  INITIAL STARTING POINT
!
              XP1=FLOAT(I)+0.5
              YP1=FLOAT(J)+0.5
!
              IC=XP1
              JC=YP1
!
              NSUM=1
              HSUM=H(IC,JC)
!
!  IOUT = 0 :  INSIDE DOMAIN
!  IOUT = 1 :  OUTSIDE DOMAIN
!
              IOUT=0
!
 40           XSTART=XP1
              YSTART=YP1
!
              CALL SONEPART(XSTART,YSTART,XP1,YP1,IOUT)
              
!
!  DETERMINE (i,j) LOCAL GRID CELL
!
               IC=XP1
               JC=YP1
!
!  SUM DEPTHS ALONG FETCH
!
               IF (IOUT.EQ.0) THEN
                 HSUM=HSUM+H(IC,JC)
                 NSUM=NSUM+1
                 GOTO 40
               ENDIF
!
!  END OF FETCH HAS BEEN REACHED
!
!  DETERMINE LENGTH OF FETCH
!
              IC=XP1    
              JC=YP1  
!
              XC=IC+0.5 
              YC=JC+0.5 
!
              XLOCAL=XP1-XC
              YLOCAL=YP1-YC
!      
!  CALC. LOCATION OF END OF FETCH
!
              XFIN=xwcor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +xwcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +xwcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +xwcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
!
              YFIN=ywcor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +ywcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +ywcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +ywcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
!
!  CALC. LOCATION OF CENTER OF ELEMENT (i,j)
!
              XSTRT=0.25*(XWCOR(I,J)+XWCOR(I+1,J)+XWCOR(I,J+1)+
     +                                          XWCOR(I+1,J+1))
!
              YSTRT=0.25*(YWCOR(I,J)+YWCOR(I+1,J)+YWCOR(I,J+1)+
     +                                          YWCOR(I+1,J+1))
!
!  CALC. LENGTH OF FETCH
!
              FETCH(N,I,J)=SQRT((XFIN-XSTRT)**2.+(YFIN-YSTRT)**2.)
!
!  CALC. MEAN DEPTH ALONG FETCH
!
              HMEAN(N,I,J)=HSUM/FLOAT(NSUM)
            ENDIF
 30     CONTINUE
 10   CONTINUE
!
!  RESET ORIGINAL VEL. VALUES (IN CASE OF HOT START)
!
      DO 105 J=1,JM
        DO 105 I=1,IM
          U(I,J,1)=UTEMP(I,J)
          V(I,J,1)=VTEMP(I,J)
 105  CONTINUE
!
!  MODIFY FETCH LENGTHS TO ACCOUNT FOR FETCH WIDTH EFFECTS
!
!  USE MODIFIED VERSION OF USACOE METHOD:  AVERAGE OVER 20 DEG ARC
!  INSTEAD OF 24 DEG ARC (3 RADIALS +/- 10 DEG INSTEAD OF 9 RADIALS
!  AT +/- 12 DEG)
!
      DO 200 N=1,36
        DO 200 J=2,JMM1
          DO 200 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
              IF (N.EQ.1.OR.N.EQ.36) THEN
                IF (N.EQ.1) THEN
                  FTEMP(1,I,J)=(FETCH(36,I,J)+FETCH(1,I,J)+
     +                              FETCH(2,I,J))/3.
                ELSE
                  FTEMP(36,I,J)=(FETCH(35,I,J)+FETCH(36,I,J)+
     +                              FETCH(1,I,J))/3.
                ENDIF
              ELSE
                FTEMP(N,I,J)=(FETCH(N-1,I,J)+FETCH(N,I,J)+
     +                              FETCH(N+1,I,J))/3.
              ENDIF
            ENDIF
 200  CONTINUE
!
      DO 210 N=1,36
        DO 210 J=2,JMM1
          DO 210 I=2,IMM1
            FETCH(N,I,J)=FTEMP(N,I,J)
 210  CONTINUE
!
      RETURN
      END
!
!**************************************************************************
!
!  PARTICLE TRACKING SUBROUTINE
!
      Subroutine sonepart(xstart,ystart,xend,yend,iout)
!
      include 'comdeck'
!
!     absolute coordinates in tranformed grid
!     for cell i:  x=[i,i+1) left close '[', right open ')' 
!     for cell j:  y=[j,j+1) left close '[', right open ')'
!
!     find cell indieces
!
      ic= xstart
      jc= ystart
!
!     coordinates of cell center
!
      xc= ic+0.5
      yc= jc+0.5
!
      kc=1
!
!---> interpolate u  
!     find i-range for u 
!
      i0u   = ic
      i1u   = ic+1
      xlocal= xstart-xc
      xdist = abs(xlocal)
!
!     find j-range for u
!
      j0u   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1u=j0u+1
      if(ylocal.le.0.) j1u=j0u-1
!
!     free-slip on land boundary
!
      if(fsm(ic,j1u).eq.0.) j1u=j0u
!
!     find k-range for u
!
      k0u   = kc
!
!     interpolate u in k0 plane
!
      ui0k0=u(i0u,j0u,k0u)+ydist*(u(i0u,j1u,k0u)-u(i0u,j0u,k0u))
      ui1k0=u(i1u,j0u,k0u)+ydist*(u(i1u,j1u,k0u)-u(i1u,j0u,k0u))
      UP   =ui0k0*(0.5-xlocal)+ui1k0*(0.5+xlocal)
!
!---> interpolate v
!
!     find i-range for v
      i0v   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1v=i0v+1
      if(xlocal.le.0.) i1v=i0v-1
!
!     free-slip on land
!
      if(fsm(i1v,jc).eq.0.) i1v=i0v
!
!     find j-range for v
!
      j0v   = jc
      j1v   = j0v+1
      ylocal= ystart-yc
      ydist = abs(ylocal)
!
!     find k-range for v
!
      k0v   = kc
!
!     interpolate v in k0 plane
!
      vj0k0=v(i0v,j0v,k0v)+xdist*(v(i1v,j0v,k0v)-v(i0v,j0v,k0v))
      vj1k0=v(i0v,j1v,k0v)+xdist*(v(i1v,j1v,k0v)-v(i0v,j1v,k0v))
      VP   =vj0k0*(0.5-ylocal)+vj1k0*(0.5+ylocal)
!
!---> interpolate h1 and h2
!
!     find i-range for h1 and h2
!
      i0h   = ic
      xlocal= xstart-xc
      xdist = abs(xlocal)
      if(xlocal.gt.0.) i1h=i0h+1
      if(xlocal.le.0.) i1h=i0h-1
!
!     use value at one element
!
      if(fsm(i1h,jc).eq.0.)i1h=i0h
!
!     find j-range for h1 and h2
!
      j0h   = jc
      ylocal= ystart-yc
      ydist = abs(ylocal)
      if(ylocal.gt.0.) j1h=j0h+1
      if(ylocal.le.0.) j1h=j0h-1
!
!     use value at one element
!
      if(fsm(ic,j1h).eq.0.)j1h=j0h
!
!xyz  this is an important change    
!
      if(fsm(i1h,j1h).eq.0.)then
      i1h=i0h
      j1h=j0h
      endif
!
!     interpolate h1
!
      h1j0   = h1(i0h,j0h) + xdist * (h1(i1h,j0h)-h1(i0h,j0h))
      h1j1   = h1(i0h,j1h) + xdist * (h1(i1h,j1h)-h1(i0h,j1h))
      h1pxyz = h1j0        + ydist * (h1j1-h1j0)
!
!     interpolate h2
!
      h2j0   = h2(i0h,j0h) + xdist * (h2(i1h,j0h)-h2(i0h,j0h))
      h2j1   = h2(i0h,j1h) + xdist * (h2(i1h,j1h)-h2(i0h,j1h))
      h2pxyz = h2j0        + ydist * (h2j1-h2j0)
!
!---> advection
!
      deltax=DTI * (Up/h1pxyz)
!
      deltay=DTI*Vp/h2pxyz
!
!---> update particle location
!
      xend=xstart + deltax
      yend=ystart + deltay
!
!     indices for new cell
!
      icnew = xend
      jcnew = yend
!
      IF (XLOCAL.GT.0.49.AND.DUM(IC+1,JC).LE.0.0) IOUT=1
      IF (XLOCAL.LT.0.-49.AND.DUM(IC,JC).LE.0.0) IOUT=1
      IF (YLOCAL.GT.0.49.AND.DVM(IC,JC+1).LE.0.0) IOUT=1
      IF (YLOCAL.LT.-0.49.AND.DVM(IC,JC).LE.0.0) IOUT=1
!
      DXNEW=XEND-XC
!
      IF (DXNEW.GT.0.0) THEN
        DXNEW1=0.5-AMOD(DXNEW,0.5)
        DXNEW1=ABS(DXNEW1)
      ELSE
        DXNEW1=ABS(DXNEW)
        DXNEW1=AMOD(DXNEW1,0.5)-0.5
        DXNEW1=ABS(DXNEW1)
      ENDIF
!
      DYNEW=YEND-YC
!
      IF (DYNEW.GT.0.0) THEN
        DYNEW1=0.5-AMOD(DYNEW,0.5)
        DYNEW1=ABS(DYNEW1)
      ELSE
        DYNEW1=ABS(DYNEW)
        DYNEW1=AMOD(DYNEW1,0.5)-0.5
        DYNEW1=ABS(DYNEW1)
      ENDIF
!
      IF (DXNEW1.LE.ABS(DELTAX)) THEN
        IF (DXNEW.GT.0.0) THEN
          ICNEW=XC
          IF (DUM(ICNEW+1,JCNEW).LE.0.0) THEN
            IOUT=1
            XEND=XC+0.5
          ENDIF
        ELSE
          ICNEW=XC
          IF (DUM(ICNEW,JCNEW).LE.0.0) THEN
            IOUT=1
            XEND=XC-0.5
          ENDIF
        ENDIF
      ENDIF
!
      IF (DYNEW1.LE.ABS(DELTAY)) THEN
        IF (DYNEW.GT.0.0) THEN
          JCNEW=YC
          IF (DVM(ICNEW,JCNEW+1).LE.0.0) THEN
            IOUT=1
            YEND=YC+0.5
          ENDIF
        ELSE
          JCNEW=YC
          IF (DVM(ICNEW,JCNEW).LE.0.0) THEN
            IOUT=1
            YEND=YC-0.5
          ENDIF
        ENDIF
      ENDIF
!
      IF (FSM(ICNEW,JCNEW).LE.0.0) IOUT=1
      IF (FSM(ICNEW,JCNEW).LE.0.0) IOUT=1
!
      RETURN
      END
