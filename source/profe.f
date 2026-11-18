      SUBROUTINE PROFE
!     VERSION(03/15/03)
!          Mixing formulation of PRITCHARD(19??) and OFFICER (1976)

!     ======>  constants to be read in by user:  prit(im,jm)

      INCLUDE 'comdeck'
      SAVE
!o
      REAL Khomo(IM,JM,KB)
      DIMENSION ZDIST(IM,JM,KB),SHEAR(IM,JM,KB)
      DIMENSION RI(IM,JM,KB),buoy(IM,JM,KB),STAB(IM,JM,KB)
      DIMENSION DH(IM,JM),PRIT(IM,JM)
!
      gee=9.806
      small=1.e-10
!
      DO 10 J=2,JMM1
      DO 10 I=2,IMM1
      PRIT(I,J)=VPRNU    ! HARDWIRED
 10   DH(I,J)=H(I,J)+ETF(I,J)
!
      DO 105 J=2,JMM1
      DO 105 K=1,KB
      DO 105 I=2,IMM1
 105  ZDIST(I,J,K)=-DH(I,J)*DH(I,J)*Z(K)*(1.+Z(K))**2
!105  ZDIST(I,J,K)=((DH(I,J)*Z(K)+ETF(I,J))+(1.+Z(K)))
!
      DO 115 J=2,JMM1
      DO 115 K=2,KBM1
      DO 115 I=2,IMM1
      Buoy(I,J,K)=-GEE*(RHO(I,J,K-1)-RHO(I,J,K))/(DZZ(K-1)*DH(I,J))
!     next line is a crude assumption, we should 
!         mix the unstable region via "convective adjustment"
      if(buoy(i,j,k).lt.0.) buoy(i,j,k)=0.0
 115  CONTINUE
!
      DO 120 J=2,JMM1
      DO 120 K=2,KBM1
      DO 120 I=2,IMM1
      SHEAR(I,J,K)=.25*
     .       ( (U(I,J,K)-U(I,J,K-1)+U(I+1,J,K)-U(I+1,J,K-1))**2
     .         +(V(I,J,K)-V(I,J,K-1)+V(I,J+1,K)-V(I,J+1,K-1))**2)
     .              /(DZZ(K-1)*DH(I,J))**2
 120  RI(I,J,K)=Buoy(I,J,K)/(SHEAR(I,J,K)+small)
!
      DO 110 K=2,KBM1
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
 110  KHOMO(I,J,K)=0.41*0.41*ZDIST(I,J,K)*sqrt(SHEAR(I,J,K))
!110  KHOMO(I,J,K)=prit(I,J)*ZDIST(I,J,K)**2*sqrt(SHEAR(I,J,K))
!
      DO 280 J=2,JMM1
      DO 280 K=2,KBM1
      DO 280 I=2,IMM1
      STAB(I,J,K)=(1.+RI(I,J,K))
      KM(I,J,K)=PRIT(I,J)*KHOMO(I,J,K)/STAB(I,J,K)
      KH(I,J,K)=PRIT(I,J)*KHOMO(I,J,K)/STAB(I,J,K)**2
 280  CONTINUE
!
      RETURN
      END

