      SUBROUTINE PROFU(DT2)
!     VERSION(09/27/90)
!
      INCLUDE 'comdeck'
      SAVE
!
      DIMENSION DH(IM,JM)
      REAL TPSI(IM,JM)
!
!-----------------------------------------------------------------------
!
!        THE FOLLOWING SECTION SOLVES THE EQUATION
!         DT2*(KM*U')'-U=-UB
!
!-----------------------------------------------------------------------
!
      DO 85 J=2,JMM1
      DO 85 I=3,IMM1
  85  DH(I,J)=.5*(H(I,J)+ETF(I,J)+H(I-1,J)+ETF(I-1,J))
!
      DO 90 K=1,KB
      DO 90 J=2,JMM1
      DO 90 I=3,IMM1
  90  C(I,J,K)=(KM(I,J,K)+KM(I-1,J,K))*.5
!
      DO 100 J=2,JMM1
      DO 100 K=2,KBM1
      DO 100 I=3,IMM1
      A(I,J,K-1)=-DT2*(C(I,J,K)+UMOL*VARUF(I,J))/
     .   (DZ(K-1)*DZZ(K-1)*DH(I,J)*DH(I,J))
      C(I,J,K)=-DT2*(C(I,J,K)+UMOL*VARUF(I,J))/
     .   (DZ(K)*DZZ(K-1)*DH(I,J)*DH(I,J))
 100  CONTINUE
!
      DO 110 J=2,JMM1
      DO 110 I=3,IMM1
      VH(I,J,1)=A(I,J,1)/(A(I,J,1)-1.)   
 110  VHP(I,J,1)=(-DT2*(.5*((100.-AREAICE(I,J))/100.*TXSURF(I,J)+      ! NG FOR ICE
     .                      (100.-AREAICE(I-1,J))/100.*TXSURF(I-1,J))  ! NG FOR ICE
     .                     +TXICE(I,J))/                               ! NG FOR ICE
     1   (-DZ(1)*DH(I,J))-UF(I,J,1))/(A(I,J,1)-1.)
!
      DO 101 K=2,KBM2
      DO 101 J=2,JMM1
      DO 101 I=3,IMM1
      VHP(I,J,K)=1./(A(I,J,K)+C(I,J,K)*(1.-VH(I,J,K-1))-1.)
      VH(I,J,K)=A(I,J,K)*VHP(I,J,K)
      VHP(I,J,K)=(C(I,J,K)*VHP(I,J,K-1)-UF(I,J,K))*VHP(I,J,K)
 101  CONTINUE
!
      DO 102 J=2,JMM1
      DO 102 I=3,IMM1
      TPS(I,J)=0.5*(CBC(I,J)+CBC(I-1,J))
     .     *SQRT(UB(I,J,KBM1)**2+(.25*(VB(I,J,KBM1)
     .     +VB(I,J+1,KBM1)+VB(I-1,J,KBM1)+VB(I-1,J+1,KBM1)))**2)
      TPSI(I,J)=0.5*(AREAICE(I,J)/100.*VARWIF(I,J)      ! NG FOR ICE
     .              +AREAICE(I-1,J)/100.*VARWIF(I-1,J)) ! NG FOR ICE
     .     *SQRT(UB(I,J,1)**2+(.25*(VB(I,J,1)
     .     +VB(I,J+1,1)+VB(I-1,J,1)+VB(I-1,J+1,1)))**2)
      UF(I,J,KBM1)=(C(I,J,KBM1)*VHP(I,J,KBM2)-UF(I,J,KBM1))/(TPS(I,J)
     . *DT2/(-DZ(KBM1)*DH(I,J))-1.-(VH(I,J,KBM2)-1.)*C(I,J,KBM1))
 102  UF(I,J,KBM1)=UF(I,J,KBM1)*DUM(I,J)*DUMWAD(I,J) ! Not in original NYHOPS
!
      DO 103 K=2,KBM1
      KI=KB-K
      DO 103 J=2,JMM1
      DO 103 I=3,IMM1
      UF(I,J,KI)=(VH(I,J,KI)*UF(I,J,KI+1)+VHP(I,J,KI))
     .          * DUM(I,J)*DUMWAD(I,J) ! Not in original NYHOPS
 103  CONTINUE     
!
      DO 120 J=2,JMM1
      DO 120 I=3,IMM1
      TXICE(I,J)=TPSI(I,J)*UF(I,J,1)*DUM(I,J)*DUMWAD(I,J) ! NG FOR ICE, NOTE SIGN AS TSURF, NOT TBOT
 120  TXBOT(I,J)=-TPS(I,J)*UF(I,J,KBM1)*DUM(I,J)*DUMWAD(I,J) ! Not in original NYHOPS
!
      RETURN
      END
