      SUBROUTINE ADVQ(QB,Q,DTI2,QF)
!     VERSION(09/27/90)
!
      INCLUDE 'comdeck'
      SAVE
!
!-----------------------------------------------------------------------
!     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS
!            FOR Q2 AND Q2L
!-----------------------------------------------------------------------
!
      DIMENSION QB(IM,JM,KB),Q(IM,JM,KB),QF(IM,JM,KB)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C)
!qa      DIMENSION AAMAX(IM,JM,KB),AAMAY(IM,JM,KB)
!qa      EQUIVALENCE (AAMAX,VH),(AAMAY,VHP)
!
      DO 10 K=1,KB
      DO 10 J=1,JM
      DO 10 I=1,IM
      XFLUX(I,J,K)=0.0
      YFLUX(I,J,K)=0.0
   10 QF(I,J,K)=0.0
!
!-------- HORIZONTAL ADVECTION -----------------------------------------
      DO 20 K=2,KBM1
      DO 20 J=2,JM
      DO 20 I=2,IM
      XFLUX(I,J,K)=0.25*(Q(I-1,J,K)+Q(I,J,K))
     .     *(XMFL3D(I,J,K-1)+XMFL3D(I,J,K))
 20   YFLUX(I,J,K)=0.25*(Q(I,J-1,K)+Q(I,J,K))
     .     *(YMFL3D(I,J,K-1)+YMFL3D(I,J,K))
!
!-------- HORIZONTAL DIFFUSION -----------------------------------------
!-------- ADD DIFFUSIVE FLUXES -----------------------------------------
      DO 30 J=2,JM
      DO 30 K=2,KBM1
      DO 30 I=2,IM
!NG05192010      AAMAX(I,J,K)=.5*(AAM(I,J,K)+AAM(I-1,J,K))
!NG05192010      AAMAY(I,J,K)=.5*(AAM(I,J,K)+AAM(I,J-1,K))
      AAMAX(I,J,K)=.5*(AAM(I,J,K)+AAM(I-1,J,K))*DUMWAD(I,J) ! Not in original NYHOPS W&D
      AAMAY(I,J,K)=.5*(AAM(I,J,K)+AAM(I,J-1,K))*DVMWAD(I,J) ! Not in original NYHOPS W&D
   30 CONTINUE
!
      DO 35 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 36 K=1,KBM1
   36         AAMAX(IC,JC,K)=0.0
            ELSE
              DO 37 K=1,KBM1
   37         AAMAX(ID,JD,K)=0.0
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 38 K=1,KBM1
   38         AAMAY(IC,JC,K)=0.0
            ELSE
              DO 39 K=1,KBM1
   39         AAMAY(ID,JD,K)=0.0
            ENDIF
      ENDIF
   35 CONTINUE
!
      DO 40 J=2,JM
      DO 40 K=2,KBM1
      DO 40 I=2,IM
      XFLUX(I,J,K)=XFLUX(I,J,K)
!NG04272010 W&D     .    -AAMAX(I,J,K)*(H(I,J)+H(I-1,J))*0.5*(H2(I,J)+H2(I-1,J))
     .    -AAMAX(I,J,K)*(H(I,J)+ETB(I,J)+H(I-1,J)+ETB(I-1,J))*0.5*
     .     (H2(I,J)+H2(I-1,J))*(QB(I,J,K)-QB(I-1,J,K))*DUM(I,J)
     .     /(H1(I,J)+H1(I-1,J))*DUMWAD(I,J) ! Not in original NYHOPS W&D
      YFLUX(I,J,K)=YFLUX(I,J,K)
!NG04272010 W&D     .    -AAMAY(I,J,K)*(H(I,J)+H(I,J-1))*0.5*(H1(I,J)+H1(I,J-1))
     .    -AAMAY(I,J,K)*(H(I,J)+ETB(I,J)+H(I,J-1)+ETB(I,J-1))*0.5*
     .     (H1(I,J)+H1(I,J-1))*(QB(I,J,K)-QB(I,J-1,K))*DVM(I,J)
     .    /(H2(I,J)+H2(I,J-1))*DVMWAD(I,J) ! Not in original NYHOPS W&D
   40 CONTINUE
!
!-------- VERTICAL ADVECTION -------------------------------------------
!-------- ADD FLUX TERMS & THEN STEP FORWARD IN TIME -------------------
      DO 50 J=2,JMM1
      DO 50 K=2,KBM1
      DO 50 I=2,IMM1
      QF(I,J,K)=(W(I,J,K-1)*Q(I,J,K-1)-W(I,J,K+1)*Q(I,J,K+1))
     .                     /(DZ(K)+DZ(K-1))*ART(I,J)
     .                      +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .                      +YFLUX(I,J+1,K)-YFLUX(I,J,K)
   50 QF(I,J,K)=((H(I,J)+ETB(I,J))*ART(I,J)*QB(I,J,K)-DTI2*QF(I,J,K))
     .             /((H(I,J)+ETF(I,J))*ART(I,J))
!
      RETURN
      END
