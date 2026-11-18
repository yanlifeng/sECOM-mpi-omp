      SUBROUTINE SMAG
!     VERSION(05/14/2012)
!
      INCLUDE 'comdeck'
      SAVE
!
!
!----------------------------------------------------------------------|
!         0.01 < horcon < 0.20                                         |
!         HOR. VISC. = CONST.*H1*H2*SQRT((DU/H1)**2+(DV/H2)**2         |
!                  +0.5*(DU/H2+DV/H1)**2)                              |
!----------------------------------------------------------------------|
!
!
      IF(TOR.NE.'BAROTROPIC') THEN
!
      DO 110 J=2,JMM1
      DO 110 K=1,KBM1
      DO 110 I=2,IMM1
      IF(FSM(I,J).EQ.0.) GOTO 110
      AAM(I,J,K)=HORCON*VARHF(I,J)*ART(I,J)*
     +  SQRT( ((U(I+1,J,K)-U(I,J,K))/H1(I,J))**2+
     +        ((V(I,J+1,K)-V(I,J,K))/H2(I,J))**2+
     +  0.5/(DVM(I,J+1)*DVMWAD(I,J+1)+DVM(I,J)*DVMWAD(I,J)+
     +       DUM(I+1,J)*DUMWAD(I+1,J)+DUM(I,J)*DUMWAD(I,J)+1.E-30)*
     +  (0.5*ABS(U(I,J+1,K)+U(I+1,J+1,K)-U(I,J,K)-U(I+1,J,K))
     +   /(0.5*H2(I,J)+0.5*H2(I,J+1))*DVM(I,J+1)*DVMWAD(I,J+1)
     +  +0.5*ABS(U(I,J,K)+U(I+1,J,K)-U(I,J-1,K)-U(I+1,J-1,K))
     +   /(0.5*H2(I,J)+0.5*H2(I,J-1))*DVM(I,J)*DVMWAD(I,J)
     +  +0.5*ABS(V(I+1,J,K)+V(I+1,J+1,K)-V(I,J,K)-V(I,J+1,K))
     +   /(0.5*H1(I,J)+0.5*H1(I+1,J))*DUM(I+1,J)*DUMWAD(I+1,J)
     +  +0.5*ABS(V(I,J,K)+V(I,J+1,K)-V(I-1,J,K)-V(I-1,J+1,K))
     +   /(0.5*H1(I,J)+0.5*H1(I-1,J))*DUM(I,J)*DUMWAD(I,J)
     +  )**2)
      AAM(I,J,K)=AMIN1(AAM(I,J,K),HORCON*VARHF(I,J)*ART(I,J)*
     .                 SQRT((25./H1(I,J))**2+(25./H2(I,J))**2)) ! CAP IT TO +/- 5m/s
 110  CONTINUE
!
      DO 120 N=1,NUMQBC
      DO 120 K=1,KBM1
      IC=IQC(N)
      JC=JQC(N)
      AAM(IC,JC,K)=0.0
  120 CONTINUE

      DO 121 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 121 K=1,KBM1
      AAM(IE,JE,K)=AAM(IC,JC,K)
  121 CONTINUE
! 
      DO 96 J=1,JM
      DO 96 I=1,IM
 96   AAM2D(I,J)=0.0
!
      DO 100 J=1,JM
      DO 100 K=1,KBM1
      DO 100 I=1,IM
      AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
 100  CONTINUE
!
      ELSE
!
      DO 130 J=2,JMM1
      DO 130 I=2,IMM1
      IF(FSM(I,J).EQ.0.) GOTO 130
      AAM2D(I,J)=HORCON*VARHF(I,J)*ART(I,J)*
     +  SQRT( ((UA(I+1,J)-UA(I,J))/H1(I,J))**2+
     +        ((VA(I,J+1)-VA(I,J))/H2(I,J))**2+
     +  0.5/(DVM(I,J+1)*DVMWAD(I,J+1)+DVM(I,J)*DVMWAD(I,J)
     +      +DUM(I+1,J)*DUMWAD(I+1,J)+DUM(I,J)*DUMWAD(I,J)+1.E-30)*
     +  (0.5*ABS(UA(I,J+1)+UA(I+1,J+1)-UA(I,J)-UA(I+1,J))
     +   /(0.5*H2(I,J)+0.5*H2(I,J+1))*DVM(I,J+1)*DVMWAD(I,J+1)
     +  +0.5*ABS(UA(I,J)+UA(I+1,J)-UA(I,J-1)-UA(I+1,J-1))
     +   /(0.5*H2(I,J)+0.5*H2(I,J-1))*DVM(I,J)*DVMWAD(I,J)
     +  +0.5*ABS(VA(I+1,J)+VA(I+1,J+1)-VA(I,J)-VA(I,J+1))
     +   /(0.5*H1(I,J)+0.5*H1(I+1,J))*DUM(I+1,J)*DUMWAD(I+1,J)
     +  +0.5*ABS(VA(I,J)+VA(I,J+1)-VA(I-1,J)-VA(I-1,J+1))
     +   /(0.5*H1(I,J)+0.5*H1(I-1,J))*DUM(I,J)*DUMWAD(I,J)
     +  )**2)
      AAM2D(I,J)=AMIN1(AAM2D(I,J),HORCON*VARHF(I,J)*ART(I,J)*
     .                 SQRT((25./H1(I,J))**2+(25./H2(I,J))**2)) ! CAP IT TO +/- 5m/s
 130  CONTINUE
!
      DO 112 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      AAM2D(IC,JC)=0.0
 112  CONTINUE
!
      DO 113 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      AAM2D(IE,JE)=AAM2D(IC,JC)
 113  CONTINUE

      ENDIF

!NG      WRITE(*,*) TIME
!NG      WRITE(*,*) (AAM2D(IETA(N),JETA(N)),N=1,NUMEBC)
!
!
      RETURN
      END
