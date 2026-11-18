      SUBROUTINE VERTVL(DTI2)
!     VERSION(09/27/90)
!
      INCLUDE 'comdeck'
      SAVE
!
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),ZMFL3D(IM,JM,KB)
      EQUIVALENCE (XFLUX,A),(YFLUX,C),(ZMFL3D,VHP)
      REAL CORRNG(IM,JM)
!
!-------- CALCULATE NEW VERTICAL VELOCITY ------------------------------
      DO 100 J=1,JM
      DO 100 I=1,IM
      CORRNG(I,J)=0.0 ! NG
 100  W(I,J,1)=(QEVAP2D(I,J)-QPREC2D(I,J))*RAMP*FSM(I,J)
!
      DO 101 K=1,KB
      DO 101 J=1,JM
      DO 101 I=1,IM
 101  ZMFL3D(I,J,K)=0.0
!-----------------------------------------------------------------------
!         IMPOSE MASS FLUX BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
      DO 120 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      DO 120 K=1,KBM1
!
        ZMFL3D(ID,JD,K)=ZMFL3D(ID,JD,K)+
     +          QDIFF(N)*RAMP*VDDIST(N,K)/100.0/DZ(K)
 120  CONTINUE

!NG---------------------------------------------------------------------
!NG04232008 IMPOSE ELEVATION POTENTIAL INTERIOR CONDITIONS
!-----------------------------------------------------------------------
!NG04232008 NOTE THAT W IS POSITIVE DOWNWARDS!
      DO N=1,NUMEPC
         I=IEPC(N)
         J=JEPC(N)
         !rrlocal=0.
         DO K=1,KBM1
!NG_ALREADY_RAMPED          ZMFL3D(I,J,K)=ZMFL3D(I,J,K) + QEPCNUDG(N)*RAMP
          ZMFL3D(I,J,K)=ZMFL3D(I,J,K) + QEPCNUDG(N) ! Q>0(incomming) -> ZQ>0 -> W=UP(negative)
          !rrlocal=rrlocal+dz(k)*
     +    !  ((-qepcnudg(n)/art(i,j))+(etf(i,j)-etb(i,j))/dti2)
          !write (*,*)k+1,-qepcnudg(n)/art(i,j),(etf(i,j)-etb(i,j))/dti2,
     +    !  rrlocal,qepcnudg(n)
         ENDDO
      ENDDO
!-----------------------------------------------------------------------

!NG 05/03/2010
!NG CALCULATE ETF CORRECTION TERM SO THAT W(KB)=0.0 (SOLVE LAST EQUATION WITH KINEMATIC BC FOR ETF)
!NG
!NG   TEMPORARILLY CHANGE FSM (Avoids setting up new mask to skip OBC and DISCHARGERS)
      DO 245 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      FSM(IC,JC)=0.0
  245 CONTINUE
      DO 246 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      FSM(IE,JE)=0.0
  246 CONTINUE
!NG      DO 247 N=1,NUMECR ! This should not be needed; it is set to 0 in setdom_corners
!NG      IE=ICHARD(N)
!NG      JE=JCHARD(N)
!NG      FSM(IE,JE)=0.0
!NG  247 CONTINUE
!NG CALCULATE AND APPLY CORRECTION TERM AT INTERNAL POINTS
      DO 140 J=2,JMM1
      DO 140 I=2,IMM1
      IF(FSM(I,J).EQ.0.) GO TO 140
      DO 141 K=1,KBM1
  141  CORRNG(I,J)=CORRNG(I,J)
     .    +DZ(K)*(-ZMFL3D(I,J,K)+XMFL3D(I+1,J,K)-XMFL3D(I,J,K)
     .            +YMFL3D(I,J+1,K)-YMFL3D(I,J,K))/ART(I,J)
      CORRNG(I,J)=CORRNG(I,J)+W(I,J,1) ! NG12202010 FURTHER CORRECTION FOR PREC/EVAP
      ETF(I,J)=ETB(I,J) - DTI2*CORRNG(I,J)
!NG      ETF(I,J)=AMAX1(ETF(I,J),-H(I,J)+1.E-3) ! DISALLOW NEGATIVE DT=ETF+H
      ETF(I,J)=AMAX1(ETF(I,J),-H(I,J)+AMAX1(WETMIN-1.E-3,1.E-3)) ! DISALLOW NEGATIVE DT=ETF+H, SAME AS IN SETDOM
  140 CONTINUE
!NG   RESTORE FSM
      DO 248 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      FSM(IC,JC)=1.0
  248 CONTINUE
      DO 249 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      FSM(IE,JE)=1.0
  249 CONTINUE
!NG      DO 250 N=1,NUMECR ! This should not be 1; it is set to 0 in setdom_corners
!NG      IE=ICHARD(N)
!NG      JE=JCHARD(N)
!NG      FSM(IE,JE)=1.0
!NG  250 CONTINUE
!
!NG      WRITE(83,*) QEVAP2D(23,39),QPREC2D(23,39),
!NG     .            QEVAP2D(23,40),QPREC2D(23,40)
!NG      WRITE(83,*) CORRNG(23,39),CORRNG(23,40)
!NG---------------------------------------------------------------------
!-----------------------------------------------------------------------
      DO 110 J=2,JMM1
      DO 110 K=1,KBM1
      DO 110 I=2,IMM1
 110  W(I,J,K+1)=W(I,J,K)
     .    +DZ(K)*((-ZMFL3D(I,J,K)+XMFL3D(I+1,J,K)-XMFL3D(I,J,K)
     .            +YMFL3D(I,J+1,K)-YMFL3D(I,J,K))/ART(I,J)
     .                        +(ETF(I,J)-ETB(I,J))/DTI2 )
!
!NG      WRITE(83,*) 
!NG     +TIME,W(10,10,KB),T(10,10,KB),S(10,10,KB),DT(10,10)
!NG      WRITE(83,*) 
!NG     +TIME,W(47,126,KB),T(47,126,KB),S(47,126,KB),DT(47,126)
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  W(I,J,KB)=0.0  
!
      RETURN
      END
