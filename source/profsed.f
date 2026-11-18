      SUBROUTINE PROFSED(F,DT2,WSET,WFBOT)

!     yang, UPDATE(Dec/29/1999), to limit the maximum deposition
!     yang, UPDATE(Dec/21/1999), to prevent from negative solution and to
!                                maintain mass balance.
!     yang, UPDATE(Apr/15/1999), for including Settling Rate and Deposition Rate
!                             in implicit scheme.
!     VERSION(09/27/90)

      INCLUDE 'comdeck'
      SAVE
!
!yang, start(Apr/15/1999)
      REAL*8 AD(IM,JM,KB),AL(IM,JM,KB),AU(IM,JM,KB)
      REAL*8 FF(IM,JM,KB),VHF(IM,JM,KB),VHPF(IM,JM,KB)
      DIMENSION WFBOT(IM,JM),WSET(IM,JM,KB)
!
      DIMENSION F(IM,JM,KB),DH(IM,JM)
      EQUIVALENCE (A,AD),(PROD,AL),(VH,VHF),(AD,FF)
      EQUIVALENCE (TPS,DH)
!yang end

      IF (TOR.EQ.'BAROTROPIC') RETURN
!
      UMOLPR=UMOL
!
!-----------------------------------------------------------------------
!
!        THE FOLLOWING SECTION SOLVES THE EQUATION
!         DT2*(KH*F')'-F=-FB
!
!-----------------------------------------------------------------------
!
      DO 90 J=2,JMM1
      DO 90 I=2,IMM1
      DH(I,J)=H(I,J)+ETF(I,J)
  90  CONTINUE
!
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
      IF(FSM(I,J).LE.0.) GO TO 100
      AL(I,J,1)=0.0
      AD(I,J,1)=1.0 + DT2*(KH(I,J,2)+UMOLPR*VARUF(I,J))/(DZZ(1)*DH(I,J)
     .     *DZ(1)*DH(I,J)) + DT2*WSET(I,J,2)/(DZ(1)*DH(I,J))
      AU(I,J,1)= - DT2*(KH(I,J,2)+UMOLPR*VARUF(I,J))/(DZZ(1)*DH(I,J)
     .            *DZ(1)*DH(I,J))

      DO 95 K=2,KBM1-1
      AL(I,J,K)= - DT2*(KH(I,J,K)+UMOLPR*VARUF(I,J))/(DZZ(K-1)*DH(I,J)
     .     *DZ(K)*DH(I,J)) - DT2*WSET(I,J,K)/(DZ(K)*DH(I,J))
      AD(I,J,K)=1.0 + DT2*(KH(I,J,K+1)+UMOLPR*VARUF(I,J))/
     .      (DZZ(K)*DH(I,J)
     .     *DZ(K)*DH(I,J)) + DT2*WSET(I,J,K+1)/(DZ(K)*DH(I,J))
     .     + DT2*(KH(I,J,K)+UMOLPR*VARUF(I,J))/
     .      (DZZ(K-1)*DH(I,J)*DZ(K)*DH(I,J)) 
      AU(I,J,K)= - DT2*(KH(I,J,K+1)+UMOLPR*VARUF(I,J))/(DZZ(K)*DH(I,J)*
     .             DZ(K)*DH(I,J))
  95  CONTINUE

      AL(I,J,KBM1)= - DT2*(KH(I,J,KBM1)+UMOLPR*VARUF(I,J))/
     .      (DZZ(KBM1-1)*DH(I,J)
     .     *DZ(KBM1)*DH(I,J)) - DT2*WSET(I,J,KBM1)/(DZ(K)*DH(I,J))
      AD(I,J,KBM1)=1.0 + DT2*(KH(I,J,KBM1)+UMOLPR*VARUF(I,J))/
     .      (DZZ(KBM1-1)*DH(I,J)
     .     *DZ(KBM1-1)*DH(I,J)) 
      AU(I,J,KBM1)=0.0

 100  CONTINUE
!
!  FORWARD SWEEP OF TRIDIAGONAL SCHEME

!-------- SURFACE BOUNDARY CONDITIONS - WFSURF -------------------------
!
!  ASSUME ZERO SURFACE FLUX FOR SEDIMENT
!
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      IF(FSM(I,J).LE.0.) GO TO 110
      VHF(I,J,1)=AD(I,J,1)
      VHPF(I,J,1)=F(I,J,1)/VHF(I,J,1)
  
      DO 101 K=2,KBM1
      VHF(I,J,K)=AD(I,J,K)-AL(I,J,K)*AU(I,J,K-1)/VHF(I,J,K-1)
      IF(K.LT.KBM1)THEN
        VHPF(I,J,K)=(F(I,J,K)-AL(I,J,K)*VHPF(I,J,K-1))/VHF(I,J,K)
      ELSE
! yang start (04/15/99), WFBOT (g/cm^2) was computed in sedflx.f or suslod.f 
! by 1 dt only.
!       VHPF(I,J,K)=(F(I,J,K)+ (DT2*2.0*WFBOT(I,J)/(DZ(K)*DH(I,J)*100.))
!    .             -AL(I,J,K)*VHPF(I,J,K-1))/VHF(I,J,K)
! yang start (dec/29/99), WFBOT (g/cm^2) was computed in sedflx.f or suslod.f 
!  CALC. MAX. DEPOSITION BASED UPON MASS AT LAYER kb-1 (UNITS: g/cm**2)
!
        IF(WFBOT(I,J).LT.0.0)THEN
          WFBOTMAX=100.*F(I,J,KBM1)*DH(I,J)*DZ(KBM1)
          WFBOT(I,J)=AMAX1(WFBOT(I,J),-WFBOTMAX)
        ENDIF
! yang end
!
        VHPF(I,J,K)=(F(I,J,K)+ (2.0*WFBOT(I,J)/(DZ(K)*DH(I,J)*100.))
     .             -AL(I,J,K)*VHPF(I,J,K-1))/VHF(I,J,K)
! yang end

! yang start (dec/21/1999): set minimum of deposition flux to prevent from 
!                           negative conc.
        IF (VHPF(I,J,K).LE. 0.0 )THEN
          IF (WFBOT(I,J).LE. 0.0 ) WFBOT(I,J)=0.0
          VHPF(I,J,K)=0.0  !! loss all mass
        ENDIF
! yang end
      ENDIF
 101  CONTINUE
 110  CONTINUE
!
      DO 120 K=1,KB
      DO 120 J=1,JM
      DO 120 I=1,IM
 120  FF(I,J,K)=F(I,J,K)
!
!  BACKWARD SWEEP FOR SOLUTION

      DO 130 J=1,JM
      DO 130 I=1,IM
        IF(FSM(I,J).LE.0.) GO TO 130
        FF(I,J,KBM1)=VHPF(I,J,KBM1)

        DO 105 K=2,KBM1
          KI=KB-K
          FF(I,J,KI)=VHPF(I,J,KI)-AU(I,J,KI)*FF(I,J,KI+1)/VHF(I,J,KI)
 105  CONTINUE
 130  CONTINUE
!
      DO 140 J=1,JM
      DO 140 K=1,KB
      DO 140 I=1,IM
!NG09232011      IF(FSM(I,J).LE.0.) GO TO 140
      IF(DH(I,J).LE.WETMIN) GOTO 140
      F(I,J,K)=FF(I,J,K)
 140  CONTINUE
!
      RETURN
      END
