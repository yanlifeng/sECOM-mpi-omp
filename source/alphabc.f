      SUBROUTINE ALPHABC(DTI2)
!C      SUBROUTINE ALPHABC
!
!     VERSION(03/05/99)
!       UPWIND FLUX AT OPEN SALT/TEMP BOUNDARY
!
      INCLUDE 'comdeck'
!
!     MODIFIED BY QUAMRUL QA ON 5/3/96 *********************************
!     ADVECTIVE/RELAXATION BOUNDARY FOR T/S 
      DO 238 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
!
      DO 239 K=1,KBM1
!
      TBDY=TBDRY(N,K)
      SBDY=SBDRY(N,K)
!
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
!
!-------- EAST SIDE ----------------------------------------------------
      VEL=U(IE,JE,K)
      IF (VEL.LE.0.0) THEN
       CPH=0.0
       ALPHAI=ALPHA
      ELSE
       ALPHAI=0.0 
       CPH=VEL*DTI*2.0/(H1(IE,JE)+H1(IE-1,JE))
      END IF
!
      UF(IE,JE,K)=TB(IE,JE,K)+CPH*(2.0*T(IE-1,JE,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*ALPHAI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0+CPH)
      VF(IE,JE,K)=SB(IE,JE,K)+CPH*(2.0*S(IE-1,JE,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*ALPHAI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0+CPH)
      GO TO 239
!
      ELSE IF (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
!
!-------- WEST SIDE ----------------------------------------------------
      VEL=U(IE+1,JE,K)
      IF (VEL.GE.0.0) THEN
      CPH=0.0
      ALPHAI=ALPHA
      ELSE
      CPH=VEL*DTI*2.0/(H1(IE,JE)+H1(IE+1,JE))
      ALPHAI=0.0 
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)-CPH*(2.0*T(IE+1,JE,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*ALPHAI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0-CPH)
      VF(IE,JE,K)=SB(IE,JE,K)-CPH*(2.0*S(IE+1,JE,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*ALPHAI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0-CPH)
      GO TO 239
      ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
!
!-------- NORTH SIDE ---------------------------------------------------
      VEL=V(IE,JE,K)
      IF (VEL.LE.0.0) THEN
      CPH=0.0
      ALPHAI=ALPHA
      ELSE
      CPH=VEL*DTI*2.0/(H2(IE,JE)+H2(IE,JE-1))
      ALPHAI=0.0 
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)+CPH*(2.0*T(IE,JE-1,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*ALPHAI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0+CPH)
      VF(IE,JE,K)=SB(IE,JE,K)+CPH*(2.0*S(IE,JE-1,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*ALPHAI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0+CPH)
      GO TO 239
      ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
!
!-------- SOUTH SIDE ---------------------------------------------------
      VEL=V(IE,JE+1,K)
      IF (VEL.GE.0.0) THEN
      CPH=0.0
      ALPHAI=ALPHA
      ELSE
      CPH=VEL*DTI*2.0/(H2(IE,JE)+H2(IE,JE+1))
      ALPHAI=0.0 
      END IF
      UF(IE,JE,K)=TB(IE,JE,K)-CPH*(2.0*T(IE,JE+1,K)-TB(IE,JE,K)) 
     1 +(TBDY-TB(IE,JE,K))*ALPHAI*DTI2 
      UF(IE,JE,K)=UF(IE,JE,K)/(1.0-CPH)
      VF(IE,JE,K)=SB(IE,JE,K)-CPH*(2.0*S(IE,JE+1,K)-SB(IE,JE,K)) 
     1 +(SBDY-SB(IE,JE,K))*ALPHAI*DTI2 
      VF(IE,JE,K)=VF(IE,JE,K)/(1.0-CPH)
!
!-------- DONE ---------------------------------------------------------
      END IF
  239 CONTINUE
  238 CONTINUE
      RETURN
      END
