      SUBROUTINE BCOND(IDX,DTI2,KS)
!
!
!       UPWIND FLUX AT OPEN SALT/TEMP BOUNDARY
!
      INCLUDE 'comdeck'
      PARAMETER (INTMAX = 100)  ! JI
      REAL SSTINT(INTMAX), WT(INTMAX), WQ(INTMAX), RLWC(INTMAX)  ! JI
*Khan 072898
      REAL TH2D(IM,JM)  ! 2D Heat Flux N.KIM 7.7.98
*Khan 072898
!
!NG 12162010
      REAL PRECTEMP(IM,JM)
      REAL AALAT,AALON
!NG 12162010
      
      CHARACTER STR*80 ! A.Jordi
      LOGICAL CLEAR_RAIN
      
      SAVE
!
!  MODIFIED BY QUAMRUL QA ON 5/3/96 **************************************
      Data SBC /2.0411E-7/     ! 4*Boltzman Constant 4*5.67E-8 Army
      Data ZERO /0.0/  
!     Following two lines are for ROSATI and Miyakoda (RANDMFLX)
!NG   Data DHR /1.0/           ! Period of Solar radiation Average in Hr
!     Data REFL /0.1/          ! REFL - FRACTION SWOBS THAT REFLECTS BACK
!
!
!===============>    S2     M2     N2     K1     P1     O1    
      DHR=AMAX1(DTI*ISKILL/3600.,0.1) ! NG Period of Solar radiation averaging in hrs based on tsepic averaging (or at least 6 minutes)
      DHR=AMIN1(DHR,1.)               ! NG And at most 1 hour
      TPI=6.283185
      PERIOD(1)=43200.
      PERIOD(2)=44712.
      PERIOD(3)=45570.
      PERIOD(4)=86164.
      PERIOD(5)=86637.
      PERIOD(6)=92950.
!
!-----------------------------------------------------------------------
!        SPECIFICATION OF OPEN BOUNDARY CONDITIONS
!   NOTE THAT AT J=2 U CALCULATION EXCLUDES ADVECTION AND DIFFUSION
!-----------------------------------------------------------------------
!        I-1           I         I+1
!
!       U(JM)=      EL(JM)      U(JM)
!
!                    V(JM)
!
!     *U(JMM1)*   *EL(JMM1)*   *U(JMM1)*
!
!                  *V(JMM1)*
!
!     *U(JMM2)*   *EL(JMM2)*   *U(JMM2)*
!                       "                                         V(IM)
!                       "
!                       "                *U(IMM1)* EL(IMM1) U(M)  EL(IM)
!                       "                              = = = = = = =
!                       "                             V(IMM1)     V(IM)
!                       "
!       *U(3)*       *EL(3)*   *U(3)*                      BC ON EL & T
!
!                    *V(3)*
!                 "
!       *U(2)*    "  *EL(2)*   *U(2)*
!                 "         "         "
!                  =  V(2)  "         "  BC ON V
!                           "         "
!          U(1)       EL(1)=    U(1) =   BC ON T
!-----------------------------------------------------------------------
!         IDX IDENTIFIES WHICH VARIABLES ARE CONSIDERED
!              1=SURFACE ELEVATION  
!              2=EXTERNAL MODE U,V 
!              3=INTERNAL MODE U,V 
!              4=TEMP,SAL
!              5=W VELOCITY
!              6=KM,KH,Q2,Q2L,L 
!              7=SURFACE FORCING AND TEMPORAL CYCLING
!              8= CONS. TRACER
!              9= SEDIMENT TRANSPORT
!             10= CHEMICAL TRANSPORT
!         
!-----------------------------------------------------------------------
!
!
!     MODIFIED BY QUAMRUL QA ON 3/10/99 ********************************
        IF(INTX.EQ.ISTART) THEN
           Open (77,file='heat.dat',FORM='FORMATTED')
        END IF
!
!     12/13: boundary entries only; an existing external-mode team
!     applies the identical dense mask immediately afterwards.
      GOTO (10,20,30,40,50,60,70,80,81,181,90,10,20), IDX
!
!     ******************************************************************
!
 10   CONTINUE
!
!-------- EXTERNAL ELEVATION BOUNDARY CCONDITIONS ----------------------
      DO 100 N=1,NUMEBC  
      II=IETA(N)
      JJ=JETA(N)
!     IF(OPTEBC(1:5).NE.'DATA ') THEN
!     FORCE=0.0
!     DO 110 I=1,6  
!110  FORCE=AMP(N,I)*COS(TPI/PERIOD(I)*(DTI*FLOAT(INTX)-PHASE(N,I)))
!    .     +FORCE
!     ELF(II,JJ)=(EMEAN(N)+FORCE)*RAMP
!     ELSE
      ELF(II,JJ)=EBDRY(N)*RAMP
!     ENDIF
 100  CONTINUE
!


      DO 120 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      ELF(IC,JC)=ELF(ID,JD)
120   CONTINUE
!
      IF(IDX.EQ.12) RETURN
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  ELF(I,J)=ELF(I,J)*FSM(I,J)
      RETURN
!
!-------- EXTERNAL VELOCITY BOUNDARY CONDITIONS-------------------------
  20  CONTINUE
      DO 200 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)  
      FRESH=QDIS(N)/(H(IC,JC)+ELF(IC,JC))*RAMP
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            UAF(IC,JC)=-FRESH/H2(IC,JC)
         ELSE
            UAF(ID,JD)=+FRESH/H2(ID,JD)
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            VAF(IC,JC)=-FRESH/H1(IC,JC)
         ELSE
            VAF(ID,JD)=+FRESH/H1(ID,JD)
         ENDIF
      ENDIF
 200  CONTINUE
!
      IF (OPTEBC(15:20).EQ.'LINEAR') THEN ! IF LINEAR OPTION ENABLED AT THE END OF OPTEBC
       DO 209 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        IC=ICON(N)
        JC=JCON(N)
        IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         VAF(IE,JE)=0.0 ! IF NEEDED LINEARIZE/PERPENDICULAR FLOW
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         VAF(IE,JE)=0.0 ! IF NEEDED LINEARIZE/PERPENDICULAR FLOW
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
         UAF(IE,JE)=0.0 ! IF NEEDED LINEARIZE/PERPENDICULAR FLOW
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
         UAF(IE,JE)=0.0 ! IF NEEDED LINEARIZE/PERPENDICULAR FLOW
        END IF
 209   CONTINUE
      ELSE
       DO 210 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        IC=ICON(N)
        JC=JCON(N)
        IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         VAF(IE,JE)=VAF(IE-1,JE)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         VAF(IE,JE)=VAF(IE+1,JE)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
         UAF(IE,JE)=UAF(IE,JE-1)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
         UAF(IE,JE)=UAF(IE,JE+1)
        END IF
 210   CONTINUE
      ENDIF
      IF(IDX.EQ.13) RETURN
      DO 135 J=1,JM
      DO 135 I=1,IM
      TXBOT(I,J)=TXBOT(I,J)*DUM(I,J) ! CNG03022009 REGARDLESS, Not in original NYHOPS
      TYBOT(I,J)=TYBOT(I,J)*DVM(I,J) ! CNG03022009 REGARDLESS, Not in original NYHOPS
      UAF(I,J)=UAF(I,J)*DUM(I,J)
 135  VAF(I,J)=VAF(I,J)*DVM(I,J)
      RETURN
!
!-------- INTERNAL VELOCITY BOUNDARY CONDITIONS ------------------------
  30  CONTINUE
      DO 320 N=1,NUMQBC  
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 320 K=1,KBM1
      FRESH=QDIS(N)/(H(IC,JC)+ELF(IC,JC))*RAMP*VQDIST(N,K)/100.
      IF(JD.EQ.JC) THEN
         IF(ID.LT.IC) THEN
            UF(IC,JC,K)=-FRESH/(H2(IC,JC)*DZ(K))
         ELSE
            UF(ID,JD,K)=+FRESH/(H2(ID,JD)*DZ(K))
         ENDIF
      ELSE
         IF(JD.LT.JC) THEN
            VF(IC,JC,K)=-FRESH/(H1(IC,JC)*DZ(K))
         ELSE
            VF(ID,JD,K)=+FRESH/(H1(ID,JD)*DZ(K))
         ENDIF
      ENDIF
  320 CONTINUE
!
      DO 321 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 321 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        VF(IE,JE,K)=VF(IE-1,JE,K)
        TYBOT(IE,JE)=TYBOT(IE-1,JE)
        TYICE(IE,JE)=TYICE(IE-1,JE) ! NG FOR ICE
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        VF(IE,JE,K)=VF(IE+1,JE,K)
        TYBOT(IE,JE)=TYBOT(IE+1,JE)
        TYICE(IE,JE)=TYICE(IE+1,JE) ! NG FOR ICE
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE-1,K)
        TXBOT(IE,JE)=TXBOT(IE,JE-1)
        TXICE(IE,JE)=TXICE(IE,JE-1) ! NG FOR ICE
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE+1,K)
        TXBOT(IE,JE)=TXBOT(IE,JE+1)
        TXICE(IE,JE)=TXICE(IE,JE+1) ! NG FOR ICE
      END IF
  321 CONTINUE
! 
      DO 360 J=1,JM
      DO 360 K=1,KBM1
      DO 360 I=1,IM
      UF(I,J,K)=UF(I,J,K)*DUM(I,J)
      VF(I,J,K)=VF(I,J,K)*DVM(I,J)
      TXBOT(I,J)=TXBOT(I,J)*DUM(I,J) ! CNG03022009 REGARDLESS, Not in original NYHOPS
      TYBOT(I,J)=TYBOT(I,J)*DVM(I,J) ! CNG03022009 REGARDLESS, Not in original NYHOPS
      TXICE(I,J)=TXICE(I,J)*DUM(I,J) ! NG FOR ICE
      TYICE(I,J)=TYICE(I,J)*DVM(I,J) ! NG FOR ICE
 360  CONTINUE
      RETURN
!
!-------- TEMPERATURE & SALINITY BOUNDARY CONDITIONS -------------------
  40  CONTINUE

!NG10062011 NEXT LOOP FOR WAD
      DO 344 J=2,JM-1 
      DO 344 K=1,KB 
      DO 344 I=2,IM-1 
       IF(FSM(I,J).EQ.0) GOTO 344
       IF(H(I,J)+ETF(I,J).LE.WETMIN) THEN 
        UF(I,J,K)=TB(I,J,K) ! NG10062011 DO NOT UPDATE DRY SEGMENTS, OR REINITIALIZE FROM A WET NEIGHBOOR BELOW TO PREPARE FOR INUNDATION
        VF(I,J,K)=SB(I,J,K) ! NG10062011 DO NOT UPDATE DRY SEGMENTS, OR REINITIALIZE FROM A WET NEIGHBOOR BELOW TO PREPARE FOR INUNDATION
        IF(H(I+1,J)+ETF(I+1,J).GT.WETMIN) THEN
         UF(I,J,K)=UF(I+1,J,K)
         T(I,J,K)=T(I+1,J,K)
         TB(I,J,K)=TB(I+1,J,K)
         VF(I,J,K)=VF(I+1,J,K)
         S(I,J,K)=S(I+1,J,K)
         SB(I,J,K)=SB(I+1,J,K)
        ENDIF
        IF(H(I-1,J)+ETF(I-1,J).GT.WETMIN) THEN
         UF(I,J,K)=UF(I-1,J,K)
         T(I,J,K)=T(I-1,J,K)
         TB(I,J,K)=TB(I-1,J,K)
         VF(I,J,K)=VF(I-1,J,K)
         S(I,J,K)=S(I-1,J,K)
         SB(I,J,K)=SB(I-1,J,K)
        ENDIF
        IF(H(I,J+1)+ETF(I,J+1).GT.WETMIN) THEN
         UF(I,J,K)=UF(I,J+1,K)
         T(I,J,K)=T(I,J+1,K)
         TB(I,J,K)=TB(I,J+1,K)
         VF(I,J,K)=VF(I,J+1,K)
         S(I,J,K)=S(I,J+1,K)
         SB(I,J,K)=SB(I,J+1,K)
        ENDIF
        IF(H(I,J-1)+ETF(I,J-1).GT.WETMIN) THEN
         UF(I,J,K)=UF(I,J-1,K)
         T(I,J,K)=T(I,J-1,K)
         TB(I,J,K)=TB(I,J-1,K)
         VF(I,J,K)=VF(I,J-1,K)
         S(I,J,K)=S(I,J-1,K)
         SB(I,J,K)=SB(I,J-1,K)
        ENDIF
       ENDIF
 344  CONTINUE

      DO 235 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 235 K=1,KBM1
      IF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).GT.0.0) THEN
      UF(IC,JC,K)=TDIS(N)
      VF(IC,JC,K)=SDIS(N)
      ELSE
      UF(IC,JC,K)=UF(ID,JD,K)
      VF(IC,JC,K)=VF(ID,JD,K)
      END IF
 235  CONTINUE

!  FOR DIFFUSER   (NKIM  5.8.97)
      Do 115 N = 1, NUMDBC
        ID = IDD(N)
        JD = JDD(N)
        Do 112 K = 1, KBM1
          If (VDDIST(N,K).NE.0.0.AND.QDIFF(N).GT.0.0) Then
            TDIF(N,K) = TDIFF(N)
            SDIF(N,K) = SDIFF(N)
            IF(H(ID,JD)+ETF(ID,JD).LE.WETMIN) THEN ! IF DIFFUSER ON DRY CELL (NG04132011)
              UF(ID,JD,K)=TDIF(N,K)       ! SET TO DIFFUSER'S SCALAR
              VF(ID,JD,K)=SDIF(N,K)       ! SET TO DIFFUSER'S SCALAR
              TB(ID,JD,K)=TDIF(N,K)       ! SET TO DIFFUSER'S SCALAR
              T(ID,JD,K)=TDIF(N,K)        ! SET TO DIFFUSER'S SCALAR
              SB(ID,JD,K)=SDIF(N,K)       ! SET TO DIFFUSER'S SCALAR
              S(ID,JD,K)=SDIF(N,K)        ! SET TO DIFFUSER'S SCALAR
            ENDIF
          Else
            TDIF(N,K) = UF(ID,JD,K)
            SDIF(N,K) = VF(ID,JD,K)
          End If
  112   Continue
  115 Continue

!  FOR ELEVATION POTENTIAL   (NG04232008)
!NG      Do 117 N = 1, NUMEPC
!NG        I = IEPC(N)
!NG        J = JEPC(N)
!NG        Do 114 K = 1, KBM1
!NG          If (QEPCNUDG(N).GT.0.0) Then
!NG2D           TQEPC(N,K) = TQEPCIN(N,1)
!NG2D           SQEPC(N,K) = SQEPCIN(N,1)
!NGUSEF           TQEPC(N,K) = TQEPCIN(N,K)
!NGUSEF           SQEPC(N,K) = SQEPCIN(N,K)
!NG          Else
!NGUSEF           TQEPC(N,K) = UF(I,J,K)
!NGUSEF           SQEPC(N,K) = VF(I,J,K)
!NG          End If
!NG  114   Continue
!NG  117 Continue
!
!     INTRODUCING ADVECTIVE B.C. FOR T AND S BY Quamrul QA 3/5/99
!
*********************************  Khan
      IF(ALPHA.GT.0.0) THEN
        CALL ALPHABC(DTI2)
      ELSE
        DO 236 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        DO 170 K=1,KBM1
        UF(IE,JE,K)= TBDRY(N,K)
        VF(IE,JE,K)= SBDRY(N,K)
  170   CONTINUE
  236   CONTINUE
      ENDIF
!
      DO 240 J=1,JM
      DO 240 K=1,KBM1
      DO 240 I=1,IM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)
 240  CONTINUE
      RETURN
************************************  Khan
!
!  FOR CONSERVATIVE TRACER
!
  80  CONTINUE
      DO 2235 N=1,NUMQBCTR
        ID=ITRQD(N)
        JD=JTRQD(N)
        IC=ITRQC(N)
        JC=JTRQC(N)
        DO 2235 K=1,KBM1
          IF (HYDTYPE.EQ.'INTERNAL') THEN
            IF (VQDIST(N,K).NE.0.0 .AND. QDIS(N).GT.0.0) THEN
              UF(IC,JC,K)=CDIS1(N)
            ELSE
              UF(IC,JC,K)=UF(ID,JD,K)
            END IF
          ELSE
            UF(IC,JC,K)=CDIS1(N)
          ENDIF
 2235 CONTINUE
      DO 2236 N=1,NUMEBCTR
        IE=ITRED(N)
        JE=JTRED(N)
        DO 2170 K=1,KBM1
          UF(IE,JE,K)= CBDRY1(N,K)
 2170  CONTINUE
 2236 CONTINUE
!  FOR DIFFUSER   (NKIM  5.8.97)
      Do 2115 N = 1, NUMDBCTR
        ID = ITRDD(N)
        JD = JTRDD(N)
        Do 2112 K = 1, KBM1
          If (VDDIST(N,K).NE.0.0.AND.QDIFF(N).GT.0.0) Then
            CDIF1(N,K) = CDIFF1(N)
          Else
            CDIF1(N,K) = CONC1(ID,JD,K)
          End If
 2112   Continue
 2115 Continue
!  END OF FIX   (NKIM  5.8.97)
      DO 2240 J=1,JM
        DO 2240 K=1,KBM1
          DO 2240 I=1,IM
            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 2240 CONTINUE
      RETURN
!
!******************************************************************
!
!  FOR SEDIMENT TRANSPORT
!
  81  CONTINUE
!
!qa   ******** RIVER SEDIMENT LOADINGS *******************
      Do 2111 N = 1, NUMQBCSE
        ID=ISEQD(N)
        JD=JSEQD(N)
        IC=ISEQC(N)
        JC=JSEQC(N)
        Do 2101 K = 1, KBM1
          IF (HYDTYPE.EQ.'INTERNAL') THEN
            If (VQDIST(N,K).NE.0.0.AND.QDIS(N).GT.0.0) Then
              UF(IC,JC,K) = CDIS(KS,N)
            Else
              UF(IC,JC,K) = UF(ID,JD,K)
            End If
          ELSE
            UF(IC,JC,K) = CDIS(KS,N)
!            UF(ID,JD,K) = CDIS(KS,N)
          ENDIF
 2101   Continue
 2111 Continue
!
!qa   ******** BOUNDARY SEDIMENT SPEC *********************
      DO 5230 N=1,NUMEBCSE
        IE=ISEED(N)
        JE=JSEED(N)
        DO 5230 K=1,KBM1
          UF(IE,JE,K)= CBDRY(KS,N,K)
 5230 CONTINUE
!
      DO 5241 J=1,JM
        DO 5241 K=1,KBM1
          DO 5241 I=1,IM
            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 5241 CONTINUE
!
!
!qa   ******** DIFFUSER SEDIMENT LOADINGS *Quamrul QA 5/12/99
!qa      IF(NUMDBCSE.GT.0)THEN
!qa       DO 7230 N=1,NUMDBCSE
!qa        ID=IDDSE(N)
!qa        JD=JDDSE(N)
!qa        DO 7230 K=1,KBM1
!qa          UF(ID,JD,K)= CSDIFF(KS,N)
!qa 7230 CONTINUE
!
!hli       DO 7241 K=1,KBM1
!hli        DO 7241 I=1,IM
!hli          DO 7241 J=1,JM
!hli            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
!hli 7241  CONTINUE
!qa      ENDIF
!
      RETURN
!
!******************************************************************
!
!  FOR CHEM TRANSPORT 
!
 181  CONTINUE
!
      Do 3111 N = 1, NUMQBCCH
        ID=ICHQD(N)
        JD=JCHQD(N)
        IC=ICHQC(N)
        JC=JCHQC(N)
        Do 3101 K = 1, KBM1
          IF (HYDTYPE.EQ.'INTERNAL') THEN
            If (VQDIST(N,K).NE.0.0.AND.QDIS(N).GT.0.0) Then
              UF(IC,JC,K) = PDIS(KS,N)
            Else
              UF(IC,JC,K) = UF(ID,JD,K)
            End If
          ELSE
            UF(IC,JC,K) = PDIS(KS,N)
          ENDIF
 3101   Continue
 3111 Continue
!
      DO 6230 N=1,NUMEBCCH
        IE=ICHED(N)
        JE=JCHED(N)
        DO 6230 K=1,KBM1
          UF(IE,JE,K)= PBDRY(KS,N,K)
 6230 CONTINUE
!
!     Diffuser Loading
!
!qa      Do 7111 N = 1, NUMDBCCH
!qa        ID=IDDCH(N)
!qa        JD=JDDCH(N)
!qa        Do 7101 K = 1, KBM1
!qa            UF(ID,JD,K) = PDIFF(KS,N)
!qa 7101   Continue
!qa 7111 Continue
!
      DO 6240 J=1,JM
        DO 6240 K=1,KBM1
          DO 6240 I=1,IM
            UF(I,J,K)=UF(I,J,K)*FSM(I,J)
 6240 CONTINUE
!
!
      RETURN
!
!******************************************************************
!
!-------- VERTICAL VELOCITY BOUNDARY CONDITIONS ------------------------
!
  50  CONTINUE
      DO 245 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 245 K=1,KBM1
      W(IC,JC,K)=W(ID,JD,K)
  245 CONTINUE
      DO 246 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 246 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        W(IE,JE,K)=W(IE-1,JE,K)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        W(IE,JE,K)=W(IE+1,JE,K)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        W(IE,JE,K)=W(IE,JE-1,K)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        W(IE,JE,K)=W(IE,JE+1,K)
      END IF
  246 CONTINUE
      DO 250 J=1,JM
      DO 250 K=1,KBM1
      DO 250 I=1,IM
      W(I,J,K)=W(I,J,K)*FSM(I,J)
      IF(H(I,J)+ETF(I,J).LE.WETMIN. ! NG10252011
     +   AND.QEVAP2D(I,J).GT.QPREC2D(I,J)) W(I,J,K)=0.0 ! NG: This needs to be checked, a)whether it should be I,J,1
                                                           ! b) whether the if statement is needed at all, and 
                                                           ! c) whether it should be ETF or ETB for example. 
 250  CONTINUE
      RETURN
!
!-------- Q2 & Q2L BOUNDARY CONDITIONS ---------------------------------
  60  CONTINUE
      DO 255 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      DO 255 K=1,KBM1
      UF(IC,JC,K)=UF(ID,JD,K)
      VF(IC,JC,K)=VF(ID,JD,K)
      L(IC,JC,K)=L(ID,JD,K)
      KM(IC,JC,K)=KM(ID,JD,K)
      KH(IC,JC,K)=KH(ID,JD,K)
      KQ(IC,JC,K)=KQ(ID,JD,K)
  255 CONTINUE
      DO 256 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      DO 256 K=1,KBM1
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        UF(IE,JE,K)=UF(IE-1,JE,K)
        VF(IE,JE,K)=VF(IE-1,JE,K)
        L(IE,JE,K)=L(IE-1,JE,K)
        KM(IE,JE,K)=KM(IE-1,JE,K)
        KH(IE,JE,K)=KH(IE-1,JE,K)
        KQ(IE,JE,K)=KQ(IE-1,JE,K)
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        UF(IE,JE,K)=UF(IE+1,JE,K)
        VF(IE,JE,K)=VF(IE+1,JE,K)
        L(IE,JE,K)=L(IE+1,JE,K)
        KM(IE,JE,K)=KM(IE+1,JE,K)
        KH(IE,JE,K)=KH(IE+1,JE,K)
        KQ(IE,JE,K)=KQ(IE+1,JE,K)
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE-1,K)
        VF(IE,JE,K)=VF(IE,JE-1,K)
        L(IE,JE,K)=L(IE,JE-1,K)
        KM(IE,JE,K)=KM(IE,JE-1,K)
        KH(IE,JE,K)=KH(IE,JE-1,K)
        KQ(IE,JE,K)=KQ(IE,JE-1,K)
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        UF(IE,JE,K)=UF(IE,JE+1,K)
        VF(IE,JE,K)=VF(IE,JE+1,K)
        L(IE,JE,K)=L(IE,JE+1,K)
        KM(IE,JE,K)=KM(IE,JE+1,K)
        KH(IE,JE,K)=KH(IE,JE+1,K)
        KQ(IE,JE,K)=KQ(IE,JE+1,K)
      END IF
  256 CONTINUE
      DO 300 J=1,JM
      DO 300 K=1,KB
      DO 300 I=1,IM
      UF(I,J,K)=UF(I,J,K)*FSM(I,J)                                      
      VF(I,J,K)=VF(I,J,K)*FSM(I,J)
      L (I,J,K)=L (I,J,K)*FSM(I,J)
      KM(I,J,K)=KM(I,J,K)*FSM(I,J)
      KH(I,J,K)=KH(I,J,K)*FSM(I,J)
      KQ(I,J,K)=KQ(I,J,K)*FSM(I,J)
      IF(H(I,J)+ETF(I,J).LE.WETMIN. ! NG10252011
     +   AND.QEVAP2D(I,J).GT.QPREC2D(I,J)) THEN ! REINITIALIZE DRY CELLS
                                                   ! only if evap>prec. This also needs to be checked more
       UF(I,J,K)=0.0
       Q2(I,J,K)=0.0
       Q2B(I,J,K)=0.0
       VF(I,J,K)=0.0
       Q2L(I,J,K)=0.0
       Q2LB(I,J,K)=0.0
       L (I,J,K)=0.0
       KM(I,J,K)=0.0
       KH(I,J,K)=0.0
       KQ(I,J,K)=0.0
      ENDIF
 300  CONTINUE
      RETURN
!
!-------- WAVE PARAMETERS AND TEMPORAL CYCLING (hli, 03/19/04)----------
  90  CONTINUE
!NG02222007      IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN      
      IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN      
        IF(NUMEBC.NE.0.AND.CREADWBC) THEN  
          IF(THOUR.GE.T2WAVE) THEN        
            T1WAVE=T2WAVE 
            DO N=1,NUMEBC
              DWHBFLX(N,1)=DWHBFLX(N,2) 
              DWPBFLX(N,1)=DWPBFLX(N,2) 
              DWDBFLX(N,1)=DWDBFLX(N,2) 
            ENDDO
            READ(IUT195,END=4021) T2WAVE
            READ(IUT195)(DWHBFLX(N,2),N=1,NUMEBC)
            READ(IUT195)(DWPBFLX(N,2),N=1,NUMEBC)
            READ(IUT195)(DWDBFLX(N,2),N=1,NUMEBC)
          ENDIF       
          FACT=(THOUR-T1WAVE)/(T2WAVE-T1WAVE)
          DO N=1,NUMEBC
!NG         Do direction averaging correctly below...
            WAVE1XTEMP=DWHBFLX(N,1)*COS(DWDBFLX(N,1)*TPI/360.)
            WAVE1YTEMP=DWHBFLX(N,1)*SIN(DWDBFLX(N,1)*TPI/360.)
            WAVE2XTEMP=DWHBFLX(N,2)*COS(DWDBFLX(N,2)*TPI/360.)
            WAVE2YTEMP=DWHBFLX(N,2)*SIN(DWDBFLX(N,2)*TPI/360.)
            WAVEXTEMP=WAVE1XTEMP+FACT*(WAVE2XTEMP-WAVE1XTEMP)
            WAVEYTEMP=WAVE1YTEMP+FACT*(WAVE2YTEMP-WAVE1YTEMP)
            WDBFLX(N)=ATAN2(WAVEYTEMP,WAVEXTEMP)/TPI*360.
            IF(WDBFLX(N).LT.0.) WDBFLX(N)=WDBFLX(N)+360
!NG         Do direction averaging correctly above...
            WHBFLX(N)=DWHBFLX(N,1)+FACT*(DWHBFLX(N,2)-DWHBFLX(N,1)) 
            WPBFLX(N)=DWPBFLX(N,1)+FACT*(DWPBFLX(N,2)-DWPBFLX(N,1)) 
            WHBFLX(N)=RAMP*WHBFLX(N)                                 
            WPBFLX(N)=RAMP*WPBFLX(N)
          ENDDO
!          write(1005,*) 'time=',time,thour,t1wave,t2wave
!          write(1005,'(8f10.5)') (WPBFLX(N),n=1,numebc)
!          write(1005,'(8f10.5)') (WPBFLX(N),n=1,numebc)
!          write(1005,'(8f10.5)') (WDBFLX(N),JETA(N)),n=1,numebc)
        ENDIF      
      ENDIF      
      RETURN
! 
!
!-------- SURFACE FORCING AND TEMPORAL CYCLING -------------------------
  70  CONTINUE
!
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 4044
! 
!-------- ELEVATION BOUNDARY CONDITIONS --------------------------------
      IF(NUMEBC.EQ.0) GOTO 4045
!     IF(OPTEBC(1:5).EQ.'DATA ') THEN  ! NO NEED TO HAVE THIS ANYMORE
       IF(THOUR.LT.T2E) GOTO 4020
       T1E=T2E 
       DO 4030 N=1,NUMEBC
 4030  DEBDRY(N,1)=DEBDRY(N,2) 
       ! A.Jordi: Read elevation boundary condition from netcdf file
       NT_T2E=NT_T2E+1
       CALL NETCDF_READ_BC_ELEV(NT_T2E,M_MAX,STR)
       ! A.Jordi: end
 4020  CONTINUE
       FACT=(THOUR-T1E)/(T2E-T1E)
       DO 4040 N=1,NUMEBC
 4040  EBDRY(N)=DEBDRY(N,1)+FACT*(DEBDRY(N,2)-DEBDRY(N,1)) 
!A    END IF    ! NO NEED TO HAVE THIS ANYMORE (N.KIM 121798)
 4045 CONTINUE      
! 
!NG04232008 ELEVATION POTENTIAL CELLS ----------------------------------
      IF(NUMEPC.EQ.0) GOTO 6545
       IF(THOUR.LT.T2EPC) GOTO 6520
       T1EPC=T2EPC 
       DO 6530 N=1,NUMEPC
 6530  DEPCOR(N,1)=DEPCOR(N,2) 
       READ (IUT65,4000,END=6510) T2EPC
       READ (IUT65,4000) (DEPCOR(N,2),N=1,NUMEPC)
 6520  CONTINUE
       FACT=(THOUR-T1EPC)/(T2EPC-T1EPC)
       DO 6540 N=1,NUMEPC
        I = IEPC(N)
        J = JEPC(N)
        IF ((DEPCOR(N,2).NE.-99.).AND.(DEPCOR(N,1).NE.-99)) THEN
         EPCOR(N)=DEPCOR(N,1)+FACT*(DEPCOR(N,2)-DEPCOR(N,1)) 
        ELSE
         EPCOR(N)=-99. ! FLAG FOR FREE SOLUTION
        ENDIF
!NG2D        TQEPCIN(N,1)= 0.0
!NG2D        SQEPCIN(N,1) = 0.0
!NG          DO K=1,KBM1
!NG2D            TQEPCIN(N,1)=TQEPCIN(N,1)+T(I,J,K)*DZ(K)
!NG2D            SQEPCIN(N,1)=SQEPCIN(N,1)+S(I,J,K)*DZ(K)
!NGUSEF            TQEPCIN(N,K)=T(I,J,K)
!NGUSEF            SQEPCIN(N,K)=S(I,J,K)
!NG          END DO
 6540  CONTINUE
 6545 CONTINUE
! 
!-------- TEMPERATURE AND SALINITY BOUNDARY CONDITIONS -----------------
! 
      IF(NUMEBC.EQ.0) GOTO 4135
      IF(THOUR.LT.T2TS) GOTO 4111
      T1TS=T2TS
      DO 4122 K=1,KBM1
      DO 4122 N=1,NUMEBC
      DTBDRY(N,K,1)=DTBDRY(N,K,2)
      DSBDRY(N,K,1)=DSBDRY(N,K,2)
4122  CONTINUE
      ! A.Jordi: Read elevation boundary condition from netcdf file
      NT_T2TS=NT_T2TS+1
      CALL NETCDF_READ_BC_TS(NT_T2TS,M_MAX,STR)
      ! A.Jordi: end
4111  CONTINUE
      FACT=(THOUR-T1TS)/(T2TS-T1TS)
      DO 4042 K=1,KBM1
      DO 4042 N=1,NUMEBC
      TBDRY(N,K)=DTBDRY(N,K,1)+FACT*(DTBDRY(N,K,2)-DTBDRY(N,K,1)) 
      SBDRY(N,K)=DSBDRY(N,K,1)+FACT*(DSBDRY(N,K,2)-DSBDRY(N,K,1)) 
 4042  CONTINUE
! 
!  DISSOLVED TRACER AT OPEN BOUNDARY
! 
 4044 IF (TRACER.EQ.'INCLUDE'.AND.NUMEBCTR.GT.0) THEN
        IF(THOUR.LT.T2CONOB) GOTO 6111
!
        T1CONOB=T2CONOB
!
        DO 6121 K=1,KBM1
          DO 6121 N=1,NUMEBCTR
            DCBDRY1(N,K,1)=DCBDRY1(N,K,2)
 6121   CONTINUE
!
        READ (IUT501,4000,END=7011) T2CONOB
        DO 6136 N=1,NUMEBCTR
          READ (IUT501,4000) (DCBDRY1(N,K,2),K=1,KBM1)
 6136   CONTINUE
!
 6111   CONTINUE
!
        FACT=(THOUR-T1CONOB)/(T2CONOB-T1CONOB)
!
        DO 6042 K=1,KBM1
          DO 6042 N=1,NUMEBCTR
            CBDRY1(N,K)=DCBDRY1(N,K,1)+FACT*(DCBDRY1(N,K,2)-
     +                                          DCBDRY1(N,K,1)) 
 6042   CONTINUE
      ENDIF
! 
!  SEDIMENT TRANSPORT AT OPEN BOUNDARY
! 
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMEBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
!
          IF(THOUR.LT.T2SEDOBM) GOTO 6010
!
          T1SEDOBM=T2SEDOBM
!
          DO 6020 K=1,KBM1
            DO 6020 N=1,NUMEBCSE
              DCBDRY(KK,N,K,1)=DCBDRY(KK,N,K,2)
 6020     CONTINUE
!
          READ (IUT502,4000,END=7012) T2SEDOBM
          DO 6030 N=1,NUMEBCSE
            READ (IUT502,4000) (DCBDRY(KK,N,K,2),K=1,KBM1)
 6030     CONTINUE
!
 6010     CONTINUE
!
          FACT=(THOUR-T1SEDOBM)/(T2SEDOBM-T1SEDOBM)
!
          DO 6040 K=1,KBM1
            DO 6040 N=1,NUMEBCSE
              CBDRY(KK,N,K)=DCBDRY(KK,N,K,1)+FACT*(DCBDRY(KK,N,K,2)-
     +                                          DCBDRY(KK,N,K,1)) 
 6040     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
!
          IF(THOUR.LT.T2SEDOBS) GOTO 6050
!
          T1SEDOBS=T2SEDOBS
!
          DO 6060 K=1,KBM1
            DO 6060 N=1,NUMEBCSE
              DCBDRY(KK,N,K,1)=DCBDRY(KK,N,K,2)
 6060     CONTINUE
!
          READ (IUT503,4000,END=7212) T2SEDOBS
          DO 6070 N=1,NUMEBCSE
            READ (IUT503,4000) (DCBDRY(KK,N,K,2),K=1,KBM1)
 6070     CONTINUE
!
 6050     CONTINUE
!
          FACT=(THOUR-T1SEDOBS)/(T2SEDOBS-T1SEDOBS)
!
          DO 6080 K=1,KBM1
            DO 6080 N=1,NUMEBCSE
              CBDRY(KK,N,K)=DCBDRY(KK,N,K,1)+FACT*(DCBDRY(KK,N,K,2)-
     +                                          DCBDRY(KK,N,K,1)) 
 6080     CONTINUE
        ENDIF
      ENDIF
! 
!  PARTICLE-BOUND TRACER AT OPEN BOUNDARY
! 
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMEBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
!
          IF(THOUR.LT.T2CHMOBM) GOTO 6110
!
          T1CHMOBM=T2CHMOBM
!
          DO 6120 K=1,KBM1
            DO 6120 N=1,NUMEBCCH
              DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 6120     CONTINUE
!
          READ (IUT504,4000,END=7013) T2CHMOBM
          DO 6130 N=1,NUMEBCCH
            READ (IUT504,4000) (DPBDRY(KK,N,K,2),K=1,KBM1)
 6130     CONTINUE
!
 6110     CONTINUE
!
          FACT=(THOUR-T1CHMOBM)/(T2CHMOBM-T1CHMOBM)
!
          DO 6140 K=1,KBM1
            DO 6140 N=1,NUMEBCCH
              PBDRY(KK,N,K)=DPBDRY(KK,N,K,1)+FACT*(DPBDRY(KK,N,K,2)-
     +                                          DPBDRY(KK,N,K,1)) 
 6140     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
!
          IF(THOUR.LT.T2CHMOBS) GOTO 6150
!
          T1CHMOBS=T2CHMOBS
!
          DO 6160 K=1,KBM1
            DO 6160 N=1,NUMEBCCH
              DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 6160     CONTINUE
!
          READ (IUT505,4000,END=7213) T2CHMOBS
          DO 6170 N=1,NUMEBCCH
            READ (IUT505,4000) (DPBDRY(KK,N,K,2),K=1,KBM1)
 6170     CONTINUE
!
 6150     CONTINUE
!
          FACT=(THOUR-T1CHMOBS)/(T2CHMOBS-T1CHMOBS)
!
          DO 6180 K=1,KBM1
            DO 6180 N=1,NUMEBCCH
              PBDRY(KK,N,K)=DPBDRY(KK,N,K,1)+FACT*(DPBDRY(KK,N,K,2)-
     +                                          DPBDRY(KK,N,K,1)) 
 6180     CONTINUE
        ENDIF
      ENDIF
!
 4135  CONTINUE
!
!-------- RIVER DISCHARGE BOUNDARY CONDITIONS --------------------------
!
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 4074
!
      IF(NUMQBC.EQ.0) GOTO 4075
!
      IF(THOUR.LT.T2Q) GOTO 4050
      T1Q=T2Q 
      DO 4060 N=1,NUMQBC
      DQDIS(N,1)=DQDIS(N,2) 
      DTDIS(N,1)=DTDIS(N,2) 
      DSDIS(N,1)=DSDIS(N,2) 
 4060 CONTINUE
!
      ! A.Jordi: Read diffuser 1 boundary condition from netcdf file
      NT_T2Q=NT_T2Q+1
      CALL NETCDF_READ_BC_DIS(NT_T2Q,NT_T2Q_MAX,STR)
      ! A.Jordi: end
!
 4050 CONTINUE
      FACT=(THOUR-T1Q)/(T2Q-T1Q)
      DO 4070 N=1,NUMQBC
      QDIS(N)=DQDIS(N,1)+FACT*(DQDIS(N,2)-DQDIS(N,1)) 
      TDIS(N)=DTDIS(N,1)+FACT*(DTDIS(N,2)-DTDIS(N,1)) 
      SDIS(N)=DSDIS(N,1)+FACT*(DSDIS(N,2)-DSDIS(N,1)) 
 4070 CONTINUE
!
!  CONSERVATIVE TRACER
!
 4074 IF (TRACER.EQ.'INCLUDE'.AND.NUMQBCTR.GT.0) THEN
        IF (THOUR.LT.T2CON) GOTO 5050
!
        T1CON=T2CON 
!
        DO 5060 N=1,NUMQBCTR
          DCDIS1(N,1)=DCDIS1(N,2)
 5060   CONTINUE
!
        READ (IUT601,4000,END=5012) T2CON
        READ (IUT601,4000) (DCDIS1(N,2),N=1,NUMQBCTR)
!
 5050   CONTINUE
!
        FACT=(THOUR-T1CON)/(T2CON-T1CON)
!
        DO 5070 N=1,NUMQBCTR
          CDIS1(N)=DCDIS1(N,1)+FACT*(DCDIS1(N,2)-DCDIS1(N,1)) 
 5070   CONTINUE
      ENDIF
!
!  SEDIMENT TRANSPORT   RIVER DISCHARGE
!
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMQBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2SEDM) GOTO 5155
          T1SEDM=T2SEDM 

          DO 5145 N=1,NUMQBCSE
            DCDIS(KK,N,1)=DCDIS(KK,N,2)
 5145     CONTINUE
!
          READ (IUT602,4000,END=5213) T2SEDM
          READ (IUT602,4000) (DCDIS(KK,N,2),N=1,NUMQBCSE)
!
 5155     CONTINUE
!
          FACT=(THOUR-T1SEDM)/(T2SEDM-T1SEDM)
!
          DO 5165 N=1,NUMQBCSE
            CDIS(KK,N)=DCDIS(KK,N,1)+FACT*(DCDIS(KK,N,2)
     +                                  -DCDIS(KK,N,1))
 5165     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
!
          IF (THOUR.LT.T2SEDS) GOTO 5156
          T1SEDS=T2SEDS 
!
          DO 5146 N=1,NUMQBCSE
            DCDIS(KK,N,1)=DCDIS(KK,N,2)
 5146     CONTINUE
!
          READ (IUT603,4000,END=5313) T2SEDS
          READ (IUT603,4000) (DCDIS(KK,N,2),N=1,NUMQBCSE)
!
 5156     CONTINUE
!
          FACT=(THOUR-T1SEDS)/(T2SEDS-T1SEDS)
!
          DO 5166 N=1,NUMQBCSE
            CDIS(KK,N)=DCDIS(KK,N,1)+FACT*(DCDIS(KK,N,2)
     +                                  -DCDIS(KK,N,1))
 5166     CONTINUE
        ENDIF
      ENDIF
!
!  PARTICLE-BOUND TRACER TRANSPORT
!
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMQBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2CHMM) GOTO 5250
          T1CHMM=T2CHMM 
!
          DO 5240 N=1,NUMQBCCH
            DPDIS(KK,N,1)=DPDIS(KK,N,2)
 5240     CONTINUE
!
          READ (IUT604,4000,END=5014) T2CHMM
          READ (IUT604,4000) (DPDIS(KK,N,2),N=1,NUMQBCCH)
!
 5250     CONTINUE
!
          FACT=(THOUR-T1CHMM)/(T2CHMM-T1CHMM)
!
          DO 5260 N=1,NUMQBCCH
            PDIS(KK,N)=DPDIS(KK,N,1)+FACT*(DPDIS(KK,N,2)
     +                                  -DPDIS(KK,N,1))
 5260     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
!
          IF (THOUR.LT.T2CHMS) GOTO 5255
          T1CHMS=T2CHMS 
          DO 5245 N=1,NUMQBCCH
            DPDIS(KK,N,1)=DPDIS(KK,N,2)
 5245     CONTINUE
          READ (IUT605,4000,END=5214) T2CHMS
          READ (IUT605,4000) (DPDIS(KK,N,2),N=1,NUMQBCCH)
 5255     CONTINUE
!
          FACT=(THOUR-T1CHMS)/(T2CHMS-T1CHMS)
!
          DO 5265 N=1,NUMQBCCH
            PDIS(KK,N)=DPDIS(KK,N,1)+FACT*(DPDIS(KK,N,2)
     +                                  -DPDIS(KK,N,1))
 5265     CONTINUE
        ENDIF
      ENDIF
!
 4075 CONTINUE
!
! add 'DIFFUSER IN LOOP' to the codes, referring to TAMPA BAY codes
! in /power6/wrai0010/MODEL/CODES/TRACER (hli, 09/28/00)
! 
!-------- DIFFUSER INTAKE/OUTFALL CONDITIONS ---------------------------
!         (REGULAR INTAKE/OUTFALLS)
! 
!PL ADD NEXT LINE 01/11/00
!
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 4104
!
      IF(NUMDBC1.EQ.0) GOTO 4105
      IF(THOUR.LT.T2D) GOTO 4080
      T1D=T2D
      DO N=1,NUMDBC1
        DQDIFF(N,1)=DQDIFF(N,2) 
        DTDIFF(N,1)=DTDIFF(N,2) 
        DSDIFF(N,1)=DSDIFF(N,2) 
      ENDDO   
!
      ! A.Jordi: Read diffuser 1 boundary condition from netcdf file
      NT_T2D=NT_T2D+1
      CALL NETCDF_READ_BC_DIFF1(NT_T2D,NT_T2D_MAX,STR)
      ! A.Jordi: end
!
 4080 CONTINUE
      FACT=(THOUR-T1D)/(T2D-T1D)
      DO N=1,NUMDBC1
        QDIFF(N)=DQDIFF(N,1)+FACT*(DQDIFF(N,2)-DQDIFF(N,1)) 
        TDIFF(N)=DTDIFF(N,1)+FACT*(DTDIFF(N,2)-DTDIFF(N,1)) 
        SDIFF(N)=DSDIFF(N,1)+FACT*(DSDIFF(N,2)-DSDIFF(N,1)) 
      ENDDO     
 4105 CONTINUE
!
!-------- DIFFUSER INTAKE/OUTFALL CONDITIONS IN LOOPS-------------------------
!           (COUPLED)  (4/30/97  NKIM)
!
      If (NUMDBC2.EQ.0) GO TO 640
      If (THOUR.GE.T2D2) Then
        T1D2 = T2D2
        DO N = NUMDBC1+1, NUMDBC1+NUMDBC2
          DQDIFF(N,1) = DQDIFF(N,2)
          DTDIFF(N,1) = DTDIFF(N,2)
          DSDIFF(N,1) = DSDIFF(N,2)
        END DO
        ! A.Jordi: Read diffuser 2 boundary condition from netcdf file
        NT_T2D2=NT_T2D2+1
        CALL NETCDF_READ_BC_DIFF2(NT_T2D2,NT_T2D2_MAX,STR)
        ! A.Jordi: end
      END IF
!
      FACT = (THOUR-T1D2) / (T2D2-T1D2)
      Do N = NUMDBC1+1, NUMDBC1+NUMDBC2,2
        ID= IDD(N)
        JD= JDD(N)
        QDIFF(N)  = DQDIFF(N,1) + FACT * (DQDIFF(N,2)-DQDIFF(N,1))
        QDIFF(N+1)= DQDIFF(N+1,1) +
     .              FACT * (DQDIFF(N+1,2)-DQDIFF(N+1,1))
        AVGTEMP= 0.0
        AVGSAL = 0.0
        VOLUME = 0.0
        DO K=1,KBM1
!qa          WRITE(55,'(20F8.3)')T(ID,JD,K),S(ID,JD,K)
          AVGTEMP=AVGTEMP+T(ID,JD,K)*DZ(K)*VDDIST(N,K)
          AVGSAL =AVGSAL +S(ID,JD,K)*DZ(K)*VDDIST(N,K)
          VOLUME=VOLUME + DZ(K)*VDDIST(N,K)
        END DO
        IF(VOLUME.EQ.0.0)GO TO  4516
        TDIFF(N)=0.0
        SDIFF(N)=0.0
        TDIFF(N+1)=AVGTEMP/VOLUME + DTDIFF(N+1,1)
     .            + FACT * (DTDIFF(N+1,2)-DTDIFF(N+1,1))
        SDIFF(N+1)=AVGSAL/VOLUME  + DSDIFF(N+1,1)
     .            + FACT * (DSDIFF(N+1,2)-DSDIFF(N+1,1))
!
      ENDDO    
!hli      IF(AMOD(THOUR,1.00).EQ.0.0)THEN  ! WRITE HOURLY SALINITY OF DIFFUSER
!hli        WRITE(55,'(20f10.3)')TIME,QDIFF(1),QDIFF(2),TDIFF(1),TDIFF(2),
!hli     .   T(IDD(1),JDD(1),1),SDIFF(1),SDIFF(2),S(IDD(1),JDD(1),1)
!hli        WRITE(55,'(20F10.3)')AVGTEMP,VOLUME,DTDIFF(2,1),DTDIFF(2,2)
!hli        WRITE(55,'(20F10.3)')AVGSAL,VOLUME,DSDIFF(2,1),DSDIFF(2,2)
!hli      END IF
!
  640 CONTINUE
!
!  TRACER FOR DIFFUSERS ARE COMBINED FOR REGULAR AND LOOP DIFFUSER CASES.
!  USER CAN SPECIFY ANY DIFFUSERS FOR TRACER INPUT (EITHER REGULAR OR LOOP).
!  THE TRACER CONC. SPECIFIED FOR LOOP DIFFUSERS ARE NOT TRACER INCREMENTS
!  AS USED IN SALINITY AND TEMPERATURE CASES.  (5/10/99)
!
! THE CODE HAS BEEN MODIFIED TO SEPARATE REGULAR DIFFUSER AND DIFFUSER-IN-LOOP
! FOR TRACERS   (NKIM 7.25.00)
!
4104  IF (TRACER.NE.'INCLUDE') GO TO 650
      IF (NUMDBCTR.EQ.0)GO TO 5102
      IF (NUMDBCTR1.EQ.0)GO TO 5101
      IF (THOUR.LT.T2DCON) GOTO 5080
!
      T1DCON=T2DCON
!
      DO N=1,NUMDBCTR1
        DCDIFF1(N,1)=DCDIFF1(N,2)
      ENDDO    
!
      READ (IUT98,4000) T2DCON
      READ (IUT98,4000) (DCDIFF1(N,2),N=1,NUMDBCTR1)
!
 5080 CONTINUE
!
      FACT=(THOUR-T1DCON)/(T2DCON-T1DCON)
!
      DO N=1,NUMDBCTR1
        CDIFF1(N)=DCDIFF1(N,1)+FACT*(DCDIFF1(N,2)-DCDIFF1(N,1))
      ENDDO    
 5101 CONTINUE
!
!  TRACER FOR DIFFUSER-IN-LOOPS  (NKIM  7/25/00)
      IF(NUMDBCTR2.EQ.0) GO TO 5102
!
      IF (THOUR.LT.T2D2CON) GOTO 580
!
      T1D2CON=T2D2CON
!
      DO N=NUMDBCTR1+1,NUMDBCTR
        DCDIFF1(N,1)=DCDIFF1(N,2)
      ENDDO 
!
      READ (IUT99,4000) T2D2CON
      READ (IUT99,4000) (DCDIFF1(N,2),N=NUMDBCTR1+1,NUMDBCTR)
!
 580  CONTINUE
!
!qa
!     Diffuser in Loop for Flows and Tracers are same so IDD,JDD are used
!     we dont use Delta Conc like Delta T but we use Average Conce 
!qa
      FACT=(THOUR-T1D2CON)/(T2D2CON-T1D2CON)
      Do N = NUMDBCTR1+1, NUMDBCTR,2
        ID= IDD(N)
        JD= JDD(N)
        AVGTR= 0.0
        VOLUME = 0.0
        DO K=1,KBM1
          AVGTR=AVGTR+CONC1(ID,JD,K)*DZ(K)*VDDIST(N,K)
          VOLUME=VOLUME + DZ(K)*VDDIST(N,K)
        END DO
        IF(VOLUME.EQ.0.0)GO TO  4516
        CDIFF1(N)=0.0    ! INTAKE CONC.
        CDIFF1(N+1)=AVGTR/VOLUME + DCDIFF1(N+1,1)
     .             + FACT * (DCDIFF1(N+1,2)-DCDIFF1(N+1,1))

      ENDDO       
!
!qa      WRITE(55,'(20E12.4)')TIME,CDIFF1(1),CDIFF1(2)
!
 5102 CONTINUE
!
!  end of adding DIFFUSER IN LOOP (hli, 09/28/00)
!
!  ******************************************************************
!  PARTICLE-BOUND TRACER LOADINGS THRU DIFFUSER By Quamrul QA 5/13/99
!
 650  IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMDBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2CHMDM) GOTO 7250
          T1CHMDM=T2CHMDM
!
          DO 7240 N=1,NUMDBCCH
            DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7240     CONTINUE
          READ (IUT704,4000,END=7014) T2CHMDM
          READ (IUT704,4000) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
 7250     CONTINUE
!
          FACT=(THOUR-T1CHMDM)/(T2CHMDM-T1CHMDM)
!
          DO 7260 N=1,NUMDBCCH
            PDIFF(KK,N)=DPDIFF(KK,N,1)+FACT*(DPDIFF(KK,N,2)
     +                                  -DPDIFF(KK,N,1))
 7260     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
!
          IF (THOUR.LT.T2CHMDS) GOTO 7255
          T1CHMDS=T2CHMDS
!
          DO 7245 N=1,NUMDBCCH
            DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7245     CONTINUE
          READ (IUT705,4000,END=5214) T2CHMDS
          READ (IUT705,4000) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
 7255     CONTINUE
!
          FACT=(THOUR-T1CHMDS)/(T2CHMDS-T1CHMDS)
          DO 7265 N=1,NUMDBCCH
            PDIFF(KK,N)=DPDIFF(KK,N,1)+FACT*(DPDIFF(KK,N,2)
     +                                  -DPDIFF(KK,N,1))
 7265     CONTINUE
        ENDIF
      ENDIF
      CONTINUE
!
!  POINT SOURCE LOADS:  CONSERVATIVE TRACER
!
      IF (TRACER.EQ.'INCLUDE'.AND.NUMPSTR.GT.0) THEN
        IF (THOUR.LT.T2PSTR) GOTO 5051
!
!  INSTANTANEOUS SPECIFICATION OF PT. SOURCE CONCENTRATION
!
        IF (OPTPSTR.EQ.'CONC') THEN
          DO 5062 N=1,NUMPSTR
            CPSTR(N)=DCPSTR(N,2)
 5062     CONTINUE
        ENDIF
!
        T1PSTR=T2PSTR 
!
        DO 5061 N=1,NUMPSTR
          DCPSTR(N,1)=DCPSTR(N,2)
 5061   CONTINUE
!
        READ (IUT701,4000,END=5012) T2PSTR
        READ (IUT701,4000) (DCPSTR(N,2),N=1,NUMPSTR)
!
 5051   CONTINUE
!
!  INSTANTANEOUS SPECIFICATION OF PT. SOURCE CONCENTRATION
!
        IF (OPTPSTR.EQ.'CONC') GOTO 5081
!
        IF (OPTPSTR.EQ.'MASS') THEN
          FACT=(THOUR-T1PSTR)/(T2PSTR-T1PSTR)
!
          DO 5071 N=1,NUMPSTR
            CPSTR(N)=DCPSTR(N,1)+FACT*(DCPSTR(N,2)-DCPSTR(N,1)) 
 5071     CONTINUE
        ENDIF
      ENDIF
!
!     *********************************************************
!     DIFFUSER SEDIMENT LOADINGS  By Quamrul QA 5/11/99
!     SEDIMENT TRANSPORT   DIFFUSER DISCHARGE
!
 5081  IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMDBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          KK=1
          IF (THOUR.LT.T2SEDDM) GOTO 7155
          T1SEDDM=T2SEDDM

          DO 7145 N=1,NUMDBCSE
            DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7145     CONTINUE
!
          READ (IUT702,4000,END=7313) T2SEDDM
          READ (IUT702,4000) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
!
 7155     CONTINUE
!
          FACT=(THOUR-T1SEDDM)/(T2SEDDM-T1SEDDM)
!
          DO 7165 N=1,NUMDBCSE
            CSDIFF(KK,N)=DCSDIFF(KK,N,1)+FACT*(DCSDIFF(KK,N,2)
     +                                  -DCSDIFF(KK,N,1))
 7165     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IF (SEDTYPE.EQ.'SAND') THEN
            KK=1
          ELSE
            KK=2
          ENDIF
!
          IF (THOUR.LT.T2SEDDS) GOTO 7156
          T1SEDDS=T2SEDDS
!
          DO 7146 N=1,NUMDBCSE
            DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7146     CONTINUE
!
          READ (IUT703,4000,END=7413) T2SEDDS
          READ (IUT703,4000) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
!
 7156     CONTINUE
!
          FACT=(THOUR-T1SEDDS)/(T2SEDDS-T1SEDDS)
!
          DO 7166 N=1,NUMDBCSE
            CSDIFF(KK,N)=DCSDIFF(KK,N,1)+FACT*(DCSDIFF(KK,N,2)
     +                                  -DCSDIFF(KK,N,1))
 7166     CONTINUE
        ENDIF
       ENDIF

!
!-------- METEOROLOGICAL BOUNDARY CONDITIONS ---------------------------
!NG09022014 First zero out 2D precipitation and evaporation.
!NG09022014 Then use synop_pet if it exists, or 1D if it does not and available
! In the uniform forcing path 4115 overwrites both fields on *all*
! cells. Keep the old clearing when PET data are used, another mode is
! active, or a forcing record must be read (including READ/END exits).
        CLEAR_RAIN=.TRUE.
        IF(OPTMBC(1:8).EQ.'AVERAGED'.AND..NOT.CREADPET) THEN
          IF(THOUR.LT.T2M) CLEAR_RAIN=.FALSE.
        ENDIF
        IF(CLEAR_RAIN) THEN
        DO J=1,JM
         DO I=1,IM
          QPREC2D(I,J)=0.0
          QEVAP2D(I,J)=0.0
         ENDDO
        ENDDO
        ENDIF
!NG09022014 Precipitation and Evapotranspiration from synop_pet... if it exists...
        IF (CREADPET) THEN
         IF (THOUR.GE.T2P) THEN
          T1P = T2P
          DO 4420 J=1,JM
           DO 4420 I=1,IM
            IF(FSM(I,J).EQ.0)GO TO 4420
            DQPREC2D(I,J,1) = DQPREC2D(I,J,2)
            DQEVAP2D(I,J,1) = DQEVAP2D(I,J,2)
 4420     CONTINUE
!   UPDATE DATA FROM 'synop_pet' data
          READ (IUT104,End=4026) T2P
          READ(IUT104) ((DQPREC2D(I,J,2),DQEVAP2D(I,J,2),I=1,IM),J=1,JM)
         END IF
         FACT = (THOUR-T1P) / (T2P-T1P)
         DO 4390 J=1,JM
          DO 4390 I=1,IM
           IF(FSM(I,J).EQ.0.0)GO TO 4390
!NG04192011 FOR PRECIPITATION AND EVAPORATION, CHECK FOR NET EVAPORATION AND TREAT AS 2D
           QPREC2D(I,J) = DQPREC2D(I,J,1) 
     .           + FACT*(DQPREC2D(I,J,2)-DQPREC2D(I,J,1))  ! PRECIPITATION
           QPREC2D(I,J)=QPREC2D(I,J)/(86400.*365.) ! NG12162010 m/year to m/sec as in bcdata
           QEVAP2D(I,J) = DQEVAP2D(I,J,1) 
     .           + FACT*(DQEVAP2D(I,J,2)-DQEVAP2D(I,J,1))  ! EVAPORATION
           QEVAP2D(I,J)=QEVAP2D(I,J)/(86400.*365.) ! NG12162010 m/year to m/sec as in bcdata
           IF (DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).LT.QEVAP2D(I,J)) THEN ! DO NOT ALLOW NET EVAP ON DRY CELL
            QPREC2D(I,J)=0.0
            QEVAP2D(I,J)=0.0
           ENDIF
 4390    CONTINUE
        ENDIF
!NG09022014 ...synop_pet
!
!-------- AVERAGED OPTION ----------------------------------------------
      IF(OPTMBC(1:8).EQ.'AVERAGED') THEN
       IF(THOUR.LT.T2M) GOTO 4110
       T1M=T2M
       DQPREC(1)=DQPREC(2) 
       DQEVAP(1)=DQEVAP(2) 
       DHFLUX(1)=DHFLUX(2)
       DTX   (1)=DTX   (2) 
       DTY   (1)=DTY   (2) 
       DWSX   (1)=DWSX   (2) 
       DWSY   (1)=DWSY   (2) 
       DWDS   (1)=DWDS   (2) 
       DWDD   (1)=DWDD   (2) 
       READ (IUT93,4000,END=4016) T2M
       READ (IUT93,4000) DHFLUX(2),DTX(2),DTY(2),
     +      DWSX(2),DWSY(2),DWDS(2),DWDD(2),DQPREC(2),DQEVAP(2)
 4110  CONTINUE
       FACT=(THOUR-T1M)/(T2M-T1M)
       QPREC=DQPREC(1)+FACT*(DQPREC(2)-DQPREC(1))
       QEVAP=DQEVAP(1)+FACT*(DQEVAP(2)-DQEVAP(1))
       HFLUX=DHFLUX(1)+FACT*(DHFLUX(2)-DHFLUX(1))
       TX   =DTX   (1)+FACT*(DTX   (2)-DTX   (1)) 
       TY   =DTY   (1)+FACT*(DTY   (2)-DTY   (1)) 
!      Quamrul QA 10/17/98
       UWIND   =DWSX   (1)+FACT*(DWSX   (2)-DWSX   (1)) 
       VWIND   =DWSY   (1)+FACT*(DWSY   (2)-DWSY   (1)) 
!      End Quamrul QA 10/17/98
!NG2014
       IF(CREADPET) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
        DO 4114 J=1,JM
         DO 4114 I=1,IM
           WSX2D(I,J) = UWIND
           WSY2D(I,J) = VWIND
4114    CONTINUE
       ELSE
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
        DO 4115 J=1,JM
         DO 4115 I=1,IM
!NG04092008           WU(I,J) = UWIND
!NG04092008           WV(I,J) = VWIND
           WSX2D(I,J) = UWIND
           WSY2D(I,J) = VWIND
!NG04192011 FOR PRECIPITATION AND EVAPORATION, CHECK FOR NET EVAPORATION AND TREAT AS 2D
           QPREC2D(I,J) = FSM(I,J)*QPREC
           QEVAP2D(I,J) = FSM(I,J)*QEVAP
           IF (DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).LT.QEVAP2D(I,J)) THEN ! DO NOT ALLOW NET EVAP ON DRY CELL
            QEVAP2D(I,J)=0.0
            QPREC2D(I,J)=0.0
           ENDIF
4115    CONTINUE
       ENDIF
       DO N=1,NUMEBC ! CNG10032011
        IE=IETA(N)
        JE=JETA(N)
        QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
        QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
       ENDDO
       DO N=1,NUMECR ! CNG2012
        IE=ICHARD(N)
        JE=JCHARD(N)
        QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
        QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
       ENDDO

!
!  following is wrong Quamrul 10/27/98
!pl       WINDSP=RAMP*(DWDS(1)-FACT*(DWDS(2)-DWDS(1)))
!qa       WINDSP=RAMP*(DWDS(1)+FACT*(DWDS(2)-DWDS(1)))
!qa       WINDIR=DWDD(1)-FACT*(DWDD(2)-DWDD(1))
!
!  CONVERT WIND DIRECTION BACK TO CORRECT ORIENTATION (180 deg ROTATION)
!
!qa       WINDIR=WINDIR+180.
!qa       WINDIR=AMOD(WINDIR,360.)
!
!  ADDED RAMPING TO WIND STRESSES AND HEAT FLUX
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K) SCHEDULE(STATIC)
       DO 4120 J=1,JM
         DO 4120 I=1,IM
!NG12162010 INCLUDE RAIN/EVAP VELOCITY
!NG12162010           WSSURF(I,J)=0.0
           IF(FSM(I,J).EQ.0) GO TO 4120 ! NG12162010
!           IF ((DQPREC(1).EQ.0.OR.INTX.EQ.ISTART).AND.DQPREC(2).GT.0.)
!     +      THEN
           PRECTEMP(I,J)=AMAX1(T(I,J,1),0.)    ! NG12162010 ASSUME RAIN IS AT THERMAL EQUILIBRIUM WITH SURFACE AT BEGINNING OF RUN OR EVENT (AND KEEP IT)
!            WRITE(IUPRT,4121) I,J,PRECTEMP(I,J) ! NG12162010 IT WOULD BE BETTER TO PROVIDE AIRT IN RUN_DATA AS WELL FOR THIS OPTION
! 4121 FORMAT(' IJ:',2I5,' STARTED RAINING WITH TEMP SET TO ',F7.3,'C')
!           ENDIF
           IF(DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).GT.QEVAP2D(I,J)) THEN ! IF RAINING ON DRY CELL (NG04132011)
            DO K=1,KBM1
             S(I,J,K)=0.0                                    ! SET ITS NOMINAL SALINITY TO ZERO
             SB(I,J,K)=0.0                                   ! SET ITS NOMINAL SALINITY TO ZERO
             T(I,J,K)=PRECTEMP(I,J)                          ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
             TB(I,J,K)=PRECTEMP(I,J)                         ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
            ENDDO
           ENDIF
           WSSURF(I,J)=RAMP*FSM(I,J)*(S(I,J,1)*QEVAP2D(I,J)- 
     .                                0.*QPREC2D(I,J))       ! NG12162010 UPWARD SURFACE SALT FLUX (RAIN HAS NO SALT)
           WTSURF(I,J)=RAMP*HFLUX/(-4.186E6)
     .                +RAMP*FSM(I,J)*(T(I,J,1)*QEVAP2D(I,J)-            ! NG12162010 UPWARD SURFACE TEMP FLUX
     .                                PRECTEMP(I,J)*QPREC2D(I,J))       ! NG12162010 (RAIN ASSUMED AT THERMAL EQUILIBRIUM WITH SURFACE, SINCE NO AIRT GIVEN)
           TXSURF(I,J)=RAMP*(-1.E-3)*
     +               (+TX*COS(ANG(I,J))+TY*SIN(ANG(I,J)))*FSM(I,J)
           TYSURF(I,J)=RAMP*
     +      (-1.E-3)*(-TX*SIN(ANG(I,J))+TY*COS(ANG(I,J)))*FSM(I,J)
 4120  CONTINUE
!
      ELSE if(OPTMBC(1:8).eq.'SYNOPTIC') then   ! JI
!-------- SYNOPTIC (PRECIPITATION, EVAPORATION & HEAT FLUX) ------------
       IF(THOUR.LT.T2M) GOTO 4210
       T1M=T2M
*Khan 072898
!C       DQPREC(1)=DQPREC(2) 
!C       DQEVAP(1)=DQEVAP(2) 
!C       DHFLUX(1)=DHFLUX(2)
!C       READ (IUT93,4000,END=4016) T2M
!C       READ (IUT93,4000) DQPREC(2),DQEVAP(2),DHFLUX(2)
       DO 4220 J=1,JM
       DO 4220 I=1,IM
 4220  DHFLX2D(I,J,1)=DHFLX2D(I,J,2)
       READ (IUT192,END=4016) T2M
       READ (IUT192) ((DHFLX2D(I,J,2),I=1,IM),J=1,JM)
 4210  CONTINUE
       FACT=(THOUR-T1M)/(T2M-T1M)
!C       QPREC=DQPREC(1)+FACT*(DQPREC(2)-DQPREC(1)) 
!C       QEVAP=DQEVAP(1)+FACT*(DQEVAP(2)-DQEVAP(1)) 
!C       HFLUX=DHFLUX(1)+FACT*(DHFLUX(2)-DHFLUX(1))
       DO 4230 J=1,JM   ! n.kim 7.7.98
       DO 4230 I=1,IM   ! n.kim 7.7.98
 4230  TH2D(I,J)=DHFLX2D(I,J,1)+FACT*(DHFLX2D(I,J,2)-DHFLX2D(I,J,1))
!      end of modification n.kim 7.7.98
*Khan 072898
!-------- SYNOPTIC (WIND STRESS) -------------------------------------
!
       IF(THOUR.LT.T2W) GOTO 4310
       T1W=T2W
       DO 4320 J=1,JM
       DO 4320 I=1,IM
       DPATM(I,J,1)=DPATM(I,J,2) 
       DTX2D(I,J,1)=DTX2D(I,J,2) 
       DTY2D(I,J,1)=DTY2D(I,J,2) 
       DWSX2D(I,J,1)=DWSX2D(I,J,2)
 4320  DWSY2D(I,J,1)=DWSY2D(I,J,2)
       READ (IUT95,END=4016) T2W
       READ (IUT95) ((DTX2D(I,J,2),DTY2D(I,J,2),DPATM(I,J,2),
     +                DWSX2D(I,J,2),DWSY2D(I,J,2), ! NG04032008
     +                I=1,IM),J=1,JM)
 4310  CONTINUE
       FACT=(THOUR-T1W)/(T2W-T1W)
       DO 4330 J=1,JM
       DO 4330 I=1,IM
       PATM(I,J)=DPATM(I,J,1)+FACT*(DPATM(I,J,2)-DPATM(I,J,1))
       TX2D(I,J)=DTX2D(I,J,1)+FACT*(DTX2D(I,J,2)-DTX2D(I,J,1))
       TY2D(I,J)=DTY2D(I,J,1)+FACT*(DTY2D(I,J,2)-DTY2D(I,J,1))
       WSX2D(I,J)=DWSX2D(I,J,1)+FACT*(DWSX2D(I,J,2)-DWSX2D(I,J,1)) ! NG04032008
 4330  WSY2D(I,J)=DWSY2D(I,J,1)+FACT*(DWSY2D(I,J,2)-DWSY2D(I,J,1)) ! NG04032008

!NG2014
       DO N=1,NUMEBC ! CNG10032011
        IE=IETA(N)
        JE=JETA(N)
        QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
        QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
       ENDDO
       DO N=1,NUMECR ! CNG2012
        IE=ICHARD(N)
        JE=JCHARD(N)
        QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
        QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
       ENDDO
!
!  ADDED RAMPING TO WIND STRESSES AND HEAT FLUX
!
       DO 4340 J=1,JM
         DO 4340 I=1,IM
           PRECTEMP(I,J)=AMAX1(T(I,J,1),0.)    ! NG12162010 ASSUME RAIN IS AT THERMAL EQUILIBRIUM WITH SURFACE AT BEGINNING OF RUN OR EVENT (AND KEEP IT)
           IF(DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).GT.QEVAP2D(I,J)) THEN ! IF RAINING ON DRY CELL (NG04132011)
            DO K=1,KBM1
             S(I,J,K)=0.0                                    ! SET ITS NOMINAL SALINITY TO ZERO
             SB(I,J,K)=0.0                                   ! SET ITS NOMINAL SALINITY TO ZERO
             T(I,J,K)=PRECTEMP(I,J)                          ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
             TB(I,J,K)=PRECTEMP(I,J)                         ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
            ENDDO
           ENDIF
           WSSURF(I,J)=RAMP*FSM(I,J)*(S(I,J,1)*QEVAP2D(I,J)- 
     .                                0.*QPREC2D(I,J))       ! NG12162010 UPWARD SURFACE SALT FLUX (RAIN HAS NO SALT)
           WTSURF(I,J)=RAMP*TH2D(I,J)/(-4.186E6)
     .                +RAMP*FSM(I,J)*(T(I,J,1)*QEVAP2D(I,J)-            ! NG12162010 UPWARD SURFACE TEMP FLUX
     .                                PRECTEMP(I,J)*QPREC2D(I,J))       ! NG12162010 (RAIN ASSUMED AT THERMAL EQUILIBRIUM WITH SURFACE, SINCE NO AIRT GIVEN)
*Khan 072898
           TXSURF(I,J)=RAMP*(-1.E-3)*(+TX2D(I,J)*COS(ANG(I,J))+
     .                      TY2D(I,J)*SIN(ANG(I,J)))*FSM(I,J)
           TYSURF(I,J)=RAMP*(-1.E-3)*(-TX2D(I,J)*SIN(ANG(I,J))+
     .                      TY2D(I,J)*COS(ANG(I,J)))*FSM(I,J)
 4340  CONTINUE
!
!---------- heatflux calculation ------------------------------------
!  MODIFIED BY QUAMRUL QA ON 3/10/99 **********************************
!
      ELSE IF (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     *          OPTMBC(1:8).EQ.'RANDMFLX'.OR.
     *          OPTMBC(1:8).EQ.'AANDBFLX'.OR.
     *          OPTMBC(1:8).EQ.'SYNOPANB'.OR.
     *          OPTMBC(1:8).EQ.'SYNOPLNP'.OR.
     *          OPTMBC(1:8).EQ.'SYNOPRNM') Then
        IF (THOUR.GE.T2M) THEN
!
          T1M = T2M
          DQPREC(1) = DQPREC(2)
          DQEVAP(1) = DQEVAP(2)
          DAIRTM(1) = DAIRTM(2)
          DRELHU(1) = DRELHU(2)
          DBAROP(1) = DBAROP(2)
          DSWOBS(1) = DSWOBS(2)
          DTX(1) = DTX(2)
          DTY(1) = DTY(2)
          dwsx(1)=dwsx(2)
          dwsy(1)=dwsy(2)  
          DWDS(1) = DWDS(2)
          DWDD(1) = DWDD(2)
!         MODIFIED BY QUAMRUL QA ON 3/10/99 **********************************
          CLOUD(1)  = CLOUD(2)
          EXTCOEF(1)  = EXTCOEF(2)
!
!       MODIFIED BY QUAMRUL QA ON 5/6/96 **********************************
           Read (IUT93,5000,End=4016) T2M
           Read (IUT93,5000) DAIRTM(2), DRELHU(2),
     *     DBAROP(2), DSWOBS(2), DTX(2), DTY(2), DWSX(2), DWSY(2),
     *     CLOUD(2), EXTCOEF(2),DQPREC(2), DQEVAP(2)
!
        END IF
        FACT = (THOUR-T1M) / (T2M-T1M)
        QPREC = DQPREC(1) + FACT * (DQPREC(2)-DQPREC(1))
        QEVAP = DQEVAP(1) + FACT * (DQEVAP(2)-DQEVAP(1))
        AIRTMP = DAIRTM(1) + FACT * (DAIRTM(2)-DAIRTM(1))
        RELHUM = DRELHU(1) + FACT * (DRELHU(2)-DRELHU(1))
        RHBY100 = RELHUM/100.  ! ROSATI and Miyakoda (RANDMFLX)
        BAROP = DBAROP(1) + FACT * (DBAROP(2)-DBAROP(1))
        BAROPNT = BAROP*100.0  !mb to Newton/m2 to use ROSA heatflux
        SWOBS = DSWOBS(1) + FACT * (DSWOBS(2)-DSWOBS(1))
!NG09222011 NOT HERE FOR NCLD!        SWOBS = SWOBS*(1.0-REFL) ! ALLOWING SWOBS TO REFLECT FROM WATER SURFACE (REFL)
        CLDFRC = CLOUD(1) + FACT * (CLOUD(2)-CLOUD(1))
        EXTC   = EXTCOEF(1) + FACT * (EXTCOEF(2)-EXTCOEF(1))
        TX = DTX(1) + FACT * (DTX(2)-DTX(1))
        TY = DTY(1) + FACT * (DTY(2)-DTY(1))
        WSX = DWSX(1) + FACT * (DWSX(2)-DWSX(1))
        WSY = DWSY(1) + FACT * (DWSY(2)-DWSY(1))
        IF(CREADPET) THEN
         DO 4116 J=1,JM
          DO 4116 I=1,IM
           WSX2D(I,J) = WSX
           WSY2D(I,J) = WSY
           TX2D(I,J) = TX
           TY2D(I,J) = TY
           PATM(I,J)=BAROPNT ! NG2011 FOR SAVING
           CLDFRC2D(I,J)=CLDFRC ! NG2011 FOR SAVING
           AIRTM2D(I,J)=AIRTMP ! NG2011 FOR SAVING
           RELHU2D(I,J)=RELHUM ! NG2011 FOR SAVING
           SWOBS2D(I,J)=SWOBS*(1.0-REFL) ! NG2011 FOR SAVING
4116     CONTINUE
        ELSE
         DO 4117 J=1,JM
          DO 4117 I=1,IM
!NG04092008            WU(I,J) = UWIND
!NG04092008            WV(I,J) = VWIND
           WSX2D(I,J) = WSX
           WSY2D(I,J) = WSY
           TX2D(I,J) = TX
           TY2D(I,J) = TY
           PATM(I,J)=BAROPNT ! NG2011 FOR SAVING
           CLDFRC2D(I,J)=CLDFRC ! NG2011 FOR SAVING
           AIRTM2D(I,J)=AIRTMP ! NG2011 FOR SAVING
           RELHU2D(I,J)=RELHUM ! NG2011 FOR SAVING
           SWOBS2D(I,J)=SWOBS*(1.0-REFL) ! NG2011 FOR SAVING
!NG04192011 FOR PRECIPITATION AND EVAPORATION, CHECK FOR NET EVAPORATION AND TREAT AS 2D
           QPREC2D(I,J) = FSM(I,J)*QPREC
           QEVAP2D(I,J) = FSM(I,J)*QEVAP
           IF (DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).LT.QEVAP2D(I,J)) THEN! DO NOT ALLOW NET EVAP ON DRY CELL
            QEVAP2D(I,J)=0.0
            QPREC2D(I,J)=0.0
           ENDIF
4117     CONTINUE
        ENDIF
        DO N=1,NUMEBC ! CNG10032011
         IE=IETA(N)
         JE=JETA(N)
         QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
         QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
        ENDDO
        DO N=1,NUMECR ! CNG2012
         IE=ICHARD(N)
         JE=JCHARD(N)
         QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
         QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
        ENDDO
!      End Quamrul QA 10/17/98

!-------- VARYING EXTINCTION COEFFICIENT -----------------------------
!
        IF(OPTEXTC.EQ.'VARI') THEN
          IF(THOUR.LT.T2EX) GOTO 7310
          T1EX=T2EX
          DO J=1,JM
            DO I=1,IM
              DEXTC(I,J,1)=DEXTC(I,J,2) 
            ENDDO
          ENDDO
          READ (IUT194,END=4016) T2EX
          READ (IUT194)((DEXTC(I,J,2),I=1,IM),J=1,JM)
7310      CONTINUE
          FACT=(THOUR-T1EX)/(T2EX-T1EX)
          DO J=1,JM
            DO I=1,IM
              EXTC1(I,J)=DEXTC(I,J,1)+FACT*(DEXTC(I,J,2)-DEXTC(I,J,1))
            ENDDO
          ENDDO 
        ENDIF
!
!
!pl        WINDSP=RAMP*(DWDS(1)-FACT*(DWDS(2)-DWDS(1)))
!qa        WINDSP=RAMP*(DWDS(1)+FACT*(DWDS(2)-DWDS(1)))
!qa        WINDIR=DWDD(1)-FACT*(DWDD(2)-DWDD(1))
!
!  CONVERT WIND DIRECTION BACK TO CORRECT ORIENTATION (180 deg ROTATION)
!
!qa        WINDIR=WINDIR+180.
!qa        WINDIR=AMOD(WINDIR,360.)
!
! find max and min sea surface temperatures
        SSTMIN = 100.
        SSTMAX = -100.
!NG10252011 Also calculate a more representative salinity for the VAPOR formula
        RMINSAL = 100.
        RMAXSAL = -100.
        Do 335 J = 1, JM
          Do 334 I = 1, IM
            If (FSM(I,J).EQ.1) Then
              SSTMAX = AMAX1(SSTMAX,T(I,J,1))
              SSTMIN = AMIN1(SSTMIN,T(I,J,1))
              RMAXSAL = AMAX1(RMAXSAL,S(I,J,1)) ! NG20111025
              RMINSAL = AMIN1(RMINSAL,S(I,J,1)) ! NG20111025
            End If
  334     Continue
  335   Continue
! expand the range slightly to avoid divide by zero problems at T=SSTMIN
! and T=SSTMAX
        SSTMIN = SSTMIN - .001
        SSTMAX = SSTMAX + .001
! Set a few constants
        UZ = SQRT(WSX*WSX+WSY*WSY)
        UZ = UZ* WNDSH                        ! Wind Sheltering Coeff WNDSH
        UZ = AMAX1(UZ,0.1) !UZ should be a non-zero number (hli,06/10/02).
        ZU = 10.
        ZT = 10.
        ZQ = 10.
        SPECHT = 1004.  ! Specific Heat (Joules/kg/deg C)
        RLHE = 2.47E3   ! Latent heat of evaporation (Joules/g)
        RHOAIR = (1.25E-03)     ! air density
!       Get Wind Direction for ROSATI and Miyakoda Heatflux
        Call NORTH(WSX,WSY,ANGLE)
!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 ****************************************
!  If CLDFRC is negative from Data i.e. missing then we wil calculate CLDFRC
!  from ncld.f, otherwise it will use cloud data ***************************
!
        IF (CLDFRC.LT.0.0) THEN
! days since start of run
         TDAY = THOUR / 24.
! Add days since start of run to Julian date of start of run
         CUJDAY = TDAY + SDAY

! subtract out Julian date of Jan 1 from current year after updating
! current year (if necessary)
         If (CUJDAY.LT.FLOAT(J2JDAY)) Then
           YDAY = CUJDAY - FLOAT(J1JDAY)
         Else
           YDAY = CUJDAY - FLOAT(J2JDAY)
           J1JDAY = J2JDAY
           J1YR = J1YR + 1
           J2YR = J1YR + 1
           Call JDAY(1,1,J2YR,J2JDAY)
         End If
! calculate cloud fraction for current day of the year, based on observed
! vs theoretical sw radiation changed CLDFRC TO CLDX IN ARGS 8/8/2001
!NG         Call NCLD(YDAY,ALAT,ALON,SWOBS,ATC,CLDX,CUTOFF,RSS)
         AALON=-ALON
         Call NCLD(YDAY,ALAT,AALON,SWOBS,ATC,CLDX,CUTOFF,RSS)
         CLDFRC = CLDX
        ENDIF
!  End MODIFIED BY QUAMRUL QA ON 3/10/99 ***********************************
!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 ********************************
!  Here We calculate ShortWave Radiation in Advance before go to the loop
!  Under 2 Circumstances 1. If OPTMBC is RANDMFLX or 2. Short Wave Data
!  is not Available i.e. Negative
!
        IF (OPTMBC(1:8).EQ.'RANDMFLX'.OR.OPTMBC(1:8).EQ.'SYNOPRNM'
     &       .OR.DSWOBS(1).LT.0.0.OR.DSWOBS(2).LT.0.) Then !if SWOBS = -999.0
          Call DATEHR(IDA,IMO,IYR,IHR,0,SHOUR)
          CHRS = SHOUR + THOUR
          Call GETDATE(IDC,IMC,IYC,IHC,IMINC,CHRS)
!NG          Call HITFLX (ALAT,ALON,IYC,IMC,IDC,IHC,BAROPNT,ZERO,
          Call HITFLX (ALAT,ALON,IYC,IMC,IDC,CHRS,BAROPNT,ZERO,
     +       AIRTMP,UZ,ANGLE,RHBY100,CLDFRC,DHR,SWR,RLW,HNOT,
     +       WENOT,UTCSHIFT)
          SWOBS = SWR
        ENDIF
!
!     ALLOWING SWOBS TO REFLECT BACK FROM WATER SURFACE (REFL)
!
!NG09222011 MOVED OUT OF PREVIOUS IF STATEMENT AND GLOBALLY UNDER NCLD SO THAT NCLD HAS FULL SWOBS
        SWOBS = SWOBS*(1.0-REFL)

! number of interpolation steps
        NINTPT = MIN0(20,INTMAX-1)   ! this value (20) can be changed
        STEP = (SSTMAX-SSTMIN) / FLOAT(NINTPT-1)
        Do 336 I = 1, NINTPT
          SSTINT(I) = SSTMIN + FLOAT(I-1) * STEP
          SSIN = SSTINT(I)
          Call BULK(UZ,ZU,AIRTMP,ZT,RELHUM,ZQ,SSIN,TOL,UW,WTOUT,WQOUT,
     *        CD,BAROP,NERR)

!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 ***************************************
!     ********************************************************************
!     NOTE: IN VAPOR ROUTINE  EA=VAPOR PRESSURE AT SSIN (SEA SURF TEMP),
!     ES=SATURATED VAPORT PRESSURE AT SSIN, IF YOU NEED TO GET SATURATED VAPOR
!     PRESSURE OR VAPOR PRESSURE AT AIR TEMPERATURE YOU SHOULD USE AIR TEMP
!     INSTEAD. ALSO NOTE EA = REL HUM*ES     QUAMRUL QA  2/22/96
!     For WHOI(LANDPFLX) ROUTINE ALL VAPOR PRESSURE ARE in mbar
!     **********************************************************************
!
!NG10252011          Call VAPOR(EA,ESS,RELHUM,SSIN,BAROP)   !For Sea Surface Temperature
          SAL=0.5*(RMINSAL+RMAXSAL) ! CNG10252011 CORRECT FOR A MORE REPRESENTATIVE SEA SURFACE SALINITY !
          Call VAPOR(EA,ESS,RELHUM,SSIN,BAROP,SAL)   !For Sea Surface Temperature
!
!  Add three new options (SYNOPANB, SYNOPLNP and SYNOPRNM)
!  (hli,04/24/02)
!
          If (OPTMBC(1:8).EQ.'AANDBFLX'.OR.
     .        OPTMBC(1:8).EQ.'SYNOPANB') Then
!
!******* Atmospheric Longwave radiation ADOPTED FROM ARMY(AANDBFLX) CODE
           T2K = 273.2+AIRTMP
           FAC = 1.0+0.17*CLDFRC**2   !CLDFRC is in 0-1 scale QA
           RAN = 1000.0/3600.0*9.37E-6*SBC*T2K**6*FAC*(1.0-0.03)
!******* Reflective Longwave radiation
           RB = 1000.0/3600.0*0.97*SBC*(SSIN+273.2)**4
!
           Call VAPOR(EAS,ES,RELHUM,AIRTMP,BAROP,0.)!  For Air Temperature
           ESEADIFF = (ESS-EAS)/1.3333           ! from mb to mmHg
!          UZ = UZ* WNDSH                        ! Wind Sheltering Coeff WNDSH
!NG11022008           fw       = 9.2 + 0.46*UZ**2           ! Wind function Edinger et al. 1974 / but based on wind at 7m! NG11022008 
           fw       = 9.2 + 0.41*UZ**2           ! Wind function Edinger et al. 1974, with wind at 10m UZ7=UZ10*(7/10)**0.15, TVA 1972, NG11022008
           WT(I)    = -0.47*fw*(SSIN-AIRTMP)     ! sensible heat flux (w/m2)
           WQ(I)    = -fw*ESEADIFF               ! Latent heat flux (w/m2)
           RLWC(I)  = RAN - RB
          Else If (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     .             OPTMBC(1:8).EQ.'SYNOPLNP') Then
           WT(I)    = -(RHOAIR*1.E3) * SPECHT * WTOUT  !sensible heat flux (w/m2)
           WQ(I)    = -RLHE * WQOUT              !latent heat flux (w/m2)
           Call LONGWAVE(SSIN,AIRTMP,EA,CLDFRC,RL1,RL2,RLW,RLWCO)
           RLWC(I) = -RLWCO
          Else If (OPTMBC(1:8).EQ.'RANDMFLX'.OR.
     .             OPTMBC(1:8).EQ.'SYNOPRNM') Then
!NG            Call HITFLX (ALAT,ALON,IYC,IMC,IDC,IHC,BAROPNT,SSIN,
            Call DATEHR(IDA,IMO,IYR,IHR,0,SHOUR)
            CHRS = SHOUR + THOUR ! CNG
            Call HITFLX (ALAT,ALON,IYC,IMC,IDC,CHRS,BAROPNT,SSIN,
     +           AIRTMP,UZ,ANGLE,RHBY100,CLDFRC,DHR,SWR,RLW,HNOT,
     +           WENOT,UTCSHIFT)
            WT(I)   = -HNOT     ! Sensible Heat
            WQ(I)   = -WENOT    ! Latent Heat
            RLWC(I) = -RLW      ! Longwave
          Endif

!
!  END MODIFIED BY QUAMRUL QA ON 3/10/99 *************************************
!
  336   Continue
!
!       Set one extra value equal to last NINTPT value to make available for
!       Interpolation purposes in next few statements.... by QA 3/25/98
!
        WT(NINTPT+1) = WT(NINTPT)
        WQ(NINTPT+1) = WQ(NINTPT)
        RLWC(NINTPT+1) = RLWC(NINTPT)
!qa
!jah    Initialization (CHECK!!!!)
!       ESSAVE = 0.0
!       EASAVE = 0.0
!       HFSAVE = 0.0
!       WT1SAVE = 0.0
!       WT2SAVE = 0.0
!       WT3SAVE = 0.0
!       SWRSAVE = 0.0


        CALL WETBULB(RELHUM,AIRTMP,BAROP,DEW,TWET) ! NG10262011 CALCULATE WET BULB TEMPERATURE FOR RAIN, MOVED HERE 2014
        Do 350 J = 1, JM
          Do 340 I = 1, IM
!**********************************************************************
!   surface heat and salt fluxes
!**********************************************************************
!NG12162010 NOT NEEDED HERE, DONE FURTHER DOWN , AND APPEARS AS A WRONG FORMULA ANYWAY.
!NG12162010             WSSURF(I,J) = 33. / 1000 * (QPREC-QEVAP) / (ART(I,J)+1.E-16)

!           interpolate heatfluxes from those calculated
!           First, find the interpolation factor

            If (FSM(I,J).EQ.1) THEN
!NG             If (T(I,J,1).eq.SSTMAX-0.001)write(*,*) I,J,ETF(I,J),H(I,J)
              ILOW = (T(I,J,1)-SSTMIN) / STEP + 1
              IHIGH = ILOW + 1
              FACTWT = (T(I,J,1)-SSTINT(ILOW)) / (SSTINT(IHIGH)-
     *            SSTINT(ILOW))

! Interpolate heat fluxes from those calculated

              WT1 = WT(ILOW) + FACTWT * (WT(IHIGH)-WT(ILOW))
              WT2 = WQ(ILOW) + FACTWT * (WQ(IHIGH)-WQ(ILOW))
              WT3 = RLWC(ILOW) + FACTWT * (RLWC(IHIGH)-RLWC(ILOW))

              HFLUX = WT1 + WT2 + WT3 + SWOBS
              HFLUXmin = WT(1) + WQ(1) + RLWC(1) + SWOBS
              HFLUXmax = WT(NINTPT) + WQ(NINTPT) + RLWC(NINTPT) + SWOBS
              SWOBS2D(I,J) = SWOBS ! SAVE TOO 12022008
!NG12162010 INCLUDE RAIN/EVAP VELOCITY
!NG12162010 THE ONE BELOW FROM ECOMSED APPEARS WRONG TO ME: QPREC,QEVAP ARE VELOCITIES, AND 33ppt APPEARS HARDWIRED, AND NO RAMP HERE
!NG12162010              WSSURF(I,J) = RAMP* 33. / 1000 * (QPREC-QEVAP) /(ART(I,J)+
!NG12162010     *            1.E-16)
              PRECTEMP(I,J)=AMAX1(TWET,0.)                      ! NG10262011 USE WET BULB TEMPERATURE FOR RAIN
              IF(DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).GT.QEVAP2D(I,J))THEN ! IF RAINING ON DRY CELL (NG04132011)
               DO K=1,KBM1
                S(I,J,K)=0.0                             ! SET ITS NOMINAL SALINITY TO ZERO
                SB(I,J,K)=0.0                            ! SET ITS NOMINAL SALINITY TO ZERO
                T(I,J,K)=PRECTEMP(I,J)                   ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
                TB(I,J,K)=PRECTEMP(I,J)                  ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
               ENDDO
              ENDIF
              WSSURF(I,J)=RAMP*FSM(I,J)*(S(I,J,1)*QEVAP2D(I,J)-
     .                                   0.*QPREC2D(I,J)) ! NG12162010 UPWARD SURFACE SALT FLUX (RAIN HAS NO SALT)
              WTSURF(I,J) = RAMP*FSM(I,J)*HFLUX / (-4.186E6)
     .                    + RAMP*FSM(I,J)*(T(I,J,1)*QEVAP2D(I,J)-        ! NG12162010 UPWARD SURFACE TEMP FLUX
     .                                   PRECTEMP(I,J)*QPREC2D(I,J))   ! NG12162010 (RAIN ASSUMED AT THERMAL EQUILIBRIUM WITH AIR)
              SWRAD(I,J)  = RAMP*FSM(I,J)*SWOBS / (-4.186E6)

!NG UNCOMMENT FOLLOWING THREE LINES FOR HEAT FLUX, IF NEEDED FOR RCA
!      IF(I.EQ.88.AND.J.EQ.538)
!     .  Write(77,777)TIME,WT1,WT2,WT3,SWOBS,HFLUX,
!     .      ESS,EAS,SSTMIN,SSTMAX,RMAXSAL,RMINSAL
!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 *****************************************
!  This is checking purposes writing heatfluxes at first location where skill
!  assessment output is written in gcm_tsr file
              IF(I.EQ.INXIE(1).AND.J.EQ.INXJE(1))THEN
                WT1SAVE = WT1SAVE + WT1*SKILLI
                WT2SAVE = WT2SAVE + WT2*SKILLI
                WT3SAVE = WT3SAVE + WT3*SKILLI
                SWRSAVE = SWRSAVE + SWOBS*SKILLI
                HFSAVE = HFSAVE + HFLUX*SKILLI
                ESSAVE = ESSAVE + ESS*SKILLI
                EASAVE = EASAVE + EAS*SKILLI
                ATSAVE = ATSAVE + AIRTMP*SKILLI
                STSAVE = STSAVE + SSIN*SKILLI
!
!      If (ISKILL.NE.0.AND.MOD(INTX,ISKILL).EQ.0) Then
!NG                If (AMOD(TIME*24,1.).EQ.0) Then ! Every hour
                IF (ISKILL.NE.0.AND.MOD(INTX,ISKILL).EQ.0) THEN
!NG            TMIDDLE = TIME - (.5*DTI*DAYI/(SKILLI*2))
!NG                 TMIDDLE = TIME 
                 TMIDDLE = TIME - (.5*DTI*DAYI/SKILLI)
!       TMIDDLE = TIME - (.5*DTI*DAYI/(SKILLI*2))
                 Write(77,777)TMIDDLE,WT1SAVE,WT2SAVE,WT3SAVE,SWRSAVE,
     .                        HFSAVE,ESSAVE,EASAVE,SSTMIN,SSTMAX,
     .                        RMAXSAL,RMINSAL
                 ESSAVE = 0.0
                 EASAVE = 0.0
                 HFSAVE = 0.0
                 WT1SAVE = 0.0
                 WT2SAVE = 0.0
                 WT3SAVE = 0.0
                 SWRSAVE = 0.0
                Endif
              ENDIF
777           FORMAT(12F12.5)
            End If
  340     Continue
  350   Continue
!-----------------------------------------------------------------
!-------- SYNOPTIC (WIND STRESSES) -------------------------------------
!NG02232007        IF(OPTMBC.EQ.'SYNOPANB'.OR.
!NG02232007     &     OPTMBC.EQ.'SYNOPLNP'.OR.
!NG02232007     &     OPTMBC.EQ.'SYNOPRNM') THEN
        IF(OPTMBC(1:8).EQ.'SYNOPANB'.OR.
     &     OPTMBC(1:8).EQ.'SYNOPLNP'.OR.
     &     OPTMBC(1:8).EQ.'SYNOPRNM') THEN
          IF(THOUR.LT.T2W) GOTO 5221
          T1W=T2W
          DO J=1,JM
            DO I=1,IM
               DPATM(I,J,1)=DPATM(I,J,2)
               DTX2D(I,J,1)=DTX2D(I,J,2)
               DTY2D(I,J,1)=DTY2D(I,J,2)
               DWSX2D(I,J,1)=DWSX2D(I,J,2)
               DWSY2D(I,J,1)=DWSY2D(I,J,2)
            ENDDO
          ENDDO
          READ (IUT95,END=4016) T2W
          READ (IUT95) ((DTX2D(I,J,2),DTY2D(I,J,2),DPATM(I,J,2),
     +                   DWSX2D(I,J,2),DWSY2D(I,J,2), ! NG04032008
     &                   I=1,IM),J=1,JM)
 5221     CONTINUE
          FACT=(THOUR-T1W)/(T2W-T1W)
          DO J=1,JM
            DO I=1,IM
             PATM(I,J)=DPATM(I,J,1)+FACT*(DPATM(I,J,2)-DPATM(I,J,1))
             TX2D(I,J)=DTX2D(I,J,1)+FACT*(DTX2D(I,J,2)-DTX2D(I,J,1))
             TY2D(I,J)=DTY2D(I,J,1)+FACT*(DTY2D(I,J,2)-DTY2D(I,J,1))
             WSX2D(I,J)=DWSX2D(I,J,1)+FACT*(DWSX2D(I,J,2)-DWSX2D(I,J,1))! NG04032008
             WSY2D(I,J)=DWSY2D(I,J,1)+FACT*(DWSY2D(I,J,2)-DWSY2D(I,J,1))! NG04032008
            ENDDO
          ENDDO
        END IF
!
! 
!  MODIFIED BY C.K.Z. ON 11/4/94
!
!  ADDED RAMPING TO WIND STRESSES AND HEAT FLUX
!
        DO 5223 J=1,JM
          DO 5223 I=1,IM
          TXSURF(I,J)=RAMP*(-1.E-3)*(+TX2D(I,J)*COS(ANG(I,J))+
     .                       TY2D(I,J)*SIN(ANG(I,J)))*FSM(I,J)
          TYSURF(I,J)=RAMP*(-1.E-3)*(-TX2D(I,J)*SIN(ANG(I,J))+
     .                       TY2D(I,J)*COS(ANG(I,J)))*FSM(I,J)
!NG04192011; IT SHOULD ONLY BE RAMPED AT FORCING TIME IN EXTERNAL AND ADV?           PATM(I,J) = PATM(I,J)*RAMP
 5223   CONTINUE
!
!  MODIFIED BY NICHOLAS KIM ON 8/11/04 **********************************
!  2D WIND AND HEATFLUX COMPUTATIONS
!
      ELSE IF (OPTMBC(1:8).EQ.'2DMETANB'.OR.
     *         OPTMBC(1:8).EQ.'2DMETLNP'.OR.
     *         OPTMBC(1:8).EQ.'2DMETRNM') Then
        IF (THOUR.GE.T2M) THEN
!
          T1M = T2M
          DO 4410 J=1,JM
           DO 4410 I=1,IM
           IF(FSM(I,J).EQ.0)GO TO 4410
           DAIRTM2D(I,J,1) = DAIRTM2D(I,J,2)
           DRELHU2D(I,J,1) = DRELHU2D(I,J,2)
           DPATM(I,J,1) = DPATM(I,J,2) ! NG04092008 DBAROP2D(I,J,1) = DBAROP2D(I,J,2)
           DSWOBS2D(I,J,1) = DSWOBS2D(I,J,2)
           DWSX2D(I,J,1)=DWSX2D(I,J,2)
           DWSY2D(I,J,1)=DWSY2D(I,J,2)  
           DCLD2D(I,J,1)  = DCLD2D(I,J,2)
           DEXTC2D(I,J,1)  = DEXTC2D(I,J,2)
 4410    CONTINUE
!
!   UPDATE DATA FROM 'synop_met' data
          ! A.Jordi: Read wind and heat flux parameters from netcdf file
          NT_T2M=NT_T2M+1
          CALL NETCDF_READ_SYNOP_MET(NT_T2M)
          ! A.Jordi: end
           ! NG12162010 NEED TO ADD DQPREC2D(I,J,2),DQEVAP2D(I,J,2)
        END IF
        FACT = (THOUR-T1M) / (T2M-T1M)
        DO 4400 J=1,JM
         DO 4400 I=1,IM
          IF(FSM(I,J).EQ.0.0)GO TO 4400
          AIRTM2D(I,J) = DAIRTM2D(I,J,1)
     .           + FACT*(DAIRTM2D(I,J,2)-DAIRTM2D(I,J,1)) ! AIR TEMP
          RELHU2D(I,J) = DRELHU2D(I,J,1)
     .           + FACT*(DRELHU2D(I,J,2)-DRELHU2D(I,J,1))  ! REL HUMIDITY
          PATM(I,J) = DPATM(I,J,1)
     .           + FACT*(DPATM(I,J,2)-DPATM(I,J,1))  ! BARO PRESSURE
          PATM(I,J) = PATM(I,J)*1.0E02 ! from mbar to N/m2 (same as in bcdata for synop_wind) NG04092008
!NG04092008          BAROP2D(I,J) = DBAROP2D(I,J,1)
!NG04092008     .           + FACT*(DBAROP2D(I,J,2)-DBAROP2D(I,J,1))  ! BARO PRESSURE
          SWOBS2D(I,J) = DSWOBS2D(I,J,1) 
     .           + FACT*(DSWOBS2D(I,J,2)-DSWOBS2D(I,J,1))  ! SHORTWAVE RAD
          CLDFRC2D(I,J) = DCLD2D(I,J,1)
     .           + FACT*(DCLD2D(I,J,2)-DCLD2D(I,J,1))    ! CLOUD COVER
          EXTC1(I,J)  = DEXTC2D(I,J,1)
     .           + FACT*(DEXTC2D(I,J,2)-DEXTC2D(I,J,1))
          WSX2D(I,J) = DWSX2D(I,J,1) +FACT*(DWSX2D(I,J,2)-DWSX2D(I,J,1)) !m/sec 
          WSY2D(I,J) = DWSY2D(I,J,1) +FACT*(DWSY2D(I,J,2)-DWSY2D(I,J,1)) 
 4400   CONTINUE
        DO N=1,NUMEBC ! CNG10032011
         IE=IETA(N)
         JE=JETA(N)
         QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
         QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON BC CELLS
        ENDDO
        DO N=1,NUMECR ! CNG2012
         IE=ICHARD(N)
         JE=JCHARD(N)
         QEVAP2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
         QPREC2D(IE,JE)=0. ! DO NOT RAIN/EVAP ON CORNER BC CELL
        ENDDO
!
!
        IF(CALCCLD.EQ.0)GO TO 4401  ! SKIPPING COMPUTING CLOUD COVER
!NG09212011 REARRANGED AND 2Ded FOR 2D-VARIABLE CLOUD COVER CALCULATION
! days since start of run
         TDAY = THOUR / 24.
! Add days since start of run to Julian date of start of run
         CUJDAY = TDAY + SDAY

! subtract out Julian date of Jan 1 from current year after updating
! current year (if necessary)
         IF (CUJDAY.LT.FLOAT(J2JDAY)) THEN
          YDAY = CUJDAY - FLOAT(J1JDAY)
         ELSE
          YDAY = CUJDAY - FLOAT(J2JDAY)
          J1JDAY = J2JDAY
          J1YR = J1YR + 1
          J2YR = J1YR + 1
          CALL JDAY(1,1,J2YR,J2JDAY)
         ENDIF
         DO 4404 J=1,JM
          DO 4404 I=1,IM
!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 ****************************************
!  If CLDFRC is negative from Data i.e. missing then we wil calculate CLDFRC
!  from ncld.f, otherwise it will use cloud data ***************************
!
           IF(FSM(I,J).EQ.0) GOTO 4404
           IF (CLDFRC2D(I,J).LT.0.0) THEN
! calculate cloud fraction for current day of the year, based on observed
! vs theoretical sw radiation changed CLDFRC TO CLDX IN ARGS 8/8/2001
!NG           Call NCLD(YDAY,ALAT,ALON,SWOBS,ATC,CLDX,CUTOFF,RSS)
            SWOBS=SWOBS2D(I,J)
            AALAT=-YGRID(I,J)
            AALON=XGRID(I,J)
            Call NCLD(YDAY,AALON,AALAT,SWOBS,ATC,CLDX,CUTOFF,RSS)
            CLDFRC2D(I,J) = CLDX
           ENDIF
!NG09222011 PUT THIS HERE, AFTER NCLD IS CALLED WITH FULL SWOBS
           SWOBS2D(I,J) = SWOBS2D(I,J)*(1.0-REFL) ! ALLOWING SWOBS TO REFLECT
                                                  ! FROM WATER SURFACE (REFL)
 4404   CONTINUE

 4401   CONTINUE
!  End MODIFIED BY QUAMRUL QA ON 3/10/99 ***********************************
!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 ********************************
!  Here We calculate ShortWave Radiation in Advance before go to the loop
!  Under 2 Circumstances 1. If OPTMBC is RANDMFLX or 2. Short Wave Data
!  is not Available i.e. Negative
!
!         If (OPTMBC(1:8).EQ.'RANDMFLX'.OR.OPTMBC(1:8).EQ.'SYNOPRNM'
!    &         .OR.DSWOBS(1).LT.0.0.OR.DSWOBS(2).LT.0.) Then
        IF(OPTMBC(1:8).EQ.'2DMETRNM'.OR.CALCSWOBS.NE.0.0)THEN !if SWOBS = -999.0
            CALL DATEHR(IDA,IMO,IYR,IHR,0,SHOUR)
            CHRS = SHOUR + THOUR
            CALL GETDATE(IDC,IMC,IYC,IHC,IMINC,CHRS)
          DO 4402 J=1,JM
           DO 4402 I=1,IM
            IF(FSM(I,J).EQ.0.0)GO TO 4402
            ALAT=YGRID(I,J)
            ALON=XGRID(I,J)
            CLDFRC=CLDFRC2D(I,J)
!NG04092008            BAROPNT=BAROP2D(I,J)*100.   ! mb to N/m2 for use in RNM heatflux
            BAROPNT=PATM(I,J)   ! Already in N/m2 for use in RNM heatflux
            AIRTMP=AIRTM2D(I,J)
            RHBY100=RELHU2D(I,J)/100.  ! FOR ROSATI AND MIYAKOTA
            WSX=WSX2D(I,J)
            WSY=WSY2D(I,J)
!           CALL NORTH(WSX,WSY,ANGLE) ! ANGLE FOR hitflx.f NO NEED
!NG            Call HITFLX (ALAT,ALON,IYC,IMC,IDC,IHC,BAROPNT,ZERO,
            Call HITFLX (ALAT,ALON,IYC,IMC,IDC,CHRS,BAROPNT,ZERO,
     +         AIRTMP,UZ,ANGLE,RHBY100,CLDFRC,DHR,SWR,RLW,HNOT,
     +         WENOT,UTCSHIFT)
!
!     ALLOWING SWOBS TO REFLECT BACK FROM WATER SURFACE (REFL)
!
            SWOBS2D(I,J) = SWR*(1.0-REFL)


!           IF(I.EQ.10.AND.J.EQ.10)
!    .    WRITE(1004,3000)TIME,I,J,SWOBS2D(I,J)
 4402     CONTINUE
        ENDIF

!   MAIN LOOP TO COMPUTE HEAT FLUX COMPONENTS
        DO 4403 J=1,JM
         DO 4403 I=1,IM
          IF(FSM(I,J).EQ.0)GO TO 4403
 
          SSIN = T(I,J,1)
          AIRTMP=AIRTM2D(I,J)
          CLDFRC=CLDFRC2D(I,J) ! CNG09212011 CORRECTION, PREVIOUSLY NOT HERE!
          BAROPNT=PATM(I,J)    ! CNG09212011 CORRECTION, PREVIOUSLY NOT HERE!
          UZ=SQRT(WSX2D(I,J)*WSX2D(I,J)+WSY2D(I,J)*WSY2D(I,J))
          UZ = UZ* WNDSH                        ! Wind Sheltering Coeff WNDSH
          UZ = AMAX1(UZ,0.1) !UZ should be a non-zero number
          RELHUM=RELHU2D(I,J)
          BAROP=PATM(I,J)/100. ! Back-convert to mbar for bulk NG04092008 BAROP2D(I,J)

          CALL BULK(UZ,ZU,AIRTMP,ZT,RELHUM,ZQ,SSIN,TOL,UW,WTOUT,WQOUT,
     *        CD,BAROP,NERR)

!
!     ********************************************************************
!     NOTE: IN VAPOR ROUTINE  EA=VAPOR PRESSURE AT SSIN (SEA SURF TEMP),
!     ES=SATURATED VAPORT PRESSURE AT SSIN, IF YOU NEED TO GET SATURATED VAPOR
!     PRESSURE OR VAPOR PRESSURE AT AIR TEMPERATURE YOU SHOULD USE AIR TEMP
!     INSTEAD. ALSO NOTE EA = REL HUM*ES     QUAMRUL QA  2/22/96
!     For WHOI(LANDPFLX) ROUTINE ALL VAPOR PRESSURE ARE in mbar
!     **********************************************************************
!
!NG10252011          Call VAPOR(EA,ESS,RELHUM,SSIN,BAROP)   !For Sea Surface Temperature
          SAL=S(I,J,1) ! CNG10252011 CORRECT FOR ACTUAL SEA SURFACE SALINITY !
          Call VAPOR(EA,ESS,RELHUM,SSIN,BAROP,SAL)   !For Sea Surface Temperature
!         WRITE(1004,*)'BULK',I,J,RELHUM,SSIN,BAROP,EA,ESS
!
          IF (OPTMBC(1:8).EQ.'2DMETANB') THEN
!
!******* Atmospheric Longwave radiation ADOPTED FROM ARMY(AANDBFLX) CODE
           T2K = 273.2+AIRTMP
           FAC = 1.0+0.17*CLDFRC**2   !CLDFRC is in 0-1 scale QA
           RAN = 1000.0/3600.0*9.37E-6*SBC*T2K**6*FAC*(1.0-0.03)
!******* Reflective Longwave radiation
           RB = 1000.0/3600.0*0.97*SBC*(SSIN+273.2)**4
!
!NG10252011           Call VAPOR(EAS,ES,RELHUM,AIRTMP,BAROP)!  For Air Temperature
           Call VAPOR(EAS,ES,RELHUM,AIRTMP,BAROP,0.)!  For Air Temperature
           ESEADIFF = (ESS-EAS)/1.3333           ! from mb to mmHg
!          UZ = UZ* WNDSH                        ! Wind Sheltering Coeff WNDSH
!NG11022008           fw       = 9.2 + 0.46*UZ**2           ! Wind function Edinger et al. 1974 / but based on wind at 7m! NG11022008 
           fw       = 9.2 + 0.41*UZ**2           ! Wind function Edinger et al. 1974, with wind at 10m UZ7=UZ10*(7/10)**0.15, TVA 1972, NG11022008
           WT2D    = -0.47*fw*(SSIN-AIRTMP)     ! sensible heat flux (w/m2)
           WQ2D    = -fw*ESEADIFF               ! Latent heat flux (w/m2)
           RLWC2D  = RAN - RB

          ELSE IF (OPTMBC(1:8).EQ.'2DMETLNP') THEN
           
           WT2D    = -(RHOAIR*1.E3)*SPECHT*WTOUT  !sensible heat flux (w/m2)
           WQ2D    = -RLHE * WQOUT              !latent heat flux (w/m2)
           CALL LONGWAVE(SSIN,AIRTMP,EA,CLDFRC,RL1,RL2,RLW,RLWCO)
           RLWC2D = -RLWCO
!          WRITE(1004,*)'WTOUT RLWC2D',WTOUT,SSIN,AIRTMP,EA,CLDFRC,
!    .                RL1,RL2,RLW,RLWCO

          ELSE IF (OPTMBC(1:8).EQ.'2DMETRNM') THEN
  
            RHBY100=RELHU2D(I,J)/100.   
            CALL DATEHR(IDA,IMO,IYR,IHR,0,SHOUR)
            CHRS = SHOUR + THOUR
!NG            CALL HITFLX (ALAT,ALON,IYC,IMC,IDC,IHC,BAROPNT,SSIN,
            CALL HITFLX (ALAT,ALON,IYC,IMC,IDC,CHRS,BAROPNT,SSIN,
     +           AIRTMP,UZ,ANGLE,RHBY100,CLDFRC,DHR,SWR,RLW,HNOT,
     +           WENOT,UTCSHIFT)
            WT2D   = -HNOT     ! Sensible Heat
            WQ2D   = -WENOT    ! Latent Heat
            RLWC2D = -RLW      ! Longwave

          END IF

!
          HFLUX = WT2D + WQ2D + RLWC2D + SWOBS2D(I,J)

!
!  MODIFIED BY QUAMRUL QA ON 3/10/99 *****************************************
!  This is checking purposes writing heatfluxes at first location where skill
!  assessment output is written in gcm_tsr file
          IF(I.EQ.INXIE(1).AND.J.EQ.INXJE(1))THEN
           WT1SAVE = WT1SAVE + WT2D*SKILLI
           WT2SAVE = WT2SAVE + WQ2D*SKILLI
           WT3SAVE = WT3SAVE + RLWC2D*SKILLI
           SWRSAVE = SWRSAVE + SWOBS2D(I,J)*SKILLI
           HFSAVE = HFSAVE + HFLUX*SKILLI
           ESSAVE = ESSAVE + ESS*SKILLI
           EASAVE = EASAVE + EAS*SKILLI
           ATSAVE = ATSAVE + AIRTMP*SKILLI
           STSAVE = STSAVE + SSIN*SKILLI
 
           IF (ISKILL.NE.0.AND.MOD(INTX,ISKILL).EQ.0) THEN
!NG            TMIDDLE = TIME - (.5*DTI*DAYI/(SKILLI*2))
            TMIDDLE = TIME - (.5*DTI*DAYI/SKILLI)
            Write(77,777)TMIDDLE,WT1SAVE,WT2SAVE,WT3SAVE,SWRSAVE,
     .                   HFSAVE,ESSAVE,EASAVE,SSTMIN,SSTMAX,
     .                   RMAXSAL,RMINSAL
            ESSAVE = 0.0
            EASAVE = 0.0
            HFSAVE = 0.0
            WT1SAVE = 0.0
            WT2SAVE = 0.0
            WT3SAVE = 0.0
            SWRSAVE = 0.0
           ENDIF
          ENDIF
!qa
!
!  ADDED RAMPING TO WIND STRESSES AND HEAT FLUX
!
!**********************************************************************
!   surface heat and salt fluxes
!**********************************************************************
!NG12162010 INCLUDE RAIN/EVAP VELOCITY
!NG12162010 THE ONE BELOW FROM ECOMSED APPEARS WRONG TO ME: QPREC,QEVAP ARE VELOCITIES, AND 33ppt APPEARS HARDWIRED, AND NO RAMP HERE
!NG12162010          WSSURF(I,J) = RAMP*
!NG12162010     *     33./1000*(QPREC2D(I,J)-QEVAP2D(I,J))/(ART(I,J)+1.E-16)
          CALL WETBULB(RELHU2D(I,J),AIRTM2D(I,J),PATM(I,J)/100,DEW,TWET) ! NG10262011 CALCULATE WET BULB TEMPERATURE FOR RAIN
          PRECTEMP(I,J)=AMAX1(TWET,0.)                      ! NG10262011 USE WET BULB TEMPERATURE FOR RAIN
          IF(DT(I,J).LE.WETMIN.AND.QPREC2D(I,J).GT.QEVAP2D(I,J)) THEN ! IF RAINING ON DRY CELL (NG04132011)
           DO K=1,KBM1
            S(I,J,K)=0.0                                    ! SET ITS NOMINAL SALINITY TO ZERO
            SB(I,J,K)=0.0                                   ! SET ITS NOMINAL SALINITY TO ZERO
            T(I,J,K)=PRECTEMP(I,J)                          ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
            TB(I,J,K)=PRECTEMP(I,J)                         ! SET ITS NOMINAL TEMPERATURE TO PRECTEMP TOO
           ENDDO
          ENDIF
          WSSURF(I,J)=RAMP*FSM(I,J)*
     .                (S(I,J,1)*QEVAP2D(I,J)-0.*QPREC2D(I,J)) ! NG12162010 UPWARD SURFACE SALT FLUX (RAIN HAS NO SALT)
          WTSURF(I,J) = RAMP* HFLUX / (-4.186E6)
     .                + RAMP*FSM(I,J)*(T(I,J,1)*QEVAP2D(I,J)-        ! NG12162010 UPWARD SURFACE TEMP FLUX
     .                                 PRECTEMP(I,J)*QPREC2D(I,J))   ! NG12162010 (RAIN ASSUMED AT THERMAL EQUILIBRIUM WITH AIR)

! CONVERT WIND SPEED TO STRESS
          WDS=SQRT(WSX2D(I,J)**2+WSY2D(I,J)**2)
          RHOA=SADENS(BAROP,AIRTMP,RELHUM) ! NG06042012 
          IF(WAVEDYN.EQ.'DONONLY '.OR.WAVEDYN.EQ.'NEGLECT') THEN ! USE LARGE AND POND
           CD=1.2E-3
           IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
           IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
          ELSE ! SIMULATION INCLUDES WAVES; USE TAYLOR-YELAND
           Z0air=3.5E-5 + 4*SIG(I,J)*1200*(4*SIG(I,J)*WN(I,J)/TPI)**4.5
           CD=0.16*ALOG(10/Z0air)**(-2)
           IF(CD.GE.3E-3) CD=3E-3
          ENDIF
          TX2D(I,J)=RHOA*CD*WSX2D(I,J)*WDS
          TY2D(I,J)=RHOA*CD*WSY2D(I,J)*WDS
          IF (WDS.GE.WDSMAX) then ! Limit internal wind stress due to local domain
           TX2D(I,J)=RHOA*CD*WSX2D(I,J)*WDSMAX/WDS*WDSMAX
           TY2D(I,J)=RHOA*CD*WSY2D(I,J)*WDSMAX/WDS*WDSMAX
          ENDIF
          TXSURF(I,J) = RAMP*(-1.E-3) * 
     +     (+TX2D(I,J)*COS(ANG(I,J))+TY2D(I,J)*SIN(ANG(I,J)))
          TYSURF(I,J) = RAMP*(-1.E-3) * 
     +     (-TX2D(I,J)*SIN(ANG(I,J))+TY2D(I,J)*COS(ANG(I,J)))
          SWRAD(I,J)  = RAMP* SWOBS2D(I,J) / (-4.186E6)
!NG04192011; IT SHOULD ONLY BE RAMPED AT FORCING TIME IN EXTERNAL AND ADV?           PATM(I,J) = PATM(I,J)*RAMP ! NG04092008 Included RAMP for consistency with synop_wind

!      IF(I.EQ.10.AND.J.EQ.10)
!    . WRITE(1004,3000)TIME,I,J,WTSURF(I,J),SWRAD(I,J)
 3000     FORMAT(F10.4,2I5,E12.4)
 4403   CONTINUE
!NG04092008        IF (ISKILL.NE.0.AND.MOD(INTX,ISKILL).EQ.0) Then
!NG04092008         WRITE(1004)TIME,SWOBS2D,WTSURF
!NG04092008        END IF
      ENDIF
!
!svv: end of paste

!NG 2011 FOR VARIABLE ICE THICKNESS AND AREA. START...
      IF (ICEFRICTION.EQ.3) THEN
        IF (THOUR.GE.T2I) THEN
!
          T1I = T2I
          DO 4510 J=1,JM
           DO 4510 I=1,IM
           IF(FSM(I,J).EQ.0)GO TO 4510
           DAREAICE(I,J,1) = DAREAICE(I,J,2)
           DTHICKICE(I,J,1) = DTHICKICE(I,J,2)
 4510    CONTINUE
!   UPDATE DATA FROM 'ifric2dt.inp' data
          ! A.Jordi: Read space-variable ice friction from netcdf file
          NT_T2I=NT_T2I+1
          CALL NETCDF_READ_IFRIC2DT(NT_T2I)
          ! A.Jordi: end
        END IF

        FACT = (THOUR-T1I) / (T2I-T1I)
        DO 4500 J=1,JM
         DO 4500 I=1,IM
          IF(FSM(I,J).EQ.0.0)GO TO 4500
          AREAICE(I,J) = DAREAICE(I,J,1) 
     .           + FACT*(DAREAICE(I,J,2)-DAREAICE(I,J,1))  ! ICE CONCENTRATION (AREA, PERCENT 0-100)
!NG          IF(AREAICE(I,J).LT.70.) AREAICE(I,J)=0.           ! IF AREA IS LESS THAN 70% DO NOT INCLUDE FRICTION
          THICKICE(I,J) = DTHICKICE(I,J,1) 
     .           + FACT*(DTHICKICE(I,J,2)-DTHICKICE(I,J,1))  ! ICE THICKNESS, CM
!          Z0ICE=0.01*THICKICE(I,J)/30 ! NIKURADSE TYPE
          Z0ICE=0.01*THICKICE(I,J)/60 ! MELLOR-KANTHA (1985) TYPE
          VARWIF(I,J)=.16/
     .      ALOG((DZ(1)*0.5)*AMAX1(DT(I,J),1.E-1,WETMIN)/Z0ICE)**2
 4500   CONTINUE
      ENDIF
!NG 2011 FOR VARIABLE ICE THICKNESS AND AREA. ...END

!
! Use SYNOPTIC wind or a local wind for wave computations (hli,03/22/04). 
!-------- SYNOPTIC (WIND SPEED FOR WAVEDYN MODEL by Quamrul QA 10/27/98)
!
!NG04092008 BELOW NOT NEEDED ANY MORE. USE WSX2D, WSY2D FOR WINDS FROM NOW ON.
!NG04092008 CNG02222007      IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN      
!NG04092008       IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN      
!NG04092008         IF(THOUR.LT.T2WV) GOTO 6310
!NG04092008         T1WV=T2WV
!NG04092008         DO J=1,JM
!NG04092008           DO I=1,IM
!NG04092008             WU1(I,J)=WU2(I,J)
!NG04092008             WV1(I,J)=WV2(I,J)
!NG04092008           ENDDO
!NG04092008         ENDDO
!NG04092008         READ (IUT193,END=4018) T2WV
!NG04092008         READ (IUT193) ((WU2(I,J),WV2(I,J),I=1,IM),J=1,JM)
!NG04092008  6310   CONTINUE
!NG04092008         FACT=(THOUR-T1WV)/(T2WV-T1WV)
!NG04092008         DO J=1,JM
!NG04092008           DO I=1,IM
!NG04092008             WU(I,J)=WU1(I,J)+FACT*(WU2(I,J)-WU1(I,J))
!NG04092008             WV(I,J)=WV1(I,J)+FACT*(WV2(I,J)-WV1(I,J))
!NG04092008           ENDDO
!NG04092008         ENDDO
!NG04092008       ENDIF
!
!  FOR WIND WAVE INPUT
!
      IF (WAVEDYN.EQ.'EXTERNAL') THEN
        IF (THOUR.LT.T2WAVE) GOTO 8010
!
        T1WAVE=T2WAVE
!
        DO 8020 J=1,JM
          DO 8020 I=1,IM
            DWVHT(I,J,1)=DWVHT(I,J,2) 
            DWVPD(I,J,1)=DWVPD(I,J,2) 
            DWVDR(I,J,1)=DWVDR(I,J,2) 
 8020   CONTINUE
!
        READ (111,END=5016) T2WAVE
        READ (111) ((DWVHT(I,J,2),I=1,IM),J=1,JM)
        READ (111) ((DWVPD(I,J,2),I=1,IM),J=1,JM)
        READ (111) ((DWVDR(I,J,2),I=1,IM),J=1,JM)
!
 8010  CONTINUE
!
        FACT=(THOUR-T1WAVE)/(T2WAVE-T1WAVE)
!
        DO 8030 J=1,JM
          DO 8030 I=1,IM
            WVHT(I,J)=DWVHT(I,J,1)+FACT*(DWVHT(I,J,2)-DWVHT(I,J,1))
            WVPD(I,J)=DWVPD(I,J,1)+FACT*(DWVPD(I,J,2)-DWVPD(I,J,1))
            WVDR(I,J)=DWVDR(I,J,1)+FACT*(DWVDR(I,J,2)-DWVDR(I,J,1))
!
            WVHT(I,J)=RAMP*WVHT(I,J)
 8030   CONTINUE
      ENDIF
!
      RETURN
!
 4000 FORMAT(8E14.7)
! Modified by Quamrul QA on 2/29/96 **********************************
 5000 Format (12E14.7)
!-----------------------------------------------------------------------
 4010 WRITE(6,4130) THOUR
 4130 FORMAT(//' THE MODEL HAS RUN OUT OF ELEVATION DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4011 WRITE(6,4132) THOUR
 4132 FORMAT(//' THE MODEL HAS RUN OUT OF TEMP-SALINITY DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT',//)
      GOTO 4140
!
 4012 WRITE(6,4013) THOUR
 4013 FORMAT(//' THE MODEL HAS RUN OUT OF DISCHARGE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4014 WRITE(6,4015) THOUR
 4015 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4016 WRITE(6,4017) THOUR
 4017 FORMAT(//' THE MODEL HAS RUN OUT OF METEOROLOGICAL DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4018 WRITE(6,4019) THOUR
 4019 FORMAT(//' THE MODEL HAS RUN OUT OF WIND DATA FOR WAVE AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4021 WRITE(6,4022) THOUR
 4022 FORMAT(//' THE MODEL HAS RUN OUT OF WAVE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4026 WRITE(6,4027) THOUR
 4027 FORMAT(//' THE MODEL HAS RUN OUT OF PRECIP/EVAPOR. DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4036 WRITE(6,4037) THOUR
 4037 FORMAT(//' THE MODEL HAS RUN OUT OF ICE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4514 WRITE(6,4515) THOUR
 4515 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER DATA IN LOOP AT TIME
     .  ',F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
 4516 WRITE(6,4517) THOUR
 4517 FORMAT(//' THE DIFFUSER IN LOOPS HAS WRONG VDDIST VALUES AT TIME',
     .   F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 5016 WRITE(6,5017) THOUR
 5017 FORMAT(//' THE MODEL HAS RUN OUT OF WIND WAVE DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 5012 WRITE(6,5113) THOUR
 5113 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** DISSOLVED TRACER INPUT ****'//)
      GOTO 4140
!
 5213 WRITE(6,5114) THOUR
 5114 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
!
 5313 WRITE(6,5124) THOUR
 5124 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
!
 5014 WRITE(6,5115) THOUR
 5115 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
!
 5214 WRITE(6,5116) THOUR
 5116 FORMAT(//' THE MODEL HAS RUN OUT OF RIVER DISCHARGE DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
!
 7011 WRITE(6,7113) THOUR
 7113 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** DISSOLVED TRACER INPUT ****'//)
      GOTO 4140
!
 7012 WRITE(6,7114) THOUR
 7114 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
!
 7212 WRITE(6,7214) THOUR
 7214 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
!
 7013 WRITE(6,7115) THOUR
 7115 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
!
 7014 WRITE(6,7116) THOUR
 7116 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER TRACER DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
!
 7213 WRITE(6,7117) THOUR
 7117 FORMAT(//' THE MODEL HAS RUN OUT OF OPEN BDY. DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE PARTICLE-BOUND TRACER INPUT ****'//)
      GOTO 4140
!
 7313 WRITE(6,7124) THOUR
 7124 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER MUD DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
!
 7413 WRITE(6,7424) THOUR
 7424 FORMAT(//' THE MODEL HAS RUN OUT OF DIFFUSER SAND DATA AT TIME',
     . 1X,F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '/
     +  5x,'**** NON-COHESIVE SEDIMENT TRANSPORT INPUT ****'//)
      GOTO 4140
!
 6510 WRITE(6,6524) THOUR
 6524 FORMAT(//' THE MODEL HAS RUN OUT OF ELEVATION POTENTIAL '//
     .         'DATA AT TIME ',
     .    F10.4,' HOURS'/,'       REVISE INPUT DECK AND RESUBMIT '//)
      GOTO 4140
!
 4140 CONTINUE
      CLOSE (IUT93)
      CALL SYSTEM('rm gcm_temp*')
!
      STOP
      END
