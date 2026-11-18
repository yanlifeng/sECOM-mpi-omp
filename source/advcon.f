      SUBROUTINE ADVCON(FB,F,FMEAN,DTI2,FF,FDIF,TSDIS,TSBDRY,PSLOAD)
!
!  FOR TRACER, CONSERVATIVE AND NON-CONSERVATIVE
!
      INCLUDE 'comdeck'
      SAVE
!
!-----------------------------------------------------------------------
!     THIS SUBROUTINE INTEGRATES CONSERVATIVE CONSTITUENT EQUATIONS
!-----------------------------------------------------------------------
!
      DIMENSION FB(IM,JM,KB),F(IM,JM,KB),FF(IM,JM,KB)
      DIMENSION FMEAN(IM,JM,KB),FBSTART(IM,JM,KB),ETA(IM,JM)
      DIMENSION XMFLUX(IM,JM,KB),YMFLUX(IM,JM,KB),ZZFLUX(IM,JM,KB)
      DIMENSION XFLUX(IM,JM,KB),YFLUX(IM,JM,KB),ZFLUX(IM,JM,KB)
!qa      DIMENSION FDIF(DBCM),TSDIS(QBCM),TSBDRY(EBCM,KBM1)
      DIMENSION FDIF(DBCM,KBM1),TSDIS(QBCM),TSBDRY(EBCM,KBM1)
      DIMENSION PSLOAD(PSCM)
      EQUIVALENCE (XFLUX,A),(YFLUX,C),(ZFLUX,PROD)
      EQUIVALENCE (FBSTART,VH),(ETA,VHP)
!
      IF(SCHEME .EQ. 'CENTRAL   ')THEN
      DO  9 J=1,JM
      DO  9 I=1,IM
      F (I,J,KB)=F (I,J,KBM1)
    9 FB(I,J,KB)=FB(I,J,KBM1)
!
!
!-------- DO ADVECTION FLUXES AS CENTER DIFFERENCES --------------------
      DO 22 K=1,KBM1
      DO 22 J=2,JM
      DO 22 I=2,IM
      XFLUX(I,J,K)=0.5*(F(I-1,J,K)+F(I,J,K)) * XMFL3D(I,J,K)
   22 YFLUX(I,J,K)=0.5*(F(I,J-1,K)+F(I,J,K)) * YMFL3D(I,J,K)
!
!
      DO 236 N=1,NUMEBCTR
      IE=ITRED(N)
      JE=JTRED(N)
      IC=ITREC(N)
      JC=JTREC(N)
!
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
!
!-------- EAST SIDE ----------------------------------------------------
      DO 21 K=1,KBM1
      UUP  =0.5*(XMFL3D(IE,JE,K)+ABS(XMFL3D(IE,JE,K)))
      UDOWN=0.5*(XMFL3D(IE,JE,K)-ABS(XMFL3D(IE,JE,K)))
   21 XFLUX(IE,JE,K)=UUP*FB(IC,JC,K)+UDOWN*FB(IE,JE,K)
      GO TO 236
!
      ELSE IF (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
!
!-------- WEST SIDE ----------------------------------------------------
      DO 222 K=1,KBM1
      UUP  =0.5*(XMFL3D(IC,JC,K)+ABS(XMFL3D(IC,JC,K)))
      UDOWN=0.5*(XMFL3D(IC,JC,K)-ABS(XMFL3D(IC,JC,K)))
  222 XFLUX(IC,JC,K)=UUP*FB(IE,JE,K)+UDOWN*FB(IC,JC,K)
      GO TO 236
!
      ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
!
!-------- NORTH SIDE ---------------------------------------------------
      DO 23 K=1,KBM1
      VUP  =0.5*(YMFL3D(IE,JE,K)+ABS(YMFL3D(IE,JE,K)))
      VDOWN=0.5*(YMFL3D(IE,JE,K)-ABS(YMFL3D(IE,JE,K)))
   23 YFLUX(IE,JE,K)=VUP*FB(IC,JC,K)+VDOWN*FB(IE,JE,K)
      GO TO 236
!
      ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
!
!-------- SOUTH SIDE ---------------------------------------------------
      DO 24 K=1,KBM1
      VUP  =0.5*(YMFL3D(IC,JC,K)+ABS(YMFL3D(IC,JC,K)))
      VDOWN=0.5*(YMFL3D(IC,JC,K)-ABS(YMFL3D(IC,JC,K)))
   24 YFLUX(IC,JC,K)=VUP*FB(IE,JE,K)+VDOWN*FB(IC,JC,K)
!
!-------- DONE ---------------------------------------------------------
      END IF
  236 CONTINUE
!
      DO 121 J=2,JMM1
      DO 121 I=2,IMM1
  121 FF(I,J,1)=-.5*DZR(1)*(F(I,J,1)+F(I,J,2))*W(I,J,2)*ART(I,J)
      DO 130 J=2,JMM1
      DO 130 K=2,KBM1
      DO 130 I=2,IMM1
  130 FF(I,J,K)=.5*DZR(K)*((F(I,J,K-1)+F(I,J,K))*W(I,J,K)
     .                    -(F(I,J,K)+F(I,J,K+1))*W(I,J,K+1))*ART(I,J)
!
!************* ADJUST FOR RIVER/WALL FLUXES *********************
!     DO 50 N=1,NUMQBCTR
!     ID=ITRQD(N)
!     JD=JTRQD(N)
!     IC=ITRQC(N)
!     JC=JTRQC(N)
      DO 50 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 60 K=1,KBM1
   60         XFLUX(IC,JC,K)=F(IC,JC,K)*XMFL3D(IC,JC,K)
            ELSE
              DO 70 K=1,KBM1
   70         XFLUX(ID,JD,K)=F(IC,JC,K)*XMFL3D(ID,JD,K)
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 80 K=1,KBM1
   80         YFLUX(IC,JC,K)=F(IC,JC,K)*YMFL3D(IC,JC,K)
            ELSE
              DO 90 K=1,KBM1
   90         YFLUX(ID,JD,K)=F(IC,JC,K)*YMFL3D(ID,JD,K)
            ENDIF
      ENDIF
   50 CONTINUE
!
!-------- ADD NET ADVECTION FLUXES & THEN STEP FORWARD IN TIME --------
!
      DO 141 J=2,JMM1
      DO 141 K=1,KBM1
      DO 141 I=2,IMM1
      FF(I,J,K)=FF(I,J,K)
     .              +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K)
      FF(I,J,K)=(FB(I,J,K)*(H(I,J)+ETB(I,J))*ART(I,J)-
     +             DTI2*FF(I,J,K))/((H(I,J)+ETF(I,J))*ART(I,J))
!
!qa     This is wrong place to put Decay rate Quamrul QA 11/12/01
!
!qa        FF(I,J,K)=FF(I,J,K)-DTI*DKRATE(I,J,K)*F(I,J,K)

 141  CONTINUE
!
      ELSE
!
      DO 10 J=1,JM
      DO 10 I=1,IM
      ETA(I,J)=ETB(I,J)
   10 FB(I,J,KB)=FB(I,J,KBM1)
      DO 11 K=1,KB
      DO 11 J=1,JM
      DO 11 I=1,IM
      FBSTART(I,J,K)=FB(I,J,K)
      XMFLUX(I,J,K)=XMFL3D(I,J,K)
      YMFLUX(I,J,K)=YMFL3D(I,J,K)
      ZZFLUX(I,J,K)=W(I,J,K)
   11 CONTINUE
!-----------------------------------------------------------------------
!            SMOLARKIEWICZ'S SCHEME  
!-----------------------------------------------------------------------
      LOOP=0
 5050 LOOP= LOOP+1
!
!-------- DO ADVECTION FLUXES, FIRST AS UPWIND  -------------
      DO 20 K=1,KBM1
      DO 20 J=2,JM
      DO 20 I=2,IM
      XFLUX(I,J,K)=AMAX1(0.,XMFLUX(I,J,K))*FBSTART(I-1,J,K)+
     .             AMIN1(0.,XMFLUX(I,J,K))*FBSTART(I,J,K)
      YFLUX(I,J,K)=AMAX1(0.,YMFLUX(I,J,K))*FBSTART(I,J-1,K)+
     .             AMIN1(0.,YMFLUX(I,J,K))*FBSTART(I,J,K)
   20 CONTINUE
      DO 120 J=2,JMM1
      DO 120 I=2,IMM1
      ZFLUX(I,J,1)=0.0
  120 ZFLUX(I,J,KB)=0.0
      DO 125 J=2,JMM1
      DO 125 K=2,KBM1
      DO 125 I=2,IMM1
      ZFLUX(I,J,K)=AMAX1(0.,ZZFLUX(I,J,K))*FBSTART(I,J,K)+
     .             AMIN1(0.,ZZFLUX(I,J,K))*FBSTART(I,J,K-1)
  125 ZFLUX(I,J,K)=ZFLUX(I,J,K)*ART(I,J)
!
!----- ADD NET ADVECTIVE FLUXES & STEP FORWARD IN TIME 
!
!hli      write(55,*) 'Before0, FF = ',FF(156,119,1),'DKRATE = ',
!hli     &   DKRATE(156,119,1)
!hli      write(55,*) 'Before0, XFLUX = ',XFLUX(157,119,1),
!hli     &   XFLUX(156,119,1)
!hli      write(55,*) 'Before0, YFLUX = ',YFLUX(156,120,1),
!hli     &   YFLUX(156,119,1)
!hli      write(55,*) 'Before0, ZFLUX = ',ZFLUX(156,119,1),
!hli     &   ZFLUX(156,119,2)
!hli      write(55,*) 'Before0, FBSTART = ',FBSTART(156,119,1),
!hli     &   ZFLUX(156,119,2)
!hli      write(55,*) 'Before0, ETA & ETF = ',ETA (156,119),
!hli     &   ETF  (156,119)
!hli      write(55,*) 'Before0, DZR,H,ART = ',DZR(1),H(156,119),
!hli     &   ART  (156,119)
      DO 140 J=2,JMM1
      DO 140 K=1,KBM1
      DO 140 I=2,IMM1
      FF(I,J,K)=     XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K)
     .      +DZR(K)*(ZFLUX(I,J,K)  -ZFLUX(I,J,K+1))
      FF(I,J,K)=(FBSTART(I,J,K)*(H(I,J)+ETA(I,J))*
     +    ART(I,J)-DTI2*FF(I,J,K))/((H(I,J)+ETF(I,J))*ART(I,J))
!
!qa     This is wrong place to put Decay rate Quamrul QA 11/12/01
!
!qa        FF(I,J,K)=FF(I,J,K)-DTI*DKRATE(I,J,K)*F(I,J,K)

 140  CONTINUE
!hli      write(55,*) 'After0, FF = ',FF(156,119,1),'DKRATE = ',
!hli     &   DKRATE(156,119,1)
!
!----- ANTIDIFFUSION VELOCITY -----------------------------------------
!
      IF(SCHEME .EQ. 'UPWIND    ' .OR. LOOP .EQ. 2) GO TO 5051
!
      CALL ANTIDIF(XMFLUX,YMFLUX,ZZFLUX,FB,FF,DTI2,TSDIS,TSBDRY)
!
      DO 146 J=1,JM
      DO 146 K=1,KB
      DO 146 I=1,IM
      ETA(I,J)=ETF(I,J)
      FBSTART(I,J,K)=FF(I,J,K)
  146 CONTINUE
      GO TO 5050
!
      END IF
!-----------------------------------------------------------------------
!            SMOLARKIEWICZ'S SCHEME END
!-----------------------------------------------------------------------
 5051 CONTINUE
!
!-------- ADD DIFFUSIVE FLUXES (HORIZONTAL) ----------------------------
!     DO 30 K=1,KB   ! SHALLOW SEAS
!     DO 30 J=1,JM
!     DO 30 I=1,IM
!  30 FB(I,J,K)=FB(I,J,K)-FMEAN(I,J,K)
      DO 40 K=1,KBM1
      DO 40 J=2,JM
      DO 40 I=2,IM
      XMFLUX(I,J,K)=.5*(AAM(I,J,K)+AAM(I-1,J,K))
      YMFLUX(I,J,K)=.5*(AAM(I,J,K)+AAM(I,J-1,K))
   40 CONTINUE
      DO 100 J=2,JM
      DO 100 K=1,KBM1
      DO 100 I=2,IM
      XFLUX(I,J,K)=
!NG04272010 W&D     .    -XMFLUX(I,J,K)/HPRNU*(H(I,J)+H(I-1,J))
     .    -XMFLUX(I,J,K)/HPRNU*(H(I,J)+ETB(I,J)+H(I-1,J)+ETB(I-1,J))
     . *(FB(I,J,K)-FB(I-1,J,K))*DUM(I,J)*DUMWAD(I,J)/(H1(I,J)+H1(I-1,J)) ! Not in original NYHOPS
     .    *0.5*(H2(I,J)+H2(I-1,J))
      YFLUX(I,J,K)=
!NG04272010 W&D     .    -YMFLUX(I,J,K)/HPRNU*(H(I,J)+H(I,J-1))
     .    -YMFLUX(I,J,K)/HPRNU*(H(I,J)+ETB(I,J)+H(I,J-1)+ETB(I,J-1))
     . *(FB(I,J,K)-FB(I,J-1,K))*DVM(I,J)*DVMWAD(I,J)/(H2(I,J)+H2(I,J-1)) ! Not in original NYHOPS
     .    *0.5*(H1(I,J)+H1(I,J-1))
  100 CONTINUE
!     DO 110 K=1,KB   ! SHALLOW SEAS
!     DO 110 J=1,JM
!     DO 110 I=1,IM
! 110 FB(I,J,K)=FB(I,J,K)+FMEAN(I,J,K)
!
!------------- ADJUST FOR RIVER/WALL FLUXES ---------------
!     DO 51 N=1,NUMQBCTR
!     ID=ITRQD(N)
!     JD=JTRQD(N)
!     IC=ITRQC(N)
!     JC=JTRQC(N)
      DO 51 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 61 K=1,KBM1
   61         XFLUX(IC,JC,K)=0.0
            ELSE
              DO 71 K=1,KBM1
   71         XFLUX(ID,JD,K)=0.0
            ENDIF
      ELSE
            IF(JC.GT.JD) THEN
              DO 81 K=1,KBM1
   81         YFLUX(IC,JC,K)=0.0
            ELSE
              DO 91 K=1,KBM1
   91         YFLUX(ID,JD,K)=0.0
            ENDIF
      ENDIF
   51 CONTINUE
!
!----- ADD NET HORIZONTAL FLUXES & STEP FORWARD IN TIME (CENTER DIFF.)
!hli      write(55,*) 'Before, FF = ',FF(156,119,1),'DKRATE = ',
!hli     &   DKRATE(156,119,1)
      DO 145 J=2,JMM1
      DO 145 K=1,KBM1
      DO 145 I=2,IMM1
      FF(I,J,K)=FF(I,J,K)-DTI2*(
     .              +XFLUX(I+1,J,K)-XFLUX(I,J,K)
     .              +YFLUX(I,J+1,K)-YFLUX(I,J,K))
     .                    /((H(I,J)+ETF(I,J))*ART(I,J))

!qa     Decay rates are now Function of T,S,Io and Vs See ecom3d.f by 
!qa     Quamrul QA 11/12/01

  145 FF(I,J,K) = FF(I,J,K) - DTI2*FB(I,J,K)*DKRATE(I,J,K)
!hli      write(55,*) 'After , FF = ',FF(156,119,1),'DKRATE = ',
!hli     &   DKRATE(156,119,1)
!-----------------------------------------------------------------------
!         IMPOSE TRACER FLUX BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
      DO 150 N=1,NUMDBCTR
      ID=IDD(N)
      JD=JDD(N)
      DO 150 K=1,KBM1
!qa      FF(ID,JD,K)=FF(ID,JD,K)+DTI2*FDIF(N)*QDIFF(N)*RAMP*
      FF(ID,JD,K)=FF(ID,JD,K)+DTI2*FDIF(N,K)*QDIFF(N)*RAMP*
     .   VDDIST(N,K)/100./DZ(K)/((H(ID,JD)+ETF(ID,JD))*ART(ID,JD))
      IF(H(ID,JD)+ETF(ID,JD).LE.WETMIN.AND.QDIFF(N).GT.0.) THEN ! IF DIFFUSER ON DRY CELL (NG04132011)
       FF(ID,JD,K)=FDIF(N,K)
      ENDIF
  150 CONTINUE
!hli      write(55,*) 'After1, FF = ',FF(156,119,1),'DKRATE = ',
!hli     &   DKRATE(156,119,1)
!
!  POINT SOURCE LOADS
!
      DO 160 N=1,NUMPSTR
        IPS=IPSTR(N)
        JPS=JPSTR(N)
        KPS=KPSTR(N)
!
!  MASS LOADING (INPUT AS kg/s, CONVERTED HERE FOR CONC IN mg/l) 
!
        IF (OPTPSTR.EQ.'MASS') THEN
          FF(IPS,JPS,KPS)=FF(IPS,JPS,KPS)+DTI2*(1000.*PSLOAD(N))*RAMP 
     .     /DZ(KPS)/((H(IPS,JPS)+ETF(IPS,JPS))*ART(IPS,JPS))
        ELSE
          IF (PSLOAD(N).GT.0.0) THEN
            FF(IPS,JPS,KPS)=PSLOAD(N)
!
!  RESET INPUT CONC. TO ZERO AFTER INSTANTANEOUS SPECIFICATION
!
            PSLOAD(N)=0.0
          ENDIF
        ENDIF
 160  CONTINUE
!
      RETURN
      END
