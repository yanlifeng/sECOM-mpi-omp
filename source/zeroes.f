      SUBROUTINE ZEROES(ADVUA,ADVVA,ADVUU,ADVVV,DRHOX,DRHOY,TRNU,TRNV)
!     VERSION(10/02/90)
!
      INCLUDE 'comdeck'
      SAVE
!
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
      DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
!
      alat=-999.9
      alon=-999.9
!NG12162010 ADDED INITIALIZATION FOR QEVAP AND QPREC TO AVOID IF STATEMENTS IN VERTVL AND EXTRNL BETWEEN 2D and 1D OPTIONS
      QPREC=0.0
      QEVAP=0.0
      DQPREC(1)=0.0
      DQPREC(2)=0.0
      DQEVAP(1)=0.0
      DQEVAP(2)=0.0
!
! Give grid rows the same static first-touch placement as the solver.
! No values are initialized that were not already set by the reference.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K)
!$OMP DO SCHEDULE(STATIC)
      DO 250 J=1,JM
      DO 250 K=1,KB
      DO 250 I=1,IM
      T(I,J,K)=0.0
      TB(I,J,K)=0.0
      S(I,J,K)=0.0
      SB(I,J,K)=0.0
      RHO(I,J,K)=0.0
      RMEAN(I,J,K)=0.0
      TMEAN(I,J,K)=0.0
      SMEAN(I,J,K)=0.0
      XMFL3D(I,J,K)=0.0
      YMFL3D(I,J,K)=0.0
      A(I,J,K)=0.0
      C(I,J,K)=0.0
      VH(I,J,K)=0.0
      VHP(I,J,K)=0.0
      PROD(I,J,K)=0.0
      DTEF(I,J,K)=0.0
      U(I,J,K)=0.0
      UB(I,J,K)=0.0
      V(I,J,K)=0.0
      VB(I,J,K)=0.0
      UF(I,J,K)=0.0
      VF(I,J,K)=0.0
      W(I,J,K)=0.0
      DRHOX(I,J,K)=0.0
      DRHOY(I,J,K)=0.0
      Q2B(I,J,K)=0.0
      Q2(I,J,K)=0.0
      Q2LB(I,J,K)=0.0
      Q2L(I,J,K)=0.0
      L(I,J,K)=0.0
      KH(I,J,K)=0.0
      KM(I,J,K)=0.0
      KQ(I,J,K)=0.0
!
      CONC1(I,J,K)=0.0
      CONC1B(I,J,K)=0.0
      ARCC1(I,J,K)=0.0
      CMEAN1(I,J,K)=0.0
        DKRATE(I,J,K)=0.0
!
      WSET1(I,J,K)=0.0
      WSET2(I,J,K)=0.0
      CSED1(I,J,K)=0.0
      CSED2(I,J,K)=0.0
      CHEM1(I,J,K)=0.0
      CHEM2(I,J,K)=0.0
      ARCSED1(I,J,K)=0.0
      ARCSED2(I,J,K)=0.0
!pl 9/15/98
      ARCTAU(I,J,K)=0.0
      TAU(I,J,K)=0.0
  250 CONTINUE
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
      DO 640 J=1,JM
      DO 640 I=1,IM
!NG04142008
      WCT1BOT(I,J)=0.0
      WCT2BOT(I,J)=0.0
      CHEMBOT1(I,J)=0.0
      CHEMBOT2(I,J)=0.0
!NG04142008
      UAF(I,J)=0.0
      UA(I,J)=0.0
      UAB(I,J)=0.0
      VAF(I,J)=0.0
      VA(I,J)=0.0
      VAB(I,J)=0.0
      ELF(I,J)=0.0
      EL(I,J)=0.0
      ELB(I,J)=0.0
      TPS(I,J)=0.0
      ETF(I,J)=0.0
      ET (I,J)=0.0
      ETB(I,J)=0.0
      UTF(I,J)=0.0
      UTB(I,J)=0.0
      VTF(I,J)=0.0
      VTB(I,J)=0.0
      EGF(I,J)=0.0
      EGB(I,J)=0.0
      PATM(I,J)=0.0
      DPATM(I,J,1)=0.0
      DPATM(I,J,2)=0.0
      WTSURF(I,J)=0.0
      WSSURF(I,J)=0.0
      TXSURF(I,J)=0.0
      TYSURF(I,J)=0.0
      TXBOT(I,J)=0.0
      TYBOT(I,J)=0.0
      VARWIF(I,J)=0.0 ! ICE-WATER FRICTION
      TXICE(I,J)=0.0 ! NG FOR ICE
      TYICE(I,J)=0.0 ! NG FOR ICE
      AREAICE(I,J)=0.0 ! NG FOR ICE
      THICKICE(I,J)=0.0 ! NG FOR ICE
      TRNU(I,J)=0.0
      TRNV(I,J)=0.0
      CURV42D(I,J)=0.0
      ADVUA(I,J)=0.0
      ADVVA(I,J)=0.0
      ADVUU(I,J)=0.0
      ADVVV(I,J)=0.0
      FLUXUA(I,J)=0.0
      FLUXVA(I,J)=0.0
!
      WCSURF(I,J)=0.0
!
!NG04092008 Added 2D winds and met
      WSX2D(I,J)=0.0
      DWSX2D(I,J,1)=0.0
      DWSX2D(I,J,2)=0.0
      WSY2D(I,J)=0.0
      DWSY2D(I,J,1)=0.0
      DWSY2D(I,J,2)=0.0
      TX2D(I,J)=0.0
      DTX2D(I,J,1)=0.0
      DTX2D(I,J,2)=0.0
      TY2D(I,J)=0.0
      DTY2D(I,J,1)=0.0
      DTY2D(I,J,2)=0.0
      SWOBS2D(I,J)=0.0
      DSWOBS2D(I,J,1)=0.0
      DSWOBS2D(I,J,2)=0.0
      AIRTM2D(I,J)=0.0
      DAIRTM2D(I,J,1)=0.0
      DAIRTM2D(I,J,2)=0.0
      RELHU2D(I,J)=0.0
      DRELHU2D(I,J,1)=0.0
      DRELHU2D(I,J,2)=0.0
      CLDFRC2D(I,J)=0.0
      DCLD2D(I,J,1)=0.0
      DCLD2D(I,J,2)=0.0
      EXTC1(I,J)=0.0
      DEXTC2D(I,J,1)=0.0
      DEXTC2D(I,J,2)=0.0
      QPREC2D(I,J)=0.0
      DQPREC2D(I,J,1)=0.0
      DQPREC2D(I,J,2)=0.0
      QEVAP2D(I,J)=0.0
      DQEVAP2D(I,J,1)=0.0
      DQEVAP2D(I,J,2)=0.0
!NG04092008 Added 2D winds and met
!NG04142008 Added waves
      WVPD(I,J)=0.0
      DWVPD(I,J,1)=0.0
      DWVPD(I,J,2)=0.0
      WVHT(I,J)=0.0
      DWVHT(I,J,1)=0.0
      DWVHT(I,J,2)=0.0
      WVDR(I,J)=0.0
      DWVDR(I,J,1)=0.0
      DWVDR(I,J,2)=0.0
!NG04142008 Added waves
!NG2011 Ice
      AREAICE(I,J)=0.0
      DAREAICE(I,J,1)=0.0
      DAREAICE(I,J,2)=0.0
      THICKICE(I,J)=0.0
      DTHICKICE(I,J,1)=0.0
      DTHICKICE(I,J,2)=0.0
!NG2011 Ice
  640 CONTINUE
!$OMP END DO
!
!$OMP DO SCHEDULE(STATIC)
      DO 100 J=1,JM
      DO 100 K=1,KB
      DO 100 I=1,IM
      ARCU (I,J,K)=0.0
      ARCV (I,J,K)=0.0
      ARCUX(I,J,K)=0.0
      ARCVX(I,J,K)=0.0
      ARCS (I,J,K)=0.0
      ARCT (I,J,K)=0.0
      ARCW (I,J,K)=0.0
      ARCKH(I,J,K)=0.0
      ARCHEM1(I,J,K)=0.0
      ARCHEM2(I,J,K)=0.0
      ARCWS1(I,J,K)=0.0
      ARCWS2(I,J,K)=0.0
      ARCSED1(I,J,K)=0.0
      ARCSED2(I,J,K)=0.0
      ARCTAU(I,J,K)=0.0
      ARCAM(I,J,K)=0.0
      ARCKM(I,J,K)=0.0
      ARCQ2(I,J,K)=0.0
         TAU(I,J,K)=0.0
 100  CONTINUE
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
      DO 110 J=1,JM
      DO 110 I=1,IM
      ARCPATM(I,J)=0.0
      ARCQPREC(I,J)=0.0
      ARCQEVAP(I,J)=0.0
      ARCCLD(I,J)=0.0
      ARCAIRT(I,J)=0.0
      ARCRHUM(I,J)=0.0
      ARCSWOBS(I,J)=0.0
      ARCWSX2D(I,J)=0.0
      ARCWSY2D(I,J)=0.0
      ARCWVPD(I,J)=0.0
      ARCWVHT(I,J)=0.0
      ARCWVDR(I,J)=0.0
 110  ARCET(I,J)=0.0
!$OMP END DO
!$OMP END PARALLEL
!
      DO 200 N=1,EPTS
      ESAVE(N)=0.0
      TXSAVE(N)=0.0
      TYSAVE(N)=0.0
      QPRECSAVE(N)=0.0
      QEVAPSAVE(N)=0.0
      PATMSAVE(N)=0.0
      CLDSAVE(N)=0.0
      AIRTSAVE(N)=0.0
      RHUMSAVE(N)=0.0
      SWOBSSAVE(N)=0.0
      WSX2DSAVE(N)=0.0
      WSY2DSAVE(N)=0.0
      WVPDSAVE(N)=0.0
      WVHTSAVE(N)=0.0
      WVDRSAVE(N)=0.0
 200  CONTINUE
      DO 210 N=1,VPTS
      DZSAVE(N)=0.0
      THSAVE(N)=0.0
      DO 210 K=1,KB
      UZSAVE(N,K)=0.0
      VZSAVE(N,K)=0.0
      SZSAVE(N,K)=0.0
      TZSAVE(N,K)=0.0
!
      WS1SAVE(N,K)=0.0
      WS2SAVE(N,K)=0.0
      C1SAVE(N,K)=0.0
      C2SAVE(N,K)=0.0
        THSAVE(N)  =0.0
      TAUSAVE(N,K)=0.0
      P1SAVE(N,K)=0.0
      P2SAVE(N,K)=0.0
!
      C1ZSAVE(N,K)=0.0
!
      KMZSAVE(N,K)=0.0
      KHZSAVE(N,K)=0.0
      AMZSAVE(N,K)=0.0
!
 210  CONTINUE
      DO 220 K=1,KB
      DO 220 N=1,FPTS
      CCFLUX(N,K)=0.0
      CS1FLUX(N,K)=0.0
      CS2FLUX(N,K)=0.0
 220  CONTINUE
      ESUM  =0.0
      TKE   =0.0
      APE   =0.0
      VOLUME=0.0
      VSTOR =0.0
      EM    =0.0
      APEC  =0.0
      TSUM  =0.0
      SSUM  =0.0
!qa   Background Extinction Coefficient for Bacterial decay rate computaion
      CLDX  =0.5  !initilaze CLDX (hli,06/12/02)
      EXTC  =1.2
!hli Moved here from wavedon.f (hli,12/04/03)
      NR=0
!hli      NWAVEOUT=0  !(hli,12/30/03)
!hli      NSTRSOUT=0  !(hli,01/12/04)
!
! ------- SET CONSPLT, CONSTSR & CONSTRN TO TRUE FOR FIRST RUN THROUGH --------
      CONSPLT=.TRUE.
      CONSTSR=.TRUE.
      CONSTRN=.TRUE.
!
!  FOR PARTICLE TRACKING
!
      IF (PARTICLE.EQ.'INCLUDE') THEN
       DO 1110 MM=1,NPART
         DO 1110 NN=1,NCONV
           DO 1110 LL=1,NSOURCE
            XP(LL,MM,NN)=0.0
            YP(LL,MM,NN)=0.0
            ZP(LL,MM,NN)=0.0
            IP(LL,MM,NN)=0
            JP(LL,MM,NN)=0
            KP(LL,MM,NN)=0
            INOUT(LL,MM,NN)=0
 1110   CONTINUE
!
        DO 1120 K=1,KBM1
          DO 1130 N=1,DBCM
            PZERO1(N,K)=0.0
 1130     CONTINUE
          DO 1140 N=1,EBCM
            PZERO3(N,K)=0.0
 1140     CONTINUE
 1120   CONTINUE
!
        DO 1150 N=1,QBCM
          PZERO2(N)=0.0
 1150   CONTINUE
      ENDIF
!
      RETURN
      END
