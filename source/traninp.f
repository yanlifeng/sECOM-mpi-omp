      SUBROUTINE TRANINP(IFLAG)
!
!  READS HYDRO INPUT FROM hqi_tran FILE AND PREPARES FOR SEDIMENT
!  TRANSPORT CALCULATIONS
!
!  FOR COUPLING WITH POM AND WET GRID OUTPUT ONLY
!
!  REVISION DATE:  JULY 30, 1996
!
!**********************************************************************
!
      INCLUDE 'comdeck'
!
!     DECLARE ARRAYS FOR WET GRID OUTPUT
!
      CHARACTER*10 RESTAR
      COMMON /COAST1/ICNT,INDX(MAXWET),JNDX(MAXWET),RESTAR
! 
      DIMENSION  WETGU(MAXWET,KBM1),  WETGV(MAXWET,KBM1),
     .           WETGW(MAXWET,KB),    
     .           WETGAAMX(MAXWET,KBM1), WETGAAMY(MAXWET,KBM1),
     .           WETGKM(MAXWET,KB),   WETGKH(MAXWET,KB),
     .           WETGES(MAXWET),WETGED(MAXWET), 
     .           WETGS(MAXWET,KBM1),WETGT(MAXWET,KBM1) 
!
!
!qa      DIMENSION AAMAX(IM,JM,KB),AAMAY(IM,JM,KB)
!qa      EQUIVALENCE (AAMAX,VH),(AAMAY,VHP)
!
      DIMENSION  KHLPF(IM,JM,KB),
     .   SLPF(IM,JM,KB),    TLPF(IM,JM,KB),
     .   ULPF(IM,JM,KB),    VLPF(IM,JM,KB),   WLPF(IM,JM,KB)
!
       REAL*8 ES(IM,JM),         ED(IM,JM)
!
      REAL KHLPF
!
!QA   CURRENTLY SEDIMENT TRANSPORT MODEL ONLY USES wet_grid TRANSPORT
!QA   INPUT STRUCTURE HAS BEEN MODIFIED BASED ON CONVENTIONAL gcm_tran
!QA   BY QUAMRUL QA 9/2/98
!
      READ (IUTRN) TMIDDLE
!qa      write (*,*)'INTX,TIME,TMIDDLE,ICNT=',INTX,TIME,TMIDDLE,ICNT
!QA   READ (IUTRN) HPRNU
!
      READ (IUTRN) ((WETGU(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGV(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGW(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGAAMX(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGAAMY(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGKH(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGKM(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) (WETGES(I),I=1,ICNT)
      READ (IUTRN) (WETGED(I),I=1,ICNT)
      READ (IUTRN) ((WETGS(I,K),I=1,ICNT),K=1,KBM1)
      READ (IUTRN) ((WETGT(I,K),I=1,ICNT),K=1,KBM1)
!
      DO 330 K=1,KB
        DO 330 I=1,ICNT
          WLPF(INDX(I),JNDX(I),K)=WETGW(I,K)
          KH(INDX(I),JNDX(I),K)=WETGKH(I,K)
          KM(INDX(I),JNDX(I),K)=WETGKM(I,K)
!qa          KM(INDX(I),JNDX(I),K)=WETGKH(I,K)
 330  CONTINUE
!
      DO 340 K=1,KBM1
        DO 340 I=1,ICNT
          ULPF(INDX(I),JNDX(I),K)=WETGU(I,K)
          VLPF(INDX(I),JNDX(I),K)=WETGV(I,K)
          AAMAX(INDX(I),JNDX(I),K)=WETGAAMX(I,K)
          AAMAY(INDX(I),JNDX(I),K)=WETGAAMY(I,K)
          T(INDX(I),JNDX(I),K)    =WETGT(I,K)
          S(INDX(I),JNDX(I),K)    =WETGS(I,K)
 340   CONTINUE
!
!
      DO 350 I=1,ICNT
        ES(INDX(I),JNDX(I))=WETGES(I)
        ED(INDX(I),JNDX(I))=WETGED(I)
 350  CONTINUE
!
!  SET ELEVATION TIME RATE OF CHANGE & INITIAL DEPTH
!
      DO 40 I=1,IM
        DO 40 J=1,JM
!
!  START OF CALC. FOR IFLAG=0
!
          IF (IFLAG.EQ.0) THEN
            ETB(I,J)=ES(I,J)
            ET(I,J)=ES(I,J)
            ETF(I,J)=ES(I,J)
            REWIND (IUTRN)
          ENDIF
!
          EL(I,J)=ES(I,J)
          DT(I,J)=H(I,J)+ES(I,J)
          DETA(I,J)=ED(I,J)
 40   CONTINUE
!
!-----------------------------------------------------------------------
!         CALCULATE HORIZONTAL MASS FLUXES, (H2*U*D) AND (H1*V*D)
!-----------------------------------------------------------------------
!
      DO 10 K=1,KBM1
        DO 10 J=2,JMM1
          DO 10 I=2,IM
            XMFL3D(I,J,K)=ULPF(I,J,K)/DZ(K)
!
!  CALC. VELOCITIES
!
            U(I,J,K)=0.0
            DBAR=0.25*(DT(I,J)+DT(I-1,J))*(H2(I,J)+H2(I-1,J))
            IF (DBAR.GT.0.0) U(I,J,K)=ULPF(I,J,K)/(DZ(K)*DBAR)
 10   CONTINUE
!
      DO 20 K=1,KBM1
        DO 20 J=2,JM
          DO 20 I=2,IMM1
            YMFL3D(I,J,K)=VLPF(I,J,K)/DZ(K)
!
!  CALC. VELOCITIES
!
            V(I,J,K)=0.0
            DBAR=0.25*(DT(I,J)+DT(I,J-1))*(H1(I,J)+H1(I,J-1))
            IF (DBAR.GT.0.0) V(I,J,K)=VLPF(I,J,K)/(DZ(K)*DBAR)
 20   CONTINUE
!
!-----------------------------------------------------------------------
!
!  SET VERTICAL VELOCITIES
!
      DO 30 J=1,JM
        DO 30 K=1,KB
          DO 30 I=1,IM
            W(I,J,K)=WLPF(I,J,K)/ART(I,J) 
 30   CONTINUE
!
!qa          write(*,*)'TMIDDLE,T,S,U,V,=',
!qa     .TMIDDLE,T(20,20,1),S(20,20,1),U(20,20,1),V(20,20,1)
      RETURN
      END
