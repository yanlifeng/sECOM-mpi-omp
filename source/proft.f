      SUBROUTINE PROFT(F,WFSURF,DT2,SW)
!     VERSION(09/27/90)
!
      INCLUDE 'comdeck'
      SAVE
!
      REAL*8 AF(IM,JM,KB),CF(IM,JM,KB),VHF(IM,JM,KB),VHPF(IM,JM,KB)
      REAL*8 FF(IM,JM,KB)
!
      DIMENSION F(IM,JM,KB),WFSURF(IM,JM),DH(IM,JM)
      DIMENSION SW(IM,JM)
      EQUIVALENCE (A,AF),(PROD,CF),(VH,VHF),(AF,FF)
      EQUIVALENCE (TPS,DH)
!
!
!     Modified by Quamrul QA on 5/6/96 *********************************
!
!
!       Modified by Quamrul QA on 2/27/96 *********************************
!       This Version of PROFT is modified to accomodate Shortwave Radiation
!       Penetrated through the Water Column. The Algorithm uses the Classi-
!       fication of Jerlov (1976), but Approximate his Methodology Such that
!       Short Wave Radiation, absorved in the first few meters, is added
!       to the surface boundary condition, and the reminder is attenuated
!       according to RAD = SW*TRC*exp(EXTC*Z)
!
!       Here: RAD is attenuated shortwave radiation
!             SW is shortwave radiation (SW=0.0 when PROFT called for Salinity)
!             TR fraction absorved in surface layer
!             TRC fraction absorved in water column (1.-TR)
!             EXTC Extinction Coefficient
!       TR and EXTC are read from MET DATA segment of run_data input file
!
!
!       FRACTION OF SHORT WAVE RAD PENETRATED THRU WATER COLUMN "TRC"
        TRC = 1. - TR
!pl         if(MOD(INTX,144).EQ.0)write(41,1212)THOUR/24.,EXTC
!         write(41,1212)FLOAT(INTX),THOUR/24.,EXTC,TR,
!     +EXTCOEF(1),EXTCOEF(2)
!pl1212   format(8F10.3)
!
!     Modified by Quamrul QA on 5/6/96 *********************************
!
        UMOLPR = UMOL
!
! Modified by Quamrul QA on 2/27/96 To Separate Total Heatflux from Shortwave
! Note: while computing Salinity SW is zero
!
        DO 25 J = 1, JM
        DO 25 I = 1, IM
            WFSURF(I,J) = WFSURF(I,J)-SW(I,J)
   25   CONTINUE

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
      DO 100 K=2,KBM1
      DO 100 I=2,IMM1
      AF(I,J,K-1)=-DT2*(KH(I,J,K)+UMOLPR*VARUF(I,J))/
     .    (DZ(K-1)*DZZ(K-1)*DH(I,J)*DH(I,J))
      CF(I,J,K)=-DT2*(KH(I,J,K)+UMOLPR*VARUF(I,J))/
     .    (DZ(K)*DZZ(K-1)*DH(I,J)*DH(I,J))
 100  CONTINUE
!
!-------- SURFACE BOUNDARY CONDITIONS - WFSURF -------------------------
      DO 110 J=2,JMM1
      DO 110 I=2,IMM1
      VHF(I,J,1)=AF(I,J,1)/(AF(I,J,1)-1.)
! Modified by Quamrul QA on 2/27/96 SW Penetration **********************
            VHPF(I,J,1) = -DT2 * (WFSURF(I,J)+SW(I,J)*(1.-TRC))
     *      /(-DZ(1)*DH(I,J)) - F(I,J,1)
! End Modified by Quamrul QA on 2/27/96 *********************************
 110  VHPF(I,J,1)=VHPF(I,J,1)/(AF(I,J,1)-1.)
!
!-------- SURFACE BOUNDARY CONDITIONS - FSURF --------------------------
!     DO 110 J=2,JMM1
!     DO 110 I=2,IMM1
!     VHF(I,J,1)=0.
!110  VHPF(I,J,1)=FSURF(I,J)
!
      DO 101 K=2,KBM2
      DO 101 J=2,JMM1
      DO 101 I=2,IMM1
      VHPF(I,J,K)=1./(AF(I,J,K)+CF(I,J,K)*(1.-VHF(I,J,K-1))-1.)
      VHF(I,J,K)=AF(I,J,K)*VHPF(I,J,K)
      VHPF(I,J,K)=(CF(I,J,K)*VHPF(I,J,K-1)-DBLE(F(I,J,K)))*VHPF(I,J,K)
 101  CONTINUE
!
      DO 130 K=1,KB
      DO 130 J=1,JM
      DO 130 I=1,IM
 130  FF(I,J,K)=F(I,J,K)
!
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
 102  FF(I,J,KBM1)=((CF(I,J,KBM1)*VHPF(I,J,KBM2)-FF(I,J,KBM1))
     .   /(CF(I,J,KBM1)*(1.-VHF(I,J,KBM2))-1.))
!
      DO 105 K=2,KBM1
      KI=KB-K
      DO 105 J=2,JMM1
      DO 105 I=2,IMM1
      FF(I,J,KI)=(VHF(I,J,KI)*FF(I,J,KI+1)+VHPF(I,J,KI))
 105  CONTINUE
!
! Modified by Quamrul QA on 5/6/96 *********************************
!       Penetrative Radiation Calculation. Any unattenuated is deposited
!       in the Bottom Layer
!
!     DO 320 K=1,KBM1
!     DO 320 J=1,JM
!     DO 320 I=1,IM
!         RADP(I,J,K)=SW(I,J)*TRC*
!    *               EXP(EXTC*Z(K)*DH(I,J))
! 320     CONTINUE
!
!     
      IF(OPTEXTC.EQ.'VARI') THEN
        DO J=1,JM
          DO K=1,KBM1
            DO I=1,IM
              RADP(I,J,K)=SW(I,J)*TRC*
     *          EXP(EXTC1(I,J)*Z(K)*DH(I,J))
            ENDDO
          ENDDO
        ENDDO     
!        write(1000,*)TIME,OPTEXC(1:7),OPTMBC(1:10)
!        write(1000,*)EXTC1(109,45),DH(109,45),SW(109,45),TRC
!        write(1000,*)Z(1),RADP(109,45,1)
      ELSE
        DO J=1,JM
          DO K=1,KBM1
            DO I=1,IM
              RADP(I,J,K)=SW(I,J)*TRC*
     *          EXP(EXTC*Z(K)*DH(I,J))
            ENDDO
          ENDDO
        ENDDO     
      ENDIF
!
      DO 420 J=1,JM
      DO 420 I=1,IM
  420 RADP(I,J,KB)=0.
!
      Do 210 J = 1, JM
        Do 220 K = 1, KBM1
          Do 200 I = 1, IM
!NG09232011            IF(DH(I,J).GT.0.0)THEN
            IF(DH(I,J).GT.WETMIN)THEN
              F(I,J,K) = FF(I,J,K)-DT2*(RADP(I,J,K)-RADP(I,J,K+1))/
     *                   (DZ(K)*DH(I,J))
            ENDIF
  200     Continue
  220   Continue
  210 Continue
!
! End Modified by Quamrul QA on 5/6/96 *********************************
!
      RETURN
      END
