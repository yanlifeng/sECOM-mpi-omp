      SUBROUTINE STRESS
!  CALC. BOTTOM SHEAR STRESS DUE TO WIND WAVES AND CURRENTS
!  USES GRANT/MADSEN/GLENN WAVE MODEL
!  REVISION DATE:  July 2, 1999
!     This Subroutine uses the CDCALC subroutine rovided by Rich Signell
!     and introduced in this code by Quamrul Ahsan, QA on 7/2/99
!     The routine Calculates effective CD due to current and waves using
!     the lowest layer velocity from ECOM and wave information
!     from WAVESMB or WAVEDON Model (PERIOD = mean wave period of energy 
!     density spectrum, UBOT   = rms value of the maxima of the orbital 
!     velocity near the bottom (m/s))
!     This routine is a simplified form of Grant & Madsen (1979)
!     which assumes colinear waves and currents.  This is described
!     in "Effect of wave-current interaction on steady wind-driven
!     circulation in narrow, shallow embayments", Signell et al,
!     JGR, 1990, v 95, n c6, 9671--9678.
!**************************************************************
!NG 04272010 Changed to make more wavedon like, and for W&D.
!     Common block REFLE needed for Grant & Madsen routines
!NG04272010      Common /REFLE/ ZR, Z0C, USTARC, USTACW, USTARW
      Common /REFLE/ ZR, Z0C, USTARC, USTACW, USTARW, PHICOL ! NG04272010 Include wave-current angle
!     USTARW = the friction velocity due to waves only
!     USTACW = the friction velocity due to combined wave and currents
      Real OMEGA, KBR
      INCLUDE 'comdeck'
      REAL FSMWSTRS(IM,JM)   !NG06092008
!NG12212010 Not used      Data IFIRST /0/
      PI = 3.141593
      TPI= 6.283185
      HPI= 1.570796
!hli
      DO J=1,JM
         DO I=1,IM
            FSMWSTRS(I,J)=FSM(I,J) 
         ENDDO
      ENDDO
      DO J=1,JMM1
         DO I=1,IMM1
            IF(DUMWAD(I,J)+DUMWAD(I+1,J)+DVMWAD(I,J)+DVMWAD(I,J+1).EQ.0)
     .        FSMWSTRS(I,J)=0.0 ! NG12212010 FOR WAD
         ENDDO
      ENDDO
      DO N=1,NUMQBC
         IC=IQC(N)
         JC=JQC(N)
         FSMWSTRS(IC,JC)=0.0
      ENDDO
!NG06082012      DO N=1,NUMEBC
!NG06082012         FSMWSTRS(IETA(N),JETA(N))=0.0 ! NG12212010 ALSO OPEN BOUNDARIES
!NG06082012      ENDDO
!hli
!
!  CALC. TOTAL BOTOM SHEAR STRESS
!
      DO 30 J=2,JMM1
        DO 30 I=2,IMM1
          IF (FSMWSTRS(I,J).EQ.0.0) GOTO 30
!
!  NO SHEAR STRESS IF D < 10 cm
!
!NG            IF (DT(I,J).LT.0.10) THEN ! CNG04272010 W&D
!NG11012013            IF (DT(I,J).LE.AMAX1(0.10,WETMIN)) THEN ! CNG04272010 W&D
            IF (DT(I,J)-WETMIN.LE.0.10) THEN ! CNG11012013 W&D
              TAU(I,J,KB)=0.0
              GOTO 30
            ENDIF
            IF(WVPD(I,J).LE.0.1.OR.WVHT(I,J).LE.0.01) THEN
!NG04272010 FOR SED TRAN AS PER PRAVI'S NOTE IN BOTTAU              TAU(I,J,KB)=10000.*CBC(I,J)*QBAR(I,J)*QBAR(I,J)
              TAU(I,J,KB)=10000.*CBC(I,J)/VARBF(I,J)*QBAR(I,J)*QBAR(I,J)
              GOTO 30
            ENDIF
!
!  CALC. WAVE FREQUENCY (1/s)
!
            OMEGA=TPI/WVPD(I,J) ! CNG04272010 MOVE IT HERE
!NG04272010C
!NG04272010C  CALC. WAVE NUMBER (1/m)
!NG04272010 WILL USE WNUMYOU DISPERSION SOLVER FROM NOW ON, CONSISTENT WITH WAVEDON
            WAVEK=WNUMYOU(OMEGA/TPI,DT(I,J),GRAV)
!NG04272010C  ASSUME OMEGAR=OMEGA
!NG04272010C  CALC. SHALLOW WATER WAVE SPEED (m/s)
!NG04272010C
!NG04272010            C0=SQRT(GRAV*DT(I,J))
!NG04272010C
!NG04272010C  ESTIMATE WAVE NO.
!NG04272010C
!NG04272010            WAVEK0=TPI/(C0*WVPD(I,J))
!NG04272010C
!NG04272010C  ITERATE FOR WAVE NO.
!NG04272010C
!NG04272010            TRAT=TPI*TPI/(GRAV*WVPD(I,J)*WVPD(I,J))
!NG04272010C
!NG04272010C  USE NEWTON'S METHOD
!NG04272010C
!NG04272010            DO 40 N=1,20
!NG04272010              TANHKH=TANH(WAVEK0*DT(I,J))
!NG04272010              FKH=WAVEK0*TANHKH-TRAT
!NG04272010              FPKH=TANHKH+WAVEK0*DT(I,J)*(1.-TANHKH*TANHKH)
!NG04272010              WAVEK=WAVEK0-FKH/FPKH
!NG04272010C
!NG04272010C  CHECK FOR CONVERGENCE
!NG04272010C
!NG04272010              DELTA=ABS(WAVEK/WAVEK0-1.)
!NG04272010              IF (DELTA.LE.0.0001) GOTO 45
!NG04272010              WAVEK0=WAVEK
!NG04272010 40         CONTINUE
!NG04272010 45         CONTINUE
!
!  CALC. BOTTOM ORBITAL AMPLITUDE (m)
!
            ABM=WVHT(I,J)/(2.*SINH(WAVEK*DT(I,J))) + 0.00001
!
!  CALC. BOTTOM ORBITAL VELOCITY (m/s)
!
            UB2=ABM*OMEGA + 0.00001
!
!  HEIGHT OF VELOCITY ABOVE BOTTOM
!
            ZR=DZ(KBM1)*DT(I,J)/2.
            IF(TOR.EQ.'BAROTROPIC') ZR=0.5*DT(I,J) ! FOR BAROTROPIC THINK OF IT AS MID-POINT?
            KBR=30.*Z0B
!
! GET CURRENTS FROM MODEL AT HEIGHT ZR ABOVE BED & CONVERT TO CGS UNITS
!
            URSPD = SQRT(UBAR(I,J,KBM1)*UBAR(I,J,KBM1)+
     +            VBAR(I,J,KBM1)*VBAR(I,J,KBM1)) * 100.
!NG04272010 WE SHOULD ADD, FOR NON-COLINEAR WAVES:
      PHICOL=((-(WVDR(I,J)-90.)-ANG(I,J))*PI/180.-
     .       ATAN2(VBAR(I,J,KBM1),UBAR(I,J,KBM1)))
!
!  CALC. SHEAR VELOCITY DUE TO WAVES & CURRENTS
!  NEAR BOTTOM CURRENT VELOCITY
!  SET MINIMUM CURRENT VELOCITY TO 1 mm/s
!
            URSPD=AMAX1(URSPD,0.1)
!
! CALL GRANT & MADSEN ROUTINE (CGS UNITS) 
! TO GET EFFECTIVE DRAG COEFFICIENT DUE TO WAVES AND CURRENTS
!
            CALL DRGCO(UB2*100,KBR*100.,OMEGA,URSPD,0.,CDE,Z0)
!NG            CALL CDNEW(UB2*100,KBR*100.,OMEGA,URSPD,CDE,Z0)
!
! DO NOT LET DRAG COEFFICIENT BE LESS THAN CBC
!
!NG06092008            CBC(I,J) = AMAX1(CDE,BFRIC)
!           IF (I.EQ.23.AND.J.EQ.26) 
!     .         write (*,*) 'TIME:',TIME,'CDE:',CDE,'CBC',CBC(I,J) ! NG LITTLE CHECK FOR DRAG
!            IF(I.EQ.93.and.J.EQ.15.AND.MOD(INTX,10).EQ.0.) THEN
!               write (1003,'(5(A,F12.9),4(A,F7.3))') 
!     .         'TIME:',TIME,' CDE:',CDE,' CBC',CBC(I,J),
!     .         ' UCUR:',URSPD/100.,' UWAV:',UB2,
!     .         ' DT:',DT(I,J),' WVPD:',WVPD(I,J),
!     .         ' WVHT:',WVHT(I,J), 'WVN:',WAVEK ! NG LITTLE CHECK FOR DRAG
!            ENDIF
            CBC(I,J)=AMIN1(AMAX1(CDE,CBC(I,J)),1.) ! NG05172010 UPDATE, BUT DO NOT ALLOW CBC TO EXCEED UNITY.
            TAU(I,J,KB)=AMAX1(USTARC,USTACW)
            TAU(I,J,KB)=TAU(I,J,KB)*TAU(I,J,KB)
 30   CONTINUE
!
!     FORCING TAU VALUE = 0 AT CELLS CONNECTED BY RIVERS (IC,JC)
!     By Quamrul QA 10/28/98
!
      DO 120 N=1,NUMQBC
        IC=IQC(N)
        JC=JQC(N)
        TAU(IC,JC,KB) = 0.0
120   CONTINUE
      RETURN
      END
!*****************************************************************
      Subroutine DRGCO(UB,KBR,OMEGA,UR,DIA,CD,Z0)
!NG I don't think we should save FW !        SAVE FW
!     ******************************************
!     SUBROUTINE DRGCO          DRAG COEFFICIENT
!     W. D. GRANT                   W.H.O.I.
!     GENERATED                    FEB 14, 1984
!     LATEST REVISION              APR  1, 1984
!     Determines the drag coefficient for mean
!     flow.
!     SUBROUTINES CALLED:
!        WFRIC
!        MOVKBR
!     ******************************************
!     If you want movable bed effects, specify KBR=100. and call
!     this routine with a realistic sand diameter (DIA).
!     If you DON'T want movable bed effects, call this routine
!     with DIA = 0.
      Real KBRAB, KBR
!NG04272010      Common /REFLE/ ZR,z0c,ustarc,ustacw,ustarw
      Common /REFLE/ ZR,z0c,ustarc,ustacw,ustarw,PHICOL ! CNG04272010 Included wave-current angle
      If (KBR.EQ.100.0) Then
        Call MOVKBR(DIA,UB,OMEGA,KBRAB)
        KBR = KBRAB * (UB/OMEGA)
      Else
!NG04272010 KBRAB is the Nikuradse bed roughness Ks over the 
!NG04272010        semi-excursion distance, UBOT(WAVE)/OMEGA
!NG04272010        where UBOT(WAVE)=pi*Ho/(Tsinh(kh)).
!NG04272010 In scientific literature, it is usually noted as (ks/A).
        KBRAB = KBR / (UB/OMEGA)
      End If
      Call WFRIC(KBRAB,FW)
!      WRITE(6,*) 'KBR ',KBR,'   KBRAB ',KBRAB,' UB ',UB,'  OMEGA ',
!     *  OMEGA,
!     * ' FW ',FW
      USTARW = SQRT(0.5*FW) * UB
!                       ------------------------------------------
!                       INITIAL DRAG COEF ESTIMATE
      DELTAW = 2.0 * 0.4 * USTARW / OMEGA
      Z0 = KBR / 30.0
      Z0C = Z0 * ((DELTAW/Z0)**0.8)
!PL PRAVI & QA 09/16/99     CD = (0.4/ALOG(ZR/Z0C)) ** 2.0   
      CD = (0.4/ALOG(ZR*100./Z0C)) ** 2.0   !ZR converted to cm
!      WRITE(6,*) 'INITIAL ESTIMATE OF CD ',CD
!      WRITE(6,*) 'Z0 ',Z0
!     WRITE(6,1000)
!1000 FORMAT(//,' I   CD    USTARC   UAUB    ALPHA  USTACW  DELTAW ',
!    *          '   B      Z0C',/)
!      write(1002,*)'SUBRO DRGCO = ',Z0
!NG      Return
!NG      Entry CDNEW(UB,KBR,OMEGA,UR,CD,Z0)
!                       ------------------------------------------
!                       UPDATE CD
      
!      WRITE(4001,*) 'INSIDE CDNEW---- ZR ',ZR
!      WRITE(4001,*) 'FW ',FW
!      WRITE(4001,*) 'Z0 ',Z0
!       WRITE(4001,*) 'UR ',UR
!       WRITE(4001,*) 'UB ',UB,'  KBR ',KBR,' OMEGA ',OMEGA
!   ----- INSERT SIMPLE VERSION OF DRAG COEFF. SUBROUTINE   -------
!      write(1002,*)'ENTRY CDNEW = ',Z0
      Do 100 I = 1, 10
        CD0 = CD
        USTARC = SQRT(CD) * UR
!        UAUB= (CD/(0.64*FW))*(UR/UB)**2.0
!        ALPHA= (1.0+UAUB)**2.0
!        USTACW= ((0.5*FW*ALPHA)**0.5)*UB
        USTAW2 = 0.5 * FW * UB ** 2
        USTAW = SQRT(USTAW2)
        USTACW = SQRT(USTAW2+USTARC**2)
!NG04272010 ADDED ANGLE FOR NON-COLINEAR WAVES
!NG04272010 TOO BIG        USTACW = SQRT(USTAW2+USTARC**2+(2.*USTAW*USTARC*COS(PHICOL)))
        !write (*,*) USTAW,USTARC,USTACW,PHICOL
        DELTAW = 2.0 * 0.4 * USTACW / OMEGA
        B = 1.0 - (USTARC/USTACW)
!NG04272010 NOW B<0 IS POSSIBLE
!NG04272010        If (B.LT.0) Then
!NG04272010          Write (6,*) 'B LESS THAN 0'
!NG04272010          Stop
!NG04272010        End If
        Z0C = Z0 * ((DELTAW/Z0)**B)
!PL PRAVI & QA 09/16/99     CD = (0.4/ALOG(ZR/Z0C)) ** 2.0   
        CD = (0.4/ALOG(ZR*100./Z0C)) ** 2.0   !ZR converted to cm
        DIFF = ABS(CD-CD0)
        If (DIFF.LE.1.E-07) Go To 110
!       WRITE(6,'(' ITERATION =',I2,'   CD =',F10.6)') I, CD
!        WRITE(6,1010) I, CD, USTARC, UAUB, ALPHA, USTACW, DELTAW,
!     *                B, Z0C
  100 Continue
! 110 WRITE(6,*) 'CD ',CD
  110 Continue
      Return
 5000 Format (I2,8(1X,F7.4))
      End
      Subroutine WFRIC(KBRAB,FW)
!     ******************************************
!     SUBROUTINE WFRIC             WAVE FRICTION
!     W. D. GRANT                   W.H.O.I.
!     GENERATED                    FEB 14, 1984
!     LATEST REVISION              AUG 21, 1986
!     Determines wave friction factor, FW, from
!     relative roughness.
!     ******************************************
      Real KBRAB
!     ------------------------------------------
!     REGION 3: KBRAB .GT. 1.0
      If (KBRAB.GT.1.0) Then ! NG If (A/ks<1)
        FW = 0.23
!     ------------------------------------------
!     REGION 2: KBRAB .LT. 0.08 ! NG If (A/ks>12.5)
      Else If (KBRAB.LT.0.08) Then
        FW = 0.13 * KBRAB ** 0.40
!                       ------------------------------------------
!     REGION 1: 0.08 .LE. KBRAB .LE. 1.0
      Else ! NG For intermediate A/ks
        FW = 0.23 * KBRAB ** 0.62
      End If
      Return
      End
      Subroutine MOVKBR(DIA,UB,OMEGA,KBRAB)
!     ******************************************
!     SUBROUTINE MOVKBR    MOVEABLE BED ROUGHNESS
!     W. D. GRANT                   W.H.O.I.
!     GENERATED                    APR  1, 1984
!     LATEST REVISION              APR  1, 1984
!     Calculates the relative roughness over
!     moveable beds approximately after Grant &
!     Madsen 1982.
!     All units are CGS.
!     This subroutine assumes quartz sand.
!     ******************************************
      Real KBRAB, KRIPAB, KBREDAB
!     ------------------------------------------
!     CONSTANTS
      S = 2.65
      G = 980.6
      VISC = 0.011
!     ------------------------------------------
!     COMMON FACTORS
      SUBWT = (S-1.0) * G * DIA
      SSTAR = (DIA/(4.0*VISC)) * SQRT(SUBWT)
      PHICR = 0.054
      AB = UB / OMEGA
      KBRAB = DIA / AB
      RE = UB * AB / VISC
      TRANRO = 0.56 / SQRT(RE)
      If (KBRAB.LT.TRANRO) Then
        FW = 0.09 * (1.0/(RE**0.2))
      Else
        Call WFRIC(KBRAB,FW)
      End If
      TAUWB = 1.0 * 0.5 * FW * (UB**2.0)
      PHIPRI = TAUWB / SUBWT
      RATIO = PHIPRI / PHICR
      BREAK = 1.8 * (SSTAR**0.6)
      If (RATIO.LE.BREAK) Then
        STEEP = 0.16 * (1.0/(RATIO**0.04))
        HTAB = 0.22 * (1.0/(RATIO**0.16))
      Else
        STEEP = 0.28 * ((SSTAR**0.6)*(1.0/RATIO))
        HTAB = 0.48 * ((SSTAR**0.8)*(1.0/(RATIO**1.5)))
      End If
      KRIPAB = 30.0 * HTAB * STEEP
      KBREDAB = 27.2 * KBRAB * ((SQRT(RATIO)-0.7)**2.0)
      KBRAB = KRIPAB + KBREDAB
      Write (6,5000) KRIPAB, KBREDAB, KBRAB
      Return
 5000 Format (/,5X,' KRIPAB =',F10.6,/,5X,' KBREDAB =',F10.6,/,5X,
     *    '   KBRAB =',F10.6)
      End
