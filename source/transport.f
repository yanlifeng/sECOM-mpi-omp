      SUBROUTINE TRANSPORT 
!     VERSION(12/18/92)
!
      INCLUDE 'comdeck'
      SAVE ICNTRAN,NHWQI ! NG03222007, bug # NG02182009 NHWQI bug for PGI7.2.
!
!**********************************************************************C
! **  THIS SUBROUTINE CALCULATES TIME AVERAGED MASS TRANSPORT 
!          FIELDS(EULERIAN+STOKES RESIDUAL FLOWS) AND OTHER
!          QUANTITIES REQUIRED BY RCA
!**********************************************************************C
!
!qa      DIMENSION AAMAX(IM,JM,KB),AAMAY(IM,JM,KB)
!qa      EQUIVALENCE (AAMAX,VH),(AAMAY,VHP)
!
      DIMENSION 
     .  KHLPF(IM,JM,KB),  AAHLPF(IM,JM,KB),
     .   SLPF(IM,JM,KB),    TLPF(IM,JM,KB),
     .   ULPF(IM,JM,KB),    VLPF(IM,JM,KB),   WLPF(IM,JM,KB),
     .   ES(IM,JM),         ED(IM,JM)
!yang 04/29/99
      REAL  CSED1LPF(IM,JM,KBM1),   CSED2LPF(IM,JM,KBM1)
     .     ,WSET1LPF(IM,JM,KBM1),   WSET2LPF(IM,JM,KBM1)
     .     ,SDEP1LPF(IM,JM),      SDEP2LPF(IM,JM)
     .     ,SRES1LPF(IM,JM),      SRES2LPF(IM,JM) 
     .     ,TAULPF(IM,JM) 
      REAL  CBED1(IM,JM),         CBED2(IM,JM)
      REAL  VRES(IM,JM),TBEDTHIK(IM,JM)
      REAL*8 BEDTHIK(IM,JM)
!yang end
 
      REAL KMLPF(IM,JM,KB),Q2LPF(IM,JM,KB),LLPF(IM,JM,KB)

!     DECLARE ARRAYS FOR WET GRID OPTION
      DIMENSION  RCAU(MAXWET,KBM1),      RCAV(MAXWET,KBM1),
     .           RCAW(MAXWET,KBM1),      
     .           RCAAMX(MAXWET,KBM1),    RCAAMY(MAXWET,KBM1),
     .           RCAKH(MAXWET,KBM1),     RCAKM(MAXWET,KBM1),
     .           RCAS(MAXWET,KBM1),      RCAT(MAXWET,KBM1),
     .           RCAES(MAXWET),          RCAED(MAXWET),
     .           INDX(MAXWET),           JNDX(MAXWET)
!yang 04/29/99
      REAL       RCACSED1(MAXWET,KBM1),  RCACSED2(MAXWET,KBM1)
     .          ,RCAWSET1(MAXWET,KBM1),  RCAWSET2(MAXWET,KBM1)
     .          ,RCASDEP1(MAXWET),       RCASDEP2(MAXWET)
     .          ,RCACBED1(MAXWET),       RCACBED2(MAXWET)
     .          ,RCAVRES(MAXWET),        RCATAU(MAXWET)
     .          ,RCATHIK(MAXWET)
!yang end

      REAL KHLPF
      DATA NHWQI/1/
!
! **  DEFINITION OF VARIABLES
!     FLTWT          FILTER WEIGHT=1/(NPLPF)
!     ES(I,J)        INITIAL ELEVATION
!     ED(I,J)        VOLUME RATE OF CHANGE
!     NHWQI          TIME STEP COUNTER CONTROLING QPERATIONS IN HWQI
!     NPLPF          NUMBER OF TIME STEPS PER FILTER CYCLE
!     KHLPF(I,J,K)   LOW PASS FILTER OF VERTICAL MASS DIFFUSIVITY
!     AAHLPF(I,J,K)  LOW PASS FILTER OF HORIZONTAL MASS DIFFUSIVITY
!     SLPF(I,J,K)    LOW PASS FILTERED SALINITY
!     TLPF(I,J,K)    LOW PASS FILTERED TEMPERATURE
!     ULPF(I,J,K)    LOW PASS FILTER OF U*(H+ET)
!     VLPF(I,J,K)    LOW PASS FILTER OF VP*(H+ET)
!     WLPF(I,J,K)    LOW PASS FILTER OF W

!yang 04/29/99
!     CSED1LPF(I,J,K)  LOW PASS FILTER OF CSED1 (mg/L), SOLID1 CONC in WC
!     CSED2LPF(I,J,K)  LOW PASS FILTER OF CSED2 (mg/L), SOLID2 CONC in WC
!     WSET1LPF(I,J,K)  LOW PASS FILTER OF WSET1 (m/day), SOLID SETTLING RATE
!     WSET2LPF(I,J,K)  LOW PASS FILTER OF WSET2 (m/day), SOLID SETTLING RATE
!     SDEP1LPF(I,J)    LOW PASS FILTER OF SOLID1 DEPOSITION RATE (m/day) 
!     SDEP2LPF(I,J)    LOW PASS FILTER OF SOLID2 DEPOSITION RATE (m/day)
!     SRES1LPF(I,J)    LOW PASS FILTER OF RESUSPENTION FLUX OF E1 (gm/m^2-day)
!     SRES2LPF(I,J)    LOW PASS FILTER OF RESUSPENTION FLUX OF E2 (gm/m^2-day)
!     TAULPF(I,J)      LOW PASS FILTER OF BOTTOM SHEAR STRESS (dynes/cm^2)
!     CBED1(I,J)       SOLID1 CONC (g/cm^3) in SEDIMENT
!     CBED2(I,J)       SOLID2 CONC (g/cm^3) in SEDIMENT
!     VRES(I,J)        RESUSPENTION RATE (cm/day) OF TOTAL SOLID 
!     BEDTHIK(I,J)     BED THICKNESS (cm) of COHESIVE SEDIMENT
!yang end
!
!
!     KMLPF(I,J,K)    LOW PASS FILTER OF KM
!     Q2LPF(I,J,K)    LOW PASS FILTER OF Q2
!     LLPF(I,J,K)     LOW PASS FILTER OF L
!     
!**********************************************************************C
!
! **  VARIABLES ARE INITIALIZED
!
      IF (NHWQI.GT.1) GO TO 110
!
!**********************************************************************C  
!
! **  INITIALIZE TIME AVERAGING VARIABLES 
!
!**********************************************************************C
!
      IF(TOR.EQ.'BAROTROPIC') THEN
        DO 120 J=1,JM
          DO 120 I=1,IM
            ES(I,J)=(ELB(I,J)+EL(I,J))/2.
            ED(I,J)=0.0
 120    CONTINUE
!
      ELSE
!
        DO 125 J=1,JM
          DO 125 I=1,IM
            ES(I,J)=(ETB(I,J)+ET(I,J))/2.
            ED(I,J)=0.0

!yang 04/30/99
            SDEP1LPF(I,J)=0.
            SDEP2LPF(I,J)=0.
            SRES1LPF(I,J)=0.
            SRES2LPF(I,J)=0.
            TAULPF(I,J)=0.
            CBED1(I,J)=0.
            CBED2(I,J)=0.
            VRES(I,J)=0.
            BEDTHIK(I,J)=0.0
!yang end

 125    CONTINUE
      END IF
!

      DO 135 K=1,KB
        DO 135 J=1,JM
          DO 135 I=1,IM
            KHLPF(I,J,K)=0.
            KMLPF(I,J,K)=0.
            Q2LPF(I,J,K)=0.
            LLPF(I,J,K)=0.
            AAHLPF(I,J,K)=0.
            SLPF(I,J,K)=0.
            TLPF(I,J,K)=0.
            ULPF(I,J,K)=0.
            VLPF(I,J,K)=0.
            WLPF(I,J,K)=0.
 135  CONTINUE 

!yang 04/29/99
      DO 137 K=1,KBM1
        DO 137 J=1,JM
          DO 137 I=1,IM
            CSED1LPF(I,J,K)=0.
            CSED2LPF(I,J,K)=0.
            WSET1LPF(I,J,K)=0.
            WSET2LPF(I,J,K)=0.
 137  CONTINUE 
!yang end
!**********************************************************************C
! **  ACCUMULATE PIECES
!
  110 CONTINUE
!**********************************************************************C
      IF(TOR.EQ.'BAROTROPIC') THEN
        DO 140 J=1,JM
          DO 140 I=1,IM
            ED(I,J)=ED(I,J)+(ELF(I,J)-ELB(I,J))/(2.*DTE)
 140    CONTINUE
!
        DO 150 J=2,JM
          DO 150 I=2,IM
            AAHLPF(I,J,1)=AAHLPF(I,J,1)+AAM2D(I,J)
            ULPF(I,J,1)=ULPF(I,J,1)+UA(I,J)*0.5*(D(I,J)+D(I-1,J))
            VLPF(I,J,1)=VLPF(I,J,1)+VA(I,J)*0.5*(D(I,J)+D(I,J-1))
 150    CONTINUE
!
      ELSE
!
        DO 145 J=1,JM
          DO 145 I=1,IM
            ED(I,J)=ED(I,J)+(ETF(I,J)-ETB(I,J))/(2.*DTI)
 145    CONTINUE

!yang 05/03/99
! convert solid deposition fluxes DD [g/cm^2] to [g/m^2-day]
! convert sediment resuspension fluxes [g/cm^2] to [g/m^2-day]
        IF (SEDTRAN.EQ.'INCLUDE') THEN
          DO 147 J=1,JM
          DO 147 I=1,IM
            SDEP1LPF(I,J)=SDEP1LPF(I,J)+10000.*DD(1,I,J)/DTI*86400.
            SDEP2LPF(I,J)=SDEP2LPF(I,J)+10000.*DD(2,I,J)/DTI*86400.
            SRES1LPF(I,J)=SRES1LPF(I,J)+10000.*E(1,I,J)/DTI*86400.
            SRES2LPF(I,J)=SRES2LPF(I,J)+10000.*E(2,I,J)/DTI*86400.
! Bottom Shear Stress
            TAULPF(I,J)=TAULPF(I,J)+TAU(I,J,KB)
! Bottom thickness
!         DO 521 J=1,JM
!           DO 521 I=1,IM
              IF (IBMSK(I,J).EQ.0) THEN
                DO 1521 LL=1,LAYMAX
                  BEDTHIK(I,J)=BEDTHIK(I,J)+TSED(LL,I,J)
 1521           CONTINUE
              ELSE
                IF (IBMSK(I,J).EQ.1) THEN
                  BEDTHIK(I,J)=BEDTHIK(I,J)+BEDTH(1,I,J)
                ELSE
                  BEDTHIK(I,J)=0.0
                ENDIF
              ENDIF
!521      CONTINUE

 147      CONTINUE

! WC solid concs CSEDx [g/cm^3]
! WC solid settling rates WSETx [m/sec]
          DO 148 K=1,KBM1
          DO 148 J=1,JM
          DO 148 I=1,IM
            CSED1LPF(I,J,K)=CSED1LPF(I,J,K)+CSED1(I,J,K)
            CSED2LPF(I,J,K)=CSED2LPF(I,J,K)+CSED2(I,J,K)
            WSET1LPF(I,J,K)=WSET1LPF(I,J,K)+WSET1(I,J,K)
            WSET2LPF(I,J,K)=WSET2LPF(I,J,K)+WSET2(I,J,K)
 148      CONTINUE
        ENDIF
!yang end
!
        DO 155 K=1,KBM1
          DO 155 J=2,JM
            DO 155 I=2,IM
              SLPF(I,J,K)=SLPF(I,J,K)+S(I,J,K)
              TLPF(I,J,K)=TLPF(I,J,K)+T(I,J,K)
!
              IF (VERTMIX.EQ.'CONSTANT  ') THEN
                KHLPF(I,J,K)=KHLPF(I,J,K)+(KH(I,J,K))*FSM(I,J)
                KMLPF(I,J,K)=KMLPF(I,J,K)+(KM(I,J,K))*FSM(I,J)
              ELSE
                KHLPF(I,J,K)=KHLPF(I,J,K)+(KH(I,J,K)+UMOL*VARUF(I,J))*
     .                       FSM(I,J)
                KMLPF(I,J,K)=KMLPF(I,J,K)+(KM(I,J,K)+UMOL*VARUF(I,J))*
     .                       FSM(I,J)
                Q2LPF(I,J,K)=Q2LPF(I,J,K)+Q2(I,J,K)*FSM(I,J)
                LLPF(I,J,K)=LLPF(I,J,K)+L(I,J,K)*FSM(I,J)
              ENDIF
              AAHLPF(I,J,K)=AAHLPF(I,J,K)+AAM(I,J,K)/HPRNU
              ULPF(I,J,K)=ULPF(I,J,K)+U(I,J,K)*0.5*(DT(I,J)+DT(I-1,J))
              VLPF(I,J,K)=VLPF(I,J,K)+V(I,J,K)*0.5*(DT(I,J)+DT(I,J-1))
              WLPF(I,J,K)=WLPF(I,J,K)+W(I,J,K)

 155    CONTINUE
      END IF
!**********************************************************************C
!     
! **  CHECK FOR END OF FILTER
!
      IF (NHWQI.LT.NPLPF) GO TO 190
!
! **  IF NHWQI=NPLPF, COMPLETE THE FILTERING
!
!**********************************************************************C
!
      DO 200 J=1,JM
        DO 200 I=1,IM
          ED(I,J)  =FLTWT*ED(I,J)
 200  CONTINUE
!yang 05/03/99
      IF (SEDTRAN.EQ.'INCLUDE') THEN
! convert WC solid concs CSEDx [g/cm^3] to [mg/L]
! convert WC solid settling rates [m/sec] to [m/day]
        DO 215 K=1,KBM1
        DO 215 J=1,JM
        DO 215 I=1,IM
          CSED1LPF(I,J,K) = FLTWT*CSED1LPF(I,J,K) *1000000.
          CSED2LPF(I,J,K) = FLTWT*CSED2LPF(I,J,K) *1000000.
          WSET1LPF(I,J,K) = FLTWT*WSET1LPF(I,J,K) *86400.
          WSET2LPF(I,J,K) = FLTWT*WSET2LPF(I,J,K) *86400.
 215    CONTINUE

! convert WC solid deposition fluxes [g/m^2-day] to deposition rate [m/day]

        DO 217 J=1,JM
        DO 217 I=1,IM
          IF(CSED1LPF(I,J,KBM1) .NE. 0.0)
     .     SDEP1LPF(I,J) = FLTWT*SDEP1LPF(I,J)/CSED1LPF(I,J,KBM1)
          IF(CSED2LPF(I,J,KBM1) .NE. 0.0)
     .     SDEP2LPF(I,J) = FLTWT*SDEP2LPF(I,J)/CSED2LPF(I,J,KBM1)
          SRES1LPF(I,J) = FLTWT*SRES1LPF(I,J)
          SRES2LPF(I,J) = FLTWT*SRES2LPF(I,J)
          TAULPF(I,J) = FLTWT*TAULPF(I,J)
! combine sediment resuspension fluxes [g/m^2-day] and convert to resuspension 
! rate [cm/day]
          IF(CBED(I,J) .NE. 0.0)
     .      VRES(I,J) = (SRES1LPF(I,J)+SRES2LPF(I,J))/CBED(I,J)/10000.

! Bed thickness: 

          IF (IBMSK(I,J).EQ.0) THEN
            TSETOT=0.0
            DO 640 LL=1,LAYMAX
             TSETOT=TSETOT+TSED0(LL,I,J)
 640       CONTINUE
!
           BEDTHIK(I,J)=BEDTHIK(I,J)-TSETOT
!
!  CONVERT FROM g/cm**2 TO cm
!
           BEDTHIK(I,J)=BEDTHIK(I,J)/CBED(I,J)*FLTWT
!
            IF (SEDTYPE.EQ.'SAND') BEDTHIK(I,J)=0.0
           ELSE

           IF (IBMSK(I,J).EQ.1) THEN
!
!onvert thickness from g/cm**2 to cm
!
            BEDTHIK(I,J)=100.*(BEDTHIK(I,J)/((CBED(I,J)/2.65)*
     +                     H1(I,J)*H2(I,J))-BEDTHI)*FLTWT
            IF (SEDTYPE.EQ.'MUD ') BEDTHIK(I,J)=0.0
           ELSE
             BEDTHIK(I,J)=0.0
           ENDIF
          ENDIF
          
! estimate sediment solid concs [g/cm^3] 
          IF(VRES(I,J) .NE. 0.0) THEN
            CBED1(I,J) = SRES1LPF(I,J)/VRES(I,J)/10000.
            CBED2(I,J) = SRES2LPF(I,J)/VRES(I,J)/10000.
          ELSE
            DO LL=1,LAYMAX
              IF(PSED1(LL,I,J).EQ.0.0 .AND. PSED2(LL,I,J).EQ.0.0) THEN
                GO TO 216
              ELSE
                CBED1(I,J) = CBED(I,J)*PSED1(LL,I,J)
                CBED2(I,J) = CBED(I,J)*PSED2(LL,I,J)
                GO TO 217
              ENDIF
 216          CONTINUE
              IF(LL.EQ.LAYMAX)THEN
                CBED1(I,J) = CBED(I,J)*0.8
                CBED2(I,J) = CBED(I,J)*0.2
              ENDIF
            ENDDO
          ENDIF
 217    CONTINUE

      ENDIF
!yang end
!
      DO 220 K=1,KBM1
        DO 220 J=1,JM
          DO 220 I=1,IM
            SLPF(I,J,K)  =FLTWT*SLPF(I,J,K)
            TLPF(I,J,K)  =FLTWT*TLPF(I,J,K)
            KHLPF(I,J,K) =FLTWT*KHLPF(I,J,K)
            KMLPF(I,J,K) =FLTWT*KMLPF(I,J,K)
            Q2LPF(I,J,K) =FLTWT*Q2LPF(I,J,K)
            LLPF(I,J,K)  =FLTWT*LLPF(I,J,K)
            AAHLPF(I,J,K)=FLTWT*AAHLPF(I,J,K)
            WLPF(I,J,K)  =FLTWT*WLPF(I,J,K)
            ULPF(I,J,K)  =FLTWT*ULPF(I,J,K)
            VLPF(I,J,K)  =FLTWT*VLPF(I,J,K)
 220  CONTINUE
!
!**********************************************************************C
! **  CALCULATE TIME AVERAGED FLOW FIELDS ACROSS CELL FACES
! **  QLRX IS STORED IN ULPF, QLRY IN VLPF, QLRZ IN WLPF .
!
      IF(TOR.NE.'BAROTROPIC') THEN
!----------------------------------------------------------------------C
!     CALCULATE THE VERTICAL COMPONENT
!----------------------------------------------------------------------C
!
        DO 250 J=2,JMM1
          DO 250 K=2,KBM1
            DO 250 I=2,IMM1
              WLPF(I,J,K)=WLPF(I,J,K)*ART(I,J)
 250    CONTINUE
      END IF
!
!----------------------------------------------------------------------C
!     CALCULATE THE HORIZONTAL COMPONENTS
!----------------------------------------------------------------------C
!
      DO 260 J=2,JM
        DO 260 K=1,KBM1
          DO 260 I=2,IM
            ULPF(I,J,K)=0.5*(H2(I,J)+H2(I-1,J))*ULPF(I,J,K)
            VLPF(I,J,K)=0.5*(H1(I,J)+H1(I,J-1))*VLPF(I,J,K)
!
            ULPF(I,J,K)=ULPF(I,J,K)*DZ(K)*DUM(I,J)
            VLPF(I,J,K)=VLPF(I,J,K)*DZ(K)*DVM(I,J)
!
 260  CONTINUE
!
      DO 340 J=2,JM
        DO 340 K=1,KBM1
          DO 340 I=2,IM
            AAMAX(I,J,K)=.5*(AAHLPF(I,J,K)+AAHLPF(I-1,J,K))*DUM(I,J)
            AAMAY(I,J,K)=.5*(AAHLPF(I,J,K)+AAHLPF(I,J-1,K))*DVM(I,J)
  340 CONTINUE
      DO 350 N=1,NUMQBC
        ID=IQD(N)
        JD=JQD(N)
        IC=IQC(N)
        JC=JQC(N)
        IF(JD.EQ.JC) THEN
            IF(IC.GT.ID) THEN
              DO 360 K=1,KBM1
                AAMAX(IC,JC,K)=0.0
 360          CONTINUE
            ELSE
              DO 370 K=1,KBM1
                AAMAX(ID,JD,K)=0.0
 370          CONTINUE
            ENDIF
       ELSE
            IF(JC.GT.JD) THEN
              DO 380 K=1,KBM1
                AAMAY(IC,JC,K)=0.0
 380          CONTINUE
            ELSE
              DO 390 K=1,KBM1
                AAMAY(ID,JD,K)=0.0
 390          CONTINUE
            ENDIF
         ENDIF
  350  CONTINUE
!**********************************************************************C
! **  OUTPUT RESULTS TO DISK  
!
! **  WRITE CONSTANTS FIRST TIME THROUGH 
      IF (CONSTRN) THEN
        DO 275 J=1,JM
          DO 275 I=1,IM
            TPS(I,J)=FSM(I,J)
 275    CONTINUE
        DO 285 N=1,NUMEBC
          IE=IETA(N)
          JE=JETA(N)
          TPS(IE,JE)=-2.
 285    CONTINUE
        DO 295 N=1,NUMQBC
          IC=IQC(N)
          JC=JQC(N)
          TPS(IC,JC)=-1.
 295    CONTINUE
!
        IF(IWET.EQ.1)THEN
*   OPEN AN FILE FOR REDUCED-SIZE GCM_TRAN FILE INFO
        ICNTRAN=0
        DO 3400 I=1,IM
          DO 3400 J=1,JM
            IF(TPS(I,J).EQ.0.0)GO TO 3400
            ICNTRAN=ICNTRAN+1
            INDX(ICNTRAN)=I
            JNDX(ICNTRAN)=J
 3400   CONTINUE

        OPEN (IUTRN,FORM='unformatted',FILE='wet_grid')
        WRITE(IUTRN)ICNTRAN
        DO I=1,ICNTRAN
         WRITE(IUTRN)INDX(I),JNDX(I)
        END DO
        CLOSE(IUTRN)
      END IF

       OPEN (IUTRN,FORM='unformatted',FILE='gcm_geom')
       WRITE(IUTRN) DZ,DZZ
       WRITE(IUTRN) H,H1,H2,TPS
       WRITE(IUTRN) ANG,NU
       CLOSE(IUTRN)
       OPEN (IUTRN,FORM='unformatted',FILE='gcm_tran')
!
!qa       IF (SEDTRAN.EQ.'INCLUDE') THEN
!qa         OPEN (UNIT=110,FILE='gcm_sedtran',FORM='UNFORMATTED')
!qa       ENDIF
!yang 04/29/99
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         OPEN (UNIT=121,FILE='gcm_sedtran',FORM='UNFORMATTED')
       ENDIF
!yang end
!
       CONSTRN=.FALSE.
      ENDIF
! **
      TMIDDLE=TIME-(.5*DTI*DAYI/FLTWT)
      WRITE(IUTRN) TMIDDLE
      IF(IWET.EQ.0)THEN
       WRITE(IUTRN) (((ULPF(I,J,K),I=1,IM),J=1,JM),K=1,KBM1)
       WRITE(IUTRN) (((VLPF(I,J,K),I=1,IM),J=1,JM),K=1,KBM1)
       WRITE(IUTRN) WLPF
       WRITE(IUTRN) (((AAMAX(I,J,K),I=1,IM),J=1,JM),K=1,KBM1)
       WRITE(IUTRN) (((AAMAY(I,J,K),I=1,IM),J=1,JM),K=1,KBM1)
       WRITE(IUTRN) KHLPF
!jjs       WRITE(IUTRN) KMLPF
       WRITE(IUTRN) ES
       WRITE(IUTRN) ED
       WRITE(IUTRN) (((SLPF(I,J,K),I=1,IM),J=1,JM),K=1,KBM1)
       WRITE(IUTRN) (((TLPF(I,J,K),I=1,IM),J=1,JM),K=1,KBM1)
!
!qa       IF (SEDTRAN.EQ.'INCLUDE') THEN
!qa         WRITE (110)KMLPF
!qa         WRITE (110)Q2LPF
!qa         WRITE (110)LLPF
!qa       ENDIF
!yang 04/29/99
          IF (SEDTRAN.EQ.'INCLUDE') THEN
!
!     TRANSFER FROM BEDTHIK(real*8) TO TBEDTHIK (real*4)
           DO J=1,JM
            DO I=1,IM
             TBEDTHIK(I,J)=BEDTHIK(I,J)
            END DO
           END DO
            WRITE(121) TMIDDLE
            WRITE (121)CSED1LPF
            WRITE (121)CSED2LPF
            WRITE (121)WSET1LPF
            WRITE (121)WSET2LPF
            WRITE (121)SDEP1LPF
            WRITE (121)SDEP2LPF
            WRITE (121)VRES
            WRITE (121)TBEDTHIK
            WRITE (121)CBED1
            WRITE (121)CBED2
            WRITE (121)TAULPF
          ENDIF
!yang end
      ELSE
! **  EXTRACT TRANSPORT INFO FOR WET GRID ONLY
      DO 4401 K=1,KBM1
      DO 4401 I=1,ICNTRAN
       RCAU(I,K)=ULPF(INDX(I),JNDX(I),K)
       RCAV(I,K)=VLPF(INDX(I),JNDX(I),K)
       RCAW(I,K)=WLPF(INDX(I),JNDX(I),K)
       RCAKH(I,K)=KHLPF(INDX(I),JNDX(I),K)
       RCAKM(I,K)=KMLPF(INDX(I),JNDX(I),K)
       RCAAMX(I,K)=AAMAX(INDX(I),JNDX(I),K)
       RCAAMY(I,K)=AAMAY(INDX(I),JNDX(I),K)
       RCAS(I,K)=SLPF(INDX(I),JNDX(I),K)
       RCAT(I,K)=TLPF(INDX(I),JNDX(I),K)
!yang 04/29/99
       RCACSED1(I,K)=CSED1LPF(INDX(I),JNDX(I),K)
       RCACSED2(I,K)=CSED2LPF(INDX(I),JNDX(I),K)
       RCAWSET1(I,K)=WSET1LPF(INDX(I),JNDX(I),K)
       RCAWSET2(I,K)=WSET2LPF(INDX(I),JNDX(I),K)
!yang end
 4401  CONTINUE
      DO 4402 I=1,ICNTRAN
       RCAES(I)=ES(INDX(I),JNDX(I))
       RCAED(I)=ED(INDX(I),JNDX(I))
!yang 04/29/99
       RCASDEP1(I)=SDEP1LPF(INDX(I),JNDX(I))
       RCASDEP2(I)=SDEP2LPF(INDX(I),JNDX(I))
       RCAVRES(I)=VRES(INDX(I),JNDX(I))
       RCATHIK(I)=BEDTHIK(INDX(I),JNDX(I))
       RCACBED1(I)=CBED1(INDX(I),JNDX(I))
       RCACBED2(I)=CBED2(INDX(I),JNDX(I))
       RCATAU(I)=TAULPF(INDX(I),JNDX(I))
!yang end
 4402 CONTINUE

       WRITE(IUTRN)((RCAU(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAV(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAW(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAAMX(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAAMY(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAKH(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAKM(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)(RCAES(I),I=1,ICNTRAN)
       WRITE(IUTRN)(RCAED(I),I=1,ICNTRAN)
       WRITE(IUTRN)((RCAS(I,K),I=1,ICNTRAN),K=1,KBM1)
       WRITE(IUTRN)((RCAT(I,K),I=1,ICNTRAN),K=1,KBM1)
!yang 04/29/99
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         WRITE(121) TMIDDLE
         WRITE (121)((RCACSED1(I,K),I=1,ICNTRAN),K=1,KBM1)
         WRITE (121)((RCACSED2(I,K),I=1,ICNTRAN),K=1,KBM1)
         WRITE (121)((RCAWSET1(I,K),I=1,ICNTRAN),K=1,KBM1)
         WRITE (121)((RCAWSET2(I,K),I=1,ICNTRAN),K=1,KBM1)
         WRITE (121)(RCASDEP1(I),I=1,ICNTRAN)
         WRITE (121)(RCASDEP2(I),I=1,ICNTRAN)
         WRITE (121)(RCAVRES(I),I=1,ICNTRAN)
         WRITE (121)(RCATHIK(I),I=1,ICNTRAN)
         WRITE (121)(RCACBED1(I),I=1,ICNTRAN)
         WRITE (121)(RCACBED2(I),I=1,ICNTRAN)
         WRITE (121)(RCATAU(I),I=1,ICNTRAN)
       ENDIF
!yang end
      END IF
!
!**********************************************************************C
!     RESET HWQI CONTROL COUNTER
!**********************************************************************C
      NHWQI=0
!**********************************************************************C
  190 CONTINUE
!**********************************************************************C
!     INCREMENT HWQI CONTROL COUNTER
!----------------------------------------------------------------------C
      NHWQI=NHWQI+1
!**********************************************************************C 
      RETURN
      END
