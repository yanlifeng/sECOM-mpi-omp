      SUBROUTINE BULK(UZ,ZU,TZ,ZT,RH,ZQ,TSFC,TOL,UW,WT,WQ,CD,BP,NERR)
!MNT  BULK ESTIMATION OF FLUXES.  FROM BILL LARGE, NCAR, OCT 84
!MNT  ESTIMATES:
!MNT   1. TOL, THE STABILITY PARAMETER Z/L AT Z=10 METERS
!MNT   2.  UW, THE KINEMATIC MOMENTUM FLUX (M/S)**2
!MNT   3.  WT, THE KINEMATIC SENSIBLE HEAT FLUX (DEGREES M/S)
!MNT   4.  WQ, THE KINEMATIC LATENT HEAT (MOISTURE) FLUX (M/S G/M**3)
!MNT   5.  CD, THE DRAG COEFFICIENT FOR THE GIVEN STABILITY, HEIGHT
!MNT
!MNT  FROM THESE INPUTS:
!MNT   1.  UZ, THE MEAN WIND SPEED AT HEIGHT ZU
!MNT   2.  TZ, THE MEAN AIR TEMPERATURE (CELSIUS) AT HEIGHT ZT
!MNT   3.  RH, THE MEAN RELATIVE HUMIDITY (PERCENT) AT HEIGHT ZQ
!MNT   4.  TSFC, THE SEA SURFACE TEMPERATURE (CELSIUS)
!MNT   5.  BP, the barometric pressure in millibars
!MNT  
!MNT  OPTIONS:
!MNT   1. IF SET ZQ=0., REL HUM ASSUMED TO BE 75 PERCENT
!MNT   2. IF SET ZT=0., UW IS FOUND ASSUMING NEUTRAL STABILITY
!MNT   3. IF SET ZU=0., 10 M STRESS IS RETURNED
!MNT
!MNT  ASSUME:
!MNT   1. A NEUTRAL 10 M DRAG COEF CDN=0.00115 U10<UCH=11M/S
!MNT                  =0.00049 + 0.000065 U10  U10>UCH=11M/S
!MNT   2. A NEUTRAL 10 M STANTON NUMBER CTN=0.00115 Z/L<0
!MNT                                       =0.00075 Z/L>0
!MNT   3. A NEUTRAL 10 M DALTON NUMBER  CEN=0.00115 ALWAYS
!MNT   4. VON KARMAN'S CONSTANT, VONK=0.040
!MNT   5. THE SMALL WT REGRESSION, WTR BELOW IF BIG=FALSE
!MNT      WITH (A=.0008 IF Z/L>0, A=.0010 IF Z/L<0)
!MNT
!**********************************************************
!        Modified November 6, 1990        M. Park Samelson
!        
!        Replace function QSAT for computing saturation
!        humidity with Buck formula (1981) for computing
!        saturation vapor pressure.  Use specific humidity
!        formula to convert saturation vapor pressure to
!        saturation humidity.  Multiply result by 1000 to 
!        convert specific humidity from KG/KG to G/KG.
!        REQUIRES BAROMETRIC PRESSURE TO BE PASSED IN !!!
!
!        Modified December 16, 1992        M. Park Samelson
!        
!        Routine returns error parameter (NERR) if iterations
!        fail to converge, rather than writing to output, and
!        exiting from subroutine.  Allows this routine to be
!        called from C routines, since we don't want to handle
!        FORTRAN I/O (all I/O should be handled in main calling
!        routine anyway).
!**********************************************************
!
!********Start program*************************************

      LOGICAL BIG
      WTR(A,U,DTH)=0.002+A*U*DTH
!MNT
!MNT
      DATA CDA,CTS,CTU,CEA,UCH/0.00115,0.00075,0.00115,0.00115,10./
      DATA RHDEFLT,VONK/75.,0.40/
!
!        added MPS - see above note.  return error value for subroutine
      nerr = 0
!MNT                         ARE TEMPS RELIABLE??
      IF (ZT .GT. 0.2) GO TO 1
      TOL=0.0
      GO TO 3
!MNT  YES, CALCULATE A BULK STABILITY PARAMETER AT 10 METERS FROM
!MNT    DELTH,SURFACE - AIR POTENTIAL TEMPERATURE DIFFERENCE
!MNT    QSFC, SURFACE ABSOLUTE HUMIDITY
!MNT    TO, LOCAL MEAN VIRTUAL TEMERATURE
1     DELTH=TSFC-TZ-0.01*ZT
!MNT  THIS IS THE SEA MINUS AIR POT. TEMP. DIFF
      IF (ZQ .LE. 0.0) RH=RHDEFLT
      QZ=RH*QSAT(TZ,BP)*.01

!        Correct for salt water at 35ppt
      SALT = 0.9815
      DELQ=SALT*QSAT(TSFC,BP)-QZ
      TKELV=273.16 + TZ
      F1=0.00000172
      T0=TKELV+QZ*F1*TKELV**2
!MNT  
!MNT  IS /WT/ LARGE ENOUGH TO USE TRANSFER COEFFICIENTS?
      UDT=UZ*DELTH
      IF ((UDT .LT. 10.) .AND. (UDT. GT. -15.)) GO TO 2
!MNT  YES
      BIG=.TRUE.
      F2=1.00
      IF (DELTH .LT. 0.0) F2=0.70
      TOL=0.0-1000.*F2/T0/UZ**2*(DELTH+F1/F2*DELQ*T0**2)
      GO TO 3
!MNT  
!MNT  NO,     DEAL WITH SMALL SENSIBLE HEAT FLUXES
2     BIG=.FALSE.
      A=0.0010
      IF (DELTH .LT. 0.0) A=0.0008
      WT=WTR(A,UZ,DELTH)
      WQ=CEA*UZ*DELQ
      U3=(SQRT(CDA)*UZ)**3
      TOL=0.0-39./U3/T0*(WT+F1*WQ*T0**2)
!MNT
!MNT     RETURN WITH ONLY TOL IF ZU NOT GREATER THAN 0.2
3     IF (ZU .LE. 0.2) RETURN
      ZUOL=TOL*ZU/10.0
      ZTOL=TOL*ZT/10.0
      ZQOL=TOL*ZQ/10.0
!MNT     EVALUATE THE STABILITY FUNCTIONS IN FZOL
      TTOL = TOL
      CALL FZOL(TTOL,P10,PG10,S10,SG10)
      TOL = TTOL
      CALL FZOL(ZUOL,PHIM,PHI,PSIM,PSI)
      CALL FZOL(ZTOL,PHI,PHIT,PSIMT,PSIT)
      CALL FZOL(ZQOL,PHI,PHIQ,PSIMQ,PSIQ)
!MNT
!MNT     SOLVE FOR U10 AND CDN ITERATIVELY
      RK=(ALOG(ZU/10.0)-PSIM+S10)/VONK
      ITER=1
      UP=UZ/(1.0+SQRT(CDA)*RK)
22    CP=CDA
      IF (UP .GT. UCH) CP=(0.49+0.065*UP)/1000.
      U10=UZ/(1.0+SQRT(CP)*RK)
      DIFF=(UP-U10)/UP
      UP=U10
      IF (DIFF .LE. 0.01) GO TO 4
      ITER=ITER+1
      IF (ITER .LE. 5) GO TO 22
!**********************************************************
!        commented out following 5 lines, added return 
!        MPS - see explanation at beginning of program
!      WRITE(6,66)ITER,U10,DIFF,CP
! 66    FORMAT('-',5X,'TOO BAD, FAILED TO CONVERGE IN 5 ITERATIONS',I6,
!     13F16.4)
!      STOP
        nerr = iter
        return
!**********************************************************
!MNT
!MNT      SOLVED FOR U10, NOW FIND CDN AND CD
!MNT      THEN GET UW FROM BULK FORMULA
4     CDN=CDA
      IF (U10 .GT. UCH) CDN=(0.49+0.065*U10)/1000.
      SQCD=SQRT(CDN)
      RCD=1.0+SQCD/VONK*(ALOG(ZU/10.0)-PSIM)
      CD=CDN/RCD**2
      UW=0.0-CD*UZ**2
!MNT
!MNT  SENSIBLE HEAT FLUX IF ZT GREATER THAN 0.02
      IF (ZT .LE. 0.02) GO TO 5
      CTN=CTU
      IF (DELTH .LT. 0.0) CTN=CTS
      RCDT=1.0+SQCD/VONK*(ALOG(ZT/10.0)-PSIMT)
      CT=CTN/RCDT/(1.0+CTN/VONK/SQCD*(ALOG(ZT/10.)-PSIT))
      WT=CT*UZ*DELTH
      IF (BIG) GO TO 5
!MNT       CALCULATE SMALL SENSIBLE HEAT FLUX
      RKT=(ALOG(ZT/10.0)-PSIT+SG10)/VONK
      DTH10=DELTH-CT/SQRT(CD)*RKT*(DELTH)
      WT=WTR(A,U10,DTH10)
!MNT
!MNT      LATENT HEAT (MOISTURE) FLUX IF ZQ SET > 0.02
5     IF (ZQ .GT. 0.02) GO TO 6
!MNT  MOD. ORIGINAL RETURNED IF ZQ .LE. 0.02 WITHOUT
!MNT  CALCULATING A NEW WQ
      WQ=CEA*UZ*DELQ
      RETURN
6     CEN=CEA
      RCDQ=1.0+SQCD/VONK*(ALOG(ZQ/10.0)-PSIMQ)
      CE=CEN/RCDQ/(1.0+CEN/VONK/SQCD*(ALOG(ZQ/10.)-PSIQ))
      WQ=CE*UZ*DELQ
      RETURN
      END
