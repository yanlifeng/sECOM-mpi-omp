      SUBROUTINE SOLAR(JD,FJD,R,DLT,ALP,SLAG,N,HANG,TAUDA,COSZ,SOLALT,
     1 PHIT)
!
!    *******************************************************************
!    *                                                                 *
!    *                            S O L A R                            *
!    *                                                                 *
!    *******************************************************************
!
!.....SOLAR COMPUTES RADIUS VECTOR, DECLINATION AND RIGHT ASCENSION OF
!.....SUN, EQUATION OF TIME, AND HOUR ANGLE OF SUN AT SUNSET FOR N
!.....EQUALLY SPACED LATITUDES GIVEN JULIAN DAY AND FRACTION.
!
      PARAMETER (JMT=1)
      DIMENSION ALAT(JMT)
!
                              D I M E N S I O N
     1   HANG( JMT ), COSZ( JMT ),  TAUDA( JMT ), SOLALT (JMT)
!
                                   D A T A
     1   TPI/6.283185/,      HPI/1.570796/,
     2   RAD/57.29578/,      CYEAR/365.25/,      CCR/1.E-6/
!
!.....TPP = DAYS BETWEEN EPOCH AND PERIHELION PASSAGE OF 1900
!.....SVT6 = DAYS BETWEEN PERIHELION PASSAGE AND MARCH EQUINOX OF 1900
!.....JDOR = JD OF EPOCH WHICH IS JANUARY 0, 1900 AT 12 HOURS UT
!
                                   D A T A
     1   TPP/1.55/,          SVT6/78.035/,       JDOR/2415020/
      DATA ITERS/100/
!
!    *******************************************************************
!
!
      PI=3.141593
      ALAT(1)=PHIT
      DAT=FLOAT(JD-JDOR)-TPP+FJD
!
!    COMPUTES TIME IN JULIAN CENTURIES AFTER EPOCH
!
      T=FLOAT(JD-JDOR)/36525.0
!
!    COMPUTES LENGTH OF ANOMALISTIC AND TROPICAL YEARS (MINUS 365 DAYS)
!
      YEAR=0.25964134+0.304E-5*T
      TYEAR=0.24219879-0.614E-5*T
!
!    COMPUTES ORBIT ECCENTRICITY AND ANGLE OF EARTH'S INCLINATION FROM T
!
      EC=0.01675104-(0.418E-4+0.126E-6*T)*T
      ANGIN=23.452294-(0.0130125+0.164E-5*T)*T
      ADOR=JDOR
      JDOE=ADOR+(SVT6*CYEAR)/(YEAR-TYEAR)
!
!    DELEQN=UPDATED SVT6 FOR CURRENT DATE
!
!
!
      DELEQN=FLOAT(JDOE-JD)*(YEAR-TYEAR)/CYEAR
!
      YEAR=YEAR+365.0
!
      SNI=SIN(ANGIN/RAD)
      TINI=1.0/TAN(ANGIN/RAD)
      ER=SQRT((1.0+EC)/(1.0-EC))
      QQ=DELEQN*TPI/YEAR
!
!    DETERMINE TRUE ANOMALY AT EQUINOX
!
      E=1.0
!
      DO 200 ITER=1,ITERS
      EP=E-(E-EC*SIN(E)-QQ)/(1.0-EC*COS(E))
      CD=ABS(E-EP)
      E=EP
      IF(CD.LE.CCR) GO TO 250
 200  CONTINUE
      PRINT 400
 250  CONTINUE
      HE=0.5*E
      EQ=2.0*ATAN(ER*TAN(HE))
!
!    DATE=DAYS SINCE LAST PERIHELION PASSAGE
!
      DATE=MOD(DAT,YEAR)
!
!    SOLVE ORBIT EQUATIONS BY NEWTON'S METHOD
!
      EM=TPI*DATE/YEAR
      E=1.0
!
      DO 300 ITER=1,ITERS
      EP=E-(E-EC*SIN(E)-EM)/(1.0-EC*COS(E))
      CR=ABS(E-EP)
      E=EP
      IF(CR.LE.CCR) GO TO 350
 300  CONTINUE
      PRINT 400
 350  CONTINUE
!
      R=1.0-EC*COS(E)
      HE=0.5*E
      W=2.0*ATAN(ER*TAN(HE))
      SIND=SNI*SIN(W-EQ)
      DLT=ASIN(SIND)
      ALP=ASIN(TAN(DLT)*TINI)
      TST=COS(W-EQ)
!
      IF(TST.LT.0.0) ALP=PI-ALP
!
      IF(ALP.LT.0.0) ALP=ALP+TPI
!
      SUN=TPI*(DATE-DELEQN)/YEAR
      IF (SUN.LT.0.0)  SUN=SUN+TPI
      SLAG=SUN-ALP-0.03255
!
!    COMPUTE HOUR ANGLE OF SUNSET AT ALL LATITUDES
!
!
      DO 1 I=1,N
      PH=ALAT(I)
      SS=SIN(PH)*SIN(DLT)
      CC=COS(PH)*COS(DLT)
!
      SOLALT(I)= SS + CC
      SOLALT(I) = ASIN(SOLALT(I)) * RAD
!
      IF(DLT.EQ.0.0) GO TO 16
!
      AP=ABS(PH)
      EPS=ABS(AP-HPI)
      IF(EPS.GT.CCR) GO TO 14
!
      H=HPI*ABS(AP/PH+ABS(DLT)/DLT)
      GO TO 5
!
   14 AR=-SS/CC
      AC=ABS(AR)
      IF(AC-1.0+CCR) 4,2,3
!
  2   H=(AC-AR)*HPI
!
      GO TO 5
!
  3   IF(AR.LT.0.0) GO TO 25
!
      H=0.0
      GO TO 5
!
 25   H=PI
      GO TO 5
!
!
 16   H=HPI
      GO TO 5
!
  4   H=ACOS(AR)
  5   HANG(I)=H
      TAUDA(I)=MAX(H/PI,0.0)
!
      IF(H.EQ.0.0) GO TO 100
!
!       write(*,*)'H in Solar=',H
      COSZ(I)=MAX((SS+CC*SIN(H)/H),0.0)
      GO TO 1
!
!
!
  100 COSZ(I)=0.0
    1 CONTINUE
!
      RETURN
 400  FORMAT('0 MAXIMUM ITERATIONS IN SUBROUTINE SOLAR')
      END
