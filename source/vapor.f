      SUBROUTINE VAPOR(EA,ES,RH,SST,BP,SAL)
!
!MNT  RH IN PERCENT, SST IN DEG C, BP IN MB
!MNT  COMPUTE VAPOR PRESSURE OF WATER FROM RELATIVE
!MNT  HUMIDITY
!MNT
!MNT  R.WELLER, 6/18/87
!MNT
!MNT  COMPUTE SATURATION VAPOR PRESSURE OF PURE WATER
!MNT  USING TABATA FORMULA
!**********************************************************
!       Modified November 6, 1990       M. Park Samelson
!       
!       Previously used Tabata formula (1973).  This
!       version uses Buck formula (1981).  Temperature
!       is in degrees C instead of degrees K.
!
!**********************************************************

!*******Tabata formula**************
!**** TK=273.16
!**** SSTK=SST+TK
!**** A=1000./SSTK
!**** E=8.42926609-(1.82717483*A)-(.071208271*A*A)
!**** ES=10.**E

!********Buck formula***************
!NG It calculates water vapor (partial) pressure of saturated moist air
!NG Original Buck (1981) used up to 2011, below:
!NG            ES = (1.0007+3.46E-6*BP)*6.1121*EXP(17.502 * SST/ (240.97 + SST))
!MNT  CORRECT FOR SALT WATER AT 35 PPT (CNG10252011 Why HARDWIRE at 35ppt?????)
!NG            EW=.9815*ES
!NGforIce      Note: this would be the following with respect to plane surface of ice:
!NGforIceUse:  ES = (1.0003+4.18E-6*BP)*6.1115*EXP(22.452 * SST/ (272.55 + SST))

!NG10252011 New Buck (1996) adopted below:
      ES = (1.0007+3.46E-6*BP)*6.1121*
     .     EXP( (18.677 - SST/234.50) * (SST/(257.14 + SST)) )
!NG10252011 Corrected for actual salinity based on, for example,
!NG10252011 List, R.J. (1984) Smithsonian Meteorological Tables. Sixth ed. Washington, 
!NG10252011 D.C.: Smithsonian Institution Press, 527 p
      EW = ES*(1 - 0.000537*SAL)
!NG10262011IfForIce      ES = (1.0003+4.18E-6*BP)*6.1121*
!NG10262011IfForIce     .     EXP( (23.036 - SST/333.70) * (SST/(279.82 + SST)) )
!NG10262011IfForIce      EW=ES ! No salinity correction (ice rejects salt)

!MNT  COMPUTE SATURATION MIXING RATIO
      B=EW/BP
      RW=(B*.62197)/(1.-B)
!MNT  COMPUTE MIXING RATIO
      R=(RH/100.)*RW
!MNT  COMPUTE VAPOR PRESSURE
      EA=(R/(.62197+R))*BP
      RETURN
      END
