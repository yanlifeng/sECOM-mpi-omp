      FUNCTION SADENS(BAROP,AIRTMP,RELHUM) ! NG12022009 
!NG
!NG  ESTIMATE SURFACE AIR DENSITY (kg/m3) BASED ON:
!NG
!NG  BAROP: Barometric Pressure, mbar
!NG  AIRTMP: Surface air temperature, C
!NG  RELHUM: Relative Humidity, %
!NG

      DATA RMOLRAT /0.62197/ ! Molecular weight ratio
      EW=6.1121*(1.0007+3.46E-6*BAROP) *
     .   EXP( (18.677 - AIRTMP/234.50) * (AIRTMP/(257.14 + AIRTMP)) ) !  saturation vapor pressure, new Buck 1996 formula
      QSAT=RMOLRAT*EW/(BAROP-(1.-RMOLRAT)*EW) ! specific humidity at saturation

      Q=(0.01*RELHUM)*QSAT ! vapor pressure (specific humidity of air)
      TT=AIRTMP+273.16
      O61=1./RMOLRAT-1. ! =0.61, moisture correction for temperature
      TV=TT*(1+O61*Q) ! virtual temperature
      SADENS=100*BAROP/(287.04*TV+1.e-30) ! 287.04 is the gas constant for dry air
      SADENS=AMAX1(AMIN1(SADENS,1.6),0.9) ! Make sure they are within "reasonable" range at earth's sea level
      RETURN
      END

