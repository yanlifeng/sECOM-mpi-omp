      FUNCTION QSAT(T,P)
!
!        T         temperature in degrees C
!        P         pressure in millibars
!        ES        saturation vapor pressure in millibars
!        QSAT      saturation humidity in g/kg 
!
      RD = 287.05
      RV = 461.5

      EPS = RD/RV

!*********Saturation vapor pressure**********
      ES = (1.0007+3.46E-6*P) * 6.1121 * EXP(17.502 * T/ (240.97 + T))

!*********Convert to saturation humidity*****
      QSAT = 1000. * EPS * ES / ( P + ES * (EPS - 1.0))

      RETURN
      END          
