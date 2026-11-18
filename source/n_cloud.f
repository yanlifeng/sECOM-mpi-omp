      Subroutine N_CLOUD(CUTOFF,ALAT,TIME,SWOBS,RAD,ATC,CLD)

!       A. Plueddemann 6-12-92
!       Modified from R. Weller's original routine CLOUD c. 6-18-87
!
!       This routine accepts:
!         ALAT = latitude of  observations, decimal degrees.
!         TIME = time in decimal yearday (example, day 1.5 is local
!              noon on the 2nd day).
!         SWOBS = observed shortwave radiation, it is assumed that 
!              SWOBS has been corrected for nighttime bias.
!         ATC = atmospheric transmission coeff., it is expected that
!              an appropriate ATC for the site has been empirically
!              derived.
!       and returns:
!         RAD = theoretical estimate of total shortwave radiation
!              reaching the earths surface
!         CLD = cloud cover fraction (0.0 to 1.0) for use in
!              correcting longwave radiation, nighttime values of
!              CLD are set to a flag value assuming that they will
!              later be interpolated.
!
!       routine to abandon attempts to compute cloud cover just prior
!       to sunset and just after dawn to avoid erroneous estmiates 
!       caused by persistent dawn/dusk differences in theoretical vs.
!       observed short wave.  The value of this cutoff must be "tuned"
!       for particular data sets by trial and error.

!       --------------------------------------------------------------

!       solar constant
      SC = 1395
!       conversion of degrees to radians
      DTOR = 6.283185 / 360.

!       compute altitude of the sun (angular elevation above the
!       horizon) from Smithsonian Met. Tables formulation.
!         SINA = sine of suns altitude
!         THTA = declination of sun 
!         HR = hour angle of the sun

      YRPH = (6.283185/365.) * (357.-TIME)
!       yrph uses 357 instead of 355.90826
      THTA = (-23.45*COS(YRPH)) * DTOR
!       thta uses 23.45 instead of 23.26658
      HR = (TIME+0.5) * 6.283185
      SINA = COS(THTA) * COS(ALAT*DTOR) * COS(HR) + SIN(THTA) * 
     *    SIN(ALAT*DTOR)

!       if low sun angle go to nighttime loop 
      If (SINA.GT.CUTOFF) Then
!       conditional cut off before R calculation (after R in old version)

!       estimate cloud cover using scheme described in Smithsonian
!       Meteorological Tables (R.J. List, 1984).  
!         R = radius vector of the earth, equal to the distance from
!              the center of the earth to the center of the sun
!         RSS = direct E solar radiation at the earths surface.
!         RDIFF = diffuse sky radiation.
!         RAD = total clear sky radiation at earths surface, can be
!              compared to SWOBS to determine ATC.
!         CLD = cloud cover fraction (0 to 1) used for correction to
!              longwave radiation (Fung et al, 1984).

        TRAD = 6.283185 * ((TIME-186.)/365.)
!       trad uses 186 instead of 186.72252
        R = 1.0002 + 0.01671 * COS(TRAD)
!       r uses 1.0002 instead of 1.0001531
        RSS = (SC/(R**2.)) * SINA * (ATC**(1./SINA))
        RDIFF = .5 * (SC/(R**2.)) * SINA * (0.91-(ATC**(1./SINA)))
        RAD = RSS + RDIFF
        CLD = (1.-(SWOBS/RAD)) * (10./7.)
!
!       dont allow cloud fraction greater than 1 or less than 0
        If (CLD.GT.1.) CLD = 1.
        If (CLD.LE.0.) CLD = 0.
      Else

!       nighttime loop, set theoretical radiation 
!       to zero.
        RAD = 0.
! redundant, but included to emphasize that the CLD passed INTO the 
! subroutine is passed OUT of the subroutine, i.e. cloud cover is
! PERSISTED at night.
        CLD = CLD
      End If

      Return
      End
