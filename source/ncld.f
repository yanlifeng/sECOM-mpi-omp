      Subroutine NCLD(YDAY,ALAT,ALON,SWOBS,ATC,CLD,CUTOFF,RAD)
!
! Subroutine to calculate percent cloud cover based on observed
! shortwave radiation and day of the year.  Converted from the
! original program to a subroutine. -hlj-        June 21, 1993.
!
! ORIGINAL COMMENTS ------------------------------------------------
! n_cldfrac.f       A.Oien/A. Plueddeman May 19, 1992
!
! This program estimates the daytime cloud fraction based on the
! ratio of observed to theoretical shortwave radiation.
! No estimate can be made at night, and the cloud fraction is set
! to a flag value.  It is assumed that the nighttime values will
! later be interpolated.  The program is similar to n_clearsky
! except that here the timebase and observed shortwave (with
! nighttime bias removed) are read in from the met data file.  
! Output data is handled more nicely, written to user-defined 
! output data file.
! Values for ATC and CUTOFF are determined empirically from
! iterations of n_clearsky and this program.  See n_clearsky and
! n_cloud for further description of ATC and CUTOFF.
! ORIGINAL COMMENTS ------------------------------------------------
!
! Since this subroutine will be run concurrently with the numerical
! circulation model, cloud cover is persisted through the night using
! the last cloud cover value obtained before cutoff was reached,
! rather than interpolated over the night. -hlj-
!
! ORIGINAL COMMENTS ------------------------------------------------
! 6-23-92
! Modified from n_clearsky to use time series from met data 
! file rather than rather than calculating from start time and
! increment.  Renamed TIME(K) array to yday(k) array to avoid 
! confusion with time variable in cloud subroutine. 
! ORIGINAL COMMENTS ------------------------------------------------

!     set lat, lon of observations, compute lon-based time correction
!     ALAT = 59.5935
!     ALON = 20.9642
      COR = ALON / 360.

!     set desired ATC and CUTOFF
!     started with ATC=.83 , CUTOFF=.05 as with n_clearsky,
!     final version uses ATC=0.77, CUTOFF=0.10 
      ATC = 0.60       !changed over phone, Ji & Signell, 12/27/93
      CUTOFF = 0.05    !Ji & Signell, 12/27/93

!   be careful of sign of COR here (alon is positive W for this convention)

      TD = YDAY - (COR)

!MNT  CORRECT LOCAL TIME TO TRUE SOLAR TIME USING
!MNT  FIT TO TIME EQUATION DATA PROVIDED ON
!MNT  SMITHSONIAN MET TABLES, PAGE 495

      TEQN = -5.5 - 9. * COS((TD-45.6)*(6.283185/182.5))
      TEQN = TEQN / (24.*60.)
      TD = TD + TEQN
      Call N_CLOUD(CUTOFF,ALAT,TD,SWOBS,RAD,ATC,CLD)
      Return
      End
