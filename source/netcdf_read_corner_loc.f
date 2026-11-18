      SUBROUTINE NETCDF_READ_CORNER_LOC
! Read corner_loc netcdf file
! A. Jordi (Jan 30, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      INTEGER START(2),COUNT(2)
      
      ! Open file
      ID_NC=NCOPN('corner_loc.nc',NCNOWRIT,IERR)
      
      ! Read variables
      START(1) = 1
      START(2) = 1
      COUNT(1) = IM+1
      COUNT(2) = JM+1
      ID_VAR=NCVID(ID_NC,'XCORNER',IERR)
      CALL NCVGT(ID_NC,ID_VAR,START,COUNT,XCORNER,IERR)
      ID_VAR=NCVID(ID_NC,'YCORNER',IERR)
      CALL NCVGT(ID_NC,ID_VAR,START,COUNT,YCORNER,IERR)

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
