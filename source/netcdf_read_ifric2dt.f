      SUBROUTINE NETCDF_READ_IFRIC2DT(M)
! Read space-variable ice friction input netcdf file
! A. Jordi (Jan 30, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      INTEGER START(3),COUNT(3)
      
      ! Open file
      ID_NC=NCOPN('ifric2dt.nc',NCNOWRIT,IERR)
      
      ! Read variables
      START(1) = M
      COUNT(1) = 1
      ID_VAR=NCVID(ID_NC,'TIME',IERR)
      CALL NCVGT(ID_NC,ID_VAR,START,COUNT,T2I,IERR)

      START(1) = 1
      START(2) = 1
      START(3) = M
      COUNT(1) = IM
      COUNT(2) = JM
      COUNT(3) = 1
      ID_VAR=NCVID(ID_NC,'DAREAICE',IERR)
      CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DAREAICE(:,:,2),IERR)
      ID_VAR=NCVID(ID_NC,'DTHICKICE',IERR)
      CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DTHICKICE(:,:,2),IERR)

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
