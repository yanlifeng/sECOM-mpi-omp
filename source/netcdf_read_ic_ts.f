      SUBROUTINE NETCDF_READ_IC_TS(IKSL,STR1,STR2)
! Read temperature and salinity initial conditions netcdf file
! A. Jordi (Feb 4, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      CHARACTER NAME*31,STR1*80,STR2*80
      INTEGER START(2),COUNT(2)
      
      ! Open file
      ID_NC=NCOPN('run_data_ic_ts.nc',NCNOWRIT,IERR)
      CALL NCAGTC(ID_NC,NCGLOBAL,'description',STR2,80,IERR)
      
      ! Read variables
      ID_DIM=NCDID(ID_NC,'KSL',IERR)
      CALL NCDINQ(ID_NC,ID_DIM,NAME,IKSL,IERR)
            
      START(1) = 1
      COUNT(1) = IKSL
      ID_VAR=NCVID(ID_NC,'DPTHSL',IERR)
      CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DPTHSL,IERR)
      CALL NCAGTC(ID_NC,ID_VAR,'description',STR1,80,IERR)
      
      START(1) = 1
      COUNT(1) = 20
      ID_VAR=NCVID(ID_NC,'OPTTSI',IERR)
      CALL NCVGTC(ID_NC,ID_VAR,START,COUNT,OPTTSI,20,IERR)
      
      IF(OPTTSI(1:5).EQ.'FIXED') THEN
        START(1) = 1
        COUNT(1) = IKSL
        ID_VAR=NCVID(ID_NC,'TSI',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,TSI,IERR)
        ID_VAR=NCVID(ID_NC,'SSI',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,SSI,IERR)
      END IF

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
