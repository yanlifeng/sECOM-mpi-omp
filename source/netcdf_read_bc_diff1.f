      SUBROUTINE NETCDF_READ_BC_DIFF1(M,M_MAX,STR)
! Read diffuser 1 boundary conditions netcdf file
! A. Jordi (Feb 3, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      CHARACTER STR*80,STR_JUNK*31
      
      INTEGER START(2),COUNT(2)
      
      ! Open file
      ID_NC=NCOPN('run_data_bc_diff1.nc',NCNOWRIT,IERR)
      
      ! Read variables
      IF (M.EQ.0) THEN
        CALL NCAGTC(ID_NC,NCGLOBAL,'description',STR,80,IERR)
        ID_DIM=NCDID(ID_NC,'NUMDBC1',IERR)
        CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,NUMDBC1,IERR)
        
        IF (NUMDBC1.NE.0) THEN
          ID_DIM=NCDID(ID_NC,'NT',IERR)
          CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,M_MAX,IERR)
       
          START(1) = 1
          COUNT(1) = NUMDBC1
          ID_VAR=NCVID(ID_NC,'IDD',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,IDD,IERR)
          ID_VAR=NCVID(ID_NC,'JDD',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JDD,IERR)
          START(1) = 1
          START(2) = 1
          COUNT(1) = NUMDBC1
          COUNT(2) = KBM1
          ID_VAR=NCVID(ID_NC,'VDDIST',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,VDDIST(1:NUMDBC1,:),IERR)
        END IF
        
      ELSE
        START(1) = M
        COUNT(1) = 1
        ID_VAR=NCVID(ID_NC,'TIME',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,T2D,IERR)

        START(1) = 1
        START(2) = M
        COUNT(1) = NUMDBC1
        COUNT(2) = 1
        ID_VAR=NCVID(ID_NC,'QDIFF',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DQDIFF(1:NUMDBC1,2),IERR)
        ID_VAR=NCVID(ID_NC,'TDIFF',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DTDIFF(1:NUMDBC1,2),IERR)
        ID_VAR=NCVID(ID_NC,'SDIFF',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DSDIFF(1:NUMDBC1,2),IERR)
      END IF

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
