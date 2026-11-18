      SUBROUTINE NETCDF_READ_BC_DIS(M,M_MAX,STR)
! Read discharge boundary conditions netcdf file
! A. Jordi (Feb 5, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      CHARACTER STR*80,STR_JUNK*31
      
      INTEGER START(2),COUNT(2)
      
      ! Open file
      ID_NC=NCOPN('run_data_bc_dis.nc',NCNOWRIT,IERR)
      
      ! Read variables
      IF (M.EQ.0) THEN
        CALL NCAGTC(ID_NC,NCGLOBAL,'description',STR,80,IERR)
        ID_DIM=NCDID(ID_NC,'NUMQBC',IERR)
        CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,NUMQBC,IERR)
        
        IF (NUMQBC.NE.0) THEN
          ID_DIM=NCDID(ID_NC,'NT',IERR)
          CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,M_MAX,IERR)
          START(1) = 1
          COUNT(1) = NUMQBC
          ID_VAR=NCVID(ID_NC,'IQD',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,IQD,IERR)
          ID_VAR=NCVID(ID_NC,'JQD',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JQD,IERR)
          ID_VAR=NCVID(ID_NC,'IQC',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,IQC,IERR)
          ID_VAR=NCVID(ID_NC,'JQC',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JQC,IERR)
          START(1) = 1
          START(2) = 1
          COUNT(1) = NUMQBC
          COUNT(2) = KBM1
          ID_VAR=NCVID(ID_NC,'VQDIST',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,VQDIST,IERR)
        END IF
        
      ELSE
        START(1) = M
        COUNT(1) = 1
        ID_VAR=NCVID(ID_NC,'TIME',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,T2D,IERR)

        START(1) = 1
        START(2) = M
        COUNT(1) = NUMQBC
        COUNT(2) = 1
        ID_VAR=NCVID(ID_NC,'QDIS',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DQDIS(:,2),IERR)
        ID_VAR=NCVID(ID_NC,'TDIS',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DTDIS(:,2),IERR)
        ID_VAR=NCVID(ID_NC,'SDIS',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DSDIS(:,2),IERR)
      END IF

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
