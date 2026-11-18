      SUBROUTINE NETCDF_READ_BC_ELEV(M,M_MAX,STR)
! Read elevation boundary conditions netcdf file
! A. Jordi (Feb 3, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      CHARACTER STR*80,STR_JUNK*31
      
      INTEGER START(2),COUNT(2)
      
      ! Open file
      ID_NC=NCOPN('run_data_bc_elev.nc',NCNOWRIT,IERR)
      
      ! Read variables
      IF (M.EQ.0) THEN
        CALL NCAGTC(ID_NC,NCGLOBAL,'description',STR,80,IERR)
        ID_DIM=NCDID(ID_NC,'NT',IERR)
        CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,M_MAX,IERR)   
        ID_DIM=NCDID(ID_NC,'NUMEBC',IERR)
        CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,NUMEBC,IERR)           
        START(1) = 1
        COUNT(1) = 20
        ID_VAR=NCVID(ID_NC,'OPTEBC',IERR)
        CALL NCVGTC(ID_NC,ID_VAR,START,COUNT,OPTEBC,20,IERR)
        START(1) = 1
        COUNT(1) = NUMEBC
        ID_VAR=NCVID(ID_NC,'IETA',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,IETA,IERR)
        ID_VAR=NCVID(ID_NC,'JETA',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JETA,IERR)
        ID_VAR=NCVID(ID_NC,'ICON',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,ICON,IERR)
        ID_VAR=NCVID(ID_NC,'JCON',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JCON,IERR)        
        IF(OPTEBC(1:4).NE.'DATA') THEN
          START(1) = 1
          COUNT(1) = NUMEBC
          ID_VAR=NCVID(ID_NC,'EMEAN',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,EMEAN,IERR)
          START(1) = 1
          START(2) = 1
          COUNT(1) = NUMEBC
          COUNT(2) = 6
          ID_VAR=NCVID(ID_NC,'AMP',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,AMP(1:NUMEBC,:),IERR)
          ID_VAR=NCVID(ID_NC,'PHASE',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,PHASE(1:NUMEBC,:),IERR)
        END IF
      
      ELSE
        START(1) = M
        COUNT(1) = 1
        ID_VAR=NCVID(ID_NC,'TIME',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,T2E,IERR)
        START(1) = 1
        START(2) = M
        COUNT(1) = NUMEBC
        COUNT(2) = 1
        ID_VAR=NCVID(ID_NC,'EBDRY',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DEBDRY(:,2),IERR)
        
      END IF

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
