      SUBROUTINE NETCDF_READ_BC_TS(M,M_MAX,STR)
! Read temperature and salinity boundary conditions netcdf file
! A. Jordi (Feb 3, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      CHARACTER STR*80,STR_JUNK*31
      
      INTEGER START(3),COUNT(3)
      
      CHARACTER TBCOPT*5
      SAVE TBCOPT
      DIMENSION ZM(KBM1),T1(KSL),T2(KBM1),S1(KSL),S2(KBM1)
      
      ! Open file
      ID_NC=NCOPN('run_data_bc_ts.nc',NCNOWRIT,IERR)
      
      ! Read variables
      IF (M.EQ.0) THEN
        CALL NCAGTC(ID_NC,NCGLOBAL,'description',STR,80,IERR)
        ID_DIM=NCDID(ID_NC,'NT',IERR)
        CALL NCDINQ(ID_NC,ID_DIM,STR_JUNK,M_MAX,IERR)      
        START(1) = 1
        COUNT(1) = 5
        ID_VAR=NCVID(ID_NC,'TBCOPT',IERR)
        CALL NCVGTC(ID_NC,ID_VAR,START,COUNT,TBCOPT,5,IERR)
      
      ELSE
        START(1) = M
        COUNT(1) = 1
        ID_VAR=NCVID(ID_NC,'TIME',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,T2TS,IERR)
      
        IF(TBCOPT.NE.'SIGMA') THEN      
          START(1) = 1
          START(2) = M
          COUNT(1) = NUMEBC
          COUNT(2) = 1
          ID_VAR=NCVID(ID_NC,'ITAS',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,ITAS,IERR)
          ID_VAR=NCVID(ID_NC,'JTAS',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JTAS,IERR)
          START(1) = 1
          START(2) = 1
          START(3) = M
          COUNT(1) = NUMEBC
          COUNT(2) = KSL
          COUNT(3) = 1
          ID_VAR=NCVID(ID_NC,'TBDRYSL',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,TBDRYSL(1:NUMEBC,:),IERR)
          ID_VAR=NCVID(ID_NC,'SBDRYSL',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,SBDRYSL(1:NUMEBC,:),IERR)
        
          DO N=1,NUMEBC
            II=ITAS(N)
            JJ=JTAS(N)
            DO K=1,KBM1
              ZM(K)=ZZ(K)*H(II,JJ)
            END DO
            DO K=1,KSL
              T1(K)=TBDRYSL(N,K)
              S1(K)=SBDRYSL(N,K)
            END DO
            CALL SINTER(DPTHSL,T1,ZM,T2,KSL,KBM1)
            CALL SINTER(DPTHSL,S1,ZM,S2,KSL,KBM1)
            DO K=1,KBM1
              DTBDRY(N,K,2)=T2(K)
              DSBDRY(N,K,2)=S2(K)
            END DO
          END DO
        
        ELSE
          START(1) = 1
          START(2) = 1
          START(3) = M
          COUNT(1) = NUMEBC
          COUNT(2) = KBM1
          COUNT(3) = 1
          ID_VAR=NCVID(ID_NC,'TBDRY',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DTBDRY(1:NUMEBC,:,2),IERR)
          ID_VAR=NCVID(ID_NC,'SBDRY',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DSBDRY(1:NUMEBC,:,2),IERR)
        
        END IF
        
      END IF

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
