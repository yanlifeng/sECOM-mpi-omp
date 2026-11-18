      SUBROUTINE NETCDF_READ_MODEL_GRID(FROM_SUBR)
! Read model_grid netcdf file
! A. Jordi (Jan 29, 2015)

      INCLUDE 'comdeck'
      INCLUDE 'netcdf.inc'
      
      CHARACTER*6 FROM_SUBR
      CHARACTER*80 STR_GRID,STR_LEVELS,STR_METRICS,STR_TDAM
      DIMENSION ISTDAM(NUMDAM),JSTDAM(NUMDAM),NTDAM(NUMDAM)
      CHARACTER DIRTDAM(NUMDAM)*4
      
      INTEGER START(2),COUNT(2)
      CHARACTER*31 STR
      
      ! Open file
      ID_NC=NCOPN('model_grid.nc',NCNOWRIT,IERR)
      
      ! Establish depth array
      CALL NCAGTC(ID_NC,NCGLOBAL,'GRID',STR_GRID,80,IERR)
      IF (FROM_SUBR.EQ.'SETDOM') WRITE(IUPRT,'(/80A/)'),STR_GRID
      
      ! Initialize sigma levels
      CALL NCAGTC(ID_NC,NCGLOBAL,'LEVELS',STR_LEVELS,80,IERR)
      WRITE(IUPRT,'(/80A/)'),STR_LEVELS
      
      ID_DIM=NCDID(ID_NC,'KB',IERR)
      CALL NCDINQ(ID_NC,ID_DIM,STR,IKB,IERR)      
      WRITE(IUPRT,'('' KB = '',I5,/)') IKB
      IF(IKB.GT.101) THEN
        WRITE(6,41) KBM1
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
   41 FORMAT(//' NUMBER OF MODEL LAYERS',I5,' (KBM1)'/
     .         '           GREATER THAN'/
     .         ' MAXIMUM NUMBER OF VERTICAL DISTRIBUTION LAYERS'/
     .         ' SPECIFIED IN run_data DISCHARGE INFORMATION'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)

      IF(IKB.NE.KB) THEN
        WRITE(6,42) IKB,KB
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
   42 FORMAT(//' NUMBER OF SIGMA LEVELS IN MODEL_GRID',I5,' (IKB)'/
     .         '           NOT EQUAL TO'/
     .         ' NUMBER OF SIGMA LEVELS IN COMDECK   ',I5,' (KB)'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)

      IF(IKB.LT.3 .AND. TOR.NE.'BAROTROPIC') THEN
        WRITE(6,43) IKB,TOR
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
   43 FORMAT(//' IKB/KB = ',I5,' NOT ALLOWED IN A ',A10,' RUN'/
     .         ' NUMBER OF SIGMA LEVELS MUST BE GREATER THAN 2'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
     
      IF (FROM_SUBR.EQ.'SETDOM') THEN
        START(1) = 1
        COUNT(1) = IKB
        ID_VAR=NCVID(ID_NC,'Z',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,Z,IERR)
        WRITE(IUPRT,'(8F10.5)') (Z(K),K=1,IKB)

        DO K=1,KBM1
          DZ(K)=Z(K)-Z(K+1)
          DZR(K)=1./DZ(K)
          ZZ(K)=.5*(Z(K)+Z(K+1))
        END DO
        DO K=1,KBM2
          DZZ(K)=ZZ(K)-ZZ(K+1)
        END DO
        DZZ(KBM1)=0.0
        DZ(KB)=0.0
      END IF
      
      ! Define the metrics of the coordinate transformation
      CALL NCAGTC(ID_NC,NCGLOBAL,'METRICS',STR_METRICS,80,IERR)
      IF (FROM_SUBR.EQ.'SETDOM') WRITE(IUPRT,'(/80A/)'),STR_METRICS
      
      ID_DIM=NCDID(ID_NC,'IM',IERR)
      CALL NCDINQ(ID_NC,ID_DIM,STR,IIX,IERR)
      ID_DIM=NCDID(ID_NC,'JM',IERR)
      CALL NCDINQ(ID_NC,ID_DIM,STR,IJY,IERR)
      WRITE(IUPRT,'(('' IM = '',I5,/'' JM = '',I5))') IIX, IJY
      IF(IIX.NE.IM) THEN
        WRITE(6,8) IIX,IM
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
    8 FORMAT (//'     MODEL_GRID I-INDEX',I5,' (IIX)',/
     .          '        DOES NOT EQUAL'/
     .          '     COMDECK    I-INDEX',I5,' (IM)'/
     .          ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
      IF(IJY.NE.JM) THEN
        WRITE(6,9) IJY,JM
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
    9 FORMAT (//'     MODEL_GRID J-INDEX',I5,' (IJY)',/
     .          '        DOES NOT EQUAL'/
     .          '     COMDECK    J-INDEX',I5,' (JM)'/
     .          ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
     
      ! Set grid parameters
      IF (FROM_SUBR.EQ.'SETDOM') THEN
        START(1) = 1
        START(2) = 1
        COUNT(1) = IM
        COUNT(2) = JM
        ID_VAR=NCVID(ID_NC,'H1',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,H1,IERR)
        ID_VAR=NCVID(ID_NC,'H2',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,H2,IERR)
        ID_VAR=NCVID(ID_NC,'H',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,H,IERR)
        ID_VAR=NCVID(ID_NC,'ANG',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,ANG,IERR)
        ID_VAR=NCVID(ID_NC,'COR',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,COR,IERR)
        ID_VAR=NCVID(ID_NC,'XGRID',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,XGRID,IERR)
        ID_VAR=NCVID(ID_NC,'DATUM',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,DATUM,IERR)
        DO I=1,IM
          DO J=1,JM
            YGRID(I,J)=COR(I,J)
            H(I,J)=H(I,J)+WETMIN ! CNG09082011 MAKE GRID DEEPER BY WETMIN
          END DO
        END DO
      
        ! Set thin dam parameter (if provided)
        CALL NCINQ(ID_NC,NDIMS,NVARS,NATTS,NDIM2,IERR)
        IF (NDIMS.GT.3) THEN
          CALL NCAGTC(ID_NC,NCGLOBAL,'TDAM',STR_TDAM,80,IERR)
          WRITE(IUPRT,'(/80A/)'),STR_TDAM
        
          ID_DIM=NCDID(ID_NC,'NUMTDAM',IERR)
          CALL NCDINQ(ID_NC,ID_DIM,STR,NUMTDAM,IERR)
          IF(NUMTDAM.GT.NUMDAM)THEN
            WRITE(IUPRT,*)'NUMBER OF THIN DAM IS GREATER THAN ',NUMDAM
            WRITE(IUPRT,*)'PLEASE INCREASE THE ARRAY SIZE IN comdeck'
            STOP
          END IF
     
          START(1) = 1
          COUNT(1) = NUMTDAM
          ID_VAR=NCVID(ID_NC,'ISTDAM',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,ISTDAM,IERR)
          ID_VAR=NCVID(ID_NC,'JSTDAM',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,JSTDAM,IERR)
          ID_VAR=NCVID(ID_NC,'NTDAM',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,NTDAM,IERR)
     
          START(1) = 1
          START(1) = 1
          COUNT(1) = NUMTDAM
          COUNT(1) = NUMTDAM
          ID_VAR=NCVID(ID_NC,'ISTDAM',IERR)
          CALL NCVGT(ID_NC,ID_VAR,START,COUNT,ISTDAM,IERR)

          ID_VAR=NCVID(ID_NC,'DIRTDAM',IERR)
          DO N=1,NUMTDAM
            START(1) = 1
            START(2) = N
            COUNT(1) = 4
            COUNT(2) = 1
            CALL NCVGTC(ID_NC,ID_VAR,START,COUNT,DIRTDAM(N),4,IERR)
          END DO
        END IF
      ELSE IF(FROM_SUBR.EQ.'CORNER') THEN
        START(1) = 1
        START(2) = 1
        COUNT(1) = IM
        COUNT(2) = JM
        ID_VAR=NCVID(ID_NC,'COR',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,YGRID,IERR)
        ID_VAR=NCVID(ID_NC,'XGRID',IERR)
        CALL NCVGT(ID_NC,ID_VAR,START,COUNT,XGRID,IERR)
      END IF

      ! Close file
      CALL NCCLOS(ID_NC,IERR)
      
      RETURN
      END
