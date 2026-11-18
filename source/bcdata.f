      SUBROUTINE BCDATA(RESTAR)
!
      INCLUDE 'comdeck'
      SAVE
!
      DIMENSION COM(80) ! CNG04142008 ,YC(EBCM)
!NG      DIMENSION COM(80),YCC(EBCM)
      DIMENSION ZM(KBM1),T1(KSL),T2(KBM1),S1(KSL),S2(KBM1)
      DIMENSION IDBC(DBCM),TEMPD(DBCM)
      CHARACTER RESTAR*10
      CHARACTER*80 WAVE_HEADER
      LOGICAL SET_ICE ! A.Jordi


*Khan 072998ƒ
      DIMENSION SHFLX(IM,JM)    ! SYNOPTIC HEAT FLUX (in Watt/M**2)
*Khan 072998
!
!  FOR CONSERVATIVE TRACER:  OPEN B.C.
!
      REAL C1D(KSL),C2D(KBM1),WTEMP1(IM,JM),WTEMP2(IM,JM),WTEMP3(IM,JM)
      REAL*4 THZERO
      REAL WDS,WDD
!   TS BOUNDARY OPTION (EITHER BLANK (KSL) OR 'SIGMA')
      CHARACTER TBCOPT*5
      CHARACTER STR1*80,STR2*80 !A.Jordi: temporary strings
!
!-----------------------------------------------------------------------
!-------- READ IN STANDARD LEVELS HERE ---------------------------------
! A.Jordi: Read standard levels and initial condition from netcdf file
      CALL NETCDF_READ_IC_TS(IKSL,STR1,STR2)
      WRITE(IUPRT,'(/1X,A80/)') STR1
      WRITE(IUPRT,'(A,I5,/)') ' KSL = ',IKSL
      IF(IKSL.NE.KSL) THEN
        WRITE(6,42) IKSL,KSL ! FIXME: write to screen (6) or IUPRT file?
        CALL SYSTEM('rm gcm_temp*') ! FIXME: need to rm gcm_temp*?
        STOP
      ENDIF
 42   FORMAT(//' NUMBER OF STANDARD LEVELS IN RUN_DATA',I5,' (IKSL)'/
     .         '           DO NOT EQUAL'/
     .         ' NUMBER OF STANDARD LEVELS IN COMDECK ',I5,' (KSL)'/
     .         ' PLEASE CORRECT THIS PROBLEM AND TRY AGAIN'//)
      WRITE(IUPRT,'(8F10.1)') (DPTHSL(K),K=1,KSL)
      DO K=1,KSL ! FIX TO H ADJUSTED FOR WETMIN NG06242013
        DPTHSL(K)=DPTHSL(K)-WETMIN
      ENDDO
 11   FORMAT(80A1)
 12   FORMAT(/1X,80A1/)
 122  FORMAT(/1X,A80/)
!
!-----------------------------------------------------------------------
!-------- INITIAL TEMPERATURE AND SALINITY DATA ------------------------
!-------- TEMPERATURE IN C, SALINITY IN PSU ----------------------------
      WRITE(IUPRT,'(/1X,A80/)') STR2
      WRITE(IUPRT,'(A20)') OPTTSI
      IF(OPTTSI(1:5).EQ.'FIXED') THEN
        WRITE(IUPRT,'(8F10.1)') (TSI(K),K=1,KSL)
        WRITE(IUPRT,'(8F10.1)') (SSI(K),K=1,KSL)
      ENDIF
! A.Jordi: end
!
!-----------------------------------------------------------------------
!-------- READ IN BOUNDARY CONDITIONS HERE -----------------------------
!-----------------------------------------------------------------------
!-------- ELEVATION BOUNDARY -------------------------------------------
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 128   ! QA 9/18/98
! A.Jordi: Read elevation boundary condition from netcdf file
      NT_T2E=0      
      CALL NETCDF_READ_BC_ELEV(NT_T2E,NT_T2E_MAX,STR1)
      WRITE(IUPRT,'(/1X,A80/)') STR1
      WRITE(IUPRT,'(I5,1X,A20)') NUMEBC,OPTEBC
! A.Jordi: end
 24   FORMAT(I5,1X,A20)
!
      IF(NUMEBC.EQ.0) GO TO 128
      IF(NUMEBC.GT.EBCM) THEN
       WRITE(IUPRT,60) NUMEBC,EBCM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  60  FORMAT(//' NUMEBC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' EBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      IF(OPTEBC(1:4).NE.'DATA') THEN
      DO 110 N=1,NUMEBC
      WRITE(IUPRT,9) IETA(N),JETA(N),ICON(N),JCON(N),EMEAN(N)
      WRITE(IUPRT,7) (AMP(N,I),I=1,6)
      WRITE(IUPRT,7) (PHASE(N,I),I=1,6)
!     DO 105 I=1,6
!105  PHASE(N,I)=60.*PHASE(N,I)
 110  CONTINUE
!BNK   CALL FOREMAN'S TIDAL PREDICTION PROGRAM 
!     NDAY=IFIX(EDAY-SDAY)      ! NKIM (11.06.2001)
      NDAY=IFIX((IEND-ISTART)*DAYI*DTI)
!hli
!      write(700,*)NDAY
!hli
      THZERO=FLOAT(INTX)*DTI/3600.0
      WRITE(IUPRT,*)'THZERO=',THZERO
!NG04142008 MOVED TO NEW NOAA STANDARD FOR TIDAL CONSTITUENTS 
!NG04142008      DO N=1,NUMEBC
!NG04142008       YCC(N)=YGRID(IETA(N),JETA(N))
!NG04142008      END DO
!NG04142008      CALL PTIDE(EBCM,NUMEBC,7,YCC,AMP,PHASE,EMEAN,IDA,IMO,
!NG04142008     .   IYR,IHOUR,NDAY+1,THZERO)
!NG04142008 MOVED TO NEW NOAA STANDARD FOR TIDAL CONSTITUENTS 
      CALL PTIDE_NOAA(EBCM,NUMEBC,7,AMP,PHASE,EMEAN,IDA,IMO,
     .   IYR,IHOUR,NDAY+1,THZERO)
 9    FORMAT(4I5,1F10.5)
 7    FORMAT(6F10.5)
      ELSE
      WRITE(IUPRT,76) (IETA(N),JETA(N),ICON(N),JCON(N),N=1,NUMEBC)
 76   FORMAT(16I5)
! A.Jordi: Read elevation boundary condition from netcdf file
      DO M=1,NT_T2E_MAX
        CALL NETCDF_READ_BC_ELEV(M,M_MAX,STR1)
        WRITE(IUPRT,'(F10.1)') T2E
        OTIME=T2E
!NG BIG HARDWIRE HERE!!!!!!!
!       DO N=1,NUMEBC
!        EBDRY(N)=EBDRY(N)+2.133-abs(TIME/24.-1.5)/1.5*2.133
!       ENDDO
!NG BIG HARDWIRE ABOVE!!!!!!!
       WRITE(IUPRT,'(8F10.1)') (DEBDRY(N,2),N=1,NUMEBC)
      END DO
! A.Jordi: end
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE ETA BC TO COMPLETE THE SIMULATIONS ' 
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
      END IF
!  CHECK THE VALIDITY OF ELEVATION BOUNDARY CELLS
      IF(RESTAR.EQ.'COLD START')THEN
      DO N=1,NUMEBC
       IF(FSM(IETA(N),JETA(N)).EQ.0.0.OR.FSM(ICON(N),JCON(N)).EQ.0.)THEN
       WRITE(IUPRT,3801)N,IETA(N),JETA(N),ICON(N),JCON(N)
        STOP
       END IF
       IF((IETA(N)-ICON(N)).EQ.0.) THEN ! CONNECTED IN J DIR
         JDIFF=JETA(N)-JCON(N)
        IF(ABS(JDIFF).GT.1.OR.JDIFF.EQ.0)THEN
         WRITE(IUPLT,3802)N,IETA(N),JETA(N),ICON(N),JCON(N)
         STOP
        END IF
       ELSE    ! CONNECTED IN I DIR
         IDIFF=IETA(N)-ICON(N)
        IF(ABS(IDIFF).GT.1.OR.IDIFF.EQ.0)THEN
         WRITE(IUPLT,3802)N,IETA(N),JETA(N),ICON(N),JCON(N)
         STOP
        END IF
       END IF
      END DO
      ENDIF
 3801 FORMAT('BOUNDARY CELL NO',I5,'SPECIFIED AT',4I5,
     .  'IS IN LAND CELL',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
 3802 FORMAT('BOUNDARY CELL NO',I5,' AT',4I5,' IS SPECIFIED WRONG',/
     .  'PLEASE FIX THIS AND RESUBMIT THE RUN')
!
!-----------------------------------------------------------------------
!NG2012    CHECK FOR POSSIBLE CORNERS THAT WILL NEED SPECIAL TREATMENT
      NUMECR=0
      DO N=2,NUMEBC
       DO M=1,N-1
        IF(ICON(N).EQ.ICON(M).AND.JCON(N).EQ.JCON(M)) THEN ! If BC cells connect to same interior cell
!         IF(FSM(ICON(N)-1,JCON(N)-1).EQ.1) THEN ! BOTTOM LEFT CORNER?  CHECK...
          IF(ICON(N)-2.EQ.1.AND.JCON(N)-2.EQ.1) THEN ! ... Bottom left corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)-1
           JCHARD(NUMECR)=JCON(N)-1
           NCHARD(NUMECR)=1
           WRITE(*,*)'Found bottom left corner at:',ICON(N)-1,JCON(N)-1
           CALL SETDOM_C(ICON(N)-1,JCON(N)-1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(N).EQ.(ICON(N)-1)).AND.(JETA(N).EQ.JCON(N)).AND.
     .            (IETA(M).EQ.ICON(M)).AND.(JETA(M).EQ.JCON(M)-1)) THEN ! ... Bottom left corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)-1
           JCHARD(NUMECR)=JCON(N)-1
           NCHARD(NUMECR)=1
           WRITE(*,*)'Found bottom left corner at:',ICON(N)-1,JCON(N)-1
           CALL SETDOM_C(ICON(N)-1,JCON(N)-1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(M).EQ.(ICON(M)-1)).AND.(JETA(M).EQ.JCON(M)).AND.
     .            (IETA(N).EQ.ICON(N)).AND.(JETA(N).EQ.JCON(N)-1)) THEN ! ... Bottom left corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)-1
           JCHARD(NUMECR)=JCON(N)-1
           NCHARD(NUMECR)=1
           WRITE(*,*)'Found bottom left corner at:',ICON(N)-1,JCON(N)-1
           CALL SETDOM_C(ICON(N)-1,JCON(N)-1) ! SET FSM TO ZERO THERE
          ENDIF
!         ENDIF
!         IF(FSM(ICON(N)+1,JCON(N)-1).EQ.1) THEN ! BOTTOM RIGHT CORNER? CHECK...
          IF(ICON(N)+2.EQ.IM.AND.JCON(N)-2.EQ.1) THEN ! ... Bottom right corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)+1
           JCHARD(NUMECR)=JCON(N)-1
           NCHARD(NUMECR)=2
           WRITE(*,*)'Found bottom right corner at:',ICON(N)+1,JCON(N)-1
           CALL SETDOM_C(ICON(N)+1,JCON(N)-1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(N).EQ.(ICON(N)+1)).AND.(JETA(N).EQ.JCON(N)).AND.
     .            (IETA(M).EQ.ICON(M)).AND.(JETA(M).EQ.JCON(M)-1)) THEN ! ... Bottom right corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)+1
           JCHARD(NUMECR)=JCON(N)-1
           NCHARD(NUMECR)=2
           WRITE(*,*)'Found bottom right corner at:',ICON(N)+1,JCON(N)-1
           CALL SETDOM_C(ICON(N)+1,JCON(N)-1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(M).EQ.(ICON(M)+1)).AND.(JETA(M).EQ.JCON(M)).AND.
     .            (IETA(N).EQ.ICON(N)).AND.(JETA(N).EQ.JCON(N)-1)) THEN ! ... Bottom right corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)+1
           JCHARD(NUMECR)=JCON(N)-1
           NCHARD(NUMECR)=2
           WRITE(*,*)'Found bottom right corner at:',ICON(N)+1,JCON(N)-1
           CALL SETDOM_C(ICON(N)+1,JCON(N)-1) ! SET FSM TO ZERO THERE
          ENDIF
!         ENDIF
!         IF(FSM(ICON(N)+1,JCON(N)+1).EQ.1) THEN ! TOP RIGHT CORNER?  CHECK ...
          IF(ICON(N)+2.EQ.IM.AND.JCON(N)+2.EQ.JM) THEN ! ... Top right corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)+1
           JCHARD(NUMECR)=JCON(N)+1
           NCHARD(NUMECR)=3
           WRITE(*,*)'Found top right corner at:',ICON(N)+1,JCON(N)+1
           CALL SETDOM_C(ICON(N)+1,JCON(N)+1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(N).EQ.(ICON(N)+1)).AND.(JETA(N).EQ.JCON(N)).AND.
     .            (IETA(M).EQ.ICON(M)).AND.(JETA(M).EQ.JCON(M)+1)) THEN ! ... Top right corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)+1
           JCHARD(NUMECR)=JCON(N)+1
           NCHARD(NUMECR)=3
           WRITE(*,*)'Found top right corner at:',ICON(N)+1,JCON(N)+1
           CALL SETDOM_C(ICON(N)+1,JCON(N)+1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(M).EQ.(ICON(M)+1)).AND.(JETA(M).EQ.JCON(M)).AND.
     .            (IETA(N).EQ.ICON(N)).AND.(JETA(N).EQ.JCON(N)+1)) THEN ! ... Top right corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)+1
           JCHARD(NUMECR)=JCON(N)+1
           NCHARD(NUMECR)=3
           WRITE(*,*)'Found top right corner at:',ICON(N)+1,JCON(N)+1
           CALL SETDOM_C(ICON(N)+1,JCON(N)+1) ! SET FSM TO ZERO THERE
          ENDIF
!         ENDIF
!         IF(FSM(ICON(N)-1,JCON(N)+1).EQ.1) THEN ! TOP LEFT CORNER? CHECK ...
          IF(ICON(N)-2.EQ.1.AND.JCON(N)+2.EQ.JM) THEN ! ... Top left corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)-1
           JCHARD(NUMECR)=JCON(N)+1
           NCHARD(NUMECR)=4
           WRITE(*,*)'Found top left corner at:',ICON(N)-1,JCON(N)+1
           CALL SETDOM_C(ICON(N)-1,JCON(N)+1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(N).EQ.(ICON(N)-1)).AND.(JETA(N).EQ.JCON(N)).AND.
     .            (IETA(M).EQ.ICON(M)).AND.(JETA(M).EQ.JCON(M)+1)) THEN ! ... Top left corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)-1
           JCHARD(NUMECR)=JCON(N)+1
           NCHARD(NUMECR)=4
           WRITE(*,*)'Found top left corner at:',ICON(N)-1,JCON(N)+1
           CALL SETDOM_C(ICON(N)-1,JCON(N)+1) ! SET FSM TO ZERO THERE
          ELSEIF ((IETA(M).EQ.(ICON(M)-1)).AND.(JETA(M).EQ.JCON(M)).AND.
     .            (IETA(N).EQ.ICON(N)).AND.(JETA(N).EQ.JCON(N)+1)) THEN ! ... Top left corner found.
           NUMECR=NUMECR+1
           ICHARD(NUMECR)=ICON(N)-1
           JCHARD(NUMECR)=JCON(N)+1
           NCHARD(NUMECR)=4
           WRITE(*,*)'Found top left corner at:',ICON(N)-1,JCON(N)+1
           CALL SETDOM_C(ICON(N)-1,JCON(N)+1) ! SET FSM TO ZERO THERE
          ENDIF
!         ENDIF
        ENDIF
       ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!-------- TEMPERATURE AND SALINITY BOUNDARY ----------------------------
! A.Jordi: Read T/S boundary condition from netcdf file
      NT_T2TS=0      
      CALL NETCDF_READ_BC_TS(NT_T2TS,NT_T2TS_MAX,STR1)
      WRITE(IUPRT,'(/1X,A80/)') STR1
      DO M=1,NT_T2TS_MAX
        CALL NETCDF_READ_BC_TS(M,M_MAX,STR1)
        WRITE(IUPRT,'(F10.1)') T2TS
        OTIME=T2TS
        DO N=1,NUMEBC
          IF(TBCOPT.NE.'SIGMA') THEN
          !  MODIFIED BY C.K.Z. ON 11/10/94
            IF ((ITAS(N).NE.IETA(N)).OR.(JTAS(N).NE.JETA(N))) THEN
              WRITE(IUPRT,5701) N,ITAS(N),JTAS(N),IETA(N),JETA(N)
 5701         FORMAT (/5X,'******  MODEL EXECUTION STOPPED  *****',
     .               //5X,'ORDER OF TEMPERATURE/SALINITY OPEN B.C.',
     .                /5X,'SPECIFIED INCORRECTLY',
     .                /5X,'B.C. NUMBER ',I4,' HAS ITAS, JTAS =',2I4,
     .                /5X,'SHOULD BE IETA, JETA =',2I4,
     .              //5X,'PLEASE CORRECT AND RESUBMIT')
              STOP
            END IF
            WRITE(IUPRT,'(2I5)') ITAS(N),JTAS(N)
            WRITE(IUPRT,'(50F5.0)') ((TBDRYSL(N,K)),K=1,KSL)
            WRITE(IUPRT,'(50F5.0)') ((SBDRYSL(N,K)),K=1,KSL)
            
          ELSE
            IF(M.EQ.1) WRITE(IUPRT,*)
     .        'SIGMA LEVEL TS BOUNDARY DATA READING'
             IF ((ITAS(N).NE.IETA(N)).OR.(JTAS(N).NE.JETA(N))) THEN
                WRITE(IUPRT,5701) N,ITAS(N),JTAS(N),IETA(N),JETA(N)
                STOP
             END IF
             WRITE(IUPRT,'(2I5)') ITAS(N),JTAS(N)
             WRITE(IUPRT,'(50F5.0)') ((TBDRY(N,K)),K=1,KBM1)
             WRITE(IUPRT,'(50F5.0)') ((SBDRY(N,K)),K=1,KBM1)            
          END IF          
        END DO
      END DO
! A.Jordi: end
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE T/S BC TO COMPLETE THE SIMULATIONS ' 
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
79    FORMAT(2I5,100F5.0)
80    FORMAT(50F5.0)
!
!  FOR CONSERVATIVE TRACER:  OPEN B.C.
!
 128  IF (TRACER.EQ.'INCLUDE') THEN
        OPEN (IUT501,FILE='gcm_temp501')
!
        READ(IUT401,11)  (COM(I),I=1,80)
        WRITE(IUPRT,12) (COM(I),I=1,80)
        READ (IUT401,24)NUMEBCTR
!
        IF (NUMEBCTR.EQ.0) GOTO 1128
!
        DO 6010 M=1,100000
          READ(IUT401,77,ERR=6020) TIME
          WRITE(IUPRT,77) TIME
          OTIME=TIME
          WRITE (IUT501,78)TIME
          DO 6030 N=1,NUMEBCTR
            READ(IUT401,179) II,JJ,IIC,JJC,(CBDRYSL1(N,K),K=1,KSL)
 179        FORMAT (4I5,100F5.0)
!
            ITRED(N)=II
            JTRED(N)=JJ
            ITREC(N)=IIC
            JTREC(N)=JJC
!
            WRITE(IUPRT,9) II,JJ,IIC,JJC
            WRITE(IUPRT,80)((CBDRYSL1(N,K)),K=1,KSL)
            DO 6050 K=1,KBM1
              ZM(K)=ZZ(K)*H(II,JJ)
 6050       CONTINUE
            DO 6060 K=1,KSL
              C1D(K)=CBDRYSL1(N,K)
 6060       CONTINUE                
            CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
            DO 6070 K=1,KBM1
              CBDRY1(N,K)=C2D(K)
 6070       CONTINUE
!
            WRITE (IUT501,78) (CBDRY1(N,K),K=1,KBM1)
!
 6030     CONTINUE
 6010   CONTINUE
 6020   BACKSPACE IUT401
      ENDIF
!
!
!*******************************************************************
!
!  SEDIMENT TRANSPORT
!
 1128 IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (IUT502,FILE='gcm_temp502')
          READ(IUT402,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT402,24)NUMEBCSE
!
          IF (NUMEBCSE.EQ.0) GOTO 2127
!
          DO 6110 M=1,100000
            READ(IUT402,77,ERR=6120) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT502,78)TIME
            DO 6130 N=1,NUMEBCSE
              READ(IUT402,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ISEED(N)=II
              JSEED(N)=JJ
              ISEEC(N)=IIC
              JSEEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
!
!  CONVERT FROM mg/l TO g/cm**3
!
              DO 6150 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000000.
 6150         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6160 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6160         CONTINUE
!
              WRITE(IUT502,78) (CBDRY(1,N,K),K=1,KBM1)
!
 6130       CONTINUE
 6110     CONTINUE
 6120     BACKSPACE IUT402
        ENDIF
!
 2127   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          OPEN (IUT503,FILE='gcm_temp503')
!
          READ(IUT403,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT403,24)NUMEBCSE
!
          IF (NUMEBCSE.EQ.0) GOTO 2128
!
          DO 6170 M=1,100000
            READ(IUT403,77,ERR=6180) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT503,78)TIME
!
            DO 6190 N=1,NUMEBCSE
              READ(IUT403,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ISEED(N)=II
              JSEED(N)=JJ
              ISEEC(N)=IIC
              JSEEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
!
!  CONVERT FROM mg/l TO g/cm**3
!
              DO 6200 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000000.
 6200         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6210 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6210         CONTINUE
!
              WRITE(IUT503,78) (CBDRY(1,N,K),K=1,KBM1)
!
 6190       CONTINUE
 6170     CONTINUE
 6180     BACKSPACE IUT403
        ENDIF
      ENDIF
!
!*******************************************************************
!
!  CHEM TRANSPORT
!
 2128 IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT504=504
          OPEN (IUT504,FILE='gcm_temp504')
!
          READ(IUT404,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT404,24)NUMEBCCH
!
          IF (NUMEBCCH.EQ.0) GOTO 2129
!
          DO 6310 M=1,100000
            READ(IUT404,77,ERR=6320) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT504,78)TIME
!
            DO 6330 N=1,NUMEBCCH
              READ(IUT404,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ICHED(N)=II
              JCHED(N)=JJ
              ICHEC(N)=IIC
              JCHEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
!
!  CONVERT FROM ug/l TO ug/cm**3
!
              DO 6350 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000.
 6350         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6360 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6360         CONTINUE
!
              WRITE(IUT504,78) (CBDRY(1,N,K),K=1,KBM1)
!
 6330       CONTINUE
 6310     CONTINUE
 6320     BACKSPACE IUT404
        ENDIF
!
 2129   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT505=505
          OPEN (IUT505,FILE='gcm_temp505')
!
          READ(IUT405,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT405,24)NUMEBCCH
!
          IF (NUMEBCCH.EQ.0) GOTO 126 
!
          DO 6370 M=1,100000
            READ(IUT405,77,ERR=6380) TIME
            WRITE(IUPRT,77) TIME
            WRITE (IUT505,78)TIME
!
            DO 6390 N=1,NUMEBCCH
              READ(IUT405,179) II,JJ,IIC,JJC,(CBDRYSL(1,N,K),K=1,KSL)
              ICHED(N)=II
              JCHED(N)=JJ
              ICHEC(N)=IIC
              JCHEC(N)=JJC
              WRITE(IUPRT,9) II,JJ,IIC,JJC
              WRITE(IUPRT,80)((CBDRYSL(1,N,K)),K=1,KSL)
!
!  CONVERT FROM ug/l TO ug/cm**3
!
              DO 6400 K=1,KSL
                C1D(K)=CBDRYSL(1,N,K)/1000.
 6400         CONTINUE                
              CALL SINTER(DPTHSL,C1D,ZM,C2D,KSL,KBM1)
              DO 6410 K=1,KBM1
                CBDRY(1,N,K)=C2D(K)
 6410         CONTINUE
!
              WRITE(IUT505,78) (CBDRY(1,N,K),K=1,KBM1)
!
 6390       CONTINUE
 6370     CONTINUE
 6380     BACKSPACE IUT405
        ENDIF
      ENDIF
!
!*******************************************************************
!
 126  CONTINUE
!
!-----------------------------------------------------------------------
!-------- ELEVATION POTENTIAL CELLS NG03242008 -------------------------
      NUMEPC=0
      OPEN (70,FILE='epcor.inp',status='old',err=99999)
      READ(70,11)  (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(70,'(I5,1X,A8,1X,F10.0)')  NUMEPC,OPTEPC,TNGEPCIN
      IF(NUMEPC.EQ.0) THEN
       CLOSE(70)
       GO TO 99999
      ENDIF
      WRITE(IUPRT,'(I5," ELEVATION POTENTIAL CELLS FOUND...")') NUMEPC
      IF (OPTEPC(1:8).NE.'CONSTANT'.AND.OPTEPC(1:8).NE.'CLAMPED '.AND.
     +    OPTEPC(1:8).NE.'PCLAMPED'.AND.OPTEPC(1:8).NE.'TIDAL   ') THEN
          OPTEPC='TIDAL   '
      ENDIF
      IF(NUMEPC.GT.EPCM) THEN
       WRITE(IUPRT,600) NUMEPC,EPCM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  600  FORMAT(//' NUMEPC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .          ' EPCM   (=',I4,') SPECIFIED IN comdeck'/
     .          ' PLEASE FIX AND RESUBMIT'//)
      WRITE(IUPRT,'(2A5,2A10)')' IEPC',' JEPC','   EPCMEAN','  TNUDGEPC'
      DO N=1,NUMEPC
        READ(70,760)  IEPC(N),JEPC(N),EPCMEAN(N)
        I=IEPC(N)
        J=JEPC(N)
        IF     (OPTEPC(1:8).EQ.'CONSTANT'.and.TNGEPCIN.NE.0.0) THEN
          TNUDGEPC(N)=TNGEPCIN ! NUDGING TIMESTEP DEFINED IN SECONDS
        ELSEIF (OPTEPC(1:8).EQ.'CLAMPED ') THEN
          TNUDGEPC(N)=DTE ! USE EXTERNAL TIME STEP FOR CLAMPED
        ELSEIF (OPTEPC(1:8).EQ.'PCLAMPED') THEN
          IF (TNGEPCIN.NE.0.0) THEN
            TNUDGEPC(N)=TNGEPCIN ! PARTIALLY CLAMPED ALPHA DEFINED IN SECONDS
          ELSE ! USE A TIDAL TNUDGE
            CPHL=DUM(I,J)/(0.5*(H1(I,J)+H1(I-1,J)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I-1,J)))
            CPHR=DUM(I+1,J)/(0.5*(H1(I,J)+H1(I+1,J)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I+1,J)))
            CPHU=DVM(I,J+1)/(0.5*(H2(I,J)+H2(I,J+1)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I,J+1)))
            CPHD=DVM(I,J)/(0.5*(H2(I,J)+H2(I,J-1)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I,J-1)))
            IF ((CPHL+CPHR+CPHU+CPHD).NE.0.) THEN
             TNUDGEPC(N)=AMAX1(DTE,0.25/AMAX1(CPHL,CPHR,CPHU,CPHD)) ! USE A QUARTER OF THE TIDAL DISTANCE
            ELSE ! 1D SIMULATION, USE CLAMPED
             TNUDGEPC(N)=DTE
             OPTEPC='CLAMPED '
            ENDIF
          ENDIF
        ELSE ! DEFAULT IS TIDAL TNUDGEPC, EVEN FOR PCLAMPED
          IF (TNGEPCIN.EQ.0.0) TNGEPCIN=0.25 ! DEFAULT 1/4 WAVELENGTH FRACTION, OTHERWISE READ IN
          CPHL=DUM(I,J)/(0.5*(H1(I,J)+H1(I-1,J)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I-1,J)))
          CPHR=DUM(I+1,J)/(0.5*(H1(I,J)+H1(I+1,J)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I+1,J)))
          CPHU=DVM(I,J+1)/(0.5*(H2(I,J)+H2(I,J+1)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I,J+1)))
          CPHD=DVM(I,J)/(0.5*(H2(I,J)+H2(I,J-1)))*
     +            SQRT(GRAV*0.5*(D(I,J)+D(I,J-1)))
          IF ((CPHL+CPHR+CPHU+CPHD).NE.0.) THEN
           TNUDGEPC(N)=AMAX1(DTE,TNGEPCIN/AMAX1(CPHL,CPHR,CPHU,CPHD))! USE TNGEPCIN OF THE TIDAL DISTANCE
          ELSE ! 1D SIMULATION, USE CLAMPED
           TNUDGEPC(N)=DTE
           OPTEPC='CLAMPED '
          ENDIF
        ENDIF
        WRITE(IUPRT,761)  I,J,EPCMEAN(N),TNUDGEPC(N)
      ENDDO
  760 FORMAT(2I5,F10.5)
  761 FORMAT(2I5,2F10.5)
      OPEN(IUT65,FILE='gcm_temp65')
      DO M=1,999999 ! Allow for up to 11 years of 6-minute data (is also max of precision for FORMAT 77)
       READ(70,77,END=9999) TIME
!NG       WRITE(IUPRT,77) TIME
       OTIME=TIME
       READ(70,77)  (EPCOR(N),N=1,NUMEPC)
!NG       WRITE(IUPRT,77) (EPCOR(N),N=1,NUMEPC)
       DO N=1,NUMEPC
        IF ((EPCOR(N).EQ.-99.).OR.(EPCMEAN(N).EQ.-99.)) THEN ! FLAG FOR BAD DATA
         EPCOR(N)=-99.
         WRITE (*,'(a,F12.3,a,3i5)') 'NO EPC OBS FOUND FOR ',TIME,
     +                               ' AND ',N,IEPC(N),JEPC(N)
        ELSE
!NG       The extensive drying check below created water. So now I just don't apply EPCOR in a dry cell. 
!         EPCOR(N)=AMAX1(EPCOR(N)-EPCMEAN(N),-H(I,J)+WETMIN-1E-3) ! Convert datum and don't allow extensive drying.
         EPCOR(N)=EPCOR(N)-EPCMEAN(N)
        ENDIF
       ENDDO
       CALL PERCELL(NUMEPC,IEPC,JEPC,EPCOR,TNUDGEPC,-99.) ! Average any elevation observations prescribed in same cells
                                                          ! Redefines NUMEPC,IEPC,JEPC,and EPCOR=mean(EPCOR~=-99.)
       WRITE(IUT65,78) TIME
       WRITE(IUT65,78) (EPCOR(N),N=1,NUMEPC)
      ENDDO
 9999 CLOSE (70)
      WRITE(IUPRT,'(A,I5,A)')'FOUND ',NUMEPC,
     +                       ' UNIQUE POTENTIAL ELEVATION CELLS'
      WRITE(IUPRT,'("THAT WILL BE USING THE OPTEPC: ",A8)') OPTEPC
      WRITE(IUPRT,'(2A5,A10)')' IEPC',' JEPC','  TNUDGEPC'
      DO N=1,NUMEPC
        I=IEPC(N)
        J=JEPC(N)
        WRITE(IUPRT,761)  I,J,TNUDGEPC(N)
      ENDDO
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,'(3(/,A,/))')'** ** WARNING ** **',
     +   'ELEVATION POTENTIAL CELLS DATA ',
     +   'DO NOT REACH THE END OF THE SIMULATION.'
        WRITE(IUPRT,'(3(/,A,F14.5/))')
     +   'WILL ASSUME FREE SOLUTION AFTER HOUR ',TIME
        WRITE(IUT65,78) TIME+1
        WRITE(IUT65,78) (-99.,N=1,NUMEPC)
        WRITE(IUT65,78) 999999.
        WRITE(IUT65,78) (-99.,N=1,NUMEPC)
      ENDIF
!  CHECK THE VALIDITY OF ELEVATION POTENTIAL CELLS
!      IF(RESTAR.EQ.'COLD START')THEN
       DO N=1,NUMEPC
        IF(FSM(IEPC(N),JEPC(N)).EQ.0.)THEN
         WRITE(IUPRT,6501)N,IEPC(N),JEPC(N)
         STOP
        END IF
       END DO
!      END IF
 6501 FORMAT('ELEVATION POTENTIAL CELL NO',I5,'SPECIFIED AT',2I5,
     .  'IS IN LAND CELL',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
99999 IF ((NUMEPC).EQ.0) WRITE (IUPRT,
     .   '("NO ELEVATION POTENTIAL CORRECTION FILE FOUND")')
!
!-----------------------------------------------------------------------
!-------- RIVER/DAM AND ONSHORE INTAKE/OUTFALL DISCHARGE BOUNDARY ------
!-----------------------------------------------------------------------
!
      IF (HYDTYPE.EQ.'EXTERNAL') THEN   ! QA 9/18/98
         NUMEBC = 0   ! QA 9/18/98
         NUMQBC=0
         GOTO 138
      ENDIF
!
! A.Jordi: Read discharge boundary condition from netcdf file
      NT_T2Q=0
      CALL NETCDF_READ_BC_DIS(NT_T2Q,NT_T2Q_MAX,STR1)
      WRITE(IUPRT,'(/1X,A80/)') STR1
      WRITE(IUPRT,'(I5)') NUMQBC
! A.Jordi: end
!
      IF(NUMQBC.EQ.0) GO TO 138
      IF(NUMQBC.GT.QBCM) THEN
       WRITE(IUPRT,61) NUMQBC,QBCM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  61  FORMAT(//' NUMQBC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' QBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      DO 321 N=1,NUMQBC
321   WRITE(IUPRT,6) IQD(N),JQD(N),IQC(N),JQC(N),(VQDIST(N,K),K=1,KBM1)
! CHECK THE VALIDITY OF RIVER INFLOW LOCATIONS
      IF(RESTAR.EQ.'COLD START')THEN
       DO N=1,NUMQBC
        IF(FSM(IQD(N),JQD(N)).NE.1.0.OR.FSM(IQC(N),JQC(N)).NE.0.0)THEN
         WRITE(IUPRT,3811)N,IQD(N),JQD(N),IQC(N),JQC(N)
         STOP
        END IF
        DO M=1,NUMQBC
         IF(M.NE.N)THEN   ! CHECKING DUPLICATE SPECIFICATION OF RIVERS
          IF(IQD(N).EQ.IQD(M).AND.IQC(N).EQ.JQD(M).AND.IQC(N).EQ.IQC(M)
     .   .AND.JQC(N).EQ.JQC(M))THEN
           WRITE(IUPRT,3812)M,N
           STOP
          END IF
         END IF
        END DO
! CHECKING WHETHER RIVER IQC AND IQD, JQC AND JQD ARE ADJACENT
         IF(JQD(N).EQ.JQC(N))THEN  ! RIVER FROM EAST OR WEST 
          IF(ABS(IQD(N)-IQC(N)).NE.1)THEN
           WRITE(IUPRT,3811)N,IQD(N),JQD(N),IQC(N),JQC(N)
           STOP
          END IF
         END IF
         IF(IQD(N).EQ.IQC(N))THEN  ! RIVER FROM NORTH OR SOUTH
          IF(ABS(JQD(N)-JQC(N)).NE.1)THEN
           WRITE(IUPRT,3811)N,IQD(N),JQD(N),IQC(N),JQC(N)
           STOP
          END IF
         END IF
! CHECKING VDIST OF RIVER INFLOW
          SUMVQDIST=0.0
         DO K=1,KBM1
          SUMVQDIST=SUMVQDIST+VQDIST(N,K)
         END DO

         IF(SUMVQDIST.NE.100.0)THEN
          WRITE(IUPRT,3813)N
          STOP
         END IF

       END DO
       END IF
 
 3811 FORMAT('RIVER INFLOW NO',I5,' SPECIFIED AT',4I5,
     .  ' IS INCORRECT',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
 3812 FORMAT('RIVER INFLOW NO',I5,' AND NO ',I5,
     .  ' ARE SPECIFIED AT THE SAME GRID CELLS',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
 3813 FORMAT('RIVER INFLOW NO',I5,' AND NO ',I5,
     .  ' HAS WRONG VQDIST: THE SUM SHOULD BE 100',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
 8    FORMAT(2I5,4F10.5)
 5    FORMAT(4I5,100F5.0)
 6    FORMAT(4I5,/,100F5.1)
! A.Jordi: Read discharge boundary condition from netcdf file
      DO M=1,NT_T2Q_MAX
        CALL NETCDF_READ_BC_DIS(M,NT_T2Q_MAX,STR1)
        WRITE(IUPRT,'(8F10.1)') T2Q
        WRITE(IUPRT,'(8F10.1)') (QDIS(N),N=1,NUMQBC)
        WRITE(IUPRT,'(8F10.1)') (TDIS(N),N=1,NUMQBC)
        WRITE(IUPRT,'(8F10.1)') (SDIS(N),N=1,NUMQBC)
        OTIME=T2Q
      END DO
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE RIVER BC TO COMPLETE THE SIMULATIONS' 
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
!
!*******************************************************************
!
!  FOR CONSERVATIVE TRACER:  RIVER DISCHARGE
!
 138  IF (TRACER.EQ.'INCLUDE') THEN
        OPEN (IUT601,FILE='gcm_temp601')
!
        READ(IUT401,11)  (COM(I),I=1,80)
        WRITE(IUPRT,12) (COM(I),I=1,80)
        READ (IUT401,24)NUMQBCTR
        WRITE (IUPRT,24)NUMQBCTR
!
        IF (NUMQBCTR.GT.NUMQBC)THEN
         WRITE(IUPRT,*)'NUMQBCTR CANNOT BE LARGER THAN NUMQBC'
         CALL SYSTEM ('rm gcm_temp*')
         STOP
        END IF
!
        IF (NUMQBCTR.EQ.0) GOTO 1129
!
        DO 6509 N=1,NUMQBCTR
          READ (IUT401,179)ITRQD(N),JTRQD(N),ITRQC(N),JTRQC(N)
          WRITE (IUPRT,179)ITRQD(N),JTRQD(N),ITRQC(N),JTRQC(N)
 6509   CONTINUE

!   FIND PROPER FLOW MATCH
          DO 6508 NQ=1,NUMQBC
           IF(ITRQD(NQ).NE.IQD(NQ).OR.JTRQD(NQ).NE.JQD(NQ).OR.
     .        ITRQC(NQ).NE.IQC(NQ).OR.JTRQC(NQ).NE.JQC(NQ)) THEN
              WRITE(IUPRT,2480) ITRQD(N),JTRQD(N)
              CALL SYSTEM ('rm gcm_temp*')
              STOP
           END IF
 6508     CONTINUE
 2480 FORMAT(' ONE OF RIVER TRACER INPUT DOES NOT MATCH WITH RIVER
     .INFLOW LOCATIONS. ITS I,J ARE',2I5)
!
        DO 6510 M=1,100000
          READ(IUT401,77,ERR=6520) TIME
          READ(IUT401,77)  (CDIS1(N),N=1,NUMQBCTR)
          OTIME=TIME
          WRITE(IUPRT,77) TIME
          WRITE (IUPRT,77)  (CDIS1(N),N=1,NUMQBCTR)
          WRITE (IUT601,78)TIME
          WRITE (IUT601,78)  (CDIS1(N),N=1,NUMQBCTR)
!
 6510   CONTINUE
 6520   BACKSPACE IUT401
        IF(OTIME/24.LE.(EDAY-SDAY))THEN
          WRITE(IUPRT,*)
     .    'NEED MORE TRACER(RIVER) BC TO COMPLETE THE SIMULATIONS'
          CALL SYSTEM ('rm gcm_temp*')
        STOP
        ENDIF
      ENDIF
!
!*******************************************************************
!
!  SEDIMENT TRANSPORT
!
 1129 IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT602=602
          OPEN (IUT602,FILE='gcm_temp602')
!
          READ(IUT402,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT402,24)NUMQBCSE
          WRITE (IUPRT,24)NUMQBCSE
!
          IF (NUMQBCSE.EQ.0) GOTO 1229
!
          DO 6529 N=1,NUMQBCSE
            READ (IUT402,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
            WRITE (IUPRT,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
 6529     CONTINUE
!
          DO 6530 M=1,100000
            READ(IUT402,77,ERR=6540) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT402,77)  (CDIS(1,N),N=1,NUMQBCSE)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCSE)
!
!  CONVERT FROM mg/l TO g/cm**3
!
            DO 6550 NN=1,NUMQBCSE
              CDIS(1,NN)=CDIS(1,NN)/1000000.
 6550       CONTINUE
!
            WRITE (IUT602,78)TIME
            WRITE (IUT602,78)  (CDIS(1,N),N=1,NUMQBCSE)
!
 6530     CONTINUE
 6540     BACKSPACE IUT402
        ENDIF
!
 1229   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT603=603
          OPEN (IUT603,FILE='gcm_temp603')
!
          READ(IUT403,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT403,24)NUMQBCSE
          WRITE (IUPRT,24)NUMQBCSE
!
          IF (NUMQBCSE.EQ.0) GOTO 1239
!
          DO 6559 N=1,NUMQBCSE
            READ (IUT403,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
            WRITE (IUPRT,179)ISEQD(N),JSEQD(N),ISEQC(N),JSEQC(N)
 6559     CONTINUE
!
          DO 6560 M=1,100000
            READ(IUT403,77,ERR=6570) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT403,77)  (CDIS(1,N),N=1,NUMQBCSE)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCSE)
!
!  CONVERT FROM mg/l TO g/cm**3
!
            DO 6580 NN=1,NUMQBCSE
              CDIS(1,NN)=CDIS(1,NN)/1000000.
 6580       CONTINUE
!
            WRITE (IUT603,78)TIME
            WRITE (IUT603,78)  (CDIS(1,N),N=1,NUMQBCSE)
!
 6560     CONTINUE
 6570     BACKSPACE IUT403
        ENDIF
      ENDIF
!
!*******************************************************************
!
!  CHEMICAL TRANSPORT
!
 1239 IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT604=604
          OPEN (IUT604,FILE='gcm_temp604')
!
          READ(IUT404,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT404,24)NUMQBCCH
          WRITE (IUPRT,24)NUMQBCCH
!
          IF (NUMQBCCH.EQ.0) GOTO 1329
!
          DO 6589 N=1,NUMQBCCH
            READ (IUT404,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
            WRITE (IUPRT,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
 6589     CONTINUE
!
          DO 6590 M=1,100000
            READ(IUT404,77,ERR=6600) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT404,77)  (CDIS(1,N),N=1,NUMQBCCH)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCCH)
!
!  CONVERT FROM ug/l TO ug/cm**3
!
            DO 6610 NN=1,NUMQBCCH
              CDIS(1,NN)=CDIS(1,NN)/1000.
 6610       CONTINUE
!
            WRITE (IUT604,78)TIME
            WRITE (IUT604,78)  (CDIS(1,N),N=1,NUMQBCCH)
!
 6590     CONTINUE
 6600     BACKSPACE IUT404
        ENDIF
!
 1329   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT605=605
          OPEN (IUT605,FILE='gcm_temp605')
!
          READ(IUT405,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT405,24)NUMQBCCH
          WRITE (IUPRT,24)NUMQBCCH
!
          IF (NUMQBCCH.EQ.0) GOTO 1330
!
          DO 6599 N=1,NUMQBCCH
            READ (IUT405,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
            WRITE (IUPRT,179)ICHQD(N),JCHQD(N),ICHQC(N),JCHQC(N)
 6599     CONTINUE
!
          DO 6620 M=1,100000
            READ(IUT405,77,ERR=6630) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT405,77)  (CDIS(1,N),N=1,NUMQBCCH)
            WRITE (IUPRT,77)  (CDIS(1,N),N=1,NUMQBCCH)
!
!  CONVERT FROM ug/l TO ug/cm**3
!
            DO 6640 NN=1,NUMQBCCH
              CDIS(1,NN)=CDIS(1,NN)/1000.
 6640       CONTINUE
!
            WRITE (IUT605,78)TIME
            WRITE (IUT605,78)  (CDIS(1,N),N=1,NUMQBCCH)
!
 6620     CONTINUE
 6630     BACKSPACE IUT405
        ENDIF
      ENDIF
!
!*******************************************************************
!
 1330 CONTINUE
      DO 140 N=1,NUMQBC
      ID=IQD(N)
      JD=JQD(N)
      IC=IQC(N)
      JC=JQC(N)
      IF(JD.EQ.JC) THEN
            H2(IC,JC)=H2(ID,JD)
            IF(IC.GT.ID) THEN
              DUM(IC,JC)=1.0
            ELSE
              DUM(ID,JD)=1.0
            ENDIF
      ELSE
            H1(IC,JC)=H1(ID,JD)
            IF(JC.GT.JD) THEN
              DVM(IC,JC)=1.0
            ELSE
              DVM(ID,JD)=1.0
            ENDIF
      ENDIF
            H (IC,JC)=H (ID,JD)
            FSM(IC,JC)=FSM(ID,JD)
!
! get H1 and H2 for exterior boundary cells of river inflow
! and this is necessary for particle tracking routines    
! (hli, 02/12/02)
!
            H1(IC,JC)=H1(ID,JD)
            H2(IC,JC)=H2(ID,JD)
            ANG(IC,JC)=ANG(ID,JD)
!
!
  140    CONTINUE
!
      IF(HORZMIX.EQ.'CONSTANT  ') THEN
      DO 220 N=1,NUMQBC
      IC=IQC(N)
      JC=JQC(N)
      AAM2D(IC,JC)=0.0
      DO 220 K=1,KBM1
      AAM(IC,JC,K)=0.0
  220 CONTINUE
      ENDIF
!
 136  CONTINUE
!
!-----------------------------------------------------------------------
!-------- OFFSHORE INTAKE/OUTFALL(DIFFUSER) BOUNDARY -------------------
      IF (HYDTYPE.EQ.'EXTERNAL') THEN   ! QA 9/18/98
         NUMDBC=0
         GOTO 3700
      ENDIF
!
! add 'DIFFUSER IN LOOP' to the codes, referring to TAMPA BAY codes 
! in /power6/wrai0010/MODEL/CODES/TRACER (hli, 09/28/00)
!
! A.Jordi: Read diffuser1 boundary condition from netcdf file
      NT_T2D=0
      CALL NETCDF_READ_BC_DIFF1(NT_T2D,NT_T2D_MAX,STR1)
      WRITE(IUPRT,'(/1X,A80/)') STR1
      WRITE(IUPRT,'(I5)') NUMDBC1
! A.Jordi: end
      IF(NUMDBC1.EQ.0) GO TO 156
      IF(NUMDBC1.GT.DBCM) THEN
       WRITE(IUPRT,62) NUMDBC1,DBCM
       CALL SYSTEM ('rm gcm_temp*')
       STOP
      ENDIF
  62  FORMAT(//' NUMDBC1 (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' DBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      DO N=1,NUMDBC1
        WRITE(IUPRT,52) IDD(N),JDD(N),(VDDIST(N,K),K=1,KBM1)
      ENDDO
! CHECK THE VALIDITY OF THE DIFFUSER INPUT CONFIGURATION
      IF(RESTAR.EQ.'COLD START')THEN
       DO N=1,NUMDBC1
        IF(FSM(IDD(N),JDD(N)).NE.1.0)THEN  ! SEE IF SPECIFIED AT LAND CELL
         WRITE(IUPRT,3821)N,IDD(N),JDD(N)
         STOP
        END IF
          SUMVDDIST=0.0
         DO K=1,KBM1
          SUMVDDIST=SUMVDDIST+VDDIST(N,K)
         END DO

         IF(SUMVDDIST.NE.100.0)THEN
          WRITE(IUPRT,3822)N
         END IF

       END DO
      END IF
 3821 FORMAT('DIFFUSER INFLOW NO',I5,' SPECIFIED AT',2I5,
     .  'LAND CELL',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
 3822 FORMAT('WARNING: DIFFUSER INFLOW NO',I5,
     .  ' HAS VDDIST NOT EQUAL TO 100',//
     . 'PLEASE BEAWARE OF THIS ')
51    FORMAT(2I5,100F5.0)
52    FORMAT(2I5,/,100F5.1)
! A.Jordi: Read diffuser1 boundary condition from netcdf file
      DO M=1,NT_T2D_MAX
        CALL NETCDF_READ_BC_DIFF1(M,NT_T2D_MAX,STR1)
        WRITE(IUPRT,'(8F10.1)') T2D
        WRITE(IUPRT,'(8F10.1)') (DQDIFF(N,2),N=1,NUMDBC1)
        WRITE(IUPRT,'(8F10.1)') (DTDIFF(N,2),N=1,NUMDBC1)
        WRITE(IUPRT,'(8F10.1)') (DSDIFF(N,2),N=1,NUMDBC1)
        OTIME=T2D
      END DO
! A.Jordi: end
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     &  'NEED MORE DIFFUSER BC TO COMPLETE THE SIMULATIONS' 
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
!
 156  CONTINUE
!
****  NKIM (4/30/97) BEGINS FROM HERE
!-------- OFFSHORE INTAKE/OUTFALL(DIFFUSER) BOUNDARY IN LOOPS------------------
!
! A.Jordi: Read diffuser2 boundary condition from netcdf file
      NT_T2D2=0
      CALL NETCDF_READ_BC_DIFF2(NT_T2D2,NT_T2D2_MAX,STR1)
      WRITE(IUPRT,'(/1X,A80/)') STR1
      NUMDBC = NUMDBC1 + NUMDBC2
      WRITE(IUPRT,'(2I5)') NUMDBC2,NUMDBC
! A.Jordi: end
      IF (NUMDBC2.EQ.0) GOTO 3700
      IF ((NUMDBC1+NUMDBC2).GT.DBCM) THEN
        WRITE (IUPRT,63) NUMDBC1+NUMDBC2, DBCM
        CALL SYSTEM('rm gcm_temp9*')
        STOP
      END IF
  63  FORMAT(//' NUMDBC (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' DBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!
      DO N = NUMDBC1+1, NUMDBC1+NUMDBC2
        WRITE (IUPRT,52) IDD(N), JDD(N), (VDDIST(N,K),K = 1,KBM1)
      ENDDO    
! CHECK THE VALIDITY OF THE DIFFUSER-IN-LOOP INPUT CONFIGURATION
      DO N=1,NUMDBC1+1,NUMDBC
        IF(FSM(IDD(N),JDD(N)).NE.1.0)THEN  ! SEE IF SPECIFIED AT LAND CELL
         WRITE(IUPRT,3831)N-NUMDBC1,IDD(N),JDD(N)
         STOP
        END IF
          SUMVDDIST=0.0
         DO K=1,KBM1
          SUMVDDIST=SUMVDDIST+VDDIST(N,K)
         END DO

         IF(SUMVDDIST.NE.100.0)THEN
          WRITE(IUPRT,3832)N-NUMDBC1
          STOP
         END IF

       END DO
 3831 FORMAT('DIFFUSER-IN-LOOP INFLOW NO',I5,' SPECIFIED AT',2I5,
     .  'LAND CELL',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
 3832 FORMAT('DIFFUSER-IN-LOOP INFLOW NO',I5,
     .  ' HAS WRONG VDDIST: THE SUM SHOULD BE 100',//
     . 'PLEASE FIX THIS AND RESUBMIT THE RUN')
! A.Jordi: Read diffuser2 boundary condition from netcdf file
      DO M=1,NT_T2D2_MAX
        CALL NETCDF_READ_BC_DIFF2(M,NT_T2D2_MAX,STR1)
        WRITE(IUPRT,'(8F10.1)') T2D2
        WRITE(IUPRT,'(8F10.1)') (DQDIFF(N,2),N = NUMDBC1+1,
     .    NUMDBC1+NUMDBC2)
        WRITE(IUPRT,'(8F10.1)') (DTDIFF(N,2),N = NUMDBC1+1,
     .    NUMDBC1+NUMDBC2)
        WRITE(IUPRT,'(8F10.1)') (DSDIFF(N,2),N = NUMDBC1+1,
     .    NUMDBC1+NUMDBC2)
        OTIME=T2D2
      END DO
! A.Jordi: end
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     .  'NEED MORE DIFFUSER IN LOOP BC TO COMPLETE THE SIMULATIONS'
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
!
 3700 CONTINUE
!
      IF (TRACER.NE.'INCLUDE') GO TO 3140  ! GO TO SED/MET INPUT DECK
      READ(IUT401,11) (COM(I),I=1,80)
      WRITE(IUPRT,12) (COM(I),I=1,80)
      READ(IUT401,179) NUMDBCTR1
      WRITE (IUPRT,179)NUMDBCTR1
      IF(NUMDBCTR1.EQ.0) GO TO 3130   ! GO TO TRACER-IN-LOOP DECK
      IF (NUMDBCTR1.GT.DBCM) THEN
        WRITE (IUPRT,64) NUMDBCTR1, DBCM
        CALL SYSTEM('rm gcm_temp9*')
        STOP
      END IF

      DO 3436 N=1,NUMDBCTR1
        READ (IUT401,179)ITRDD(N),JTRDD(N)
        WRITE (IUPRT,179)ITRDD(N),JTRDD(N)
3436  CONTINUE
!
!   FIND PROPER DIFFUSER FLOW MATCH
!   2 OPTIONS: IF NUMDBCTR1=NUMDBC1, USER WANTS TO USE SAME ORDER OF DIFFUSERS
!               SPECIFIED IN DIFFUSER INFLOW CONDITION  (THIS OPTION MUST BE
!               USED IF THERE ARE MULTIPLE DIFFUSERS ARE ASSIGNED IN A
!               GRID CELL.  (ORIGINAL CODING)
!              IF NUMDBCTR1 .NE. NUMDBC1 THEN, PROGRAM ASSUMES THAT USER IS
!                ASSIGNING DIFFUSER TRACER INPUT SELECTIVELY.  (NKIM 5/10/99)
!
        IF(NUMDBCTR1.EQ.NUMDBC1) THEN

         DO 3434 ND = 1,NUMDBCTR1
          IDBC(ND)=ND
          IF(ITRDD(ND).NE.IDD(ND).OR.JTRDD(ND).NE.JDD(ND))THEN
            WRITE(IUPRT,2481) ITRDD(ND),JTRDD(ND)
            CALL SYSTEM ('rm gcm_temp*')
            STOP
          END IF
 3434    CONTINUE     

        ELSE IF (NUMDBCTR1.GT.NUMDBC1) THEN
            WRITE(IUPRT,2482) NUMDBCTR1,NUMDBC1
            CALL SYSTEM ('rm gcm_temp*')
            STOP
        ELSE
         WRITE(IUPRT,2484)
         DO 3435 ND = 1,NUMDBC1
            WRITE(IUPRT,2483)ITRDD(ND),JTRDD(ND),IDD(ND),JDD(ND)
 3435    CONTINUE     
        END IF
 3438    CONTINUE     

 2479 FORMAT(' ONE OF DIFFUSER-LOOP TRACER INPUT DOES NOT MATCH WITH
     .DIFFUSER LOOP  LOCATIONS. ITS I,J ARE',2I5)
 2481 FORMAT(' ONE OF DIFFUSER TRACER INPUT DOES NOT MATCH WITH
     .DIFFUSER INFLOW LOCATIONS. ITS I,J ARE',2I5)
 2482 FORMAT(' NUMBER OF DIFFUSER (TRACER) =',I5,' IS GREATER THAN
     .NUMBER OF DIFFUSER (FLOW)=',I5)
 2483 FORMAT('ITRDD JTRDD IDD JDD =',4I5)
 2484 FORMAT(' CAUTION NUMBER OF DIFFUSER (TRACER) =',I5,' IS LESS THAN
     .NUMBER OF DIFFUSER (FLOW)=',I5)


      DO N=1,NUMDBC
        CDIFF1(N)=0.0
      END DO
      OPEN(IUT98,FILE='gcm_temp98')
!
      DO 3437 M=1,100000
        READ(IUT401,77,ERR=3128) TIME
        READ(IUT401,77)  (TEMPD(N),N=1,NUMDBCTR1)
        DO N=1,NUMDBCTR1
          CDIFF1(IDBC(N))=TEMPD(N)
        END DO
!
        WRITE(IUPRT,77) TIME
        OTIME=TIME
        WRITE (IUPRT,77)  (CDIFF1(N),N=1,NUMDBCTR1)
!
        WRITE (IUT98,78) TIME
        WRITE (IUT98,78)  (CDIFF1(N),N=1,NUMDBCTR1)
!
 3437 CONTINUE

 3128 BACKSPACE IUT401
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     .  'NEED MORE DIFFUSER TRACER BC TO COMPLETE THE RUN'
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
***********  NKIM (4/30/97) ENDS HERE
 3130 CONTINUE
!
! FOR TRACER LOADINGS IN DIFFUSER-IN-LOOP  (NKIM 07/25/00)
!-------- OFFSHORE INTAKE/OUTFALL(DIFFUSER IN LOOPS) TRACER BOUNDARY------------------
! CAUTION !!! DIFFUSER-IN-LOOP TRACER INPUTS SHOULD MATCH WITH THE
!             DIFFUSER-IN-LOOP THAT Is ID =ITRDD and JD=JTRDD
!             THEREFORE NO NEED TO READ ITRDD and JTRDD
!
      READ (IUT401,11) (COM(I),I = 1,80)
      WRITE(IUPRT,12) (COM(I),I = 1,80)
      READ (IUT401,8) NUMDBCTR2
      IF(NUMDBCTR2.EQ.0) GO TO 3140   ! EXIT TRACER-IN-LOOP DECK ! NG 08052011 
!
      IF(NUMDBCTR2.NE.NUMDBC2)THEN
        WRITE(IUPRT,*)'THE NUMBER OF DIFFUSER-IN-LOOP TRACER INPUT
     .  IS NOT EQUAL TO NUMDBC2'
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      END IF
!
      NUMDBCTR = NUMDBCTR1 + NUMDBCTR2
      WRITE (IUPRT,8200) NUMDBCTR2,NUMDBCTR
 8200 FORMAT('NUMBER OF DIFFUSER IN LOOP TRACER LOADING IS',
     . I6,/'TOTAL NUMBER OF DIFFUSER TRACER LOADING IS',I6)
      IF (NUMDBCTR2.EQ.0) GOTO 3140

      OPEN(IUT99,FILE='gcm_temp99')

      DO 3439 N=NUMDBCTR1+1,NUMDBCTR
        READ (IUT401,179)ITRDD(N),JTRDD(N)
        WRITE (IUPRT,179)ITRDD(N),JTRDD(N)
3439  CONTINUE

      IF (NUMDBCTR2.EQ.NUMDBC2) THEN

         DO 3432 ND = NUMDBCTR1+1,NUMDBCTR
          IDBC(ND)=ND
          IF(ITRDD(ND).NE.IDD(ND).OR.JTRDD(ND).NE.JDD(ND))THEN
            WRITE(IUPRT,2479) ITRDD(ND),JTRDD(ND)
            CALL SYSTEM ('rm gcm_temp*')
            STOP
          END IF
 3432    CONTINUE

      ELSE IF (NUMDBCTR2.GT.NUMDBC2) THEN
            WRITE(IUPRT,64) NUMDBCTR,DBCM,NUMDBCTR2,NUMDBC2
            CALL SYSTEM ('rm gcm_temp*')
            STOP
      ELSE
         WRITE(IUPRT,2484)
         DO 3440 ND = 1,NUMDBC1
            WRITE(IUPRT,2483)ITRDD(ND),JTRDD(ND),IDD(ND),JDD(ND)
 3440    CONTINUE
      END IF
  64  FORMAT(//' NUMDBCTR (=',I4,') MUST BE LESS THAN OR EQUAL TO'/
     .         ' DBCM   (=',I4,') SPECIFIED IN comdeck'/
     .         ' PLEASE FIX AND RESUBMIT'//
     .         ' OR OR OR'/
     .         ' NUMDBCTR2 (=',I4,') MUST BE EQUAL TO'/
     .         ' NUMDBC2   (=',I4,') SPECIFIED IN run_data'/
     .         ' PLEASE FIX AND RESUBMIT'//)
!


      DO N=NUMDBCTR1+1,NUMDBCTR
        CDIFF1(N)=0.0
      ENDDO

      DO M = 1, 100000
        READ (IUT401,77,ERR=3133) TIME
        WRITE (IUPRT,77) TIME
        OTIME=TIME
        READ (IUT401,77) (CDIFF1(N),N = NUMDBCTR1+1,NUMDBCTR1+NUMDBCTR2)
        WRITE(IUPRT,77) (CDIFF1(N),N = NUMDBCTR1+1,NUMDBCTR1+NUMDBCTR2)
        WRITE (IUT99,78) TIME
        WRITE (IUT99,78) (CDIFF1(N),N = NUMDBCTR1+1,NUMDBCTR1+NUMDBCTR2)
      ENDDO     
!
 3133 BACKSPACE IUT401
      IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     .  'NEED MORE DIFFUSER IN LOOP BC TO COMPLETE THE SIMULATIONS'
        CALL SYSTEM ('rm gcm_temp*')
        STOP
      ENDIF
!
 3140 CONTINUE
!
!  end of adding DIFFUSER IN LOOP (hli, 09/28/00)
! 
!
!*******************************************************************
!
!  SEDIMENT TRANSPORT By Quamrul 5/13/99
!
      IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT702=702
          OPEN (IUT702,FILE='gcm_temp702')
!
          READ(IUT402,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT402,24)NUMDBCSE
          WRITE (IUPRT,24)NUMDBCSE
!qa          write (*,*)'NUMDBCSE=',NUMDBCSE
!
          IF (NUMDBCSE.EQ.0) GOTO 157
!
          IF (NUMDBCSE.GT.DBCM) THEN
           WRITE(IUPRT,5402)NUMDBCSE,DBCM
           CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
!
          DO 7529 N=1,NUMDBCSE
            READ(IUT402,51)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
 7529     CONTINUE
!
          DO 7560 M=1,100000
            READ(IUT402,77,END=7570) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT402,771)  (CDIFF(1,N),N=1,NUMDBCSE)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCSE)
!
!  Quamrul/Pravi May 25, 1999
!  Diffuser Loadings are in Kg/day. We need to compute sediment concentration
!  in g/cm3. Therefore Loading should be divided by a factor (86.4*1.0e6)
!
            DO 7580 NN=1,NUMDBCSE
              CDIFF(1,NN)=CDIFF(1,NN)/(86.4*1.0e06)
 7580       CONTINUE
!
            WRITE (IUT702,78)TIME
            WRITE (IUT702,78)  (CDIFF(1,N),N=1,NUMDBCSE)
!
 7560     CONTINUE
 7570     CONTINUE
        ENDIF

 157    IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT703=703
          OPEN (IUT703,FILE='gcm_temp703')
!
          READ(IUT403,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT403,24)NUMDBCSE
          WRITE (IUPRT,24)NUMDBCSE
!
          IF (NUMDBCSE.EQ.0) GOTO 158
!
          IF (NUMDBCSE.GT.DBCM) THEN
           WRITE(IUPRT,5402)NUMDBCSE,DBCM
           CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
!
          DO 7559 N=1,NUMDBCSE
            READ(IUT403,51)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDSE(N),JDDSE(N),(VDDISTSE(N,K),K=1,KBM1)
 7559     CONTINUE
!
!
          DO 7561 M=1,100000
            READ(IUT403,77,END=7571) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT403,771)  (CDIFF(1,N),N=1,NUMDBCSE)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCSE)
!
!
!  Quamrul/Pravi May 25, 1999
!  Diffuser Loadings are in Kg/day. We need to compute sediment concentration
!  in g/cm3. Therefore Loading should be divided by a factor (86.4*1.0e6)
!
            DO 7581 NN=1,NUMDBCSE
              CDIFF(1,NN)=CDIFF(1,NN)/(86.4*1.0e06)
 7581       CONTINUE
!
            WRITE (IUT703,78)TIME
            WRITE (IUT703,78)  (CDIFF(1,N),N=1,NUMDBCSE)
!
 7561     CONTINUE
 7571     CONTINUE
        ENDIF
      ENDIF

!*******************************************************************
!
!  CHEMICAL TRANSPORT
!
  158 IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT704=704
          OPEN (IUT704,FILE='gcm_temp704')
!
          READ(IUT404,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT404,24)NUMDBCCH
          WRITE (IUPRT,24)NUMDBCCH
!
          IF (NUMDBCCH.EQ.0) GOTO 1588
          IF (NUMDBCCH.GT.DBCM) THEN
           WRITE(IUPRT,7402)NUMDBCCH,DBCM
           CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
!
          DO 7589 N=1,NUMDBCCH
            READ (IUT404,51)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
 7589     CONTINUE
!
!
          DO 7590 M=1,100000
            READ(IUT404,77,END=7600) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT404,771)  (CDIFF(1,N),N=1,NUMDBCCH)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCCH)
!
!  Quamrul/Pravi May 25, 1999
!  Diffuser Loadings are in Kg/day. We need to compute chemical concentration
!  in ug/cm3. Therefore Loading should be divided by a factor (86.4)
!
            DO 7610 NN=1,NUMDBCCH
              CDIFF(1,NN)=CDIFF(1,NN)/(86.4)
 7610       CONTINUE
!
            WRITE (IUT704,78)TIME
            WRITE (IUT704,78)  (CDIFF(1,N),N=1,NUMDBCCH)
!
 7590     CONTINUE
 7600     CONTINUE
        ENDIF
!
 1588   IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          IUT705=705
          OPEN (IUT705,FILE='gcm_temp705')
!
          READ(IUT405,11)  (COM(I),I=1,80)
          WRITE(IUPRT,12) (COM(I),I=1,80)
          READ (IUT405,24)NUMDBCCH
          WRITE (IUPRT,24)NUMDBCCH
!
          IF (NUMDBCCH.EQ.0) GOTO 1589
          IF (NUMDBCCH.GT.DBCM) THEN
           WRITE(IUPRT,7402)NUMDBCCH,DBCM
           CALL SYSTEM ('rm gcm_temp*')
           STOP
          ENDIF
!
          DO 7599 N=1,NUMDBCCH
            READ (IUT405,51)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
            WRITE (IUPRT,52)IDDCH(N),JDDCH(N),(VDDISTCH(N,K),K=1,KBM1)
 7599     CONTINUE
!
!
          DO 7620 M=1,100000
            READ(IUT405,77,END=7630) TIME
            WRITE(IUPRT,77) TIME
            READ(IUT405,771)  (CDIFF(1,N),N=1,NUMDBCCH)
            WRITE (IUPRT,771)  (CDIFF(1,N),N=1,NUMDBCCH)
!
!  Quamrul/Pravi May 25, 1999
!  Diffuser Loadings are in Kg/day. We need to compute chemical concentration
!  in ug/cm3. Therefore Loading should be divided by a factor (86.4)
!
            DO 7640 NN=1,NUMDBCCH
              CDIFF(1,NN)=CDIFF(1,NN)/86.4
 7640       CONTINUE
!
            WRITE (IUT705,78)TIME
            WRITE (IUT705,78)  (CDIFF(1,N),N=1,NUMDBCCH)
!
 7620     CONTINUE
 7630     CONTINUE
        ENDIF
      ENDIF

!
!*******************************************************************
!
!  POINT SOURCE LOADS:  CONSERVATIVE TRACER
!
 1589  IF (TRACER.EQ.'INCLUDE') THEN
        OPEN (IUT701,FILE='gcm_temp701')
!
        READ(IUT401,11)  (COM(I),I=1,80)
        WRITE(IUPRT,12) (COM(I),I=1,80)
!
        READ (IUT401,1120)NUMPSTR,OPTPSTR 
        WRITE (IUPRT,1120)NUMPSTR,OPTPSTR 
 1120   FORMAT (I5,6X,A4)
!
        IF (NUMPSTR.EQ.0) GOTO 3129
!qa        IF (NUMPSTR.EQ.0) GOTO 1123
!
        DO 1121 N=1,NUMPSTR
          READ (IUT401,179)IPSTR(N),JPSTR(N),KPSTR(N)
          WRITE (IUPRT,179)IPSTR(N),JPSTR(N),KPSTR(N)
 1121   CONTINUE
!
        DO 1122 M=1,100000
          READ(IUT401,77,END=1123) TIME
          READ(IUT401,77) (CDIS1(N),N=1,NUMPSTR)
!
          WRITE(IUPRT,77) TIME
          OTIME=TIME
          WRITE (IUPRT,77) (CDIS1(N),N=1,NUMPSTR)
!
          WRITE (IUT701,78)TIME
          WRITE (IUT701,78) (CDIS1(N),N=1,NUMPSTR)
!
 1122   CONTINUE
!
 1123   CONTINUE 
       IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)
     . 'NEED MORE POINT SOURCE TRACER IN LOOP BC TO COMPLETE THE RUN' 
         CALL SYSTEM ('rm gcm_temp*')
         STOP
       ENDIF
!
      ENDIF
!
3129  CONTINUE
!
!-----------------------------------------------------------------------
!-------- METEOROLOGICAL BOUNDARY CONDITIONS ---------------------------
      IF (HYDTYPE.EQ.'EXTERNAL') GOTO 1341   ! QA 9/18/98
!
! Open a temporary file for wave computations (hli, 03/22/04).
!NG IUT193 This is used to store wind speeds.
!
!NG02222007      IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008 BELOW NOT NEEDED ANY MORE. USE WSX2D, WSY2D FOR WINDS FROM NOW ON.
!NG04092008      IF(WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008        OPEN(IUT193,FILE='gcm_temp193',FORM='unformatted')
!NG04092008      ENDIF
!
! A.Jordi: Read meteorological boundary condition from netcdf file
      NT_T2M=1
      CALL NETCDF_READ_BC_MET(NT_T2M,NT_T2M_MAX,STR1,WDS,WDD)
      WRITE(IUPRT,'(/1X,A80/)') STR1
!svv080326 added 2d heatflux
      WRITE(IUPRT,7734) OPTMBC,ALAT,ALON,TR,WNDSH,REFL,OPTEXTC,
     .  CALCCLD,CALCSWOBS
! A.Jordi: end   
 131  FORMAT(A10,5F10.2,6X,A4,2F5.0)
!  MODIFIED BY QUAMRUL QA ON 5/3/96 for fraction of SWRAD, TR (0.-1.),
!  Absorbed in surface layer of water colum, rest (1-TR) is penetrated
!  Modified By Quamrul QA to add Wind Sheltering Coefficient WNDSH
!
7734  FORMAT('METEOROLOGICAL FORCING OPTION IS ',A10,
     +     /'Latitude= ', f10.4,'     Longitude= ', f10.4/
     +     5x,'Fraction of SWRad Absorbed in Surface = ', F10.3/
     +     5x,'Wind Sheltering Coeffient Used =', F10.3/
     +     5x,'Reflection Coeffient Used =', F10.3/
     +     5x,'EXTINCTION COEFFICIENT OPTION IS ',A10/
     +     5x,'Cloud Cover Request =', F10.3/
     +     5x,'Shortwave Radiation Request =', F10.3)
!
!  MODIFIED BY QUAMRUL QA ON 5/6/96 *************************************
!  Add three new options (SYNOPANB, SYNOPLNP and SYNOPRNM) 
!  (hli,04/24/02)
!
!svv080306 added new options 2DMETANB 2DMETLP 2DMETRNM
      If (OPTMBC(1:8).NE.'AVERAGED'.AND.OPTMBC(1:8).NE.'LANDPFLX'.
     *AND.OPTMBC(1:8).NE.'AANDBFLX'.AND.OPTMBC(1:8).NE.'RANDMFLX'.
     *AND.OPTMBC(1:8).NE.'SYNOPTIC'.AND.OPTMBC(1:8).NE.'SYNOPANB'.
     *AND.OPTMBC(1:8).NE.'SYNOPLNP'.AND.OPTMBC(1:8).NE.'SYNOPRNM'.
     *AND.OPTMBC(1:8).NE.'2DMETANB'.AND.OPTMBC(1:8).NE.'2DMETLNP'.
     *AND.OPTMBC(1:8).NE.'2DMETRNM') then
        Write (IUPRT,5401) OPTMBC
        Call SYSTEM('rm gcm_temp9*')
        Stop
      End If

!NG2014 2D Precipitation and Evaporation
! TODO: synop_pet file (IUT104)
      INQUIRE(FILE='synop_pet',EXIST=CREADPET) ! CNG2014 For 2D precipitation and evaporation
      IF(CREADPET) OPEN(IUT104,FILE='synop_pet',FORM='unformatted') ! CNG2014
!NG2014
!
!-------- AVERAGED METEOROLOGICAL BOUNDARY CONDITIONS ------------------
!-------- PRECIPITATION(m/yr), EVAPORATION(m/yr) & HEAT FLUX(w/m2) -----
!-------- WIND FROM DIRECTION WIND IS BLOWING --------------------------
      IF(OPTMBC(1:8).EQ.'AVERAGED') THEN
! A.Jordi: Read meteorological boundary condition from netcdf file
        DO M = 1,NT_T2M_MAX
          CALL NETCDF_READ_BC_MET(M,M_MAX,STR1,WDS,WDD)
          WRITE(IUPRT,'(8F10.1)') T2M
          OTIME=T2M
          WRITE(IUPRT,'(8F10.1)') WDS,WDD,HFLUX,QPREC,QEVAP
! A.Jordi: end
       QPREC=QPREC/(86400.*365.)
       QEVAP=QEVAP/(86400.*365.)
!-------- WIND STRESS --------------------------------------------------
       WDD=180.+WDD
       WDD=AMOD(WDD,360.)
       VWIND=WDS*COS(6.283185*WDD/360.)
       UWIND=WDS*SIN(6.283185*WDD/360.)
!
! Save wind speed to IUT193 for wave computations (hli, 03/22/04).
!
!NG02222007      IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008 BELOW NOT NEEDED ANY MORE. READ UWIND AND VWIND FROM IUT93
!NG04092008       IF(WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008         DO I=2,IM1
!NG04092008           DO J=2,JM1
!NG04092008             WSX2D(I,J)=UWIND ! NG04092008 TX2D(I,J)=UWIND
!NG04092008             WSY2D(I,J)=VWIND ! NG04092008 TY2D(I,J)=VWIND
!NG04092008           ENDDO
!NG04092008         ENDDO
!NG04092008         WRITE(IUT193) TIME
!NG04092008         WRITE(IUT193) ((WSX2D(I,J),WSY2D(I,J),I=1,IM),J=1,JM)
!NG04092008       ENDIF
!   
       IF(WAVEDYN.EQ.'DONONLY '.OR.WAVEDYN.EQ.'NEGLECT') THEN ! USE LARGE AND POND
        CD=1.2E-3
        IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
        IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
       ELSE ! SIMULATION INCLUDES WAVES; USE TAYLOR-YELAND
        Z0air=3.5E-5 + 4*SIG(I,J)*1200*(4*SIG(I,J)*WN(I,J)/TPI)**4.5
        CD=0.16*ALOG(10/Z0air)**(-2)
        IF(CD.GE.3E-3) CD=3E-3
       ENDIF
       TX=1.2*CD*UWIND*WDS
       TY=1.2*CD*VWIND*WDS
       IF (WDS.GE.WDSMAX) then ! Limit internal wind stress due to local domain
        TX=1.2*CD*UWIND*WDSMAX/WDS*WDSMAX
        TY=1.2*CD*VWIND*WDSMAX/WDS*WDSMAX
       ENDIF
       WRITE(IUT93,78) TIME
       WRITE(IUT93,78) HFLUX,TX,TY,UWIND,VWIND,WDS,WDD,QPREC,QEVAP
       END DO
       IF(OTIME/24.LE.(EDAY-SDAY))THEN
         WRITE(IUPRT,*)'NEED MORE MET DATA TO COMPLETE THE RUN' 
         CALL SYSTEM ('rm gcm_temp*')
         STOP
       ENDIF
!
!hli Move the "SYNOPTIC" input to wave section (03/18/04)
!
!------ HEATFLUX CALCULATION --------------------------------------
!START, MODIFIED FOR HEATFLUX CALCULATION BY jeff JI, 8/30/94
!       ELSE IF (OPTMBC(1:8).EQ.'HEATFLUX') THEN
!
!  MODIFIED BY QUAMRUL QA ON 5/6/96 ********************************
!
      ELSE IF (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     *         OPTMBC(1:8).EQ.'AANDBFLX'.OR.
     *         OPTMBC(1:8).EQ.'RANDMFLX'.OR.
     .         OPTMBC(1:8).EQ.'SYNOPANB'.OR.
     .         OPTMBC(1:8).EQ.'SYNOPLNP'.OR.
     .         OPTMBC(1:8).EQ.'SYNOPRNM') Then
        IF(ALAT.LE.-999.0) THEN
!
!NG THIS IS A WEIRD LEGACY OPTION, ALAT IS ALSO NEEDED...
         ALON=-ALON  ! TO FOLLOW THE CONVENTION IN NCLD.F
        ENDIF
!
        AIRLOW=-999.9 ! DO ALLOW NEGATIVE TEMPERATURES, ICE, NG02152007
        AIRHIGH=999.9
        WRITE(IUPRT,7735) AIRLOW,AIRHIGH
7735  FORMAT(' THE LOWEST AIR TEMPERATURE ALLOWED: ',f10.2/,
     &       ' THE HIGHEST AIR TEMPERATURE ALLOWED: ',f10.2)
! A.Jordi: Read meteorological boundary condition from netcdf file
        DO M = 1,NT_T2M_MAX
          CALL NETCDF_READ_BC_MET(M,M_MAX,STR1,WDS,WDD)
          WRITE(IUPRT,'(8F10.1)') T2M
          OTIME=T2M
! A.Jordi: end
!
!       MODIFIED BY Quamrul QA 5/6/96 **********************************
!       Adding Extinction Coefficient **********************************
!
!---       
          IF(AIRTMP.LT.AIRLOW) AIRTMP=AIRLOW
          IF(AIRTMP.GT.AIRHIGH) AIRTMP=AIRHIGH
          WRITE (IUPRT,7750) WDS,WDD,SWOBS,AIRTMP,RELHUM,BAROP,
     &      CLD,EXTC,QPREC,QEVAP
! CHANGE EVAP AND PRECIP UNITS FROM m/year TO m/s.
          QPREC = QPREC / (86400.*365.)
          QEVAP = QEVAP / (86400.*365.)
!-------- WIND STRESS (large & pond, newtons/m2)
          WDD = 180. + WDD
          WDD = AMOD(WDD,360.)
          VWIND = WDS * COS(6.283185*WDD/360.)
          UWIND = WDS * SIN(6.283185*WDD/360.)
!
! Save wind speed to IUT193 for wave computations (hli, 03/22/04).
!
!NG02222007      IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008 BELOW NOT NEEDED ANY MORE. READ UWIND AND VWIND FROM IUT93
!NG04092008          IF(WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008            IF(OPTMBC(1:8).EQ.'LANDPFLX'.OR.
!NG04092008     *         OPTMBC(1:8).EQ.'AANDBFLX'.OR.
!NG04092008     *         OPTMBC(1:8).EQ.'RANDMFLX') THEN
!NG04092008              DO I=2,IM1
!NG04092008                DO J=2,JM1
!NG04092008                  WSX2D(I,J)=UWIND ! NG04092008 TX2D(I,J)=UWIND
!NG04092008                  WSY2D(I,J)=VWIND ! NG04092008 TY2D(I,J)=VWIND
!NG04092008                ENDDO
!NG04092008              ENDDO
!NG04092008              WRITE(IUT193) TIME
!NG04092008              WRITE(IUT193) ((WSX2D(I,J),WSY2D(I,J),I=1,IM),J=1,JM)
!NG04092008            ENDIF
!NG04092008          ENDIF
!  

          RHOA=SADENS(BAROP,AIRTMP,RELHUM) ! NG06042012 
          IF(WAVEDYN.EQ.'DONONLY '.OR.WAVEDYN.EQ.'NEGLECT') THEN ! USE LARGE AND POND
           CD=1.2E-3
           IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
           IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
          ELSE ! SIMULATION INCLUDES WAVES; USE TAYLOR-YELAND
           Z0air=3.5E-5 + 4*SIG(I,J)*1200*(4*SIG(I,J)*WN(I,J)/TPI)**4.5
           CD=0.16*ALOG(10/Z0air)**(-2)
           IF(CD.GE.3E-3) CD=3E-3
          ENDIF
          TX = RHOA * CD * UWIND * WDS
          TY = RHOA * CD * VWIND * WDS
          IF (WDS.GE.WDSMAX) then ! Limit internal wind stress due to local domain
           TX=RHOA*CD*UWIND*WDSMAX/WDS*WDSMAX
           TY=RHOA*CD*VWIND*WDSMAX/WDS*WDSMAX
          ENDIF
          WRITE (IUT93,78) TIME
          WRITE (IUT93,7850) AIRTMP, RELHUM, BAROP,
     &     SWOBS, TX, TY, UWIND, VWIND, CLD, EXTC, QPREC, QEVAP
!
        END DO
        IF(OTIME/24.LE.(EDAY-SDAY))THEN
          WRITE(IUPRT,*)'NEED MORE MET DATA TO COMPLETE THE RUN' 
          CALL SYSTEM ('rm gcm_temp*')
          STOP
        ENDIF
! END,  MODIFIED FOR HEARFLUX CALCULATION BY jeff ji, 8/30/94
      ENDIF
!-------- Spatial and temporal varying extinction coefficient ----------
!-------- (hli, 10/02/03) ---------------------------------------------
!svv080326: added 2dmet options
      IF(OPTMBC(1:8).NE.'2DMETLNP'
     .  .AND.OPTMBC(1:8).NE.'2DMETANB'
     .  .AND.OPTMBC(1:8).NE.'2DMETRNM'
     .  .AND.OPTEXTC.EQ.'VARI') THEN

        OPEN(IUT801,FILE='extc2d.inp',STATUS='old',FORM='unformatted')
        OPEN(IUT194,FILE='gcm_temp194',FORM='unformatted')
*
        DO M=1,100000
          READ(IUT801,END=1175) TIME
          READ(IUT801)((EXTC1(I,J),I=1,IM),J=1,JM) 
          OTIME=TIME
          WRITE(IUT194)TIME
          WRITE(IUT194)((EXTC1(I,J),I=1,IM),J=1,JM)
        ENDDO
 1175   CONTINUE
        IF(OTIME/24.LE.(EDAY-SDAY))THEN
          WRITE(IUPRT,*)
     .      'NEED MORE EXTC COEFF DATA TO COMPLETE THE RUN' 
          CALL SYSTEM ('rm gcm_temp*')
          STOP
        ENDIF
      ENDIF

      IF (OPTMBC(1:8).EQ.'SYNOPANB'.OR.
     &    OPTMBC(1:8).EQ.'SYNOPLNP'.OR.
     &    OPTMBC(1:8).EQ.'SYNOPRNM'.OR.
     &    OPTMBC(1:8).EQ.'SYNOPTIC') THEN
!-------- FROM U AND V WIND COMPONENTS
!---------------------------------
       IF(OPTMBC(1:8).EQ.'SYNOPTIC') THEN
         RHOA=1.2 ! NG 2014 TO USE FOR SYNOPTIC, OTHERWISE USES RHOA FROM ABOVE
         OPEN(IUT191,FILE='synop_hflx',FORM='unformatted',STATUS='old')
         OPEN(IUT192,FILE='gcm_temp192',FORM='unformatted')
*
         DO M=1,100000
           READ(IUT191,END=175) TIME
           READ(IUT191) ((SHFLX(I,J),I=1,IM),J=1,JM)
           OTIME=TIME
           WRITE(IUT192)TIME
           WRITE(IUT192)((SHFLX(I,J),I=1,IM),J=1,JM)
         ENDDO
 175     CONTINUE
       ENDIF
!-------- FROM U AND V WIND COMPONENTS
!---------------------------------
       OPEN(IUUAV,FILE='synop_wind',FORM='unformatted')
       OPEN(IUT95,FILE='gcm_temp95',FORM='unformatted')
!hli(11/26/03)
!
       DO 281 M=1,100000
       READ(IUUAV,END=285) TIME
!NG04092008       READ(IUUAV) ((TX2D(I,J),TY2D(I,J),PATM(I,J),I=1,IM),J=1,JM)
       READ(IUUAV) ((WSX2D(I,J),WSY2D(I,J),PATM(I,J),I=1,IM),J=1,JM)
!hli
!NG02222007      IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008 BELOW NOT NEEDED ANY MORE. READ UWIND AND VWIND FROM IUT95
!NG04092008       IF(WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG04092008          WRITE(IUT193) TIME
!NG04092008          WRITE(IUT193) ((WSX2D(I,J),WSY2D(I,J),I=1,IM),J=1,JM)
!NG04092008          print *,time,WSX2D(65,138),WSY2D(65,138) ! Print Castle Point wind
!NG04092008       ENDIF
!NG       print *,time,WSX2D(65,138),WSY2D(65,138) ! Print Castle Point wind
!hli
       DO 282 J=1,JM
       DO 282 I=1,IM
       WDS=SQRT(WSX2D(I,J)**2+WSY2D(I,J)**2)
       IF(WAVEDYN.EQ.'DONONLY '.OR.WAVEDYN.EQ.'NEGLECT') THEN ! USE LARGE AND POND
        CD=1.2E-3
        IF(WDS.GE.11.) CD=(0.49+0.065*WDS)*1.E-3
        IF(WDS.GE.25.) CD=(0.49+0.065*25.)*1.E-3
       ELSE ! SIMULATION INCLUDES WAVES; USE TAYLOR-YELAND
        Z0air=3.5E-5 + 4*SIG(I,J)*1200*(4*SIG(I,J)*WN(I,J)/TPI)**4.5
        CD=0.16*ALOG(10/Z0air)**(-2)
        IF(CD.GE.3E-3) CD=3E-3
       ENDIF
       TX2D(I,J)=RHOA*CD*WSX2D(I,J)*WDS
       TY2D(I,J)=RHOA*CD*WSY2D(I,J)*WDS
       IF (WDS.GE.WDSMAX) then ! Limit internal wind stress due to local domain
        TX2D(I,J)=RHOA*CD*WSX2D(I,J)*WDSMAX/WDS*WDSMAX
        TY2D(I,J)=RHOA*CD*WSY2D(I,J)*WDSMAX/WDS*WDSMAX
       ENDIF
       PATM(I,J) = PATM(I,J)*1.0E02 ! from mbar to N/m2
282    CONTINUE

       WRITE(IUT95) TIME
!NG04032008 281    WRITE(IUT95) ((TX2D(I,J),TY2D(I,J),PATM(I,J),I=1,IM),J=1,JM)
281    WRITE(IUT95) ((TX2D(I,J),TY2D(I,J),PATM(I,J),
     +                WSX2D(I,J),WSY2D(I,J),I=1,IM),J=1,JM)
285    CONTINUE
       IF(OTIME/24.LE.(EDAY-SDAY))THEN
        WRITE(IUPRT,*)'NEED MORE SYNOPTIC WIND DATA TO COMPLETE THE RUN'
        CALL SYSTEM ('rm gcm_temp*')
        STOP
       ENDIF
       CLOSE(IUUAV)
      END IF
!
!Based on sfan's (SIT) input, specify wave parameters on open boundaries
!(hli, 03/19/04).
!     WHBDRY: SIGNIFICANT WAVE HEIGHT OF WIND WAVES AND SWELL (meters)
!     WPBDRY: MEAN WAVE PERIOD (seconds)
!     WDBDRY: MEAN WAVE DIRECTION. The units are degrees CW true north
!             in oceanographic convention (pointing towards, not from)
!     WHBFLX: MOMENTUM FLUX
!     WDBFLX: MOMENTUM FLUX DIRECTION, W.R.T. ZETA1-DIRECTION AS
!             DONELAN AND SCHWAB'S CONVENTION. Wave is moving towards the
!             direction not from, counter clock is positive.
!
!NG02222007      IF(WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
!NG05172010 FOR LAKES   IF(WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL')THEN
      IF((WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL')
     .   .AND.NUMEBC.GT.0) THEN
!NG10302013 INQUIRE FOR wave.inp. If boundary conditions don't exist,
!NG10302013 assign from internal solution. (No-gradient boundary condition).
       INQUIRE(FILE='wave.inp',EXIST=CREADWBC)
       IF(CREADWBC) THEN
        OPEN(IUT802,FILE='wave.inp',STATUS='old',FORM='unformatted')
        OPEN(IUT195,FILE='gcm_temp195',FORM='unformatted')
!
        READ(IUT802)WAVE_HEADER        
        WRITE(IUPRT,122)WAVE_HEADER    
        DO M=1,100000
          READ(IUT802,END=286) TIME
          write(iuprt,*) 'wave boundary condition'
          WRITE(IUPRT,77) TIME
!          write(*,*) time,numebc
          READ(IUT802)(WHBDRY(N),N=1,NUMEBC)
          READ(IUT802)(WPBDRY(N),N=1,NUMEBC)
          READ(IUT802)(WDBDRY(N),N=1,NUMEBC)
          WRITE(IUPRT,'(8F10.3)')(WHBDRY(N),N=1,NUMEBC)
          WRITE(IUPRT,'(8F10.3)')(WPBDRY(N),N=1,NUMEBC)
          WRITE(IUPRT,'(8F10.3)')(WDBDRY(N),N=1,NUMEBC)
          WRITE(IUT195)TIME
          OTIME=TIME
          DO N=1,NUMEBC
            WHBFLX(N)=WHBDRY(N)/4.0
!NG06042012 WAVE DIRECTIONS ARE IN OCEANOGRAPHIC CONVENTION THROUGHOUT
            WDBFLX(N)=90.0-WDBDRY(N)-ANG(IETA(N),JETA(N))*180/
     &                3.141593
            WDBFLX(N)=AMOD(WDBFLX(N),360.0)
            IF(WDBFLX(N).LT.0.0) WDBFLX(N)=360+WDBFLX(N)
            WPBFLX(N)=WPBDRY(N)
          ENDDO
          WRITE(IUT195)(WHBFLX(N),N=1,NUMEBC)
          WRITE(IUT195)(WPBFLX(N),N=1,NUMEBC)
          WRITE(IUT195)(WDBFLX(N),N=1,NUMEBC)
        ENDDO
286     CONTINUE
        IF(OTIME/24.LE.(EDAY-SDAY))THEN
          WRITE(IUPRT,*)'NEED MORE WAVE DATA TO COMPLETE THE RUN'
          CALL SYSTEM ('rm gcm_temp*')
          STOP
        ENDIF
        CLOSE(IUT802)
       ELSE
        WRITE(IUPRT,*)' ** No wave.inp found. Using open wave BC. **'
       ENDIF
      ENDIF
!svv080326
      IF (OPTMBC(1:8).EQ.'2DMETANB'.OR.
     &    OPTMBC(1:8).EQ.'2DMETLNP'.OR.
     &    OPTMBC(1:8).EQ.'2DMETRNM') THEN
!-------- 2D WIND AND HEAT FLUX PARAMETERS
!---------------------------------
       OPTEXTC='VARI'
      END IF
!End (hli, 03/19/04)
!
 77    FORMAT(8F10.1)
 771   FORMAT(5E14.7)
 7750  FORMAT(10F10.5)
 78    FORMAT(8E14.7)
 7850  FORMAT(12E14.7)

!NG2011 SPACE-VARIABLE ICE DRAG OR SPACE-TIME VARIABLE ICE THICKNESS AND ICE CONCENTRATION. START...
!NG02182011   SPACE-VARIABLE ICE FRICTION, BUT NOT TIME-VARIABLE
      ICEFRICTION=-1
      OPEN (IUT196,FILE='ifric2d.inp',STATUS='old',ERR=39) 
      ICEFRICTION=2 ! Space-Variable (2D) ice friction drag enabled
 39   IF (ICEFRICTION.LE.0) THEN 
        CLOSE(IUT196)
      ELSE ! Read space-variable ice friction drag
        WRITE(IUPRT,*)
     + 'Space-variable ice-water friction enabled'
        READ(IUT196,11)  (COM(I),I=1,80)
        WRITE(IUPRT,11)  (COM(I),I=1,80)
        READ(IUT196,'(I5,E10.4)')NVARIF,WIFRIC
        WRITE(IUPRT,'(I5,E10.4)')NVARIF,WIFRIC
        IF(NVARIF.GT.MAXWET)THEN
          WRITE(IUPRT,56)NVARIF,MAXWET
          CALL SYSTEM ('rm gcm_temp*')
          STOP
        ENDIF
        READ(IUT196,11)  (COM(I),I=1,80)
        WRITE(IUPRT,11)  (COM(I),I=1,80)
        DO N=1,NVARIF
         READ(IUT196,55)II,JJ,VARWIF(II,JJ)
         VARWIF(II,JJ)=VARWIF(II,JJ)*WIFRIC
         AREAICE(II,JJ)=100. ! AND ALL ICE
         WRITE(IUPRT,55)II,JJ,VARWIF(II,JJ)
        END DO
        CLOSE(IUT196)
      ENDIF
 55   FORMAT(2I5,E10.3)
 56   FORMAT(/' NUMBER OF CELLS OF VARIABLE IFRIC ',I5,' IS GREATER THAN
     * MAXIMUM ALLOWED.',I5// ' PLEASE FIX THIS AND RESUBMIT'//)
!NG02182011   SPACE-AND-TIME VARIABLE ICE FRICTION
! A.Jordi: Set space-variable ice friction flag
      INQUIRE(FILE='ifric2dt.nc',EXIST=SET_ICE)
      IF (SET_ICE) ICEFRICTION=3 ! Space-and-time-Variable (3D) ice friction drag enabled
! A.Jordi: end
 40   IF (ICEFRICTION.LE.0) THEN 
        WRITE(IUPRT,*)
     + 'No Ice-Water Friction'
        CLOSE(IUT196)
      ELSEIF (ICEFRICTION.EQ.3) THEN 
        WRITE(IUPRT,*)
     + 'Space-and-time-variable ice-water friction enabled'
      ENDIF
!NG2011 SPACE-VARIABLE ICE DRAG OR SPACE-TIME VARIABLE ICE THICKNESS AND ICE CONCENTRATION. ...END
!
!  FOR INPUTTING WAM MODEL OUTPUT FOR WIND WAVE CALCULATIONS
!
1341  IF (WAVEDYN.EQ.'EXTERNAL') THEN
!
        OPEN (UNIT=110,FILE='wave_input',FORM='UNFORMATTED')
        OPEN (UNIT=111,FILE='gcm_temp111',FORM='UNFORMATTED')
!
        DO 8010 M=1,1000000
          READ (110,END=8100)TIME
          READ (110) WTEMP1   ! WAVE HEIGHT
          READ (110) WTEMP2   ! WAVE PERIOD
          READ (110) WTEMP3   ! WAVE DIRECTION
!
          WRITE (111)TIME
          WRITE (111) ((WTEMP1(I,J),I=1,IM),J=1,JM)
          WRITE (111) ((WTEMP2(I,J),I=1,IM),J=1,JM)
          WRITE (111) ((WTEMP3(I,J),I=1,IM),J=1,JM)
 8010   CONTINUE 
 8100   CLOSE (110)
      ENDIF
!
 5401 Format (//'OPTMBC HAS BEEN SELECTED AS ',A20//'THIS OPTION IS
     * NOT CORRECT, PLEASE FIX AND RESUBMIT'//)
 5402 Format (//'NUMDBCSE IS EQUAL TO ',I20//'THIS VALUE IS
     * NOT LESS OR EQUAL TO DBCM in comdeck',I20,//,
     * 'PLEASE FIX AND RESUBMIT'//)
 7402 Format (//'NUMDBCCH IS GIVEN EQUAL TO ',I20//'THIS VALUE IS
     * NOT LESS OR EQUAL TO DBCM in comdeck',I20,//,
     * 'PLEASE FIX AND RESUBMIT'//)
      RETURN
      END
