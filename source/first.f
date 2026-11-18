      SUBROUTINE FIRST
!
!     VERSION(10/02/90)
! 
      INCLUDE 'comdeck'
      
      CHARACTER STR*80 ! A.Jordi
      
      SAVE

! A.Jordi: Read elevation boundary data from netcdf file
      IF(NUMEBC.NE.0) THEN
!     IF(OPTEBC(1:5).EQ.'DATA ') THEN !NO NEED TO HAVE THIS ANYMORE(NKIM 121798)
        ! Read initial data
        NT_T2E=1
        CALL NETCDF_READ_BC_ELEV(NT_T2E,M_MAX,STR)
        DO WHILE (THOUR.GE.T2E)
          T1E=T2E 
          DO N=1,NUMEBC
            DEBDRY(N,1)=DEBDRY(N,2)
          END DO
          ! Read data at next time
          NT_T2E=NT_T2E+1
          CALL NETCDF_READ_BC_ELEV(NT_T2E,M_MAX,STR)
        END DO
      END IF
! A.Jordi: end
! 
! A.Jordi: change if order (originally NUMEBC was 2nd and HYDTYPE 1st)
      IF (NUMEBC.EQ.0) GOTO 135

! A.Jordi: Read elevation boundary data from netcdf file
      IF (HYDTYPE.NE.'EXTERNAL') THEN
        ! Read initial data
        NT_T2TS=1
        CALL NETCDF_READ_BC_TS(NT_T2TS,M_MAX,STR)
        DO WHILE (THOUR.GE.T2TS)
          T1TS=T2TS 
          DO K=1,KBM1
            DO N=1,NUMEBC
              DTBDRY(N,K,1)=DTBDRY(N,K,2)
              DSBDRY(N,K,1)=DSBDRY(N,K,2)
            END DO
          END DO
          ! Read data at next time
          NT_T2TS=NT_T2TS+1
          CALL NETCDF_READ_BC_TS(NT_T2TS,M_MAX,STR)
        END DO
      END IF
! A.Jordi: end
!
!  DISSOLVED TRACER AT OPEN BOUNDARY
!
      IF (TRACER.EQ.'INCLUDE'.AND.NUMEBCTR.GT.0) THEN
        REWIND (IUT501)
        DO 1010 M=1,100000
          IFLAG=3
          READ (IUT501,10,ERR=130) T2CONOB
!
          DO 1020 N=1,NUMEBCTR
            READ (IUT501,10) (DCBDRY1(N,K,2),K=1,KBM1)
 1020     CONTINUE
!
          IF(THOUR.LT.T2CONOB) GOTO 1030
!
          T1CONOB=T2CONOB
          DO 1040 K=1,KBM1
            DO 1040 N=1,NUMEBCTR
              DCBDRY1(N,K,1)=DCBDRY1(N,K,2)
 1040     CONTINUE
 1010     CONTINUE
 1030   CONTINUE
      ENDIF
!
!  SEDIMENT TRANSPORT AT OPEN BOUNDARY
!
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMEBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT502)
          DO 1110 M=1,100000
            IFLAG=4
            READ (IUT502,10,ERR=130) T2SEDOBM
!
            DO 1120 N=1,NUMEBCSE
              READ (IUT502,10) (DCBDRY(1,N,K,2),K=1,KBM1)
 1120       CONTINUE
!
            IF(THOUR.LT.T2SEDOBM) GOTO 1130
!
            T1SEDOBM=T2SEDOBM
            DO 1140 K=1,KBM1
              DO 1140 N=1,NUMEBCSE
                DCBDRY(1,N,K,1)=DCBDRY(1,N,K,2)
 1140       CONTINUE
 1110       CONTINUE
 1130     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT503)
          DO 1150 M=1,100000
            IFLAG=5
            READ (IUT503,10,ERR=130) T2SEDOBS
!
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
!
            DO 1160 N=1,NUMEBCSE
              READ (IUT503,10) (DCBDRY(KK,N,K,2),K=1,KBM1)
 1160       CONTINUE
!
            IF(THOUR.LT.T2SEDOBS) GOTO 1180
!
            T1SEDOBS=T2SEDOBS
            DO 1170 K=1,KBM1
              DO 1170 N=1,NUMEBCSE
                DCBDRY(KK,N,K,1)=DCBDRY(KK,N,K,2)
 1170       CONTINUE
 1150       CONTINUE
 1180     CONTINUE
        ENDIF
      ENDIF
!
!  PARTICLE-BOUND TRACER AT OPEN BOUNDARY
!
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMEBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT504)
          DO 1210 M=1,100000
            IFLAG=6
            READ (IUT504,10,ERR=130) T2CHMOBM
!
            KK=1
!
            DO 1220 N=1,NUMEBCCH
              READ (IUT504,10) (DPBDRY(KK,N,K,2),K=1,KBM1)
 1220       CONTINUE
!
            IF(THOUR.LT.T2CHMOBM) GOTO 1230
!
            T1CHMOBM=T2CHMOBM
            DO 1240 K=1,KBM1
              DO 1240 N=1,NUMEBCCH
                DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 1240       CONTINUE
 1210     CONTINUE
 1230     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT505)
          DO 1250 M=1,100000
            IFLAG=7
            READ (IUT505,10,ERR=130) T2CHMOBS
!
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
!
            DO 1260 N=1,NUMEBCCH
              READ (IUT505,10) (DPBDRY(KK,N,K,2),K=1,KBM1)
 1260       CONTINUE
!
            IF(THOUR.LT.T2CHMOBS) GOTO 1270
!
            T1CHMOBS=T2CHMOBS
            DO 1280 K=1,KBM1
              DO 1280 N=1,NUMEBCCH
                DPBDRY(KK,N,K,1)=DPBDRY(KK,N,K,2)
 1280       CONTINUE
 1250     CONTINUE
 1270     CONTINUE
        ENDIF
      ENDIF
!
 135  CONTINUE
!NG04232008 ELEVATION POTENTIAL CELLS
      IF(NUMEPC.EQ.0) GO TO 6540
       REWIND IUT65
       DO 6500 M=1,100000
       IFLAG=65
       READ (IUT65,10,ERR=130) T2EPC
       READ (IUT65,10) (DEPCOR(N,2),N=1,NUMEPC)
       IF(THOUR.LT.T2EPC) GO TO 6510
       T1EPC=T2EPC
       DO N=1,NUMEPC
        QEPCNUDF(N)=0.0 ! INITIALIZE
        IF (DEPCOR(N,2).NE.-99.) THEN
         DEPCOR(N,1)=DEPCOR(N,2)
        ELSE
         DEPCOR(N,1)=-99. ! FLAG FOR FREE SOLUTION
        ENDIF
       ENDDO
6500   CONTINUE
6510   CONTINUE
6540  CONTINUE
!
!-------- WAVE PARAMETERS BOUNDARY (hli, 03/19/04) -------------
!
!NG02222007      IF (WAVEDYN.EQ.'DONMODEL'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
      IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
!
        IF (NUMEBC.NE.0.AND.CREADWBC) THEN       
          REWIND IUT195
!
          DO M=1,100000
            IFLAG=16
            READ (IUT195,ERR=130) T2WAVE
            READ (IUT195)(DWHBFLX(N,2),N=1,NUMEBC)
            READ (IUT195)(DWPBFLX(N,2),N=1,NUMEBC)
            READ (IUT195)(DWDBFLX(N,2),N=1,NUMEBC)
            IF(THOUR.LT.T2WAVE) GOTO 1235
            T1WAVE=T2WAVE
            DO N=1,NUMEBC
              DWHBFLX(N,1)=DWHBFLX(N,2)
              DWPBFLX(N,1)=DWPBFLX(N,2)
              DWDBFLX(N,1)=DWDBFLX(N,2)
            ENDDO     
          ENDDO      
        ENDIF   
      ENDIF   
1235  CONTINUE
!
!  RIVER DISCHARGE BOUNDARY
!
! A.Jordi: Read discharge boundary condition from netcdf file 
      IF (NUMQBC.NE.0) THEN
        IF (HYDTYPE.NE.'EXTERNAL') THEN   
          NT_T2D=1
          CALL NETCDF_READ_BC_DIS(NT_T2Q,NT_T2Q_MAX,STR)
          DO WHILE (THOUR.GE.T2Q)
            T1Q=T2Q
            DO N=1,NUMQBC
              DQDIS(N,1)=DQDIS(N,2)
              DTDIS(N,1)=DTDIS(N,2)
              DSDIS(N,1)=DSDIS(N,2)
            END DO 
            ! Read data at next time
            NT_T2D=NT_T2D+1
            CALL NETCDF_READ_BC_DIS(NT_T2Q,NT_T2Q_MAX,STR)
          END DO
        END IF
      END IF
! A.Jordi: end
!
!  DISSOLVED TRACER AT RIVER DISCHARGE BOUNDARY
!
      IF (TRACER.EQ.'INCLUDE'.AND.NUMQBCTR.GT.0) THEN
        REWIND (IUT601)
        DO 1275 M=1,100000
          IFLAG=9
          READ (IUT601,10,ERR=130) T2CON
          READ (IUT601,10) (DCDIS1(N,2),N=1,NUMQBCTR)
!
          IF(THOUR.LT.T2CON) GO TO 1285
          T1CON=T2CON 
!
          DO 1295 N=1,NUMQBCTR
            DCDIS1(N,1)=DCDIS1(N,2) 
 1295     CONTINUE
 1275   CONTINUE
 1285   CONTINUE
      ENDIF
!
!  SEDIMENT TRANSPORT AT RIVER DISCHARGE BOUNDARY
!
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMQBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT602)
          DO 1300 M=1,100000
            IFLAG=10
            READ (IUT602,10,ERR=130) T2SEDM
            KK=1
            READ (IUT602,10) (DCDIS(KK,N,2),N=1,NUMQBCSE)
            IF(THOUR.LT.T2SEDM) GO TO 1310
            T1SEDM=T2SEDM 
!
            DO 1320 N=1,NUMQBCSE
              DCDIS(KK,N,1)=DCDIS(KK,N,2) 
 1320       CONTINUE
 1300     CONTINUE
 1310     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT603)
          DO 1330 M=1,100000
            IFLAG=11
            READ (IUT603,10,ERR=130) T2SEDS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT603,10) (DCDIS(KK,N,2),N=1,NUMQBCSE)
            IF(THOUR.LT.T2SEDS) GO TO 1340
            T1SEDS=T2SEDS 
!
            DO 1350 N=1,NUMQBCSE
              DCDIS(KK,N,1)=DCDIS(KK,N,2) 
 1350       CONTINUE
 1330     CONTINUE
 1340     CONTINUE
        ENDIF
      ENDIF
!
!  PARTICLE-BOUND TRACER AT RIVER DISCHARGE BOUNDARY
!
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMQBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT604)
          DO 1400 M=1,100000
            IFLAG=12
            READ (IUT604,10,ERR=130) T2CHMM
            KK=1
            READ (IUT604,10) (DPDIS(KK,N,2),N=1,NUMQBCCH)
            IF(THOUR.LT.T2CHMM) GO TO 1410
            T1CHMM=T2CHMM 
!
            DO 1420 N=1,NUMQBCCH
              DPDIS(KK,N,1)=DPDIS(KK,N,2) 
 1420       CONTINUE
 1400     CONTINUE
 1410     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT605)
          DO 1430 M=1,100000
            IFLAG=13
            READ (IUT605,10,ERR=130) T2CHMS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT605,10) (DPDIS(KK,N,2),N=1,NUMQBCCH)
            IF(THOUR.LT.T2CHMS) GO TO 1440
            T1CHMS=T2CHMS 
!
            DO 1450 N=1,NUMQBCCH
              DPDIS(KK,N,1)=DPDIS(KK,N,2) 
 1450       CONTINUE
 1430     CONTINUE
 1440     CONTINUE
        ENDIF
      ENDIF
! 
! add 'DIFFUSER IN LOOP' to the codes, referring to TAMPA BAY codes
! in /power6/wrai0010/MODEL/CODES/TRACER (hli, 09/28/00)
!
!  ************************************************************
!  DIFFUSER BOUNDARY
!
! A.Jordi: change if order (originally NUMEBC was 2nd and HYDTYPE 1st)
      IF(NUMDBC1.EQ.0) GO TO 1510
      
! A.Jordi: Read diffuser 1 boundary condition from netcdf file 
        IF (HYDTYPE.NE.'EXTERNAL') THEN     
          NT_T2D=1
          CALL NETCDF_READ_BC_DIFF1(NT_T2D,NT_T2D_MAX,STR)
          DO WHILE (THOUR.GE.T2D)
            T1D=T2D
            DO N=1,NUMDBC1
              DQDIFF(N,1)=DQDIFF(N,2)
              DTDIFF(N,1)=DTDIFF(N,2)
              DSDIFF(N,1)=DSDIFF(N,2)
            END DO 
            ! Read data at next time
            NT_T2D=NT_T2D+1
            CALL NETCDF_READ_BC_DIFF1(NT_T2D,NT_T2D_MAX,STR)
          END DO
        END IF
! A.Jordi: end
!
!  DISSOLVED TRACER AT DIFFUSER
!
      IF (TRACER.NE.'INCLUDE'.OR.NUMDBCTR1.EQ.0) GO TO 1510
      REWIND IUT98
!
      DO M=1,100000
        IFLAG=15
        READ (IUT98,10,ERR=130) T2DCON
        READ (IUT98,10) (DCDIFF1(N,2),N=1,NUMDBCTR1)
!
        IF(THOUR.LT.T2DCON) GO TO 1510
!
        T1DCON=T2DCON
!
        DO N=1,NUMDBCTR1
          DCDIFF1(N,1)=DCDIFF1(N,2)
        ENDDO     
      ENDDO     
 1510 CONTINUE
!
!  DIFFUSER IN LOOP BOUNDARY
!
      IF(NUMDBC2.EQ.0) GO TO 1511
      
! A.Jordi: Read diffuser 2 boundary condition from netcdf file 
        NT_T2D2=1
        CALL NETCDF_READ_BC_DIFF2(NT_T2D2,NT_T2D2_MAX,STR)
        DO WHILE (THOUR.GE.T2D2)
          T1D2=T2D2
          DO N=NUMDBC1+1,NUMDBC1+NUMDBC2
            DQDIFF(N,1)=DQDIFF(N,2)
            DTDIFF(N,1)=DTDIFF(N,2)
            DSDIFF(N,1)=DSDIFF(N,2)
          END DO
          ! Read data at next time
          NT_T2D2=NT_T2D2+1
          CALL NETCDF_READ_BC_DIFF2(NT_T2D2,NT_T2D2_MAX,STR)
        END DO
! A.Jordi: end
!
!  DISSOLVED TRACER AT DIFFUSER IN LOOPS
!        GO TO 1511   ! DO NOT USE THIS FOR NOW
      IF (TRACER.NE.'INCLUDE'.OR.NUMDBCTR2.EQ.0) GO TO 1511
      REWIND IUT99
      DO M=1,100000
        IFLAG=150
        READ (IUT99,10,ERR=130) T2D2CON
        READ (IUT99,10) (DCDIFF1(N,2),N=NUMDBCTR1+1,NUMDBCTR)
!
        IF(THOUR.LT.T2D2CON) GO TO 1511
!
        T1D2CON=T2D2CON
!
        DO N=NUMDBCTR1+1,NUMDBCTR
          DCDIFF1(N,1)=DCDIFF1(N,2)
        ENDDO   
      ENDDO   
 1511 CONTINUE
!
!  end of adding DIFFUSER IN LOOP (hli, 09/28/00)
!
!  POINT SOURCE LOADS:  DISSOLVED TRACER
!
      IF (TRACER.EQ.'INCLUDE'.AND.NUMPSTR.GT.0) THEN
        REWIND (IUT701)
        DO 4275 M=1,100000
          IFLAG=91
          READ (IUT701,10,ERR=130) T2PSTR
          READ (IUT701,10) (DCPSTR(N,2),N=1,NUMPSTR)
!
          IF(THOUR.LT.T2PSTR) GO TO 4285
!
          T1PSTR=T2PSTR 
!
          DO 4295 N=1,NUMPSTR
            DCPSTR(N,1)=DCPSTR(N,2) 
 4295     CONTINUE
 4275   CONTINUE
 4285   CONTINUE
      ENDIF
!
!  SEDIMENT LOADS AT DIFFUSER ***** Quamrul QA 5/12/99
!
      IF (SEDTRAN.EQ.'INCLUDE'.AND.NUMDBCSE.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT702)
          DO 7300 M=1,100000
            IFLAG=702
            READ (IUT702,10,ERR=730) T2SEDDM
            KK=1
            READ (IUT702,10) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
            IF(THOUR.LT.T2SEDDM) GO TO 7310
            T1SEDDM=T2SEDDM
!
            DO 7320 N=1,NUMDBCSE
              DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7320       CONTINUE
 7300     CONTINUE
 7310     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT703)
          DO 7330 M=1,100000
            IFLAG=703
            READ (IUT703,10,ERR=730) T2SEDDS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT703,10) (DCSDIFF(KK,N,2),N=1,NUMDBCSE)
            IF(THOUR.LT.T2SEDDS) GO TO 7340
            T1SEDDS=T2SEDDS
!
            DO 7350 N=1,NUMDBCSE
              DCSDIFF(KK,N,1)=DCSDIFF(KK,N,2)
 7350       CONTINUE
 7330     CONTINUE
 7340     CONTINUE
        ENDIF
      ENDIF
!
!  PARTICLE-BOUND TRACER AT DIFFUSER BOUNDARY By Quamrul QA 5/13/99
!
      IF (CHEMTRAN.EQ.'INCLUDE'.AND.NUMDBCCH.GT.0) THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT704)
          DO 7400 M=1,100000
            IFLAG=704
            READ (IUT704,10,ERR=130) T2CHMDM
            KK=1
            READ (IUT704,10) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
            IF(THOUR.LT.T2CHMDM) GO TO 7410
            T1CHMDM=T2CHMDM
!
            DO 7420 N=1,NUMDBCCH
              DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7420       CONTINUE
 7400     CONTINUE
 7410     CONTINUE
        ENDIF
!
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          REWIND (IUT705)
          DO 7430 M=1,100000
            IFLAG=705
            READ (IUT705,10,ERR=130) T2CHMDS
            IF (SEDTYPE.EQ.'SAND') THEN
              KK=1
            ELSE
              KK=2
            ENDIF
            READ (IUT705,10) (DPDIFF(KK,N,2),N=1,NUMDBCCH)
            IF(THOUR.LT.T2CHMDS) GO TO 7440
            T1CHMDS=T2CHMDS
!
            DO 7450 N=1,NUMDBCCH
              DPDIFF(KK,N,1)=DPDIFF(KK,N,2)
 7450       CONTINUE
 7430     CONTINUE
 7440     CONTINUE
        ENDIF
      ENDIF
!
! 
      IF(OPTMBC(1:8).EQ.'AVERAGED') THEN
       REWIND IUT93
       DO 400 M=1,100000
            IFLAG=93
       READ (IUT93,10,ERR=130) T2M
       READ (IUT93,10) DHFLUX(2),DTX(2),DTY(2),
     +   DWSX(2),DWSY(2),DWDS(2),DWDD(2),DQPREC(2),DQEVAP(2)
       IF(THOUR.LT.T2M) GO TO 410
       T1M=T2M
       DQPREC(1)=DQPREC(2) 
       DQEVAP(1)=DQEVAP(2) 
       DHFLUX(1)=DHFLUX(2)
       DTX   (1)=DTX   (2) 
       DTY   (1)=DTY   (2) 
       DWSX  (1)=DWSX  (2) 
       DWSY  (1)=DWSY  (2) 
       DWDS  (1)=DWDS  (2) 
       DWDD  (1)=DWDD  (2) 
 400   CONTINUE
 410   CONTINUE
! MOVE SYNOPTIC DOWN
!hli      ELSE IF(OPTMBC(1:8).EQ.'SYNOPTIC') THEN    ! JI
!C       REWIND IUT93
!hli       REWIND IUT192
!hli       DO 500 M=1,100000
!hli            IFLAG=92
*Khan 072898
!C       READ (IUT93,10,ERR=130) T2M
!C       READ (IUT93,10) DQPREC(2),DQEVAP(2),DHFLUX(2)
!hli       READ (IUT192,ERR=130) T2M
!hli       READ (IUT192) ((DHFLX2D(I,J,2),I=1,IM),J=1,JM)
*Khan 072898
!hli       IF(THOUR.LT.T2M) GO TO 510
!hli       T1M=T2M
*Khan 072898
!C       DQPREC(1)=DQPREC(2) 
!C       DQEVAP(1)=DQEVAP(2) 
!C       DHFLUX(1)=DHFLUX(2)
!hli       DO 511 J=1,JM       ! n.kim 7.7.98
!hli       DO 511 I=1,IM
!hli        DHFLX2D(I,J,1)=DHFLX2D(I,J,2)
!hli 511   CONTINUE
*Khan 072898
!hli 500   CONTINUE
!hli 510   CONTINUE
!hli
!------  WIND STRESS  ----------------------------------------
!hli       REWIND IUT95
!hli       DO 520 M=1,100000
!hli            IFLAG=95
!hli       READ (IUT95,ERR=130) T2W
!hli       READ (IUT95) ((DTX2D(I,J,2),DTY2D(I,J,2),DPATM(I,J,2),
!hli     +                I=1,IM),J=1,JM)
!hli       IF(THOUR.LT.T2W) GO TO 530
!hli       T1W=T2W
!hli       DO 540 J=1,JM
!hli       DO 540 I=1,IM
!hli       DPATM(I,J,1)=DPATM(I,J,2) 
!hli       DTX2D(I,J,1)=DTX2D(I,J,2) 
!hli 540   DTY2D(I,J,1)=DTY2D(I,J,2) 
!hli 520   CONTINUE
!hli 530   CONTINUE
!hli

!qa        ELSE IF (OPTMBC(1:8).EQ.'HEATFLUX') THEN
!
!       MODIFIED BY QUAMRUL QA ON 3/10/99 ***************************************
!  Add three new options (SYNOPANB, SYNOPLNP and SYNOPRNM) 
!  (hli,04/24/02)
!
       else If (OPTMBC(1:8).EQ.'LANDPFLX'.OR.
     *          OPTMBC(1:8).EQ.'RANDMFLX'.OR.
     *          OPTMBC(1:8).EQ.'AANDBFLX'.OR.
     *          OPTMBC(1:8).EQ.'SYNOPANB'.OR.
     *          OPTMBC(1:8).EQ.'SYNOPLNP'.OR.
     *          OPTMBC(1:8).EQ.'SYNOPRNM') Then
!qa
        REWIND IUT93
        DO 150 M = 1, 100000
          READ (IUT93,810,ERR=130) T2M
          READ (IUT93,810) DAIRTM(2), DRELHU(2),
     *        DBAROP(2), DSWOBS(2), DTX(2), DTY(2), DWSX(2), DWSY(2),
     *     CLOUD(2), EXTCOEF(2), DQPREC(2), DQEVAP(2)
810     format(12e14.7)
!
          IF (THOUR.LT.T2M) GO TO 160
          T1M = T2M
          DQPREC(1) = DQPREC(2)
          DQEVAP(1) = DQEVAP(2)
          DAIRTM(1) = DAIRTM(2)
          DRELHU(1) = DRELHU(2)
          DBAROP(1) = DBAROP(2)
          DSWOBS(1) = DSWOBS(2)
!  MODIFIED BY QUAMRUL QA ON 5/3/96 *******************************************
          CLOUD(1)  = CLOUD(2)
          EXTCOEF(1)  = EXTCOEF(2)
          DTX(1) = DTX(2)
          DTY(1) = DTY(2)
          DWSX(1) = DWSX(2)
          DWSY(1) = DWSY(2)
          DWDS  (1)=DWDS  (2) 
          DWDD  (1)=DWDD  (2) 
  150   CONTINUE
  160   CONTINUE

!------  VARYING EXTINCTION COEFFICIENT ----------------------
!
       IF(OPTEXTC.EQ.'VARI') THEN
         REWIND IUT194
         DO M=1,100000
           IFLAG=194
           READ (IUT194,ERR=130) T2EX
           READ (IUT194)((DEXTC(I,J,2),I=1,IM),J=1,JM)
           IF(THOUR.LT.T2EX) GO TO 1530
           T1EX=T2EX
           DO J=1,JM
             DO I=1,IM
               DEXTC(I,J,1)=DEXTC(I,J,2) 
             ENDDO
           ENDDO                  
         ENDDO
1530     CONTINUE
       ENDIF
!
      ENDIF
!
       IF ( OPTMBC(1:8).EQ.'SYNOPANB'.OR.
     &      OPTMBC(1:8).EQ.'SYNOPLNP'.OR.
     &      OPTMBC(1:8).EQ.'SYNOPRNM'.OR.
     &      OPTMBC(1:8).EQ.'SYNOPTIC') Then
        IF(OPTMBC(1:8).EQ.'SYNOPTIC') Then
          REWIND IUT192
          DO M=1,100000
            IFLAG=92
            READ (IUT192,ERR=130) T2M
            READ (IUT192) ((DHFLX2D(I,J,2),I=1,IM),J=1,JM)
            IF(THOUR.LT.T2M) GOTO 510   
            T1M=T2M
            DO J=1,JM       ! n.kim 7.7.98
              DO I=1,IM
                DHFLX2D(I,J,1)=DHFLX2D(I,J,2)
              ENDDO
            ENDDO
          ENDDO
510       CONTINUE
        ENDIF

        REWIND IUT95
        DO M=1,100000
          IFLAG=18
          READ (IUT95,ERR=130) T2W
          READ (IUT95) ((DTX2D(I,J,2),DTY2D(I,J,2),DPATM(I,J,2),
     .                   DWSX2D(I,J,2),DWSY2D(I,J,2), ! NG04032008
     .                   I=1,IM),J=1,JM)
          IF(THOUR.LT.T2W) GO TO 630
          T1W=T2W
          DO J=1,JM ! NG04102008
            DO I=1,IM
              DPATM(I,J,1)=DPATM(I,J,2)
              DTX2D(I,J,1)=DTX2D(I,J,2)
              DTY2D(I,J,1)=DTY2D(I,J,2)
              DWSX2D(I,J,1)=DWSX2D(I,J,2) ! NG04032008
              DWSY2D(I,J,1)=DWSY2D(I,J,2) ! NG04032008
            ENDDO
          ENDDO
        ENDDO
 630    CONTINUE
!
      ENDIF  
!                
!NG04092008 BELOW NOT NEEDED ANY MORE. USE WSX2D, WSY2D FOR WINDS FROM NOW ON.
!NG04092008      IF (WAVEDYN(1:3).EQ.'DON'.OR.WAVEDYN.EQ.'SMBMODEL') THEN
!NG04092008        REWIND IUT193
!NG04092008        DO M=1,100000
!NG04092008          IFLAG=193
!NG04092008          READ (IUT193,ERR=130) T2WV
!NG04092008          READ (IUT193) ((WU2(I,J),WV2(I,J),I=1,IM),J=1,JM)
!NG04092008          IF(THOUR.LT.T2WV) GO TO 5311
!NG04092008          T1WV=T2WV
!NG04092008          DO J=1,JM
!NG04092008            DO I=1,IM
!NG04092008              WU1(I,J)=WU2(I,J)
!NG04092008              WV1(I,J)=WV2(I,J)
!NG04092008            ENDDO
!NG04092008           ENDDO
!NG04092008        ENDDO
!NG04092008 5311   CONTINUE
!NG04092008      ENDIF

!svv080326 inserted
!  READ 2D WIND AND HEAT FLUX PARAMETERS
      IF (OPTMBC(1:8).EQ.'2DMETANB'.OR.
     .    OPTMBC(1:8).EQ.'2DMETLNP'.OR.
     .    OPTMBC(1:8).EQ.'2DMETRNM') THEN
! A.Jordi: Read wind and heat flux parameters from netcdf file
        ! Read initial data
        NT_T2M=1
        CALL NETCDF_READ_SYNOP_MET(NT_T2M)
        WRITE(IUPRT,*) T2M
        DO WHILE (THOUR.GE.T2M)
          T1M=T2M
          DO I=1,IM
            DO J=1,JM
              DAIRTM2D(I,J,1)=DAIRTM2D(I,J,2)
              DRELHU2D(I,J,1)=DRELHU2D(I,J,2)
              DPATM(I,J,1)=DPATM(I,J,2)
              DSWOBS2D(I,J,1)=DSWOBS2D(I,J,2)
              DWSX2D(I,J,1)=DWSX2D(I,J,2)
              DWSY2D(I,J,1)=DWSY2D(I,J,2)  
              DCLD2D(I,J,1)=DCLD2D(I,J,2)
              DEXTC2D(I,J,1)=DEXTC2D(I,J,2)
            END DO
          END DO
          ! Read data at next time
          NT_T2M=NT_T2M+1
          CALL NETCDF_READ_SYNOP_MET(NT_T2M)
          WRITE(IUPRT,*) T2M
        END DO
      END IF
! A.Jordi: end

!NG 2011 FOR VARIABLE ICE THICKNESS AND AREA. START...
      IF (ICEFRICTION.EQ.3) THEN
! A.Jordi: Read space-variable ice friction from netcdf file
        NT_T2I=1
        CALL NETCDF_READ_IFRIC2DT(NT_T2I)
        WRITE(IUPRT,*) T2I
        DO WHILE (THOUR.GE.T2I)
          T1I=T2I
          DO I=1,IM
            DO J=1,JM
              DAREAICE(I,J,1) = DAREAICE(I,J,2)
              DTHICKICE(I,J,1) = DTHICKICE(I,J,2)
            END DO
          END DO
          ! Read data at next time
          NT_T2I=NT_T2I+1
          CALL NETCDF_READ_IFRIC2DT(NT_T2I)
          WRITE(IUPRT,*) T2I
        END DO
! A.Jordi: end
      END IF
!NG 2011 FOR VARIABLE ICE THICKNESS AND AREA. ...END

!NG09022014 For synop_pet, if it exists...
      IF(CREADPET) THEN
       REWIND IUT104
       DO 721 M=1,100000
        IFLAG=104
        READ (IUT104,ERR=130) T2P
        write(IUPRT,*)T2P
        READ(IUT104) ((DQPREC2D(I,J,2),DQEVAP2D(I,J,2),I=1,IM),J=1,JM)
        IF(THOUR.LT.T2P) GO TO 731
        T1P=T2P
        DO 741 J=1,JM
         DO 741 I=1,IM
          DQPREC2D(I,J,1) = DQPREC2D(I,J,2) ! NG2014 Use data
 741      DQEVAP2D(I,J,1) = DQEVAP2D(I,J,2) ! NG2014 Use data
 721   CONTINUE
 731   CONTINUE
      END IF
!NG09022014 ...For synop_pet, if it exists.
!
!hli
!
!  FOR WIND WAVE INPUT
!
      IF (WAVEDYN.EQ.'EXTERNAL') THEN
        REWIND (111)
        DO 8010 M=1,100000
          READ (111,ERR=8020) T2WAVE
          READ (111) ((DWVHT(I,J,2),I=1,IM),J=1,JM)
          READ (111) ((DWVPD(I,J,2),I=1,IM),J=1,JM)
          READ (111) ((DWVDR(I,J,2),I=1,IM),J=1,JM)
!
          IF (THOUR.LT.T2WAVE) GO TO 8020
!
          T1WAVE=T2WAVE
!
          DO 8030 J=1,JM
            DO 8030 I=1,IM
              DWVHT(I,J,1)=DWVHT(I,J,2) 
              DWVPD(I,J,1)=DWVPD(I,J,2) 
              DWVDR(I,J,1)=DWVDR(I,J,2) 
 8030     CONTINUE
 8010   CONTINUE
 8020   CONTINUE
      ENDIF
!
      RETURN
!
 130  WRITE(6,20) IFLAG
 20   FORMAT(//' THERE IS INSUFFICIENT TEMPORAL DATA FOR THIS RUN'/,
     .         '        REVISE INPUT DECK AND RESUBMIT '
     +     // 'IFLAG =',I4)
!
 730  WRITE(6,70) IFLAG
 70   FORMAT(//' THERE IS INSUFFICIENT TEMPORAL DATA FOR THIS RUN'/,
     .         '        REVISE INPUT DECK AND RESUBMIT '
     +     // 'IFLAG =',I4)
 10   FORMAT(8E14.7)
      STOP
      END
