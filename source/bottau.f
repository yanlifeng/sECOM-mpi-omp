      SUBROUTINE BOTTAU
!
!  CALCULATES BOTTOM SHEAR STRESSES FOR USE IN SEDIMENT TRANSPORT
!  (SUBROUTINES sedflx AND suslod)
!
!  REVISION DATE :  AUGUST 6, 1996
!
!**********************************************************************
!
      INCLUDE 'comdeck'
      INTEGER JOTLO,JOTHI
      SAVE
!
!  CALC. BOTTOM SHEAR STRESS
!
      JOTLO=2
      JOTHI=JMM1
      IF(SMPI_OUTER) THEN
        IF(TOR.NE.'BAROTROPIC'.OR.WAVEDYN.NE.'NEGLECT'.OR.
     &     NUMQBC.NE.0) STOP 97
        JOTLO=SMPI_FIRST(2)
        JOTHI=SMPI_LAST(JMM1)
! Step-entry VA is the owner from the last ROTATE, not an output
! snapshot. VBAR on the last owned row reads the next owner's VA.
        CALL SMPI_HALO_RAW(VA,IM,JM)
      END IF
      IF(TOR.EQ.'BAROTROPIC') THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
       DO 10 J=JOTLO,JOTHI
        DO 10 I=2,IMM1
!         UBAR(I,J,KBM1)=(DUMWAD(I,J)*UA(I,J)+DUMWAD(I+1,J)*UA(I+1,J))/
!     +                  (DUMWAD(I,J)+DUMWAD(I+1,J)+1.0E-30)
!         VBAR(I,J,KBM1)=(DVMWAD(I,J)*VA(I,J)+DVMWAD(I,J+1)*VA(I,J+1))/
!     +                  (DVMWAD(I,J)+DVMWAD(I,J+1)+1.0E-30)
             UBAR(I,J,KBM1)=0.5*(UA(I,J)+UA(I+1,J))
             VBAR(I,J,KBM1)=0.5*(VA(I,J)+VA(I,J+1))
             DT(I,J)=D(I,J) ! NG11012013, for stress that needs DT, use D.
 10    CONTINUE
      ELSE
       DO 11 K=1,KBM1
        DO 11 J=2,JMM1
         DO 11 I=2,IMM1
!          UBAR(I,J,K)=(DUMWAD(I,J)*U(I,J,K)+DUMWAD(I+1,J)*U(I+1,J,K))/
!     +                (DUMWAD(I,J)+DUMWAD(I+1,J)+1.0E-30)
!          VBAR(I,J,K)=(DVMWAD(I,J)*V(I,J,K)+DVMWAD(I,J+1)*V(I,J+1,K))/
!     +                (DVMWAD(I,J)+DVMWAD(I,J+1)+1.0E-30)
               UBAR(I,J,K)=0.5*(U(I,J,K)+U(I+1,J,K))
               VBAR(I,J,K)=0.5*(V(I,J,K)+V(I,J+1,K))
 11    CONTINUE
      ENDIF
!
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
      DO 20 J=JOTLO,JOTHI
        DO 20 I=2,IMM1
           QBAR(I,J)=SQRT(UBAR(I,J,KBM1)*UBAR(I,J,KBM1)+
     +                    VBAR(I,J,KBM1)*VBAR(I,J,KBM1))
 20   CONTINUE
!
!-------- VELOCITY BOUNDARY CONDITION ----------------------------------
!
!qa      Do 21 N = 1, NUMQBC
!qa        ID = IQD(N)
!qa        JD = JQD(N)
!qa        IC = IQC(N)
!qa        JC = JQC(N)
!
!qa        Do 22 K = 1, KBM1
!qa          FRESH = QDIS(N) * RAMP * VQDIST(N,K) / 100.
!qa          If (JD.EQ.JC) Then
!qa            If (ID.LT.IC) Then
!qa              UBAR(IC,JC,K) = U(ID,JD,K)
!qa            Else
!qa              UBAR(ID,JD,K) = U(ID+1,JD,K)
!qa            End If
!qa          Else
!qa            If (JD.LT.JC) Then
!qa              VBAR(IC,JC,K) = -FRESH/(DZ(K)*DT(I,J)*H1(ID,JD))
!qa            Else
!qa              VBAR(ID,JD,K) = FRESH/(DZ(K)*DT(I,J)*H1(ID,JD))
!qa            End If
!qa          End If
!qa   22   Continue
!qa   21 Continue
!
!qa      IF (HYDTYPE.EQ.'INTERNAL') THEN
!qa        DO 30 K=1,KB
!qa          DO 30 J=2,JMM1
!qa            DO 30 I=2,IMM1
!qa              TAU(I,J,K)=0.0
!qa 30     CONTINUE
!qa      ELSE
!qa        DO 35 J=2,JMM1
!qa          DO 35 I=2,IMM1
!qa            TAU(I,J,KB)=0.0
!qa 35     CONTINUE
!qa      ENDIF
!
!  NO WIND WAVES
!
!NG02222007
      IF (WAVEDYN.EQ.'DONONLY ') THEN
!NG02222007
!NG12222010 ! INCORPORATED NWAVE AND BCOND CALL HERE LIKE FOR DONMODEL
        NWAVECNT=NWAVECNT+1
        IF (NWAVECNT.EQ.1) THEN
          CALL BCOND(11,DTI2,0)
          CALL WAVEDON
        ENDIF
        IF (NWAVECNT.EQ.NWAVE) NWAVECNT=0
!NG12222010
      ELSEIF (WAVEDYN.EQ.'NEGLECT ') THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
        DO 40 J=JOTLO,JOTHI
          DO 40 I=2,IMM1
            IF (FSM(I,J).GT.0.0) THEN
!
!  BOTTOM STRESS (AT K= KB)
!
!  IN SED TRANSPORT, THE CBC DOES NOT NEED TO BE SCALED. OTHERWISE, IT WOULD
!  RESUSPEND EXCESSIVE SEDIMENT (PER PRAVI, NKIM 10.23.02)
!             TAU(I,J,KB)=10000.*CBC(I,J)*QBAR(I,J)*QBAR(I,J)
             TAU(I,J,KB)=10000.*CBC(I,J)/VARBF(I,J)*QBAR(I,J)*QBAR(I,J)
            ENDIF
 40     CONTINUE
!
!     FORCING TAU VALUE = 0 AT CELLS CONNECTED BY RIVERS (IC,JC)
!     By Quamrul QA 10/28/98
!
        DO 120 N=1,NUMQBC
          IC=IQC(N)
          JC=JQC(N)
          TAU(IC,JC,KB) = 0.0
120     CONTINUE

      ELSE
!
!  INCLUDE WIND WAVE EFFECTS
!
!NKIM11232010        IF (WAVEDYN.EQ.'SMBMODEL') CALL WAVESMB
!NKIM11232010        IF (WAVEDYN.EQ.'DONMODEL') CALL WAVEDON
!NKIM11232010 BRING NWAVECNT BELOW HERE
!NG12022010 NEED TO BE CALLED IN INITIAL POINT
!NG12022010        IF (WAVEDYN.EQ.'SMBMODEL'.OR.
!NG12022010     +      WAVEDYN.EQ.'DONMODEL'.OR.
!NG12022010     +      WAVEDYN.EQ.'EXTERNAL') THEN
!NG12022010            NWAVECNT=NWAVECNT+1
!NG12022010            IF (NWAVECNT.EQ.NWAVE) THEN
!NG12022010              NWAVECNT=0
!NG12022010              CALL BCOND(11,DTI2,0) ! UPDATE WAVE BC (NKIM 11.23.2010)
!NG12022010              IF (WAVEDYN.EQ.'SMBMODEL') CALL WAVESMB
!NG12022010              IF (WAVEDYN.EQ.'DONMODEL') CALL WAVEDON
!NG12022010            ENDIF
!NG12022010        ENDIF
        IF (WAVEDYN.EQ.'SMBMODEL'.OR.
     +      WAVEDYN.EQ.'DONMODEL'.OR.
     +      WAVEDYN.EQ.'EXTERNAL') THEN
            NWAVECNT=NWAVECNT+1
            IF (NWAVECNT.EQ.1) THEN
              CALL BCOND(11,DTI2,0) ! UPDATE WAVE BC (NKIM 11.23.2010)
              IF (WAVEDYN.EQ.'SMBMODEL') CALL WAVESMB
              IF (WAVEDYN.EQ.'DONMODEL') CALL WAVEDON
            ENDIF
            IF (NWAVECNT.EQ.NWAVE) NWAVECNT=0
        ENDIF
!NKIM11232010 BRING NWAVECNT ABOVE HERE
        CALL STRESS
      ENDIF
!
      RETURN
      END
