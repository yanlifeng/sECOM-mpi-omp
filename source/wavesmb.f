      SUBROUTINE WAVESMB
!
!  CALC. BOTTOM SHEAR STRESS DUE TO WIND WAVES AND CURRENTS
!  USES GRANT/MADSEN/GLENN WAVE MODEL
!  REVISION DATE:  SEPTEMBER 26, 1997
!**************************************************************
!
      INCLUDE 'comdeck'
      SAVE
!
!sfan moved to comdeck      dimension FSMW(IM,JM)

      PI=  3.141593
      TPI= 6.283185
      HPI= 1.570796
!
!    The EXTERNAL should be taken care of later
!qa      IF (WAVEDYN.EQ.'EXTERNAL') GOTO 15   
!
!  CALM WIND ( < 1 MPH = 0.45 m/s) ==> PURE CURRENTS
!
!PL      IF (WINDSP.LT.0.45) GOTO 99
!
!  CALC. SIGNIFICANT WAVE HEIGHT AND PERIOD
!  USE SHALLOW WATER SMB THEORY
!
!  WVHT = SIG. WAVE HEIGHT (m)
!  WVPD = SIG. WAVE PERIOD (s)
!
!  WINDSP = WIND SPEED (m/s)
!  WINDIR = WIND DIRECTION (degrees)
!  HMEAN(N) = MEAN DEPTH ALONG FETCH FOR DIRECTION n (m)
!  FETCH(N) = FETCH FOR DIRECTION n (m)
!
!  Previously WINDSP and AINDIR was passed into this subroutine calculated
!  in bcond.f was wrong. We introduce spatially variable wind field and
!  incorporate correct interpolation to compute WSX2D and WSY2D at time THOUR
!  by Quamrul QA 10/27/98
!
!hli
      DO J=1,JM
         DO I=1,IM
            FSMW(I,J)=FSM(I,J)
         ENDDO
      ENDDO
      DO N=1,NUMQBC
         IC=IQC(N)
         JC=JQC(N)
         FSMW(IC,JC)=0.0
      ENDDO
!hli
      DO 10 J=2,JM-1
        DO 10 I=2,IM-1
          IF (FSMW(I,J).GT.0.0) THEN
!
!PL 4/20/99 ADD SMALL VALUE
!
           UWIND=(WSX2D(I,J)*COS(ANG(I,J))+WSY2D(I,J)*SIN(ANG(I,J)))*
     &        FSMW(I,J)
           VWIND=(-WSX2D(I,J)*SIN(ANG(I,J))+WSY2D(I,J)*COS(ANG(I,J)))*
     &        FSMW(I,J)
           WSPDSQ=UWIND*UWIND+VWIND*VWIND
           WINDSP=SQRT(WSPDSQ)+1.E-10
!
!     WINDIR is w.r.t. zeta1-direction
!
           WINDIR=ATAN2(VWIND,UWIND)*180./PI
           IF(WINDIR .LT. 0.) WINDIR = WINDIR + 360
           IF(WINDIR .GT. 360.) WINDIR = WINDIR - 360.
!
!     Direction with respect to east (x-direction)
!
           WINDIR = WINDIR + ANG(I,J)*180./PI
           WINDIR = AMOD(WINDIR,360.)
!
!     IN stress.f and in SMB Model WAVE DIR IS IN w.r.t North Clock Wise
!
           WINDIR= 90. - WINDIR
           IF(WINDIR.LE.0.) WINDIR = WINDIR + 360.
           NDIR=NINT(WINDIR/10.)
           IF (NDIR.EQ.0) NDIR=36
           WRAT0=GRAV/(WINDSP*WINDSP)
           HRAT0=0.283/WRAT0
           TRAT0=2.4*PI*WINDSP/GRAV
           WRAT1=WRAT0*HMEAN(NDIR,I,J)
           WRAT2=WRAT0*FETCH(NDIR,I,J)
           HA1=TANH(0.53*(WRAT1**0.75))
           TA1=TANH(0.833*(WRAT1**0.375))
           WVHT(I,J)=HRAT0*HA1*TANH(0.0125*(WRAT2**0.42)/HA1)
           WVPD(I,J)=TRAT0*TA1*TANH(0.077*(WRAT2**0.25)/TA1)
           WVDR(I,J)=WINDIR
!
!  LIMIT WAVE HEIGHT TO BREAKING HEIGHT
!
            TBREAK=TIME*24.0
            HBREAK=0.78*DT(I,J)
            WVHT(I,J)=AMIN1(HBREAK,WVHT(I,J))
          ENDIF
 10   CONTINUE
!
!  OUTPUT WAVE PARAMETERS IF WAVEDYN IS SMBMODEL
!
      WRITE(1001)TIME
      WRITE(1001)WVHT
      WRITE(1001)WVPD
      WRITE(1001)WVDR
!
      RETURN
      END
