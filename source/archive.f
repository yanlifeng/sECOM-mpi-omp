      SUBROUTINE ARCHIVE
!     VERSION(08/07/91)
!
      INCLUDE 'comdeck'
!
!  FOR PARTICLE TRACKING
!
      REAL XOUTP(NPARTM),YOUTP(NPARTM),ZOUTP(NPARTM)
      REAL ENER_E(IM,JM),ENER_A(IM,JM),ENER_K(IM,JM)
!
      SAVE
!
!-----------------------------------------------------------------------
!-------- COMPUTATIONAL HISTORY WRITES AND ACCUMULATION ----------------
!-----------------------------------------------------------------------
      IF ((INTX.EQ.0).AND.(IAVGE.EQ.1)) GOTO 8350 ! NG03182010, NEED TO SKIP PRINTING OF COLD START IC

      DO 8300 JHIST=1,JHM
      IF(INTX.LT.IHIST(JHIST,1).OR.INTX.GT.IHIST(JHIST,2)) GOTO 8300
      IF ((INTX.EQ.0).AND.(JHIST.NE.1)) GOTO 8300 ! NG03182010, FOR COLD START PLT AVERAGING
      IF(TOR.EQ.'PROGNOSTIC' .OR. TOR.EQ.'DIAGNOSTIC'.OR.
     .   TOR.EQ.'TEMP_ONLY ' .OR. TOR.EQ.'SALT_ONLY ') THEN
!$doacross local(i,j,k)
       DO 500 K=1,KB
       DO 500 J=2,JMM1
       DO 500 I=2,IMM1
       ARCU (I,J,K)=ARCU (I,J,K)+U (I,J,K)*DEI
       ARCV (I,J,K)=ARCV (I,J,K)+V (I,J,K)*DEI
!NG03312009       ARCUX(I,J,K)=ARCUX(I,J,K)+U (I,J,K)*DEI*
!NG03312009     2    0.5*(H(I,J)+ET(I,J)+H(I-1,J)+ET(I-1,J))*DZ(K)
!NG03312009       ARCVX(I,J,K)=ARCVX(I,J,K)+V (I,J,K)*DEI*
!NG03312009     2    0.5*(H(I,J)+ET(I,J)+H(I,J-1)+ET(I,J-1))*DZ(K)
       ARCUX(I,J,K)=ARCUX(I,J,K)+U (I,J,K)*DEI*(H2(I,J)+H2(I-1,J))*
     2    0.25*(H(I,J)+ET(I,J)+H(I-1,J)+ET(I-1,J))*DZ(K)
       ARCVX(I,J,K)=ARCVX(I,J,K)+V (I,J,K)*DEI*(H1(I,J)+H1(I,J-1))*
     2    0.25*(H(I,J)+ET(I,J)+H(I,J-1)+ET(I,J-1))*DZ(K)
       ARCS (I,J,K)=ARCS (I,J,K)+S (I,J,K)*DEI
       ARCT (I,J,K)=ARCT (I,J,K)+T (I,J,K)*DEI
       ARCW (I,J,K)=ARCW (I,J,K)+WR(I,J,K)*DEI
       ARCKH(I,J,K)=ARCKH(I,J,K)+KH(I,J,K)*DEI
       ARCKM(I,J,K)=ARCKM(I,J,K)+KM(I,J,K)*DEI
       ARCAM(I,J,K)=ARCAM(I,J,K)+AAM(I,J,K)*DEI
!
!  CONSERVATIVE TRACER
!
       ARCC1(I,J,K)=ARCC1(I,J,K)+CONC1(I,J,K)*DEI
!
!  FOR SEDIMENT TRANSPORT
!
        ARCWS1(I,J,K)=ARCWS1(I,J,K)+WSET1(I,J,K)*DEI
        ARCWS2(I,J,K)=ARCWS2(I,J,K)+WSET2(I,J,K)*DEI
        ARCSED1(I,J,K)=ARCSED1(I,J,K)+CSED1(I,J,K)*DEI
        ARCSED2(I,J,K)=ARCSED2(I,J,K)+CSED2(I,J,K)*DEI
        ARCTAU(I,J,K)=ARCTAU(I,J,K)+TAU(I,J,K)*DEI
!
!  FOR CHEMICAL TRANSPORT
!
        ARCHEM1(I,J,K)=ARCHEM1(I,J,K)+CHEM1(I,J,K)*DEI
        ARCHEM2(I,J,K)=ARCHEM2(I,J,K)+CHEM2(I,J,K)*DEI
!
 500   CONTINUE
       DO 520 J=1,JM
         DO 520 I=1,IM
           ARCET (I,J)=ARCET (I,J)+ET (I,J)*DEI
!NG04092008           ARCWU (I,J)=ARCWU (I,J)+WUSURF (I,J)*DEI
!NG04092008           ARCWV (I,J)=ARCWV (I,J)+WVSURF (I,J)*DEI
           ARCPATM (I,J)=ARCPATM (I,J)+PATM(I,J)/100*DEI  ! mbar
!NG 04092008 ADDED 2D WINDS AND MET
           ARCWSX2D (I,J)=ARCWSX2D (I,J)+WSX2D(I,J)*DEI
           ARCWSY2D (I,J)=ARCWSY2D (I,J)+WSY2D(I,J)*DEI
           ARCCLD (I,J)=ARCCLD (I,J)+CLDFRC2D(I,J)*DEI
           ARCAIRT (I,J)=ARCAIRT (I,J)+AIRTM2D(I,J)*DEI
           ARCRHUM (I,J)=ARCRHUM (I,J)+RELHU2D(I,J)*DEI
           ARCSWOBS (I,J)=ARCSWOBS (I,J)+SWOBS2D(I,J)*DEI
           ARCQPREC (I,J)=ARCQPREC (I,J)+QPREC2D(I,J)*1000.*3600.*DEI ! mm/hr
           ARCQEVAP (I,J)=ARCQEVAP (I,J)+QEVAP2D(I,J)*86400.*365.*DEI ! m/yr
!NG 04142008 INCLUDED WAVES
           ARCWVPD (I,J)=ARCWVPD (I,J)+WVPD(I,J)*DEI
           ARCWVHT (I,J)=ARCWVHT (I,J)+WVHT(I,J)*DEI
           ARCWVDR (I,J)=ARCWVDR (I,J)+WVDR(I,J)*DEI
 520   CONTINUE
!
        IF (SEDTRAN.EQ.'INCLUDE') THEN
          DO 521 J=1,JM
            DO 521 I=1,IM
              IF (IBMSK(I,J).EQ.0) THEN
                DO 1521 LL=1,LAYMAX
                  ARCTHIK(I,J)=ARCTHIK(I,J)+TSED(LL,I,J)*DEI
 1521           CONTINUE
              ELSE
                IF (IBMSK(I,J).EQ.1) THEN
                  ARCTHIK(I,J)=ARCTHIK(I,J)+BEDTH(1,I,J)*DEI
                ELSE
                  ARCTHIK(I,J)=0.0
                ENDIF
              ENDIF
 521      CONTINUE
!
!
          IF (CHEMTRAN.EQ.'INCLUDE') THEN
            DO 522 J=1,JM
              DO 522 I=1,IM
                DO 522 LL=1,NCHEMLAY
                    ARCPBED(LL,I,J)=ARCPBED(LL,I,J)
     +                                  +CBEDCHEM(LL,I,J)*DEI
 522        CONTINUE
          ENDIF
        ENDIF
!
      ELSE
! Distinct history cells; preserve the time/layer summation order.
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(J,I,LL)
       DO 530 J=1,JM
         DO 530 I=1,IM
!NG03312009           ARCUX(I,J,1)=ARCUX(I,J,1)+UAF(I,J)*DEI*
!NG03312009     2       0.5*(H(I,J)+EL(I,J)+H(I-1,J)+EL(I-1,J))
!NG03312009           ARCVX(I,J,1)=ARCVX(I,J,1)+VAF(I,J)*DEI*
!NG03312009     2       0.5*(H(I,J)+EL(I,J)+H(I,J-1)+EL(I,J-1))
           ARCAM(I,J,1)=ARCAM(I,J,1)+AAM2D(I,J)*DEI
           IF(I.GT.1) THEN
           ARCUX(I,J,1)=ARCUX(I,J,1)+UAF(I,J)*DEI*(H2(I,J)+H2(I-1,J))*
     2       0.25*(H(I,J)+EL(I,J)+H(I-1,J)+EL(I-1,J))
           ENDIF
           IF(J.GT.1) THEN
           ARCVX(I,J,1)=ARCVX(I,J,1)+VAF(I,J)*DEI*(H1(I,J)+H1(I,J-1))*
     2       0.25*(H(I,J)+EL(I,J)+H(I,J-1)+EL(I,J-1))
           ENDIF
           ARCU (I,J,1)=ARCU (I,J,1)+UAF(I,J)*DEI 
           ARCV (I,J,1)=ARCV (I,J,1)+VAF(I,J)*DEI 
           ARCET (I,J)=ARCET (I,J)+EL (I,J)*DEI
!NG04092008           ARCWU (I,J)=ARCWU (I,J)+WUSURF (I,J)*DEI
!NG04092008           ARCWV (I,J)=ARCWV (I,J)+WVSURF (I,J)*DEI
           ARCPATM (I,J)=ARCPATM (I,J)+PATM(I,J)/100*DEI ! mbar
!NG 04092008 ADDED 2D WINDS AND MET
           ARCWSX2D (I,J)=ARCWSX2D (I,J)+WSX2D(I,J)*DEI
           ARCWSY2D (I,J)=ARCWSY2D (I,J)+WSY2D(I,J)*DEI
           ARCCLD (I,J)=ARCCLD (I,J)+CLDFRC2D(I,J)*DEI
           ARCAIRT (I,J)=ARCAIRT (I,J)+AIRTM2D(I,J)*DEI
           ARCRHUM (I,J)=ARCRHUM (I,J)+RELHU2D(I,J)*DEI
           ARCSWOBS (I,J)=ARCSWOBS (I,J)+SWOBS2D(I,J)*DEI
           ARCQPREC (I,J)=ARCQPREC (I,J)+QPREC2D(I,J)*1000.*3600.*DEI ! mm/hr
           ARCQEVAP (I,J)=ARCQEVAP (I,J)+QEVAP2D(I,J)*86400.*365.*DEI ! m/yr
!NG 04142008 INCLUDED WAVES
           ARCWVPD (I,J)=ARCWVPD (I,J)+WVPD(I,J)*DEI
           ARCWVHT (I,J)=ARCWVHT (I,J)+WVHT(I,J)*DEI
           ARCWVDR (I,J)=ARCWVDR (I,J)+WVDR(I,J)*DEI
!
!  FOR SEDIMENT TRANSPORT
!
           ARCWS1(I,J,1)=ARCWS1(I,J,1)+WSET1(I,J,1)*DEI
           ARCWS2(I,J,1)=ARCWS2(I,J,1)+WSET2(I,J,1)*DEI
           ARCSED1(I,J,1)=ARCSED1(I,J,1)+CSED1(I,J,1)*DEI
           ARCSED2(I,J,1)=ARCSED2(I,J,1)+CSED2(I,J,1)*DEI
           ARCTAU(I,J,1)=ARCTAU(I,J,1)+TAU(I,J,1)*DEI
!
           IF (IBMSK(I,J).EQ.0) THEN
             DO 61 LL=1,LAYMAX
               ARCTHIK(I,J)=ARCTHIK(I,J)+TSED(LL,I,J)*DEI
 61          CONTINUE
           ELSE
             IF (IBMSK(I,J).EQ.1) THEN
               ARCTHIK(I,J)=ARCTHIK(I,J)+BEDTH(1,I,J)*DEI
             ELSE
               ARCTHIK(I,J)=0.0
             ENDIF
           ENDIF
!
!  FOR CHEMICAL TRANSPORT
!
           ARCHEM1(I,J,1)=ARCHEM1(I,J,1)+CHEM1(I,J,1)*DEI
           ARCHEM2(I,J,1)=ARCHEM2(I,J,1)+CHEM2(I,J,1)*DEI
           DO 62 LL=1,NCHEMLAY
             ARCPBED(LL,I,J)=ARCPBED(LL,I,J)
     +                           +CBEDCHEM(LL,I,J)*DEI
 62        CONTINUE
!
 530  CONTINUE
!$OMP END PARALLEL DO
      ENDIF
!
      IF(INTX.NE.IHIST(JHIST,2)) GO TO 560
!
!-------- WRITE CONSTANTS FIRST TIME THROUGH ---------------------------
      IF (CONSPLT.AND.(.NOT.CDFOUT)) THEN
       OPEN (IUPLT,FORM='unformatted',FILE='gcmplt')
!qa     + CONVERT='BIG_ENDIAN')
!
       WRITE(IUPLT) IM,JM,KB
!-------- NOT NECESSARY TO SPECIFY EBCM & QBCM -------------------------
!-------- PER JAH INSTRUCTIONS - LEAVE IN FOR NOW (06-28-91) -----------
       WRITE(IUPLT) EBCM,QBCM,NCHEMLAY
       WRITE(IUPLT) DTI,GRAV,UMOL,TOR,TRACER,SEDTRAN,CHEMTRAN,SEDTYPE
       WRITE(IUPLT) NUMEBC
       IF(NUMEBC.GT.0)
     . WRITE(IUPLT) (IETA(I),JETA(I),ICON(I),JCON(I),I=1,NUMEBC)
       WRITE(IUPLT) NUMQBC
       IF(NUMQBC.GT.0)
     . WRITE(IUPLT) (IQC(I),JQC(I),I=1,NUMQBC)
       WRITE(IUPLT) H
       WRITE(IUPLT) H1
       WRITE(IUPLT) H2
       WRITE(IUPLT) ANG
       WRITE(IUPLT) DUM
       WRITE(IUPLT) DVM
       WRITE(IUPLT) FSM
       CONSPLT=.FALSE.
      ENDIF
!
!  CONVERT THICKNESS FROM g/cm**2 TO cm
!
            IF (SEDTRAN.EQ.'INCLUDE') THEN
!
              DO 650 J=1,JM
                DO 650 I=1,IM
                  IF (FSM(I,J).GT.0.0) THEN
                    IF (IBMSK(I,J).EQ.0) THEN
                      TSETOT=0.0
                      DO 640 LL=1,LAYMAX
                        TSETOT=TSETOT+TSED0(LL,I,J)
 640                  CONTINUE
!
                      ARCTHIK(I,J)=ARCTHIK(I,J)-TSETOT
!
!  CONVERT FROM g/cm**2 TO cm
!
                      ARCTHIK(I,J)=ARCTHIK(I,J)/CBED(I,J)
!
                      IF (SEDTYPE.EQ.'SAND') ARCTHIK(I,J)=0.0
                    ELSE
                      IF (IBMSK(I,J).EQ.1) THEN
!
!  IN cm
!
                         ARCTHIK(I,J)=100.*(ARCTHIK(I,J)/
     +                          ((CBED(I,J)/2.65)*
     +                          H1(I,J)*H2(I,J))-BEDTHI)
                         IF (SEDTYPE.EQ.'MUD ') ARCTHIK(I,J)=0.0
                      ELSE
                         ARCTHIK(I,J)=0.0
                      ENDIF
                    ENDIF
                  ENDIF
!
!  CONVERT FROM g/cm**3 TO mg/l
!
                  DO 660 K=1,KBM1
                    ARCSED1(I,J,K)=1000000.*ARCSED1(I,J,K)
                    ARCSED2(I,J,K)=1000000.*ARCSED2(I,J,K)
 660              CONTINUE
 650          CONTINUE
            ENDIF
!
            IF (CHEMTRAN.EQ.'INCLUDE') THEN
!
!  CONVERT FROM ug CHEM/cm**3 TO ug CHEM/l
!
              DO 1660 K=1,KBM1
                DO 1660 J=1,JM
                  DO 1660 I=1,IM
                    ARCHEM1(I,J,K)=1000.*ARCHEM1(I,J,K)               
                    ARCHEM2(I,J,K)=1000.*ARCHEM2(I,J,K)               
 1660         CONTINUE
!
           ENDIF
!
      TMIDDLE=TIME-(.5*DTI*DAYI/DEI)
      IF (CDFOUT) THEN 
       CALL PUTCDF(TMIDDLE)
       icdfcheck=1
       GO TO 5600
      END IF

      WRITE(IUPLT) TMIDDLE
      WRITE(IUPLT) ARCET
!
      IF(TOR.EQ.'BAROTROPIC') THEN
       WRITE(IUPLT) ((ARCU  (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT) ((ARCV  (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT) ((ARCUX (I,J,1),I=1,IM),J=1,JM)
       WRITE(IUPLT) ((ARCVX (I,J,1),I=1,IM),J=1,JM)
!
       IF (TRACER.EQ.'INCLUDE') 
     +             WRITE (IUPLT) ((ARCC1(I,J,1),I=1,IM),J=1,JM)
!
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
           WRITE (IUPLT) ((ARCWS1(I,J,1),I=1,IM),J=1,JM)
           WRITE (IUPLT) ((ARCSED1(I,J,1),I=1,IM),J=1,JM)
         END IF
         IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
           WRITE (IUPLT) ((ARCWS2(I,J,1),I=1,IM),J=1,JM)
           WRITE (IUPLT) ((ARCSED2(I,J,1),I=1,IM),J=1,JM)
         END IF
         WRITE (IUPLT) ((ARCTAU(I,J,1),I=1,IM),J=1,JM)
!
!  TRANSFER FROM ARCHTHIK (real*8) TO ARCET (real*4) FOR OUTPUT
!
          DO 1665 J=1,JM
            DO 1665 I=1,IM
              ARCET(I,J)=ARCTHIK(I,J)
 1665     CONTINUE
!
         WRITE (IUPLT) ((ARCET(I,J),I=1,IM),J=1,JM)
       ENDIF
!
       IF (CHEMTRAN.EQ.'INCLUDE') THEN
         WRITE (IUPLT) ((ARCHEM1(I,J,1),I=1,IM),J=1,JM)
         WRITE (IUPLT) ((ARCHEM2(I,J,1),I=1,IM),J=1,JM)
         WRITE (IUPLT) (((ARCPBED(N,I,J),I=1,IM),J=1,JM),
     +                      N=1,NCHEMLAY)
       ENDIF
!
      ELSE
       WRITE(IUPLT) Z
       WRITE(IUPLT) ZZ
       WRITE(IUPLT) DZ
       WRITE(IUPLT) ARCU
       WRITE(IUPLT) ARCV
       WRITE(IUPLT) ARCUX
       WRITE(IUPLT) ARCVX
       WRITE(IUPLT) ARCT
       WRITE(IUPLT) ARCS
       WRITE(IUPLT) ARCW
       WRITE(IUPLT) ARCKH
!
!  FOR CONSERVATIVE TRACER
!
       IF (TRACER.EQ.'INCLUDE') WRITE(IUPLT) ARCC1
!
        IF (SEDTRAN.EQ.'INCLUDE') THEN
         IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
          WRITE (IUPLT) ARCWS1
          WRITE (IUPLT) ARCSED1
         END IF
         IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
          WRITE (IUPLT) ARCWS2
          WRITE (IUPLT) ARCSED2
         END IF
         WRITE (IUPLT) ARCTAU
!
!  TRANSFER FROM ARCHTHIK (real*8) TO ARCET (real*4) FOR OUTPUT
!
          DO 1670 J=1,JM
            DO 1670 I=1,IM
              ARCET(I,J)=ARCTHIK(I,J)
 1670     CONTINUE
!
          WRITE (IUPLT) ARCET
          WRITE (IUPLT) ARCTAU
        ENDIF
!
        IF (CHEMTRAN.EQ.'INCLUDE') THEN
          WRITE (IUPLT) ARCHEM1
          WRITE (IUPLT) ARCHEM2
          WRITE (IUPLT) (((ARCPBED(N,I,J),I=1,IM),J=1,JM),
     +                      N=1,NCHEMLAY)
        ENDIF
!
      ENDIF
!
 5600 CONTINUE
      IF(PARTICLE.NE.'INCLUDE')GO TO 5700
!
!*************************************************************
!
!  PARTICLE TRACKING OUTPUT
!
!---> convert to physical domain
!
        CPI=3.141593/180.
!
        do LL=1,nsource
          do NN=1,ngradeloop
!
!  ICHK = 0 :  NO PARTICLES IN THIS CLASS - NO OUTPUT
!  ICHK = 1 :  PARTICLES IN THIS CLASS - OUTPUT
!
            ICHK=0
!
            do MM=1,NPART
!
              XOUTP(MM)=0.0
              YOUTP(MM)=0.0
              ZOUTP(MM)=0.0
!
              iout=inout(LL,MM,NN)
!             if (iout.eq.1) goto 4444                       !n.kim 4/21/98
              if (iout.eq.1)then                             !n.kim 4/21/98
               write(*,*)'skipping 4444 jump',TIME,LL,MM,NN  !n.kim 4/21/98
               XOUTP(MM)=-999.0    !locations are reset for outgoing
               YOUTP(MM)=-999.0    !particles (hli,02/14/02)
               ZOUTP(MM)=-999.0
               goto 4444                                     !hli   2/14/02
              end if                                         !n.kim 4/21/98

!
              ICHK=1
!
              xstart=XP(LL,MM,NN)
              ystart=YP(LL,MM,NN)
              zstart=ZP(LL,MM,NN)
!
              ic=xstart
              jc=ystart
!
              do k=1,kbm1
                if (zstart.le.z(k).and.zstart.ge.z(k+1)) goto 4445
              enddo
              WRITE (*,*) 'trouble !',zstart,INTX
 4445         kc=k        
!
              xc= ic+0.5
              yc= jc+0.5
              zc= zz(kc)
!
!---> interpolate D
!     find i-range for D
!
              i0h   = ic
              xlocal= xstart-xc
              xdist = abs(xlocal)
              if(xlocal.gt.0.) i1h=i0h+1
              if(xlocal.le.0.) i1h=i0h-1
!
!     use value at one element
!
              if(fsm(i1h,jc).eq.0.)i1h=i0h
!
!     find j-range for D
!
              j0h   = jc
              ylocal= ystart-yc
              ydist = abs(ylocal)
              if(ylocal.gt.0.) j1h=j0h+1
              if(ylocal.le.0.) j1h=j0h-1
!
!     use value at one element
!
              if(fsm(ic,j1h).eq.0.)j1h=j0h
!
!xyz  important change
!
              if(fsm(i1h,j1h).eq.0) then
                i1h=i0h
                j1h=j0h
              endif
!
!     interpolate D
!
              Dj0= D(i0h,j0h) + xdist * (D(i1h,j0h)-D(i0h,j0h))
              Dj1= D(i0h,j1h) + xdist * (D(i1h,j1h)-D(i0h,j1h))
              Dp = Dj0        + ydist * (Dj1-Dj0)
!
!     interpolate EL
!
              ELj0= EL(i0h,j0h) + xdist * (EL(i1h,j0h)-EL(i0h,j0h))
              ELj1= EL(i0h,j1h) + xdist * (EL(i1h,j1h)-EL(i0h,j1h))
              ELp = ELj0        + ydist * (ELj1-ELj0)
!      
              XOUTP(MM)=xcor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +xcor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +xcor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +xcor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
!
              YOUTP(MM)=ycor(ic,jc)     * (0.5-xlocal)*(0.5-ylocal)
     +         +ycor(ic+1,jc)   * (0.5+xlocal)*(0.5-ylocal) 
     +         +ycor(ic+1,jc+1) * (0.5+xlocal)*(0.5+ylocal)    
     +         +ycor(ic,jc+1)   * (0.5-xlocal)*(0.5+ylocal)
!
              ZOUTP(MM)=ELp + zstart * Dp
!
 4444         continue
             enddo
!
!  OUTPUT PARTICLE POSITIONS FOR PARTICLES IN CLASS nn
!
!               IF (ICHK.EQ.0) GOTO 900     !! Original Code (n.kim commented out)
!
               WRITE (32,930)LL,NN,NPART,TIME   
 930           FORMAT (3I8,F10.4)
!
               WRITE (32,940) (XOUTP(NP),NP=1,NPART)
               WRITE (32,940) (YOUTP(NP),NP=1,NPART)
               WRITE (32,940) (ZOUTP(NP),NP=1,NPART)
 940           FORMAT (20F10.1)
!hli 900           CONTINUE
           enddo
         enddo
!
!*****************************************************************
!
 5700 CONTINUE
!
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I,J,K,LL)
!$OMP DO SCHEDULE(STATIC)
      DO 570 J=1,JM
        DO 570 I=1,IM
          ARCET (I,J)=0.0
!NG04092008          ARCWU (I,J)=0.0
!NG04092008          ARCWV (I,J)=0.0
          ARCPATM (I,J)=0.0
          ARCTHIK(I,J) = 0.0
!NG 04092008 ADDED 2D WINDS AND MET
          ARCWSX2D (I,J)=0.0
          ARCWSY2D (I,J)=0.0
          ARCCLD (I,J)=0.0
          ARCAIRT (I,J)=0.0
          ARCRHUM (I,J)=0.0
          ARCSWOBS (I,J)=0.0
          ARCQPREC (I,J)=0.0
          ARCQEVAP (I,J)=0.0
!NG 04142008 INCLUDED WAVES
          ARCWVPD (I,J)=0.0
          ARCWVHT (I,J)=0.0
          ARCWVDR (I,J)=0.0
          DO 570 LL=1,NCHEMLAY
            ARCPBED(LL,I,J) = 0.0
 570  CONTINUE
!$OMP END DO NOWAIT
!
!$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
      DO 580 K=1,KB
      DO 580 J=1,JM
      DO 580 I=1,IM
      ARCU (I,J,K)=0.0
      ARCV (I,J,K)=0.0
      ARCUX(I,J,K)=0.0
      ARCVX(I,J,K)=0.0
      ARCT (I,J,K)=0.0
      ARCS (I,J,K)=0.0
      ARCW (I,J,K)=0.0
      ARCKH(I,J,K)=0.0
      ARCKM(I,J,K)=0.0
      ARCAM(I,J,K)=0.0
!
      ARCC1(I,J,K)=0.0
!
      ARCWS1(I,J,K) = 0.0
      ARCWS2(I,J,K) = 0.0
      ARCSED1(I,J,K) = 0.0
      ARCSED2(I,J,K) = 0.0
      ARCTAU(I,J,K)  = 0.0
!
      ARCHEM1(I,J,K) = 0.0
      ARCHEM2(I,J,K) = 0.0
 580  CONTINUE
!$OMP END DO
!$OMP END PARALLEL
 560  CONTINUE
!
 8300 CONTINUE


 8350 CONTINUE ! NG03182010, NEED TO SKIP PRINTING OF COLD START IC

!
!-----------------------------------------------------------------------
!------- TIME SERIES WRITES AND ACCUMULATIONS --------------------------
      IF ((INTX.EQ.0).AND.(ISKILL.EQ.1)) GOTO 8450 ! NG03182010, NEED TO SKIP PRINTING OF COLD START IC

      IF(TOR.EQ.'BAROTROPIC') THEN
       DO 140 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
       ESAVE(N)=ESAVE(N)+EL(II,JJ)*SKILLI
!NG04092008       WUSAVE(N)=WUSAVE(N)+WUSURF(II,JJ)*SKILLI
!NG04092008       WVSAVE(N)=WVSAVE(N)+WVSURF(II,JJ)*SKILLI
       TXSAVE(N)=TXSAVE(N)+TXSURF(II,JJ)*SKILLI
       TYSAVE(N)=TYSAVE(N)+TYSURF(II,JJ)*SKILLI
       PATMSAVE(N)=PATMSAVE(N)+PATM(II,JJ)/100*SKILLI ! mbar
!NG 04092008 ADDED 2D WINDS AND MET
       WSX2DSAVE(N)=WSX2DSAVE(N)+WSX2D(II,JJ)*SKILLI
       WSY2DSAVE(N)=WSY2DSAVE(N)+WSY2D(II,JJ)*SKILLI
       CLDSAVE(N)=CLDSAVE(N)+CLDFRC2D(II,JJ)*SKILLI
       AIRTSAVE(N)=AIRTSAVE(N)+AIRTM2D(II,JJ)*SKILLI
       RHUMSAVE(N)=RHUMSAVE(N)+RELHU2D(II,JJ)*SKILLI
       SWOBSSAVE(N)=SWOBSSAVE(N)+SWOBS2D(II,JJ)*SKILLI
       QPRECSAVE(N)=QPRECSAVE(N)+QPREC2D(II,JJ)*1000.*3600.*SKILLI ! mm/hr
       QEVAPSAVE(N)=QEVAPSAVE(N)+QEVAP2D(II,JJ)*86400.*365.*SKILLI ! m/yr
!NG 04142008 INCLUDED WAVES
       WVPDSAVE(N)=WVPDSAVE(N)+WVPD(II,JJ)*SKILLI
       WVHTSAVE(N)=WVHTSAVE(N)+WVHT(II,JJ)*SKILLI
       WVDRSAVE(N)=WVDRSAVE(N)+WVDR(II,JJ)*SKILLI
       

140    CONTINUE
       DO 150 N=1,VPTS
         II=INXIV(N)
         JJ=INXJV(N)
         UZSAVE(N,1)=UZSAVE(N,1)+.5*(UA(II,JJ)+UA(II+1,JJ))*SKILLI
         VZSAVE(N,1)=VZSAVE(N,1)+.5*(VA(II,JJ)+VA(II,JJ+1))*SKILLI
         AMZSAVE(N,1)=AMZSAVE(N,1)+AAM2D(II,JJ)*SKILLI
!
!  FOR SEDIMENT TRANSPORT
!
         WS1SAVE(N,1)=WS1SAVE(N,1)+WSET1(II,JJ,1)*SKILLI
         WS2SAVE(N,1)=WS2SAVE(N,1)+WSET2(II,JJ,1)*SKILLI
         C1SAVE(N,1)=C1SAVE(N,1)+CSED1(II,JJ,1)*SKILLI
         C2SAVE(N,1)=C2SAVE(N,1)+CSED2(II,JJ,1)*SKILLI
         TAUSAVE(N,1)=TAUSAVE(N,1)+TAU(II,JJ,1)*SKILLI
!
         IF (IBMSK(II,JJ).EQ.0) THEN
           DO 151 LL=1,LAYMAX
             THSAVE(N)=THSAVE(N)+TSED(LL,II,JJ)*SKILLI
 151       CONTINUE
         ELSE
           IF (IBMSK(II,JJ).EQ.1) THEN
             THSAVE(N)=THSAVE(N)+BEDTH(1,II,JJ)*SKILLI
           ELSE
             THSAVE(N)=0.0
           ENDIF
         ENDIF
!
!  FOR CHEMICAL TRANSPORT
!
          P1SAVE(N,1)=P1SAVE(N,1)+CHEM1(II,JJ,1)*SKILLI
          P2SAVE(N,1)=P2SAVE(N,1)+CHEM2(II,JJ,1)*SKILLI
          DO 152 LL=1,NCHEMLAY
            PBEDSAVE(N,LL)=PBEDSAVE(N,LL)+CBEDCHEM(LL,I,J)*SKILLI
 152      CONTINUE
!
150    CONTINUE
!
!-------- COMPUTE CROSS SECTIONAL FLUXES -------------------------------
       DO 180 N=1,FPTS
       IS=ISFLX(N)
       JS=JSFLX(N)
       IF(DIRFLX(N).EQ.'IDIR') THEN
!       IE=IS+NFLXE(N)-1
        DO 181 I=IS,IS+NFLXE(N)-1
!PL 04/09/2002 MODIFY TO INCLUDE SEDIMENT CROSS-SECTION FLUX - NEXT 14 LINES
!PL
!PL 181    CCFLUX(N,1)=CCFLUX(N,1)+VAF(I,JS)*SKILLI*0.25*
!PL     2     (H(I,JS)+EL(I,JS)+H(I,JS-1)+EL(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
        CCFLUX(N,1)=CCFLUX(N,1)+VAF(I,JS)*SKILLI*0.25*
     2   (H(I,JS)+EL(I,JS)+H(I,JS-1)+EL(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
        IF (SEDTRAN.EQ.'INCLUDE') THEN
        CS1FLUX(N,1)=CS1FLUX(N,1)+VAF(I,JS)*SKILLI*0.25*
     +   (H(I,JS)+EL(I,JS)+H(I,JS-1)+EL(I,JS-1))*(H1(I,JS)+H1(I,JS-1))*
     +   (CSED1(I,JS,1)+CSED1(I,JS-1,1))*0.5*1000.0 !IN kg/s 
        CS2FLUX(N,1)=CS2FLUX(N,1)+VAF(I,JS)*SKILLI*0.25*
     +   (H(I,JS)+EL(I,JS)+H(I,JS-1)+EL(I,JS-1))*(H1(I,JS)+H1(I,JS-1))*
     +   (CSED2(I,JS,1)+CSED2(I,JS-1,1))*0.5*1000.0 !IN kg/s 
        END IF
 181    CONTINUE
       ELSE
!        DIRFLX(N).EQ.'JDIR'
!       JE=JS+NFLXE(N)-1
        DO 182 J=JS,JS+NFLXE(N)-1
!PL 04/09/2002 MODIFY TO INCLUDE SEDIMENT CROSS-SECTION FLUX - NEXT 14 LINES
!PL
!PL 182    CCFLUX(N,1)=CCFLUX(N,1)+UAF(IS,J)*SKILLI*0.25*
!PL     2     (H(IS,J)+EL(IS,J)+H(IS-1,J)+EL(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
        CCFLUX(N,1)=CCFLUX(N,1)+UAF(IS,J)*SKILLI*0.25*
     2   (H(IS,J)+EL(IS,J)+H(IS-1,J)+EL(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
        IF (SEDTRAN.EQ.'INCLUDE') THEN
        CS1FLUX(N,1)=CS1FLUX(N,1)+UAF(IS,J)*SKILLI*0.25*
     +   (H(IS,J)+EL(IS,J)+H(IS-1,J)+EL(IS-1,J))*(H2(IS,J)+H2(IS-1,J))*
     +   (CSED1(IS,J,1)+CSED1(IS-1,J,1))*0.5*1000.0 !IN kg/s 
        CS2FLUX(N,1)=CS2FLUX(N,1)+UAF(IS,J)*SKILLI*0.25*
     +   (H(IS,J)+EL(IS,J)+H(IS-1,J)+EL(IS-1,J))*(H1(IS,J)+H1(IS-1,J))*
     +   (CSED2(IS,J,1)+CSED2(IS-1,J,1))*0.5*1000.0 !IN kg/s 
        END IF
 182    CONTINUE
       ENDIF
 180   CONTINUE
!
!---- TOR = PROGNOSTIC OR DIAGNOSTIC OR TEMP_ONLY OR SALT_ONLY ----
      ELSE
       DO 160 N=1,EPTS
       II=INXIE(N)
       JJ=INXJE(N)
       ESAVE(N)=ESAVE(N)+ET(II,JJ)*SKILLI
!NG04092008       WUSAVE(N)=WUSAVE(N)+WUSURF(II,JJ)*SKILLI
!NG04092008       WVSAVE(N)=WVSAVE(N)+WVSURF(II,JJ)*SKILLI
       TXSAVE(N)=TXSAVE(N)+TXSURF(II,JJ)*SKILLI
       TYSAVE(N)=TYSAVE(N)+TYSURF(II,JJ)*SKILLI
       PATMSAVE(N)=PATMSAVE(N)+PATM(II,JJ)/100*SKILLI ! mbar
!NG 04092008 ADDED 2D WINDS AND MET
       WSX2DSAVE(N)=WSX2DSAVE(N)+WSX2D(II,JJ)*SKILLI
       WSY2DSAVE(N)=WSY2DSAVE(N)+WSY2D(II,JJ)*SKILLI
       CLDSAVE(N)=CLDSAVE(N)+CLDFRC2D(II,JJ)*SKILLI
       AIRTSAVE(N)=AIRTSAVE(N)+AIRTM2D(II,JJ)*SKILLI
       RHUMSAVE(N)=RHUMSAVE(N)+RELHU2D(II,JJ)*SKILLI
       SWOBSSAVE(N)=SWOBSSAVE(N)+SWOBS2D(II,JJ)*SKILLI
       QPRECSAVE(N)=QPRECSAVE(N)+QPREC2D(II,JJ)*1000.*3600.*SKILLI ! mm/hr
       QEVAPSAVE(N)=QEVAPSAVE(N)+QEVAP2D(II,JJ)*86400.*365.*SKILLI ! m/yr
!NG 04142008 INCLUDED WAVES
       WVPDSAVE(N)=WVPDSAVE(N)+WVPD(II,JJ)*SKILLI
       WVHTSAVE(N)=WVHTSAVE(N)+WVHT(II,JJ)*SKILLI
       WVDRSAVE(N)=WVDRSAVE(N)+WVDR(II,JJ)*SKILLI
160    CONTINUE
       DO 170 N=1,VPTS
       II=INXIV(N)
       JJ=INXJV(N)
       DZSAVE(N)=DZSAVE(N)+DT(II,JJ)*SKILLI
!
!  FOR SEDIMENT TRANSPORT
!
        IF (IBMSK(II,JJ).EQ.0) THEN
          DO 201 LL=1,LAYMAX
            THSAVE(N)=THSAVE(N)+TSED(LL,II,JJ)*SKILLI
 201      CONTINUE
        ELSE 
          IF (IBMSK(II,JJ).EQ.1) THEN
            THSAVE(N)=THSAVE(N)+BEDTH(1,II,JJ)*SKILLI
          ELSE
            THSAVE(N)=0.0
          ENDIF
        ENDIF
!
!  FOR CHEMICAL TRANSPORT
!
          DO 202 LL=1,NCHEMLAY
            PBEDSAVE(N,LL)=PBEDSAVE(N,LL)+CBEDCHEM(LL,II,JJ)*SKILLI
 202      CONTINUE
!
       DO 170 K=1,KB
       UZSAVE(N,K)=UZSAVE(N,K)+.5*(U(II,JJ,K)+U(II+1,JJ,K))*SKILLI
       VZSAVE(N,K)=VZSAVE(N,K)+.5*(V(II,JJ,K)+V(II,JJ+1,K))*SKILLI
       WZSAVE(N,K)=WZSAVE(N,K)+WR(II,JJ,K)*SKILLI
       SZSAVE(N,K)=SZSAVE(N,K)+S(II,JJ,K)*SKILLI
       TZSAVE(N,K)=TZSAVE(N,K)+T(II,JJ,K)*SKILLI
       KHZSAVE(N,K)=KHZSAVE(N,K)+KH(II,JJ,K)*SKILLI
       KMZSAVE(N,K)=KMZSAVE(N,K)+KM(II,JJ,K)*SKILLI
       AMZSAVE(N,K)=AMZSAVE(N,K)+AAM(II,JJ,K)*SKILLI
!
!  FOR CONSERVATIVE TRACER
!
       C1ZSAVE(N,K)=C1ZSAVE(N,K)+CONC1(II,JJ,K)*SKILLI
!
!  FOR SEDIMENT TRANSPORT
!
         WS1SAVE(N,K)=WS1SAVE(N,K)+WSET1(II,JJ,K)*SKILLI
         WS2SAVE(N,K)=WS2SAVE(N,K)+WSET2(II,JJ,K)*SKILLI
         C1SAVE(N,K)=C1SAVE(N,K)+CSED1(II,JJ,K)*SKILLI
         C2SAVE(N,K)=C2SAVE(N,K)+CSED2(II,JJ,K)*SKILLI
         TAUSAVE(N,K)=TAUSAVE(N,K)+TAU(II,JJ,K)*SKILLI
!
!  FOR CHEMICAL TRANSPORT
!
         P1SAVE(N,K)=P1SAVE(N,K)+CHEM1(II,JJ,K)*SKILLI
         P2SAVE(N,K)=P2SAVE(N,K)+CHEM2(II,JJ,K)*SKILLI
!
 170   CONTINUE
!
!-------- COMPUTE CROSS SECTIONAL FLUXES -------------------------------
       DO 250 N=1,FPTS
       IS=ISFLX(N)
       JS=JSFLX(N)
       IF(DIRFLX(N).EQ.'IDIR') THEN
!       IE=IS+NFLXE(N)-1
        DO 251 K=1,KBM1
        DO 251 I=IS,IS+NFLXE(N)-1
!PL 04/09/2002 MODIFY TO INCLUDE SEDIMENT CROSS-SECTION FLUX - NEXT 14 LINES
!PL
!PL 251    CCFLUX(N,K)=CCFLUX(N,K)+V(I,JS,K)*DZ(K)*SKILLI*0.25*
!PL     2     (H(I,JS)+ET(I,JS)+H(I,JS-1)+ET(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
        CCFLUX(N,K)=CCFLUX(N,K)+V(I,JS,K)*DZ(K)*SKILLI*0.25*
     2   (H(I,JS)+ET(I,JS)+H(I,JS-1)+ET(I,JS-1))*(H1(I,JS)+H1(I,JS-1))
        IF (SEDTRAN.EQ.'INCLUDE') THEN
        CS1FLUX(N,K)=CS1FLUX(N,K)+V(I,JS,K)*DZ(K)*SKILLI*0.25*
     +   (H(I,JS)+ET(I,JS)+H(I,JS-1)+ET(I,JS-1))*(H1(I,JS)+H1(I,JS-1))*
     +   (CSED1(I,JS,K)+CSED1(I,JS-1,K))*0.5*1000.0 !IN kg/s 
        CS2FLUX(N,K)=CS2FLUX(N,K)+V(I,JS,K)*DZ(K)*SKILLI*0.25*
     +   (H(I,JS)+ET(I,JS)+H(I,JS-1)+ET(I,JS-1))*(H1(I,JS)+H1(I,JS-1))*
     +   (CSED2(I,JS,K)+CSED2(I,JS-1,K))*0.5*1000.0 !IN kg/s 
        END IF
 251    CONTINUE
       ELSE
!        DIRFLX(N).EQ.'JDIR'
!       JE=JS+NFLXE(N)-1
        DO 252 K=1,KBM1
        DO 252 J=JS,JS+NFLXE(N)-1
!PL 04/09/2002 MODIFY TO INCLUDE SEDIMENT CROSS-SECTION FLUX - NEXT 14 LINES
!PL
!PL 252    CCFLUX(N,K)=CCFLUX(N,K)+U(IS,J,K)*DZ(K)*SKILLI*0.25*
!PL     2     (H(IS,J)+ET(IS,J)+H(IS-1,J)+ET(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
        CCFLUX(N,K)=CCFLUX(N,K)+U(IS,J,K)*DZ(K)*SKILLI*0.25*
     2   (H(IS,J)+ET(IS,J)+H(IS-1,J)+ET(IS-1,J))*(H2(IS,J)+H2(IS-1,J))
        IF (SEDTRAN.EQ.'INCLUDE') THEN
        CS1FLUX(N,K)=CS1FLUX(N,K)+U(IS,J,K)*DZ(K)*SKILLI*0.25*
     +   (H(IS,J)+ET(IS,J)+H(IS-1,J)+ET(IS-1,J))*(H2(IS,J)+H2(IS-1,J))*
     +   (CSED1(IS,J,K)+CSED1(IS-1,J,K))*0.5*1000.0 !IN kg/s 
        CS2FLUX(N,K)=CS2FLUX(N,K)+U(IS,J,K)*DZ(K)*SKILLI*0.25*
     +   (H(IS,J)+ET(IS,J)+H(IS-1,J)+ET(IS-1,J))*(H2(IS,J)+H2(IS-1,J))*
     +   (CSED2(IS,J,K)+CSED2(IS-1,J,K))*0.5*1000.0 !IN kg/s 
        END IF
 252    CONTINUE
       ENDIF
 250   CONTINUE
      ENDIF
!
      IF(TOR.EQ.'BAROTROPIC') THEN
!$doacross local(i,j),reduction(esum,ape,tke)
!jah   initialize variables
       ESUM=0.
       APE=0.
       TKE=0.
       IF(INTX.GT.0) THEN
       IF(.NOT.(AREA.GT.0.)) STOP 95
!      These three LOCAL diagnostics are reset on every ARCHIVE call;
!      they are not the running station/history averages. Their only
!      consumer in this branch is the skill-record writer below.
!      Keep the input/area guard above on every step. All sampled sums
!      retain both the arithmetic and original serial fold order.
       IF(ISKILL.EQ.0) GOTO 7100
       IF(MOD(INTX,ISKILL).NE.0) GOTO 7100
!      Evaluate independent contributions in parallel. Fold original
!      single-precision accumulators in the original i-fastest order.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J) SCHEDULE(STATIC)
       DO J=1,JM
        DO I=1,IM
         ENER_E(I,J)=EL(I,J)*ART(I,J)*FSM(I,J)/AREA*SKILLI
         ENER_A(I,J)=.5*GRAV*EL(I,J)*EL(I,J)*FSM(I,J)*ART(I,J)/
     2      AREA*SKILLI
         IF(I.LT.IM.AND.J.LT.JM)
     1    ENER_K(I,J)=0.125*D(I,J)*((UA(I,J)+UA(I+1,J))**2+
     2     (VA(I,J)+VA(I,J+1))**2)*FSM(I,J)*ART(I,J)/
     3     AREA*SKILLI
        ENDDO
       ENDDO
!$OMP END PARALLEL DO
       DO 100 J=1,JM
       DO 100 I=1,IM
       ESUM=ESUM+ENER_E(I,J)
       APE=APE+ENER_A(I,J)
       IF(I.LT.IM.AND.J.LT.JM) THEN
       TKE=TKE+ENER_K(I,J)
       ENDIF
 100   CONTINUE
 7100  CONTINUE
       ENDIF
      ELSE
!$doacross local(i,j),reduction(volume)
!jah   INITIALIZE VARIABLES
       VOLUME=0.
       ESUM=0.
       APE=0.
       EM=0.
       DO 200 J=1,JM
       DO 200 I=1,IM
       VOLUME=VOLUME+DT(I,J)*ART(I,J)*FSM(I,J)*SKILLI
  200  ESUM=ESUM+ET(I,J)*ART(I,J)*FSM(I,J)*SKILLI
!$doacross local(i,j,k),reduction(tsum,ssum)
!jah   INITIALIZE VARIABLES
       TSUM=0.
       SSUM=0.
       DO 210 K=1,KBM1
       DO 210 J=1,JM
       DO 210 I=1,IM
       TSUM=TSUM+T(I,J,K)*DT(I,J)*DZ(K)*ART(I,J)*FSM(I,J)*SKILLI
  210  SSUM=SSUM+S(I,J,K)*DT(I,J)*DZ(K)*ART(I,J)*FSM(I,J)*SKILLI
!$doacross local(i,j,k,trho,vol),reduction(em,ape)
       DO 220 K=1,KBM1
       DO 220 J=1,JM
       DO 220 I=1,IM
       TRHO=(RHO(I,J,K)+1.)*1000.
       VOL=DT(I,J)*ART(I,J)*DZ(K)
       EM=EM+TRHO*VOL*FSM(I,J)
  220  APE=APE+GRAV*TRHO*ZZ(K)*DT(I,J)*VOL*FSM(I,J)
      ENDIF
!
      TMIDDLE=TIME-(0.5*DTI*DAYI/SKILLI)
      IF(ISKILL.EQ.0) GO TO 8450
      IF(MOD(INTX,ISKILL).NE.0.OR.INTX.EQ.0) GO TO 8450
!NG02192010      IF (CDFOUT) THEN
      IF ((CDFOUT).AND.((EPTS+VPTS+FPTS).NE.0)) THEN ! NG02192010 DO NOT CREATE tsepic.cdf IF NONE REQUESTED
!         print *, 'tmiddle:', tmiddle
       CALL TSCDF(TMIDDLE,icdfcheck)
       GO TO 2300
      ENDIF
!
!-------- WRITE CONSTANTS THE FIRST TIME THROUGH -----------------------
      IF (CONSTSR) THEN
       OPEN (IUTSR,FORM='unformatted',FILE='gcmtsr')
!qa     + CONVERT='BIG_ENDIAN')
       WRITE(IUTSR) TOR,TRACER,SEDTRAN,CHEMTRAN,SEDTYPE
       WRITE(IUTSR) KBM1,NCHEMLAY
       WRITE(IUTSR) EPTS
       IF(EPTS.NE.0)
     . WRITE(IUTSR) (INXIE(N),INXJE(N),N=1,EPTS)
       WRITE(IUTSR) VPTS
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (INXIV(N),INXJV(N),N=1,VPTS)
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (ANG(INXIV(N),INXJV(N)),N=1,VPTS)
       WRITE(IUTSR) FPTS
       IF(FPTS.NE.0)
     . WRITE(IUTSR) (ISFLX(N),JSFLX(N),DIRFLX(N),NFLXE(N),N=1,FPTS)
       CONSTSR=.FALSE.
      ENDIF
!
      IF (SEDTRAN.EQ.'INCLUDE') THEN
!
            DO 680 N = 1, VPTS
              II = INXIV(N)
              JJ = INXJV(N)
!
              IF (IBMSK(II,JJ).EQ.0) THEN
                TSETOT=0.0
                DO 670 LL=1,LAYMAX
                  TSETOT=TSETOT+TSED0(LL,II,JJ)
 670            CONTINUE
!
                THSAVE(N)=THSAVE(N)-TSETOT
!
                IF (SEDTYPE.EQ.'SAND') THSAVE(N)=0.0
!
!  CONVERT FROM g/cm**2 TO cm
!
                THSAVE(N)=THSAVE(N)/CBED(II,JJ)
              ELSE
                THSAVE(N)=THSAVE(N)/((CBED(II,JJ)/2.65)*H1(II,JJ)
     +                                *H2(II,JJ))-BEDTHI
!
!  CONVERT FROM m TO cm
!
                THSAVE(N)=100.*THSAVE(N)
!
                IF (SEDTYPE.EQ.'MUD ') THSAVE(N)=0.0
              ENDIF
!
!  CONVERT FROM g/cm**3 TO mg/l
!
              DO 690 K=1,KBM1
                C1SAVE(N,K)=1000000.*C1SAVE(N,K)
                C2SAVE(N,K)=1000000.*C2SAVE(N,K)
 690          CONTINUE
!
              IF (CHEMTRAN.EQ.'INCLUDE') THEN
!
!  CONVERT FROM ug CHEM/cm**3 TO ug CHEM/l
!
                DO 691 K=1,KBM1
                  P1SAVE(N,K)=1000.*P1SAVE(N,K)
                  P2SAVE(N,K)=1000.*P2SAVE(N,K)
 691            CONTINUE
              ENDIF
!
 680        CONTINUE
      ENDIF
!
      TMIDDLE=TIME-(.5*DTI*DAYI/SKILLI)
      WRITE(IUTSR) TMIDDLE
!
      IF(EPTS.NE.0)
     .WRITE(IUTSR) (ESAVE(N),N=1,EPTS)
!NG04092008      WRITE(IUTSR) (WUSAVE(N),N=1,EPTS)
!NG04092008      WRITE(IUTSR) (WVSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (TXSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (TYSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (PATMSAVE(N),N=1,EPTS)
!NG 04092008 ADDED 2D WINDS AND MET
      WRITE(IUTSR) (WSX2DSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (WSY2DSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (CLDSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (AIRTSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (RHUMSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (SWOBSSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (QPRECSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (QEVAPSAVE(N),N=1,EPTS)
!NG 04142008 INCLUDED WAVES
      WRITE(IUTSR) (WVPDSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (WVHTSAVE(N),N=1,EPTS)
      WRITE(IUTSR) (WVDRSAVE(N),N=1,EPTS)

      IF(TOR.EQ.'BAROTROPIC') THEN
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (UZSAVE(N,1),VZSAVE(N,1),N=1,VPTS)
!
       IF (TRACER.EQ.'INCLUDE') THEN
         WRITE(IUTSR) (C1ZSAVE(N,1),N=1,VPTS)
       ENDIF
!
       IF (SEDTRAN.EQ.'INCLUDE') THEN
         IF (VPTS.NE.0) THEN
           IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
            WRITE (IUTSR) (WS1SAVE(N,1),N=1,VPTS)
            WRITE (IUTSR) (C1SAVE(N,1),N=1,VPTS)
           END IF
           IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
            WRITE (IUTSR) (WS2SAVE(N,1),N=1,VPTS)
            WRITE (IUTSR) (C2SAVE(N,1),N=1,VPTS)
           END IF
           WRITE (IUTSR) (TAUSAVE(N,1),N=1,VPTS)
           WRITE (IUTSR) (THSAVE(N),N=1,VPTS)
         ENDIF
       ENDIF
!
       IF (CHEMTRAN.EQ.'INCLUDE') THEN
         IF (VPTS.NE.0) THEN
           WRITE (IUTSR) (P1SAVE(N,1),P2SAVE(N,1),N=1,VPTS)
           WRITE (IUTSR) ((PBEDSAVE(N,LL),N=1,VPTS),LL=1,NCHEMLAY)
         ENDIF
       ENDIF
!
!PL 04/09/2002 MODIFY TO INCLUDE SEDIMENT CROSS-SECTION FLUX - NEXT 12 LINES
!PL
!PL       IF(FPTS.NE.0)
!PL     . WRITE(IUTSR) (CCFLUX(N,1),N=1,FPTS)
!PL       WRITE(IUTSR)  ESUM,TKE,APE
       IF(FPTS.NE.0) THEN
         WRITE(IUTSR) (CCFLUX(N,1),N=1,FPTS)
         IF(SEDTRAN.EQ.'INCLUDE')THEN
           WRITE(IUTSR) (CS1FLUX(N,1),N=1,FPTS)
           WRITE(IUTSR) (CS2FLUX(N,1),N=1,FPTS)
         END IF
       END IF
       WRITE(IUTSR)  ESUM,TKE,APE
!
      ELSE
       IF(VPTS.NE.0)
     . WRITE(IUTSR) (DZSAVE(N),N=1,VPTS)
       IF(VPTS.NE.0)
     . WRITE(IUTSR) ((UZSAVE(N,K),VZSAVE(N,K),SZSAVE(N,K),
     .                TZSAVE(N,K),N=1,VPTS),K=1,KBM1)
!
       IF (TRACER.EQ.'INCLUDE') THEN
         WRITE(IUTSR) ((C1ZSAVE(N,K),N=1,VPTS),K=1,KBM1)
       ENDIF
!
       IF (SEDTRAN.EQ.'INCLUDE') THEN
        IF (SEDTYPE.EQ.'MUD '.OR.SEDTYPE.EQ.'BOTH') THEN
         WRITE (IUTSR) ((WS1SAVE(N,K),N=1,VPTS),K=1,KBM1)
         WRITE (IUTSR) ((C1SAVE(N,K), N=1,VPTS),K=1,KBM1)
        END IF
        IF (SEDTYPE.EQ.'SAND'.OR.SEDTYPE.EQ.'BOTH') THEN
         WRITE (IUTSR) ((WS2SAVE(N,K),N=1,VPTS),K=1,KBM1)
         WRITE (IUTSR) ((C2SAVE(N,K), N=1,VPTS),K=1,KBM1)
        END IF
         WRITE (IUTSR) ((TAUSAVE(N,K),N=1,VPTS),K=1,KB)
         WRITE (IUTSR) (THSAVE(N),N=1,VPTS)
           WRITE (IUTSR) (TAUSAVE(N,KB),N=1,VPTS)
       ENDIF
!
       IF (CHEMTRAN.EQ.'INCLUDE') THEN
        IF (VPTS.NE.0) THEN
          WRITE (IUTSR) ((P1SAVE(N,K),P2SAVE(N,K),N=1,VPTS),
     +                            K=1,KBM1)
          WRITE (IUTSR) ((PBEDSAVE(N,LL),N=1,VPTS),LL=1,NCHEMLAY)
        ENDIF
       ENDIF
!
!PL 04/09/2002 MODIFY TO INCLUDE SEDIMENT CROSS-SECTION FLUX - NEXT 10 LINES
!PL
!PL    IF(FPTS.NE.0)
!PL  . WRITE(IUTSR) ((CCFLUX(N,K),N=1,FPTS),K=1,KBM1)
       IF(FPTS.NE.0) THEN
         WRITE(IUTSR) ((CCFLUX(N,K),N=1,FPTS),K=1,KBM1)
         IF(SEDTRAN.EQ.'INCLUDE')THEN
           WRITE(IUTSR) ((CS1FLUX(N,K),N=1,FPTS),K=1,KBM1)
           WRITE(IUTSR) ((CS2FLUX(N,K),N=1,FPTS),K=1,KBM1)
         END IF
       END IF
!
       VSTOR=ESUM
       EM=(EM-EMI)*SKILLI
       APEC=(APE-APEI)*SKILLI
       TSUM=TSUM/VOLUME
       SSUM=SSUM/VOLUME
!
       WRITE(IUTSR) VSTOR,EM,APEC,TSUM,SSUM
      ENDIF
 2300 CONTINUE
!
      DO N=1,EPTS
         ESAVE(N)=0.0
!NG04092008         WUSAVE(N)=0.0
!NG04092008         WVSAVE(N)=0.0
         TXSAVE(N)=0.0
         TYSAVE(N)=0.0
         PATMSAVE(N)=0.0
!NG 04092008 ADDED 2D WINDS AND MET
         WSX2DSAVE(N)=0.0
         WSY2DSAVE(N)=0.0
         CLDSAVE(N)=0.0
         AIRTSAVE(N)=0.0
         RHUMSAVE(N)=0.0
         SWOBSSAVE(N)=0.0
         QPRECSAVE(N)=0.0
         QEVAPSAVE(N)=0.0
!NG 04142008 INCLUDED WAVES
         WVPDSAVE(N)=0.0
         WVHTSAVE(N)=0.0
         WVDRSAVE(N)=0.0
      ENDDO
      DO 235 N=1,VPTS
      DZSAVE(N)=0.0
      THSAVE(N)=0.0
 235  CONTINUE
      DO 240 K=1,KB
      DO 240 N=1,VPTS
      UZSAVE(N,K)=0.0
      VZSAVE(N,K)=0.0
      WZSAVE(N,K)=0.0
      SZSAVE(N,K)=0.0
      TZSAVE(N,K)=0.0
      KHZSAVE(N,K)=0.0
      KMZSAVE(N,K)=0.0
      AMZSAVE(N,K)=0.0
      C1ZSAVE(N,K)=0.0
      WS1SAVE(N,K)= 0.0
      WS2SAVE(N,K)= 0.0
      TAUSAVE(N,K)= 0.0
      C1SAVE(N,K) = 0.0
      C2SAVE(N,K) = 0.0
      THSAVE(N)   = 0.0
      P1SAVE(N,K) = 0.0
      P2SAVE(N,K) = 0.0
 240  CONTINUE
!
      DO 241 LL=1,NCHEMLAY
      DO 241 N=1,VPTS
      PBEDSAVE(N,LL)=0.0
 241  CONTINUE
!
      DO 245 K=1,KB
      DO 245 N=1,FPTS
      CCFLUX(N,K)=0.0
      CS1FLUX(N,K)=0.0
      CS2FLUX(N,K)=0.0
 245  CONTINUE
      ESUM  =0.0
      TKE   =0.0
      APE   =0.0
      VSTOR =0.0
      VOLUME=0.0
      EM    =0.0
      APEC  =0.0
      TSUM  =0.0
      SSUM  =0.0
 8450 CONTINUE
!
      RETURN
      END

