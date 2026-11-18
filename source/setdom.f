      SUBROUTINE SETDOM
!
!     VERSION(06/27/91)
!
      INCLUDE 'comdeck'
!
!  FOR WET GRID INPUT FROM POM
!        
      CHARACTER*10 RESTAR
      COMMON /COAST1/ICNT,INDX(MAXWET),JNDX(MAXWET),RESTAR
!
      SAVE
!
!-----------------------------------------------------------------------
!     INITIAL CONDITION PREPARATION PROGRAM FOR A GENERAL CIRCULATION
!        MODEL IN ORTHOGONAL CURVILINEAR COORDINATES
!-----------------------------------------------------------------------
!
!     DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),ADVUU(IM,JM),ADVVV(IM,JM)
!     DIMENSION DRHOX(IM,JM,KB),DRHOY(IM,JM,KB),TRNU(IM,JM),TRNV(IM,JM)
      DIMENSION IVAR(IM,JM),PRT(IM,KB)
      EQUIVALENCE (IVAR,TPS),(PRT,A)
      DIMENSION COM(80)
      DIMENSION ISTDAM(NUMDAM),JSTDAM(NUMDAM),NTDAM(NUMDAM)
      CHARACTER DIRTDAM(NUMDAM)*4
! 
!-------- ALL UNITS IN M.K.S. SYSTEM -----------------------------------
      ISTART=1
      IEND=1
      INTX=0 
!
!  FOR EXTERNAL HYDRODYNAMICS RUN
!
!  ASSUME POM OUTPUT FOR WET GRID
!
      IF (HYDTYPE.EQ.'EXTERNAL') THEN
        IF(RESTAR.EQ.'HOT START '.AND.IWET.EQ.1) GOTO 153
!qa
        OPEN (IUTRN,FILE='gcm_geom',FORM='UNFORMATTED')
          READ (IUTRN)DZ,DZZ
          READ (IUTRN)H,H1,H2,TPS
          READ (IUTRN) ANG,NU
!         CONVERTED TO DEGREES FROM RADIANS BY QUAMRUL QA 9/15/98
          DO J=1,JM
           DO I=1,IM
             ANG(I,J)=ANG(I,J)*360.0/(2.*3.141593)
           ENDDO
          ENDDO
        CLOSE(IUTRN)
!
        OPEN (IUTRN,FORM='unformatted',FILE='wet_grid')
          READ (IUTRN)ICNT
        DO I=1,ICNT
          READ(IUTRN)INDX(I),JNDX(I)
        END DO
        CLOSE (IUTRN)
!
        DO 151 K=1,KBM1
          DZR(K)=1./DZ(K)
 151    CONTINUE
!
!  RESET LAND ELEMENTS TO ZERO DEPTH (SET TO 1.000 m IN POM)
!
        DO 152 J=1,JM
          DO 152 I=1,IM
            IF (H(I,J).LT.1.002) H(I,J)=0.0
 152    CONTINUE
        GOTO 999
!
!qa
!qa     Reads wet_grid when "HOT START" and IWET=1 i.e. IWET on
!qa
 153    CONTINUE
        OPEN (IUTRN,FORM='unformatted',FILE='wet_grid')
          READ (IUTRN)ICNT
        DO I=1,ICNT
          READ(IUTRN)INDX(I),JNDX(I)
        END DO
        CLOSE (IUTRN)
        RETURN
      ENDIF
!
      WRITE(IUPRT,51)            
 51   FORMAT(/'...... MODEL STARTING UP FROM INITAL CONDITIONS .......')
      
      CALL NETCDF_READ_MODEL_GRID('SETDOM') ! A.Jordi
!
      call RDCORLOC !csvv080327
!
!-------- DEFINE MASK FOR FREE SURFACE HEIGHT = FSM --------------------
!-------- DEFINE MASK FOR (U,V) VELOCITY = (DUM,DVM) -------------------
!-------- IF DEPTH <= 0.0, DEPTH AND MASK SET TO 0 ---------------------
 999  DO 400 J=1,JM
      DO 400 I=1,IM
      FSM(I,J)=1.0
      DUM(I,J)=1.0
      DVM(I,J)=1.0
      IF(H1(I,J).LE.0.0) THEN ! BASE IT ON H1, NOT H, W&D NG04212010
!      IF(H(I,J).LE.0.0) THEN ! BASE IT ON H1, NOT H, W&D NG04212010
!PL 8/16/99      H(I,J)=0.0
!pl      H(I,J)=0.0+1.0E-3
       H(I,J)=0.0+1.0E-10
       FSM(I,J)=0.0
       DUM(I,J)=0.0
       DVM(I,J)=0.0
!      ELSEIF(H(I,J).EQ.-999.99) THEN ! FLAG FOR OPEN BC CORNERS
!       H(I,J)=0.0+1.0E-10
!       FSM(I,J)=0.0
!       DUM(I,J)=0.0
!       DVM(I,J)=0.0
      ELSEIF(H(I,J).EQ.0.0) THEN
       H(I,J)=0.0+1.0E-10
!      ELSEIF(H(I,J).LT.-20.0) THEN
!       H(I,J)=0.0+1.0E-10-1.*float(I*J)/float(IM*JM)
      ENDIF
 400  CONTINUE 
!PL 8/16/99 400  H(I,J)=H(I,J)+1.0E-03
!
!      DO 450 I=1,IM
!       DUM(I,3)=0.0 !!!!!!!!!!!!!!!!!!!!!!!!!!! BIG HARDWIRE HERE
!  450 CONTINUE
      DO 420 J=1,JMM1
      DO 420 I=1,IM
      IF(FSM(I,J).EQ.0.0.AND.FSM(I,J+1).NE.0.0) DVM(I,J+1)=0.0
  420 CONTINUE
      DO 440 J=1,JM
      DO 440 I=1,IMM1
      IF(FSM(I,J).EQ.0.0.AND.FSM(I+1,J).NE.0.0) DUM(I+1,J)=0.0
  440 CONTINUE
! FOR THIN DAM  
      DO 180 N=1,NUMTDAM
       IS=ISTDAM(N)
       JS=JSTDAM(N)
       IF(DIRTDAM(N).EQ.'IDIR') THEN
        DO 181 I=IS,IS+NTDAM(N)-1
         DVM(I,JS)=0.0
  181   CONTINUE
       ELSE    ! JDIR
        DO 182 J=JS,JS+NTDAM(N)-1
         DUM(IS,J)=0.0
  182   CONTINUE
       END IF
  180  CONTINUE

!
      DO 160 J=1,JM
      DO 160 I=1,IM
      DUMWAD(I,J)=DUM(I,J) ! NG05192010
      DVMWAD(I,J)=DVM(I,J) ! NG05192010
      H1(I,J)=H1(I,J)+1.E-10
      H2(I,J)=H2(I,J)+1.E-10
      ANG(I,J)=ANG(I,J)*2.*3.141593/360.
      COR(I,J)=2.*7.292E-5*SIN(COR(I,J)*2.*3.141593/360.)*FSM(I,J)
!NG      COR(I,J)=0.0 ! HARDWIRE
  160 CONTINUE
!
      !FLOODIT=AMAX1(AMIN1(WETEPS,WETMIN/2.),1.E-3) ! W&D NG04212010 MAKE SURE D>0
!NG      FLOODIT=AMAX1(WETMIN-1.E-5,1.E-3) ! W&D NG04212010 MAKE SURE D>0
      FLOODIT=AMAX1(WETMIN-1.E-3,1.E-3) ! W&D NG04212010 MAKE SURE D>0
      DO 260 J=1,JM
      DO 260 I=1,IM
      IF(FSM(I,J).GT.0.0) THEN ! W&D NG05182010 IF ORIGINALLY ABOVE GEOID, OR JUST BELOW, MAKE SURE D>0=FLOODIT
       D(I,J)=AMAX1(FLOODIT,H(I,J))
       DT(I,J)=D(I,J)
       ELB(I,J)=D(I,J)-H(I,J)
       EL(I,J)=D(I,J)-H(I,J)
       ETB(I,J)=D(I,J)-H(I,J)
       ET(I,J)=D(I,J)-H(I,J)
       EGB(I,J)=D(I,J)-H(I,J)
!NG       ELF(I,J)=D(I,J)-H(I,J)
!NG       ETF(I,J)=D(I,J)-H(I,J)
       IF (ELB(I,J).GT.0.0) WRITE(IUPRT,170) I,J,D(I,J)
      ELSE ! NG, JUST FOR TIME STEP CALCULATION
       D(I,J)=H(I,J)
       DT(I,J)=D(I,J)
      ENDIF
  260 CONTINUE
  170 FORMAT(' IJ : ',2I5,' BELOW GEOID. INITIAL D SET TO ',F7.3)
!
      DO 340 J=2,JMM1
      DO 340 I=2,IMM1
      ART(I,J)=H1(I,J)*H2(I,J)
      ARU(I,J)=.25E0*(H1(I,J)+H1(I-1,J))*(H2(I,J)+H2(I-1,J))
      ARV(I,J)=.25E0*(H1(I,J)+H1(I,J-1))*(H2(I,J)+H2(I,J-1))
  340 CONTINUE
      DO 310 J=1,JM
      ARU(IM,J)=ARU(IMM1,J)
  310 ARU(1,J)=ARU(2,J)
      DO 320 I=1,IM
      ARV(I,JM)=ARV(I,JMM1)
  320 ARV(I,1)=ARV(I,2)
!
      DO 410 J=1,JM
      DO 410 I=1,IM
 410  TPS(I,J)=SQRT(H1(I,J)**2+H2(I,J)**2)
! 
      CALL MAXMIN(TPS,FSM,IM,JM,FMAX,FMIN)
! 
      IF(TOR.NE.'BAROTROPIC') THEN
      DO 350 K=1,KBM1
      DO 350 J=1,JM
      DO 350 I=1,IM
      AAM(I,J,K)=HORCON*VARHF(I,J)*TPS(I,J)/FMIN*FSM(I,J)
  350 CONTINUE
      DO 375 J=1,JM
      DO 375 I=1,IM
 375  AAM2D(I,J)=0.0
      DO 380 K=1,KBM1
      DO 380 J=1,JM
      DO 380 I=1,IM
      AAM2D(I,J)=AAM2D(I,J)+AAM(I,J,K)*DZ(K)
 380  CONTINUE
!
      ELSE
!
      DO 355 J=1,JM
      DO 355 I=1,IM
      AAM2D(I,J)=HORCON*VARHF(I,J)*TPS(I,J)/FMIN*FSM(I,J)
  355 CONTINUE
! 
      ENDIF
!
      DO 360 J=1,JM
      DO 360 I=1,IM
      KM(I,J,1)=0.0
      KM(I,J,KB)=0.0
      KH(I,J,1)=0.0
      KH(I,J,KB)=0.0
      KQ(I,J,1)=0.0
      KQ(I,J,KB)=0.0
      L (I,J,1)=0.0
      L (I,J,KB)=0.0
      Q2(I,J,1)=0.0
      Q2(I,J,KB)=0.0
      Q2B(I,J,1)=0.0
      Q2B(I,J,KB)=0.0
      Q2L(I,J,1)=0.0
      Q2L(I,J,KB)=0.0
      Q2LB(I,J,1)=0.0
      Q2LB(I,J,KB)=0.0
  360 CONTINUE
!
      DO 370 K=2,KBM1
      DO 370 J=1,JM
      DO 370 I=1,IM
      L(I,J,K)   =1.*FSM(I,J)
      Q2B(I,J,K) =1.E-12*FSM(I,J)
      Q2(I,J,K)  =Q2B(I,J,K)
      Q2LB(I,J,K)=Q2B(I,J,K)*L(I,J,K)
      Q2L(I,J,K) =Q2LB(I,J,K)
      KM(I,J,K)  =UMOL*VARUF(I,J)*FSM(I,J)
      KQ(I,J,K)  =UMOL*VARUF(I,J)*FSM(I,J)
      KH(I,J,K)  =UMOL*VARUF(I,J)/VPRNU*FSM(I,J)
  370 CONTINUE
!
      IF(VERTMIX.EQ.'CONSTANT  ') UMOL=0.0
!
!*********************************************************************
!
!  FOR PARTICLE TRACKING
!
      IF (PARTICLE.EQ.'INCLUDE') THEN

!
!  OPEN CORNER LOCATIONS FILE
!
!  CORNER LOCATION CONVENTION:  XCOR(I,J) = x(i-1/2, j-1/2) 
!     (LOWER LEFT-HAND CORNER)  YCOR(I,J) = y(i-1/2, j-1/2)
!
        OPEN (UNIT=33,FILE='corner_loc',FORM='FORMATTED')
!
        DO 1300 N=1,1000000
          READ (33,*,END=1310)I,J,XCOR(I,J),YCOR(I,J)
 1300   CONTINUE
!
!  CALC. H1 AND H2 AT ELEMENT INTERFACES
!
 1310   DO 1320 J=1,JM
          DO 1320 I=1,IM
            H1P(I,J)=SQRT((XCOR(I+1,J)-XCOR(I,J))**2.+
     +                        (YCOR(I+1,J)-YCOR(I,J))**2.)
!
            H2P(I,J)=SQRT((XCOR(I,J+1)-XCOR(I,J))**2.+
     +                        (YCOR(I,J+1)-YCOR(I,J))**2.)
 1320   CONTINUE
      ENDIF
!
!*********************************************************************
!
!
      RETURN
      END 
