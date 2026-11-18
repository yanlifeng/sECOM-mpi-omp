      SUBROUTINE EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2)
      INCLUDE 'comdeck'
      REAL ADVUA(IM,JM),ADVVA(IM,JM),TRNU(IM,JM),TRNV(IM,JM)
!$OMP PARALLEL DEFAULT(SHARED)
      CALL EXTRNL_TEAM(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,.FALSE.)
!$OMP END PARALLEL
      END
      RECURSIVE SUBROUTINE EXTRNL_TEAM(ADVUA,ADVVA,TRNU,TRNV,DTE2,DTI2,ROTATE)
!pl      SUBROUTINE EXTRNL(ADVUA,ADVVA,TRNU,TRNV,DTE2,BFRIC,DTI2)
!     VERSION(11/16/90)
!
      INCLUDE 'comdeck'
      SAVE GEOMX,GEOMY
      REAL*8 QEPC
!
!-----------------------------------------------------------------------
!     CALCULATE FREE SURFACE HEIGHT
!-----------------------------------------------------------------------
!
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM),TRNU(IM,JM),TRNV(IM,JM)
      REAL GEOMX(IM,JM),GEOMY(IM,JM)
      LOGICAL ROTATE
      REAL CORBASE
!
!$OMP DO SCHEDULE(STATIC)
      DO 400 J=2,JMM1
      DO 400 I=2,IM
 400  FLUXUA(I,J)=.25*(D(I,J)+D(I-1,J))*(H2(I,J)+H2(I-1,J))*UA(I,J)
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO 405 J=2,JM
      DO 405 I=2,IMM1
 405  FLUXVA(I,J)=.25*(D(I,J)+D(I,J-1))*(H1(I,J)+H1(I,J-1))*VA(I,J)
!
!wad      DO J=2,JMM1
!wad         DO I=2,IMM1
!wad            IF(D(I,J).LE.WETMIN) THEN
!wad               ELF(I,J)=EL(I,J)-0.5*DTE2*
!wad     .         (FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))/
!wad     &         ART(I,J)
!wad            ELSE
!wad               ELF(I,J)=ELB(I,J)-DTE2*
!wad     .         (FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))/
!wad     &         ART(I,J)
!wad            ENDIF
!wad         ENDDO
!wad      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
      DO 410 J=2,JMM1
      DO 410 I=2,IMM1
 410  ELF(I,J)=ELB(I,J)-DTE2*
     . (FLUXUA(I+1,J)-FLUXUA(I,J)+FLUXVA(I,J+1)-FLUXVA(I,J))/ART(I,J)
!NG12162010 INCLUDE RAIN/EVAP FLUX
     . -DTE2*(QEVAP2D(I,J)-QPREC2D(I,J))*RAMP*FSM(I,J)
!-----------------------------------------------------------------------
!         IMPOSE MASS FLUX BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
!$OMP END DO
!$OMP SINGLE
      DO 406 N=1,NUMDBC
      ID=IDD(N)
      JD=JDD(N)
      ELF(ID,JD)=ELF(ID,JD)+ DTE2*QDIFF(N)*RAMP /ART(ID,JD)
  406 CONTINUE
!
!-----------------------------------------------------------------------
!NG04232008 IMPOSE ELEVATION POTENTIAL INTERIOR CONDITIONS
!-----------------------------------------------------------------------
      IF (OPTEPC(1:8).EQ.'PCLAMPED') THEN
       DO N=1,NUMEPC
        I=IEPC(N)
        J=JEPC(N)
        IF (EPCOR(N).NE.-99.AND.D(I,J).GT.WETMIN) THEN
         CPHL=DUM(I,J)*DTE2/(0.5*(H1(I,J)+H1(I-1,J)))*
     +      SQRT(GRAV*0.5*(D(I,J)+D(I-1,J)))
         CPHR=DUM(I+1,J)*DTE2/(0.5*(H1(I,J)+H1(I+1,J)))*
     +      SQRT(GRAV*0.5*(D(I,J)+D(I+1,J)))
         CPHU=DVM(I,J+1)*DTE2/(0.5*(H2(I,J)+H2(I,J+1)))*
     +      SQRT(GRAV*0.5*(D(I,J)+D(I,J+1)))
         CPHD=DVM(I,J)*DTE2/(0.5*(H2(I,J)+H2(I,J-1)))*
     +      SQRT(GRAV*0.5*(D(I,J)+D(I,J-1)))
         ELFTEMP=(ELB(I,J)+CPHL*(2.*EL(I-1,J)-ELB(I,J))+DTE2*
     +        (EPCOR(N)-ELB(I,J))/TNUDGEPC(N))/(1.+CPHL) +
     +           (ELB(I,J)+CPHR*(2.*EL(I+1,J)-ELB(I,J))+DTE2*
     +        (EPCOR(N)-ELB(I,J))/TNUDGEPC(N))/(1.+CPHR) +
     +           (ELB(I,J)+CPHU*(2.*EL(I,J+1)-ELB(I,J))+DTE2*
     +        (EPCOR(N)-ELB(I,J))/TNUDGEPC(N))/(1.+CPHU) +
     +           (ELB(I,J)+CPHD*(2.*EL(I,J-1)-ELB(I,J))+DTE2*
     +        (EPCOR(N)-ELB(I,J))/TNUDGEPC(N))/(1.+CPHD)
         ELFTEMP=0.25*ELFTEMP
         QEPC=RAMP*ART(I,J)*(ELFTEMP-ELF(I,J))/DTE2 ! Q is positive incoming (if ELFTEMP>ELF)
         QEPCNUDF(N)=QEPCNUDF(N)+QEPC*ISPI ! AVERAGE OVER THE EXTERNAL MODE
         ELF(I,J)=ELF(I,J) + DTE2*QEPC/ART(I,J)
        ELSE
         QEPCNUDF(N)=0.0 ! ALLOW NO EPC FLUX
        ENDIF
       ENDDO
      ELSE ! FOR NUDGING
       DO N=1,NUMEPC
        I=IEPC(N)
        J=JEPC(N)
        IF (EPCOR(N).NE.-99.AND.D(I,J).GT.WETMIN) THEN
!NG     Note: I tested the one below as well, (epcor applied on a dry
!NG     element, if positive). It also worked, but it was not as smooth as the
!NG     above (only allowing lateral, not EPCOR-related wetting).
!NG        IF ( (EPCOR(N).NE.-99).AND.
!NG     .       (D(I,J).GT.WETMIN.OR.EPCOR(N).GT.ELB(I,J)) ) THEN
         QEPC=RAMP*ART(I,J)*(EPCOR(N)-ELB(I,J))/TNUDGEPC(N) ! Q is positive incoming (if EPCOR(N)>ELB)
         QEPCNUDF(N)=QEPCNUDF(N)+QEPC*ISPI ! AVERAGE OVER THE EXTERNAL MODE
         ELF(I,J)=ELF(I,J) + DTE2*QEPC/ART(I,J)
        ELSE
         QEPCNUDF(N)=0.0 ! ALLOW NO EPC FLUX
        ENDIF
       ENDDO
      ENDIF
! 
      CALL BCOND(12,DTI2,0)
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
      DO J=1,JM
        DO I=1,IM
          ELF(I,J)=ELF(I,J)*FSM(I,J)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP SINGLE
!
!     Introduced Inverted Reid and Bodine (IRANDB) and Optimized Clamped 
!     and MIXED TYPE BOUNDARY CONDITIONS In Barotropic or external Mode
!
!  MODIFIED BY Quamrul QA 01/03/99
!
!qa        IF (BCTYPE.EQ.'RANDB  '.OR.BCTYPE.EQ.'SHULMAN') THEN
!qa            CALL SHULMAN
!qa        ELSE IF (BCTYPE.EQ.'PCLAMP ') THEN
!qa            CALL PCLAMP
!qa        ENDIF

!     Modified by Quamrul QA adopted from JKL 8/23/00

      IF (BCTYPE.EQ.'OCLAMP '.OR.BCTYPE.EQ.'IRANDB '.OR.
     +    BCTYPE.EQ.'MIXED  ') CALL OCLAMP(CLAM)
      IF (BCTYPE.EQ.'PCLAMP ') CALL PCLAMP
!
!$OMP END SINGLE
! ADVAVE uses the existing EXTRNL team; all workers must call it.
      CALL ADVAVE(ADVUA,ADVVA)
!
! Nonbarotropic callers keep the original complete Coriolis field.
      IF(TOR.NE.'BAROTROPIC') THEN
!$OMP DO SCHEDULE(STATIC)
       DO 95 J=1,JM
       DO 95 I=1,IM
 95    CURV42D(I,J)=.25*COR(I,J)
!$OMP END DO
       GO TO 5000
      ENDIF
!
!     The reference's permanent grid metrics and FSM cannot change
!     within this external subcycle. Dynamic wet masks are NOT cached.
!     Refresh on every outer step (also works on hot starts).
      IF(IEXT.EQ.1) THEN
!$OMP DO SCHEDULE(STATIC)
       DO J=2,JMM1
        DO I=2,IMM1
         GEOMX(I,J)=(H2(I+1,J)*FSM(I+1,J)+H2(I,J)*FSM(I,J))/
     .          (FSM(I+1,J)+FSM(I,J)+1.E-30)-
     .      (H2(I,J)*FSM(I,J)+H2(I-1,J)*FSM(I-1,J))/
     .          (FSM(I,J)+FSM(I-1,J)+1.E-30)
         GEOMY(I,J)=(H1(I,J+1)*FSM(I,J+1)+H1(I,J)*FSM(I,J))/
     .          (FSM(I,J+1)+FSM(I,J)+1.E-30)-
     .      (H1(I,J)*FSM(I,J)+H1(I,J-1)*FSM(I,J-1))/
     .          (FSM(I,J)+FSM(I,J-1)+1.E-30)
        ENDDO
       ENDDO
!$OMP END DO
      ENDIF
! Both edges and interior are owned by one row workshare. Each cell
! gets one store; the separate REAL base preserves the original
! rounding before the unchanged curvature correction. No halo read of
! CURV42D occurs until after the subsequent stress barrier.
!$OMP DO SCHEDULE(STATIC) PRIVATE(CORBASE)
      DO J=1,JM
       IF(J.EQ.1.OR.J.EQ.JM) THEN
        DO I=1,IM
         CURV42D(I,J)=.25*COR(I,J)
        ENDDO
       ELSE
        CURV42D(1,J)=.25*COR(1,J)
        CURV42D(IM,J)=.25*COR(IM,J)
        DO I=2,IMM1
         CORBASE=.25*COR(I,J)
         CURV42D(I,J)= CORBASE
     .  +.125*((VA(I,J+1)+VA(I,J))*GEOMX(I,J)
     .    -(UA(I+1,J)+UA(I,J))*GEOMY(I,J))
     .                /ART(I,J)
        ENDDO
       ENDIF
      ENDDO
!
! Neither stress direction reads CURV42D. The barrier after TYICE
! still completes all three independent fields before boundary edits.
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO 100 J=2,JMM1
      DO 100 I=2,IMM1
!NG04232008      TXBOT(I,J)=-BFRIC
      TXBOT(I,J)=-0.5*(CBC(I,J)+CBC(I-1,J))
     .     *SQRT(UAB(I,J)**2+(.25*(VAB(I,J)
     .     +VAB(I,J+1)+VAB(I-1,J)+VAB(I-1,J+1)))**2)*UAB(I,J)
      TXICE(I,J)=0.5*(AREAICE(I,J)/100.*VARWIF(I,J)      ! NG FOR ICE
     .               +AREAICE(I-1,J)/100.*VARWIF(I-1,J)) ! NG FOR ICE
     .     *SQRT(UAB(I,J)**2+(.25*(VAB(I,J)
     .     +VAB(I,J+1)+VAB(I-1,J)+VAB(I-1,J+1)))**2)*UAB(I,J)
 100  CONTINUE
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO 102 J=2,JMM1
      DO 102 I=2,IMM1
!NG04232008      TYBOT(I,J)=-BFRIC
      TYBOT(I,J)=-0.5*(CBC(I,J)+CBC(I,J-1))
     .     *SQRT((.25*(UAB(I,J)+UAB(I+1,J)
     .     +UAB(I,J-1)+UAB(I+1,J-1)))**2+VAB(I,J)**2)*VAB(I,J)
      TYICE(I,J)=-0.5*(AREAICE(I,J)/100.*VARWIF(I,J)      ! NG FOR ICE
     .                +AREAICE(I,J-1)/100.*VARWIF(I,J-1)) ! NG FOR ICE
     .     *SQRT((.25*(UAB(I,J)+UAB(I+1,J)
     .     +UAB(I,J-1)+UAB(I+1,J-1)))**2+VAB(I,J)**2)*VAB(I,J)
 102  CONTINUE

!NG INCLUDE SIMILARLY TO 3D TREATMENT 
!NG WHETHER THIS IS A GOOD IDEA NEEDS FURTHER INVESTIGATION!
!$OMP END DO
!$OMP SINGLE
      DO 310 N=1,NUMEBC
      IE=IETA(N)
      JE=JETA(N)
      IC=ICON(N)
      JC=JCON(N)
      IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        CURV42D(IE,JE)=CURV42D(IE-1,JE)
        TYBOT(IE,JE)=TYBOT(IE-1,JE) ! NG same as internal bcond, cosmetic
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
        CURV42D(IE,JE)=CURV42D(IE+1,JE)
        TYBOT(IE,JE)=TYBOT(IE+1,JE) ! NG same as internal bcond, cosmetic
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
        CURV42D(IE,JE)=CURV42D(IE,JE-1)
        TXBOT(IE,JE)=TXBOT(IE,JE-1) ! NG same as internal bcond, cosmetic
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
        CURV42D(IE,JE)=CURV42D(IE,JE+1)
        TXBOT(IE,JE)=TXBOT(IE,JE+1) ! NG same as internal bcond, cosmetic
      ENDIF
 310  CONTINUE
!$OMP END SINGLE

 5000 CONTINUE

!
!-----------------------------------------------------------------------
!     CALCULATE U COMPONENT OF VELOCITY
!-----------------------------------------------------------------------
!

!NG2012 CORNERS FIX - NO NEED TO HARDWIRE ANYMORE
!$OMP SINGLE
      DO N=1,NUMECR
       IF(NCHARD(N).EQ.1) THEN ! Bottom Left Corner
        Va(ICHARD(N),JCHARD(N)+1) =
     +   .5*(Va(ICHARD(N)+1,JCHARD(N)+1)+Va(ICHARD(N),JCHARD(N)+2))
        Ua(ICHARD(N)+1,JCHARD(N)) =
     +   .5*(Ua(ICHARD(N)+1,JCHARD(N)+1)+Ua(ICHARD(N)+2,JCHARD(N)))
       ELSEIF(NCHARD(N).EQ.2) THEN ! Bottom Right Corner
        Va(ICHARD(N),JCHARD(N)+1) =
     +   .5*(Va(ICHARD(N)-1,JCHARD(N)+1)+Va(ICHARD(N),JCHARD(N)+2))
        Ua(ICHARD(N),JCHARD(N)) =
     +   .5*(Ua(ICHARD(N)-1,JCHARD(N))+Ua(ICHARD(N),JCHARD(N)+1))
       ELSEIF(NCHARD(N).EQ.3) THEN ! Top Right Corner
        Va(ICHARD(N),JCHARD(N)) =
     +   .5*(Va(ICHARD(N)-1,JCHARD(N))+Va(ICHARD(N),JCHARD(N)-1))
        Ua(ICHARD(N),JCHARD(N)) =
     +   .5*(Ua(ICHARD(N)-1,JCHARD(N))+Ua(ICHARD(N),JCHARD(N)-1))
       ELSE ! Top Left Corner
        Va(ICHARD(N),JCHARD(N)) =
     +   .5*(Va(ICHARD(N)+1,JCHARD(N))+Va(ICHARD(N),JCHARD(N)-1))
        Ua(ICHARD(N)+1,JCHARD(N)) =
     +   .5*(Ua(ICHARD(N)+2,JCHARD(N))+Ua(ICHARD(N)+1,JCHARD(N)-1))
       ENDIF
      ENDDO

!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC)
      DO 425 J=2,JMM1
      DO 425 I=3,IMM1
      UAF(I,J)=ADVUA(I,J)
     .    -ARU(I,J)*( CURV42D(I,J)*D(I,J)*(VA(I,J+1)+VA(I,J))
     .              +CURV42D(I-1,J)*D(I-1,J)*(VA(I-1,J+1)+VA(I-1,J)) )
     .         +.25*GRAV*(H2(I,J)+H2(I-1,J))*(D(I,J)+D(I-1,J))
     .             *( (1.-2.*THETA)*(EL(I,J)-EL(I-1,J))
     .            +THETA*(ELB(I,J)-ELB(I-1,J)+ELF(I,J)-ELF(I-1,J)) )
     .         +.25*GRAV*(H2(I,J)+H2(I-1,J))*(D(I,J)+D(I-1,J))
     .             *(DATUM(I,J)-DATUM(I-1,J))*RAMP
     .         +RAMP*TRNU(I,J)
!NG     .      -ARU(I,J)*(-.5*(TXSURF(I,J)+TXSURF(I-1,J))+TXBOT(I,J)  )
     .      -ARU(I,J)*
     .      ( -.5*((100.-AREAICE(I,J))/100.*TXSURF(I,J)+      ! NG FOR ICE
     .             (100.-AREAICE(I-1,J))/100.*TXSURF(I-1,J) ) ! NG FOR ICE
     .                -TXICE(I,J)                             ! NG FOR ICE
     .                +TXBOT(I,J)  )
!qa         Add Atmospheric Pressure Terms by Quamrul QA 4/20/00
     *            +ARU(I,J)*(D(I,J)+D(I-1,J))*(PATM(I,J)-PATM(I-1,J))
!NG     *             /RHO0/(H1(I,J)+H1(I-1,J))
     *            *RAMP*PATMOPT/RHO0/(H1(I,J)+H1(I-1,J))

 425  UAF(I,J)=
     .         ((H(I,J)+ELB(I,J)+H(I-1,J)+ELB(I-1,J))*ARU(I,J)*UAB(I,J)
     .                -4.*DTE*UAF(I,J))
     .        /((H(I,J)+ELF(I,J)+H(I-1,J)+ELF(I-1,J))*ARU(I,J))
!
!-----------------------------------------------------------------------
!     CALCULATE V COMPONENT OF VELOCITY
!-----------------------------------------------------------------------
!
! VAF reads UA/VA, not UAF. Both updates finish before BCOND(13).
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO 435 J=3,JMM1
      DO 435 I=2,IMM1
      VAF(I,J)=ADVVA(I,J)  
     .    +ARV(I,J)*( CURV42D(I,J)*D(I,J)*(UA(I+1,J)+UA(I,J))
     .                +CURV42D(I,J-1)*D(I,J-1)*(UA(I+1,J-1)+UA(I,J-1)) )
     .           +.25*GRAV*(H1(I,J)+H1(I,J-1))*(D(I,J)+D(I,J-1))
     .                *( (1.-2.*THETA)*(EL(I,J)-EL(I,J-1))
     .                +THETA*(ELB(I,J)-ELB(I,J-1)+ELF(I,J)-ELF(I,J-1)) )
     .           +.25*GRAV*(H1(I,J)+H1(I,J-1))*(D(I,J)+D(I,J-1))
     .                *(DATUM(I,J)-DATUM(I,J-1))*RAMP
     .           +RAMP*TRNV(I,J)
!NG     .    + ARV(I,J)*( .5*(TYSURF(I,J)+TYSURF(I,J-1))-TYBOT(I,J)   )
     .    + ARV(I,J)*
     .      (  .5*((100.-AREAICE(I,J))/100.*TYSURF(I,J)+      ! NG FOR ICE
     .             (100.-AREAICE(I,J-1))/100.*TYSURF(I,J-1) ) ! NG FOR ICE
     .               +TYICE(I,J)                              ! NG FOR ICE
     .               -TYBOT(I,J)   )
!qa         Add Atmospheric Pressure Terms by Quamrul QA 4/20/00
     *            +ARV(I,J)*(D(I,J)+D(I,J-1))*(PATM(I,J)-PATM(I,J-1))
!NG     *             /RHO0/(H2(I,J)+H2(I,J-1))
     *            *RAMP*PATMOPT/RHO0/(H2(I,J)+H2(I,J-1))
 435  VAF(I,J)=
     .        ((H(I,J)+ELB(I,J)+H(I,J-1)+ELB(I,J-1))*VAB(I,J)*ARV(I,J)
     .              -4.*DTE*VAF(I,J))
     .       /((H(I,J)+ELF(I,J)+H(I,J-1)+ELF(I,J-1))*ARV(I,J))
!
!JKL
!$OMP END DO
!$OMP SINGLE
      IF(BCTYPE.EQ.'RANDB  ') CALL RANDB
      IF(BCTYPE.EQ.'MIXED  ') CALL MIXED(1)
!JKL

      CALL BCOND(13,DTI2,0)
!$OMP END SINGLE
! On the no-transport barotropic path filtering reads the masked
! results only at the very same cell. Rotate before those values
! leave cache; retain the old interface and timing for other paths.
      IF(ROTATE) THEN
!$OMP DO SCHEDULE(STATIC)
        DO J=1,JM
          DO I=1,IM
            TXBOT(I,J)=TXBOT(I,J)*DUM(I,J)
            TYBOT(I,J)=TYBOT(I,J)*DVM(I,J)
            UAF(I,J)=UAF(I,J)*DUM(I,J)
            VAF(I,J)=VAF(I,J)*DVM(I,J)
            UAB(I,J)=UA(I,J)+.5*NU*(UAB(I,J)-2.*UA(I,J)+UAF(I,J))
            VAB(I,J)=VA(I,J)+.5*NU*(VAB(I,J)-2.*VA(I,J)+VAF(I,J))
            ELB(I,J)=EL(I,J)+.5*NU*(ELB(I,J)-2.*EL(I,J)+ELF(I,J))
            EL(I,J)=ELF(I,J)
            D(I,J)=H(I,J)+EL(I,J)
            UA(I,J)=UAF(I,J)
            VA(I,J)=VAF(I,J)
          ENDDO
        ENDDO
!$OMP END DO
      ELSE
!$OMP DO SCHEDULE(STATIC)
      DO J=1,JM
        DO I=1,IM
          TXBOT(I,J)=TXBOT(I,J)*DUM(I,J)
          TYBOT(I,J)=TYBOT(I,J)*DVM(I,J)
          UAF(I,J)=UAF(I,J)*DUM(I,J)
          VAF(I,J)=VAF(I,J)*DVM(I,J)
        ENDDO
      ENDDO
!$OMP END DO
      ENDIF
!
!hli      print *,'returning from EXTRNL'
      RETURN
      END
