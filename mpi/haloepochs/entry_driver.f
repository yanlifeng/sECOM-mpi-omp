! Isolated pair of consecutive outer steps, including real complete EXTRNL,
! ADC/V flux stages, WAD and CHECK. Raw owner-only captures; no solver run.
      PROGRAM EPOCH_DRIVER
      INCLUDE 'comdeck'
      REAL AX(IM,JM),AY(IM,JM),RX(IM,JM),RY(IM,JM),RV(JM),RE(JM)
      CHARACTER*8 ARG
      INTEGER LAYOUT,ISEQ,IXOUT,IIT,JIT
      LOGICAL READY
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) LAYOUT
      CALL SMPI_INIT(IM,JM)
      CALL GET_COMMAND_ARGUMENT(2,ARG)
      READ(ARG,*) NRGRID
      TOR='BAROTROPIC'
      HYDTYPE='INTERNAL'
      BCTYPE='CLAMPED'
      OPTEBC='               LINEAR'
      OPTEBC(15:20)='LINEAR'
      OPTEPC='PCLAMPED'
      SMPI_OUTER=LAYOUT.EQ.1
      SMPI_LOCALBC=LAYOUT.NE.2
      NUMQBC=0
      NUMEPC=0
      NUMDBC=0
      NUMEBC=1
      IETA(1)=3
      JETA(1)=2
      ICON(1)=4
      JCON(1)=2
      EBDRY(1)=.031
      NUMECR=1
! Top corner at a slab end; its two required rows are rank-local.
      JE=JM/NRGRID+MIN(1,MOD(JM,NRGRID))
      IF(JE.LT.2) NUMECR=0
      ICHARD(1)=7
      JCHARD(1)=MAX(2,JE)
      NCHARD(1)=3
      WETMIN=.3
      WETEPS=.017
      DTE=.015
      THETA=.151
      NU=.11
      RAMP=.83
      RHO0=1000.
      PATMOPT=.71
      GRAV=9.81
      DO J=1,JM
       DO I=1,IM
       X=.013*MOD(I,17)+.027*J
       H(I,J)=1.2+X
       IF(MOD(I+2*J,7).EQ.1) H(I,J)=.30+.009*MOD(I+J,3)
       H1(I,J)=1.+X*.05
       H2(I,J)=1.+X*.08
       ART(I,J)=H1(I,J)*H2(I,J)
       ARU(I,J)=ART(I,J)*1.1
       ARV(I,J)=ART(I,J)*1.2
       FSM(I,J)=1.
       IF(MOD(I*2+J,19).EQ.0) FSM(I,J)=0.
       DUM(I,J)=1.
       DVM(I,J)=1.
       IF(MOD(I+3*J,11).EQ.4) DUM(I,J)=0.
       IF(MOD(2*I+J,9).EQ.3) DVM(I,J)=0.
       COR(I,J)=.00017*X
       CBC(I,J)=.016+X*.001
       AAM2D(I,J)=.039+X*.01
       DATUM(I,J)=.00031*X
       PATM(I,J)=1000.+17.*X
       VARWIF(I,J)=1.+X*.1
       AREAICE(I,J)=.1*X
       QEVAP2D(I,J)=X*.00001
       QPREC2D(I,J)=X*.00002
       TXSURF(I,J)=.03*X
       TYSURF(I,J)=.021*X
       ETB(I,J)=.001*J-.002*I
       IF(MOD(I+J,7).EQ.0) ETB(I,J)=-0.
       ET(I,J)=ETB(I,J)*.7
       ETF(I,J)=ET(I,J)*1.1
       EL(I,J)=.019*MOD(I+J,3)-.025
       ELB(I,J)=EL(I,J)*1.2
       ELF(I,J)=EL(I,J)*.97
       D(I,J)=H(I,J)+EL(I,J)
       UA(I,J)=.041*MOD(I,17)-.037*J
       VA(I,J)=-.027*MOD(I,17)+.043*J
       UAF(I,J)=UA(I,J)*.89
       VAF(I,J)=VA(I,J)*.91
       UAB(I,J)=UA(I,J)*.73
       VAB(I,J)=VA(I,J)*(-.61)
       RX(I,J)=.023*X
       RY(I,J)=-.017*X
       FLUXUA(I,J)=X*.01
       FLUXVA(I,J)=-X*.02
       AX(I,J)=0.
       AY(I,J)=0.
       ENDDO
      ENDDO
      FSM(4,2)=0.
      OPEN(22,FILE='owner.raw',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='NEW')
      VMAX=0.
      IEXT=0
      EMIN=124.
      IIT=-44
      JIT=77
      DO IXOUT=1,2
       ADVECT='NON-LINEAR'
       IF(IXOUT.EQ.2) ADVECT='LINEAR'
       DUMWAD=DUM
       DVMWAD=DVM
       INTX=IXOUT
       DO ISEQ=1,4
! Same row-wise freshness as the production protocol. Poison both sides
! of all transmitted fields, and all unused ghost sides. Preserve only
! both D ghosts, south EL and north VA at substep 2+ from actual WAD.
        DO J=1,JM
         IF(J.GE.SMPI_JLO.AND.J.LE.SMPI_JHI) CYCLE
         X=20.+.3*J+IEXT*.01
         IF(ISEQ.EQ.1.OR.J.LT.SMPI_JLO-1.OR.J.GT.SMPI_JHI+1) D(:,J)=X
         IF(ISEQ.EQ.1.OR.J.NE.SMPI_JLO-1) EL(:,J)=-X
         IF(ISEQ.EQ.1.OR.J.NE.SMPI_JHI+1) VA(:,J)=X*1.13
         UA(:,J)=-3.*X
         UAB(:,J)=2.*X
         VAB(:,J)=-2.25*X
         ELB(:,J)=.37*X
         DUMWAD(:,J)=.63*X
         DVMWAD(:,J)=.73*X
        ENDDO
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP MASTER
        IEXT=ISEQ
!$OMP END MASTER
        CALL EXTRNL_TEAM(AX,AY,RX,RY,.03,.06,.TRUE.)
!@@WAD@@
        CALL SECOM_CHECK_TEAM(VMAX,EMIN,IIT,JIT,RV,RE)
!$OMP END PARALLEL
        WRITE(22) INTX,IEXT,VMAX,EMIN,IIT,JIT
!@@CAPTURE@@
       ENDDO
      ENDDO
      CLOSE(22)
      CALL SMPI_FINALIZE()
      END
      SUBROUTINE OCLAMP(X)
      STOP 73
      END
      SUBROUTINE PCLAMP
      STOP 73
      END
      SUBROUTINE RANDB
      STOP 73
      END
      SUBROUTINE MIXED(I)
      STOP 73
      END
