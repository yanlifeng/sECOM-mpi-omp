      PROGRAM Q34
      USE OMP_LIB
      INCLUDE 'comdeck'
      REAL AX(IM,JM),AY(IM,JM),XN(IM,JM),YN(IM,JM)
      LOGICAL QROT
      INTEGER QFLAG,QNT
      COMMON /Q34FLAGS/ QFLAG
      CHARACTER*8 SARG
      CALL GET_COMMAND_ARGUMENT(1,SARG)
      READ(SARG,*) QNT
      OPEN(91,FILE='fields.bin',ACCESS='STREAM',FORM='UNFORMATTED',STATUS='NEW')
      RHO0=1025.
      PATMOPT=1.
      DTE=.02
      GRAV=9.81
      RAMP=.8
      NU=.105
      THETA=.21
      NUMDBC=0
      NUMEPC=0
      OPTEPC=' '
      NUMEBC=7
      IETA(1:7)=(/10,11,20,30,10,IM-1,50/)
      JETA(1:7)=(/10,10,20,30,10,50,JM-1/)
      ICON(1:7)=(/9,10,20,30,9,IM,50/)
      JCON(1:7)=(/10,10,21,29,10,50,JM/)
      NUMECR=5
      NCHARD(1:5)=(/1,2,3,4,1/)
      ICHARD(1:5)=(/18,20,21,19,18/)
      JCHARD(1:5)=(/14,14,14,13,14/)
      DO IQ=0,23
        TOR='BAROTROPIC'
        IF(MOD(IQ/4,3).EQ.2) TOR='PROGNOSTIC'
        QROT=MOD(IQ/4,3).NE.1
        BCTYPE='CLAMPED'
        IF(IQ.GE.12) BCTYPE='PCLAMP '
        ADVECT='LINEAR'
        IF(MOD(IQ,2).EQ.0) ADVECT='NON-LINEAR'
        OPTEBC='DATA'
        IF(MOD(IQ/2,2).EQ.0) OPTEBC(15:20)='LINEAR'
        FSM=1.
        FSM(11,10)=0.
        FSM(9,10)=0.
        FSM(12,10)=0.
        FSM(20,19)=0.
        FSM(30,31)=0.
        FSM(IM-2,50)=0.
        FSM(50,JM-2)=0.
        DO J=1,JM
          DO I=1,IM
            H(I,J)=4.+MOD(I+J,9)*.03
            H1(I,J)=9.+MOD(I+J,5)*.27
            H2(I,J)=10.+MOD(I*J,9)*.31
            ART(I,J)=H1(I,J)*H2(I,J)
            ARU(I,J)=ART(I,J)*.91
            ARV(I,J)=ART(I,J)*.93
            QEVAP2D(I,J)=MOD(I+J,5)*.0003
            QPREC2D(I,J)=MOD(I+3*J,11)*.0002
            PATM(I,J)=MOD(I*J,9)*.15
            DATUM(I,J)=MOD(I+J,4)*.025
            CBC(I,J)=.003+MOD(I*J,7)*.001
            COR(I,J)=.00012+MOD(I+J,9)*.00001
            VARWIF(I,J)=.12
            AREAICE(I,J)=MOD(I+J,9)*.15
            AAM2D(I,J)=1.+MOD(I*J,4)*.29
            DUM(I,J)=1.
            IF(MOD(I+3*J,19).EQ.0) DUM(I,J)=0.
            DVM(I,J)=1.
            IF(MOD(2*I+J,17).EQ.0) DVM(I,J)=0.
            DUMWAD(I,J)=1.
            IF(MOD(I+J,5).EQ.0) DUMWAD(I,J)=0.
            DVMWAD(I,J)=1.
            IF(MOD(I+J,7).EQ.0) DVMWAD(I,J)=0.
            UA(I,J)=(MOD(I+2*J,7)-3)*.113
            VA(I,J)=(MOD(4*I+J,9)-4)*.087
            UAB(I,J)=VA(I,J)
            VAB(I,J)=UA(I,J)
            UAF(I,J)=VAB(I,J)
            VAF(I,J)=UAB(I,J)
            XN(I,J)=.018*VA(I,J)
            YN(I,J)=.016*UA(I,J)
            ELB(I,J)=.1+MOD(I+J,7)*.02
            EL(I,J)=.2+MOD(I+J,13)*.03
            ELF(I,J)=EL(I,J)
            D(I,J)=H(I,J)+EL(I,J)
            TXBOT(I,J)=MOD(I+3*J,15)*.032
            TYBOT(I,J)=-MOD(I+3*J,15)*.071
            TXSURF(I,J)=.03*UA(I,J)
            TYSURF(I,J)=.03*VA(I,J)
            TXICE(I,J)=.1
            TYICE(I,J)=.2
          ENDDO
        ENDDO
        DO KQSTEP=1,2
          IEXT=KQSTEP
          QFLAG=0
!$OMP PARALLEL DEFAULT(SHARED)
          IF(OMP_GET_NUM_THREADS().NE.QNT) STOP 92
          CALL EXTRNL_TEAM(AX,AY,XN,YN,.04,.4,QROT)
!$OMP END PARALLEL
          WRITE(91) AX,AY,FLUXUA,FLUXVA,TPS,CURV42D,TXBOT,TYBOT,TXICE,TYICE,
     &              UA,VA,UAF,VAF,UAB,VAB,ELB,ELF,EL,D
        ENDDO
      ENDDO
      CLOSE(91)
      END
! Only unrelated BC library calls are replaced by deterministic fixture hooks.
      SUBROUTINE BCOND(IDX,X,K)
      INCLUDE 'comdeck'
      IF(IDX.EQ.12) ELF(2,JM-1)=ELF(3,JM-1)+.031
      IF(IDX.EQ.13) VAF(2,JM-1)=UAF(3,JM-1)-.015
      END
      SUBROUTINE PCLAMP
      INCLUDE 'comdeck'
      ELF(3,4)=ELF(3,5)*.921
      END
      SUBROUTINE OCLAMP(X)
      STOP 81
      END
      SUBROUTINE RANDB
      STOP 82
      END
      SUBROUTINE MIXED(K)
      STOP 83
      END
      SUBROUTINE DELAY34
      USE OMP_LIB
      DOUBLE PRECISION T
      IF(OMP_GET_NUM_THREADS().EQ.1) RETURN
      T=OMP_GET_WTIME()
      DO WHILE(OMP_GET_WTIME()-T.LT..025D0)
      ENDDO
      END
      SUBROUTINE PUBLISH34
      INTEGER QFLAG
      COMMON /Q34FLAGS/ QFLAG
!$OMP ATOMIC WRITE
      QFLAG=1
      END
      SUBROUTINE REQUIRE34
      INTEGER QFLAG,V
      COMMON /Q34FLAGS/ QFLAG
!$OMP ATOMIC READ
      V=QFLAG
      IF(V.NE.1) STOP 94
      END
