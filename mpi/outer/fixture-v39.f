      PROGRAM CHECK_ARCHIVE
      USE ISO_FORTRAN_ENV, ONLY: INT32,INT64
      INCLUDE 'comdeck'
      LOGICAL PROOF,SECOM_HIST_CERT0
      CHARACTER*10 ARG
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) MODE
      TOR='BAROTROPIC'
      SEDTRAN='NEGLECT'
      CHEMTRAN='NEGLECT'
      PARTICLE='NEGLECT'
      NCHEMLAY=0
      EPTS=0
      FPTS=0
      VPTS=0
      ISKILL=10000
      SKILLI=1.0/FLOAT(ISKILL)
      AREA=1.0
      IF (MODE.EQ.9) AREA=0.
      JHM=3
      IAVGE=3
      DEI=1./FLOAT(IAVGE)
      DO IH=1,JHM
        IHIST(IH,1)=1+(IH-1)*IAVGE
        IHIST(IH,2)=IH*IAVGE
      ENDDO
      DAYI=1./86400.
      DTI=1.
      CDFOUT=.TRUE.
      CONSPLT=.FALSE.
      IBMSK=0
      DO J=1,JM
      DO I=1,IM
        IBMSK(I,J)=MOD(I+3*J,3)
        H(I,J)=I/7.+J/11.
        H1(I,J)=I/3.+J/5.
        H2(I,J)=I/11.+J/13.
      ENDDO
      ENDDO
      WSET1=0.
      WSET2=0.
      CSED1=0.
      CSED2=0.
      CHEM1=0.
      CHEM2=0.
      TSED=0.D0
      BEDTH=0.D0
      ARCWS1=0.
      ARCWS2=0.
      ARCSED1=0.
      ARCSED2=0.
      ARCHEM1=0.
      ARCHEM2=0.
      ARCTHIK=0.D0
      IF (MODE.EQ.1) THEN
        CSED2(IM,JM,1)=1.25
        TSED(7,3,JM)=9.25D0
        BEDTH(1,1,JM)=-0.125D0
      ENDIF
      IF (MODE.EQ.2) WSET1(3,JM,1)=TRANSFER(IBSET(0_INT32,31),0.)
      IF (MODE.EQ.3) CHEM2(IM,JM,1)=TRANSFER(INT(Z'7F800123',INT32),0.)
      IF (MODE.EQ.4) THEN
        DEI=1.
        WSET2(1,1,1)=HUGE(DEI)
        ARCWS2(1,1,1)=HUGE(DEI)
      ENDIF
      IF (MODE.EQ.5) THEN
! Unused opposite branches contain toxic bits; certificate must not
! touch TSED for 1/other masks nor BEDTH for 0/other masks.
        TSED(7,1,JM)=TRANSFER(INT(Z'7FF0000000000001',INT64),0.D0)
        BEDTH(1,3,JM)=TRANSFER(INT(Z'7FF0000000000001',INT64),0.D0)
      ENDIF
      IF (MODE.EQ.6) ARCHEM1(IM,JM,1)=TRANSFER(IBSET(0_INT32,31),0.)
      IF (MODE.EQ.7) CSED2(7,8,1)=7.5
      IF (MODE.EQ.8) DEI=-DEI
      PROOF=.FALSE.
!@@CERT@@
      OPEN(23,FILE='archive.raw',STATUS='NEW',FORM='UNFORMATTED',ACCESS='STREAM')
      DO INTX=1,9
        TIME=DAYI*INTX
        DO J=1,JM
        DO I=1,IM
          EL(I,J)=I*.01+J*.0003+INTX*.1
          UAF(I,J)=I*.2-J*.3+INTX*.4
          VAF(I,J)=I*.3-J*.2-INTX*.5
          AAM2D(I,J)=1.+EL(I,J)
          TAU(I,J,1)=EL(I,J)*7.
          WSX2D(I,J)=UAF(I,J)
          RELHU2D(I,J)=VAF(I,J)
          PATM(I,J)=REAL(I+J)
        ENDDO
        ENDDO
!@@CALL@@
        CALL SNAP39
      ENDDO
      CLOSE(23)
      END
      SUBROUTINE SNAP39
      INCLUDE 'comdeck'
      WRITE(23) INTX,ARCU,ARCV,ARCUX,ARCVX,ARCET,ARCAM,ARCPATM,
     &   ARCWSX2D,ARCWSY2D,ARCCLD,ARCRHUM,ARCAIRT,ARCSWOBS,ARCQPREC,ARCQEVAP,
     &   ARCWVPD,ARCWVHT,ARCWVDR,ARCWS1,ARCWS2,ARCSED1,ARCSED2,
     &   ARCTAU,ARCTHIK,ARCHEM1,ARCHEM2
      END
      SUBROUTINE PUTCDF(TMIDDLE)
      WRITE(23) TMIDDLE
      CALL SNAP39
      END
      SUBROUTINE TSCDF(TMIDDLE,I)
      STOP 94
      END
