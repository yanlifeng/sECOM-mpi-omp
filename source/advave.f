      SUBROUTINE ADVAVE(ADVUA,ADVVA)
!     VERSION(09/27/90)
!
      INCLUDE 'comdeck'
      SAVE
!
!-----------------------------------------------------------------------
!      THIS SUBROUTINE CALCULATES ADVECTION & DIFFUSION
!-----------------------------------------------------------------------
!
      DIMENSION ADVUA(IM,JM),ADVVA(IM,JM)
!
      REAL UFX(IM,JM),UFY(IM,JM),XU0,YU0,XV0,YV0
! Build both components at the same cell. U fluxes have their own
! scratch: the legacy shared FLUX arrays retain the final V component.
! Only divergence needs neighboring rows and hence a full-team barrier.
!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J,XU0,YU0,XV0,YV0)
      DO J=1,JM
       DO I=1,IM
        XU0=0.0
        YU0=0.0
        XV0=0.0
        YV0=0.0
        TPS(I,J)=0.0
        IF(I.GE.2.AND.J.GE.2) THEN
         IF(I.LE.IMM1) THEN
          IF(J.LE.JMM1) THEN
           IF(ADVECT.EQ.'NON-LINEAR') THEN
      XU0=.125*((D(I+1,J)+D(I,J))*UA(I+1,J)
     .                 +(D(I,J)+D(I-1,J))*UA(I,J))
     .                  *(UA(I+1,J)+UA(I,J))
      YV0=.125*((D(I,J+1)+D(I,J))
     .       *VA(I,J+1)+(D(I,J)+D(I,J-1))*VA(I,J))
     .      *(VA(I,J+1)+VA(I,J))
           ENDIF
      XU0=XU0
     .         -D(I,J)*2.*AAM2D(I,J)*(UAB(I+1,J)-UAB(I,J))/H1(I,J)
      YV0=YV0
     .        -D(I,J)*2.*AAM2D(I,J)*(VAB(I,J+1)-VAB(I,J))/H2(I,J)
          ENDIF
          IF(ADVECT.EQ.'NON-LINEAR') THEN
      YU0=.125*((D(I,J)+D(I,J-1))*VA(I,J)
     .                 +(D(I-1,J)+D(I-1,J-1))*VA(I-1,J))
     .                 *(UA(I,J)+UA(I,J-1))
          ENDIF
      TPS(I,J)=.25*(D(I,J)+D(I-1,J)+D(I,J-1)+D(I-1,J-1))
     .            *(AAM2D(I,J)+AAM2D(I,J-1)+AAM2D(I-1,J)+AAM2D(I-1,J-1))
     .                *((UAB(I,J)-UAB(I,J-1))
     .                /(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
     .                 +(VAB(I,J)-VAB(I-1,J))
     .               /(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1)) )
      XU0=XU0*H2(I,J)
      YU0=(YU0-TPS(I,J)*DUM(I,J)*DUM(I,J-1)
     .            *DUMWAD(I,J)*DUMWAD(I,J-1)*DVMWAD(I,J)*DVMWAD(I-1,J) ! Not in original NYHOPS W&D
     .                                 *DVM(I,J)*DVM(I-1,J))! THIN DAMS NG08192008
     .            *.25*(H1(I,J)+H1(I-1,J)+H1(I,J-1)+H1(I-1,J-1))
         ENDIF
         IF(J.LE.JMM1) THEN
          IF(J.GE.3.AND.ADVECT.EQ.'NON-LINEAR') THEN
      XV0=.125*((D(I,J)+D(I-1,J))*UA(I,J)
     .         +(D(I,J-1)+D(I-1,J-1))*UA(I,J-1))*
     .                        (VA(I-1,J)+VA(I,J))
          ENDIF
      YV0=YV0*H1(I,J)
      XV0=(XV0-TPS(I,J)*DVM(I,J)*DVM(I-1,J)
     .            *DVMWAD(I,J)*DVMWAD(I-1,J)*DUMWAD(I,J)*DUMWAD(I,J-1) ! Not in original NYHOPS W&D
     .                                 *DUM(I,J)*DUM(I,J-1))! THIN DAMS NG08192008
     .             *.25*(H2(I,J)+H2(I-1,J)+H2(I,J-1)+H2(I-1,J-1))
         ENDIF
        ENDIF
        UFX(I,J)=XU0
        UFY(I,J)=YU0
        FLUXUA(I,J)=XV0
        FLUXVA(I,J)=YV0
       ENDDO
      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC) PRIVATE(I,J)
      DO J=1,JM
       DO I=1,IM
        ADVUA(I,J)=0.0
        ADVVA(I,J)=0.0
        IF(J.GE.2.AND.J.LE.JMM1.AND.I.GE.3.AND.I.LE.IMM1) THEN
      ADVUA(I,J)=UFX(I,J)-UFX(I-1,J)
     .           +UFY(I,J+1)-UFY(I,J)
        ENDIF
        IF(J.GE.3.AND.J.LE.JMM1.AND.I.GE.2.AND.I.LE.IMM1) THEN
      ADVVA(I,J)=FLUXUA(I+1,J)-FLUXUA(I,J)
     .          +FLUXVA(I,J)-FLUXVA(I,J-1)
        ENDIF
       ENDDO
      ENDDO
!$OMP END DO
!
! Keep potentially overlapping boundary entries in original order.
!$OMP SINGLE
      IF (OPTEBC(15:20).EQ.'LINEAR') THEN ! IF LINEAR OPTION ENABLED AT THE END OF OPTEBC
       DO 309 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        IC=ICON(N)
        JC=JCON(N)
        IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE,JE)=0.0
         ADVVA(IE,JE)=0.0   ! LINEARIZE THIS AND BCOND IF NEEDED
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE+1,JE)=0.0
         ADVVA(IE+1,JE)=0.0 ! LINEARIZE THIS AND BCOND IF NEEDED
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE)=0.0
         ADVUA(IE,JE)=0.0   ! LINEARIZE THIS AND BCOND IF NEEDED
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE+1)=0.0
         ADVUA(IE,JE+1)=0.0 ! LINEARIZE THIS AND BCOND IF NEEDED
        ENDIF
 309   CONTINUE
      ELSE
       DO 310 N=1,NUMEBC
        IE=IETA(N)
        JE=JETA(N)
        IC=ICON(N)
        JC=JCON(N)
        IF(FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE,JE)=0.0
        ELSE IF(FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) THEN
         ADVUA(IE+1,JE)=0.0
        ELSE IF(FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE)=0.0
        ELSE IF(FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) THEN
         ADVVA(IE,JE+1)=0.0
        ENDIF
 310   CONTINUE
      ENDIF
!$OMP END SINGLE
!
      RETURN
      END
