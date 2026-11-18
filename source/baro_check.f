      RECURSIVE SUBROUTINE SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
! First reduce only values, making the per-cell loop vectorizable.
! Recover indices only on the two winning rows, preserving the legacy
! tie order and speed-before-minimum abort behavior (including all-land).
!$OMP DO SCHEDULE(STATIC)
      DO J=1,JMM1
        VROWMAX=-1.E10
!$OMP SIMD REDUCTION(MAX:VROWMAX) PRIVATE(VMAXTMP)
        DO I=1,IMM1
      VMAXTMP=SQRT((0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2)
          IF(VMAXTMP.GT.VROWMAX) VROWMAX=VMAXTMP
        ENDDO
        ROWVMAX(J)=VROWMAX
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO J=1,JM
        ZROWMIN=1.E10
!$OMP SIMD REDUCTION(MIN:ZROWMIN)
        DO I=1,IM
          IF(FSM(I,J).NE.0..AND.ELF(I,J).LT.ZROWMIN)
     &      ZROWMIN=ELF(I,J)
        ENDDO
        ROWEMIN(J)=ZROWMIN
      ENDDO
!$OMP END DO
!$OMP SINGLE
      VAMAX=-1.E10
      JVMAX=0
      DO J=1,JMM1
        IF(ROWVMAX(J).GE.VAMAX.AND.ROWVMAX(J).GE.0.) THEN
          JVMAX=J
          VAMAX=ROWVMAX(J)
        ENDIF
      ENDDO
      IF(JVMAX.NE.0) THEN
        J=JVMAX
        DO I=IMM1,1,-1
      VMAXTMP=SQRT((0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2)
          IF(VMAXTMP.EQ.VAMAX) THEN
            IAMAX=I
            JAMAX=J
            EXIT
          ENDIF
        ENDDO
      ENDIF
      IF(.NOT.(VAMAX.GT.100.)) THEN
      ELMINX=1.E10
      JWATER=0
      DO J=1,JM
        IF(ROWEMIN(J).LT.ELMINX) THEN
          JWATER=J
          ELMINX=ROWEMIN(J)
        ENDIF
      ENDDO
      IF(JWATER.NE.0) THEN
        J=JWATER
        DO I=1,IM
          IF(FSM(I,J).NE.0..AND.ELF(I,J).EQ.ELMINX) THEN
            IAMAX=I
            JAMAX=J
            ELMINX=ELF(I,J)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      ENDIF
!$OMP END SINGLE
      END
