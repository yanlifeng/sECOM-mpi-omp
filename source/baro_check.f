      RECURSIVE SUBROUTINE SECOM_CHECK_TEAM(VAMAX,ELMINX,IAMAX,JAMAX,ROWVMAX,ROWEMIN)
      INCLUDE 'comdeck'
      REAL ROWVMAX(JM),ROWEMIN(JM)
! Wet/dry owns and commits all face changes first. Only these two
! dynamic J+1 inputs cross a row owner for the speed calculation.
      IF(SMPI_SIZE.GT.1) THEN
        CALL SMPI_HALO(DVMWAD,IM,JM)
        CALL SMPI_HALO(VAF,IM,JM)
      ENDIF
! First reduce only values, making the per-cell loop vectorizable.
! Recover indices only on the two winning rows, preserving the legacy
! tie order and speed-before-minimum abort behavior (including all-land).
!$OMP DO SCHEDULE(STATIC)
      DO J=SMPI_FIRST(1),SMPI_LAST(JMM1)
        VROWMAX=-1.E10
!$OMP SIMD REDUCTION(MAX:VROWMAX) PRIVATE(VMAXTMP)
        DO I=1,IMM1
! SQRT is monotone on a finite nonnegative sum of squares. Reduce the
! original rounded squared norm, then round the root once per row.
! The final scan below still compares rounded roots: distinct squared
! norms may tie after SQRT, even in different rows or at the abort limit.
      VMAXTMP=(0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2
          IF(VMAXTMP.GT.VROWMAX) VROWMAX=VMAXTMP
        ENDDO
        ROWVMAX(J)=VROWMAX
        IF(VROWMAX.GE.0.) ROWVMAX(J)=SQRT(VROWMAX)
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
      DO J=SMPI_FIRST(1),SMPI_LAST(JM)
! SMPI_ROWS transfers JM entries of each vector; JM has no speed cell.
        IF(SMPI_SIZE.GT.1.AND.J.EQ.JM) ROWVMAX(J)=-1.E10
        ZROWMIN=1.E10
!$OMP SIMD REDUCTION(MIN:ZROWMIN)
        DO I=1,IM
          IF(FSM(I,J).NE.0..AND.ELF(I,J).LT.ZROWMIN)
     &      ZROWMIN=ELF(I,J)
        ENDDO
        ROWEMIN(J)=ZROWMIN
      ENDDO
!$OMP END DO
      IF(SMPI_SIZE.GT.1) CALL SMPI_ROWS(ROWVMAX,ROWEMIN,JM)
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
        IPICK=0
        IF(SMPI_OWNS(J)) THEN
        DO I=IMM1,1,-1
      VMAXTMP=SQRT((0.5*DUMWAD(I,J)*UAF(I,J)+
     .              0.5*DUMWAD(I+1,J)*UAF(I+1,J))**2+
     .             (0.5*DVMWAD(I,J)*VAF(I,J)+
     .              0.5*DVMWAD(I,J+1)*VAF(I,J+1))**2)
          IF(VMAXTMP.EQ.VAMAX) THEN
            IPICK=I
            EXIT
          ENDIF
        ENDDO
        ENDIF
! A remote row's dense data need not exist locally. Only the owner
! reports its descending-I winner. Keep old coordinates if no match.
        IF(SMPI_SIZE.GT.1) CALL SMPI_INDEX_RAW(IPICK)
        IF(IPICK.GT.0) THEN
          IAMAX=IPICK
          JAMAX=JVMAX
        ENDIF
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
        IPICK=0
        IF(SMPI_OWNS(J)) THEN
        DO I=1,IM
          IF(FSM(I,J).NE.0..AND.ELF(I,J).EQ.ELMINX) THEN
            IPICK=I
            ELMINX=ELF(I,J)
            EXIT
          ENDIF
        ENDDO
        ENDIF
        IF(SMPI_SIZE.GT.1) CALL SMPI_INDEX_RAW(IPICK)
        IF(IPICK.GT.0) THEN
          IAMAX=IPICK
          JAMAX=JWATER
! Preserve the *matched cell*, including its original sign of zero.
          IF(SMPI_SIZE.GT.1) CALL SMPI_RBCAST_ROW_RAW(ELMINX,JWATER)
        ENDIF
      ENDIF
      ENDIF
!$OMP END SINGLE
      END
