      SUBROUTINE WANDD_KKY(ELEV,TD)
      INCLUDE 'comdeck'
      REAL ELEV(IM,JM),TD(IM,JM)
!$OMP PARALLEL DEFAULT(SHARED)
      CALL WANDD_TEAM(ELEV,TD)
!$OMP END PARALLEL
      END
      RECURSIVE SUBROUTINE WANDD_TEAM(ELEV,TD)
      INCLUDE 'comdeck'
      REAL ELEV(IM,JM),TD(IM,JM),TMPV41(1),TMPE41(1)
! Public and nonpersistent callers never request a diagnostic.
      CALL WANDD_ROWS41(ELEV,TD,.FALSE.,TMPV41,TMPE41)
      END
      RECURSIVE SUBROUTINE WANDD_ROWS41(ELEV,TD,SCAN41,RV41,RE41)
!
!     INTRODUCE WET/DRY SCHEME BASED ON KKY's DIISERTATION (hli, 08/03/05).
!NGMAR2009   SIMPLIFIED AND ADDED AN EXTRA CONDITION (NO FLOW OUT OF "DRY" CELL).
!
!
      INCLUDE 'comdeck'
!
      DIMENSION ELEV(IM,JM),TD(IM,JM)
      LOGICAL SCAN41
      REAL RV41(*),RE41(*)
!
!     Owner (i,j) handles both tests of its east and north faces. Its
!     ET/ETF scatter may also touch row j+1. Alternate row colors with
!     a barrier so there are no concurrent writes (even equal writes).
!     The sole caller supplies EL,D, neither aliases the written fields.
      DO ICOLOR=1,2
!$OMP DO SCHEDULE(STATIC)
      DO J = ICOLOR, JMM1, 2
         DO I = 1, IMM1
               IF(TD(I,J).LE.WETMIN.AND.FSM(I,J).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I+1,J)).LE.WETMIN).OR.
     &            ((ELEV(I+1,J)-ELEV(I,J)).LE.WETEPS).OR.
     *             (UA(I+1,J).GT.0).OR.DUMWAD(I+1,J).EQ.0.0) THEN
                     UA(I+1,J)=0.0
                     UAF(I+1,J)=0.0
                     DUMWAD(I+1,J)=0.0
                     ET(I,J)=ETB(I,J) ! NG04012011
                     ETF(I,J)=ETB(I,J) ! NG04012011
                  ENDIF
                  IF((0.5*(TD(I,J)+TD(I,J+1)).LE.WETMIN).OR.
     &            ((ELEV(I,J+1)-ELEV(I,J)).LE.WETEPS).OR.
     *             (VA(I,J+1).GT.0).OR.DVMWAD(I,J+1).EQ.0.0) THEN
                     VA(I,J+1)=0.0
                     VAF(I,J+1)=0.0
                     DVMWAD(I,J+1)=0.0
                     ET(I,J)=ETB(I,J) ! NG04012011
                     ETF(I,J)=ETB(I,J) ! NG04012011
                  ENDIF
!                  write(*,*)'W&D1 ',time,I,J,TD(I,J),UA(I+1,J),VA(I,J+1)
               ENDIF
!
               IF(TD(I+1,J).LE.WETMIN.AND.DUM(I+1,J).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I+1,J)).LE.WETMIN).OR.
     &            ((ELEV(I,J)-ELEV(I+1,J)).LE.WETEPS).OR.
     *             (UA(I+1,J).LT.0).OR.DUMWAD(I+1,J).EQ.0.0) THEN
                     UA(I+1,J)=0.0
                     UAF(I+1,J)=0.0
                     DUMWAD(I+1,J)=0.0
                     ET(I+1,J)=ETB(I+1,J) ! NG04012011
                     ETF(I+1,J)=ETB(I+1,J) ! NG04012011
!                  write(*,*)'W&D2 ',time,I+1,J,TD(I+1,J),UA(I+1,J)
                  ENDIF
               ENDIF
!
               IF(TD(I,J+1).LE.WETMIN.AND.DVM(I,J+1).NE.0) THEN
                  IF((0.5*(TD(I,J)+TD(I,J+1)).LE.WETMIN).OR.
     &            ((ELEV(I,J)-ELEV(I,J+1)).LE.WETEPS).OR.
     *             (VA(I,J+1).LT.0).OR.DVMWAD(I,J+1).EQ.0.0) THEN
                     VA(I,J+1)=0.0
                     VAF(I,J+1)=0.0
                     DVMWAD(I,J+1)=0.0
                     ET(I,J+1)=ETB(I,J+1) ! NG04012011
                     ETF(I,J+1)=ETB(I,J+1) ! NG04012011
!                  write(*,*)'W&D3 ',time,I,J+1,TD(I,J+1),VA(I,J+1)
                  ENDIF
               ENDIF
         ENDDO
! U(:,J) and V(:,J+1) can only be written by this face-owner row.
! First odd row uses an untouched physical south face. After odd
! publication, even J finalizes the last faces for centers J and J+1.
         IF(SCAN41) THEN
          IF(J.EQ.1) CALL SECOM_PAIR41(J,J,RV41,RE41)
          IF(ICOLOR.EQ.2) CALL SECOM_PAIR41(J,MIN(J+1,JMM1),RV41,RE41)
          IF(J.EQ.JMM1) CALL SECOM_WATER41(JM,RE41)
         ENDIF
      ENDDO

!$OMP END DO
      ENDDO
      RETURN
      END
!
! Only called with fully finalized face data for this one row / adjacent
! pair. Evaluate every norm. Never replace a threshold test by a bound.
      RECURSIVE SUBROUTINE SECOM_PAIR41(J1,J2,RV41,RE41)
      INCLUDE 'comdeck'
      REAL RV41(*),RE41(*),R1,R2,V1,V2,VMID
      R1=-1.E10
      IF(J1.EQ.J2) THEN
!$OMP SIMD REDUCTION(MAX:R1) PRIVATE(V1)
       DO I=1,IMM1
        V1=(0.5*DUMWAD(I,J1)*UAF(I,J1)+
     .      0.5*DUMWAD(I+1,J1)*UAF(I+1,J1))**2+
     .     (0.5*DVMWAD(I,J1)*VAF(I,J1)+
     .      0.5*DVMWAD(I,J1+1)*VAF(I,J1+1))**2
        IF(V1.GT.R1) R1=V1
       ENDDO
      ELSE
       R2=-1.E10
!$OMP SIMD REDUCTION(MAX:R1,R2) PRIVATE(V1,V2,VMID)
       DO I=1,IMM1
        VMID=0.5*DVMWAD(I,J2)*VAF(I,J2)
        V1=(0.5*DUMWAD(I,J1)*UAF(I,J1)+
     .      0.5*DUMWAD(I+1,J1)*UAF(I+1,J1))**2+
     .     (0.5*DVMWAD(I,J1)*VAF(I,J1)+VMID)**2
        V2=(0.5*DUMWAD(I,J2)*UAF(I,J2)+
     .      0.5*DUMWAD(I+1,J2)*UAF(I+1,J2))**2+
     .     (VMID+0.5*DVMWAD(I,J2+1)*VAF(I,J2+1))**2
        IF(V1.GT.R1) R1=V1
        IF(V2.GT.R2) R2=V2
       ENDDO
       RV41(J2)=R2
       IF(R2.GE.0.) RV41(J2)=SQRT(R2)
       CALL SECOM_WATER41(J2,RE41)
      ENDIF
      RV41(J1)=R1
      IF(R1.GE.0.) RV41(J1)=SQRT(R1)
      CALL SECOM_WATER41(J1,RE41)
      END
      RECURSIVE SUBROUTINE SECOM_WATER41(J,RE41)
      INCLUDE 'comdeck'
      REAL RE41(*)
      ZROWMIN=1.E10
!$OMP SIMD REDUCTION(MIN:ZROWMIN)
      DO I=1,IM
        IF(FSM(I,J).NE.0..AND.ELF(I,J).LT.ZROWMIN)
     &    ZROWMIN=ELF(I,J)
      ENDDO
      RE41(J)=ZROWMIN
      END
