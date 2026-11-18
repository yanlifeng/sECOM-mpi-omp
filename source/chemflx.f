      SUBROUTINE CHEMFLX
!
!  CALCULATES DEPOSITION AND DEPOSITION FLUXES OF CHEM AT
!  SEDIMENT-WATER INTERFACE
!
!  REVISION DATE :  APRIL 18, 1996
!
!**********************************************************************
!
      INCLUDE 'comdeck'
!
!
!  CLASS 1 & 2 SEDIMENTS:  FLUX IN ug CHEM/cm**2
!
!  CBEDCHEM = BED CONC. IN LAYER 1 (ug CHEM/g SOLIDS)
!
!  COHESIVE ELEMENTS
!
      DO 10 J=2,JMM1
        DO 10 I=2,IMM1
          CHEMBOT1(I,J)=0.0
          CHEMBOT2(I,J)=0.0
!
          IF (FSM(I,J).GT.0.0) THEN
!
!  EROSION FROM LAYER 1
!
            CHEMBOT1(I,J)=CBEDCHEM(1,I,J)*E(1,I,J)
            CHEMBOT2(I,J)=CBEDCHEM(1,I,J)*E(2,I,J)
!
!  DEPOSITION
!
            IF (CSED1(I,J,KBM1).GT.0.0) THEN
              CHEMBOT1(I,J)=CHEMBOT1(I,J)-
     +                   DD(1,I,J)*CHEM1(I,J,KBM1)/CSED1(I,J,KBM1)
            ENDIF
!
            IF (CSED2(I,J,KBM1).GT.0.0) THEN
              CHEMBOT2(I,J)=CHEMBOT2(I,J)-
     +                   DD(2,I,J)*CHEM2(I,J,KBM1)/CSED2(I,J,KBM1)
            ENDIF
          ENDIF
 10   CONTINUE
!
      RETURN
      END
