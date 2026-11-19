      SUBROUTINE UPDATE_CBCM
!NG2014 This version if for barotropic only with Manning's n.
      INCLUDE 'comdeck'
      INTEGER JOCLO,JOCHI
      SAVE

      JOCLO=1
      JOCHI=JM
      IF(SMPI_OUTER) THEN
        JOCLO=SMPI_JLO
        JOCHI=SMPI_JHI
      END IF
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,CBCMIN) SCHEDULE(STATIC)
      Do 6030 J = JOCLO, JOCHI
        Do 6020 I = 1, IM
          If (FSM(I,J).NE.0.0) Then
              CBCMIN = BFRIC
              CBC(I,J)=AMAX1(CBCMIN,GRAV*VARMN(I,J)*VARMN(I,J)/
     +                  (AMAX1(DT(I,J),WETMIN)**(1./3.)))
          Else
           CBC(I,J)=BFRIC !Initialize CBC even at FSM=0 cells
          END IF
          CBC(I,J) = AMIN1(CBC(I,J),1.) ! DO NOT ALLOW CBC TO EXCEED UNITY.
 6020   CONTINUE
 6030 CONTINUE
! EXTRNL_TEAM computes TYBOT using CBC(I,J-1); all remaining local
! consumers of freshly computed CBC use the same J.
      IF(SMPI_OUTER) CALL SMPI_HALO_RAW(CBC,IM,JM)

      RETURN
      END
