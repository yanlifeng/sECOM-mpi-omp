      SUBROUTINE UPDATE_CBC
!NG04272010 UPDATE CBC BASED ON INSTANTANEOUS DEPTH DT(I,J)
      INCLUDE 'comdeck'
      SAVE

      Do 6030 J = 1, JM
        Do 6020 I = 1, IM
          If (FSM(I,J).NE.0.0) Then
!
!  VARIABLE BOTTOM FRICTION FOR SEDIMENT TRANSPORT
!
            IF (SEDTRAN.EQ.'INCLUDE') THEN
!
!  COHESIVE ELEMENTS
!
              IF (IBMSK(I,J).EQ.0) THEN
                Z0 = Z0BCOH
                CBCMIN = BFCOH
                IF (TOR.EQ.'BAROTROPIC') THEN
                  CBC(I,J)=BFCOH*FSM(I,J)
                ELSE
                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
     +             AMAX1(DT(I,J),Z0+5.E-1,WETMIN)/Z0)**2)*FSM(I,J)
                ENDIF
!
!  NON-COHESIVE ELEMENTS
!
              ELSE
                Z0 = Z0B
                CBCMIN = BFRIC
                IF (TOR.EQ.'BAROTROPIC') THEN
                  CBC(I,J)=BFRIC*FSM(I,J)
                ELSE
                  CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
     +             AMAX1(DT(I,J),Z0+5.E-1,WETMIN)/Z0)**2)*FSM(I,J)
                 ENDIF
              ENDIF
            ELSE
              Z0 = Z0B
              CBCMIN = BFRIC
              IF (TOR.EQ.'BAROTROPIC') THEN
                CBC(I,J)=BFRIC*FSM(I,J)
              ELSE
                CBC(I,J)=AMAX1(CBCMIN,.16/ALOG((DZ(KBM1)*0.5)*
     +             AMAX1(DT(I,J),Z0+5.E-1,WETMIN)/Z0)**2)*FSM(I,J)
              ENDIF
            ENDIF
          Else
           CBC(I,J)=BFRIC !Initialize CBC even at FSM=0 cells 05302012
          END IF
          CBC(I,J) = AMIN1(CBC(I,J)*VARBF(I,J),1.) ! FACTOR IN 2D VARIABLE BFRIC, AND DO NOT ALLOW CBC TO EXCEED UNITY.
 6020   CONTINUE
 6030 CONTINUE

      RETURN
      END
