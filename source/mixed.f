!        SUBROUTINE TO FOR A MIXED OBC: ZERO VELOCITY GRADIENTS
!        OUT TO WATER DEPTHS OF 50 M, OPTIMIZED CLAMPED IN WATER
!        DEPTHS OF > 100 M, AND A GRADUAL TRANSITION
!        BETWEEN THE TWO OBC'S FROM 50 TO 100 M DEPTHS
      SUBROUTINE MIXED(NN)
      INCLUDE 'comdeck'
      Dimension DH(IM,JM)
!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
      GO TO (10,20) NN
 10   CONTINUE
!        GO THROUGH OPEN BOUNDARIES, CALC AMULT, RECALCULATE UAF
!        AND VAF
      DO N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        DEP=(H(IE,JE)+H(IC,JC))/2.
        AMULT=SQRT(AMAX1(0.,AMIN1(DEP-50.,50.)/50.))
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
! right
          HL=(H(IE-2,JE)+ELF(IE-2,JE)+
     1        H(IC,JC)+ELF(IC,JC))/2.
          HR=(H(IC,JC)+ELF(IC,JC)+H(IE,JE)+ELF(IE,JE))/2.
          UAF(IE,JE) = AMULT * UAF(IE,JE) + (1.-AMULT) *
     1      UAF(IC,JC)*(H2(IE-2,JE)+H2(IC,JC))*HL*
     *      DUM(IC,JC)/(HR*(H2(IE,JE)+H2(IC,JC)))
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
! left
          HL=(H(IC,JC)+ELF(IC,JC)+H(IE,JE)+ELF(IE,JE))/2.
          HR=(H(IC+1,JC)+ELF(IC+1,JC)+
     1        H(IC,JC)+ELF(IC,JC))/2.
          UAF(IC,JC) = AMULT * UAF(IC,JC) + (1.-AMULT) *
     1      UAF(IC+1,JC)*(H2(IC,JC)+H2(IC+1,JC))*HR*
     *      DUM(IC+1,JC)/(HL*(H2(IE,JE)+H2(IC,JC)))
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then 
! top
          HL=(H(IE,JE-2)+ELF(IE,JE-2)+
     1        H(IC,JC)+ELF(IC,JC))/2.
          HR=(H(IC,JC)+ELF(IC,JC)+H(IE,JE)+ELF(IE,JE))/2.
          VAF(IE,JE) = AMULT * VAF(IE,JE) + (1.-AMULT) *
     1      VAF(IC,JC)*(H1(IE,JE-2)+H1(IC,JC))*HL*
     *      DVM(IC,JC)/(HR*(H1(IE,JE)+H1(IC,JC)))
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then 
! bottom
          HL=(H(IC,JC)+ELF(IC,JC)+H(IE,JE)+ELF(IE,JE))/2.
          HR=(H(IC,JC+1)+ELF(IC,JC+1)+
     1        H(IC,JC)+ELF(IC,JC))/2.
          VAF(IC,JC) = AMULT * VAF(IC,JC) + (1.-AMULT) *
     1      VAF(IC,JC+1)*(H1(IC,JC)+H1(IC,JC+1))*HR*
     *      DVM(IC,JC+1)/(HL*(H1(IE,JE)+H1(IC,JC)))
        End If
      END DO
      RETURN
!        HANDLE BAROCLINIC CURRENTS
 20   CONTINUE
        PRINT *,'MIXED: NEED TO MODIFY FOR BAROCLINIC CURRENTS'
!hli        PAUSE
!        GO THROUGH OPEN BOUNDARIES, CALC AMULT, RECALCULATE UF
!        AND VF
      DO N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        DEP=(H(IE,JE)+H(IC,JC))/2.
        AMULT=SQRT(AMAX1(0.,AMIN1(DEP-50.,50.)/50.))
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
! right
          DO K=1,KBM1
            HL=(H(IE-2,JE)+ELF(IE-2,JE)+
     1        H(IC,JC)+ELF(IC,JC))/2.
            HR=H(IC,JC)+ELF(IC,JC)
            UF(IE,JE,K)=AMULT * UF(IE,JE,K) + (1.-AMULT) *
     *        UF(IC,JC,K)*(H2(IE-2,JE)+H2(IC,JC))*HL
     *        /(HR*(H2(IE,JE)+H2(IC,JC)))
          END DO
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
! left
          DO K=1,KBM1
            HL=H(IC,JC)+ELF(IC,JC)
            HR=(H(IC+1,JC)+ELF(IC+1,JC)+
     1        H(IC,JC)+ELF(IC,JC))/2.
            UF(IC,JC,K) = AMULT * UF(IC,JC,K)+(1.-AMULT)*
     *        UF(IC+1,JC,K)*(H2(IC,JC)+H2(IC+1,JC))*HR
     *        /(HL*(H2(IE,JE)+H2(IC,JC)))
          END DO
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then 
! top
          DO K=1,KBM1
            HL=(H(IE,JE-2)+ELF(IE,JE-2)+
     1        H(IC,JC)+ELF(IC,JC))/2.
            HR=H(IC,JC)+ELF(IC,JC)
            VF(IE,JE,K)=AMULT*VF(IE,JE,K)+(1.-AMULT)*
     *        VF(IC,JC,K)*(H1(IE,JE-2)+H1(IC,JC))*HL
     *        /(HR*(H1(IE,JE)+H1(IC,JC)))
          END DO
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then 
! bottom
          DO K=1,KBM1
            HL=H(IC,JC)+ELF(IC,JC)
            HR=(H(IC,JC+1)+ELF(IC,JC+1)+
     1        H(IC,JC)+ELF(IC,JC))/2.
            VF(IC,JC,K)=AMULT*VF(IC,JC,K)+(1.-AMULT)*
     *          VF(IC,JC+1,K)*(H1(IC,JC)+H1(IC,JC+1))*HR
     *          /(HL*(H1(IE,JE)+H1(IC,JC)))
          END DO
        End If
      END DO
      RETURN
      END
