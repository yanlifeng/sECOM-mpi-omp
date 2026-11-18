      SUBROUTINE OCLAMP(CLAM)
!     Optimized open boundary conditions for ECOM3D
!     Open boundaries cross the location of the sea surface elevation
!     Version from Shulman et al. (1998)
!     Otimized Boundary Conditions and Data Assimilation with
!     Application to the M2 Tide in the Yellow Sea
!     By Igor Shulman, J.K. Lewis, A.F. Blumberg and N. Kim
!     Journal of Atmospheric and Oceanic Technology
!     *********************************************
! ************************************
!     CLAM is lambda in eqs. (4) -(5)
!     FT* are estimates of energy fluxes through the open boundary (Pt)
!     FOT* are estimates of integrals in numerator of (5)
!     DIFF3* are estimates of integrals in denominator of (5)
! ************************************
      INCLUDE 'comdeck'
      CLAM=1.0
      IF(BCTYPE.EQ.'IRANDB ') GO TO 2000
      CLAM=0.0
      FT=0.0
      DIFF3=0.0
      DO 1000 N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then
! right
          SIGN=1.0
          DEP=(D(IE,JE)+D(IC,JC))/2.
          DEL=(H2(IE,JE)+H2(IC,JC))/2.
          UVM=UA(IC,JC)
!          UVA=UA(IE,JE)+SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=UA(IE,JE)
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then
! left
          SIGN=-1.0
          DEP=(D(IE,JE)+D(IC,JC))/2.
          DEL=(H2(IC,JC)+H2(IE,JE))/2.
          UVM=UA(IC+1,JC)
!          UVA=UA(IC,JC)-SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=UA(IC,JC)
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then
! top
          SIGN=1.0
          DEP=(D(IC,JC)+D(IE,JE))/2.
          DEL=(H1(IE,JE)+H1(IC,JC))/2.
          UVM=VA(IC,JC)
!          UVA=VA(IE,JE)+SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=VA(IE,JE)
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then
! bottom
          SIGN=-1.0
          DEP=(D(IC,JC)+D(IE,JE))/2.
          DEL=(H1(IC,JC)+H1(IE,JE))/2.
          UVM=VA(IC,JC+1)
!          UVA=VA(IC,JC)-SQRT(GRAV/DEP)*(EL(IC,JC)-EL(IE,JE))
          UVA=VA(IC,JC)
        End If
        FT=FT+SIGN*GRAV*DEP*DEL*ELF(IE,JE)*(UVM-UVA)
        DIFF3=DIFF3+DEP*SQRT(GRAV*DEP)*DEL*UVM*UVM
1000  CONTINUE
      GAMMA=DIFF3
      IF(TIME.GT.0.05) CLAM=-FT/(DIFF3+GAMMA)
2000  CONTINUE
!     UPDATING BOUNDARY ETA FOLLOWING Eq (4)
      DO N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        DEPTH=D(IE,JE)
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
! right
          UVA=UA(IE,JE)*DUM(IE,JE)
          SIGN=1.
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then 
! left
          UVA=UA(IC,JC)*DUM(IC,JC)
          SIGN=-1.
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then 
! top
          UVA=VA(IE,JE)*DVM(IE,JE)
          SIGN=1.
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then 
! bottom
          UVA=VA(IC,JC)*DVM(IC,JC)
          SIGN=-1.
        End If
        ELF(IE,JE)=ELF(IE,JE)+SIGN*UVA*SQRT(DEPTH/GRAV)*CLAM
      ENDDO
      RETURN
      END
