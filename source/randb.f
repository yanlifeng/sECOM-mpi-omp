      SUBROUTINE RANDB
      INCLUDE 'comdeck'
      DO 1000 N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
        If (FSM(IE+1,JE).EQ.0.0.AND.JE.EQ.JC) Then
! right
          DEPTH=(D(IE,JE)+D(IC,JC))/2.
          UAF(IE,JE)=SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        Else If (FSM(IE-1,JE).EQ.0.0.AND.JE.EQ.JC) Then
! left
          DEPTH=(D(IE,JE)+D(IC,JC))/2.
          UAF(IC,JC)=-SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        Else If (FSM(IE,JE+1).EQ.0.0.AND.IE.EQ.IC) Then
! top
          DEPTH=(D(IC,JC)+D(IE,JE))/2.
!         DEPTH=(h(IC,JC)+h(IE,JE))/2.
          VAF(IE,JE)=SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        Else If (FSM(IE,JE-1).EQ.0.0.AND.IE.EQ.IC) Then
! bottom
          DEPTH=(D(IC,JC)+D(IE,JE))/2.
          VAF(IC,JC)=-SQRT(GRAV/DEPTH)*(EL(IC,JC)-EL(IE,JE))
        End If
1000  CONTINUE
      RETURN
      END
