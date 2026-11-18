      SUBROUTINE PCLAMP
!     VERSION(12/31/98)
!NG     VERSION(10/05/11) WITH GRAV, AVERAGED H + EL, CORRECTED ORLANSKI AND SIMPLIFIED
!
      INCLUDE 'comdeck'
!
!  MODIFIED BY Quamrul QA 12/31/98 
!
        DO 100 N = 1, NUMEBC
        IE = IETA(N)
        JE = JETA(N)
        IC = ICON(N)
        JC = JCON(N)
!
! --- NORTH OR SOUTH BOUNDARY (RADIATION & FORCING) --------------------
!
      IF (IE.EQ.IC) Then
        CPH=SQRT(GRAV*0.5*(H(IE,JE)+EL(IE,JE)+H(IC,JC)+EL(IC,JC)))*2.
     .      *DTE/(H2(IE,JE)+H2(IC,JC))
!
! --- EAST OR WEST BOUNDARY (RADIATION & FORCING) --------------------
!
      ELSE
        CPH=SQRT(GRAV*0.5*(H(IE,JE)+EL(IE,JE)+H(IC,JC)+EL(IC,JC)))*2.
     .      *DTE/(H1(IE,JE)+H1(IC,JC))
      ENDIF
!NG-ORLANSKI-IMPLICIT      ELF(IE,JE)=(ELB(IE,JE)+CPH*(2.*EL(IC,JC)-ELB(IC,JC))+2.*DTE*
      ELF(IE,JE)=(ELB(IE,JE)+CPH*(2.*EL(IC,JC)-ELB(IE,JE))+2.*DTE*
     .           (ELF(IE,JE)-ELB(IE,JE))*TLAG)/(1.+CPH)
 100  CONTINUE
!
      RETURN
      END
