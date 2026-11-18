! ------------------------------------------------------------------
      SUBROUTINE SMOOTH(A,FSM,IM,JM,NITS)
!-------------------------------------------------------------------
!     THIS ROUTINE SMOOTHS DATA WITH A FIVE POINT LAPLACIAN FILTER.
!-------------------------------------------------------------------
      DIMENSION A(IM,JM),FSM(IM,JM)
      DO 100 N=1,NITS
!      DO 200 J=2,JM-1
!      DO 200 I=2,IM-1
      DO 200 J=JM-1,2,-1
      DO 200 I=IM-1,2,-1
!
      SMFAC=FSM(I+1,J)+FSM(I,J-1)+FSM(I-1,J)+FSM(I,J+1)+1.E-10
      A(I,J)=A(I,J)+(.5/SMFAC)
     1               *(A(I+1,J)*FSM(I+1,J)+A(I,J-1)*FSM(I,J-1)
     2                +A(I-1,J)*FSM(I-1,J)+A(I,J+1)*FSM(I,J+1)
     3                -SMFAC*A(I,J))
      A(I,J)=A(I,J)*FSM(I,J)
  200 CONTINUE
  100 CONTINUE
      RETURN
      END
