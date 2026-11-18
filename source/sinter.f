      SUBROUTINE SINTER(X,A,Y,B,M,N)
      SAVE
!     VERSION(09/27/90)
!
!----------------------------------------------------------------------|
!       A SPECIAL CASE OF INTERP ....NO EXTRAPOLATION BELOW & ABOVE DATA
!       THIS ROUTINE LINEARLY INTERPOLATES AN ARRAY B                  |
!       X(M) MUST BE DESCENDING                                        |
!       A(X) GIVEN FUNCTION                                            |
!       B(Y) FOUND BY LINEAR INTERPOLATION AND EXTRAPOLATION           |
!       Y(N) THE DESIRED DEPTHS                                        |
!       M    THE NUMBER OF POINTS IN X AND A                           |
!       N    THE NUMBER OF POINTS IN Y AND B                           |
!----------------------------------------------------------------------|
!
      DIMENSION X(M),A(M),Y(N),B(N)
!
!-------- BOUNDARIES ---------------------------------------------------
      DO 30 I=1,N
      IF(Y(I).GT.X(1)) B(I)=A(1)
      IF(Y(I).LT.X(M)) B(I)=A(M)
 30   CONTINUE
!
!-------- INTERPOLATION CASES ------------------------------------------
      NM=M-1
      DO 10 I=1,N
      DO 20 J=1,NM
      IF(Y(I).LE.X(J).AND.Y(I).GE.X(J+1))
     +  B(I)=A(J) -(A(J)-A(J+1))*(X(J)-Y(I))/(X(J)-X(J+1))
 20   CONTINUE
 10   CONTINUE
!
      RETURN
      END
