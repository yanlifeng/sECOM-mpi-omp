      Subroutine JDAY(D,M,Y,J)
!     ************************************************************
!**   WHERE...
!        D.... DAY (1-31) COMPONENT OF GREGORIAN DATE
!        M.... MONTH (1-12) COMPONENT
!        Y.... YEAR (E.G., 1979) COMPONENT
!        J.... CORRESPONDING JULIAN DAY NUMBER
!     ************************************************************
      Integer*4 C, J, MO, YR
!     INTEGER*2 D,M,Y
      Integer*4 D, M, Y
!
      YR = Y
      If (M.GT.2) Then
        MO = M - 3
      Else
        MO = M + 9
        YR = YR - 1
      End If
!
      C = YR / 100
      YR = YR - C * 100
      J = (146097*C) / 4 + (1461*YR) / 4 + (153*MO+2) / 5 + D + 1721119
!
!     convert to days since May 23, 1968 
!       (which in True Julian day is 2440000)
      J = J - 2440000
      Return
      End
