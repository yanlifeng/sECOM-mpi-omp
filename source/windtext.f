      RECURSIVE SUBROUTINE SECOM_WINDTEXT(IUNIT,NX,NY,FM,E,WX,WY,SUN,ICEFR,ICEA)
      USE ISO_FORTRAN_ENV, ONLY: INT32
      IMPLICIT NONE
      INTEGER IUNIT,NX,NY,I,J,K,N,P
      REAL FM(NX,NY),E(NX,NY),WX(NX,NY),WY(NX,NY),SUN(NX,NY),ICEFR(NX,NY),ICEA(NX,NY)
      INTEGER COL(NX)
      INTEGER(INT32) BASE(5),HERE(5)
      REAL VALUES(5)
      LOGICAL SAME
      CHARACTER*50 SUFFIX
      CHARACTER*32 RND,SGN,DEC
! Match the connected unit's edit modes, even on an application that
! opens the text output with nondefault rounding, sign or decimal.
      INQUIRE(UNIT=IUNIT,ROUND=RND,SIGN=SGN,DECIMAL=DEC)
      DO J=1,NY
        N=0
        SAME=STORAGE_SIZE(E(1,J)).EQ.32
        DO I=1,NX
          IF(FM(I,J).NE.0.) THEN
            N=N+1
            COL(N)=I
            IF(SAME) THEN
              VALUES=(/WX(I,J),WY(I,J),SUN(I,J),ICEFR(I,J)*1000.,ICEA(I,J)/)
              HERE=TRANSFER(VALUES,HERE,SIZE=5)
              IF(N.EQ.1) THEN
                BASE=HERE
              ELSE
                SAME=ALL(HERE.EQ.BASE)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
! No WRITE for all-dry rows: an empty implied-do could emit a blank line.
        IF(N.EQ.0) CYCLE
        IF(SAME) THEN
          P=COL(1)
          WRITE(SUFFIX,'(5F10.3)',ROUND=RND,SIGN=SGN,DECIMAL=DEC)
     &       WX(P,J),WY(P,J),SUN(P,J),ICEFR(P,J)*1000.,ICEA(P,J)
! The outer format reverts once per wet grid cell, exactly one record.
          WRITE(IUNIT,'(2I5,F10.3,A50)') (COL(K),J,E(COL(K),J),SUFFIX,K=1,N)
        ELSE
          WRITE(IUNIT,'(2I5,6F10.3)') (COL(K),J,E(COL(K),J),
     &      WX(COL(K),J),WY(COL(K),J),SUN(COL(K),J),ICEFR(COL(K),J)*1000.,ICEA(COL(K),J),K=1,N)
        ENDIF
      ENDDO
      END
