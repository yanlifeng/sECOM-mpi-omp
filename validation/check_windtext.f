      PROGRAM TEST_WIND
      IMPLICIT NONE
      INTEGER,PARAMETER::NX=29,NY=11
      INTEGER I,J,M,IU,JU
      REAL FM(NX,NY),E(NX,NY),X(NX,NY),Y(NX,NY),S(NX,NY),F(NX,NY),A(NX,NY)
      CHARACTER*24 ANAME,BNAME
      CHARACTER*32 RND(4),SIG(4),DEC(4)
      DATA RND /'NEAREST','UP','DOWN','ZERO'/
      DATA SIG /'PROCESSOR_DEFINED','PLUS','SUPPRESS','PLUS'/
      DATA DEC /'POINT','POINT','COMMA','COMMA'/
      DO J=1,NY
        DO I=1,NX
          FM(I,J)=1.
          IF(MOD(I+J,4).EQ.2.OR.J.EQ.3) FM(I,J)=0.
          E(I,J)=.0625*(I-J)
          X(I,J)=.0625
          Y(I,J)=1.17
          S(I,J)=13.
          F(I,J)=0.
          A(I,J)=0.
          IF(J.EQ.4.AND.MOD(I,3).EQ.1) X(I,J)=-0.
          IF(J.EQ.4.AND.MOD(I,3).NE.1) X(I,J)=0.
          IF(J.EQ.5) Y(I,J)=1.E25
          IF(J.EQ.6) S(I,J)=-121.678432
          IF(J.EQ.7.AND.I.EQ.NX) F(I,J)=.0174
          IF(J.EQ.8) A(I,J)=I*.0625
          IF(J.EQ.9.AND.I.NE.20) FM(I,J)=0.
        ENDDO
      ENDDO
      DO M=1,4
        WRITE(ANAME,'("reference-",I1,".dat")') M
        WRITE(BNAME,'("candidate-",I1,".dat")') M
        OPEN(NEWUNIT=IU,FILE=ANAME,STATUS='NEW',ROUND=RND(M),SIGN=SIG(M),DECIMAL=DEC(M))
        OPEN(NEWUNIT=JU,FILE=BNAME,STATUS='NEW',ROUND=RND(M),SIGN=SIG(M),DECIMAL=DEC(M))
        DO J=1,NY
          DO I=1,NX
            IF(FM(I,J).NE.0.) WRITE(IU,'(2I5,6F10.3)') I,J,E(I,J),X(I,J),Y(I,J),S(I,J),F(I,J)*1000.,A(I,J)
          ENDDO
        ENDDO
        CALL SECOM_WINDTEXT(JU,NX,NY,FM,E,X,Y,S,F,A)
        CLOSE(IU)
        CLOSE(JU)
      ENDDO
      END
