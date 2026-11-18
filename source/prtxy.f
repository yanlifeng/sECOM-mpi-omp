      SUBROUTINE PRTXY(FD,FMASK,IM,JM,KB,KT,IARRAY,SCALE,IUNIT,DEV)
!     VERSION(05/24/91)
!
      SAVE
      CHARACTER*3 DEV
      CHARACTER*7 LEVEL(2)
      DIMENSION FD(IM,JM,KB),IARRAY(IM,JM),FMASK(IM,JM)
      DIMENSION KP(2)
      KP(1)=KT
      KP(2)=KB-1
      LEVEL(1)='SURFACE'
      LEVEL(2)='BOTTOM '
!
!-----------------------------------------------------------------------
!        THIS SUBROUTINE WRITES HORIZONTAL LAYERS OF A 3-D FIELD
!-----------------------------------------------------------------------
!
      DO 10 KM=1,2
      K=KP(KM)
      WRITE(IUNIT,20) LEVEL(KM)
 20   FORMAT(//A7,' LAYER')
      CALL PRINT(FD(1,1,K),FMASK,IM,JM,IARRAY,SCALE,IUNIT,DEV)
 10   CONTINUE
!
      RETURN
      END
