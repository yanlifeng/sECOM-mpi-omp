      SUBROUTINE SETDOM_C(INGX,INGY)
!
!     VERSION(04/23/2013)
!
      INCLUDE 'comdeck'
!
!  TAKE OUT FOUND BOUNDARY CORNER AS IN SETDOM
!        
!
      SAVE

      FSM(INGX,INGY)=0.0
      DUM(INGX,INGY)=0.0
      DVM(INGX,INGY)=0.0
      DUM(INGX+1,INGY)=0.0
      DVM(INGX,INGY+1)=0.0
      DUMWAD(INGX,INGY)=0.0
      DVMWAD(INGX,INGY)=0.0
      DUMWAD(INGX+1,INGY)=0.0
      DVMWAD(INGX,INGY+1)=0.0
      COR(INGX,INGY)=0.0
! 
      IF(TOR.NE.'BAROTROPIC') THEN
      DO 350 K=1,KBM1
      AAM(INGX,INGY,K)=0.0
  350 CONTINUE
      AAM2D(INGX,INGY)=0.0
      ELSE
      AAM2D(INGX,INGY)=0.0
      ENDIF
!
      DO 370 K=2,KBM1
      L(INGX,INGY,K)   =0.0
      Q2B(INGX,INGY,K) =0.0
      Q2(INGX,INGY,K)  =0.0
      Q2LB(INGX,INGY,K)=0.0
      Q2L(INGX,INGY,K) =0.0
      KM(INGX,INGY,K)  =0.0
      KQ(INGX,INGY,K)  =0.0
      KH(INGX,INGY,K)  =0.0
  370 CONTINUE
!
!*********************************************************************
!
!
      RETURN
      END 
