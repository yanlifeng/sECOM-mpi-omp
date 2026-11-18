        SUBROUTINE NORTH(WX,WY,dir)
!
!       This program was updated by Ray Sheppard. After a couple of bugs
!       are removed, it is incorporated into ECOMSED (hli, 03/06/01)
! 
!       WX=wind in x-dir, WY=wind in y-dir
!       DIRECTION  GIVES ANGLE WITH X-
!       AXIS (NORTH) POSITIVE- CLOCKWISE, NEGATIVE- ANTI CLOCKWISE
!       this gives the direction of wind blow from (meteorological conv)
!       
        data RAD/57.29577/
        IF (WX .EQ. 0.) THEN
            IF (WY .GE. 0.) THEN
               dir = 180. 
            ELSE
               dir = 0.0
            ENDIF
        ELSE
            GRAD = WY/WX
            dir = ATAN(GRAD) * RAD
            IF (WY .GE. 0.) THEN
                IF (WX .GT. 0.) then
                   dir = 270. - dir
                ELSE
                   dir = abs(dir) + 90.
                ENDIF
            ELSE
                IF (WX .LT. 0.) THEN
                   dir=90-abs(dir)
                ELSE
                   dir=abs(dir)+270
                ENDIF
            ENDIF
        ENDIF
        RETURN
        END
