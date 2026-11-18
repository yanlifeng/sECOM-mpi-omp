      SUBROUTINE RDCORLOC
! This program reads corner_loc a file with the corners of each cell
! stores corners in XCORNER,YCORNER, and cell centroids in
! XGRID,YGRID (by simple division by four for now).
!
!  XCORNER,YCORNER are by definition IM+1,JM+1 ARRAYS.
!  CORNER LOCATION CONVENTION:  XCORNER(I,J) = XGRID(i-1/2, j-1/2) 
!     (LOWER LEFT-HAND CORNER)  YCORNER(I,J) = YGRID(i-1/2, j-1/2)
!
!  NOTE: IF corner_loc is not provided, corners will be created
!        for the model space in the above convention, and centroids
!        will be read in from model_grid.
!
! (Stevens Institute of Technology)
!                                     March 26, 2008

! Modified by A.Jordi on Jan 30, 2015

      INCLUDE 'comdeck'
      LOGICAL IS_CORNER_LOC,IS_MODEL_GRID
      
!NG   Read corners    
      INQUIRE(FILE='corner_loc.nc',EXIST=IS_CORNER_LOC)
      IF (IS_CORNER_LOC) THEN ! Read from corner_loc.nc file
        WRITE(IUPRT,*) '***********************************'
        WRITE(IUPRT,*) 'Reading corners from corner_loc.   '
        WRITE(IUPRT,*) '***********************************'
        CALL NETCDF_READ_CORNER_LOC
!NG     Create simple centroids
        DO J=1,JM
          DO I=1,IM
            IF (FSM(I,J).NE.0) THEN
              XGRID(I,J)=(XCORNER(I,J)+
     .                    XCORNER(I+1,J)+
     .                    XCORNER(I+1,J+1)+
     .                    XCORNER(I,J+1))/4.
              YGRID(I,J)=(YCORNER(I,J)+
     .                    YCORNER(I+1,J)+
     .                    YCORNER(I+1,J+1)+
     .                    YCORNER(I,J+1))/4.
            ELSE
              XGRID(I,J)=-99999. ! NG03312009 PUTCDF FIX
              YGRID(I,J)=-99999. ! TO FLAG LAND VALUES
            ENDIF
          END DO
        END DO
        
        ELSE ! No corner_loc.nc file found
        WRITE(IUPRT,*) '**********************************************'
        WRITE(IUPRT,*) 'No corner_loc file found. Set to model space. '
        WRITE(IUPRT,*) 'For centers, the model_grid file will be used.'
        WRITE(IUPRT,*) '**********************************************'
        DO J=1,JM+1
          DO I=1,IM+1
            XCORNER(I,J)=FLOAT(I)-0.5
            YCORNER(I,J)=FLOAT(J)-0.5
          END DO
        END DO
        DO J=1,JM
          DO I=1,IM
            XGRID(I,J)=0.0
            YGRID(I,J)=0.0
          END DO
        END DO
        INQUIRE(FILE='model_grid.nc',EXIST=IS_MODEL_GRID)
        IF (IS_MODEL_GRID) THEN ! Read XGRID & YGRID from model_grid.nc
          CALL NETCDF_READ_MODEL_GRID('CORNER')
        ELSE ! No model_grid file found
          WRITE(IUPRT,*) '*********************************************'
          WRITE(IUPRT,*) 'NO MODEL_GRID FILE FOUND. PLEASE READ BELOW: '
          WRITE(IUPRT,*) 'Centroids of cells are used to compute tidal '
          WRITE(IUPRT,*) 'constituents and in 2D solar heat flux and   '
          WRITE(IUPRT,*) 'cloud cover computations. If you use any of  '
          WRITE(IUPRT,*) 'these options, the model_grid or corner_loc  '
          WRITE(IUPRT,*) 'files are required for physically correct    '
          WRITE(IUPRT,*) 'results.                                     '
          WRITE(IUPRT,*) '*********************************************'
        END IF
      END IF
      
      RETURN
      END
