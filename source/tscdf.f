      Subroutine TSCDF(TMIDDLE,icdfcheck)
!
! write EPIC compliant netCDF time series
!
!NG03312009 Brought to compliance with the latest
!NG03312009 NetCDF Climate and Forecast (CF)
!NG03312009 MetaData Conventions Version 1.4, 27 February, 2009
!NG03312009
!NG03132009 http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/cf-conventions.pdf
!NG02192010
!NG02192010 Included independent definitions for 2D (elevation) versus 3D (velocity)
!NG02192010 stations.
!NG03182010 Fixed the issue of not writting the first requested timestamp (unless it's the
!NG03182010 first cold start).
      Include 'netcdf.inc'
      Include 'comdeck'
      Integer IRET, BASE_DATE(4)
      Character*60 CDFFILE
      CHARACTER*1 PTC
!NG03312009
      CHARACTER modelname*80
      CHARACTER NOWDATE*9,BASEDATE*26
* netCDF id
      Integer CDFID
* dimension ids
!NG02192010      Integer STADIM, ZPOSDIM, TIMEDIM, DIMDIM, DUM1DIM, DUM2DIM, 
!NG02192010     *    DUM3DIM
      Integer STA2DDIM, STA3DDIM, ZPOSDIM, TIMEDIM,  ! NG02192010 POSSIBLY DIFFERENT FOR EPTS,VPTS (NOW 2D, 3D)
     *    DIMDIM, DUM1DIM, DUM2DIM, DUM3DIM
     *   ,NFLUXDIM ! NG04022008 store fluxes too if needed
* variable ids
      Integer LOC2DID, LOC3DID, SIGMAID, DEPTH2DID, DEPTH3DID,TID, ELID
     *   ,UID, VID ,SALID, TMPID, FLXID, STA2DID, STA3DID,DUM1ID
     *   ,DUM2ID, DUM3ID, TXID, TYID, DIMID, WID, KHID, KMID, AMID, HFID
     *   ,CID,X2DID,Y2DID,X3DID,Y3DID,PATMID,QPRECID,QEVAPID
     *   ,FLUXID ! NG04022008 store fluxes too
     *   ,AIRTID,CLDID,RHUMID,SWOBSID ! NG04092008 store 2d met too
     *   ,WVPDID,WVHTID,WVDRID ! NG04142008 store waves too
     *   ,UWNDID,VWNDID, DEPTHTID ! NG03312009 for wind and instant depth
* variable shapes
      Integer DIMS(4), NDUM(EPTSM), NDUM2(2)
* corners and edge lengths
      Integer CORNER(4), EDGES(4)
* attribute vectors
      Real FLOATVAL(1)
!NG      Parameter (NWORK = EPTSM*KBM1)
      Parameter (NWORK = EPTSM*KB)      ! NG 11/05/02
      Real WA1(NWORK)
!  temporary declaration of variables (NKIM 1.12.00) clean up later
      REAL CZSAVE(VPTSM,KB)    ! NG 11/07/02
      REAL XCENTROID(EPTSM),YCENTROID(EPTSM),BATHY(EPTSM)     ! NG 08242007
      
* count the number of times this routine is called
      Data ICDF /1/ ! NG03182010
      Save
      If (ICDF.EQ.1) Then ! NG03182010
        BASE_DATE(1) = IYR ! CNG04172008
        BASE_DATE(2) = IMO ! CNG04172008
        BASE_DATE(3) = IDA ! CNG04172008
        BASE_DATE(4) = IHR ! CNG04172008
        BASEDATE="YYYY-0M-0D 0H:00:00 -00:00"
        write (BASEDATE(1:4),'(i4)') IYR
        if (IMO.GT.9) THEN
         write (BASEDATE(6:7),'(i2)') IMO
        else
         write (BASEDATE(7:7),'(i1)') IMO
        endif
        if (IDA.GT.9) THEN
         write (BASEDATE(9:10),'(i2)') IDA
        else
         write (BASEDATE(10:10),'(i1)') IDA
        endif
        if (IHR.GT.9) THEN
         write (BASEDATE(12:13),'(i2)') IHR
        else
         write (BASEDATE(13:13),'(i1)') IHR
        endif
        ISNG1=IFIX(UTCSHIFT)
        ISNG2=IFIX(ABS((UTCSHIFT-AINT(UTCSHIFT))*60.))
        if (ISNG1.GT.9) THEN
         write (BASEDATE(21:21),'(a1)') '+'
         write (BASEDATE(22:23),'(i2)') ISNG1
        elseif (ISNG1.GE.0) THEN
         write (BASEDATE(21:21),'(a1)') '+'
         write (BASEDATE(23:23),'(i1)') ISNG1
        elseif (ISNG1.GT.-10) THEN
         write (BASEDATE(23:23),'(i1)') IABS(ISNG1)
        else
         write (BASEDATE(22:23),'(i2)') IABS(ISNG1)
        endif
        if (ISGN2.GT.9) THEN
         write (BASEDATE(25:26),'(i2)') ISNG2
        else
         write (BASEDATE(26:26),'(i1)') ISNG2
        endif
!        write (*,*) 'Basedate',BASEDATE
!NG        CDFFILE = 'tsepic.nc' ! NG03312009 changed name to .nc as per CF
        CDFFILE = 'tsepic.cdf' ! For NYHOPS, too much of a pain otherwise.
* HARDWIRE heat flux save
!NG HARDWIRE heat flux and the rest as in putcdf.f
        PTHF = 'n' ! Not implemented yet
        PTQ2 = 'n' ! Not implemented yet
* REMOVED HARDWIRE: don't save mixing parameters
!NG        PTKH = 'y'
!NG        PTKM = 'y'
!NG        PTAM = 'y'
!       PTU = 'y'
!       PTV = 'y'
!       PTW = 'y'
!       PTS = 'y'
!       PTT = 'y'
!NG End of putcdf.f-like HARDWIRE ! NG 11/08/02
* Don't need T, S, W, if run is barotropic (2-D) (Overwrite ! NG 11/08/02)
        If (TOR.EQ.'BAROTROPIC') Then ! Moved here as in putcdf.f ! NG 11/08/02
          PTS = 'n'
          PTT = 'n'
          PTS = 'n'
          PTW = 'n'
          PTKH = 'n'
          PTKM = 'n'
        End If
* FOR CONSERVATIVE TRACER
       IF(TRACER.EQ.'INCLUDE')THEN
        PTC = 'y'
       ELSE
        PTC = 'n'
       END IF

* enter define mode
        CDFID = NCCRE(CDFFILE,NCCLOB,IRET)
        IFILL = NCSFIL(CDFID,NCNOFILL,IRET)
* define dimensions
!NG02192010 Included cases when EPTS or VPTS might be zero, and VPTS i
        If (EPTS.GT.0) Then
         STA2DDIM = NCDDEF(CDFID,'2Dstation',EPTS,IRET)
        Endif
        If (VPTS.GT.0) Then
         STA3DDIM = NCDDEF(CDFID,'3Dstation',VPTS,IRET)
        Endif
        If (FPTS.GT.0) Then
         NFLUXDIM = NCDDEF(CDFID,'transects',FPTS,IRET) ! NG 04022008
        Endif
!NG        ZPOSDIM = NCDDEF(CDFID,'sigma',KBM1,IRET)
        ZPOSDIM = NCDDEF(CDFID,'sigma',KB,IRET) ! NG 11/05/02
        TIMEDIM = NCDDEF(CDFID,'time',NCUNLIM,IRET)
        DIMDIM = NCDDEF(CDFID,'dim',2,IRET)
        DUM1DIM = NCDDEF(CDFID,'dum1',1,IRET)
        DUM2DIM = NCDDEF(CDFID,'dum2',1,IRET)
        DUM3DIM = NCDDEF(CDFID,'dum3',1,IRET)
* define variables
* 1D Vars 
* time
        DIMS(1) = TIMEDIM
        TID = NCVDEF(CDFID,'time',NCLONG,1,DIMS,IRET)
* z (sigma levels)
        DIMS(1) = ZPOSDIM
        SIGMAID = NCVDEF(CDFID,'sigma',NCFLOAT,1,DIMS,IRET)
        Call NCAPT(CDFID,SIGMAID,'epic_code',NCLONG,1,1,IRET)
* 2D Stations
        IF(EPTS.GT.0) THEN ! NG02192010
        DIMS(1) = STA2DDIM
* x (stations)
        STA2DID = NCVDEF(CDFID,'2Dstation',NCLONG,1,DIMS,IRET)
        Call NCAPTC(CDFID,STA2DID,'long_name',NCCHAR,13,'2D station id',
     *       IRET)
        Call NCAPT(CDFID,STA2DID,'epic_code',NCLONG,1,501,IRET)
        X2DID = NCVDEF(CDFID,'2Dlon',NCFLOAT,1,DIMS,IRET)
        Call NCAPTC(CDFID,X2DID,'long_name',NCCHAR,9,'longitude',IRET)
        Call NCAPTC(CDFID,X2DID,'units',NCCHAR,12,'degrees_east',IRET)
        Call NCAPTC(CDFID,X2DID,'standard_name',NCCHAR,9,
     *              'longitude',IRET)
        Y2DID = NCVDEF(CDFID,'2Dlat',NCFLOAT,1,DIMS,IRET)
        Call NCAPTC(CDFID,Y2DID,'long_name',NCCHAR,8,'latitude',IRET)
        Call NCAPTC(CDFID,Y2DID,'units',NCCHAR,13,'degrees_north',IRET)
        Call NCAPTC(CDFID,Y2DID,'standard_name',NCCHAR,8,
     *              'latitude',IRET)
        DEPTH2DID = NCVDEF(CDFID,'2Ddepth',NCFLOAT,1,DIMS,IRET)
        Call NCAPTC(CDFID,DEPTH2DID,'long_name',NCCHAR,10,'Bathymetry',
     *      IRET)
        Call NCAPTC(CDFID,DEPTH2DID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,DEPTH2DID,'standard_name',NCCHAR,5,
     *              'depth',IRET)
        Call NCAPTC(CDFID,DEPTH2DID,'positive',NCCHAR,4,
     *              'down',IRET)
        Call NCAPTC(CDFID,DEPTH2DID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ENDIF ! NG02192010
* 3D Stations
        IF(VPTS.GT.0) THEN ! NG02192010
        DIMS(1) = STA3DDIM
* x (stations)
        STA3DID = NCVDEF(CDFID,'3Dstation',NCLONG,1,DIMS,IRET)
        Call NCAPTC(CDFID,STA3DID,'long_name',NCCHAR,13,'3D station id',
     *       IRET)
        Call NCAPT(CDFID,STA3DID,'epic_code',NCLONG,1,501,IRET)
        X3DID = NCVDEF(CDFID,'3Dlon',NCFLOAT,1,DIMS,IRET)
        Call NCAPTC(CDFID,X3DID,'long_name',NCCHAR,9,'longitude',IRET)
        Call NCAPTC(CDFID,X3DID,'units',NCCHAR,12,'degrees_east',IRET)
        Call NCAPTC(CDFID,X3DID,'standard_name',NCCHAR,9,
     *              'longitude',IRET)
        Y3DID = NCVDEF(CDFID,'3Dlat',NCFLOAT,1,DIMS,IRET)
        Call NCAPTC(CDFID,Y3DID,'long_name',NCCHAR,8,'latitude',IRET)
        Call NCAPTC(CDFID,Y3DID,'units',NCCHAR,13,'degrees_north',IRET)
        Call NCAPTC(CDFID,Y3DID,'standard_name',NCCHAR,8,
     *              'latitude',IRET)
        DEPTH3DID = NCVDEF(CDFID,'3Ddepth',NCFLOAT,1,DIMS,IRET)
        Call NCAPTC(CDFID,DEPTH3DID,'long_name',NCCHAR,10,'Bathymetry',
     *      IRET)
        Call NCAPTC(CDFID,DEPTH3DID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,DEPTH3DID,'standard_name',NCCHAR,5,
     *              'depth',IRET)
        Call NCAPTC(CDFID,DEPTH3DID,'positive',NCCHAR,4,
     *              'down',IRET)
        Call NCAPTC(CDFID,DEPTH3DID,'coordinates',NCCHAR,11,
     *              '3Dlat 3Dlon',IRET)
        ENDIF ! NG02192010
! NG08242007
* psuedo x for wind 
        DUM3ID = NCVDEF(CDFID,'dum3',NCLONG,1,DUM3DIM,IRET)
        Call NCAPTC(CDFID,DUM3ID,'long_name',NCCHAR,6,'dummy3',IRET)
        Call NCAPT(CDFID,DUM3ID,'epic_code',NCLONG,1,501,IRET)
* psuedo y for all variables
        DUM1ID = NCVDEF(CDFID,'dum1',NCLONG,1,DUM1DIM,IRET)
        Call NCAPTC(CDFID,DUM1ID,'long_name',NCCHAR,6,'dummy1',IRET)
        Call NCAPT(CDFID,DUM1ID,'epic_code',NCLONG,1,500,IRET)
* psuedo z for elevation
        DUM2ID = NCVDEF(CDFID,'dum2',NCLONG,1,DUM2DIM,IRET)
        Call NCAPTC(CDFID,DUM2ID,'long_name',NCCHAR,6,'dummy2',IRET)
        Call NCAPT(CDFID,DUM2ID,'epic_code',NCLONG,1,1,IRET)
* stations
        DIMID = NCVDEF(CDFID,'dim',NCLONG,1,DIMDIM,IRET)
        Call NCAPTC(CDFID,DIMID,'long_name',NCCHAR,6,'dummy4',IRET)
*   2D Vars
        IF(EPTS.GT.0) THEN ! NG02192010
        DIMS(1) = STA2DDIM
        DIMS(2) = DUM1DIM
        DIMS(3) = DUM2DIM
        DIMS(4) = TIMEDIM
!NG04092008 Unfortunatelly, for legacy purposes, we store stress here as wu,wv.
!NG04092008 It should really have been wtaux, wtauy. One could change for
!NG04092008 a new model with new graphics. wu, wv should have been wind from
!NG04092008 WSX/Y2DSAVE.
!NG03312009 FIXED AND ADDED WINDS
        TXID = NCVDEF(CDFID,'wtaux',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,TXID,'long_name',NCCHAR,15,'E-W Wind Stress',
     *      IRET)
        Call NCAPTC(CDFID,TXID,'units',NCCHAR,6,'pascal',IRET)
        Call NCAPTC(CDFID,TXID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        TYID = NCVDEF(CDFID,'wtauy',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,TYID,'long_name',NCCHAR,15,'N-S Wind Stress',
     *      IRET)
        Call NCAPTC(CDFID,TYID,'units',NCCHAR,6,'pascal',IRET)
        Call NCAPTC(CDFID,TYID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)

        UWNDID = NCVDEF(CDFID,'wu',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,UWNDID,'long_name',NCCHAR,8,'E-W Wind',
     *      IRET)
        Call NCAPTC(CDFID,UWNDID,'units',NCCHAR,5,'m s-1',IRET)
        Call NCAPTC(CDFID,UWNDID,'standard_name',NCCHAR,13,
     *            'eastward_wind',IRET)
        Call NCAPTC(CDFID,UWNDID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        VWNDID = NCVDEF(CDFID,'wv',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,VWNDID,'long_name',NCCHAR,8,'N-S Wind',
     *      IRET)
        Call NCAPTC(CDFID,VWNDID,'units',NCCHAR,5,'m s-1',IRET)
        Call NCAPTC(CDFID,VWNDID,'standard_name',NCCHAR,14,
     *            'northward_wind',IRET)
        Call NCAPTC(CDFID,VWNDID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
!NG03312009
        ELID = NCVDEF(CDFID,'elev',NCFLOAT,4,DIMS,IRET)
        ! NG01172008 Included Atmospheric Pressure
        PATMID = NCVDEF(CDFID,'patm',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,PATMID,'long_name',NCCHAR,20,
     *     'Atmospheric Pressure',IRET)
        Call NCAPTC(CDFID,PATMID,'units',NCCHAR,4,'mbar',IRET)
        Call NCAPTC(CDFID,PATMID,'standard_name',NCCHAR,25,
     *             'air_pressure_at_sea_level',IRET)
        Call NCAPTC(CDFID,PATMID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! NG09212011 Included Precipitation
        QPRECID = NCVDEF(CDFID,'qprec',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,QPRECID,'long_name',NCCHAR,13,
     *     'Rainfall Rate',IRET)
        Call NCAPTC(CDFID,QPRECID,'units',NCCHAR,7,'mm hr-1',IRET)
        Call NCAPTC(CDFID,QPRECID,'standard_name',NCCHAR,13,
     *             'rainfall_rate',IRET)
        Call NCAPTC(CDFID,QPRECID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! NG09212011 Included Evaporation
        QEVAPID = NCVDEF(CDFID,'qevap',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,QEVAPID,'long_name',NCCHAR,16,
     *     'Evaporation Rate',IRET)
        Call NCAPTC(CDFID,QEVAPID,'units',NCCHAR,6,'m yr-1',IRET)
        Call NCAPTC(CDFID,QEVAPID,'standard_name',NCCHAR,26,
     *             'lwe_water_evaporation_rate',IRET)
        Call NCAPTC(CDFID,QEVAPID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! NG04092008 Included other 2d met variables
        ! Air Temperature
        AIRTID = NCVDEF(CDFID,'airt',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,AIRTID,'long_name',NCCHAR,32,
     *     '2m-above-surface Air Temperature',IRET)
        Call NCAPTC(CDFID,AIRTID,'units',NCCHAR,7,'Celsius',IRET)
        Call NCAPTC(CDFID,AIRTID,'standard_name',NCCHAR,19,
     *              'surface_temperature',IRET)
        Call NCAPTC(CDFID,AIRTID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! Cloud Cover
        CLDID = NCVDEF(CDFID,'cld',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,CLDID,'long_name',NCCHAR,20,
     *     'Cloud Cover fraction',IRET)
        Call NCAPTC(CDFID,CLDID,'units',NCCHAR,1,'1',IRET)
        Call NCAPTC(CDFID,CLDID,'standard_name',NCCHAR,19,
     *             'cloud_area_fraction',IRET)
        Call NCAPTC(CDFID,CLDID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! Relative Humidity
        RHUMID = NCVDEF(CDFID,'rhum',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,RHUMID,'long_name',NCCHAR,34,
     *      '2m-above-surface Relative Humidity',IRET)
        Call NCAPTC(CDFID,RHUMID,'units',NCCHAR,7,'percent',IRET)
        Call NCAPTC(CDFID,RHUMID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! Surface-Penetrating Shortwave Radiation
        SWOBSID = NCVDEF(CDFID,'swobs',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,SWOBSID,'long_name',NCCHAR,27,
     *     'Shortwave Radiation in W/m2',IRET)
        Call NCAPTC(CDFID,SWOBSID,'units',NCCHAR,5,'W m-2',IRET)
        Call NCAPTC(CDFID,SWOBSID,'standard_name',NCCHAR,41,
     *              'surface_downwelling_shortwave_flux_in_air',IRET)
        Call NCAPTC(CDFID,SWOBSID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! Wave Period
        WVPDID = NCVDEF(CDFID,'wp',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,WVPDID,'long_name',NCCHAR,11,
     *     'Wave Period',IRET)
        Call NCAPTC(CDFID,WVPDID,'units',NCCHAR,3,'sec',IRET)
        Call NCAPTC(CDFID,WVPDID,'standard_name',NCCHAR,28,
     *             'sea_surface_wind_wave_period',IRET)
        Call NCAPTC(CDFID,WVPDID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! Significant Wave Height
        WVHTID = NCVDEF(CDFID,'wh',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,WVHTID,'long_name',NCCHAR,23,
     *     'Significant Wave Height',IRET)
        Call NCAPTC(CDFID,WVHTID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,WVHTID,'standard_name',NCCHAR,40,
     *             'sea_surface_wind_wave_significant_height',IRET)
        Call NCAPTC(CDFID,WVHTID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        ! Significant Wave Direction
        WVDRID = NCVDEF(CDFID,'wd',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,WVDRID,'long_name',NCCHAR,54,
     *    'Significant Wave Direction in Oceanographic Convention',IRET)
        Call NCAPTC(CDFID,WVDRID,'units',NCCHAR,13,'degrees_north',IRET)
        Call NCAPTC(CDFID,WVDRID,'standard_name',NCCHAR,34,
     *             'sea_surface_wind_wave_to_direction',IRET)
        Call NCAPTC(CDFID,WVDRID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)

        If (PTHF.EQ.'Y'.OR.PTHF.EQ.'y') Then
          HFID = NCVDEF(CDFID,'heat_flux',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,HFID,'long_name',NCCHAR,17,
     *         'heat flux in W/m2',IRET)
          Call NCAPTC(CDFID,HFID,'units',NCCHAR,5,'W m-2',IRET)
          Call NCAPTC(CDFID,HFID,'standard_name',NCCHAR,39,
     *              'surface_downward_heat_flux_in_sea_water',IRET)
          Call NCAPTC(CDFID,HFID,'coordinates',NCCHAR,11,
     *                '2Dlat 2Dlon',IRET)
        endif
        ENDIF ! NG02192010

*   1D, 2D, and 3D attributes
        Call NCAPTC(CDFID,SIGMAID,'long_name',NCCHAR,36,
     *      'stretched vertical coordinate levels',IRET)
!NGCF CAREFUL IF COARDS        Call NCAPTC(CDFID,SIGMAID,'units',NCCHAR,11,'sigma_level',
!NGCF CAREFUL IF COARDS     *      IRET)
        Call NCAPTC(CDFID,SIGMAID,'units',NCCHAR,1,'1',
     *      IRET)
        Call NCAPTC(CDFID,SIGMAID,'positive',NCCHAR,2,
     *              'up',IRET)
        Call NCAPTC(CDFID,SIGMAID,'standard_name',NCCHAR,22,
     *              'ocean_sigma_coordinate',IRET)
        Call NCAPTC(CDFID,SIGMAID,'formula_terms',NCCHAR,35,
     *              'sigma: sigma eta: elev depth: depth',IRET)
!NG03312009 The following are IDV specific (?)
        Call NCAPTC(CDFID,SIGMAID,'_CoordinateAxisType',NCCHAR,4,
     *              'GeoZ',IRET)
        Call NCAPTC(CDFID,SIGMAID,'_CoordinateZisPositive',NCCHAR,2,
     *              'up',IRET)
        Call NCAPTC(CDFID,SIGMAID,'_CoordinateTransformType',NCCHAR,8,
     *              'Vertical',IRET)
        Call NCAPTC(CDFID,SIGMAID,'_CoordinateAxes',NCCHAR,5,
     *              'sigma',IRET)

        Call NCAPTC(CDFID,TID,'long_name',NCCHAR,4,'time',IRET)
        Call NCAPT(CDFID,TID,'epic_code',NCLONG,1,625,IRET)
        Call NCAPTC(CDFID,TID,'units',NCCHAR,40,
     *              'seconds since '//BASEDATE,IRET)
        Call NCAPTC(CDFID,TID,'standard_name',NCCHAR,4,
     *              'time',IRET)
        Call NCAPTC(CDFID,ELID,'long_name',NCCHAR,9,'elevation',IRET)
        Call NCAPTC(CDFID,ELID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,ELID,'standard_name',NCCHAR,34,
     *              'sea_surface_height_above_sea_level',IRET)
        Call NCAPTC(CDFID,ELID,'coordinates',NCCHAR,11,
     *              '2Dlat 2Dlon',IRET)
        Call NCAPTC(CDFID,ELID,'positive',NCCHAR,2,
     *              'up',IRET)
        ! NG04262007 moved DEPTHT to a 2D variable as it should be
        ! (=H+ET)

*   3D Vars
        IF(VPTS.GT.0) THEN ! NG02192010
        DIMS(1) = STA3DDIM
        DIMS(2) = DUM1DIM
        DIMS(3) = DUM2DIM
        DIMS(4) = TIMEDIM
        ! NG04262007 moved DEPTHT to a 2D variable as it should be
        ! (=H+ET)
        DEPTHTID = NCVDEF(CDFID,'depth_t',NCFLOAT,4,DIMS,IRET)
        Call NCAPTC(CDFID,DEPTHTID,'long_name',NCCHAR,11,'water_depth',
     *      IRET)
        Call NCAPTC(CDFID,DEPTHTID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,DEPTHTID,'coordinates',NCCHAR,11,
     *              '3Dlat 3Dlon',IRET) ! NG02192010 But on the 3d lat-long

        DIMS(1) = STA3DDIM
        DIMS(2) = DUM1DIM
        DIMS(3) = ZPOSDIM
        DIMS(4) = TIMEDIM
!
        If (PTU.EQ.'Y'.OR.PTU.EQ.'y') Then
!   u velocity
          UID = NCVDEF(CDFID,'u',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,UID,'long_name',NCCHAR,10,'U velocity',IRET)
          Call NCAPTC(CDFID,UID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,UID,'standard_name',NCCHAR,27,
     *              'eastward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,UID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
!
        If (PTV.EQ.'Y'.OR.PTV.EQ.'y') Then
!   v velocity
          VID = NCVDEF(CDFID,'v',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,VID,'long_name',NCCHAR,10,'V velocity',IRET)
          Call NCAPTC(CDFID,VID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,VID,'standard_name',NCCHAR,28,
     *              'northward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,VID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
!
        If (PTW.EQ.'Y'.OR.PTW.EQ.'y') Then
!   w velocity
          WID = NCVDEF(CDFID,'w',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,WID,'long_name',NCCHAR,10,'W velocity',IRET)
          Call NCAPTC(CDFID,WID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,WID,'standard_name',NCCHAR,25,
     *              'upward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,WID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        If (PTC.EQ.'Y'.OR.PTC.EQ.'y') Then
!   nonconservative tracer
          CID = NCVDEF(CDFID,'conc',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,CID,'long_name',NCCHAR,13,'Concentration',
     *        IRET)
          Call NCAPTC(CDFID,CID,'units',NCCHAR,4,'none',IRET)
          Call NCAPTC(CDFID,CID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        If (PTS.EQ.'Y'.OR.PTS.EQ.'y') Then
!   salinity
          SALID = NCVDEF(CDFID,'salt',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,SALID,'long_name',NCCHAR,15,
     *         'Salinity in psu',IRET)
          Call NCAPTC(CDFID,SALID,'units',NCCHAR,4,'1e-3',IRET)
          Call NCAPTC(CDFID,SALID,'standard_name',NCCHAR,18,
     *              'sea_water_salinity',IRET)
          Call NCAPTC(CDFID,SALID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        If (PTKH.EQ.'Y'.OR.PTKH.EQ.'y') Then
!   vertical diffusivity
          KHID = NCVDEF(CDFID,'kh',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,KHID,'long_name',NCCHAR,20,
     *        'Vertical Diffusivity',IRET)
          Call NCAPTC(CDFID,KHID,'units',NCCHAR,6,'m2 s-1',IRET)
          Call NCAPTC(CDFID,KHID,'standard_name',NCCHAR,31,
     *              'ocean_vertical_heat_diffusivity',IRET)
          Call NCAPTC(CDFID,KHID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        If (PTKM.EQ.'Y'.OR.PTKM.EQ.'y') Then
!   vertical viscosity
          KMID = NCVDEF(CDFID,'km',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,KMID,'long_name',NCCHAR,18,
     *        'Vertical Viscosity',IRET)
          Call NCAPTC(CDFID,KMID,'units',NCCHAR,6,'m2 s-1',IRET)
          Call NCAPTC(CDFID,KMID,'standard_name',NCCHAR,35,
     *              'ocean_vertical_momentum_diffusivity',IRET)
          Call NCAPTC(CDFID,KMID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        If (PTAM.EQ.'Y'.OR.PTAM.EQ.'y') Then
!   horizontal viscosity
          AMID = NCVDEF(CDFID,'am',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,AMID,'long_name',NCCHAR,20,
     *        'Horizontal Viscosity',IRET)
          Call NCAPTC(CDFID,AMID,'units',NCCHAR,6,'m2 s-1',IRET)
          Call NCAPTC(CDFID,AMID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        If (PTT.EQ.'Y'.OR.PTT.EQ.'y') Then
!   temperature
          TMPID = NCVDEF(CDFID,'temp',NCFLOAT,4,DIMS,IRET)
          Call NCAPTC(CDFID,TMPID,'long_name',NCCHAR,17,
     *        'Water Temperature',IRET)
          Call NCAPTC(CDFID,TMPID,'units',NCCHAR,7,'Celsius',IRET)
          Call NCAPTC(CDFID,TMPID,'standard_name',NCCHAR,21,
     *              'sea_water_temperature',IRET)
          Call NCAPTC(CDFID,TMPID,'coordinates',NCCHAR,17,
     *                '3Dlat 3Dlon sigma',IRET)
        End If
        ENDIF ! NG02192010

*   Fluxes NG04022008
        If (FPTS.GT.0) Then
         DIMS(1) = NFLUXDIM
         DIMS(2) = DUM1DIM
         DIMS(3) = ZPOSDIM
         DIMS(4) = TIMEDIM
         FLUXID = NCVDEF(CDFID,'Q',NCFLOAT,4,DIMS,IRET)
         Call NCAPTC(CDFID,FLUXID,'long_name',NCCHAR,40,
     *              'Flux across transect defined in run_data',IRET)
         Call NCAPTC(CDFID,FLUXID,'units',NCCHAR,6,'m3 s-1',IRET)
        End If

        IF(EPTS.GT.0) THEN ! NG02192010
        DIMS(1) = STA2DDIM
        DIMS(2) = DIMDIM
        LOC2DID = NCVDEF(CDFID,'2Dloc',NCLONG,2,DIMS,IRET)
        Call NCAPTC(CDFID,LOC2DID,'long_name',NCCHAR,42,
     *              'I J locations of 2D stations on model grid',IRET)
        ENDIF ! NG02192010

        IF(VPTS.GT.0) THEN ! NG02192010
        DIMS(1) = STA3DDIM
        DIMS(2) = DIMDIM
        LOC3DID = NCVDEF(CDFID,'3Dloc',NCLONG,2,DIMS,IRET)
        Call NCAPTC(CDFID,LOC3DID,'long_name',NCCHAR,42,
     *              'I J locations of 3D stations on model grid',IRET)
        ENDIF ! NG02192010

!
! Global attributes
!NG03312009 CF Vs 1.4 and COARDS CONVENTIONS GLOBAL ATTRIBUTES
!NG03312009 CF stands for Climate and Forecast; Metadata Convention
        Call NCAPTC(CDFID,NCGLOBAL,'file_type',NCCHAR,7,
     *      "Station",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'cdm_datatype',NCCHAR,5,
     *      "Point",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'Conventions',NCCHAR,6,
     *      "CF-1.4",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'Conventions_comment',NCCHAR,21,
     *      "The version is CF-1.4",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'grid_type',NCCHAR,11,
     *      "Curvilinear",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'z_type',NCCHAR,5,
     *      "Sigma",IRET)
        modelname='New York Harbor Observing and Prediction System '
     *   //'(NYHOPS)'
        Call NCAPTC(CDFID,NCGLOBAL,'model',NCCHAR,80,
     *      modelname,IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'title',NCCHAR,30,
     *      "NYHOPS Forecast Model Results",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'comment',NCCHAR,42,
     *      "Once-daily 24hr Hindcast and 48hr Forecast",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'source',NCCHAR,10,
     *      "S::POM_vs3",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'institution',NCCHAR,31,
     *      "Stevens Institute of Technology",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'history',NCCHAR,29,
     *      "ECOMSED/Princeton Ocean Model",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'references',NCCHAR,55,
     *      "Alan.Blumberg@stevens.edu, Nickitas.Georgas@stevens.edu"
     *      ,IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'related_url',NCCHAR,39,
     *      "http://www.stevens.edu/maritimeforecast"
     *      ,IRET)
        CALL DATE_AND_TIME(NOWDATE)
        Call NCAPTC(CDFID,NCGLOBAL,'creation_date',NCCHAR,9,
     *      NOWDATE,IRET)
!NG        write (*,*) 'Basedate',IYR,IMO,IDA,IHR
        Call NCAPT(CDFID,NCGLOBAL,'base_date',NCLONG,4,BASE_DATE,IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'type_of_run',NCCHAR,10,TOR,IRET)
        Call NCAPT(CDFID,NCGLOBAL,'ishapiro',NCLONG,1,ISHAPIRO,IRET)
        Call NCAPT(CDFID,NCGLOBAL,'imeanprof',NCLONG,1,IMEANPROF,IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'hormix_type',NCCHAR,10,HORZMIX,
     *     IRET)
        Call NCAPT(CDFID,NCGLOBAL,'hormix_value',NCFLOAT,1,HORCON,
     *     IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'vrtmix_type',NCCHAR,10,VERTMIX,
     *     IRET)
        Call NCAPT(CDFID,NCGLOBAL,'vrtmix_value',NCFLOAT,1,UMOL,
     *     IRET)
        Call NCAPT(CDFID,NCGLOBAL,'tlag',NCFLOAT,1,1./tlag,IRET)
        Call NCAPT(CDFID,NCGLOBAL,'tlax_w',NCFLOAT,1,TRLXW,IRET)
        Call NCAPT(CDFID,NCGLOBAL,'tlax_s',NCFLOAT,1,TRLXS,IRET)
        Call NCAPT(CDFID,NCGLOBAL,'tlax_e',NCFLOAT,1,TRLXE,IRET)
        Call NCAPT(CDFID,NCGLOBAL,'tlax_n',NCFLOAT,1,TRLXN,IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'DATA_TYPE',NCCHAR,5,'MODEL',IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'COORD_SYSTEM',NCCHAR,3,'UTC',IRET)
        Call NCENDF(CDFID,IRET)
!=========================================================================
!         END OF DEFINE MODE
!=========================================================================
      Call NCCLOS(CDFID,IRET)
      Endif ! ICDF=1 NG03182010

      CDFID = NCOPN(CDFFILE,NCWRITE,IRET)
      IFILL = NCSFIL(CDFID,NCNOFILL,IRET)

      IF (ICDF.EQ.1) THEN ! NG03182010

!NG     Elevation stations ! NG 11/07/02
!NG     Previously used to store everything.

         IF (EPTS.GT.0) THEN ! NG02192010

         CORNER(1) = 1
         EDGES(1) = EPTS
!NG08242007 Store x and y coordinates of centroids
!svv         OPEN(502,FILE='corner_loc',STATUS='OLD',ERR=899)
!         DO K=1,(IM+1)*(JM+1) ! NG07252007
!            READ(502,'(2i5,2f12.5)',END=899)I,J,RJUNKX,RJUNKY
!            IF ((I.LE.IM).AND.(J.LE.JM)) THEN
!             XGRID(I,J)=RJUNKX
!             YGRID(I,J)=RJUNKY
!            ENDIF
!         END DO
! 899     CLOSE (502)
         
      
         Do I = 1, EPTS
          XCENTROID(I)=(XGRID(INXIE(I),INXJE(I))) !+
!svv     +                  XGRID(INXIE(I)+1,INXJE(I))+
!     +                  XGRID(INXIE(I)+1,INXJE(I)+1)+
!     +                  XGRID(INXIE(I),INXJE(I)+1))/4.
          YCENTROID(I)=(YGRID(INXIE(I),INXJE(I))) !+
!     +                  YGRID(INXIE(I)+1,INXJE(I))+
!     +                  YGRID(INXIE(I)+1,INXJE(I)+1)+
!     +                  YGRID(INXIE(I),INXJE(I)+1))/4.
          BATHY(I)=H(INXIE(I),INXJE(I))-WETMIN
         End Do
         Call NCVPT(CDFID,X2DID,CORNER,EDGES,XCENTROID,IRET)
         Call NCVPT(CDFID,Y2DID,CORNER,EDGES,YCENTROID,IRET)
         Call NCVPT(CDFID,DEPTH2DID,CORNER,EDGES,BATHY,IRET)
!NG08242007 End
         Do I = 1, EPTS
           NDUM(I) = I
         End Do
!NG04262007        Call NCVPT(CDFID,DEPTHID,CORNER,EDGES,DZSAVE,IRET)
         Call NCVPT1(CDFID,DUM1ID,CORNER,1,IRET)
         Call NCVPT1(CDFID,DUM2ID,CORNER,1,IRET)
         Call NCVPT1(CDFID,DUM3ID,CORNER,1,IRET)
         NDUM2(1) = 1
         NDUM2(2) = 2
         Call NCVPT(CDFID,DIMID,CORNER,2,NDUM2,IRET)
         Call NCVPT(CDFID,STA2DID,CORNER,EDGES,NDUM,IRET)
! store i coordinates for elevation stations
         CORNER(2) = 1
         EDGES(2) = 1
         Call NCVPT(CDFID,LOC2DID,CORNER,EDGES,INXIE,IRET)
! store j coordinates for elevation stations
         CORNER(2) = 2
         Call NCVPT(CDFID,LOC2DID,CORNER,EDGES,INXJE,IRET)
         ENDIF ! NG02192010

         IF (VPTS.GT.0) THEN ! NG02192010
!NG     Velocity, T, S, etc. stations ! NG 11/07/02
         CORNER(1) = 1
         EDGES(1) = VPTS
         Do I = 1, VPTS
          XCENTROID(I)=(XGRID(INXIV(I),INXJV(I)))
          YCENTROID(I)=(YGRID(INXIV(I),INXJV(I)))
          BATHY(I)=H(INXIV(I),INXJV(I))-WETMIN
         End Do
         Call NCVPT(CDFID,X3DID,CORNER,EDGES,XCENTROID,IRET)
         Call NCVPT(CDFID,Y3DID,CORNER,EDGES,YCENTROID,IRET)
         Call NCVPT(CDFID,DEPTH3DID,CORNER,EDGES,BATHY,IRET)
         Do I = 1, VPTS
           NDUM(I) = I
         End Do
         Call NCVPT1(CDFID,DUM1ID,CORNER,1,IRET)
         Call NCVPT1(CDFID,DUM2ID,CORNER,1,IRET)
         Call NCVPT1(CDFID,DUM3ID,CORNER,1,IRET)
         NDUM2(1) = 1
         NDUM2(2) = 2
         Call NCVPT(CDFID,DIMID,CORNER,2,NDUM2,IRET)
         Call NCVPT(CDFID,STA3DID,CORNER,EDGES,NDUM,IRET)
! store i coordinates for velocity etc. stations
         CORNER(2) = 1
         EDGES(2) = 1
         Call NCVPT(CDFID,LOC3DID,CORNER,EDGES,INXIV,IRET)
! store j coordinates for velocity etc. stations
         CORNER(2) = 2
         Call NCVPT(CDFID,LOC3DID,CORNER,EDGES,INXJV,IRET)
!NG     End of ! NG 11/07/02 Addition
         ENDIF ! NG02192010
! store stretched vertical coordinate levels (sigma-levels)
         CORNER(1) = 1
         EDGES(1) = KB ! NG 11/05/02
         Call NCVPT(CDFID,SIGMAID,CORNER,EDGES,Z,IRET)
         
      ENDIF ! ICDF=1 NG03182010

! store time
      CORNER(1) = ICDF
      TID = NCVID(CDFID,'time',IRET)
      Call NCVPT1(CDFID,TID,CORNER,NINT(TMIDDLE*3600*24),IRET)

!NG02192010 2D VARS
      IF(EPTS.GT.0) THEN ! NG02192010
! store wind stress
        CORNER(1) = 1
        CORNER(2) = 1
        CORNER(3) = 1
        CORNER(4) = ICDF
        EDGES(1) = EPTS
        EDGES(2) = 1
        EDGES(3) = 1
        EDGES(4) = 1
        Do I = 1, EPTS  
           RADS = ANG(INXIV(I),INXJV(I))
           WA1(I) = COS(RADS) * TXSAVE(I) - SIN(RADS) * TYSAVE(I)
        End Do
        TXID = ncvid(cdfid,'wtaux',iret)
        Call NCVPT(CDFID,TXID,CORNER,EDGES,WA1,IRET)
        Do I = 1, EPTS  
           RADS = ANG(INXIV(I),INXJV(I))
           WA1(I) = SIN(RADS) * TXSAVE(I) + COS(RADS) * TYSAVE(I)
        End Do
        TYID = ncvid(cdfid,'wtauy',iret)
        Call NCVPT(CDFID,TYID,CORNER,EDGES,WA1,IRET)

! store wind ! NG03312009
        CORNER(1) = 1
        CORNER(2) = 1
        CORNER(3) = 1
        CORNER(4) = ICDF
        EDGES(1) = EPTS
        EDGES(2) = 1
        EDGES(3) = 1
        EDGES(4) = 1
        Do I = 1, EPTS  
!NG           RADS = ANG(INXIV(I),INXJV(I))
!NG           WA1(I) = COS(RADS) * WSX2DSAVE(I) - SIN(RADS) * WSY2DSAVE(I)
           WA1(I) = WSX2DSAVE(I)
        End Do
        UWNDID = ncvid(cdfid,'wu',iret)
        Call NCVPT(CDFID,UWNDID,CORNER,EDGES,WA1,IRET)
        Do I = 1, EPTS  
!NG           RADS = ANG(INXIV(I),INXJV(I))
!NG           WA1(I) = SIN(RADS) * WSX2DSAVE(I) + COS(RADS) * WSY2DSAVE(I)
           WA1(I) = WSY2DSAVE(I)
        End Do
        VWNDID = ncvid(cdfid,'wv',iret)
        Call NCVPT(CDFID,VWNDID,CORNER,EDGES,WA1,IRET)

! store atmospheric pressure ! NG01172008
        CORNER(1) = 1
        CORNER(2) = 1
        CORNER(3) = 1
        CORNER(4) = ICDF
        EDGES(1) = EPTS
        EDGES(2) = 1
        EDGES(3) = 1
        EDGES(4) = 1
        Do I = 1, EPTS  
         WA1(I) = PATMSAVE(I)
        End Do
        PATMID = ncvid(cdfid,'patm',iret)
!      print *, 'writing patm ', patmsave(1)
        Call NCVPT(CDFID,PATMID,CORNER,EDGES,WA1,IRET)
! store rainfall rate ! NG09212011
        Do I = 1, EPTS  
         WA1(I) = QPRECSAVE(I)
        End Do
        QPRECID = ncvid(cdfid,'qprec',iret)
!      print *, 'writing qprec ', qprecsave(1)
        Call NCVPT(CDFID,QPRECID,CORNER,EDGES,WA1,IRET)
! store evaporation rate ! NG09212011
        Do I = 1, EPTS  
         WA1(I) = QEVAPSAVE(I)
        End Do
        QEVAPID = ncvid(cdfid,'qevap',iret)
!      print *, 'writing qevap ', qevapsave(1)
        Call NCVPT(CDFID,QEVAPID,CORNER,EDGES,WA1,IRET)
!NG04092008 Store other 2D Met variables
! store air temperature
        Do I = 1, EPTS  
         WA1(I) = AIRTSAVE(I)
        End Do
        AIRTID = ncvid(cdfid,'airt',iret)
!      print *, 'writing airt ', airtsave(1)
        Call NCVPT(CDFID,AIRTID,CORNER,EDGES,WA1,IRET)
! store cloud cover
        Do I = 1, EPTS  
         WA1(I) = CLDSAVE(I)
        End Do
        CLDID = ncvid(cdfid,'cld',iret)
!      print *, 'writing cld ', cldsave(1)
        Call NCVPT(CDFID,CLDID,CORNER,EDGES,WA1,IRET)
! store relative humidity
        Do I = 1, EPTS  
         WA1(I) = RHUMSAVE(I)
        End Do
        RHUMID = ncvid(cdfid,'rhum',iret)
!      print *, 'writing rhum ', rhumsave(1)
        Call NCVPT(CDFID,RHUMID,CORNER,EDGES,WA1,IRET)
! store surface-penetrating shortwave radiation
        Do I = 1, EPTS  
         WA1(I) = SWOBSSAVE(I)
        End Do
        SWOBSID = ncvid(cdfid,'swobs',iret)
!      print *, 'writing swobs ', swobssave(1)
        Call NCVPT(CDFID,SWOBSID,CORNER,EDGES,WA1,IRET)
! store wave period
        Do I = 1, EPTS  
         WA1(I) = WVPDSAVE(I)
        End Do
        WVPDID = ncvid(cdfid,'wp',iret)
!      print *, 'writing wp ', wvpdsave(1)
        Call NCVPT(CDFID,WVPDID,CORNER,EDGES,WA1,IRET)
! store significant wave height
        Do I = 1, EPTS  
         WA1(I) = WVHTSAVE(I)
        End Do
        WVHTID = ncvid(cdfid,'wh',iret)
!      print *, 'writing wh ', wvhtsave(1)
        Call NCVPT(CDFID,WVHTID,CORNER,EDGES,WA1,IRET)
! store significant wave direction
        Do I = 1, EPTS  
         WA1(I) = WVDRSAVE(I)
        End Do
        WVDRID = ncvid(cdfid,'wd',iret)
!      print *, 'writing wd ', wvdrsave(1)
        Call NCVPT(CDFID,WVDRID,CORNER,EDGES,WA1,IRET)

! store elevation
        CORNER(1) = 1
        CORNER(2) = 1
        CORNER(3) = 1
        CORNER(4) = ICDF
        EDGES(1) = EPTS
        EDGES(2) = 1
        EDGES(3) = 1
        EDGES(4) = 1
        ELID = NCVID(CDFID,'elev',IRET)
        Call NCVPT(CDFID,ELID,CORNER,EDGES,ESAVE,IRET)
        If (PTHF.EQ.'Y'.OR.PTHF.EQ.'y') Then
          HFID = NCVID(CDFID,'heat_flux',IRET)
          Call NCVPT(CDFID,HFID,CORNER,EDGES,HF2SAVE,IRET)
        Endif
      ENDIF ! NG02192010

      IF(VPTS.GT.0) THEN ! NG02192010
        CORNER(1) = 1
        CORNER(2) = 1
        CORNER(3) = 1
        CORNER(4) = ICDF
        EDGES(1) = VPTS   ! NG 11/05/02
! store velocity station depths here NG04262007
        EDGES(2) = 1
        EDGES(3) = 1
        EDGES(4) = 1
        DEPTHTID = NCVID(CDFID,'depth_t',IRET)
        Call NCVPT(CDFID,DEPTHTID,CORNER,EDGES,DZSAVE,IRET)
        EDGES(2) = 1
        EDGES(3) = KB ! NG 11/05/02
        EDGES(4) = 1
! store salinity
        If (PTS.EQ.'Y'.OR.PTS.EQ.'y') Then
           IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              WA1(IC) = SZSAVE(I,J)
            End Do
          End Do
          SALID = NCVID(CDFID,'salt',IRET)
          Call NCVPT(CDFID,SALID,CORNER,EDGES,WA1,IRET)
        End If
! store concentration
        If (PTC.EQ.'Y'.OR.PTC.EQ.'y') Then
           IC=0
          Do J = 1, KB        ! NG 11/08/02
          Do I = 1, VPTS  ! NG 11/05/02
            IC=IC+1
            WA1(IC) = CZSAVE(I,J)
          End Do
        End Do
        CID = NCVID(CDFID,'conc',IRET)
        Call NCVPT(CDFID,CID,CORNER,EDGES,WA1,IRET)
        END IF
! store temperature
        If (PTT.EQ.'Y'.OR.PTT.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              WA1(IC) = TZSAVE(I,J)
            End Do
          End Do
          TMPID = NCVID(CDFID,'temp',IRET)
          Call NCVPT(CDFID,TMPID,CORNER,EDGES,WA1,IRET)
        End If
! store vertical diffusivity
        If (PTKH.EQ.'Y'.OR.PTKH.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              WA1(IC) = KHZSAVE(I,J)
            End Do
          End Do
          KHID = NCVID(CDFID,'kh',IRET)
          Call NCVPT(CDFID,KHID,CORNER,EDGES,WA1,IRET)
        End If
! store vertical viscosity
        If (PTKM.EQ.'Y'.OR.PTKM.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              WA1(IC) = KMZSAVE(I,J)
            End Do
          End Do
          KMID = NCVID(CDFID,'km',IRET)
          Call NCVPT(CDFID,KMID,CORNER,EDGES,WA1,IRET)
        End If
! store horizontal viscosity
        If (PTAM.EQ.'Y'.OR.PTAM.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              WA1(IC) = AMZSAVE(I,J)
            End Do
          End Do
          AMID = NCVID(CDFID,'am',IRET)
          Call NCVPT(CDFID,AMID,CORNER,EDGES,WA1,IRET)
        End If
! store velocity
        If (PTU.EQ.'Y'.OR.PTU.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              RADS = ANG(INXIV(I),INXJV(I))       ! NG 11/07/02
            WA1(IC) = COS(RADS) * UZSAVE(I,J) - SIN(RADS) * VZSAVE(I,J)
            End Do
          End Do
          UID = NCVID(CDFID,'u',IRET)
          Call NCVPT(CDFID,UID,CORNER,EDGES,WA1,IRET)
        End If
        If (PTV.EQ.'Y'.OR.PTV.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
              RADS = ANG(INXIV(I),INXJV(I))       ! NG 11/07/02
            WA1(IC) = SIN(RADS) * UZSAVE(I,J) + COS(RADS) * VZSAVE(I,J)
            End Do
          End Do
          VID = NCVID(CDFID,'v',IRET)
          Call NCVPT(CDFID,VID,CORNER,EDGES,WA1,IRET)
        End If
        If (PTW.EQ.'Y'.OR.PTW.EQ.'y') Then
          IC=0
          Do J = 1, KB        ! NG 11/08/02
            Do I = 1, VPTS        ! NG 11/05/02
              IC=IC+1
!NG03312009              WA1(IC) = WZSAVE(I,J)  ! WZSAVE IS NOT DEFINED !!!!!
              IF(J.EQ.KB) THEN
               WA1(IC) = WZSAVE(I,J)
              ELSE
               WA1(IC) = 0.5*(WZSAVE(I,J)+WZSAVE(I,J+1))
              ENDIF
            End Do
          End Do
!NG03312009          WID = NCVID(CDFID,'w',IRET) ! NG CRASHES HERE
          WID = NCVID(CDFID,'w',IRET)
          Call NCVPT(CDFID,WID,CORNER,EDGES,WA1,IRET)
        End If
      ENDIF ! NG02192010
! store cross-sectional fluxes NG04022008
      If (FPTS.GT.0) Then
        EDGES(1) = FPTS
        EDGES(2) = 1
        EDGES(3) = KB
        EDGES(4) = 1
        IC=0
        Do J = 1, KB
          Do I = 1, FPTS
            IC=IC+1
            WA1(IC) = CCFLUX(I,J)
          End Do
        End Do
        FLUXID = NCVID(CDFID,'Q',IRET)
        Call NCVPT(CDFID,FLUXID,CORNER,EDGES,WA1,IRET)
      End If

      Call NCCLOS(CDFID,IRET)
!      write (*,*) TMIDDLE,ICDF
      ICDF = ICDF + 1
      RETURN ! NG 11/05/02
      END
