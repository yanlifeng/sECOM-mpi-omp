      SUBROUTINE PUTCDF(TMID)
! This program is a modified version of 'putcdf.f' which was used
! to generate netCDF format output of ECOM. USED FOR MC-02, 2002
!     Store time-averaged variables at all grid points in
!     netCDF data file
!     (Stevens Institute of Technology)
!                                     June 25, 2003
!
!NG03312009 Brought to compliance with the latest
!NG03312009 NetCDF Climate and Forecast (CF)
!NG03312009 MetaData Conventions Version 1.4, 27 February, 2009
!NG03312009
!NG03132009 http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.4/cf-conventions.pdf
!NG03312009

      Parameter (ELMIN = -10.,ELMAX = 10.)
      Parameter (SALMIN = -5.,SALMAX = 80.)
      Parameter (TMPMIN = -10.,TMPMAX = 50.)
      Parameter (VELMIN = -6.,VELMAX = 6.)
      Parameter (WTMIN = -.1,WTMAX = 0.1)
      Parameter (FLXMIN = -7200,FLXMAX = 7200.)
      Parameter (RKHMIN = 0.,RKHMAX = 0.2)
      Parameter (RKMMIN = 0.,RKMMAX = 0.2)
      Parameter (AMMIN = 0.,AMMAX = 1000.)
      Parameter (CTMIN = 0., CTMAX = 1.)
      Parameter (BFMIN = 0.,BFMAX = 10000.) ! Bottom Friction x 1e4
      Parameter (WPMIN = -10.,WPMAX = 7200.) ! WAVE PERIODup to 2 hours (sec)
      Parameter (WDMIN = -10.,WDMAX = 360.)  ! WAVE DIRECTION (deg)
      Parameter (WHMIN =  -1.,WHMAX = 30.)   ! WAVE HEIGHT (m)
      Parameter (WNDMIN = -150.,WNDMAX = 150.) ! WIND MAG (m/s)
      Parameter (PATMMIN = 800.,PATMMAX = 1100.) ! MSL ATMOSPHERIC PRESSURE MAG (mbar) ! NG01172008
      Parameter (QPRECMIN = 0.,QPRECMAX = 250.) ! MSL ATMOSPHERIC PRESSURE MAG (mbar) ! NG01172008
      Parameter (AIRTMIN = -90.,AIRTMAX = 60.) ! AIR TEMPERATURE (C) ! NG04092008
      Parameter (CLDMIN = 0.,CLDMAX = 1.) ! CLOUD COVER (FRACTION) ! NG04092008
      Parameter (RHUMMIN = 0.,RHUMMAX = 100.) ! RELATIVE HUMIDITY (%) ! NG04092008
      Parameter (SWOBSMIN = -50.,SWOBSMAX = 2500.) ! SURFACE-PENETRATING SHORT WAVE RADIATION (W/m2) ! NG04092008
      Parameter (HFMIN = -4000.,HFMAX = 4000.) ! Careful: NOT IMPLEMENTED YET
      Parameter (Q2MIN = 0,Q2MAX = 1E-2) ! Careful: NOT IMPLEMENTED YET
      Parameter (SHORTRANGE = 2.*32765.) ! NG03312009 CF 1.4 CONVENTIONS, 2 OFF FROM FILL_Values of -32767
      Parameter (IFILLVAL = -32767) ! NG06202013
      Parameter (RFILLVAL = -99999.) ! NG06202013
      Include 'comdeck'
      Include 'netcdf.inc'
      integer i1,j1,k1,i2,j2,k2
      Integer IRET, BASE_DATE(4)
* netCDF id
      Integer CDFID
* dimension ids
      Integer XPOSDIM, YPOSDIM, ZPOSDIM, TIMEDIM, ZSDVDIM, CBNDDIM ! NG03312009
* variable ids
      Integer XID, YID, SIGMAID, DEPTHID, TID, ELID, UID, VID, WID, 
     *    SALID, TMPID, UFLXID, VFLXID, KHID, ANGID,BTYPEID,
     *     H1ID, H2ID, TXID, TYID, KMID, AMID, HFID, CCID,BFID,
     *     Q2ID,DATEID,WPID,WDID,WHID,UWNDID,VWNDID, FSMID,PATMID
     *    ,QPRECID ! NG09212011
     *    ,XCID,YCID ! Store centers too. NG04022008
     *    ,AIRTID,CLDID,RHUMID,SWOBSID ! 2D MET NG04092008
     *    ,DUMID,DVMID ! FOR CELL CONNECTIVITY (E.G. THIN DAMS) NG08072008
     *    ,XPOSID, YPOSID ! CNG03312009 FOR CARTESIAN CF 1.4
     *    ,LNBNDID, LTBNDID ! CNG03312009 FOR CELL BOUNDS
     *    ,EXPLTID, DEPTHTID ! CNG03312009 FOR EXPLICIT DEPTH MAPPING
     *    ,ZSDVID ! CNG05292009 FOR SDV SIGMA LEVELS
* variable shapes
      Integer DIMS(4)
* corners and edge lengths
      Integer CORNER(4), EDGES(4)
* attribute vectors
      Real FLOATVAL(1)
!NG11022008 CF convention      INTEGER*2 VRNG(2)
      Real VRNG(2)
      Integer*2 IRNG(2)
* short integer workspace arrays
      Integer*2 WS0(IM,JM)
      Integer*2 WA1(IM,JM,KB)
!   corners file
      REAL XC(IM,JM),YC(IM,JM)
      REAL CNTII(IM),CNTJJ(JM) ! CNG03312009 CF 1.4 CARTESIAN COORDS
      REAL LONBOUNDS(4,IM,JM),LATBOUNDS(4,IM,JM) ! CORNERS FOR CELL BOUNDS CF1.7
!   DUM,DVM NG08072008
      REAL DUM_IMJM(IM,JM),DVM_IMJM(IM,JM)
!   bottom types
      REAL BTYPE(IM,JM)
! Temporary local arrays holding 3D parameters
      REAL TCDF(IM,JM,KB),SCDF(IM,JM,KB)
!NG05292009 FOR SIGMA AT MID-POINTS
      REAL SIGMEAN(KB)
! long integer 
!NG04172008      Integer*4 DATEHOUR,IMONTH,IYY,IDAY,IHH,IMN
      REAL*8 DATEHOUR,DMONTH,DYY,DDAY,DHH,DMN
* scale and offset
      Real KHSC, KHOF, KMSC, KMOF
      CHARACTER*1 PTC,PTWND,PTPATM,PTQPREC,PTD
      CHARACTER*1 PT2DMET ! NG04092008 2D MET
! arrays for netCDF output: but we are not going to use this
! anyway  (NKIM 1.12.00) clean up later!
      REAL ARCC(IM,JM,KB),ARCHF(IM,JM)
!NG02232007     .     ARCAM(IM,JM,KB) ! moved to comdeck
      Integer E1
      real dx,dy,costheta,sintheta
!NG03312009
      CHARACTER modelname*80
      CHARACTER NOWDATE*9,BASEDATE*26
!NG03222010
      REAL UTEMP,VTEMP,VELU(IM,JM,KB),VELV(IM,JM,KB)
!NG06202013
      REAL HOUT(IM,JM)
* count the number of times this routine is called
      Data ICDF /0/
      Save
* first time through...
      If (ICDF.EQ.0) Then
!NG05292009 FOR SIGMA AT MID-POINTS
        DO K=1,KBM1
         SIGMEAN(K)=ZZ(K)
        ENDDO
        SIGMEAN(KB)=Z(KB) ! CLOSURE AT -1
!NG03312009 and 06202013
        HMIN=12000.
        HMAX=-9000.
        DO J=1,JM
         DO I=1,IM
          IF (FSM(I,J).EQ.0) THEN
           HOUT(I,J)=RFILLVAL ! TAKE OUT DEPTHS AT ALWAYS DRY POINTS
          ELSE
           HOUT(I,J)=H(I,J)-WETMIN ! TAKE BACK OUT WETMIN TO ORIGINAL DEPTH
           HMIN=AMIN1(HOUT(I,J),HMIN)
           HMAX=AMAX1(HOUT(I,J),HMAX)
          ENDIF 
         ENDDO
        ENDDO
        HMIN=HMIN+ELMIN
        HMAX=HMAX+ELMAX
        DPMIN=-HMAX
        DPMAX=-HMIN

!NG03312009
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
* Don't need T, S, W, and KH if run is barotropic (2-D)
        If (TOR.EQ.'BAROTROPIC') Then ! Always overwrite ! NG 11/08/02
          PTD = 'n'
          PTS = 'n'
          PTT = 'n'
          PTW = 'n'
          PTKH = 'n'
          PTKM = 'n'
          kbot = 1
       else
          PTD = 'n' ! NG03312009, set to y for explicit mapping, e.g.  for IDV
          kbot = kb
        End If
* calculate scale factors and offsets based on expected range of data
        DPOF = (DPMAX+DPMIN) * .5
        ELOF = (ELMAX+ELMIN) * .5
        HFOF = (HFMAX+HFMIN) * .5
        BFOF = (BFMAX+BFMIN) * .5
        SALOF = (SALMAX+SALMIN) * .5
        TMPOF = (TMPMAX+TMPMIN) * .5
        CTOF = (CTMAX+CTMIN) * .5
        VELOF = (VELMAX+VELMIN) * .5
        WTOF = (WTMAX+WTMIN) * .5
        FLXOF = (FLXMAX+FLXMIN) * .5
        KHOF = (RKHMAX+RKHMIN) * .5
        KMOF = (RKMMAX+RKMMIN) * .5
        AMOF = (AMMAX+AMMIN) * .5
        Q2OF = (Q2MAX+Q2MIN) * .5
        WPOF = (WPMAX+WPMIN) * .5
        WDOF = (WDMAX+WDMIN) * .5
        WHOF = (WHMAX+WHMIN) * .5
        WNDOF = (WNDMAX+WNDMIN) * .5
        PATMOF = (PATMMAX+PATMMIN) * .5
        QPRECOF = (QPRECMAX+QPRECMIN) * .5
        AIRTOF = (AIRTMAX+AIRTMIN) * .5
        CLDOF = (CLDMAX+CLDMIN) * .5
        RHUMOF = (RHUMMAX+RHUMMIN) * .5
        SWOBSOF = (SWOBSMAX+SWOBSMIN) * .5
        DPSC = (DPMAX-DPMIN) / SHORTRANGE
        ELSC = (ELMAX-ELMIN) / SHORTRANGE
        HFSC = (HFMAX-HFMIN) / SHORTRANGE
        BFSC = (BFMAX-BFMIN) / SHORTRANGE
        SALSC = (SALMAX-SALMIN) / SHORTRANGE
        TMPSC = (TMPMAX-TMPMIN) / SHORTRANGE
        CTSC = (CTMAX-CTMIN) / SHORTRANGE
        VELSC = (VELMAX-VELMIN) / SHORTRANGE
        WTSC = (WTMAX-WTMIN) / SHORTRANGE
        FLXSC = (FLXMAX-FLXMIN) / SHORTRANGE
        KHSC = (RKHMAX-RKHMIN) / SHORTRANGE
        KMSC = (RKMMAX-RKMMIN) / SHORTRANGE
        AMSC = (AMMAX-AMMIN) / SHORTRANGE
        Q2SC = (Q2MAX-Q2MIN) / SHORTRANGE
        WPSC = (WPMAX-WPMIN) / SHORTRANGE
        WDSC = (WDMAX-WDMIN) / SHORTRANGE
        WHSC = (WHMAX-WHMIN) / SHORTRANGE
        WNDSC = (WNDMAX-WNDMIN) / SHORTRANGE
        PATMSC = (PATMMAX-PATMMIN) / SHORTRANGE
        QPRECSC = (QPRECMAX-QPRECMIN) / SHORTRANGE
        AIRTSC = (AIRTMAX-AIRTMIN) / SHORTRANGE
        CLDSC = (CLDMAX-CLDMIN) / SHORTRANGE
        RHUMSC = (RHUMMAX-RHUMMIN) / SHORTRANGE
        SWOBSSC = (SWOBSMAX-SWOBSMIN) / SHORTRANGE
* HARDWIRE: if one wants to save wind
        PTWND='y'
* HARDWIRE: if one wants to save atmospheric pressure
        PTPATM='y'
* HARDWIRE: if one wants to save precipitation rate
        PTQPREC='n'
* if 2DMET then save other 2D met variables NG04092008
        IF (OPTMBC(1:5).EQ.'2DMET') THEN
         PT2DMET='y'
        ENDIF
**1* HARDWIRE: don't save fluxes, but you could
        PTUF = 'n'
        PTVF = 'n'
* REMOVED HARDWIRE: don't save mixing parameters
!NG        PTKH = 'n'
!NG        PTKM = 'n'
!NG        PTAM = 'n'
* HARDWIRE: if one wants to save turbulence kinetic energy
        PTQ2 = 'n'
* HARDWIRE: if one wants to save heat flux
        PTHF = 'n'
* HARDWIRE: if one wants to save drag coefficient
        PTCD='y'
*  save concentration
!NG02222007
       IF(TRACER.EQ.'INCLUDE'.AND.WAVEDYN.NE.'DONONLY ')THEN
        PTC = 'y'
       ELSE
        PTC = 'n'
       END IF
* enter define mode
!NG        CDFID = NCCRE('gcmplt.nc',NCCLOB,IRET) ! NG03312009 changed name to .nc as per CF
        CDFID = NCCRE('gcmplt.cdf',NCCLOB,IRET) ! This is for NYHOPS, too much of a pain otherwise
        IFILL = NCSFIL(CDFID,NCNOFILL,IRET)
* define dimensions
        XPOSDIM = NCDDEF(CDFID,'xpos',IM,IRET)
        YPOSDIM = NCDDEF(CDFID,'ypos',JM,IRET)
        ZSDVDIM = NCDDEF(CDFID,'zpos',KBot,IRET) ! CNG03312009 FOR BACKCOMPATIBILITY WITH SDV
        ZPOSDIM = NCDDEF(CDFID,'sigma',KBot,IRET) ! CNG03312009 FOR CF  1.4 CONVENTIONS
        TIMEDIM = NCDDEF(CDFID,'time',NCUNLIM,IRET)
        CBNDDIM = NCDDEF(CDFID,'vertices',4,IRET) ! CNG03312009 4 cell-bounding certices 
* define variables
*   1D Vars 
!        print *, "timedim=", timedim
!NG03312009 DEFINE CARTESIAN COORDINATES TOO, CF1.4 CONVENTIONS
        DIMS(1) = XPOSDIM
        XPOSID = NCVDEF(CDFID,'xpos',NCFLOAT,1,DIMS,IRET)
        DIMS(1) = YPOSDIM
        YPOSID = NCVDEF(CDFID,'ypos',NCFLOAT,1,DIMS,IRET)

        DIMS(1) = TIMEDIM
        TID = NCVDEF(CDFID,'time',NCFLOAT,1,DIMS,IRET)
        DATEID = NCVDEF(CDFID,'date',NCDOUBLE,1,DIMS,IRET)
!       TXID = NCVDEF(CDFID,'taux',NCFLOAT,1,DIMS,IRET)  !NKIM 1.12.00
!       TYID = NCVDEF(CDFID,'tauy',NCFLOAT,1,DIMS,IRET)  !NKIM 1.12.00
!NG02222007
!NG03092007        IF(WAVEDYN.NE.'DONONLY ')THEN
         DIMS(1) = ZSDVDIM
         ZSDVID = NCVDEF(CDFID,'layer_bnds',NCFLOAT,1,DIMS,IRET)
         DIMS(1) = ZPOSDIM
         SIGMAID = NCVDEF(CDFID,'sigma',NCFLOAT,1,DIMS,IRET)
!NG03092007        ENDIF
*   2D Vars
        DIMS(2) = YPOSDIM
        DIMS(1) = XPOSDIM
        XID = NCVDEF(CDFID,'x',NCFLOAT,2,DIMS,IRET)
        YID = NCVDEF(CDFID,'y',NCFLOAT,2,DIMS,IRET)
        XCID = NCVDEF(CDFID,'lon',NCFLOAT,2,DIMS,IRET) ! NG04022008
        YCID = NCVDEF(CDFID,'lat',NCFLOAT,2,DIMS,IRET) ! NG04022008
        H1ID = NCVDEF(CDFID,'h1',NCFLOAT,2,DIMS,IRET)
        H2ID = NCVDEF(CDFID,'h2',NCFLOAT,2,DIMS,IRET)
        DEPTHID = NCVDEF(CDFID,'depth',NCFLOAT,2,DIMS,IRET)
        ANGID = NCVDEF(CDFID,'ang',NCFLOAT,2,DIMS,IRET)
        FSMID = NCVDEF(CDFID,'FSM',NCFLOAT,2,DIMS,IRET)
        DUMID = NCVDEF(CDFID,'DUM',NCFLOAT,2,DIMS,IRET) ! NG08072008
        DVMID = NCVDEF(CDFID,'DVM',NCFLOAT,2,DIMS,IRET) ! NG08072008
!       BTYPEID = NCVDEF(CDFID,'btype',NCFLOAT,2,DIMS,IRET)
*   3D Vars
!NG03312009 NEW VARS FOR BOUNDS
        DIMS(3) = YPOSDIM
        DIMS(2) = XPOSDIM
        DIMS(1) = CBNDDIM ! This needs to be the fastest-varying dimension in CF-1.7, NG 04/23/2014

        LNBNDID = NCVDEF(CDFID,'lon_bnds',NCFLOAT,3,DIMS,IRET)
        Call NCAPTC(CDFID,LNBNDID,'long_name',NCCHAR,16,
     *             'Vertex longitude',IRET)
        Call NCAPTC(CDFID,LNBNDID,'units',NCCHAR,12,'degrees_east',IRET)
        VRNG(1) = -361.
        VRNG(2) = 361.
        Call NCAPT(CDFID,LNBNDID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,LNBNDID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)

        LTBNDID = NCVDEF(CDFID,'lat_bnds',NCFLOAT,3,DIMS,IRET)
        Call NCAPTC(CDFID,LTBNDID,'long_name',NCCHAR,15,
     *             'Vertex latitude',IRET)
        Call NCAPTC(CDFID,LTBNDID,'units',NCCHAR,13,'degrees_north'
     *              ,IRET)
        VRNG(1) = -361.
        VRNG(2) = 361.
        Call NCAPT(CDFID,LTBNDID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,LTBNDID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)

        DIMS(3) = TIMEDIM
        DIMS(2) = YPOSDIM
        DIMS(1) = XPOSDIM
!NG02222007
        IF(WAVEDYN.NE.'DONONLY ')THEN
         ELID = NCVDEF(CDFID,'elev',NCSHORT,3,DIMS,IRET)
        ENDIF
*   1D, 2D, and 3D attributes
        Call NCAPTC(CDFID,XID,'long_name',NCCHAR,16,
     *             'Corner longitude',IRET)
        Call NCAPTC(CDFID,XID,'units',NCCHAR,12,'degrees_east',IRET)
        VRNG(1) = -361.
        VRNG(2) = 361.
        Call NCAPT(CDFID,XID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,XID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)

        Call NCAPTC(CDFID,XCID,'long_name',NCCHAR,11,'center long',IRET)! NG04022008
        Call NCAPTC(CDFID,XCID,'units',NCCHAR,12,'degrees_east',IRET) ! NG04022008
        Call NCAPTC(CDFID,XCID,'standard_name',NCCHAR,9,
     *              'longitude',IRET)
        Call NCAPTC(CDFID,XCID,'bounds',NCCHAR,8,
     *              'lon_bnds',IRET)
        VRNG(1) = -361.
        VRNG(2) = 361.
        Call NCAPT(CDFID,XCID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,XCID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)

        Call NCAPTC(CDFID,H1ID,'long_name',NCCHAR,9,'dx metric',IRET)
        Call NCAPTC(CDFID,H1ID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,H1ID,'coordinates',NCCHAR,7,
     *              'lat lon',IRET)

        Call NCAPTC(CDFID,YID,'long_name',NCCHAR,15,
     *             'Corner latitude',IRET)
        Call NCAPTC(CDFID,YID,'units',NCCHAR,13,'degrees_north',IRET)
        VRNG(1) = -361.
        VRNG(2) = 361.
        Call NCAPT(CDFID,YID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,YID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)

        Call NCAPTC(CDFID,YCID,'long_name',NCCHAR,10,'center lat',IRET) ! NG04022008
        Call NCAPTC(CDFID,YCID,'units',NCCHAR,13,'degrees_north',IRET) ! NG04022008
        Call NCAPTC(CDFID,YCID,'standard_name',NCCHAR,8,
     *              'latitude',IRET)
        Call NCAPTC(CDFID,YCID,'bounds',NCCHAR,8,
     *              'lat_bnds',IRET)
        VRNG(1) = -361.
        VRNG(2) = 361.
        Call NCAPT(CDFID,YCID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,YCID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)

        Call NCAPTC(CDFID,H2ID,'long_name',NCCHAR,9,'dy metric',IRET)
        Call NCAPTC(CDFID,H2ID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,H2ID,'coordinates',NCCHAR,7,
     *              'lat lon',IRET)
!NG02222007
!NG03092007        IF(WAVEDYN.NE.'DONONLY ')THEN

        Call NCAPTC(CDFID,XPOSID,'long_name',NCCHAR,32,
     *              'X-coordinate in Cartesian system',IRET)
        Call NCAPTC(CDFID,XPOSID,'units',NCCHAR,1,'m',IRET)
        Call NCAPTC(CDFID,XPOSID,'axis',NCCHAR,1,'X',IRET)
        Call NCAPTC(CDFID,YPOSID,'long_name',NCCHAR,32,
     *              'Y-coordinate in Cartesian system',IRET)
        Call NCAPTC(CDFID,YPOSID,'units',NCCHAR,1,'m',IRET)
        Call NCAPTC(CDFID,YPOSID,'axis',NCCHAR,1,'Y',IRET)

!NG05292009 Store Sigma Levels (bounds) also for back-compatibility
        Call NCAPTC(CDFID,ZSDVID,'long_name',NCCHAR,46,
     *       'bounds of stretched vertical coordinate levels',IRET)
        Call NCAPTC(CDFID,ZSDVID,'units',NCCHAR,1,'1',
     *      IRET)
        Call NCAPTC(CDFID,ZSDVID,'positive',NCCHAR,2,
     *              'up',IRET)

!NGCF CAREFUL IF COARDS        Call NCAPTC(CDFID,SIGMAID,'units',NCCHAR,11,'sigma_level',
!NGCF CAREFUL IF COARDS     *      IRET)
        Call NCAPTC(CDFID,SIGMAID,'long_name',NCCHAR,36,
     *       'stretched vertical coordinate levels',IRET)
        Call NCAPTC(CDFID,SIGMAID,'units',NCCHAR,1,'1',
     *      IRET)
        Call NCAPTC(CDFID,SIGMAID,'positive',NCCHAR,2,
     *              'up',IRET)
        Call NCAPTC(CDFID,SIGMAID,'standard_name',NCCHAR,22,
     *              'ocean_sigma_coordinate',IRET)
        Call NCAPTC(CDFID,SIGMAID,'formula_terms',NCCHAR,35,
     *              'sigma: sigma eta: elev depth: depth',IRET)
        Call NCAPTC(CDFID,SIGMAID,'axis',NCCHAR,1,'Z',IRET)
!NG03312009 The following are IDV specific (?)
!NG        Call NCAPTC(CDFID,SIGMAID,'_CoordinateAxisType',NCCHAR,4,
!NG     *              'GeoZ',IRET)
!NG        Call NCAPTC(CDFID,SIGMAID,'_CoordinateZisPositive',NCCHAR,2,
!NG     *              'up',IRET)
!NG        Call NCAPTC(CDFID,SIGMAID,'_CoordinateTransformType',NCCHAR,8,
!NG     *              'vertical',IRET)
!NG        Call NCAPTC(CDFID,SIGMAID,'_CoordinateAxes',NCCHAR,5,
!NG     *              'sigma',IRET)

!NG03092007        ENDIF
        Call NCAPTC(CDFID,DEPTHID,'long_name',NCCHAR,10,'Bathymetry',
     *      IRET)
        Call NCAPTC(CDFID,DEPTHID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,DEPTHID,'standard_name',NCCHAR,5,
     *              'depth',IRET)
        Call NCAPTC(CDFID,DEPTHID,'positive',NCCHAR,4,
     *              'down',IRET)
        Call NCAPTC(CDFID,DEPTHID,'coordinates',NCCHAR,7,
     *              'lat lon',IRET)
        Call NCAPTC(CDFID,ANGID,'long_name',NCCHAR,10,'grid angle',IRET)
        Call NCAPTC(CDFID,ANGID,'units',NCCHAR,7,'radians',IRET)
        Call NCAPTC(CDFID,ANGID,'standard_name',NCCHAR,32,
     *              'angle_of_rotation_from_east_to_x',IRET)
        Call NCAPTC(CDFID,ANGID,'coordinates',NCCHAR,7,
     *              'lat lon',IRET)
        Call NCAPTC(CDFID,FSMID,'long_name',NCCHAR,17,
     *       'free surface mask',IRET)
        Call NCAPTC(CDFID,FSMID,'coordinates',NCCHAR,7,
     *              'lat lon',IRET)
        Call NCAPTC(CDFID,DUMID,'long_name',NCCHAR,17,
     *       'U1-direction mask',IRET) ! NG08072008
        Call NCAPTC(CDFID,DVMID,'long_name',NCCHAR,17,
     *       'V1-direction mask',IRET) ! NG08072008
!       Call NCAPTC(CDFID,BTYPEID,'long_name',NCCHAR,12,'bottom types',
!    *      IRET)
!       Call NCAPTC(CDFID,BTYPEID,'units',NCCHAR,9,'see_table',IRET)
        VRNG(1) = HMIN ! CAREFUL HERE, ONLY FOR DEPTH NG03312009
        VRNG(2) = HMAX ! CAREFUL HERE, ONLY FOR DEPTH NG03312009
!        Write (*,*) 'hmin,hmax', HMIN, HMAX
        Call NCAPT(CDFID,DEPTHID,'valid_range',NCFLOAT,2,VRNG,IRET)
        Call NCAPT(CDFID,DEPTHID,'_FillValue',NCFLOAT,1,RFILLVAL,IRET)
        Call NCAPTC(CDFID,TID,'long_name',NCCHAR,4,'Time',IRET)
        Call NCAPT(CDFID,TID,'base_date',NCLONG,4,BASE_DATE,IRET)
        Call NCAPTC(CDFID,TID,'units',NCCHAR,37,
     *              'days since '//BASEDATE,IRET)
        Call NCAPTC(CDFID,TID,'standard_name',NCCHAR,4,
     *              'time',IRET)
        Call NCAPTC(CDFID,TID,'axis',NCCHAR,1,'T',IRET)
        Call NCAPTC(CDFID,DATEID,'long_name',NCCHAR,9,'Date_Time',IRET)
!       Call NCAPTC(CDFID,TXID,'long_name',NCCHAR,16,'East Wind Stress',
!    *      IRET)
!       Call NCAPTC(CDFID,TXID,'units',NCCHAR,6,'pascal',IRET)
!       Call NCAPTC(CDFID,TYID,'long_name',NCCHAR,17,
!    *      'North Wind Stress',IRET)
!       Call NCAPTC(CDFID,TYID,'units',NCCHAR,6,'pascal',IRET)
!       VRNG(1) = ELMIN
!       VRNG(2) = ELMAX
!        Write (*,*) 'elmin,elmax', ELMIN, ELMAX

!NG02222007
        IF(WAVEDYN.NE.'DONONLY ')THEN
        IRNG(1) = NINT((ELMIN-ELOF)/ELSC)
        IRNG(2) = NINT((ELMAX-ELOF)/ELSC)
!        Write (*,*) 'elmin,elmax', ELMIN, ELMAX
        Call NCAPT(CDFID,ELID,'valid_range',NCSHORT,2,IRNG,IRET) 
        Call NCAPT(CDFID,ELID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
        Call NCAPTC(CDFID,ELID,'long_name',NCCHAR,9,'Elevation',IRET)
        Call NCAPTC(CDFID,ELID,'units',NCCHAR,6,'meters',IRET)
        Call NCAPTC(CDFID,ELID,'standard_name',NCCHAR,34,
     *              'sea_surface_height_above_sea_level',IRET)
        Call NCAPTC(CDFID,ELID,'positive',NCCHAR,2,
     *              'up',IRET)
        Call NCAPTC(CDFID,ELID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
        FLOATVAL(1) = ELSC
        Call NCAPT(CDFID,ELID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
        FLOATVAL(1) = ELOF
        Call NCAPT(CDFID,ELID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        If (PTHF.EQ.'y'.OR.PTHF.EQ.'Y') Then
* define heatflux
          DIMS(3) = TIMEDIM
          DIMS(2) = YPOSDIM
          DIMS(1) = XPOSDIM
          HFID = NCVDEF(CDFID,'heat_flux',NCSHORT,3,DIMS,IRET)
!         Write (*,*) 'hfmin,hfmax', HFMIN, HFMAX
          IRNG(1) = NINT((HFMIN-HFOF)/HFSC)
          IRNG(2) = NINT((HFMAX-HFOF)/HFSC)
          Call NCAPT(CDFID,HFID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,HFID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,HFID,'long_name',NCCHAR,17,
     *         'heat_flux in W/m2',IRET)
          Call NCAPTC(CDFID,HFID,'units',NCCHAR,5,'W m-2',IRET)
          Call NCAPTC(CDFID,HFID,'standard_name',NCCHAR,39,
     *              'surface_downward_heat_flux_in_sea_water',IRET)
          Call NCAPTC(CDFID,HFID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = HFSC
          Call NCAPT(CDFID,HFID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = HFOF
          Call NCAPT(CDFID,HFID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTCD.EQ.'y'.OR.PTCD.EQ.'Y') Then
* define bottom drag
          DIMS(3) = TIMEDIM
          DIMS(2) = YPOSDIM
          DIMS(1) = XPOSDIM
          BFID = NCVDEF(CDFID,'cbc',NCSHORT,3,DIMS,IRET)
          IRNG(1) = NINT((BFMIN-BFOF)/BFSC)
          IRNG(2) = NINT((BFMAX-BFOF)/BFSC)
!         Write (*,*) 'cdmin,cdmax', BFMIN, BFMAX
          Call NCAPT(CDFID,BFID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,BFID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,BFID,'long_name',NCCHAR,24,
     &       'Bottom Drag Coefficient',IRET)
          Call NCAPTC(CDFID,BFID,'units',NCCHAR,4,'none',IRET)
          Call NCAPTC(CDFID,BFID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = BFSC
          Call NCAPT(CDFID,BFID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = BFOF
          Call NCAPT(CDFID,BFID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If

        ENDIF

!   Wind parameter
        If (PTWND.EQ.'y'.OR.PTWND.EQ.'Y') Then
!   First, EW wind component
          UWNDID = NCVDEF(CDFID,'wu',NCSHORT,3,DIMS,IRET)
          IRNG(1) = NINT((WNDMIN-WNDOF)/WNDSC)
          IRNG(2) = NINT((WNDMAX-WNDOF)/WNDSC)
          Call NCAPT(CDFID,UWNDID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,UWNDID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,UWNDID,'long_name',NCCHAR,12,
     *        'EW Wind Comp',IRET)
          Call NCAPTC(CDFID,UWNDID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,UWNDID,'standard_name',NCCHAR,13,
     *              'eastward_wind',IRET)
          Call NCAPTC(CDFID,UWNDID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = WNDSC
         Call NCAPT(CDFID,UWNDID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = WNDOF
          Call NCAPT(CDFID,UWNDID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
!   Second, NS wind component
          VWNDID = NCVDEF(CDFID,'wv',NCSHORT,3,DIMS,IRET)
          IRNG(1) = NINT((WNDMIN-WNDOF)/WNDSC)
          IRNG(2) = NINT((WNDMAX-WNDOF)/WNDSC)
          Call NCAPT(CDFID,VWNDID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,VWNDID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,VWNDID,'long_name',NCCHAR,12,
     *        'NS Wind Comp',IRET)
          Call NCAPTC(CDFID,VWNDID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,VWNDID,'standard_name',NCCHAR,14,
     *              'northward_wind',IRET)
          Call NCAPTC(CDFID,VWNDID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = WNDSC
         Call NCAPT(CDFID,VWNDID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = WNDOF
          Call NCAPT(CDFID,VWNDID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End if
!   Atmospheric Pressure NG01172008
        If ((PTPATM.EQ.'y'.OR.PTPATM.EQ.'Y').AND.WAVEDYN.NE.'DONONLY ')
     +   Then
         PATMID = NCVDEF(CDFID,'patm',NCSHORT,3,DIMS,IRET)
         IRNG(1) = NINT((PATMMIN-PATMOF)/PATMSC)
         IRNG(2) = NINT((PATMMAX-PATMOF)/PATMSC)
         Call NCAPT(CDFID,PATMID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,PATMID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,PATMID,'long_name',NCCHAR,20,
     *       'Atmospheric pressure',IRET)
         Call NCAPTC(CDFID,PATMID,'units',NCCHAR,4,'mbar',IRET)
         Call NCAPTC(CDFID,PATMID,'standard_name',NCCHAR,25,
     *              'air_pressure_at_sea_level',IRET)
         Call NCAPTC(CDFID,PATMID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
         FLOATVAL(1) = PATMSC
         Call NCAPT(CDFID,PATMID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
         FLOATVAL(1) = PATMOF
         Call NCAPT(CDFID,PATMID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End if
!   Rainfall rate NG09212011
       If ((PTQPREC.EQ.'y'.OR.PTQPREC.EQ.'Y').AND.WAVEDYN.NE.'DONONLY ')
     +   Then
         QPRECID = NCVDEF(CDFID,'qprec',NCSHORT,3,DIMS,IRET)
         IRNG(1) = NINT((QPRECMIN-QPRECOF)/QPRECSC)
         IRNG(2) = NINT((QPRECMAX-QPRECOF)/QPRECSC)
         Call NCAPT(CDFID,QPRECID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,QPRECID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,QPRECID,'long_name',NCCHAR,13,
     *       'Rainfall rate',IRET)
         Call NCAPTC(CDFID,QPRECID,'units',NCCHAR,7,'mm hr-1',IRET)
         Call NCAPTC(CDFID,QPRECID,'standard_name',NCCHAR,13,
     *             'rainfall_rate',IRET)
         Call NCAPTC(CDFID,QPRECID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
         FLOATVAL(1) = QPRECSC
        Call NCAPT(CDFID,QPRECID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
         FLOATVAL(1) = QPRECOF
         Call NCAPT(CDFID,QPRECID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End if
!   Other 2D Met variables NG04072008
        If((PT2DMET.EQ.'y'.OR.PT2DMET.EQ.'Y').AND.WAVEDYN.NE.'DONONLY ')
     +   Then
         AIRTID = NCVDEF(CDFID,'airt',NCSHORT,3,DIMS,IRET)
         IRNG(1) = NINT((AIRTMIN-AIRTOF)/AIRTSC)
         IRNG(2) = NINT((AIRTMAX-AIRTOF)/AIRTSC)
         Call NCAPT(CDFID,AIRTID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,AIRTID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,AIRTID,'long_name',NCCHAR,32,
     *       '2m-above-surface Air Temperature',IRET)
         Call NCAPTC(CDFID,AIRTID,'units',NCCHAR,7,'Celsius',IRET)
         Call NCAPTC(CDFID,AIRTID,'standard_name',NCCHAR,19,
     *              'surface_temperature',IRET)
         Call NCAPTC(CDFID,AIRTID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
         FLOATVAL(1) = AIRTSC
         Call NCAPT(CDFID,AIRTID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
         FLOATVAL(1) = AIRTOF
         Call NCAPT(CDFID,AIRTID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
         CLDID = NCVDEF(CDFID,'cld',NCSHORT,3,DIMS,IRET)
         IRNG(1) = NINT((CLDMIN-CLDOF)/CLDSC)
         IRNG(2) = NINT((CLDMAX-CLDOF)/CLDSC)
         Call NCAPT(CDFID,CLDID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,CLDID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,CLDID,'long_name',NCCHAR,20,
     *       'Cloud Cover fraction',IRET)
         Call NCAPTC(CDFID,CLDID,'units',NCCHAR,1,'1',IRET)
         Call NCAPTC(CDFID,CLDID,'standard_name',NCCHAR,19,
     *              'cloud_area_fraction',IRET)
         Call NCAPTC(CDFID,CLDID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
         FLOATVAL(1) = CLDSC
         Call NCAPT(CDFID,CLDID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
         FLOATVAL(1) = CLDOF
         Call NCAPT(CDFID,CLDID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
         RHUMID = NCVDEF(CDFID,'rhum',NCSHORT,3,DIMS,IRET)
         IRNG(1) = NINT((RHUMMIN-RHUMOF)/RHUMSC)
         IRNG(2) = NINT((RHUMMAX-RHUMOF)/RHUMSC)
         Call NCAPT(CDFID,RHUMID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,RHUMID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,RHUMID,'long_name',NCCHAR,34,
     *       '2m-above-surface Relative Humidity',IRET)
         Call NCAPTC(CDFID,RHUMID,'units',NCCHAR,7,'percent',IRET)
         Call NCAPTC(CDFID,RHUMID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
         FLOATVAL(1) = RHUMSC
         Call NCAPT(CDFID,RHUMID,'scale_factor',
     +              NCFLOAT,1,FLOATVAL,IRET)
         FLOATVAL(1) = RHUMOF
         Call NCAPT(CDFID,RHUMID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
         SWOBSID = NCVDEF(CDFID,'swobs',NCSHORT,3,DIMS,IRET)
         IRNG(1) = NINT((SWOBSMIN-SWOBSOF)/SWOBSSC)
         IRNG(2) = NINT((SWOBSMAX-SWOBSOF)/SWOBSSC)
         Call NCAPT(CDFID,SWOBSID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,SWOBSID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,SWOBSID,'long_name',NCCHAR,27,
     *       'Shortwave Radiation in W/m2',IRET)
         Call NCAPTC(CDFID,SWOBSID,'units',NCCHAR,5,'W m-2',IRET)
         Call NCAPTC(CDFID,SWOBSID,'standard_name',NCCHAR,41,
     *               'surface_downwelling_shortwave_flux_in_air',IRET)
         Call NCAPTC(CDFID,SWOBSID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
         FLOATVAL(1) = SWOBSSC
         Call NCAPT(CDFID,SWOBSID,'scale_factor',
     +              NCFLOAT,1,FLOATVAL,IRET)
         FLOATVAL(1) = SWOBSOF
         Call NCAPT(CDFID,SWOBSID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End if
!   Wave parameters
        If (PTWAVE.EQ.'Y'.OR.PTWAVE.EQ.'y') Then
!    First, Wave Heights
          WHID = NCVDEF(CDFID,'wh',NCSHORT,3,DIMS,IRET)
          IRNG(1) = NINT((WHMIN-WHOF)/WHSC)
          IRNG(2) = NINT((WHMAX-WHOF)/WHSC)
          Call NCAPT(CDFID,WHID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,WHID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,WHID,'long_name',NCCHAR,12,
     *        'wave heights',IRET)
          Call NCAPTC(CDFID,WHID,'units',NCCHAR,6,'meters',IRET)
          Call NCAPTC(CDFID,WHID,'standard_name',NCCHAR,40,
     *               'sea_surface_wind_wave_significant_height',IRET)
          Call NCAPTC(CDFID,WHID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = WHSC
          Call NCAPT(CDFID,WHID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = WHOF
          Call NCAPT(CDFID,WHID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
!    Second, Wave Periods
          WPID = NCVDEF(CDFID,'wp',NCSHORT,3,DIMS,IRET)
          IRNG(1) = NINT((WPMIN-WPOF)/WPSC)
          IRNG(2) = NINT((WPMAX-WPOF)/WPSC)
          Call NCAPT(CDFID,WPID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,WPID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,WPID,'long_name',NCCHAR,12,
     *        'wave periods',IRET)
          Call NCAPTC(CDFID,WPID,'units',NCCHAR,3,'sec',IRET)
          Call NCAPTC(CDFID,WPID,'standard_name',NCCHAR,28,
     *               'sea_surface_wind_wave_period',IRET)
          Call NCAPTC(CDFID,WPID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = WPSC
          Call NCAPT(CDFID,WPID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = WPOF
          Call NCAPT(CDFID,WPID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
!    Third, Wave Direction
          WDID = NCVDEF(CDFID,'wd',NCSHORT,3,DIMS,IRET)
          IRNG(1) = NINT((WDMIN-WDOF)/WDSC)
          IRNG(2) = NINT((WDMAX-WDOF)/WDSC)
          Call NCAPT(CDFID,WDID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,WDID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,WDID,'long_name',NCCHAR,54,
     *    'Significant Wave Direction in Oceanographic Convention',IRET)
          Call NCAPTC(CDFID,WDID,'units',NCCHAR,13,'degrees_north',IRET)
          Call NCAPTC(CDFID,WDID,'standard_name',NCCHAR,34,
     *               'sea_surface_wind_wave_to_direction',IRET)
          Call NCAPTC(CDFID,WDID,'coordinates',NCCHAR,12,
     *              'lat lon time',IRET)
          FLOATVAL(1) = WDSC
          Call NCAPT(CDFID,WDID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = WDOF
          Call NCAPT(CDFID,WDID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
*   4D Vars
        
!NG02222007
        IF(WAVEDYN.NE.'DONONLY ')THEN
      
        DIMS(4) = TIMEDIM
        DIMS(3) = ZPOSDIM

!NG03312009 The following are IDV specific for depth mapping
!NG03312009 IDV does not recognize the CF sigma.positive.up convention.
        If (PTD.EQ.'y') Then
         DEPTHTID = NCVDEF(CDFID,'depth_t',NCSHORT,4,DIMS,IRET)
         IRNG(1) = NINT((DPMIN-DPOF)/DPSC)
         IRNG(2) = NINT((DPMAX-DPOF)/DPSC)
         Call NCAPT(CDFID,DEPTHTID,'valid_range',NCSHORT,2,IRNG,IRET)
         Call NCAPT(CDFID,DEPTHTID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
         Call NCAPTC(CDFID,DEPTHTID,'long_name',NCCHAR,17,
     *       'vertical_position',IRET)
         Call NCAPTC(CDFID,DEPTHTID,'units',NCCHAR,6,'meters',IRET)
         Call NCAPTC(CDFID,DEPTHTID,'positive',NCCHAR,2,
     *             'up',IRET)
         Call NCAPTC(CDFID,DEPTHTID,'coordinates',NCCHAR,18,
     *               'lat lon sigma time',IRET)
         FLOATVAL(1) = DPSC
         Call NCAPT(CDFID,DEPTHTID,'scale_factor',NCFLOAT,1
     *             ,FLOATVAL,IRET)
         FLOATVAL(1) = DPOF
         Call NCAPT(CDFID,DEPTHTID,'add_offset',NCFLOAT,1
     *             ,FLOATVAL,IRET)
         EXPLTID= NCVDEF(CDFID,'ExplicitField',NCCHAR,0,0,IRET)
         Call NCAPTC(CDFID,EXPLTID,'standard_name',NCCHAR,14,
     *               'explicit_field',IRET)
         Call NCAPTC(CDFID,EXPLTID,'existingDataField',NCCHAR,7,
     *               'depth_t',IRET)
         Call NCAPTC(CDFID,EXPLTID,'_CoordinateTransformType',NCCHAR
     *               ,8,'vertical',IRET)
         Call NCAPTC(CDFID,EXPLTID,'_CoordinateAxes',NCCHAR,5,
     *               'sigma',IRET)
        Endif
    
        If (PTU.EQ.'Y'.OR.PTU.EQ.'y') Then
!   u velocity
          UID = NCVDEF(CDFID,'u',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((VELMIN-VELOF)/VELSC)
          IRNG(2) = NINT((VELMAX-VELOF)/VELSC)
!           Write (*,*) 'umin,umax', VELMIN, VELMAX, VRNG
          Call NCAPT(CDFID,UID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,UID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,UID,'long_name',NCCHAR,27,
     *              'eastward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,UID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,UID,'standard_name',NCCHAR,27,
     *              'eastward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,UID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = VELSC
          Call NCAPT(CDFID,UID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = VELOF
          Call NCAPT(CDFID,UID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTUF.EQ.'Y'.OR.PTUF.EQ.'y') Then
!   u flux
          UFLXID = NCVDEF(CDFID,'uflux',NCSHORT,4,DIMS,IRET)
          Call NCAPTC(CDFID,UFLXID,'type',NCCHAR,6,'scalar',IRET)
          IRNG(1) = NINT((FLXMIN-FLXOF)/FLXSC)
          IRNG(2) = NINT((FLXMAX-FLXOF)/FLXSC)
!          Write (*,*) 'flxmin,flxmax', FLXMIN, FLXMAX, VRNG
          Call NCAPT(CDFID,UFLXID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,UFLXID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,UFLXID,'long_name',NCCHAR,38
     *                ,'U1 advective transport per sigma layer',IRET)
          Call NCAPTC(CDFID,UFLXID,'units',NCCHAR,6,'m3 s-1',IRET)
          Call NCAPTC(CDFID,UFLXID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = FLXSC
          Call NCAPT(CDFID,UFLXID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET
     *        )
          FLOATVAL(1) = FLXOF
          Call NCAPT(CDFID,UFLXID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTV.EQ.'Y'.OR.PTV.EQ.'y') Then
!   v velocity
          VID = NCVDEF(CDFID,'v',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((VELMIN-VELOF)/VELSC)
          IRNG(2) = NINT((VELMAX-VELOF)/VELSC)
!          Write (*,*) 'vmin,vmax', VELMIN, VELMAX, VRNG
          Call NCAPT(CDFID,VID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,VID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,VID,'long_name',NCCHAR,28,
     *              'northward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,VID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,VID,'standard_name',NCCHAR,28,
     *              'northward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,VID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = VELSC
          Call NCAPT(CDFID,VID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = VELOF
          Call NCAPT(CDFID,VID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTVF.EQ.'Y'.OR.PTVF.EQ.'y') Then
!   v flux
          VFLXID = NCVDEF(CDFID,'vflux',NCSHORT,4,DIMS,IRET)
          Call NCAPTC(CDFID,VFLXID,'type',NCCHAR,6,'scalar',IRET)
          IRNG(1) = NINT((FLXMIN-FLXOF)/FLXSC)
          IRNG(2) = NINT((FLXMAX-FLXOF)/FLXSC)
!         Write (*,*) 'flxmin,flxmax', FLXMIN, FLXMAX, VRNG
          Call NCAPT(CDFID,VFLXID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,VFLXID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,VFLXID,'long_name',NCCHAR,38
     *               ,'V1 advective transport per sigma layer',IRET)
          Call NCAPTC(CDFID,VFLXID,'units',NCCHAR,6,'m3 s-1',IRET)
          Call NCAPTC(CDFID,VFLXID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = FLXSC
          Call NCAPT(CDFID,VFLXID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET
     *        )
          FLOATVAL(1) = FLXOF
          Call NCAPT(CDFID,VFLXID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTW.EQ.'Y'.OR.PTW.EQ.'y') Then
!   w velocity
          WID = NCVDEF(CDFID,'w',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((WTMIN-WTOF)/WTSC)
          IRNG(2) = NINT((WTMAX-WTOF)/WTSC)
!          Write (*,*) 'wmin,wmax', WTMIN, WTMAX, VRNG
          Call NCAPT(CDFID,WID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,WID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,WID,'long_name',NCCHAR,10,'W velocity',IRET)
          Call NCAPTC(CDFID,WID,'units',NCCHAR,5,'m s-1',IRET)
          Call NCAPTC(CDFID,WID,'standard_name',NCCHAR,25,
     *              'upward_sea_water_velocity',IRET)
          Call NCAPTC(CDFID,WID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = WTSC
          Call NCAPT(CDFID,WID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = WTOF
          Call NCAPT(CDFID,WID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTC.EQ.'Y'.OR.PTC.EQ.'y') Then
!   concentration
          CCID = NCVDEF(CDFID,'conc',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((CTMIN-CTOF)/CTSC)
          IRNG(2) = NINT((CTMAX-CTOF)/CTSC)
!         Write (*,*) 'cmin,cmax', CTMIN, CTMAX, VRNG
          Call NCAPT(CDFID,CCID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,CCID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,CCID,'long_name',NCCHAR,13,'Concentration',
     *        IRET)
          Call NCAPTC(CDFID,CCID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = CTSC
          Call NCAPT(CDFID,CCID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = CTOF
          Call NCAPT(CDFID,CCID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
          Call NCAPTC(CDFID,CCID,'units',NCCHAR,4,'none',IRET)
        End If
        If (PTS.EQ.'Y'.OR.PTS.EQ.'y') Then
!   salinity
          SALID = NCVDEF(CDFID,'salt',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((SALMIN-SALOF)/SALSC)
          IRNG(2) = NINT((SALMAX-SALOF)/SALSC)
!          Write (*,*) 'salmin,salmax', SALMIN, SALMAX, VRNG
          Call NCAPT(CDFID,SALID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,SALID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,SALID,'long_name',NCCHAR,15,
     *         'Salinity in psu',IRET)
          Call NCAPTC(CDFID,SALID,'standard_name',NCCHAR,18,
     *              'sea_water_salinity',IRET)
          Call NCAPTC(CDFID,SALID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = SALSC
          Call NCAPT(CDFID,SALID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = SALOF
          Call NCAPT(CDFID,SALID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
          Call NCAPTC(CDFID,SALID,'units',NCCHAR,4,'1e-3',IRET)
        End If
        If (PTT.EQ.'Y'.OR.PTT.EQ.'y') Then
!   temperature
          TMPID = NCVDEF(CDFID,'temp',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((TMPMIN-TMPOF)/TMPSC)
          IRNG(2) = NINT((TMPMAX-TMPOF)/TMPSC)
!          Write (*,*) 'tmpmin,tmpmax', TMPMIN, TMPMAX, VRNG
          Call NCAPT(CDFID,TMPID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,TMPID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,TMPID,'long_name',NCCHAR,17,
     *        'Water Temperature',IRET)
          Call NCAPTC(CDFID,TMPID,'units',NCCHAR,7,'Celsius',IRET)
          Call NCAPTC(CDFID,TMPID,'standard_name',NCCHAR,21,
     *              'sea_water_temperature',IRET)
          Call NCAPTC(CDFID,TMPID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = TMPSC
          Call NCAPT(CDFID,TMPID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = TMPOF
          Call NCAPT(CDFID,TMPID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTKH.EQ.'Y'.OR.PTKH.EQ.'y') Then
!    vertical diffusivity
          KHID = NCVDEF(CDFID,'kh',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((RKHMIN-KHOF)/KHSC)
          IRNG(2) = NINT((RKHMAX-KHOF)/KHSC)
!         Write (*,*) 'khmin,khmax', RKHMIN, RKHMAX, VRNG
          Call NCAPT(CDFID,KHID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,KHID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,KHID,'long_name',NCCHAR,20,
     *        'vertical diffusivity',IRET)
          Call NCAPTC(CDFID,KHID,'units',NCCHAR,6,'m2 s-1',IRET)
          Call NCAPTC(CDFID,KHID,'standard_name',NCCHAR,31,
     *              'ocean_vertical_heat_diffusivity',IRET)
          Call NCAPTC(CDFID,KHID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = KHSC
          Call NCAPT(CDFID,KHID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = KHOF
          Call NCAPT(CDFID,KHID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTKM.EQ.'Y'.OR.PTKM.EQ.'y') Then
!    vertical viscosity
          KMID = NCVDEF(CDFID,'km',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((RKMMIN-KMOF)/KMSC)
          IRNG(2) = NINT((RKMMAX-KMOF)/KMSC)
!         Write (*,*) 'kmmin,kmmax', RKMMIN, RKMMAX, VRNG
          Call NCAPT(CDFID,KMID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,KMID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,KMID,'long_name',NCCHAR,18,
     *        'vertical viscosity',IRET)
          Call NCAPTC(CDFID,KMID,'units',NCCHAR,6,'m2 s-1',IRET)
          Call NCAPTC(CDFID,KMID,'standard_name',NCCHAR,35,
     *              'ocean_vertical_momentum_diffusivity',IRET)
          Call NCAPTC(CDFID,KMID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = KMSC
          Call NCAPT(CDFID,KMID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = KMOF
          Call NCAPT(CDFID,KMID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTAM.EQ.'Y'.OR.PTAM.EQ.'y') Then
!    horizontal eddy viscosity
          AMID = NCVDEF(CDFID,'am',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((AMMIN-AMOF)/AMSC)
          IRNG(2) = NINT((AMMAX-AMOF)/AMSC)
!         Write (*,*) 'ammin,ammax', AMMIN, AMMAX
          Call NCAPT(CDFID,AMID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,AMID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,AMID,'long_name',NCCHAR,25,
     *        'horizontal eddy viscosity',IRET)
          Call NCAPTC(CDFID,AMID,'units',NCCHAR,6,'m2 s-1',IRET)
          Call NCAPTC(CDFID,AMID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          FLOATVAL(1) = AMSC
          Call NCAPT(CDFID,AMID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = AMOF
          Call NCAPT(CDFID,AMID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        If (PTQ2.EQ.'Y'.OR.PTQ2.EQ.'y') Then
!    turbulence kinetic energy
          Q2ID = NCVDEF(CDFID,'q2',NCSHORT,4,DIMS,IRET)
          IRNG(1) = NINT((Q2MIN-Q2OF)/Q2SC)
          IRNG(2) = NINT((Q2MAX-Q2OF)/Q2SC)
!         Write (*,*) 'q2min,q2max', Q2MIN, Q2MAX
          Call NCAPT(CDFID,Q2ID,'valid_range',NCSHORT,2,IRNG,IRET)
          Call NCAPT(CDFID,Q2ID,'_FillValue',NCSHORT,1,IFILLVAL,IRET)
          Call NCAPTC(CDFID,Q2ID,'long_name',NCCHAR,25,
     *        'turbulence kinetic energy',IRET)
          Call NCAPTC(CDFID,Q2ID,'coordinates',NCCHAR,18,
     *                'lat lon sigma time',IRET)
          Call NCAPTC(CDFID,Q2ID,'units',NCCHAR,6,'m2 s-1',IRET)
          FLOATVAL(1) = Q2SC
          Call NCAPT(CDFID,Q2ID,'scale_factor',NCFLOAT,1,FLOATVAL,IRET)
          FLOATVAL(1) = Q2OF
          Call NCAPT(CDFID,Q2ID,'add_offset',NCFLOAT,1,FLOATVAL,IRET)
        End If
        
        ENDIF

! Global attributes
!NG03312009 CF Vs 1.4 and COARDS CONVENTIONS GLOBAL ATTRIBUTES
!NG03312009 CF stands for Climate and Forecast; Metadata Convention
        Call NCAPTC(CDFID,NCGLOBAL,'file_type',NCCHAR,9,
     *      "Full_Grid",IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'Conventions',NCCHAR,6,
     *      "CF-1.0",IRET)
!NG     *      "CF-1.0",IRET) ! FOR IDV AT SIGMA DEPTHS
!NG     *      "CF-1.4",IRET) ! FOR IDV VECTORS, AS FLAT, NO VERTICAL RECOGNIZED
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
     *      "S_ECOM_vs8",IRET)
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
        Call NCAPTC(CDFID,NCGLOBAL,'hormix_type',NCCHAR,10,HORZMIX,
     *     IRET)
        Call NCAPT(CDFID,NCGLOBAL,'hormix_value',NCFLOAT,1,HORCON,
     *     IRET)
        Call NCAPTC(CDFID,NCGLOBAL,'vrtmix_type',NCCHAR,10,VERTMIX,
     *     IRET)
        Call NCAPT(CDFID,NCGLOBAL,'vrtmix_value',NCFLOAT,1,UMOL,
     *     IRET)      
        Call NCAPT(CDFID,NCGLOBAL,'refl_coeff',NCFLOAT,1,REFL, 
     *     IRET) ! NG04092008      
        Call NCENDF(CDFID,IRET)
!=========================================================================
!         END OF DEFINE MODE
!=========================================================================
! store easting
        call RDCORLOC
!svv creating dummy arrays to reduse size from IM+1 to IM
        do j=1,JM
           do i=1,IM
              XC(i,j) = XCORNER(i,j)
              YC(i,j) = YCORNER(i,j)
              DUM_IMJM(i,j) = DUM(i,j)
              DVM_IMJM(i,j) = DVM(i,j)
!NG03312009   CF 1.7 BOUNDS
              LONBOUNDS(1,I,J)=XCORNER(I,J)
              LONBOUNDS(2,I,J)=XCORNER(I+1,J)
              LONBOUNDS(3,I,J)=XCORNER(I+1,J+1)
              LONBOUNDS(4,I,J)=XCORNER(I,J+1)
              LATBOUNDS(1,I,J)=YCORNER(I,J)
              LATBOUNDS(2,I,J)=YCORNER(I+1,J)
              LATBOUNDS(3,I,J)=YCORNER(I+1,J+1)
              LATBOUNDS(4,I,J)=YCORNER(I,J+1)
           enddo
        enddo
        CORNER(1) = 1
        CORNER(2) = 1
        CORNER(3) = 1
        EDGES(1) = 4
        EDGES(2) = IM
        EDGES(3) = JM
!       Write (*,*) 'storing vertices'
        Call NCVPT(CDFID,LNBNDID,CORNER,EDGES,LONBOUNDS,IRET)
        Call NCVPT(CDFID,LTBNDID,CORNER,EDGES,LATBOUNDS,IRET)

        CORNER(1) = 1
        CORNER(2) = 1
        EDGES(1) = IM
        EDGES(2) = JM
!       Write (*,*) 'storing corners'
        Call NCVPT(CDFID,XID,CORNER,EDGES,XC,IRET)
        Call NCVPT(CDFID,YID,CORNER,EDGES,YC,IRET)
!       Write (*,*) 'storing centers' ! NG04022008
        Call NCVPT(CDFID,XCID,CORNER,EDGES,XGRID,IRET) ! NG04022008
        Call NCVPT(CDFID,YCID,CORNER,EDGES,YGRID,IRET) ! NG04022008
! store h1
        Call NCVPT(CDFID,H1ID,CORNER,EDGES,H1,IRET)
!       Write (*,*) 'storing h1'
! store h2
        Call NCVPT(CDFID,H2ID,CORNER,EDGES,H2,IRET)
!       Write (*,*) 'storing h2'
! store depth
!       Write (*,*) 'storing depth'
        Call NCVPT(CDFID,DEPTHID,CORNER,EDGES,HOUT,IRET)
! store grid angles
!       Write (*,*) 'storing grid angles'
        Call NCVPT(CDFID,ANGID,CORNER,EDGES,ANG,IRET)
! store FSM
!       Write (*,*) 'storing FSM'
        Call NCVPT(CDFID,FSMID,CORNER,EDGES,FSM,IRET)
! store DUM and DVM NG08072008
!       Write (*,*) 'storing DUM'
        Call NCVPT(CDFID,DUMID,CORNER,EDGES,DUM_IMJM,IRET)
!       Write (*,*) 'storing DVM'
        Call NCVPT(CDFID,DVMID,CORNER,EDGES,DVM_IMJM,IRET)
! store bottom types
!       INITIALIZE THE BTYPE WITH 999 (land)
!       DO J=1,JM
!        DO I=1,IM
!         BTYPE(I,J)=999
!        END DO
!       END DO      
        
!       open(41,FILE='bot_types.txt',STATUS='OLD')
!       READ(41,*)
!50     READ(41,*,END=51)IDUM,IIJJ,IBT
!       II=IIJJ/1000
!       JJ=MOD(IIJJ,1000)
!       write(*,*)iijj,ii,jj
!       BTYPE(II,JJ)=IBT
!       GO TO 50
!51     WRITE(*,*)'Last record',IIJJ,BT
      
!        write(*,*)' storing bottom types'
!       Call NCVPT(CDFID,BTYPEID,CORNER,EDGES,BTYPE,IRET)
! store stretched vertical coordinate levels (sigma-levels)

!NG03312009 STORING CARTESIAN VARS
        DO I=1,IM
         CNTII(I)=FLOAT(I)
        ENDDO
        DO J=1,JM
         CNTJJ(J)=FLOAT(J)
        ENDDO
        CORNER(1) = 1
        EDGES(1) = IM
        Call NCVPT(CDFID,XPOSID,CORNER,EDGES,CNTII,IRET)
        CORNER(1) = 1
        EDGES(1) = JM
        Call NCVPT(CDFID,YPOSID,CORNER,EDGES,CNTJJ,IRET)
!       Write (*,*) 'storing sigma levels'
        CORNER(1) = 1
        EDGES(1) = KBot
!NG05292008        Call NCVPT(CDFID,SIGMAID,CORNER,EDGES,Z,IRET)
        Call NCVPT(CDFID,SIGMAID,CORNER,EDGES,SIGMEAN,IRET) ! FOR SIGMA AT MID-POINTS
        Call NCVPT(CDFID,ZSDVID,CORNER,EDGES,Z,IRET) ! FOR SIGMA LEVEL BOUNDS
!NG03092007        ENDIF
        
      Else
!NG        CDFID = NCOPN('gcmplt.nc',NCWRITE,IRET) ! NG03312009 changed name to .nc as per CF
        CDFID = NCOPN('gcmplt.cdf',NCWRITE,IRET) ! This is for NYHOPS, too much of a pain otherwise
        IFILL = NCSFIL(CDFID,NCNOFILL,IRET)
      End If
      ICDF = ICDF + 1
! store time
      CORNER(1) = ICDF
!      Write (*,*) 'storing time'
      TID = NCVID(CDFID,'time',IRET)
      Call NCVPT1(CDFID,TID,CORNER,TMID,IRET)

!NG04172008 All in UTC and self-consistent.
!NG      write (*,*) 'start',IYR,IMO,IDA,IHR
      CALL DATEHR(IDA,IMO,IYR,IHR,0,SHOURPC)
      CHRSPC = SHOURPC + TMID*24.
      CALL GETDATE(IDCPC,IMCPC,IYCPC,IHCPC,IMINCPC,CHRSPC)
!NG      write (*,*) 'new',IYCPC,IMCPC,IDCPC,IHCPC,IMINCPC
      DYY=DFLOAT(IYCPC)*1.D8
      DMONTH=DFLOAT(IMCPC)*1.D6
      DDAY=DFLOAT(IDCPC)*1.D4
      DHH=DFLOAT(IHCPC)*1.D2
      DMN=DFLOAT(IMINCPC)*1.D0
      DATEHOUR=DYY+DMONTH+DDAY+DHH+DMN
!NG      write (*,*) DATEHOUR
!      write(*,2000)tmiddle,lyear,lmonth,lday,ihour
!      write(*,1000)datehour,iyy,imonth,iday,ihh
 1000 format(5I15)
 2000 format(F10.5,4I15)
      DATEID = NCVID(CDFID,'date',IRET)
      Call NCVPT1(CDFID,DATEID,CORNER,DATEHOUR,IRET)

!  2D Parameters
      CORNER(1) = 1
      CORNER(2) = 1
      CORNER(3) = ICDF
      EDGES(1) = IM
      EDGES(2) = JM
      EDGES(3) = 1
!NG02222007
      IF(WAVEDYN.NE.'DONONLY ')THEN
! store elevation
      Do J = 1, JM
        Do I = 1, IM
          WS0(I,J) = NINT((ARCET(I,J)-ELOF)/ELSC)
!          IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5)) WS0(I,J)=IFILLVAL ! to short out
          IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WS0(I,J)=IFILLVAL ! to short out
        End Do
      End Do
!      Write (*,*) 'storing elevations'
      ELID = NCVID(CDFID,'elev',IRET)
      Call NCVPT(CDFID,ELID,CORNER,EDGES,WS0,IRET)
! store heatflux
      If (PTHF.EQ.'y'.OR.PTHF.EQ.'Y') Then
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCHF(I,J)*(-4.186E6)-HFOF)/HFSC)
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing heat flux'
        HFID = NCVID(CDFID,'heat_flux',IRET)
        Call NCVPT(CDFID,HFID,CORNER,EDGES,WS0,IRET)
      End If
! store drag coefficient
      If (PTCD.EQ.'y'.OR.PTCD.EQ.'Y') Then
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((CBC(I,J)*1E4-BFOF)/BFSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing drag coefficient'
        BFID = NCVID(CDFID,'cbc',IRET)
        Call NCVPT(CDFID,BFID,CORNER,EDGES,WS0,IRET)
      End If

      ENDIF

!  Wave parameters
      If (PTWAVE.EQ.'Y'.OR.PTWAVE.EQ.'y') Then
! First, Wave heights
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCWVHT(I,J)-WHOF)/WHSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF(ARCWVHT(I,J).LT.0.) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!       Write (*,*) 'storing Wave Heights'
        WHID = NCVID(CDFID,'wh',IRET)
        Call NCVPT(CDFID,WHID,CORNER,EDGES,WS0,IRET)
! Second, Wave periods
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCWVPD(I,J)-WPOF)/WPSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF(ARCWVHT(I,J).LT.0.) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!       Write (*,*) 'storing Wave periods'
        WPID = NCVID(CDFID,'wp',IRET)
        Call NCVPT(CDFID,WPID,CORNER,EDGES,WS0,IRET)
! Third, Wave direction
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCWVDR(I,J)-WDOF)/WDSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WS0(I,J)=IFILLVAL ! to short out
            IF(ARCWVHT(I,J).LT.0.) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!       Write (*,*) 'storing Wave Direction'
        WDID = NCVID(CDFID,'wd',IRET)
        Call NCVPT(CDFID,WDID,CORNER,EDGES,WS0,IRET)
      End If
!  Wind parameters
      If (PTWND.EQ.'Y'.OR.PTWND.EQ.'y') Then
! First, EW wind
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCWSX2D(I,J)-WNDOF)/WNDSC)
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing EW Winds'
        UWNDID = NCVID(CDFID,'wu',IRET)
        Call NCVPT(CDFID,UWNDID,CORNER,EDGES,WS0,IRET)
! Second, NS Wind
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCWSY2D(I,J)-WNDOF)/WNDSC)
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing NS Winds'
        VWNDID = NCVID(CDFID,'wv',IRET)
        Call NCVPT(CDFID,VWNDID,CORNER,EDGES,WS0,IRET)
!  temporary output of wind and surface solar radiation
        WRITE(1002,'(f10.3,1x,f14.0)')TMID,DATEHOUR
        WRITE(1006,'(f10.3,1x,f14.0,i10)')TMID,DATEHOUR,NUMEPC
       DO J=1,JM
        DO I=1,IM
         IF (FSM(I,J).NE.0) 
     +   write(1002,'(2i5,6f10.3)')I,J,ARCET(I,J)
     +                         ,ARCWSX2D(I,J),ARCWSY2D(I,J)
     +                         ,ARCSWOBS(I,J),VARWIF(I,J)*1000.
     +                         ,AREAICE(I,J)
        END DO
       END DO
       DO J=39,41
        DO I=22,24
         write(1006,'(2i5,13E10.3)')I,J,ARCET(I,J)
     +                    ,ARCWSX2D(I,J),ARCWSY2D(I,J),ARCPATM(I,J)
     +                    ,QPREC2D(I,J),QEVAP2D(I,J)
     +                    ,DUMWAD(I,J),DVMWAD(I,J)
     +                    ,DUM(I,J),DVM(I,J)
     +                    ,AAM2D(I,J),AAM(I,J,2),KM(I,J,2)
        END DO
       END DO
      End if
!  Atmospheric Pressure NG01172008
      If ((PTPATM.EQ.'Y'.OR.PTPATM.EQ.'y').AND.WAVEDYN.NE.'DONONLY ')
     + Then
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCPATM(I,J)-PATMOF)/PATMSC) 
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing atmospheric pressure'
        PATMID = NCVID(CDFID,'patm',IRET)
        Call NCVPT(CDFID,PATMID,CORNER,EDGES,WS0,IRET)
      Endif
!  Rainfall Rate NG09212011
      If ((PTQPREC.EQ.'Y'.OR.PTQPREC.EQ.'y').AND.WAVEDYN.NE.'DONONLY ')
     + Then
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCQPREC(I,J)-QPRECOF)/QPRECSC) 
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing atmospheric pressure'
        QPRECID = NCVID(CDFID,'qprec',IRET)
        Call NCVPT(CDFID,QPRECID,CORNER,EDGES,WS0,IRET)
      Endif
      If ((PT2DMET.EQ.'Y'.OR.PT2DMET.EQ.'y').AND.WAVEDYN.NE.'DONONLY ')
     + Then
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCAIRT(I,J)-AIRTOF)/AIRTSC) 
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing air temperature'
        AIRTID = NCVID(CDFID,'airt',IRET)
        Call NCVPT(CDFID,AIRTID,CORNER,EDGES,WS0,IRET)
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCCLD(I,J)-CLDOF)/CLDSC) 
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing cloud cover'
        CLDID = NCVID(CDFID,'cld',IRET)
        Call NCVPT(CDFID,CLDID,CORNER,EDGES,WS0,IRET)
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCRHUM(I,J)-RHUMOF)/RHUMSC) 
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing relative humidity'
        RHUMID = NCVID(CDFID,'rhum',IRET)
        Call NCVPT(CDFID,RHUMID,CORNER,EDGES,WS0,IRET)
        Do J = 1, JM
          Do I = 1, IM
            WS0(I,J) = NINT((ARCSWOBS(I,J)-SWOBSOF)/SWOBSSC) 
            IF(FSM(I,J).EQ.0) WS0(I,J)=IFILLVAL ! to short out
          End Do
        End Do
!      Write (*,*) 'storing shortwave radiation'
        SWOBSID = NCVID(CDFID,'swobs',IRET)
        Call NCVPT(CDFID,SWOBSID,CORNER,EDGES,WS0,IRET)
      Endif

!NG02222007
      IF(WAVEDYN.NE.'DONONLY ')THEN

!  3D Parameters
      CORNER(3) = 1
      CORNER(4) = ICDF
      EDGES(3) = KBot
      EDGES(4) = 1
! store instantanous depth CNG03312009
      If (PTD.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              IF(K.EQ.KBot) THEN
!                 WA1(I,J,K) = NINT((ARCET(I,J)+Z(K)*(H(I,J)+ARCET(I,J))
               WA1(I,J,K) = NINT((ARCET(I,J)+Z(K)*(HOUT(I,J)+ARCET(I,J))
     *              -DPOF)/DPSC)
              ELSE
!                 WA1(I,J,K) = NINT((ARCET(I,J)+ZZ(K)*(H(I,J)+ARCET(I,J))
              WA1(I,J,K) = NINT((ARCET(I,J)+ZZ(K)*(HOUT(I,J)+ARCET(I,J))
     *              -DPOF)/DPSC)
              ENDIF
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!       Write (*,*) 'storing salinity'
        DEPTHTID = NCVID(CDFID,'depth_t',IRET)
        Call NCVPT(CDFID,DEPTHTID,CORNER,EDGES,WA1,IRET)
      End If
! store salinity
      If (PTS.EQ.'Y'.OR.PTS.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCS(I,J,K)-SALOF)/SALSC)
               IF(K.EQ.KBot)
     .          WA1(I,J,K) = NINT((ARCS(I,J,K-1)-SALOF)/SALSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!       Write (*,*) 'storing salinity'
        SALID = NCVID(CDFID,'salt',IRET)
        Call NCVPT(CDFID,SALID,CORNER,EDGES,WA1,IRET)
      End If
! store temperature
      If (PTT.EQ.'Y'.OR.PTT.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCT(I,J,K)-TMPOF)/TMPSC)
               IF(K.EQ.KBot)
     .          WA1(I,J,K) = NINT((ARCT(I,J,K-1)-TMPOF)/TMPSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!       Write (*,*) 'storing temperature'
        TMPID = NCVID(CDFID,'temp',IRET)
        Call NCVPT(CDFID,TMPID,CORNER,EDGES,WA1,IRET)
      End If
! store vertical diffusivity
      If (PTKH.EQ.'Y'.OR.PTKH.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCKH(I,J,K)-KHOF)/KHSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!      Write (*,*) 'storing vertical diffusivity'
        KHID = NCVID(CDFID,'kh',IRET)
        Call NCVPT(CDFID,KHID,CORNER,EDGES,WA1,IRET)
      End If
! store vertical viscosity
      If (PTKM.EQ.'Y'.OR.PTKM.EQ.'y') Then
       Do J = 1, JM
        Do K = 1, KBot
         Do I = 1, IM
              WA1(I,J,K) = NINT((ARCKM(I,J,K)-KMOF)/KMSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!      Write (*,*) 'storing vertical viscosity'
        KMID = NCVID(CDFID,'km',IRET)
        Call NCVPT(CDFID,KMID,CORNER,EDGES,WA1,IRET)
      End If
! store horizontal viscosity
      If (PTAM.EQ.'Y'.OR.PTAM.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCAM(I,J,K)-AMOF)/AMSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!      Write (*,*) 'storing horizontal viscosity'
        AMID = NCVID(CDFID,'am',IRET)
        Call NCVPT(CDFID,AMID,CORNER,EDGES,WA1,IRET)
      End If
! store turbulence kinetic energy
      If (PTQ2.EQ.'Y'.OR.PTQ2.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCQ2(I,J,K)-Q2OF)/Q2SC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
        Q2ID = NCVID(CDFID,'q2',IRET)
        Call NCVPT(CDFID,Q2ID,CORNER,EDGES,WA1,IRET)
      End If
! store concentration
      If (PTC.EQ.'Y'.OR.PTC.EQ.'y') Then
        Do J = 1, JM
          Do K = 1, KBot
            Do I = 1, IM
              WA1(I,J,K) = NINT((MAX(ARCC1(I,J,K),0.)-CTOF)/CTSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!        Write (*,*) 'storing concentration'
        CCID = NCVID(CDFID,'conc',IRET)
        Call NCVPT(CDFID,CCID,CORNER,EDGES,WA1,IRET)
      End If
! store flux
      If (PTUF.EQ.'Y'.OR.PTUF.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCUX(I,J,K)-FLXOF)/FLXSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!      Write (*,*) 'storing U component of flux'
        UFLXID = NCVID(CDFID,'uflux',IRET)
        Call NCVPT(CDFID,UFLXID,CORNER,EDGES,WA1,IRET)
      End If
      If (PTVF.EQ.'Y'.OR.PTVF.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              WA1(I,J,K) = NINT((ARCVX(I,J,K)-FLXOF)/FLXSC)
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!      Write (*,*) 'storing V component of flux'
! FIXME: order of next 2 lines should be reversed (vflux), A.Jordi
        Call NCVPT(CDFID,VFLXID,CORNER,EDGES,WA1,IRET)
        VFLXID = NCVID(CDFID,'vflux',IRET)
      End If
! store velocity
!NG03222010 CENTER AND ROTATE TO E-W, N-S FOR GODIVA
      Do J = 1, JM-1
        Do K = 1, KBot
          Do I = 1, IM-1
            UTEMP = (DUM(I,J)*ARCU(I,J,K)+
     .               DUM(I+1,J)*ARCU(I+1,J,K))
     .              /(DUM(I,J)+DUM(I+1,J)+1.E-8) 
            VTEMP = (DVM(I,J)*ARCV(I,J,K)+
     .               DVM(I,J+1)*ARCV(I,J+1,K))
     .              /(DVM(I,J)+DVM(I,J+1)+1.E-8) 
            VELU(I,J,K)=((UTEMP*COS(ANG(I,J)))
     .           -(VTEMP*SIN(ANG(I,J))))
            VELV(I,J,K)=((UTEMP*SIN(ANG(I,J)))
     .           +(VTEMP*COS(ANG(I,J))))
          End Do
        End Do
      End Do
      If (PTU.EQ.'Y'.OR.PTU.EQ.'y') Then
        WA1 = IFILLVAL
        Do K = 1, KBot
          Do J = 1, JM-1
            Do I = 1, IM-1
              WA1(I,J,K) = NINT((VELU(I,J,K)-VELOF)/VELSC) ! AT THE CENTER OF GRID 
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
            WA1(IM,J,K)=IFILLVAL ! to short out
          End Do
        End Do
        UID = NCVID(CDFID,'u',IRET)
        Call NCVPT(CDFID,UID,CORNER,EDGES,WA1,IRET)
      End If
      If (PTV.EQ.'Y'.OR.PTV.EQ.'y') Then
        WA1 = IFILLVAL
        Do K = 1, KBot
          Do J = 1, JM-1
            Do I = 1, IM-1
              WA1(I,J,K) = NINT((VELV(I,J,K)-VELOF)/VELSC) ! AT THE CENTER OF GRID 
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
! Entire outer frame initialized to fill above.
          End Do
        End Do
        VID = NCVID(CDFID,'v',IRET)
        Call NCVPT(CDFID,VID,CORNER,EDGES,WA1,IRET)
      End If
      If (PTW.EQ.'Y'.OR.PTW.EQ.'y') Then
        Do I = 1, IM
          Do J = 1, JM
            Do K = 1, KBot
              IF(K.EQ.KBot) THEN
               WA1(I,J,K) = NINT((0.-WTOF)/WTSC)
              ELSE
               WA1(I,J,K) = NINT((0.5*(ARCW(I,J,K)+ARCW(I,J,K+1))
     .                    -WTOF)/WTSC) ! AT THE CENTER OF SIGMA CELL
              ENDIF
!            IF((H(I,J)+ARCET(I,J)).LE.(WETMIN+1.E-5))WA1(I,J,K)=IFILLVAL ! to short out
            IF((HOUT(I,J)+ARCET(I,J)).LE.(1.E-5)) WA1(I,J,K)=IFILLVAL ! to short out
            End Do
          End Do
        End Do
!       Write (*,*) 'storing W component of velocity'
        WID = NCVID(CDFID,'w',IRET)
        Call NCVPT(CDFID,WID,CORNER,EDGES,WA1,IRET)
      End If

      ENDIF
      
      Call NCCLOS(CDFID,IRET)
      RETURN ! NG 11/05/02
      End
