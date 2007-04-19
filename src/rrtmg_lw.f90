!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            * 
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                       email:  miacono@aer.com                            *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

       subroutine rrtmg_lw &
            (play    ,plev    ,tlay    ,tlev    ,tsfc    ,h2ovmr  , &
             o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,cfc11vmr,cfc12vmr, &
             emis    ,inflglw ,iceflglw,liqflglw,cldfmcl ,taucmcl , &
             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc)

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! Note: The call to RRTMG_LW_INIT should be moved to the GCM initialization
!  area, since this has to be called only once. 
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    b) calls CLDPRMC to set cloud optical depth for McICA based 
!       on input cloud properties 
!    c) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation using McICA, the Monte-Carlo 
!       Independent Column Approximation, to represent sub-grid scale 
!       cloud variability
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! This call to RRTMG_LW must be preceeded by a call to the module
!     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprop.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflglw and liqflglw
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
! set of g-points for application to GCMs.  
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.  
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM3.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability. 
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrtm, only : mxlay, nbands, ngpt, maxxsec, mxmol, nlon
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave
      use rrlw_vsn

      implicit none

! ------- Declarations -------

! ----- Input -----
! Add necessary GCM include and parameter statements here (e.g. longitude, layer dimensions)

! Input arrays from GCM
      real(kind=jprb), intent(in) :: play(nlon,mxlay)           ! Layer pressures (hPa, mb)
      real(kind=jprb), intent(in) :: plev(nlon,mxlay+1)         ! Interface pressures (hPa, mb)
      real(kind=jprb), intent(in) :: tlay(nlon,mxlay)           ! Layer temperatures (K)
      real(kind=jprb), intent(in) :: tlev(nlon,mxlay+1)         ! Interface temperatures (K)
      real(kind=jprb), intent(in) :: tsfc(nlon)                 ! Surface temperature (K)
      real(kind=jprb), intent(in) :: h2ovmr(nlon,mxlay)         ! H2O volume mixing ratio
      real(kind=jprb), intent(in) :: o3vmr(nlon,mxlay)          ! O3 volume mixing ratio
      real(kind=jprb), intent(in) :: co2vmr                     ! CO2 volume mixing ratio
      real(kind=jprb), intent(in) :: ch4vmr(nlon,mxlay)         ! Methane volume mixing ratio
      real(kind=jprb), intent(in) :: n2ovmr(nlon,mxlay)         ! Nitrous oxide volume mixing ratio
      real(kind=jprb), intent(in) :: cfc11vmr(nlon,mxlay)       ! CFC11 volume mixing ratio
      real(kind=jprb), intent(in) :: cfc12vmr(nlon,mxlay)       ! CFC12 volume mixing ratio
      real(kind=jprb), intent(in) :: emis(nlon,nbands)          ! Surface emissivity

      integer(kind=jpim), intent(in) :: inflglw                 ! Flag for cloud optical properties
      integer(kind=jpim), intent(in) :: iceflglw                ! Flag for ice particle specification
      integer(kind=jpim), intent(in) :: liqflglw                ! Flag for liquid droplet specification

      real(kind=jprb), intent(in) :: cldfmcl(ngpt,nlon,mxlay)   ! Cloud fraction
      real(kind=jprb), intent(in) :: ciwpmcl(ngpt,nlon,mxlay)   ! Cloud ice water path (g/m2)
      real(kind=jprb), intent(in) :: clwpmcl(ngpt,nlon,mxlay)   ! Cloud liquid water path (g/m2)
      real(kind=jprb), intent(in) :: reicmcl(nlon,mxlay)        ! Cloud ice effective radius (microns)
      real(kind=jprb), intent(in) :: relqmcl(nlon,mxlay)        ! Cloud water drop effective radius (microns)
      real(kind=jprb), intent(in) :: taucmcl(ngpt,nlon,mxlay)   ! Cloud optical depth
!      real(kind=jprb), intent(in) :: ssacmcl(ngpt,nlon,mxlay)   ! Cloud single scattering albedo
                                                                !   for future expansion
                                                                !   lw scattering not yet available
!      real(kind=jprb), intent(in) :: asmcmcl(ngpt,nlon,mxlay)   ! Cloud asymmetry parameter
                                                                !   for future expansion
                                                                !   lw scattering not yet available
!      real(kind=jprb), intent(in) :: tauaer(nlon,mxlay,nbands)   ! aerosol optical depth
                                                                  ! for future expansion 
                                                                  !   (lw aerosols/scattering not yet available)
!      real(kind=jprb), intent(in) :: ssaaer(nlon,mxlay,nbands)   ! aerosol single scattering albedo
                                                                  ! for future expansion 
                                                                  !   (lw aerosols/scattering not yet available)
!      real(kind=jprb), intent(in) :: asmaer(nlon,mxlay,nbands)   ! aerosol asymmetry parameter
                                                                  ! for future expansion 
                                                                  !   (lw aerosols/scattering not yet available)

! ----- Output -----

      real(kind=jprb), intent(out) :: uflx(nlon,mxlay+1)        ! Total sky longwave upward flux (W/m2)
      real(kind=jprb), intent(out) :: dflx(nlon,mxlay+1)        ! Total sky longwave downward flux (W/m2)
      real(kind=jprb), intent(out) :: hr(nlon,mxlay)            ! Total sky longwave radiative heating rate (K/d)
      real(kind=jprb), intent(out) :: uflxc(nlon,mxlay+1)       ! Clear sky longwave upward flux (W/m2)
      real(kind=jprb), intent(out) :: dflxc(nlon,mxlay+1)       ! Clear sky longwave downward flux (W/m2)
      real(kind=jprb), intent(out) :: hrc(nlon,mxlay)           ! Clear sky longwave radiative heating rate (K/d)

! ----- Local -----

! Control
      integer(kind=jpim) :: nlayers             ! total number of layers
      integer(kind=jpim) :: istart              ! beginning band of calculation
      integer(kind=jpim) :: iend                ! ending band of calculation
      integer(kind=jpim) :: icld                ! clear/cloud flag
      integer(kind=jpim) :: iout                ! output option flag (inactive)
      integer(kind=jpim) :: iplon               ! column loop index
      integer(kind=jpim) :: imca                ! flag for mcica [0=off, 1=on]
      integer(kind=jpim) :: ims                 ! value for changing mcica permute seed
      integer(kind=jpim) :: k                   ! layer loop index

! Atmosphere
      real(kind=jprb) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=jprb) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=jprb) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=jprb) :: tbound                 ! surface temperature (K)
      real(kind=jprb) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=jprb) :: coldry(mxlay)          ! dry air column density (mol/cm2)
      real(kind=jprb) :: wbrodl(mxlay)          ! broadening gas column density (mol/cm2)
      real(kind=jprb) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)
      real(kind=jprb) :: wx(maxxsec,mxlay)      ! cross-section amounts (mol/cm-2)
      real(kind=jprb) :: pwvcm                  ! precipitable water vapor (cm)
      real(kind=jprb) :: semiss(nbands)         ! lw surface emissivity
      real(kind=jprb) :: fracs(mxlay,ngpt)      ! 
      real(kind=jprb) :: taug(mxlay,ngpt)       ! gaseous optical depths

!      real(kind=jprb) :: taua(mxlay,nbands)    ! aerosol optical depth
                                                ! for future expansion 
                                                !   (lw aerosols/scattering not yet available)
!      real(kind=jprb) :: ssaa(mxlay,nbands)    ! aerosol single scattering albedo
                                                ! for future expansion 
                                                !   (lw aerosols/scattering not yet available)
!      real(kind=jprb) :: asma(mxlay,nbands)    ! aerosol asymmetry parameter
                                                ! for future expansion 
                                                !   (lw aerosols/scattering not yet available)

! Atmosphere - setcoef
      integer(kind=jpim) :: laytrop            ! tropopause layer index
      integer(kind=jpim) :: jp(mxlay)          ! lookup table index 
      integer(kind=jpim) :: jt(mxlay)          ! lookup table index 
      integer(kind=jpim) :: jt1(mxlay)         ! lookup table index 
      real(kind=jprb) :: planklay(mxlay,nbands) ! 
      real(kind=jprb) :: planklev(0:mxlay,nbands) ! 
      real(kind=jprb) :: plankbnd(nbands)       ! 

      real(kind=jprb) :: colh2o(mxlay)         ! column amount (h2o)
      real(kind=jprb) :: colco2(mxlay)         ! column amount (co2)
      real(kind=jprb) :: colo3(mxlay)          ! column amount (o3)
      real(kind=jprb) :: coln2o(mxlay)         ! column amount (n2o)
      real(kind=jprb) :: colco(mxlay)          ! column amount (co)
      real(kind=jprb) :: colch4(mxlay)         ! column amount (ch4)
      real(kind=jprb) :: colo2(mxlay)          ! column amount (o2)
      real(kind=jprb) :: colbrd(mxlay)         ! column amount (broadening gases)

      integer(kind=jpim) :: indself(mxlay)
      integer(kind=jpim) :: indfor(mxlay)
      real(kind=jprb) :: selffac(mxlay)
      real(kind=jprb) :: selffrac(mxlay)
      real(kind=jprb) :: forfac(mxlay)
      real(kind=jprb) :: forfrac(mxlay)

      integer(kind=jpim) :: indminor(mxlay)
      real(kind=jprb) :: minorfrac(mxlay)
      real(kind=jprb) :: scaleminor(mxlay)
      real(kind=jprb) :: scaleminorn2(mxlay)

      real(kind=jprb) :: &                     !
                         fac00(mxlay), fac01(mxlay), &
                         fac10(mxlay), fac11(mxlay) 
      real(kind=jprb) :: &                     !
                         rat_h2oco2(mxlay),rat_h2oco2_1(mxlay), &
                         rat_h2oo3(mxlay),rat_h2oo3_1(mxlay), &
                         rat_h2on2o(mxlay),rat_h2on2o_1(mxlay), &
                         rat_h2och4(mxlay),rat_h2och4_1(mxlay), &
                         rat_n2oco2(mxlay),rat_n2oco2_1(mxlay), &
                         rat_o3co2(mxlay),rat_o3co2_1(mxlay)

! Atmosphere/clouds - cldprop
      integer(kind=jpim) :: ncbands             ! number of cloud spectral bands
      integer(kind=jpim) :: inflag              ! flag for cloud property method
      integer(kind=jpim) :: iceflag             ! flag for ice cloud properties
      integer(kind=jpim) :: liqflag             ! flag for liquid cloud properties

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=jprb) :: cldfmc(ngpt,mxlay)     ! cloud fraction [mcica]
      real(kind=jprb) :: ciwpmc(ngpt,mxlay)     ! cloud ice water path [mcica]
      real(kind=jprb) :: clwpmc(ngpt,mxlay)     ! cloud liquid water path [mcica]
      real(kind=jprb) :: relqmc(mxlay)          ! liquid particle size (microns)
      real(kind=jprb) :: reicmc(mxlay)          ! ice partcle size (microns)
      real(kind=jprb) :: taucmc(ngpt,mxlay)     ! cloud optical depth [mcica]
!      real(kind=jprb) :: ssacmc(ngpt,mxlay)     ! cloud single scattering albedo [mcica]
                                                ! for future expansion 
                                                !   (lw scattering not yet available)
!      real(kind=jprb) :: asmcmc(ngpt,mxlay)     ! cloud asymmetry parameter [mcica]
                                                ! for future expansion 
                                                !   (lw scattering not yet available)

! Output
      real(kind=jprb) :: totuflux(0:mxlay)      ! upward longwave flux (w/m2)
      real(kind=jprb) :: totdflux(0:mxlay)      ! downward longwave flux (w/m2)
      real(kind=jprb) :: fnet(0:mxlay)          ! net longwave flux (w/m2)
      real(kind=jprb) :: htr(0:mxlay)           ! longwave heating rate (k/day)
      real(kind=jprb) :: totuclfl(0:mxlay)      ! clear sky upward longwave flux (w/m2)
      real(kind=jprb) :: totdclfl(0:mxlay)      ! clear sky downward longwave flux (w/m2)
      real(kind=jprb) :: fnetc(0:mxlay)         ! clear sky net longwave flux (w/m2)
      real(kind=jprb) :: htrc(0:mxlay)          ! clear sky longwave heating rate (k/day)

!
! Initializations

      oneminus = 1._jprb - 1.e-6_jprb
      pi = 2._jprb*asin(1._jprb)
      fluxfac = pi * 2.e4_jprb                ! orig:   fluxfac = pi * 2.d4  
      istart = 1
      iend = 16
      iout = 0
      ims = 1

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud overlap method  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      icld = 2 


! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption coefficient 
! data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  

      call rrtmg_lw_init

!  This is the main longitude/column loop within RRTMG.
      do iplon = 1, nlon

!  Prepare atmospheric profile from GCM for use in RRTMG, and define
!  other input parameters.  

         call inatm (iplon, icld, play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, n2ovmr, cfc11vmr, cfc12vmr, emis, &
              inflglw, iceflglw, liqflglw, &
              cldfmcl, taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
              nlayers, pavel, pz, pdp, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, relqmc)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

         call cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, ciwpmc, &
                      clwpmc, reicmc, relqmc, ncbands, taucmc)

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

         call setcoef(nlayers, istart, pavel, tavel, tz, tbound, semiss, &
                      coldry, wkl, wbrodl, &
                      laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                      colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                      colbrd, fac00, fac01, fac10, fac11, &
                      rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                      rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                      rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                      selffac, selffrac, indself, forfac, forfrac, indfor, &
                      minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

         call taumol(nlayers, pavel, wx, coldry, &
                     laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                     colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                     colbrd, fac00, fac01, fac10, fac11, &
                     rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                     rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                     rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     minorfrac, scaleminor, scaleminorn2, indminor, &
                     fracs, taug)

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, RTRNMC is called for clear and cloudy calculations.

         call rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                     cldfmc, taucmc, planklay, planklev, plankbnd, &
                     pwvcm, fracs, taug, &
                     totuflux, totdflux, fnet, htr, &
                     totuclfl, totdclfl, fnetc, htrc )

!  Transfer up and down fluxes and heating rate to output arrays.
!  Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

         do k = 0, nlayers
            uflx(iplon,k+1) = totuflux(k)
            dflx(iplon,k+1) = totdflux(k)
            uflxc(iplon,k+1) = totuclfl(k)
            dflxc(iplon,k+1) = totdclfl(k)
         enddo
         do k = 0, nlayers-1
            hr(iplon,k+1) = htr(k)
            hrc(iplon,k+1) = htrc(k)
         enddo

      enddo

      return
      end

!***************************************************************************
      subroutine lwdatinit
!***************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrtm, only : mxlay, nbands, maxxsec, maxinpx
      use rrlw_con, only: heatfac, grav, planck, boltz, &
                          clight, avogad, alosmt, gascon, radcn1, radcn2 
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave, &
                          nxmol, ixindx
      use rrlw_vsn

      implicit none
      save 
 
! Longwave spectral band limits (wavenumbers)
      wavenum1(:) = (/ 10., 350., 500., 630., 700., 820., 980.,1080., &
                     1180.,1390.,1480.,1800.,2080.,2250.,2390.,2600./)
      wavenum2(:) = (/350., 500., 630., 700., 820., 980.,1080.,1180., &
                     1390.,1480.,1800.,2080.,2250.,2390.,2600.,3250./)
      delwave(:) =  (/340., 150., 130.,  70., 120., 160., 100., 100., &
                      210.,  90., 320., 280., 170., 130., 220., 650./)

! Spectral band information
      ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
      nspa(:) = (/1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)
      nspb(:) = (/1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)

!     Heatfac is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)
      heatfac = 8.4391_jprb

!     Modified values for consistency with CAM3:
!        =  (9.80616)(86400)(1e-5)/(1.00464)
!      heatfac = 8.43339130434_jprb

!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
      nxmol = 4
      ixindx(1) = 1
      ixindx(2) = 2
      ixindx(3) = 3
      ixindx(4) = 4
      ixindx(5:maxinpx) = 0

!    Constants from NIST 01/11/2002

      grav = 9.8066_jprb
      planck = 6.62606876e-27_jprb
      boltz = 1.3806503e-16_jprb
      clight = 2.99792458e+10_jprb
      avogad = 6.02214199e+23_jprb
      alosmt = 2.6867775e+19_jprb
      gascon = 8.31447200e+07_jprb
      radcn1 = 1.191042722e-12_jprb
      radcn2 = 1.4387752_jprb
!
!     units are generally cgs
!
!     The first and second radiation constants are taken from NIST.
!     They were previously obtained from the relations:
!          radcn1 = 2.*planck*clight*clight*1.e-07
!          radcn2 = planck*clight/boltz

      return
      end

!***************************************************************************
      subroutine inatm (iplon, icld, play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, n2ovmr, cfc11vmr, cfc12vmr, emis, &
              inflglw, iceflglw, liqflglw, &
              cldfmcl, taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
              nlayers, pavel, pz, pdp, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, relqmc)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_LW.
!  Set other RRTMG_LW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrtm, only : mxlay, nbands, ngpt, nmol, maxxsec, mxmol, nlon
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave, ixindx
      use rrlw_vsn

      implicit none

! ------- Declarations -------

! ----- Input -----
! Add necessary GCM include and parameter statements here (e.g. longitude, layer dimensions)

! Input arrays from GCM
      integer(kind=jpim), intent(in) :: iplon                   ! column loop index
      integer(kind=jpim), intent(in) :: icld                    ! clear/cloud flag

      real(kind=jprb), intent(in) :: play(nlon,mxlay)           ! Layer pressures (hPa, mb)
      real(kind=jprb), intent(in) :: plev(nlon,mxlay+1)         ! Interface pressures (hPa, mb)
      real(kind=jprb), intent(in) :: tlay(nlon,mxlay)           ! Layer temperatures (K)
      real(kind=jprb), intent(in) :: tlev(nlon,mxlay+1)         ! Interface temperatures (K)
      real(kind=jprb), intent(in) :: tsfc(nlon)                 ! Surface temperature (K)
      real(kind=jprb), intent(in) :: h2ovmr(nlon,mxlay)         ! H2O volume mixing ratio
      real(kind=jprb), intent(in) :: o3vmr(nlon,mxlay)          ! O3 volume mixing ratio
      real(kind=jprb), intent(in) :: co2vmr                     ! CO2 volume mixing ratio
      real(kind=jprb), intent(in) :: ch4vmr(nlon,mxlay)         ! Methane volume mixing ratio
      real(kind=jprb), intent(in) :: n2ovmr(nlon,mxlay)         ! Nitrous oxide volume mixing ratio
      real(kind=jprb), intent(in) :: cfc11vmr(nlon,mxlay)       ! CFC11 volume mixing ratio
      real(kind=jprb), intent(in) :: cfc12vmr(nlon,mxlay)       ! CFC12 volume mixing ratio
      real(kind=jprb), intent(in) :: emis(nlon,nbands)          ! Surface emissivity

      integer(kind=jpim), intent(in) :: inflglw                 ! Flag for cloud optical properties
      integer(kind=jpim), intent(in) :: iceflglw                ! Flag for ice particle specification
      integer(kind=jpim), intent(in) :: liqflglw                ! Flag for liquid droplet specification

      real(kind=jprb), intent(in) :: cldfmcl(ngpt,nlon,mxlay)   ! Cloud fraction
      real(kind=jprb), intent(in) :: ciwpmcl(ngpt,nlon,mxlay)   ! Cloud ice water path (g/m2)
      real(kind=jprb), intent(in) :: clwpmcl(ngpt,nlon,mxlay)   ! Cloud liquid water path (g/m2)
      real(kind=jprb), intent(in) :: reicmcl(nlon,mxlay)        ! Cloud ice effective radius (microns)
      real(kind=jprb), intent(in) :: relqmcl(nlon,mxlay)        ! Cloud water drop effective radius (microns)
      real(kind=jprb), intent(in) :: taucmcl(ngpt,nlon,mxlay)   ! Cloud optical depth

! Atmosphere
      integer(kind=jpim), intent(out) :: nlayers             ! number of layers

      real(kind=jprb), intent(out) :: pavel(mxlay)           ! layer pressures (mb) 
      real(kind=jprb), intent(out) :: tavel(mxlay)           ! layer temperatures (K)
      real(kind=jprb), intent(out) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb), intent(out) :: tz(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=jprb), intent(out) :: tbound                 ! surface temperature (K)
      real(kind=jprb), intent(out) :: pdp(mxlay)             ! layer pressure thickness (hPa, mb)
      real(kind=jprb), intent(out) :: coldry(mxlay)          ! dry air column density (mol/cm2)
      real(kind=jprb), intent(out) :: wbrodl(mxlay)          ! broadening gas column density (mol/cm2)
      real(kind=jprb), intent(out) :: wkl(mxmol,mxlay)       ! molecular amounts (mol/cm-2)
      real(kind=jprb), intent(out) :: wx(maxxsec,mxlay)      ! cross-section amounts (mol/cm-2)
      real(kind=jprb), intent(out) :: pwvcm                  ! precipitable water vapor (cm)
      real(kind=jprb), intent(out) :: semiss(nbands)         ! lw surface emissivity

! Atmosphere/clouds - cldprop
      integer(kind=jpim), intent(out) :: inflag              ! flag for cloud property method
      integer(kind=jpim), intent(out) :: iceflag             ! flag for ice cloud properties
      integer(kind=jpim), intent(out) :: liqflag             ! flag for liquid cloud properties

      real(kind=jprb), intent(out) :: cldfmc(ngpt,mxlay)     ! cloud fraction [mcica]
      real(kind=jprb), intent(out) :: ciwpmc(ngpt,mxlay)     ! cloud ice water path [mcica]
      real(kind=jprb), intent(out) :: clwpmc(ngpt,mxlay)     ! cloud liquid water path [mcica]
      real(kind=jprb), intent(out) :: relqmc(mxlay)          ! liquid particle size (microns)
      real(kind=jprb), intent(out) :: reicmc(mxlay)          ! ice partcle size (microns)
      real(kind=jprb), intent(out) :: taucmc(ngpt,mxlay)     ! cloud optical depth [mcica]


! ----- Local -----
      real(kind=jprb), parameter :: amd = 28.9660_jprb       ! Effective molecular weight of dry air (g/mol)
      real(kind=jprb), parameter :: amw = 18.0160_jprb       ! Molecular weight of water vapor (g/mol)
!      real(kind=jprb), parameter :: amc = 44.0098_jprb       ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=jprb), parameter :: amo = 47.9998_jprb       ! Molecular weight of ozone (g/mol)
!      real(kind=jprb), parameter :: amo2 = 31.9999_jprb      ! Molecular weight of oxygen (g/mol)
!      real(kind=jprb), parameter :: amch4 = 16.0430_jprb     ! Molecular weight of methane (g/mol)
!      real(kind=jprb), parameter :: amn2o = 44.0128_jprb     ! Molecular weight of nitrous oxide (g/mol)
!      real(kind=jprb), parameter :: amc11 = 137.3684_jprb    ! Molecular weight of CFC11 (g/mol) - CFCL3
!      real(kind=jprb), parameter :: amc12 = 120.9138_jprb    ! Molecular weight of CFC12 (g/mol) - CF2CL2
! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=jprb), parameter :: amdw = 1.607793_jprb      ! Molecular weight of dry air / water vapor
      real(kind=jprb), parameter :: amdc = 0.658114_jprb      ! Molecular weight of dry air / carbon dioxide
      real(kind=jprb), parameter :: amdo = 0.603428_jprb      ! Molecular weight of dry air / ozone
      real(kind=jprb), parameter :: amdm = 1.805423_jprb      ! Molecular weight of dry air / methane
      real(kind=jprb), parameter :: amdn = 0.658090_jprb      ! Molecular weight of dry air / nitrous oxide
      real(kind=jprb), parameter :: amdo2 = 0.905140_jprb     ! Molecular weight of dry air / oxygen
      real(kind=jprb), parameter :: amdc1 = 0.210852_jprb     ! Molecular weight of dry air / CFC11
      real(kind=jprb), parameter :: amdc2 = 0.239546_jprb     ! Molecular weight of dry air / CFC12

      real(kind=jprb), parameter :: sbc = 5.67e-08_jprb       ! Stefan-Boltzmann constant (W/m2K4)
      real(kind=jprb), parameter :: o2mmr = 0.23143_jprb      ! o2 mass mixing ratio

      integer(kind=jpim) :: isp, l, ix, n, imol, ig           ! Loop indices
      real(kind=jprb) :: amm, amttl, wvttl, wvsh, summol      ! 

      nlayers = mxlay

!  Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
!  into RRTM arrays below.

      wkl(:,:) = 0.0_jprb
      wx(:,:) = 0.0_jprb
      cldfmc(:,:) = 0.0_jprb
      taucmc(:,:) = 0.0_jprb
      ciwpmc(:,:) = 0.0_jprb
      clwpmc(:,:) = 0.0_jprb
      reicmc(:) = 0.0_jprb
      relqmc(:) = 0.0_jprb
      amttl = 0.0_jprb
      wvttl = 0.0_jprb
 
!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_LW arrays for pressure, temperature,
!  and molecular amounts.  
!  Pressures are input in mb, or are converted to mb here.
!  Molecular amounts are input in volume mixing ratio, or are converted from 
!  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
!  here. These are then converted to molecular amount (molec/cm2) below.  
!  The dry air column COLDRY (in molec/cm2) is calculated from the level 
!  pressures, pz (in mb), based on the hydrostatic equation and includes a 
!  correction to account for h2o in the layer.  The molecular weight of moist 
!  air (amm) is calculated for each layer.  
!  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
!  assumes GCM input fields are also bottom to top. Input layer indexing
!  from GCM fields should be reversed here if necessary.

      pz(0) = plev(iplon,1)
      tz(0) = tlev(iplon,1)
      do l = 1, nlayers
         pavel(l) = play(iplon,l)
         tavel(l) = tlay(iplon,l)
         pz(l) = plev(iplon,l+1)
         tz(l) = tlev(iplon,l+1)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,l)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,l)/(1._jprb - h2o(iplon,l)))*amdw
         wkl(2,l) = co2vmr
         wkl(3,l) = o3vmr(iplon,l)
         wkl(4,l) = n2ovmr(iplon,l)
         wkl(6,l) = ch4vmr(iplon,l)
         wkl(7,l) = o2mmr*amdo2
         amm = (1._jprb - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_jprb * avogad / &
                     (1.e2_jprb * grav * amm * (1._jprb + wkl(1,l)))
      enddo

! Set cross section molecule amounts from input; convert to vmr
      do l=1, nlayers
         wx(2,l) = cfc11vmr(iplon,l)
         wx(3,l) = cfc12vmr(iplon,l)
      enddo      

! The following section can be used to set values for an additional layer (from
! the GCM top level to 0. mb) for improved calculation of TOA fluxes. 
! Temperature and molecular amounts in the extra model layer are set to 
! their values in the top GCM model layer, though these can be modified
! here if necessary. 
! If this feature is utilized, increase nlayers by one above, limit the two
! loops above to (nlayers-1), and set the top most (extra) layer values here. 

!      pavel(nlayers) = 0.5_jprb * pz(nlayers-1)
!      tavel(nlayers) = tavel(nlayers-1)
!      pz(nlayers) = 1.e-4jprb
!      tz(nlayers-1) = 0.5_jprb * (tavel(nlayers)+tavel(nlayers-1))
!      tz(nlayers) = tz(nlayers-1)
!      wkl(1,nlayers) = wkl(1,nlayers-1)
!      wkl(2,nlayers) = wkl(2,nlayers-1)
!      wkl(3,nlayers) = wkl(3,nlayers-1)
!      wkl(4,nlayers) = wkl(4,nlayers-1)
!      wkl(6,nlayers) = wkl(6,nlayers-1)
!      wkl(7,nlayers) = wkl(7,nlayers-1)
!      amm = (1._jprb - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
!      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_jprb * avogad / &
!                        (1.e2_jprb * grav * amm * (1._jprb + wkl(1,nlayers-1)))
!      wx(2,nlayers) = wx(2,nlayers-1)
!      wx(3,nlayers) = wx(3,nlayers-1)

! At this point all moleculular amounts in wkl and wx are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  also, compute precipitable
! water vapor for diffusivity angle adjustments in rtrn and rtrnmr.

      do l = 1, nlayers
         summol = 0.0_jprb
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         wbrodl(l) = coldry(l) * (1._jprb - summol)
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
         do ix = 1,maxxsec
            if (ixindx(ix) .ne. 0) then
               wx(ixindx(ix),l) = coldry(l) * wx(ix,l) * 1.e-20_jprb
            endif
         enddo
      enddo

      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_jprb * pz(0)) / (1.e2_jprb * grav)

! Set spectral surface emissivity for each longwave band.  

      do n=1,nbands
         semiss(n) = emis(iplon,n)
!          semiss(n) = 1.0
      enddo

! Transfer cloud fraction and cloud optical properties to RRTM variables,
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflglw
         iceflag = iceflglw
         liqflag = liqflglw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is in effective radius; for Fu parameterization (iceflag = 3)
! convert effective radius to generalized effective size using method of Mitchell, JAS, 2002:

         do l = 1, nlayers
            do ig = 1, ngpt
               cldfmc(ig,l) = cldfmcl(ig,iplon,l)
               taucmc(ig,l) = taucmcl(ig,iplon,l)
               ciwpmc(ig,l) = ciwpmcl(ig,iplon,l)
               clwpmc(ig,l) = clwpmcl(ig,iplon,l)
            enddo
            reicmc(l) = reicmcl(iplon,l)
            if (iceflag .eq. 3) then
               reicmc(l) = 1.5396_jprb * reicmcl(iplon,l)
            endif
            relqmc(l) = relqmcl(iplon,l)
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

!         cldfmc(:,nlayers) = 0.0_jprb
!         taucmc(:,nlayers) = 0.0_jprb
!         ciwpmc(:,nlayers) = 0.0_jprb
!         clwpmc(:,nlayers) = 0.0_jprb
!         reicmc(nlayers) = 0.0_jprb
!         relqmc(nlayers) = 0.0_jprb

      endif
      
      return
      end 










