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

      program rrtmg_lw

! -------- Description --------

! This program is the driver for RRTMG_LW, the AER LW radiation model for 
! application to GCMs, that has been adapted from RRTM_LW for improved 
! efficiency.
!
! This routine:
!    a) calls RRTMG_LW_INIT to intialize data and to perform
!       the g-point interval reduction from 256 to 140. 
!    b) calls READPROF to read in the atmospheric profile;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPROP to set cloud optical depth based on input
!       cloud properties, or CLDPRMC to set cloud optical depth
!       for McICA
!    d) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    f) calls RTRNMR (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation with a maximum-random cloud
!       overlap method, or calls RTRN to use random cloud overlap,
!       or calls RTRNMC to use McICA, the Monte-Carlo Independent 
!       Column Approximation, to represent sub-grid scale cloud 
!       variability
!    g) writes out the calculated fluxes and cooling rates
!
! Two modes of operation are possible:
!     The mode is chosen by setting flag imca below.  
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       For single column calculations, this method also requires setting flag
!       nmca below to the sample size of the Monte Carlo calculation; 
!       (nmca = 200 is recommended)
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_lw_instructions
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
      use mcica_subcol_gen_lw, only: mcica_subcol_lw

      implicit none

! ------- Declarations -------

! ----- Local -----

! Control
      integer(kind=jpim) :: nlayers             ! total number of layers
      integer(kind=jpim) :: istart              ! beginning band of calculation
      integer(kind=jpim) :: iend                ! ending band of calculation
      integer(kind=jpim) :: icld                ! clear/cloud flag
      integer(kind=jpim) :: iflag               ! control flag
      integer(kind=jpim) :: iout                ! output option flag
      integer(kind=jpim) :: ird                 ! input unit
      integer(kind=jpim) :: iwr                 ! output unit
      integer(kind=jpim) :: i                   ! output index
      integer(kind=jpim) :: iplon               ! column loop index
      integer(kind=jpim) :: ims                 ! mcica statistical loop index
      integer(kind=jpim) :: imca                ! flag for mcica [0=off, 1=on]
      integer(kind=jpim) :: nmca                ! number of mcica samples (mcica mode)
      character page 

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
      real(kind=jprb) :: taug(mxlay,ngpt)       ! 

!      real(kind=jprb) :: tauaer(mxlay,nbands)   ! aerosol optical depth
                                                ! for future expansion 
                                                !   (lw aerosols/scattering not yet available)
!      real(kind=jprb) :: ssaaer(mxlay,nbands)   ! aerosol single scattering albedo
                                                ! for future expansion 
                                                !   (lw aerosols/scattering not yet available)
!      real(kind=jprb) :: asmaer(mxlay,nbands)   ! aerosol asymmetry parameter
                                                ! for future expansion 
                                                !   (lw aerosols/scattering not yet available)

! Atmosphere - setcoef
      integer(kind=jpim) :: laytrop            ! tropopause layer index
      integer(kind=jpim) :: jp(mxlay)          ! 
      integer(kind=jpim) :: jt(mxlay)          !
      integer(kind=jpim) :: jt1(mxlay)         !
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

      real(kind=jprb) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=jprb) :: tauc(nbands,mxlay)     ! cloud optical depth (non-delta scaled)
!      real(kind=jprb) :: ssac(nbands,mxlay)     ! cloud single scattering albedo (non-delta scaled)
                                                ! for future expansion 
                                                !   (lw scattering not yet available)
!      real(kind=jprb) :: asmc(nbands,mxlay)     ! cloud asymmetry parameter (non-delta scaled)
                                                ! for future expansion 
                                                !   (lw scattering not yet available)
      real(kind=jprb) :: ciwp(mxlay)            ! cloud ice water path
      real(kind=jprb) :: clwp(mxlay)            ! cloud liquid water path
      real(kind=jprb) :: fice(mxlay)            ! cloud ice fraction
      real(kind=jprb) :: rei(mxlay)             ! cloud ice particle size
      real(kind=jprb) :: rel(mxlay)             ! cloud liquid particle size

      real(kind=jprb) :: taucloud(mxlay,nbands) ! cloud optical depth; delta scaled
!      real(kind=jprb) :: ssacloud(mxlay,nbands) ! cloud single scattering albedo; delta scaled
                                                ! for future expansion 
                                                !   (lw scattering not yet available)
!      real(kind=jprb) :: asmcloud(mxlay,nbands) ! cloud asymmetry parameter; delta scaled
                                                ! for future expansion 
                                                !   (lw scattering not yet available)

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
! Output (mean output for McICA calculation)
      real(kind=jprb) :: uflxsum(0:mxlay)       ! upward longwave flux (w/m2)
      real(kind=jprb) :: dflxsum(0:mxlay)       ! downward longwave flux (w/m2)
      real(kind=jprb) :: fnetsum(0:mxlay)       ! net longwave flux (w/m2)
      real(kind=jprb) :: htrsum(0:mxlay)        ! longwave heating rate (k/day)

!
! Initializations

      hvrrtm = '$Revision$'
      hvrini = 'NOT USED'
      hvrcld = 'NOT USED'
      hvrclc = 'NOT USED'
      hvrrtr = 'NOT USED'
      hvrrtx = 'NOT USED'
      hvrrtc = 'NOT USED'
      hvrset = 'NOT USED'
      hvrtau = 'NOT USED'
      hvrkg  = '$Revision$'
      hvratm = 'NOT USED'
      hvrutl = 'NOT USED'
      hvrext = 'NOT USED'

      hnamrtm = '  rrtmg_lw.1col.f90:'
      hnamini = '  rrtmg_lw_init.f90:'
      hnamcld = 'rrtmglw_cldprop.f90:'
      hnamclc = 'rrtmglw_cldprmc.f90:'
      hnamrtr = '  rrtmg_lw_rtrn.f90:'
      hnamrtx = 'rrtmg_lw_rtrnmr.f90:'
      hnamrtc = 'rrtmg_lw_rtrnmc.f90:'
      hnamset = 'rrtmglw_setcoef.f90:'
      hnamtau = 'rrtmg_lw_taumol.f90:'
      hnamkg  = '   rrtmg_lw_k_g.f90:'
      hnamatm = '           rrtatm.f:'
      hnamutl = '         util_xxx.f:'
      hnamext = '            extra.f:'

      oneminus = 1._jprb - 1.e-6_jprb
      pi = 2._jprb * asin(1._jprb)
      fluxfac = pi * 2.e4_jprb                ! orig:   fluxfac = pi * 2.d4  
      ird = 9
      iwr = 10
      page = char(12)

      uflxsum(0:mxlay) = 0._jprb
      dflxsum(0:mxlay) = 0._jprb
      fnetsum(0:mxlay) = 0._jprb
      htrsum(0:mxlay) = 0._jprb

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type
!   (input by subroutine readprof from input file INPUT_RRTM):
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!      imca = 0
!      imca = 1

! Set icld to select of clear or cloud calculation and cloud overlap method
!   (read by subroutine readprof from input file INPUT_RRTM):  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)


! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  

      call rrtmg_lw_init

! Open the input set of atmospheres
      open (ird,file='INPUT_RRTM',form='formatted')
! Open the output file
      open (iwr,file='OUTPUT_RRTM',form='formatted')
      
! This is the main longitude/column loop within rrtmg.

      do iplon = 1, nlon

! Input atmospheric profile from INPUT_RRTM.
         call readprof(ird, nlayers, iout, imca, icld, pdp, &
                       pavel, tavel, pz, tz, tbound, semiss, &
                       coldry, wkl, wbrodl, wx, pwvcm, &
                       inflag, iceflag, liqflag, cldfrac, &
                       tauc, fice, ciwp, clwp, rei, rel)

         istart = 1
         iend = 16
         iflag = iout

! Set nmca to sample size for Monte Carlo calculation
         if (imca.eq.0) nmca = 1
         if (imca.eq.1) nmca = 200

! Return here for multiple band output
 1000    continue
         if (iflag .gt. 0 .and. iflag .le. 40) then
            istart = iflag
            iend = iflag
         endif

! This is the statistical sampling loop for McICA

         do ims = 1, nmca

! Call sub-colum cloud generator for McICA calculations
! Output will be written for all nmca samples.  This will be excessive if
! band output (iout=99) option is selected. 

            if (imca.eq.1) then
               call mcica_subcol_lw(iplon, nlayers, icld, ims, pavel, &
                          cldfrac, ciwp, clwp, rei, rel, tauc, cldfmc, &
                          ciwpmc, clwpmc, reicmc, relqmc, taucmc)
            endif

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprop.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprop.  

!  If McICA is requested use cloud fraction and cloud physical properties 
!  generated by sub-column cloud generator above. 

            if (imca.eq.0) then
               call cldprop(nlayers, inflag, iceflag, liqflag, cldfrac, tauc, &
                            ciwp, clwp, rei, rel, ncbands, taucloud)
            else
               call cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, &
                            ciwpmc, clwpmc, reicmc, relqmc, ncbands, taucmc)
            endif

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
! For McICA, only RTRN is called for clear and cloudy calculations.

            if (imca .eq. 0) then
               if (icld .eq. 1) then
                  call rtrn(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                            cldfrac, taucloud, planklay, planklev, plankbnd, &
                            pwvcm, fracs, taug, &
                            totuflux, totdflux, fnet, htr, &
                            totuclfl, totdclfl, fnetc, htrc ) 
               else
                  call rtrnmr(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                              cldfrac, taucloud, planklay, planklev, plankbnd, &
                              pwvcm, fracs, taug, &
                              totuflux, totdflux, fnet, htr, &
                              totuclfl, totdclfl, fnetc, htrc ) 
               endif
            elseif (imca .eq. 1) then
               call rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                           cldfmc, taucmc, planklay, planklev, plankbnd, &
                           pwvcm, fracs, taug, &
                           totuflux, totdflux, fnet, htr, &
                           totuclfl, totdclfl, fnetc, htrc )
            endif

! Process output.
            if (iout .lt. 0) goto 2000

            if (imca .eq. 0) then
     
               write(iwr,9899)wavenum1(istart),wavenum2(iend),iplon
               write(iwr,9900)
               write(iwr,9901)

               do i = nlayers, 0, -1
                  if (pz(i) .lt. 1.e-2_jprb) then
                     write(iwr,9952) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  elseif (pz(i) .lt. 1.e-1) then
                     write(iwr,9953) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  elseif (pz(i) .lt. 1.) then
                     write(iwr,9954) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  elseif (pz(i) .lt. 10.) then
                     write(iwr,9955) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  elseif (pz(i) .lt. 100.) then
                     write(iwr,9956) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  elseif (pz(i) .lt. 1000.) then
                     write(iwr,9957) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  else
                     write(iwr,9958) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
                  endif
               enddo
               write(iwr,9903)page

            elseif (imca .eq. 1) then 

               do i = nlayers, 0, -1
                  uflxsum(i) = uflxsum(i) + totuflux(i)
                  dflxsum(i) = dflxsum(i) + totdflux(i)
                  fnetsum(i) = fnetsum(i) + fnet(i)
                  htrsum(i) = htrsum(i) + htr(i)
               enddo

! Output average over samples when last sample reached.  Comment this if-check
! and swap flux output write statements below to output all McICA samples.
               if (ims .eq. nmca) then

                  do i = nlayers, 0, -1
                     uflxsum(i) = uflxsum(i)/nmca
                     dflxsum(i) = dflxsum(i)/nmca
                     fnetsum(i) = fnetsum(i)/nmca
                     htrsum(i) = htrsum(i)/nmca
                  enddo

                  write(iwr,9899) wavenum1(istart),wavenum2(iend),iplon
                  write(iwr,9900)
                  write(iwr,9901)

                  do i = nlayers, 0, -1
                     if (pz(i) .lt. 1.e-2_jprb) then
                        write(iwr,9952) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     elseif (pz(i) .lt. 1.e-1) then
                        write(iwr,9953) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     elseif (pz(i) .lt. 1.) then
                        write(iwr,9954) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     elseif (pz(i) .lt. 10.) then
                        write(iwr,9955) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     elseif (pz(i) .lt. 100.) then
                        write(iwr,9956) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     elseif (pz(i) .lt. 1000.) then
                        write(iwr,9957) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     else
                        write(iwr,9958) i,pz(i),uflxsum(i),dflxsum(i),fnetsum(i),htrsum(i)
                     endif
                  enddo
                  write(iwr,9903)page

!               do i = nlayers, 0, -1
!                  if (pz(i) .lt. 1.e-2_jprb) then
!                     write(iwr,9952) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  elseif (pz(i) .lt. 1.e-1) then
!                     write(iwr,9953) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  elseif (pz(i) .lt. 1.) then
!                     write(iwr,9954) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  elseif (pz(i) .lt. 10.) then
!                     write(iwr,9955) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  elseif (pz(i) .lt. 100.) then
!                     write(iwr,9956) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  elseif (pz(i) .lt. 1000.) then
!                     write(iwr,9957) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  else
!                     write(iwr,9958) i,pz(i),totuflux(i),totdflux(i),fnet(i),htr(i)
!                  endif
!               enddo
!               write(iwr,9903)page

               endif

            endif

 2000       continue
! End statistical loop for McICA
         enddo

         if (iout .le. 40 .or. iflag .eq. 16) goto 3500
         if (iflag .eq. 99) then
            iflag = 1
         elseif (iout .eq. 99) then
            iflag = iflag + 1
         endif
         goto 1000

 3500    continue

!
! Output module version numbers
!
         write(iwr,9910) hnamrtm,hvrrtm,hnamini,hvrini,hnamcld,hvrcld, &
           hnamclc,hvrclc,hnamrtr,hvrrtr,hnamrtx,hvrrtx,hnamrtc,hvrrtc, &
           hnamset,hvrset,hnamtau,hvrtau,hnamatm,hvratm,hnamutl,hvrutl, &
           hnamext,hvrext,hnamkg,hvrkg

 4000    continue
! 4444    continue

! End longitude/column loop
      enddo

      close(ird)
      close(iwr)

 9952 format(1x,i3,9x,f7.6,3x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9953 format(1x,i3,9x,f6.5,4x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9954 format(1x,i3,8x,f6.4,5x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9955 format(1x,i3,7x,f6.3,6x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9956 format(1x,i3,6x,f6.2,7x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9957 format(1x,i3,5x,f6.1,8x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9958 format(1x,i3,5x,f6.1,8x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9899 format(1x,'Wavenumbers: ',f6.1,' - ',f6.1,' cm-1, ATM ',i6)
 9900 format(1x,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET FLUX       HEATING RATE')
 9901 format(1x,'            mb          W/m2          W/m2           W/m2          degree/day')
 9902 format(1x,i3,3x,f11.6,4x,1p,2(g12.6,2x),g13.6,3x,g16.9,0p)
 9903 format(a)
 9910 format('  Modules and versions used in this calculation:',/,/, &
              7(5x,a20,2x,a18,10x,a20,2x,a18,/))

      stop
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
      wavenum1(:) = (/ 10._jprb, 350._jprb, 500._jprb, 630._jprb, 700._jprb, 820._jprb, &
                      980._jprb,1080._jprb,1180._jprb,1390._jprb,1480._jprb,1800._jprb, &
                     2080._jprb,2250._jprb,2390._jprb,2600._jprb/)
      wavenum2(:) = (/350._jprb, 500._jprb, 630._jprb, 700._jprb, 820._jprb, 980._jprb, &
                     1080._jprb,1180._jprb,1390._jprb,1480._jprb,1800._jprb,2080._jprb, &
                     2250._jprb,2390._jprb,2600._jprb,3250._jprb/)
      delwave(:) =  (/340._jprb, 150._jprb, 130._jprb,  70._jprb, 120._jprb, 160._jprb, &
                      100._jprb, 100._jprb, 210._jprb,  90._jprb, 320._jprb, 280._jprb, &
                      170._jprb, 130._jprb, 220._jprb, 650._jprb/)

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


!************************************************************************
      subroutine readprof(ird_in, nlayers_out, iout_out, imca, icld_out, &
          pdp, pavel_out, tavel_out, pz_out, tz_out, tbound_out, semiss_out, &
          coldry_out, wkl_out, wbrodl_out, wx_out, pwvcm_out, &
          inflag_out, iceflag_out, liqflag_out, cldfrac_out, &
          tauc, fice, ciwp, clwp, rei, rel)
!************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use rrlw_con, only: pi, planck, boltz, clight, avogad, alosmt, gascon, &
                          radcn1, radcn2 

! Note: COMMON blocks are left in this routine for array passing with 
!       rrtatm.f for reading input profiles in single column mode.
!       Scalars and arrays in subroutine call are renamed to avoid conflict
!       with COMMON blocks.
!
! Purpose: Read in atmospheric profile.

! ------- Parameters -------
      parameter (mxlay=203, mxmol=38)
      parameter (nbands = 16)
      parameter (maxinpx=mxmol)
      parameter (maxxsec=4)
      parameter (maxprod = mxlay*maxxsec)
      parameter (amd = 28.9660_jprb)        ! Effective molecular weight of dry air (g/mol)
      parameter (amw = 18.0160_jprb)        ! Molecular weight of water vapor (g/mol)
      parameter (grav = 9.8066_jprb)        ! Acceleration of gravity (m/s2)

      dimension altz(0:mxlay),ixtrans(14),semis(nbands)

      common /consts/   pic,planckc,boltzc,clightc,avogadc,alosmtc,gasconc, &
                        radcn1c,radcn2c
      common /control/  numangs, iout, istart, iend, icld
      common /profile/  nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay)
      common /surface/  tbound,ireflect,semiss(nbands)
      common /species/  coldry(mxlay),wkl(mxmol,mxlay),wbrodl(mxlay),colmol(mxlay),nmol
      common /ifil/     ird,ipr,ipu,idum(15)
      common /xsecctrl/ nxmol,ixindx(maxinpx)
      common /xsec/     wx(maxxsec,mxlay)
      common /pathx/    ixmax,nxmol0,ixindx0(maxinpx),wx0(maxinpx,mxlay)    
      common /xrrtatm/  ixsect

      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay),taucloud(mxlay,nbands)

      character*80 form1(0:1),form2(0:1),form3(0:1)
      character*1 ctest, cdollar, cprcnt,cdum

! Dimensions for transfer to rrtmg
      integer(kind=jpim), intent(in) :: ird_in                  ! input file unit
      integer(kind=jpim), intent(out) :: nlayers_out             ! total number of layers
      integer(kind=jpim), intent(out) :: icld_out                ! clear/cloud flag
      integer(kind=jpim), intent(out) :: imca                    ! McICA on/off flag (1 = use McICA)
      integer(kind=jpim), intent(out) :: iout_out                ! output option flag

      real(kind=jprb), intent(out) :: pavel_out(mxlay)           ! layer pressures (mb) 
      real(kind=jprb), intent(out) :: tavel_out(mxlay)           ! layer temperatures (K)
      real(kind=jprb), intent(out) :: pz_out(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb), intent(out) :: tz_out(0:mxlay)            ! level (interface) temperatures (K)
      real(kind=jprb), intent(out) :: pdp(mxlay)                 ! layer pressure thickness (hPa, mb)
      real(kind=jprb), intent(out) :: tbound_out                 ! surface temperature (K)
      real(kind=jprb), intent(out) :: coldry_out(mxlay)          ! dry air column density (mol/cm2)
      real(kind=jprb), intent(out) :: wbrodl_out(mxlay)          ! broadening gas column density (mol/cm2)
      real(kind=jprb), intent(out) :: wkl_out(mxmol,mxlay)       ! molecular amounts (mol/cm-2)
      real(kind=jprb), intent(out) :: wx_out(maxxsec,mxlay)      ! cross-section amounts (mol/cm-2)
      real(kind=jprb), intent(out) :: pwvcm_out                  ! precipitable water vapor (cm)
      real(kind=jprb), intent(out) :: semiss_out(nbands)         ! lw surface emissivity

      integer(kind=jpim), intent(out) :: inflag_out              ! cloud property option flag
      integer(kind=jpim), intent(out) :: iceflag_out             ! ice cloud property flag
      integer(kind=jpim), intent(out) :: liqflag_out             ! liquid cloud property flag

      real(kind=jprb), intent(out) :: cldfrac_out(mxlay)         ! cloud fraction
      real(kind=jprb), intent(out) :: tauc(nbands,mxlay)         ! cloud optical depth
!      real(kind=jprb), intent(out) :: ssac(nbands,mxlay)         ! cloud single scattering albedo
                                                                 !   for future expansion
!      real(kind=jprb), intent(out) :: asmc(nbands,mxlay)         ! cloud asymmetry parameter
                                                                 !   for future expansion
      real(kind=jprb), intent(out) :: ciwp(mxlay)                ! cloud ice water path
      real(kind=jprb), intent(out) :: clwp(mxlay)                ! cloud liquid water path
      real(kind=jprb), intent(out) :: fice(mxlay)                ! cloud ice fraction
      real(kind=jprb), intent(out) :: rei(mxlay)                 ! cloud ice particle size
      real(kind=jprb), intent(out) :: rel(mxlay)                 ! cloud liquid particle size
!      real(kind=jprb), intent(out) :: tauaer_out(mxlay,nbands)   ! aerosol optical depth
                                                                 !   for future expansion
!      real(kind=jprb), intent(out) :: ssaaer_out(mxlay,nbands)   ! aerosol single scattering albedo
                                                                 !   for future expansion
!      real(kind=jprb), intent(out) :: asmaer_out(mxlay,nbands)   ! aerosol asymmetry parameter
                                                                 !   for future expansion
!

! Initializations

      data cdollar /'$'/
      data cprcnt /'%'/
      data ixtrans /0,0,0,1,2,3,0,0,0,0,0,4,0,0/
!      data wx /maxprod*0.0/

      form1(0) = '(3f10.4,a3,i2,1x,2(f7.2,f8.3,f7.2))'
      form2(0) = '(3f10.4,a3,i2,23x,(f7.2,f8.3,f7.2))'
      form3(0) = '(8e10.3)'
      form1(1) = '(g15.7,g10.4,g10.4,a3,i2,1x,2(g7.2,g8.3,g7.2))'
      form2(1) = '(g15.7,g10.4,g10.4,a3,i2,23x,(g7.2,g8.3,g7.2))'
      form3(1) = '(8g15.7)'

! Pass constants to common block names for rrtatm
      pic = pi
      planckc = planck
      boltzc = boltz
      clightc = clight
      avogadc = avogad
      alosmtc = alosmt
      gasconc = gascon
      radcn1c = radcn1
      radcn2c = radcn2 

      ixmax = maxinpx

      ird = ird_in

      do l = 1, mxlay
         do ix = 1, maxxsec
           wx(ix,l) = 0.0_jprb
         enddo
      enddo

! Top of read input loop
 1000 continue
      read (ird,9010,end=8800) ctest
      if (ctest .eq. cprcnt) goto 8900 
      if (ctest .ne. cdollar) goto 1000

      read (ird,9011) iatm, ixsect, numangs, iout, imca, icld

!  If numangs set to -1, reset to default rt code for
!  backwards compatibility with original rrtm
      if (numangs .eq. -1) numangs = 0

!  If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      if (icld .ge. 1) call readcld

!  Read in surface information.
      read (ird,9012) tbound,iemiss,ireflect,(semis(i),i=1,16)

      do iband = 1, nbands
         semiss(iband) = 1.0_jprb
         if (iemiss .eq. 1 .and. semis(1) .ne. 0._jprb) then
            semiss(iband) = semis(1)
         elseif (iemiss .eq. 2) then
            if (semis(iband) .ne. 0._jprb) then
               semiss(iband) = semis(iband)
            endif
         endif
      enddo

      if (iatm .eq. 0) then
         read (ird,9013) iform,nlayers,nmol
         if (nmol.eq.0) nmol = 7                                    
         read (ird,form1(iform)) pavel(1),tavel(1),secntk,cinp, &
              ipthak,altz(0),pz(0),tz(0),altz(1),pz(1),tz(1)
         read (ird,form3(iform)) (wkl(m,1),m=1,7), wbrodl(1)
         if(nmol .gt. 7) read (ird,form3(iform)) (wkl(m,1),m=8,nmol)

         do l = 2, nlayers
            read (ird,form2(iform)) pavel(l),tavel(l),secntk,cinp, &
                 ipthrk,altz(l),pz(l),tz(l)
            read (ird,form3(iform)) (wkl(m,l),m=1,7), wbrodl(l)
            if(nmol .gt. 7) read (ird,form3(iform)) (wkl(m,l),m=8,nmol)
         enddo
           
         if (ixsect .eq. 1) then                                 
            read (ird,9300) nxmol0
            nxmol = nxmol0
            call xsident(ird)
            read (ird,9301) iformx
     
            do l = 1, nlayers       
               read (ird,9010) cdum
               read (ird, form3(iformx)) (wx0(m,l),m=1,7),wbrodx    
               if (nxmol0 .gt. 7) read (ird,form3(iformx)) (wx0(m,l),m=8,nxmol0)
            enddo
         endif
      else
         ipu = 7
         ipr = 66
         open(unit=ipr,file='TAPE6',status='unknown')
         call rrtatm
         if (ixsect .eq. 1) then
            do mx = 1, nxmol0
               ixindx(mx) = ixtrans(ixindx0(mx))
            enddo
         endif
      endif
      if (tbound .lt. 0) tbound = tz(0)

!  Test for mixing ratio input.
      imix = 1
      do m = 1, nmol
         if (wkl(m,1) .gt. 1.0_jprb) then
            imix = 0
            goto 3600
         endif
      enddo
 3600 continue

      if (ixsect .eq. 1) then
         imixx = 0
         if (wx0(1,1) .le. 1.0_jprb) imixx = 1
      endif
      do l = 1, nlayers
         summol = 0.0_jprb
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         if (imix .eq. 1) then
            coldry(l) = wbrodl(l) / (1._jprb - summol)
            do imol = 1, nmol
               wkl(imol,l) = coldry(l) * wkl(imol,l)
            enddo
         else
            coldry(l) = wbrodl(l) + summol
         endif
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
         if (ixsect .eq. 1) then
            do ix = 1, nxmol0
               if (ixindx(ix) .ne. 0) then
                  if (imixx .eq. 1) then
                     wx(ixindx(ix),l) = coldry(l) * wx0(ix,l) * 1.e-20_jprb
                  else
                     wx(ixindx(ix),l) = wx0(ix,l) * 1.e-20_jprb
                  endif
               endif
            enddo
         endif
      enddo

!  Calculate total precipitable water 
      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_jprb * pz(0)) / (1.e2_jprb * grav)

! Pass output arrays to new variables for transfer to rrtmg through subroutine call.
      nlayers_out = nlayers
      iout_out = iout
      icld_out = icld
      tbound_out = tbound
      pwvcm_out = pwvcm
      inflag_out = inflag
      iceflag_out = iceflag
      liqflag_out = liqflag

      pz_out(0) = pz(0)
      tz_out(0) = tz(0)
      do l = 1, mxlay
         pavel_out(l) = pavel(l)
         tavel_out(l) = tavel(l)
         pz_out(l) = pz(l)
         tz_out(l) = tz(l)
         pdp(l) = (pz(l-1) - pz(l)) 
         coldry_out(l) = coldry(l)
         wbrodl_out(l) = wbrodl(l)
         cldfrac_out(l) = cldfrac(l)
         do m = 1, mxmol
            wkl_out(m,l) = wkl(m,l)
         enddo  
         do ix = 1, maxxsec
            wx_out(ix,l) = wx(ix,l)
         enddo  
      enddo
      do n = 1, nbands
         semiss_out(n) = semiss(n)
      enddo

      do l = 1, mxlay
         if (inflag.eq.0) then
            do n = 1, nbands
               tauc(n,l) = clddat1(l)
!               ssac(n,l) = 1._jprb
!               asmc(n,l) = 1._jprb
            enddo
            ciwp(l) = 0._jprb
            clwp(l) = 0._jprb
            fice(l) = 0._jprb
            rei(l) = 0._jprb
            rel(l) = 0._jprb
         else
            do n = 1, nbands
               tauc(n,l) = 0._jprb
!               ssac(n,l) = 1._jprb
!               asmc(n,l) = 1._jprb
            enddo
            cwp = clddat1(l)
            fice(l) = clddat2(l)
            ciwp(l) = cwp * fice(l)
            clwp(l) = cwp * (1. - fice(l))
            rei(l) = clddat3(l)
            rel(l) = clddat4(l)
         endif 
      enddo

!      do l = 1, mxlay
!         do n = 1, nbands
!            tauaer_out(l,n) = 0._jprb
!            ssaaer_out(l,n) = 1._jprb
!            asmaer_out(l,n) = 1._jprb
!         enddo
!      enddo

      goto 9000

 8800 continue
 8900 if (ctest.eq.'%') stop 'END OF INPUT FILE'
 9000 continue

 9010 format (a1)
 9011 format (49x,i1,19x,i1,13x,i2,2x,i3,3x,i1,i1)
 9012 format (e10.3,1x,i1,2x,i1,16e5.3)
 9013 format (1x,i1,i3,i5)                                     
 9300 format (i5)
 9301 format (1x,i1)

      return
      end 

!*************************************************************************
      subroutine readcld
!*************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 

! Purpose:  To read in IN_CLD_RRTM, the file that contains input 
!           cloud properties.

! ------- Parameters -------
      parameter (mxlay=203)
      parameter (nbands = 16)

      common /profile/   nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay)
      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay),taucloud(mxlay,nbands)

      character*1 ctest, cpercent

      data cpercent /'%'/
      irdcld = 11

      open(irdcld,file='IN_CLD_RRTM',form='formatted')

! Read in cloud input option.  
      read(irdcld,9050) inflag, iceflag, liqflag

      do lay = 1, nlayers
         cldfrac(lay) = 0._jprb
      enddo

! Top of read input loop
 1000 continue

!  For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
!  DAT1 are pertinent.  If CTEST = '%', then there are no more 
!  cloudy layers to process.
      read (irdcld,9100,end=9000) ctest,lay,frac,dat1,dat2,dat3,dat4
      if (ctest .eq. cpercent) goto 9000
      cldfrac(lay) = frac
      clddat1(lay) = dat1
      clddat2(lay) = dat2
      clddat3(lay) = dat3
      clddat4(lay) = dat4
      goto 1000

 9000 continue
      close(irdcld)

 9050 format (3x,i2,4x,i1,4x,i1)
 9100 format (a1,1x,i3,5e10.5)

      return
      end

!**********************************************************************
      subroutine xsident(ird)
!**********************************************************************

! Purpose:  This subroutine identifies which cross-sections are to be used.

! ------- Parameters -------
      parameter (maxinpx=38)
      parameter (maxxsec=4)

      common /xsecctrl/ nxmol,ixindx(maxinpx)

!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
!                                                                         
!     xsname=names, alias=aliases of the cross-section molecules          
!                                                                         
      character*10 xsname(maxinpx),alias(maxxsec,4),blank               
                                                                         
      data (alias(1,i),i=1,4)/ &
         'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/ 
      data (alias(2,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/         
      data (alias(3,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/         
      data (alias(4,i),i=1,4)/ &
          ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/        

      data blank / '          '/                                          
                                                                         
      do i = 1, nxmol                                                 
         xsname(i) = blank                                                
      enddo
                                                                         
! Read in the names of the molecules                                  
                                                                         
      if (nxmol.gt.7) then                                               
         read (ird,'(7a10)') (xsname(i),i=1,7)                            
         read (ird,'(8a10)') (xsname(i),i=8,nxmol)                       
      else                                                                
         read (ird,'(7a10)') (xsname(i),i=1,nxmol)                       
      endif                                                               

!  Match the names read in against the names stored in alias           
!  and determine the index value.  
      ixmax = 4                                                          
      do i = 1, nxmol                                                 
!  Left-justify all inputed names.                                      
         call cljust (xsname(i),10)
         ixindx(i) = 0
         do j = 1, ixmax
            if ((xsname(i).eq.alias(1,j)) .or. &
                (xsname(i).eq.alias(2,j)) .or. &
                (xsname(i).eq.alias(3,j)) .or. &
                (xsname(i).eq.alias(4,j))) then                           
               ixindx(i) = j                                              
            endif                                                         
         enddo
      enddo

      return
      end

