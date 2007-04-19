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

!-----------------------------------------------------------------------------
      subroutine rtrn(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                      cldfrac, taucloud, planklay, planklev, plankbnd, &
                      pwvcm, fracs, taug, & 
                      totuflux, totdflux, fnet, htr, &
                      totuclfl, totdclfl, fnetc, htrc ) 
!-----------------------------------------------------------------------------
!
!  RRTM Longwave Radiative Transfer Model
!  Atmospheric and Environmental Research, Inc., Cambridge, MA
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with a random cloud overlap method.
!***************************************************************************

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrtm, only : mxlay, mg, nbands, ngpt
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, ngs
      use rrlw_tbl, only: tblint, bpade, tautbl, trans, tf
      use rrlw_vsn, only: hvrrtr, hnamrtr

      implicit none

! ------- Declarations -------

! Control
      integer(kind=jpim), intent(in) :: istart              ! beginning band of calculation
      integer(kind=jpim), intent(in) :: iend                ! ending band of calculation
      integer(kind=jpim), intent(in) :: iout                ! output option flag

! Atmosphere
      integer(kind=jpim), intent(in) :: nlayers             ! total number of layers
      real(kind=jprb), intent(in) :: pz(0:mxlay)            ! level (interface) pressures (hPa, mb)
      real(kind=jprb), intent(in) :: pwvcm                  ! precipitable water vapor (cm)
      real(kind=jprb), intent(in) :: semiss(nbands)         ! lw surface emissivity
      real(kind=jprb), intent(in) :: planklay(mxlay,nbands) ! 
      real(kind=jprb), intent(in) :: planklev(0:mxlay,nbands) ! 
      real(kind=jprb), intent(in) :: plankbnd(nbands)       ! 
      real(kind=jprb), intent(in) :: fracs(mxlay,ngpt)      ! 
      real(kind=jprb), intent(in) :: taug(mxlay,ngpt)       ! 

! Clouds
      integer(kind=jpim), intent(in) :: ncbands             ! number of cloud spectral bands
      real(kind=jprb), intent(in) :: cldfrac(mxlay)         ! layer cloud fraction
      real(kind=jprb), intent(in) :: taucloud(mxlay,nbands) ! layer cloud optical depth

! Output
      real(kind=jprb), intent(out) :: totuflux(0:mxlay)      ! upward longwave flux (w/m2)
      real(kind=jprb), intent(out) :: totdflux(0:mxlay)      ! downward longwave flux (w/m2)
      real(kind=jprb), intent(out) :: fnet(0:mxlay)          ! net longwave flux (w/m2)
      real(kind=jprb), intent(out) :: htr(0:mxlay)           ! longwave heating rate (k/day)
      real(kind=jprb), intent(out) :: totuclfl(0:mxlay)      ! clear sky upward longwave flux (w/m2)
      real(kind=jprb), intent(out) :: totdclfl(0:mxlay)      ! clear sky downward longwave flux (w/m2)
      real(kind=jprb), intent(out) :: fnetc(0:mxlay)         ! clear sky net longwave flux (w/m2)
      real(kind=jprb), intent(out) :: htrc(0:mxlay)          ! clear sky longwave heating rate (k/day)

! Local
! Declarations for radiative transfer
      real(kind=jprb) :: abscld(mxlay,nbands)
      real(kind=jprb) :: atot(mxlay)
      real(kind=jprb) :: atrans(mxlay)
      real(kind=jprb) :: bbugas(mxlay)
      real(kind=jprb) :: bbutot(mxlay)
      real(kind=jprb) :: clrurad(0:mxlay)
      real(kind=jprb) :: clrdrad(0:mxlay)
      real(kind=jprb) :: efclfrac(mxlay,nbands)
      real(kind=jprb) :: uflux(0:mxlay)
      real(kind=jprb) :: dflux(0:mxlay)
      real(kind=jprb) :: urad(0:mxlay)
      real(kind=jprb) :: drad(0:mxlay)
      real(kind=jprb) :: uclfl(0:mxlay)
      real(kind=jprb) :: dclfl(0:mxlay)
      real(kind=jprb) :: odcld(mxlay,nbands)


      real(kind=jprb) :: secdiff(nbands)                 ! secant of diffusivity angle
      real(kind=jprb) :: a0(nbands),a1(nbands),a2(nbands)! diffusivity angle adjustment coefficients
      real(kind=jprb) :: wtdiff, rec_6
      real(kind=jprb) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=jprb) :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real(kind=jprb) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=jprb) :: rad0, reflect, radlu, radclru

      integer(kind=jpim) :: icldlyr(mxlay)               ! flag for cloud in layer
      integer(kind=jpim) :: ibnd, ib, iband, lay, lev, l ! loop indices
      integer(kind=jpim) :: igc                          ! g-point interval counter
      integer(kind=jpim) :: iclddn                       ! flag for cloud in down path
      integer(kind=jpim) :: ittot, itgas, itr            ! lookup table indices
      integer(kind=jpim) :: ipat(16,0:2)


! ------- Definitions -------
! input
!    mxlay                        ! maximum number of model layers
!    ngpt                         ! total number of g-point subintervals
!    nbands                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    nlayers                      ! number of model layers (plev+1)
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tautbl                       ! clear sky optical depth look-up table
!    tf                           ! tau transition function look-up table
!    trans                        ! clear sky transmittance look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer

! output
!    totuflux(0:mxlay)            ! upward longwave flux (w/m2)
!    totdflux(0:mxlay)            ! downward longwave flux (w/m2)
!    fnet(0:mxlay)                ! net longwave flux (w/m2)
!    htr(0:mxlay)                 ! longwave heating rate (k/day)
!    totuclfl(0:mxlay)            ! clear sky upward longwave flux (w/m2)
!    totdclfl(0:mxlay)            ! clear sky downward longwave flux (w/m2)
!    fnetc(0:mxlay)               ! clear sky net longwave flux (w/m2)
!    htrc(0:mxlay)                ! clear sky longwave heating rate (k/day)

! These arrays indicate the spectral 'region' (used in the 
! calculation of ice cloud optical depths) corresponding
! to each spectral band.  See cldprop.f for more details.
      data ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
                 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

! This secant and weight corresponds to the standard diffusivity 
! angle.  This initial value is redefined below for some bands.
      data wtdiff /0.5_jprb/
      data rec_6 /0.166667_jprb/

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.
      data a0 / 1.66_jprb,  1.55_jprb,  1.58_jprb,  1.66_jprb, &
                1.54_jprb, 1.454_jprb,  1.89_jprb,  1.33_jprb, &
               1.668_jprb,  1.66_jprb,  1.66_jprb,  1.66_jprb, &
                1.66_jprb,  1.66_jprb,  1.66_jprb,  1.66_jprb /
      data a1 / 0.00_jprb,  0.25_jprb,  0.22_jprb,  0.00_jprb, &
                0.13_jprb, 0.446_jprb, -0.10_jprb,  0.40_jprb, &
              -0.006_jprb,  0.00_jprb,  0.00_jprb,  0.00_jprb, &
                0.00_jprb,  0.00_jprb,  0.00_jprb,  0.00_jprb /
      data a2 / 0.00_jprb, -12.0_jprb, -11.7_jprb,  0.00_jprb, &
               -0.72_jprb,-0.243_jprb,  0.19_jprb,-0.062_jprb, &
               0.414_jprb,  0.00_jprb,  0.00_jprb,  0.00_jprb, &
                0.00_jprb,  0.00_jprb,  0.00_jprb,  0.00_jprb /

      hvrrtr = '$Revision$'

      do ibnd = 1,nbands
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff(ibnd) = 1.66_jprb
         else
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm)
         endif
      enddo
      if (pwvcm.lt.1.0) secdiff(6) = 1.80_jprb
      if (pwvcm.gt.7.1) secdiff(7) = 1.50_jprb

      urad(0) = 0.0_jprb
      drad(0) = 0.0_jprb
      totuflux(0) = 0.0_jprb
      totdflux(0) = 0.0_jprb
      clrurad(0) = 0.0_jprb
      clrdrad(0) = 0.0_jprb
      totuclfl(0) = 0.0_jprb
      totdclfl(0) = 0.0_jprb

      do lay = 1, nlayers
         urad(lay) = 0.0_jprb
         drad(lay) = 0.0_jprb
         totuflux(lay) = 0.0_jprb
         totdflux(lay) = 0.0_jprb
         clrurad(lay) = 0.0_jprb
         clrdrad(lay) = 0.0_jprb
         totuclfl(lay) = 0.0_jprb
         totdclfl(lay) = 0.0_jprb

         do ib = 1, ncbands
            if (cldfrac(lay) .ge. 1.e-6_jprb) then
               odcld(lay,ib) = secdiff(ib) * taucloud(lay,ib)
               transcld = exp(-odcld(lay,ib))
               abscld(lay,ib) = 1. - transcld
               efclfrac(lay,ib) = abscld(lay,ib) * cldfrac(lay)
               icldlyr(lay) = 1
            else
               odcld(lay,ib) = 0.0_jprb
               abscld(lay,ib) = 0.0_jprb
               efclfrac(lay,ib) = 0.0_jprb
               icldlyr(lay) = 0
            endif
         enddo
      enddo

      igc = 1
! Loop over frequency bands.
      do iband = istart, iend

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1
         if (ncbands .eq. 1) then
            ib = ipat(iband,0)
         elseif (ncbands .eq.  5) then
            ib = ipat(iband,1)
         elseif (ncbands .eq. 16) then
            ib = ipat(iband,2)
         endif

! Loop over g-channels.
 1000    continue

! Radiative transfer starts here.
         radld = 0._jprb
         radclrd = 0._jprb
         iclddn = 0

! Downward radiative transfer loop.  

         do lev = nlayers, 1, -1
               plfrac = fracs(lev,igc)
               blay = planklay(lev,iband)
               dplankup = planklev(lev,iband) - blay
               dplankdn = planklev(lev-1,iband) - blay
               odepth = secdiff(iband) * taug(lev,igc)
               if (odepth .lt. 0.0_jprb) odepth = 0.0_jprb
! Cloudy layer
               if (icldlyr(lev).eq.1) then
                  iclddn = 1
                  odtot = odepth + odcld(lev,ib)
                  if (odtot .lt. 0.06_jprb) then
                     atrans(lev) = odepth - 0.5_jprb*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     atot(lev) =  odtot - 0.5_jprb*odtot*odtot
                     odtot_rec = rec_6*odtot
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(lev,ib) * (1. - atrans(lev))) + &
                         gassrc + cldfrac(lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld
                  
                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)

                  elseif (odepth .le. 0.06_jprb) then
                     atrans(lev) = odepth - 0.5_jprb*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     odtot = odepth + odcld(lev,ib)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_jprb
                     tfactot = tf(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     atot(lev) = 1._jprb - trans(ittot)

                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(lev,ib) * (1._jprb - atrans(lev))) + &
                         gassrc + cldfrac(lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld

                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)

                  else

                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5_jprb
                     odepth = tautbl(itgas)
                     atrans(lev) = 1._jprb - trans(itgas)
                     tfacgas = tf(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odcld(lev,ib)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_jprb
                     tfactot = tf(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1._jprb - trans(ittot)

                  radld = radld - radld * (atrans(lev) + &
                     efclfrac(lev,ib) * (1._jprb - atrans(lev))) + &
                     gassrc + cldfrac(lev) * &
                     (bbdtot * atot(lev) - gassrc)
                  drad(lev-1) = drad(lev-1) + radld
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  endif
! Clear layer
               else
                  if (odepth .le. 0.06_jprb) then
                     atrans(lev) = odepth-0.5_jprb*odepth*odepth
                     odepth = rec_6*odepth
                     bbd = plfrac*(blay+dplankdn*odepth)
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)
                  else
                     tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5_jprb
                     transc = trans(itr)
                     atrans(lev) = 1._jprb-transc
                     tausfac = tf(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                  endif   
                  radld = radld + (bbd-radld)*atrans(lev)
                  drad(lev-1) = drad(lev-1) + radld
               endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                  if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev) 
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd
                  else
                     radclrd = radld
                     clrdrad(lev-1) = drad(lev-1)
                  endif
            enddo

! Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.

         rad0 = fracs(1,igc) * plankbnd(iband)
!  Add in specular reflection of surface downward radiance.
         reflect = 1._jprb - semiss(iband)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd


! Upward radiative transfer loop.
         urad(0) = urad(0) + radlu
         clrurad(0) = clrurad(0) + radclru

         do lev = 1, nlayers
! Cloudy layer
            if (icldlyr(lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               radlu = radlu - radlu * (atrans(lev) + &
                   efclfrac(lev,ib) * (1._jprb - atrans(lev))) + &
                   gassrc + cldfrac(lev) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               urad(lev) = urad(lev) + radlu
! Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               urad(lev) = urad(lev) + radlu
            endif
!  Set clear sky stream to total sky stream as long as all layers
!  are clear (iclddn=0).  Streams must be calculated separately at 
!  all layers when a cloud is present (iclddn=1), because surface 
!  reflectance is different for each stream.
               if (iclddn.eq.1) then
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev) 
                  clrurad(lev) = clrurad(lev) + radclru
               else
                  radclru = radlu
                  clrurad(lev) = urad(lev)
               endif
         enddo

! Increment g-point counter
         igc = igc + 1
! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs(iband)) go to 1000

! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.
         do lev = nlayers, 0, -1
            uflux(lev) = urad(lev)*wtdiff
            dflux(lev) = drad(lev)*wtdiff
            urad(lev) = 0.0_jprb
            drad(lev) = 0.0_jprb
            totuflux(lev) = totuflux(lev) + uflux(lev) * delwave(iband)
            totdflux(lev) = totdflux(lev) + dflux(lev) * delwave(iband)
            uclfl(lev) = clrurad(lev)*wtdiff
            dclfl(lev) = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_jprb
            clrdrad(lev) = 0.0_jprb
            totuclfl(lev) = totuclfl(lev) + uclfl(lev) * delwave(iband)
            totdclfl(lev) = totdclfl(lev) + dclfl(lev) * delwave(iband)
         enddo

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(0) = totuflux(0) * fluxfac
      totdflux(0) = totdflux(0) * fluxfac
      fnet(0) = totuflux(0) - totdflux(0)
      totuclfl(0) = totuclfl(0) * fluxfac
      totdclfl(0) = totdclfl(0) * fluxfac
      fnetc(0) = totuclfl(0) - totdclfl(0)

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(lev) = totuflux(lev) * fluxfac
         totdflux(lev) = totdflux(lev) * fluxfac
         fnet(lev) = totuflux(lev) - totdflux(lev)
         totuclfl(lev) = totuclfl(lev) * fluxfac
         totdclfl(lev) = totdclfl(lev) * fluxfac
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         l = lev - 1

! Calculate heating rates at model layers
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev)) 
      enddo

! Set heating rate to zero in top layer
      htr(nlayers) = 0.0_jprb
      htrc(nlayers) = 0.0_jprb

      return
      end   

