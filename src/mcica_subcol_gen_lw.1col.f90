!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!

      module mcica_subcol_gen_lw

!----------------------------------------------------------------------------
! Copyright (c) 2002-2020, Atmospheric & Environmental Research, Inc. (AER)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC., 
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
! THE POSSIBILITY OF SUCH DAMAGE.
!                        (http://www.rtweb.aer.com/)                        
!----------------------------------------------------------------------------

! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.
! Two options are possible:
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0; longwave scattering is not
!    yet available, ssac and asmc are for future expansion)

! --------- Modules ----------

      use parkind, only : im => kind_im, rb => kind_rb 
      use parrrtm, only : nbndlw, ngptlw
      use rrlw_con, only: grav
      use rrlw_wvn, only: ngb
      use rrlw_vsn

      implicit none

! public interfaces/functions/subroutines
      public :: get_alpha, mcica_subcol_lw, generate_stochastic_clouds 

      contains

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine get_alpha(iplon, nlayers, icld, idcor, decorr_con, &
                           dz, lat, juldat, cldfrac, alpha)

! Subroutine alpha calculates the alpha parameter required for the 
! exponential and exponential-random cloud overlap methods. Calling 
! this subroutine is therefore only required when icld = 4 or 5.  
!
! This subroutine calculates the exponential transition, alpha, from 
! maximum to random overlap required to define the fractional cloud
! vertical correlations for the exponential or exponential-random
! cloud overlap options. For exponential, the transition from maximum
! to random with distance through model layers occurs without regard
! to the configuration of clear and cloudy layers. For the exponential-
! random method, each block of adjacent cloudy layers is treated with a
! separate transition from maximum to random, and blocks of cloudy 
! layers separated by one or more clear layers are correlated randomly. 

! Inputs
      integer(kind=im), intent(in) :: iplon           ! column/longitude dimension
      integer(kind=im), intent(in) :: nlayers         ! number of model layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
                                                      !  0 = clear
                                                      !  1 = random
                                                      !  2 = maximum-random
                                                      !  3 = maximum
                                                      !  4 = exponential
                                                      !  5 = exponential-random
      integer(kind=im), intent(in) :: idcor           ! decorrelation length method
                                                      !  0 = constant
                                                      !  1 = latitude-varying (Oreopolous et al, 2012)
      integer(kind=im), intent(in) :: juldat          ! Julian day of year

      real(kind=rb), intent(in) :: decorr_con         ! decorrelation length, constant (m) 

      real(kind=rb), intent(in) :: dz(:)              ! layer thickness (m) 
                                                      ! Dimensions: (nlayers)
      real(kind=rb), intent(in) :: lat                ! latitude (degrees)

      real(kind=rb), intent(in) :: cldfrac(:)         ! layer cloud fraction
                                                      ! Dimensions: (nlayers)
! Outputs
      real(kind=rb), intent(out) :: alpha(:)          ! vertical cloud fraction correlation parameter
                                                      ! Dimensions: (nlayers)

! Local
      integer :: k
      real(kind=rb) :: decorr_lat                     ! decorrelation length, latitude-varying (m)
      real(kind=rb) :: decorr_len                     ! final decorrelation length (m)
      real(kind=rb) :: decorr_inv                     ! 1. / decorr_len

! Constants for latitude and day-of-year depenendent decorrelation length (Oreopolous et al., 2012)
! when idcor = 1
      real(kind=rb), parameter :: am1 = 1.4315
      real(kind=rb), parameter :: am2 = 2.1219
      real(kind=rb), parameter :: am4 = -25.584
      real(kind=rb), parameter :: amr = 7.0
      real(kind=rb) :: am3

      real(kind=rb), parameter :: f_one = 1.0

! If exponential or exponential-random cloud overlap is used, 
! Derive latitude-varying decorrelation length if requested; 
! otherwise use the provided constant decorrelation length, decorr_con
      decorr_inv = f_one
      if (icld .eq. 4 .or. icld .eq. 5) then
         if (idcor .eq. 1) then 
            if (juldat .gt. 181) then
               am3 = -4._rb * amr / 365._rb * (juldat - 272) 
            else
               am3 = 4._rb * amr / 365._rb * (juldat - 91)
            endif
! latitude in degrees, decorr_lat in km
            decorr_lat = am1 + am2 * exp( -(lat - am3)**2 / am4**2)
            decorr_len = decorr_lat * 1.e3_rb
         else
            decorr_len = decorr_con
         endif
         if (decorr_len .ge. 0.0_rb) decorr_inv = f_one / decorr_len
      endif

! Atmosphere data defined from sfc to toa; define alpha from sfc to toa
! Exponential cloud overlap
      if (icld .eq. 4) then 
         alpha(1) = 0.0_rb
         do k = 2, nlayers
            alpha(k) = exp( -(0.5_rb * (dz(k) + dz(k-1))) * decorr_inv)
         enddo
      endif
! Exponential-random cloud overlap
      if (icld .eq. 5) then 
         alpha(1) = 0.0_rb
         do k = 2, nlayers
            alpha(k) = exp( -(0.5_rb * (dz(k) + dz(k-1))) * decorr_inv)
      ! Decorrelate layers when a clear layer follows a cloudy layer to enforce
      ! random correlation between non-adjacent blocks of cloudy layers
            if (cldfrac(k) .eq. 0.0_rb .and. cldfrac(k-1) .gt. 0.0_rb) then
               alpha(k) = 0.0_rb
            endif
         enddo
      endif

      end subroutine get_alpha

!-------------------------------------------------------------------------------------------------
      subroutine mcica_subcol_lw(iplon, nlayers, icld, ims, irng, play, &
                       cldfrac, ciwp, clwp, rei, rel, tauc, alpha, cldfmc, &
                       ciwpmc, clwpmc, reicmc, relqmc, taucmc)

! Control
      integer(kind=im), intent(in) :: iplon           ! column/longitude dimension
      integer(kind=im), intent(in) :: nlayers         ! number of model layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(in) :: ims             ! mcica statistical loop index; also
                                                      ! value for changing mcica permute seed
      integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister

! Atmosphere
      real(kind=rb), intent(in) :: play(:)            ! layer pressures (mb) 
                                                      !    Dimensions: (nlayers)

! Atmosphere/clouds - cldprop
      real(kind=rb), intent(in) :: cldfrac(:)         ! layer cloud fraction
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: tauc(:,:)          ! in-cloud optical depth
                                                      !    Dimensions: (nbndlw,nlayers)
!      real(kind=rb), intent(in) :: ssac(:,:)         ! in-cloud single scattering albedo
                                                      !    Dimensions: (nbndlw,nlayers)
!      real(kind=rb), intent(in) :: asmc(:,:)         ! in-cloud asymmetry parameter
                                                      !    Dimensions: (nbndlw,nlayers)
      real(kind=rb), intent(in) :: ciwp(:)            ! in-cloud ice water path
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: clwp(:)            ! in-cloud liquid water path
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: rei(:)             ! cloud ice particle size
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: rel(:)             ! cloud liquid particle size
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: alpha(:)           ! vertical cloud fraction correlation parameter
                                                      ! Dimensions: (nlayers)

! ----- Output -----
! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb), intent(out) :: cldfmc(:,:)       ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(out) :: ciwpmc(:,:)       ! in-cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(out) :: clwpmc(:,:)       ! in-cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(out) :: relqmc(:)         ! liquid particle size (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(out) :: reicmc(:)         ! ice partcle size (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(out) :: taucmc(:,:)       ! in-cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
!      real(kind=rb), intent(out) :: ssacmc(:,:)      ! in-cloud single scattering albedo [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
!      real(kind=rb), intent(out) :: asmcmc(:,:)      ! in-cloud asymmetry parameter [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)

! ----- Local -----

! Stochastic cloud generator variables [mcica]
      integer(kind=im), parameter :: nsubclw = ngptlw ! number of sub-columns (g-point intervals)
      integer(kind=im) :: permuteseed                 ! if the cloud generator is called multiple times, 
                                                      ! permute the seed between each call.    
      integer(kind=im) :: ilev                        ! loop index

      real(kind=rb) :: pmid(nlayers)                  ! layer pressures (Pa) 
!      real(kind=rb) :: pdel(nlayers)                 ! layer pressure thickness (Pa) 
!      real(kind=rb) :: qi(nlayers)                   ! ice water (specific humidity)
!      real(kind=rb) :: ql(nlayers)                   ! liq water (specific humidity)


! Return if clear sky; or stop if icld out of range
      if (icld.eq.0) return
      if (icld.lt.0.or.icld.gt.5) then 
         stop 'MCICA_SUBCOL_LW: INVALID ICLD'
      endif 

! Set permute seed for random number generator
! For single column mode, permuteseed must be different for each (ims) sample performed

      permuteseed = ims*nsubclw

! Pass particle sizes to new arrays, no subcolumns for these properties yet
! Convert pressures from mb to Pa

      reicmc(:nlayers) = rei(:nlayers)
      relqmc(:nlayers) = rel(:nlayers)
      pmid(:nlayers) = play(:nlayers)*1.e2_rb

! Convert input ice and liquid cloud water paths to specific humidity ice and liquid components 

!      cwp =  (q * pdel * 1000.) / gravit)
!           = (kg/kg * kg m-1 s-2 *1000.) / m s-2
!           = (g m-2)
!
!      q  = (cwp * gravit) / (pdel *1000.)
!         = (g m-2 * m s-2) / (kg m-1 s-2 * 1000.)
!         =  kg/kg

!      do ilev = 1, nlayers
!         qi(ilev) = (ciwp(ilev) * grav) / (pdel(ilev) * 1000._rb)
!         ql(ilev) = (clwp(ilev) * grav) / (pdel(ilev) * 1000._rb)
!      enddo

!  Generate the stochastic subcolumns of cloud optical properties for the longwave;
      call generate_stochastic_clouds (nlayers, icld, irng, pmid, cldfrac, clwp, ciwp, &
                                   alpha, tauc, &
                                   cldfmc, clwpmc, ciwpmc, taucmc, permuteseed)

      end subroutine mcica_subcol_lw


!-------------------------------------------------------------------------------------------------
      subroutine generate_stochastic_clouds(nlayers, icld, irng, pmid, cld, clwp, ciwp, &
                                   alpha, tauc, &
                                   cld_stoch, clwp_stoch, ciwp_stoch, tauc_stoch, changeSeed) 
!-------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------------------
  ! ---------------------
  ! Contact: Cecile Hannay (hannay@ucar.edu)
  ! 
  ! Original code: Based on Raisanen et al., QJRMS, 2004.
  ! 
  ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
  !   random number generator, which can be changed to the optional kissvec random number generator
  !   with flag 'irng'. Some extra functionality has been commented or removed.  
  !   Michael J. Iacono, AER, Inc., February 2007
  !
  ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
  ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
  ! and uniform cloud liquid and cloud ice concentration.
  ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
  ! and obeys an overlap assumption in the vertical.   
  ! 
  ! Overlap assumption:
  !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
  !  The default option is maximum-random (option 3)
  !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
  !  This is set with the variable "overlap" 
  !mji - Exponential overlap option (overlap=4) has been deactivated in this version
  !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
  ! 
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are different.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  !
  ! PDF assumption:
  !  We can use arbitrary complicated PDFS. 
  !  In the present version, we produce homogeneuous clouds (the simplest case).  
  !  Future developments include using the PDF scheme of Ben Johnson. 
  !
  ! History file:
  !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
  !  nsubcol = number of subcolumns
  !  overlap = overlap type (1-3)
  !  Zo = length scale 
  !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
  !  CLDLIQ_S = mean of the subcolumn cloud water
  !  CLDICE_S = mean of the subcolumn cloud ice 
  !
  ! Note:
  !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
  !   i.e we only have cloud condensate when the cell is cloudy. 
  !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
  !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
  !   without cloud condensate or the opposite).
  !---------------------------------------------------------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      use MersenneTwister, only: randomNumberSequence, &   
                                 new_RandomNumberSequence, getRandomReal

      type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer(kind=im), intent(in) :: nlayers         ! number of layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister
      integer(kind=im), optional, intent(in) :: changeSeed     ! allows permuting the seed

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real(kind=rb), intent(in) :: pmid(:)            ! layer pressure (Pa)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: cld(:)             ! cloud fraction 
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: clwp(:)            ! in-cloud liquid water path
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: ciwp(:)            ! in-cloud ice water path
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: tauc(:,:)          ! in-cloud optical depth
                                                      !    Dimensions: (nbndlw,nlayers)
!      real(kind=rb), intent(in) :: ssac(:,:)         ! in-cloud single scattering albedo
                                                      !    Dimensions: (nbndlw,nlayers)
                                                      !   inactive - for future expansion
!      real(kind=rb), intent(in) :: asmc(:,:)         ! in-cloud asymmetry parameter
                                                      !    Dimensions: (nbndlw,nlayers)
                                                      !   inactive - for future expansion
      real(kind=rb), intent(in) :: alpha(:)           ! vertical cloud fraction correlation parameter
                                                      !    Dimensions: (nlayers)

      real(kind=rb), intent(out) :: cld_stoch(:,:)    ! subcolumn cloud fraction 
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(out) :: clwp_stoch(:,:)   ! subcolumn in-cloud liquid water path
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(out) :: ciwp_stoch(:,:)   ! subcolumn in-cloud ice water path
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(out) :: tauc_stoch(:,:)   ! subcolumn in-cloud optical depth
                                                      !    Dimensions: (ngptlw,nlayers)
!      real(kind=rb), intent(out) :: ssac_stoch(:,:)  ! subcolumn in-cloud single scattering albedo
                                                      !    Dimensions: (ngptlw,nlayers)
                                                      !   inactive - for future expansion
!      real(kind=rb), intent(out) :: asmc_stoch(:,:)  ! subcolumn in-cloud asymmetry parameter
                                                      !    Dimensions: (ngptlw,nlayers)
                                                      !   inactive - for future expansion

! -- Local variables
      integer(kind=im), parameter :: nsubcol = ngptlw ! number of sub-columns (g-point intervals)
      real(kind=rb) :: cldf(nlayers)                  ! cloud fraction 

! Mean over the subcolumns (cloud fraction, cloud water , cloud ice) - inactive
!      real(kind=rb) :: mean_cld_stoch(nlayers)       ! cloud fraction 
!      real(kind=rb) :: mean_clwp_stoch(nlayers)      ! cloud water
!      real(kind=rb) :: mean_ciwp_stoch(nlayers)      ! cloud ice
!      real(kind=rb) :: mean_tauc_stoch(nlayers)      ! cloud optical depth
!      real(kind=rb) :: mean_ssac_stoch(nlayers)      ! cloud single scattering albedo
!      real(kind=rb) :: mean_asmc_stoch(nlayers)      ! cloud asymmetry parameter

! Set overlap
      integer(kind=im) :: overlap                     ! 1 = random overlap, 2 = maximum_random,
                                                      ! 3 = maximum overlap, 4 = exponential overlap
                                                      ! 5 = exponential_random overlap

! Constants (min value for cloud fraction and cloud water and ice)
      real(kind=rb), parameter :: cldmin = 1.0e-20_rb ! min cloud fraction
!      real(kind=rb), parameter :: qmin   = 1.0e-10_rb   ! min cloud water and cloud ice (not used)

! Variables related to random number and seed 
      real(kind=rb), dimension(nsubcol, nlayers) :: CDF, CDF2 ! random numbers
      integer(kind=im) :: seed1, seed2, seed3, seed4          ! seed to create random number (kissvec)
      real(kind=rb) :: rand_num                       ! random number (kissvec)
      integer(kind=im) :: iseed                       ! seed to create random number (Mersenne Twister)
      real(kind=rb) :: rand_num_mt                    ! random number (Mersenne Twister)

! Flag to identify cloud fraction in subcolumns
      logical,  dimension(nsubcol, nlayers) :: iscloudy ! flag that says whether a gridbox is cloudy

! Indices
      integer(kind=im) :: ilev, isubcol, i, n         ! indices

!------------------------------------------------------------------------------------------ 

! Check that irng is in bounds; if not, set to default
      if (irng .ne. 0) irng = 1

! Pass input cloud overlap settings to local variables
      overlap = icld

! Ensure that cloud fractions are in bounds 
      do ilev = 1, nlayers
         cldf(ilev) = cld(ilev)
         if (cldf(ilev) < cldmin) then
            cldf(ilev) = 0._rb
         endif
      enddo

! ----- Create seed  --------

! Advance randum number generator by changeseed values
      if (irng.eq.0) then   
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
! Must use pmid from bottom four layers. 
         if (pmid(1).lt.pmid(2)) then 
            stop 'MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.'
         endif 
         seed1 = (pmid(1) - int(pmid(1)))  * 1000000000_im
         seed2 = (pmid(2) - int(pmid(2)))  * 1000000000_im
         seed3 = (pmid(3) - int(pmid(3)))  * 1000000000_im
         seed4 = (pmid(4) - int(pmid(4)))  * 1000000000_im
         do i=1,changeSeed
            call kissvec(seed1, seed2, seed3, seed4, rand_num)
         enddo
      elseif (irng.eq.1) then
         randomNumbers = new_RandomNumberSequence(seed = changeSeed)
      endif 

! ------ Apply overlap assumption --------

! generate the random numbers  

      select case (overlap)

      case(1) 
! Random overlap
! i) pick a random value at every level
  
         if (irng.eq.0) then 
            do isubcol = 1,nsubcol
               do ilev = 1,nlayers
                  call kissvec(seed1, seed2, seed3, seed4, rand_num)  ! we get different random number for each level
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               do ilev = 1, nlayers
                  rand_num_mt = getRandomReal(randomNumbers)
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

      case(2) 
! Maximum_Random overlap
! i) pick a random number for top layer.
! ii) walk down the column: 
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number 

         if (irng.eq.0) then 
            do isubcol = 1,nsubcol
               do ilev = 1,nlayers
                  call kissvec(seed1, seed2, seed3, seed4, rand_num) 
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               do ilev = 1, nlayers
                  rand_num_mt = getRandomReal(randomNumbers)
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

         do ilev = 2,nlayers
            do isubcol = 1, nsubcol
               if (CDF(isubcol, ilev-1) > 1._rb - cldf(ilev-1) ) then
                  CDF(isubcol,ilev) = CDF(isubcol,ilev-1) 
               else
                  CDF(isubcol,ilev) = CDF(isubcol,ilev) * (1._rb - cldf(ilev-1)) 
               endif
            enddo
         enddo
       
      case(3) 
! Maximum overlap
! i) pick the same random number at every level  

         if (irng.eq.0) then 
            do isubcol = 1,nsubcol
               call kissvec(seed1, seed2, seed3, seed4, rand_num) 
               do ilev = 1,nlayers
                  CDF(isubcol, ilev) = rand_num
               enddo
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               rand_num_mt = getRandomReal(randomNumbers)
               do ilev = 1, nlayers
                  CDF(isubcol,ilev) = rand_num_mt
               enddo
             enddo
         endif

     case(4)
       ! Exponential overlap: transition from maximum to random cloud overlap increases 
       ! exponentially with layer thickness and distance through layers
       !
       ! generate 2 streams of random numbers
       ! CDF2 is used to select which sub-columns are vertically correlated relative to alpha
       ! CDF  is used to select which sub-columns are treated as cloudy relative to cloud fraction
       if (irng.eq.0) then 
          do isubcol = 1,nsubcol
             do ilev = 1,nlayers
                call kissvec(seed1, seed2, seed3, seed4, rand_num)
                CDF(isubcol, ilev) = rand_num
                call kissvec(seed1, seed2, seed3, seed4, rand_num)
                CDF2(isubcol, ilev) = rand_num
             end do
          end do
       elseif (irng.eq.1) then
          do isubcol = 1, nsubcol
             do ilev = 1,nlayers
                rand_num_mt = getRandomReal(randomNumbers)
                CDF(isubcol,ilev) = rand_num_mt
                rand_num_mt = getRandomReal(randomNumbers)
                CDF2(isubcol,ilev) = rand_num_mt
             enddo
          enddo
       endif

       ! generate vertical correlations in random number arrays: bottom to top
       do ilev = 2,nlayers
          where (CDF2(:, ilev) < spread(alpha (ilev), dim=1, nCopies=nsubcol) )
             CDF(:,ilev) = CDF(:,ilev-1) 
          end where
       end do

     case(5)
       ! Exponential_Random overlap: transition from maximum to random cloud overlap increases 
       ! exponentially with layer thickness and with distance through adjacent cloudy layers. 
       ! Non-adjacent blocks of clouds are treated randomly, and each block begins a new 
       ! exponential transition from maximum to random. 
       !
       ! generate 2 streams of random numbers
       ! CDF2 is used to select which sub-columns are vertically correlated relative to alpha
       ! CDF  is used to select which sub-columns are treated as cloudy relative to cloud fraction
       if (irng.eq.0) then 
          do isubcol = 1,nsubcol
             do ilev = 1,nlayers
                call kissvec(seed1, seed2, seed3, seed4, rand_num)
                CDF(isubcol, ilev) = rand_num
                call kissvec(seed1, seed2, seed3, seed4, rand_num)
                CDF2(isubcol, ilev) = rand_num
             end do
          end do
       elseif (irng.eq.1) then
          do isubcol = 1, nsubcol
             do ilev = 1,nlayers
                rand_num_mt = getRandomReal(randomNumbers)
                CDF(isubcol,ilev) = rand_num_mt
                rand_num_mt = getRandomReal(randomNumbers)
                CDF2(isubcol,ilev) = rand_num_mt
             enddo
          enddo
       endif

       ! generate vertical correlations in random number arrays - bottom to top
       do ilev = 2,nlayers
          where (CDF2(:, ilev) < spread(alpha (ilev), dim=1, nCopies=nsubcol) )
             CDF(:,ilev) = CDF(:,ilev-1) 
          end where
       end do

     end select

! -- generate logical cloud mask across subcolumns including vertical correlations --

      do ilev = 1,nlayers
         iscloudy(:,ilev) = (CDF(:,ilev) >= 1._rb - spread(cldf(ilev), dim=1, nCopies=nsubcol) )
      enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
! where there is a cloud, define the subcolumn cloud properties, 
! otherwise set these to zero

      do ilev = 1,nlayers
         do isubcol = 1, nsubcol
            if ( iscloudy(isubcol,ilev) ) then
               cld_stoch(isubcol,ilev) = 1._rb
               clwp_stoch(isubcol,ilev) = clwp(ilev)
               ciwp_stoch(isubcol,ilev) = ciwp(ilev)
               n = ngb(isubcol)
               tauc_stoch(isubcol,ilev) = tauc(n,ilev)
!               ssac_stoch(isubcol,ilev) = ssac(n,ilev)
!               asmc_stoch(isubcol,ilev) = asmc(n,ilev)
            else
               cld_stoch(isubcol,ilev) = 0._rb
               clwp_stoch(isubcol,ilev) = 0._rb
               ciwp_stoch(isubcol,ilev) = 0._rb
               tauc_stoch(isubcol,ilev) = 0._rb
!               ssac_stoch(isubcol,ilev) = 1._rb
!               asmc_stoch(isubcol,ilev) = 1._rb
            endif
         enddo
      enddo

! -- compute the means of the subcolumns ---
!      mean_cld_stoch(:) = 0._rb
!      mean_clwp_stoch(:) = 0._rb
!      mean_ciwp_stoch(:) = 0._rb
!      mean_tauc_stoch(:) = 0._rb
!      mean_ssac_stoch(:) = 0._rb
!      mean_asmc_stoch(:) = 0._rb
!      do i = 1, nsubcol
!         mean_cld_stoch(:) =  cld_stoch(i,:) + mean_cld_stoch(:) 
!         mean_clwp_stoch( :) =  clwp_stoch( i,:) + mean_clwp_stoch( :) 
!         mean_ciwp_stoch( :) =  ciwp_stoch( i,:) + mean_ciwp_stoch( :) 
!         mean_tauc_stoch( :) =  tauc_stoch( i,:) + mean_tauc_stoch( :) 
!         mean_ssac_stoch( :) =  ssac_stoch( i,:) + mean_ssac_stoch( :) 
!         mean_asmc_stoch( :) =  asmc_stoch( i,:) + mean_asmc_stoch( :) 
!      end do
!      mean_cld_stoch(:) = mean_cld_stoch(:) / nsubcol
!      mean_clwp_stoch( :) = mean_clwp_stoch( :) / nsubcol
!      mean_ciwp_stoch( :) = mean_ciwp_stoch( :) / nsubcol
!      mean_tauc_stoch( :) = mean_tauc_stoch( :) / nsubcol
!      mean_ssac_stoch( :) = mean_ssac_stoch( :) / nsubcol
!      mean_asmc_stoch( :) = mean_asmc_stoch( :) / nsubcol

      end subroutine generate_stochastic_clouds


!------------------------------------------------------------------
! Private subroutines
!------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------- 
      subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)
!-------------------------------------------------------------------------------------------------- 

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
      real(kind=rb), intent(inout)  :: ran_arr
      integer(kind=im), intent(inout) :: seed1,seed2,seed3,seed4
!      integer(kind=im) :: i,sz,kiss
      integer(kind=im) :: kiss
      integer(kind=im) :: m, k, n

! inline function 
      m(k, n) = ieor (k, ishft (k, n) )

!      sz = size(ran_arr)
!      do i = 1, sz
         seed1 = 69069_im * seed1 + 1327217885_im
         seed2 = m (m (m (seed2, 13_im), - 17_im), 5_im)
         seed3 = 18000_im * iand (seed3, 65535_im) + ishft (seed3, - 16_im)
         seed4 = 30903_im * iand (seed4, 65535_im) + ishft (seed4, - 16_im)
         kiss = seed1 + seed2 + ishft (seed3, 16_im) + seed4
         ran_arr = kiss*2.328306e-10_rb + 0.5_rb
!     end do
    
      end subroutine kissvec

      end module mcica_subcol_gen_lw

