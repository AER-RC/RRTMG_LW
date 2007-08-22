!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
      module rrtmg_lw_cldprop

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

      use parkind, only : jpim, jprb 
      use parrrtm, only : nbndlw
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3
      use rrlw_vsn, only: hvrcld, hnamcld

      implicit none

      contains

! ------------------------------------------------------------------------------
      subroutine cldprop(nlayers, inflag, iceflag, liqflag, cldfrac, tauc, &
                         ciwp, clwp, rei, rel, ncbands, taucloud)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

      integer(kind=jpim), intent(in) :: nlayers         ! total number of layers
      integer(kind=jpim), intent(in) :: inflag          ! see definitions
      integer(kind=jpim), intent(in) :: iceflag         ! see definitions
      integer(kind=jpim), intent(in) :: liqflag         ! see definitions

      real(kind=jprb), intent(in) :: cldfrac(:)         ! cloud fraction
                                                        !    Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: ciwp(:)            ! cloud ice water path
                                                        !    Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: clwp(:)            ! cloud liquid water path
                                                        !    Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: rei(:)             ! cloud ice particle size
                                                        !    Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: rel(:)             ! cloud liquid particle size
                                                        !    Dimensions: (nlayers)
      real(kind=jprb), intent(in) :: tauc(:,:)          ! cloud optical depth
                                                        !    Dimensions: (nbndlw,nlayers)

! ------- Output -------

      integer(kind=jpim), intent(out) :: ncbands        ! number of cloud spectral bands
      real(kind=jprb), intent(out) :: taucloud(:,:)     ! cloud optical depth
                                                        !    Dimensions: (nlayers,nbndlw)

! ------- Local -------

      integer(kind=jpim) :: lay                 ! Layer index
      integer(kind=jpim) :: ib                  ! spectral band index
      integer(kind=jpim) :: index 
      integer(kind=jpim) :: icepat
      integer(kind=jpim) :: liqpat
      integer(kind=jpim) :: ipat(16,0:2)

      real(kind=jprb) :: abscoice(nbndlw)       ! ice absorption coefficients
      real(kind=jprb) :: abscoliq(nbndlw)       ! liquid absorption coefficients
      real(kind=jprb) :: cwp                    ! cloud water path
      real(kind=jprb) :: radliq                 ! cloud liquid droplet radius (microns)
      real(kind=jprb) :: radice                 ! cloud ice effective radius (microns)
      real(kind=jprb) :: dgeice                 ! cloud ice generalized effective size
      real(kind=jprb) :: factor                 ! 
      real(kind=jprb) :: fint                   ! 
      real(kind=jprb) :: tauctot(nlayers)       ! band integrated cloud optical depth
      real(kind=jprb), parameter :: eps = 1.e-6 ! epsilon

! ------- Definitions -------

!     Explanation of the method for each value of INFLAG.  Values of
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
!     INFLAG = 2 does distinguish between liquid and ice clouds, and
!     requires further user input to specify the method to be used to 
!     compute the aborption due to each.
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
!                  optical depth are input.  
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
!                  water path (g/m2) are input.  The (gray) cloud optical 
!                  depth is computed as in CCM2.
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
!                  water path (g/m2), and cloud ice fraction are input.
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in CCM3.
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in 
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
!                     spectral regions in this work have been matched with
!                     the spectral bands in RRTM to as great an extent 
!                     as possible:  
!                     E&C 1      IB = 5      RRTM bands 9-16
!                     E&C 2      IB = 4      RRTM bands 6-8
!                     E&C 3      IB = 3      RRTM bands 3-5
!                     E&C 4      IB = 2      RRTM band 2
!                     E&C 5      IB = 1      RRTM band 1
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer 
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CCM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!                     and the optical depths due to water clouds are computed 
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a 
!                     range of effective radii by an averaging procedure 
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption 
!                     coefficients for the input effective radius.

      data ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
                 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

      hvrcld = '$Revision$'

      ncbands = 1
      tauctot(:) = 0._jprb

      do lay = 1, nlayers
         do ib = 1, nbndlw
            taucloud(lay,ib) = 0.0_jprb
            tauctot(lay) = tauctot(lay) + tauc(ib,lay)
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers
         cwp = ciwp(lay) + clwp(lay)
         if (cldfrac(lay) .ge. eps .and. &
            (cwp .ge. eps .or. tauctot(lay) .ge. eps)) then

! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
               ncbands = 16
               do ib = 1, ncbands
                  taucloud(lay,ib) = tauc(ib,lay)
               end do

            elseif (inflag .eq. 1) then
               ncbands = 16
               do ib = 1, ncbands
                  taucloud(lay,ib) = abscld1 * cwp
               end do

! Separate treatement of ice clouds and water clouds.
            elseif (inflag .eq. 2) then
               radice = rei(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwp(lay) .eq. 0.0_jprb) then
                  abscoice(1) = 0.0_jprb
                  icepat = 0

               elseif (iceflag .eq. 0) then
                  if (radice .lt. 10.0_jprb) stop 'ICE RADIUS TOO SMALL'
                  abscoice(1) = absice0(1) + absice0(2)/radice
                  icepat = 0

               elseif (iceflag .eq. 1) then
                  if (radice .lt. 13.0_jprb .or. radice .gt. 130._jprb) stop &
                       'ICE RADIUS OUT OF BOUNDS'
                  ncbands = 5
                  do ib = 1, ncbands
                     abscoice(ib) = absice1(1,ib) + absice1(2,ib)/radice
                  enddo
                  icepat = 1

! For iceflag=2 option, combine with iceflag=0 option to handle out of bounds 
! particle sizes.
! Use iceflag=2 option for ice particle effective radii from 5.0 and 131.0 microns
! and use iceflag=0 option for ice particles greater than 131.0 microns.
! *** NOTE: Transition between two methods has not been smoothed. 

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0_jprb) stop 'ICE RADIUS OUT OF BOUNDS'
                  if (radice .ge. 5.0_jprb .and. radice .le. 131._jprb) then
                     ncbands = 16
                     factor = (radice - 2._jprb)/3._jprb
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - float(index)
                     do ib = 1, ncbands
                        abscoice(ib) = &
                            absice2(index,ib) + fint * &
                            (absice2(index+1,ib) - (absice2(index,ib)))
                     enddo
                     icepat = 2
                  elseif (radice .gt. 131._jprb) then
                     do ib = 1, ncbands
                       abscoice(ib) = absice0(1) + absice0(2)/radice
                       icepat = 0
                     enddo
                  endif

! For iceflag=3 option, combine with iceflag=0 option to handle large particle sizes.
! Use iceflag=3 option for ice particle effective radii from 5.0 and 91.0 microns
! (generalized effective size, dge, from 8 to 140 microns), and use iceflag=0 option
! for ice particle effective radii greater than 91.0 microns (dge = 140 microns).
! *** NOTE: Fu parameterization requires particle size in generalized effective size.
! *** NOTE: Transition between two methods has not been smoothed. 

               elseif (iceflag .eq. 3) then
                   dgeice = radice
                  if (dgeice .lt. 5.0_jprb) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                  if (dgeice .ge. 5.0_jprb .and. dgeice .le. 140._jprb) then
                     ncbands = 16
                     factor = (dgeice - 2._jprb)/3._jprb
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - float(index)
                     do ib = 1, ncbands
                        abscoice(ib) = &
                          absice3(index,ib) + fint * &
                          (absice3(index+1,ib) - (absice3(index,ib)))
                     enddo
                     icepat = 2
                  elseif (dgeice .gt. 140._jprb) then
                     do ib = 1, ncbands
                       abscoice(ib) = absice0(1) + absice0(2)/radice
                       icepat = 0
                     enddo
                  endif
   
               endif
                  
! Calculation of absorption coefficients due to water clouds.
               if (clwp(lay) .eq. 0.0_jprb) then
                  abscoliq(1) = 0.0_jprb
                  liqpat = 0
                  if (icepat .eq. 1) icepat = 2

               elseif (liqflag .eq. 0) then
                  abscoliq(1) = absliq0
                  liqpat = 0
                  if (icepat .eq. 1) icepat = 2

               elseif (liqflag .eq. 1) then
                  radliq = rel(lay)
                  if (radliq .lt. 1.5_jprb .or. radliq .gt. 60._jprb) stop &
                       'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                  index = radliq - 1.5_jprb
                  if (index .eq. 58) index = 57
                  if (index .eq. 0) index = 1
                  fint = radliq - 1.5_jprb - index
                  ncbands = 16
                  do ib = 1, ncbands
                     abscoliq(ib) = &
                         absliq1(index,ib) + fint * &
                         (absliq1(index+1,ib) - (absliq1(index,ib)))
                  enddo
                  liqpat = 2
               endif

               do ib = 1, ncbands
                  taucloud(lay,ib) = ciwp(lay) * abscoice(ipat(ib,icepat)) + &
                                     clwp(lay) * abscoliq(ipat(ib,liqpat))
               enddo
            endif
         endif
      enddo

      end subroutine cldprop

!***************************************************************************
      subroutine lwcldpr
!***************************************************************************

      save

! ABSCLDn is the liquid water absorption coefficient (m2/g). 
! For INFLAG = 1.
      abscld1 = 0.0602410_jprb
!  
! Everything below is for INFLAG = 2.

! ABSICEn(J,IB) are the parameters needed to compute the liquid water 
! absorption coefficient in spectral region IB for ICEFLAG=n.  The units
! of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units (microns (m2/g)).
! For ICEFLAG = 0.

      absice0(:)= (/0.005_jprb,  1.0_jprb/)

! For ICEFLAG = 1.
      absice1(1,:) = (/0.0036_jprb, 0.0068_jprb, 0.0003_jprb, 0.0016_jprb, 0.0020_jprb/)
      absice1(2,:) = (/1.136_jprb , 0.600_jprb , 1.338_jprb , 1.166_jprb , 1.118_jprb /)

! For ICEFLAG = 2.  In each band, the absorption
! coefficients are listed for a range of effective radii from 5.0
! to 131.0 microns in increments of 3.0 microns.
! Spherical Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice2(:,1) = (/ &
! band 1
       7.798999e-02_jprb,6.340479e-02_jprb,5.417973e-02_jprb,4.766245e-02_jprb,4.272663e-02_jprb, &
       3.880939e-02_jprb,3.559544e-02_jprb,3.289241e-02_jprb,3.057511e-02_jprb,2.855800e-02_jprb, &
       2.678022e-02_jprb,2.519712e-02_jprb,2.377505e-02_jprb,2.248806e-02_jprb,2.131578e-02_jprb, &
       2.024194e-02_jprb,1.925337e-02_jprb,1.833926e-02_jprb,1.749067e-02_jprb,1.670007e-02_jprb, &
       1.596113e-02_jprb,1.526845e-02_jprb,1.461739e-02_jprb,1.400394e-02_jprb,1.342462e-02_jprb, &
       1.287639e-02_jprb,1.235656e-02_jprb,1.186279e-02_jprb,1.139297e-02_jprb,1.094524e-02_jprb, &
       1.051794e-02_jprb,1.010956e-02_jprb,9.718755e-03_jprb,9.344316e-03_jprb,8.985139e-03_jprb, &
       8.640223e-03_jprb,8.308656e-03_jprb,7.989606e-03_jprb,7.682312e-03_jprb,7.386076e-03_jprb, &
       7.100255e-03_jprb,6.824258e-03_jprb,6.557540e-03_jprb/)
      absice2(:,2) = (/ &
! band 2
       2.784879e-02_jprb,2.709863e-02_jprb,2.619165e-02_jprb,2.529230e-02_jprb,2.443225e-02_jprb, &
       2.361575e-02_jprb,2.284021e-02_jprb,2.210150e-02_jprb,2.139548e-02_jprb,2.071840e-02_jprb, &
       2.006702e-02_jprb,1.943856e-02_jprb,1.883064e-02_jprb,1.824120e-02_jprb,1.766849e-02_jprb, &
       1.711099e-02_jprb,1.656737e-02_jprb,1.603647e-02_jprb,1.551727e-02_jprb,1.500886e-02_jprb, &
       1.451045e-02_jprb,1.402132e-02_jprb,1.354084e-02_jprb,1.306842e-02_jprb,1.260355e-02_jprb, &
       1.214575e-02_jprb,1.169460e-02_jprb,1.124971e-02_jprb,1.081072e-02_jprb,1.037731e-02_jprb, &
       9.949167e-03_jprb,9.526021e-03_jprb,9.107615e-03_jprb,8.693714e-03_jprb,8.284096e-03_jprb, &
       7.878558e-03_jprb,7.476910e-03_jprb,7.078974e-03_jprb,6.684586e-03_jprb,6.293589e-03_jprb, &
       5.905839e-03_jprb,5.521200e-03_jprb,5.139543e-03_jprb/)
      absice2(:,3) = (/ &
! band 3
       1.065397e-01_jprb,8.005726e-02_jprb,6.546428e-02_jprb,5.589131e-02_jprb,4.898681e-02_jprb, &
       4.369932e-02_jprb,3.947901e-02_jprb,3.600676e-02_jprb,3.308299e-02_jprb,3.057561e-02_jprb, &
       2.839325e-02_jprb,2.647040e-02_jprb,2.475872e-02_jprb,2.322164e-02_jprb,2.183091e-02_jprb, &
       2.056430e-02_jprb,1.940407e-02_jprb,1.833586e-02_jprb,1.734787e-02_jprb,1.643034e-02_jprb, &
       1.557512e-02_jprb,1.477530e-02_jprb,1.402501e-02_jprb,1.331924e-02_jprb,1.265364e-02_jprb, &
       1.202445e-02_jprb,1.142838e-02_jprb,1.086257e-02_jprb,1.032445e-02_jprb,9.811791e-03_jprb, &
       9.322587e-03_jprb,8.855053e-03_jprb,8.407591e-03_jprb,7.978763e-03_jprb,7.567273e-03_jprb, &
       7.171949e-03_jprb,6.791728e-03_jprb,6.425642e-03_jprb,6.072809e-03_jprb,5.732424e-03_jprb, &
       5.403748e-03_jprb,5.086103e-03_jprb,4.778865e-03_jprb/)
      absice2(:,4) = (/ &
! band 4
       1.804566e-01_jprb,1.168987e-01_jprb,8.680442e-02_jprb,6.910060e-02_jprb,5.738174e-02_jprb, &
       4.902332e-02_jprb,4.274585e-02_jprb,3.784923e-02_jprb,3.391734e-02_jprb,3.068690e-02_jprb, &
       2.798301e-02_jprb,2.568480e-02_jprb,2.370600e-02_jprb,2.198337e-02_jprb,2.046940e-02_jprb, &
       1.912777e-02_jprb,1.793016e-02_jprb,1.685420e-02_jprb,1.588193e-02_jprb,1.499882e-02_jprb, &
       1.419293e-02_jprb,1.345440e-02_jprb,1.277496e-02_jprb,1.214769e-02_jprb,1.156669e-02_jprb, &
       1.102694e-02_jprb,1.052412e-02_jprb,1.005451e-02_jprb,9.614854e-03_jprb,9.202335e-03_jprb, &
       8.814470e-03_jprb,8.449077e-03_jprb,8.104223e-03_jprb,7.778195e-03_jprb,7.469466e-03_jprb, &
       7.176671e-03_jprb,6.898588e-03_jprb,6.634117e-03_jprb,6.382264e-03_jprb,6.142134e-03_jprb, &
       5.912913e-03_jprb,5.693862e-03_jprb,5.484308e-03_jprb/)
      absice2(:,5) = (/ &
! band 5
       2.131806e-01_jprb,1.311372e-01_jprb,9.407171e-02_jprb,7.299442e-02_jprb,5.941273e-02_jprb, &
       4.994043e-02_jprb,4.296242e-02_jprb,3.761113e-02_jprb,3.337910e-02_jprb,2.994978e-02_jprb, &
       2.711556e-02_jprb,2.473461e-02_jprb,2.270681e-02_jprb,2.095943e-02_jprb,1.943839e-02_jprb, &
       1.810267e-02_jprb,1.692057e-02_jprb,1.586719e-02_jprb,1.492275e-02_jprb,1.407132e-02_jprb, &
       1.329989e-02_jprb,1.259780e-02_jprb,1.195618e-02_jprb,1.136761e-02_jprb,1.082583e-02_jprb, &
       1.032552e-02_jprb,9.862158e-03_jprb,9.431827e-03_jprb,9.031157e-03_jprb,8.657217e-03_jprb, &
       8.307449e-03_jprb,7.979609e-03_jprb,7.671724e-03_jprb,7.382048e-03_jprb,7.109032e-03_jprb, &
       6.851298e-03_jprb,6.607615e-03_jprb,6.376881e-03_jprb,6.158105e-03_jprb,5.950394e-03_jprb, &
       5.752942e-03_jprb,5.565019e-03_jprb,5.385963e-03_jprb/)
      absice2(:,6) = (/ &
! band 6
       1.546177e-01_jprb,1.039251e-01_jprb,7.910347e-02_jprb,6.412429e-02_jprb,5.399997e-02_jprb, &
       4.664937e-02_jprb,4.104237e-02_jprb,3.660781e-02_jprb,3.300218e-02_jprb,3.000586e-02_jprb, &
       2.747148e-02_jprb,2.529633e-02_jprb,2.340647e-02_jprb,2.174723e-02_jprb,2.027731e-02_jprb, &
       1.896487e-02_jprb,1.778492e-02_jprb,1.671761e-02_jprb,1.574692e-02_jprb,1.485978e-02_jprb, &
       1.404543e-02_jprb,1.329489e-02_jprb,1.260066e-02_jprb,1.195636e-02_jprb,1.135657e-02_jprb, &
       1.079664e-02_jprb,1.027257e-02_jprb,9.780871e-03_jprb,9.318505e-03_jprb,8.882815e-03_jprb, &
       8.471458e-03_jprb,8.082364e-03_jprb,7.713696e-03_jprb,7.363817e-03_jprb,7.031264e-03_jprb, &
       6.714725e-03_jprb,6.413021e-03_jprb,6.125086e-03_jprb,5.849958e-03_jprb,5.586764e-03_jprb, &
       5.334707e-03_jprb,5.093066e-03_jprb,4.861179e-03_jprb/)
      absice2(:,7) = (/ &
! band 7
       7.583404e-02_jprb,6.181558e-02_jprb,5.312027e-02_jprb,4.696039e-02_jprb,4.225986e-02_jprb, &
       3.849735e-02_jprb,3.538340e-02_jprb,3.274182e-02_jprb,3.045798e-02_jprb,2.845343e-02_jprb, &
       2.667231e-02_jprb,2.507353e-02_jprb,2.362606e-02_jprb,2.230595e-02_jprb,2.109435e-02_jprb, &
       1.997617e-02_jprb,1.893916e-02_jprb,1.797328e-02_jprb,1.707016e-02_jprb,1.622279e-02_jprb, &
       1.542523e-02_jprb,1.467241e-02_jprb,1.395997e-02_jprb,1.328414e-02_jprb,1.264164e-02_jprb, &
       1.202958e-02_jprb,1.144544e-02_jprb,1.088697e-02_jprb,1.035218e-02_jprb,9.839297e-03_jprb, &
       9.346733e-03_jprb,8.873057e-03_jprb,8.416980e-03_jprb,7.977335e-03_jprb,7.553066e-03_jprb, &
       7.143210e-03_jprb,6.746888e-03_jprb,6.363297e-03_jprb,5.991700e-03_jprb,5.631422e-03_jprb, &
       5.281840e-03_jprb,4.942378e-03_jprb,4.612505e-03_jprb/)
      absice2(:,8) = (/ &
! band 8
       9.022185e-02_jprb,6.922700e-02_jprb,5.710674e-02_jprb,4.898377e-02_jprb,4.305946e-02_jprb, &
       3.849553e-02_jprb,3.484183e-02_jprb,3.183220e-02_jprb,2.929794e-02_jprb,2.712627e-02_jprb, &
       2.523856e-02_jprb,2.357810e-02_jprb,2.210286e-02_jprb,2.078089e-02_jprb,1.958747e-02_jprb, &
       1.850310e-02_jprb,1.751218e-02_jprb,1.660205e-02_jprb,1.576232e-02_jprb,1.498440e-02_jprb, &
       1.426107e-02_jprb,1.358624e-02_jprb,1.295474e-02_jprb,1.236212e-02_jprb,1.180456e-02_jprb, &
       1.127874e-02_jprb,1.078175e-02_jprb,1.031106e-02_jprb,9.864433e-03_jprb,9.439878e-03_jprb, &
       9.035637e-03_jprb,8.650140e-03_jprb,8.281981e-03_jprb,7.929895e-03_jprb,7.592746e-03_jprb, &
       7.269505e-03_jprb,6.959238e-03_jprb,6.661100e-03_jprb,6.374317e-03_jprb,6.098185e-03_jprb, &
       5.832059e-03_jprb,5.575347e-03_jprb,5.327504e-03_jprb/)
      absice2(:,9) = (/ &
! band 9
       1.294087e-01_jprb,8.788217e-02_jprb,6.728288e-02_jprb,5.479720e-02_jprb,4.635049e-02_jprb, &
       4.022253e-02_jprb,3.555576e-02_jprb,3.187259e-02_jprb,2.888498e-02_jprb,2.640843e-02_jprb, &
       2.431904e-02_jprb,2.253038e-02_jprb,2.098024e-02_jprb,1.962267e-02_jprb,1.842293e-02_jprb, &
       1.735426e-02_jprb,1.639571e-02_jprb,1.553060e-02_jprb,1.474552e-02_jprb,1.402953e-02_jprb, &
       1.337363e-02_jprb,1.277033e-02_jprb,1.221336e-02_jprb,1.169741e-02_jprb,1.121797e-02_jprb, &
       1.077117e-02_jprb,1.035369e-02_jprb,9.962643e-03_jprb,9.595509e-03_jprb,9.250088e-03_jprb, &
       8.924447e-03_jprb,8.616876e-03_jprb,8.325862e-03_jprb,8.050057e-03_jprb,7.788258e-03_jprb, &
       7.539388e-03_jprb,7.302478e-03_jprb,7.076656e-03_jprb,6.861134e-03_jprb,6.655197e-03_jprb, &
       6.458197e-03_jprb,6.269543e-03_jprb,6.088697e-03_jprb/)
      absice2(:,10) = (/ &
! band 10
       1.593628e-01_jprb,1.014552e-01_jprb,7.458955e-02_jprb,5.903571e-02_jprb,4.887582e-02_jprb, &
       4.171159e-02_jprb,3.638480e-02_jprb,3.226692e-02_jprb,2.898717e-02_jprb,2.631256e-02_jprb, &
       2.408925e-02_jprb,2.221156e-02_jprb,2.060448e-02_jprb,1.921325e-02_jprb,1.799699e-02_jprb, &
       1.692456e-02_jprb,1.597177e-02_jprb,1.511961e-02_jprb,1.435289e-02_jprb,1.365933e-02_jprb, &
       1.302890e-02_jprb,1.245334e-02_jprb,1.192576e-02_jprb,1.144037e-02_jprb,1.099230e-02_jprb, &
       1.057739e-02_jprb,1.019208e-02_jprb,9.833302e-03_jprb,9.498395e-03_jprb,9.185047e-03_jprb, &
       8.891237e-03_jprb,8.615185e-03_jprb,8.355325e-03_jprb,8.110267e-03_jprb,7.878778e-03_jprb, &
       7.659759e-03_jprb,7.452224e-03_jprb,7.255291e-03_jprb,7.068166e-03_jprb,6.890130e-03_jprb, &
       6.720536e-03_jprb,6.558794e-03_jprb,6.404371e-03_jprb/)
      absice2(:,11) = (/ &
! band 11
       1.656227e-01_jprb,1.032129e-01_jprb,7.487359e-02_jprb,5.871431e-02_jprb,4.828355e-02_jprb, &
       4.099989e-02_jprb,3.562924e-02_jprb,3.150755e-02_jprb,2.824593e-02_jprb,2.560156e-02_jprb, &
       2.341503e-02_jprb,2.157740e-02_jprb,2.001169e-02_jprb,1.866199e-02_jprb,1.748669e-02_jprb, &
       1.645421e-02_jprb,1.554015e-02_jprb,1.472535e-02_jprb,1.399457e-02_jprb,1.333553e-02_jprb, &
       1.273821e-02_jprb,1.219440e-02_jprb,1.169725e-02_jprb,1.124104e-02_jprb,1.082096e-02_jprb, &
       1.043290e-02_jprb,1.007336e-02_jprb,9.739338e-03_jprb,9.428223e-03_jprb,9.137756e-03_jprb, &
       8.865964e-03_jprb,8.611115e-03_jprb,8.371686e-03_jprb,8.146330e-03_jprb,7.933852e-03_jprb, &
       7.733187e-03_jprb,7.543386e-03_jprb,7.363597e-03_jprb,7.193056e-03_jprb,7.031072e-03_jprb, &
       6.877024e-03_jprb,6.730348e-03_jprb,6.590531e-03_jprb/)
      absice2(:,12) = (/ &
! band 12
       9.194591e-02_jprb,6.446867e-02_jprb,4.962034e-02_jprb,4.042061e-02_jprb,3.418456e-02_jprb, &
       2.968856e-02_jprb,2.629900e-02_jprb,2.365572e-02_jprb,2.153915e-02_jprb,1.980791e-02_jprb, &
       1.836689e-02_jprb,1.714979e-02_jprb,1.610900e-02_jprb,1.520946e-02_jprb,1.442476e-02_jprb, &
       1.373468e-02_jprb,1.312345e-02_jprb,1.257858e-02_jprb,1.209010e-02_jprb,1.164990e-02_jprb, &
       1.125136e-02_jprb,1.088901e-02_jprb,1.055827e-02_jprb,1.025531e-02_jprb,9.976896e-03_jprb, &
       9.720255e-03_jprb,9.483022e-03_jprb,9.263160e-03_jprb,9.058902e-03_jprb,8.868710e-03_jprb, &
       8.691240e-03_jprb,8.525312e-03_jprb,8.369886e-03_jprb,8.224042e-03_jprb,8.086961e-03_jprb, &
       7.957917e-03_jprb,7.836258e-03_jprb,7.721400e-03_jprb,7.612821e-03_jprb,7.510045e-03_jprb, &
       7.412648e-03_jprb,7.320242e-03_jprb,7.232476e-03_jprb/)
      absice2(:,13) = (/ &
! band 13
       1.437021e-01_jprb,8.872535e-02_jprb,6.392420e-02_jprb,4.991833e-02_jprb,4.096790e-02_jprb, &
       3.477881e-02_jprb,3.025782e-02_jprb,2.681909e-02_jprb,2.412102e-02_jprb,2.195132e-02_jprb, &
       2.017124e-02_jprb,1.868641e-02_jprb,1.743044e-02_jprb,1.635529e-02_jprb,1.542540e-02_jprb, &
       1.461388e-02_jprb,1.390003e-02_jprb,1.326766e-02_jprb,1.270395e-02_jprb,1.219860e-02_jprb, &
       1.174326e-02_jprb,1.133107e-02_jprb,1.095637e-02_jprb,1.061442e-02_jprb,1.030126e-02_jprb, &
       1.001352e-02_jprb,9.748340e-03_jprb,9.503256e-03_jprb,9.276155e-03_jprb,9.065205e-03_jprb, &
       8.868808e-03_jprb,8.685571e-03_jprb,8.514268e-03_jprb,8.353820e-03_jprb,8.203272e-03_jprb, &
       8.061776e-03_jprb,7.928578e-03_jprb,7.803001e-03_jprb,7.684443e-03_jprb,7.572358e-03_jprb, &
       7.466258e-03_jprb,7.365701e-03_jprb,7.270286e-03_jprb/)
      absice2(:,14) = (/ &
! band 14
       1.288870e-01_jprb,8.160295e-02_jprb,5.964745e-02_jprb,4.703790e-02_jprb,3.888637e-02_jprb, &
       3.320115e-02_jprb,2.902017e-02_jprb,2.582259e-02_jprb,2.330224e-02_jprb,2.126754e-02_jprb, &
       1.959258e-02_jprb,1.819130e-02_jprb,1.700289e-02_jprb,1.598320e-02_jprb,1.509942e-02_jprb, &
       1.432666e-02_jprb,1.364572e-02_jprb,1.304156e-02_jprb,1.250220e-02_jprb,1.201803e-02_jprb, &
       1.158123e-02_jprb,1.118537e-02_jprb,1.082513e-02_jprb,1.049605e-02_jprb,1.019440e-02_jprb, &
       9.916989e-03_jprb,9.661116e-03_jprb,9.424457e-03_jprb,9.205005e-03_jprb,9.001022e-03_jprb, &
       8.810992e-03_jprb,8.633588e-03_jprb,8.467646e-03_jprb,8.312137e-03_jprb,8.166151e-03_jprb, &
       8.028878e-03_jprb,7.899597e-03_jprb,7.777663e-03_jprb,7.662498e-03_jprb,7.553581e-03_jprb, &
       7.450444e-03_jprb,7.352662e-03_jprb,7.259851e-03_jprb/)
      absice2(:,15) = (/ &
! band 15
       8.254229e-02_jprb,5.808787e-02_jprb,4.492166e-02_jprb,3.675028e-02_jprb,3.119623e-02_jprb, &
       2.718045e-02_jprb,2.414450e-02_jprb,2.177073e-02_jprb,1.986526e-02_jprb,1.830306e-02_jprb, &
       1.699991e-02_jprb,1.589698e-02_jprb,1.495199e-02_jprb,1.413374e-02_jprb,1.341870e-02_jprb, &
       1.278883e-02_jprb,1.223002e-02_jprb,1.173114e-02_jprb,1.128322e-02_jprb,1.087900e-02_jprb, &
       1.051254e-02_jprb,1.017890e-02_jprb,9.873991e-03_jprb,9.594347e-03_jprb,9.337044e-03_jprb, &
       9.099589e-03_jprb,8.879842e-03_jprb,8.675960e-03_jprb,8.486341e-03_jprb,8.309594e-03_jprb, &
       8.144500e-03_jprb,7.989986e-03_jprb,7.845109e-03_jprb,7.709031e-03_jprb,7.581007e-03_jprb, &
       7.460376e-03_jprb,7.346544e-03_jprb,7.238978e-03_jprb,7.137201e-03_jprb,7.040780e-03_jprb, &
       6.949325e-03_jprb,6.862483e-03_jprb,6.779931e-03_jprb/)
      absice2(:,16) = (/ &
! band 16
       1.382062e-01_jprb,8.643227e-02_jprb,6.282935e-02_jprb,4.934783e-02_jprb,4.063891e-02_jprb, &
       3.455591e-02_jprb,3.007059e-02_jprb,2.662897e-02_jprb,2.390631e-02_jprb,2.169972e-02_jprb, &
       1.987596e-02_jprb,1.834393e-02_jprb,1.703924e-02_jprb,1.591513e-02_jprb,1.493679e-02_jprb, &
       1.407780e-02_jprb,1.331775e-02_jprb,1.264061e-02_jprb,1.203364e-02_jprb,1.148655e-02_jprb, &
       1.099099e-02_jprb,1.054006e-02_jprb,1.012807e-02_jprb,9.750215e-03_jprb,9.402477e-03_jprb, &
       9.081428e-03_jprb,8.784143e-03_jprb,8.508107e-03_jprb,8.251146e-03_jprb,8.011373e-03_jprb, &
       7.787140e-03_jprb,7.577002e-03_jprb,7.379687e-03_jprb,7.194071e-03_jprb,7.019158e-03_jprb, &
       6.854061e-03_jprb,6.697986e-03_jprb,6.550224e-03_jprb,6.410138e-03_jprb,6.277153e-03_jprb, &
       6.150751e-03_jprb,6.030462e-03_jprb,5.915860e-03_jprb/)

! ICEFLAG = 3; Fu parameterization. Particle size 5 - 140 micron in 
! increments of 3 microns.
! units = m2/g
! Hexagonal Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice3(:,1) = (/ &
! band 1
       3.110649e-03_jprb,4.666352e-02_jprb,6.606447e-02_jprb,6.531678e-02_jprb,6.012598e-02_jprb, &
       5.437494e-02_jprb,4.906411e-02_jprb,4.441146e-02_jprb,4.040585e-02_jprb,3.697334e-02_jprb, &
       3.403027e-02_jprb,3.149979e-02_jprb,2.931596e-02_jprb,2.742365e-02_jprb,2.577721e-02_jprb, &
       2.433888e-02_jprb,2.307732e-02_jprb,2.196644e-02_jprb,2.098437e-02_jprb,2.011264e-02_jprb, &
       1.933561e-02_jprb,1.863992e-02_jprb,1.801407e-02_jprb,1.744812e-02_jprb,1.693346e-02_jprb, &
       1.646252e-02_jprb,1.602866e-02_jprb,1.562600e-02_jprb,1.524933e-02_jprb,1.489399e-02_jprb, &
       1.455580e-02_jprb,1.423098e-02_jprb,1.391612e-02_jprb,1.360812e-02_jprb,1.330413e-02_jprb, &
       1.300156e-02_jprb,1.269801e-02_jprb,1.239127e-02_jprb,1.207928e-02_jprb,1.176014e-02_jprb, &
       1.143204e-02_jprb,1.109334e-02_jprb,1.074243e-02_jprb,1.037786e-02_jprb,9.998198e-03_jprb, &
       9.602126e-03_jprb/)
      absice3(:,2) = (/ &
! band 2
       3.984966e-04_jprb,1.681097e-02_jprb,2.627680e-02_jprb,2.767465e-02_jprb,2.700722e-02_jprb, &
       2.579180e-02_jprb,2.448677e-02_jprb,2.323890e-02_jprb,2.209096e-02_jprb,2.104882e-02_jprb, &
       2.010547e-02_jprb,1.925003e-02_jprb,1.847128e-02_jprb,1.775883e-02_jprb,1.710358e-02_jprb, &
       1.649769e-02_jprb,1.593449e-02_jprb,1.540829e-02_jprb,1.491429e-02_jprb,1.444837e-02_jprb, &
       1.400704e-02_jprb,1.358729e-02_jprb,1.318654e-02_jprb,1.280258e-02_jprb,1.243346e-02_jprb, &
       1.207750e-02_jprb,1.173325e-02_jprb,1.139941e-02_jprb,1.107487e-02_jprb,1.075861e-02_jprb, &
       1.044975e-02_jprb,1.014753e-02_jprb,9.851229e-03_jprb,9.560240e-03_jprb,9.274003e-03_jprb, &
       8.992020e-03_jprb,8.713845e-03_jprb,8.439074e-03_jprb,8.167346e-03_jprb,7.898331e-03_jprb, &
       7.631734e-03_jprb,7.367286e-03_jprb,7.104742e-03_jprb,6.843882e-03_jprb,6.584504e-03_jprb, &
       6.326424e-03_jprb/)
      absice3(:,3) = (/ &
! band 3
       6.933163e-02_jprb,8.540475e-02_jprb,7.701816e-02_jprb,6.771158e-02_jprb,5.986953e-02_jprb, &
       5.348120e-02_jprb,4.824962e-02_jprb,4.390563e-02_jprb,4.024411e-02_jprb,3.711404e-02_jprb, &
       3.440426e-02_jprb,3.203200e-02_jprb,2.993478e-02_jprb,2.806474e-02_jprb,2.638464e-02_jprb, &
       2.486516e-02_jprb,2.348288e-02_jprb,2.221890e-02_jprb,2.105780e-02_jprb,1.998687e-02_jprb, &
       1.899552e-02_jprb,1.807490e-02_jprb,1.721750e-02_jprb,1.641693e-02_jprb,1.566773e-02_jprb, &
       1.496515e-02_jprb,1.430509e-02_jprb,1.368398e-02_jprb,1.309865e-02_jprb,1.254634e-02_jprb, &
       1.202456e-02_jprb,1.153114e-02_jprb,1.106409e-02_jprb,1.062166e-02_jprb,1.020224e-02_jprb, &
       9.804381e-03_jprb,9.426771e-03_jprb,9.068205e-03_jprb,8.727578e-03_jprb,8.403876e-03_jprb, &
       8.096160e-03_jprb,7.803564e-03_jprb,7.525281e-03_jprb,7.260560e-03_jprb,7.008697e-03_jprb, &
       6.769036e-03_jprb/)
      absice3(:,4) = (/ &
! band 4
       1.765735e-01_jprb,1.382700e-01_jprb,1.095129e-01_jprb,8.987475e-02_jprb,7.591185e-02_jprb, &
       6.554169e-02_jprb,5.755500e-02_jprb,5.122083e-02_jprb,4.607610e-02_jprb,4.181475e-02_jprb, &
       3.822697e-02_jprb,3.516432e-02_jprb,3.251897e-02_jprb,3.021073e-02_jprb,2.817876e-02_jprb, &
       2.637607e-02_jprb,2.476582e-02_jprb,2.331871e-02_jprb,2.201113e-02_jprb,2.082388e-02_jprb, &
       1.974115e-02_jprb,1.874983e-02_jprb,1.783894e-02_jprb,1.699922e-02_jprb,1.622280e-02_jprb, &
       1.550296e-02_jprb,1.483390e-02_jprb,1.421064e-02_jprb,1.362880e-02_jprb,1.308460e-02_jprb, &
       1.257468e-02_jprb,1.209611e-02_jprb,1.164628e-02_jprb,1.122287e-02_jprb,1.082381e-02_jprb, &
       1.044725e-02_jprb,1.009154e-02_jprb,9.755166e-03_jprb,9.436783e-03_jprb,9.135163e-03_jprb, &
       8.849193e-03_jprb,8.577856e-03_jprb,8.320225e-03_jprb,8.075451e-03_jprb,7.842755e-03_jprb, &
       7.621418e-03_jprb/)
      absice3(:,5) = (/ &
! band 5
       2.339673e-01_jprb,1.692124e-01_jprb,1.291656e-01_jprb,1.033837e-01_jprb,8.562949e-02_jprb, &
       7.273526e-02_jprb,6.298262e-02_jprb,5.537015e-02_jprb,4.927787e-02_jprb,4.430246e-02_jprb, &
       4.017061e-02_jprb,3.669072e-02_jprb,3.372455e-02_jprb,3.116995e-02_jprb,2.894977e-02_jprb, &
       2.700471e-02_jprb,2.528842e-02_jprb,2.376420e-02_jprb,2.240256e-02_jprb,2.117959e-02_jprb, &
       2.007567e-02_jprb,1.907456e-02_jprb,1.816271e-02_jprb,1.732874e-02_jprb,1.656300e-02_jprb, &
       1.585725e-02_jprb,1.520445e-02_jprb,1.459852e-02_jprb,1.403419e-02_jprb,1.350689e-02_jprb, &
       1.301260e-02_jprb,1.254781e-02_jprb,1.210941e-02_jprb,1.169468e-02_jprb,1.130118e-02_jprb, &
       1.092675e-02_jprb,1.056945e-02_jprb,1.022757e-02_jprb,9.899560e-03_jprb,9.584021e-03_jprb, &
       9.279705e-03_jprb,8.985479e-03_jprb,8.700322e-03_jprb,8.423306e-03_jprb,8.153590e-03_jprb, &
       7.890412e-03_jprb/)
      absice3(:,6) = (/ &
! band 6
       1.145369e-01_jprb,1.174566e-01_jprb,9.917866e-02_jprb,8.332990e-02_jprb,7.104263e-02_jprb, &
       6.153370e-02_jprb,5.405472e-02_jprb,4.806281e-02_jprb,4.317918e-02_jprb,3.913795e-02_jprb, &
       3.574916e-02_jprb,3.287437e-02_jprb,3.041067e-02_jprb,2.828017e-02_jprb,2.642292e-02_jprb, &
       2.479206e-02_jprb,2.335051e-02_jprb,2.206851e-02_jprb,2.092195e-02_jprb,1.989108e-02_jprb, &
       1.895958e-02_jprb,1.811385e-02_jprb,1.734245e-02_jprb,1.663573e-02_jprb,1.598545e-02_jprb, &
       1.538456e-02_jprb,1.482700e-02_jprb,1.430750e-02_jprb,1.382150e-02_jprb,1.336499e-02_jprb, &
       1.293447e-02_jprb,1.252685e-02_jprb,1.213939e-02_jprb,1.176968e-02_jprb,1.141555e-02_jprb, &
       1.107508e-02_jprb,1.074655e-02_jprb,1.042839e-02_jprb,1.011923e-02_jprb,9.817799e-03_jprb, &
       9.522962e-03_jprb,9.233688e-03_jprb,8.949041e-03_jprb,8.668171e-03_jprb,8.390301e-03_jprb, &
       8.114723e-03_jprb/)
      absice3(:,7) = (/ &
! band 7
       1.222345e-02_jprb,5.344230e-02_jprb,5.523465e-02_jprb,5.128759e-02_jprb,4.676925e-02_jprb, &
       4.266150e-02_jprb,3.910561e-02_jprb,3.605479e-02_jprb,3.342843e-02_jprb,3.115052e-02_jprb, &
       2.915776e-02_jprb,2.739935e-02_jprb,2.583499e-02_jprb,2.443266e-02_jprb,2.316681e-02_jprb, &
       2.201687e-02_jprb,2.096619e-02_jprb,2.000112e-02_jprb,1.911044e-02_jprb,1.828481e-02_jprb, &
       1.751641e-02_jprb,1.679866e-02_jprb,1.612598e-02_jprb,1.549360e-02_jprb,1.489742e-02_jprb, &
       1.433392e-02_jprb,1.380002e-02_jprb,1.329305e-02_jprb,1.281068e-02_jprb,1.235084e-02_jprb, &
       1.191172e-02_jprb,1.149171e-02_jprb,1.108936e-02_jprb,1.070341e-02_jprb,1.033271e-02_jprb, &
       9.976220e-03_jprb,9.633021e-03_jprb,9.302273e-03_jprb,8.983216e-03_jprb,8.675161e-03_jprb, &
       8.377478e-03_jprb,8.089595e-03_jprb,7.810986e-03_jprb,7.541170e-03_jprb,7.279706e-03_jprb, &
       7.026186e-03_jprb/)
      absice3(:,8) = (/ &
! band 8
       6.711058e-02_jprb,6.918198e-02_jprb,6.127484e-02_jprb,5.411944e-02_jprb,4.836902e-02_jprb, &
       4.375293e-02_jprb,3.998077e-02_jprb,3.683587e-02_jprb,3.416508e-02_jprb,3.186003e-02_jprb, &
       2.984290e-02_jprb,2.805671e-02_jprb,2.645895e-02_jprb,2.501733e-02_jprb,2.370689e-02_jprb, &
       2.250808e-02_jprb,2.140532e-02_jprb,2.038609e-02_jprb,1.944018e-02_jprb,1.855918e-02_jprb, &
       1.773609e-02_jprb,1.696504e-02_jprb,1.624106e-02_jprb,1.555990e-02_jprb,1.491793e-02_jprb, &
       1.431197e-02_jprb,1.373928e-02_jprb,1.319743e-02_jprb,1.268430e-02_jprb,1.219799e-02_jprb, &
       1.173682e-02_jprb,1.129925e-02_jprb,1.088393e-02_jprb,1.048961e-02_jprb,1.011516e-02_jprb, &
       9.759543e-03_jprb,9.421813e-03_jprb,9.101089e-03_jprb,8.796559e-03_jprb,8.507464e-03_jprb, &
       8.233098e-03_jprb,7.972798e-03_jprb,7.725942e-03_jprb,7.491940e-03_jprb,7.270238e-03_jprb, &
       7.060305e-03_jprb/)
      absice3(:,9) = (/ &
! band 9
       1.236780e-01_jprb,9.222386e-02_jprb,7.383997e-02_jprb,6.204072e-02_jprb,5.381029e-02_jprb, &
       4.770678e-02_jprb,4.296928e-02_jprb,3.916131e-02_jprb,3.601540e-02_jprb,3.335878e-02_jprb, &
       3.107493e-02_jprb,2.908247e-02_jprb,2.732282e-02_jprb,2.575276e-02_jprb,2.433968e-02_jprb, &
       2.305852e-02_jprb,2.188966e-02_jprb,2.081757e-02_jprb,1.982974e-02_jprb,1.891599e-02_jprb, &
       1.806794e-02_jprb,1.727865e-02_jprb,1.654227e-02_jprb,1.585387e-02_jprb,1.520924e-02_jprb, &
       1.460476e-02_jprb,1.403730e-02_jprb,1.350416e-02_jprb,1.300293e-02_jprb,1.253153e-02_jprb, &
       1.208808e-02_jprb,1.167094e-02_jprb,1.127862e-02_jprb,1.090979e-02_jprb,1.056323e-02_jprb, &
       1.023786e-02_jprb,9.932665e-03_jprb,9.646744e-03_jprb,9.379250e-03_jprb,9.129409e-03_jprb, &
       8.896500e-03_jprb,8.679856e-03_jprb,8.478852e-03_jprb,8.292904e-03_jprb,8.121463e-03_jprb, &
       7.964013e-03_jprb/)
      absice3(:,10) = (/ &
! band 10
       1.655966e-01_jprb,1.134205e-01_jprb,8.714344e-02_jprb,7.129241e-02_jprb,6.063739e-02_jprb, &
       5.294203e-02_jprb,4.709309e-02_jprb,4.247476e-02_jprb,3.871892e-02_jprb,3.559206e-02_jprb, &
       3.293893e-02_jprb,3.065226e-02_jprb,2.865558e-02_jprb,2.689288e-02_jprb,2.532221e-02_jprb, &
       2.391150e-02_jprb,2.263582e-02_jprb,2.147549e-02_jprb,2.041476e-02_jprb,1.944089e-02_jprb, &
       1.854342e-02_jprb,1.771371e-02_jprb,1.694456e-02_jprb,1.622989e-02_jprb,1.556456e-02_jprb, &
       1.494415e-02_jprb,1.436491e-02_jprb,1.382354e-02_jprb,1.331719e-02_jprb,1.284339e-02_jprb, &
       1.239992e-02_jprb,1.198486e-02_jprb,1.159647e-02_jprb,1.123323e-02_jprb,1.089375e-02_jprb, &
       1.057679e-02_jprb,1.028124e-02_jprb,1.000607e-02_jprb,9.750376e-03_jprb,9.513303e-03_jprb, &
       9.294082e-03_jprb,9.092003e-03_jprb,8.906412e-03_jprb,8.736702e-03_jprb,8.582314e-03_jprb, &
       8.442725e-03_jprb/)
      absice3(:,11) = (/ &
! band 11
       1.775615e-01_jprb,1.180046e-01_jprb,8.929607e-02_jprb,7.233500e-02_jprb,6.108333e-02_jprb, &
       5.303642e-02_jprb,4.696927e-02_jprb,4.221206e-02_jprb,3.836768e-02_jprb,3.518576e-02_jprb, &
       3.250063e-02_jprb,3.019825e-02_jprb,2.819758e-02_jprb,2.643943e-02_jprb,2.487953e-02_jprb, &
       2.348414e-02_jprb,2.222705e-02_jprb,2.108762e-02_jprb,2.004936e-02_jprb,1.909892e-02_jprb, &
       1.822539e-02_jprb,1.741975e-02_jprb,1.667449e-02_jprb,1.598330e-02_jprb,1.534084e-02_jprb, &
       1.474253e-02_jprb,1.418446e-02_jprb,1.366325e-02_jprb,1.317597e-02_jprb,1.272004e-02_jprb, &
       1.229321e-02_jprb,1.189350e-02_jprb,1.151915e-02_jprb,1.116859e-02_jprb,1.084042e-02_jprb, &
       1.053338e-02_jprb,1.024636e-02_jprb,9.978326e-03_jprb,9.728357e-03_jprb,9.495613e-03_jprb, &
       9.279327e-03_jprb,9.078798e-03_jprb,8.893383e-03_jprb,8.722488e-03_jprb,8.565568e-03_jprb, &
       8.422115e-03_jprb/)
      absice3(:,12) = (/ &
! band 12
       9.465447e-02_jprb,6.432047e-02_jprb,5.060973e-02_jprb,4.267283e-02_jprb,3.741843e-02_jprb, &
       3.363096e-02_jprb,3.073531e-02_jprb,2.842405e-02_jprb,2.651789e-02_jprb,2.490518e-02_jprb, &
       2.351273e-02_jprb,2.229056e-02_jprb,2.120335e-02_jprb,2.022541e-02_jprb,1.933763e-02_jprb, &
       1.852546e-02_jprb,1.777763e-02_jprb,1.708528e-02_jprb,1.644134e-02_jprb,1.584009e-02_jprb, &
       1.527684e-02_jprb,1.474774e-02_jprb,1.424955e-02_jprb,1.377957e-02_jprb,1.333549e-02_jprb, &
       1.291534e-02_jprb,1.251743e-02_jprb,1.214029e-02_jprb,1.178265e-02_jprb,1.144337e-02_jprb, &
       1.112148e-02_jprb,1.081609e-02_jprb,1.052642e-02_jprb,1.025178e-02_jprb,9.991540e-03_jprb, &
       9.745130e-03_jprb,9.512038e-03_jprb,9.291797e-03_jprb,9.083980e-03_jprb,8.888195e-03_jprb, &
       8.704081e-03_jprb,8.531306e-03_jprb,8.369560e-03_jprb,8.218558e-03_jprb,8.078032e-03_jprb, &
       7.947730e-03_jprb/)
      absice3(:,13) = (/ &
! band 13
       1.560311e-01_jprb,9.961097e-02_jprb,7.502949e-02_jprb,6.115022e-02_jprb,5.214952e-02_jprb, &
       4.578149e-02_jprb,4.099731e-02_jprb,3.724174e-02_jprb,3.419343e-02_jprb,3.165356e-02_jprb, &
       2.949251e-02_jprb,2.762222e-02_jprb,2.598073e-02_jprb,2.452322e-02_jprb,2.321642e-02_jprb, &
       2.203516e-02_jprb,2.096002e-02_jprb,1.997579e-02_jprb,1.907036e-02_jprb,1.823401e-02_jprb, &
       1.745879e-02_jprb,1.673819e-02_jprb,1.606678e-02_jprb,1.544003e-02_jprb,1.485411e-02_jprb, &
       1.430574e-02_jprb,1.379215e-02_jprb,1.331092e-02_jprb,1.285996e-02_jprb,1.243746e-02_jprb, &
       1.204183e-02_jprb,1.167164e-02_jprb,1.132567e-02_jprb,1.100281e-02_jprb,1.070207e-02_jprb, &
       1.042258e-02_jprb,1.016352e-02_jprb,9.924197e-03_jprb,9.703953e-03_jprb,9.502199e-03_jprb, &
       9.318400e-03_jprb,9.152066e-03_jprb,9.002749e-03_jprb,8.870038e-03_jprb,8.753555e-03_jprb, &
       8.652951e-03_jprb/)
      absice3(:,14) = (/ &
! band 14
       1.559547e-01_jprb,9.896700e-02_jprb,7.441231e-02_jprb,6.061469e-02_jprb,5.168730e-02_jprb, &
       4.537821e-02_jprb,4.064106e-02_jprb,3.692367e-02_jprb,3.390714e-02_jprb,3.139438e-02_jprb, &
       2.925702e-02_jprb,2.740783e-02_jprb,2.578547e-02_jprb,2.434552e-02_jprb,2.305506e-02_jprb, &
       2.188910e-02_jprb,2.082842e-02_jprb,1.985789e-02_jprb,1.896553e-02_jprb,1.814165e-02_jprb, &
       1.737839e-02_jprb,1.666927e-02_jprb,1.600891e-02_jprb,1.539279e-02_jprb,1.481712e-02_jprb, &
       1.427865e-02_jprb,1.377463e-02_jprb,1.330266e-02_jprb,1.286068e-02_jprb,1.244689e-02_jprb, &
       1.205973e-02_jprb,1.169780e-02_jprb,1.135989e-02_jprb,1.104492e-02_jprb,1.075192e-02_jprb, &
       1.048004e-02_jprb,1.022850e-02_jprb,9.996611e-03_jprb,9.783753e-03_jprb,9.589361e-03_jprb, &
       9.412924e-03_jprb,9.253977e-03_jprb,9.112098e-03_jprb,8.986903e-03_jprb,8.878039e-03_jprb, &
       8.785184e-03_jprb/)
      absice3(:,15) = (/ &
! band 15
       1.102926e-01_jprb,7.176622e-02_jprb,5.530316e-02_jprb,4.606056e-02_jprb,4.006116e-02_jprb, &
       3.579628e-02_jprb,3.256909e-02_jprb,3.001360e-02_jprb,2.791920e-02_jprb,2.615617e-02_jprb, &
       2.464023e-02_jprb,2.331426e-02_jprb,2.213817e-02_jprb,2.108301e-02_jprb,2.012733e-02_jprb, &
       1.925493e-02_jprb,1.845331e-02_jprb,1.771269e-02_jprb,1.702531e-02_jprb,1.638493e-02_jprb, &
       1.578648e-02_jprb,1.522579e-02_jprb,1.469940e-02_jprb,1.420442e-02_jprb,1.373841e-02_jprb, &
       1.329931e-02_jprb,1.288535e-02_jprb,1.249502e-02_jprb,1.212700e-02_jprb,1.178015e-02_jprb, &
       1.145348e-02_jprb,1.114612e-02_jprb,1.085730e-02_jprb,1.058633e-02_jprb,1.033263e-02_jprb, &
       1.009564e-02_jprb,9.874895e-03_jprb,9.669960e-03_jprb,9.480449e-03_jprb,9.306014e-03_jprb, &
       9.146339e-03_jprb,9.001138e-03_jprb,8.870154e-03_jprb,8.753148e-03_jprb,8.649907e-03_jprb, &
       8.560232e-03_jprb/)
      absice3(:,16) = (/ &
! band 16
       1.688344e-01_jprb,1.077072e-01_jprb,7.994467e-02_jprb,6.403862e-02_jprb,5.369850e-02_jprb, &
       4.641582e-02_jprb,4.099331e-02_jprb,3.678724e-02_jprb,3.342069e-02_jprb,3.065831e-02_jprb, &
       2.834557e-02_jprb,2.637680e-02_jprb,2.467733e-02_jprb,2.319286e-02_jprb,2.188299e-02_jprb, &
       2.071701e-02_jprb,1.967121e-02_jprb,1.872692e-02_jprb,1.786931e-02_jprb,1.708641e-02_jprb, &
       1.636846e-02_jprb,1.570743e-02_jprb,1.509665e-02_jprb,1.453052e-02_jprb,1.400433e-02_jprb, &
       1.351407e-02_jprb,1.305631e-02_jprb,1.262810e-02_jprb,1.222688e-02_jprb,1.185044e-02_jprb, &
       1.149683e-02_jprb,1.116436e-02_jprb,1.085153e-02_jprb,1.055701e-02_jprb,1.027961e-02_jprb, &
       1.001831e-02_jprb,9.772141e-03_jprb,9.540280e-03_jprb,9.321966e-03_jprb,9.116517e-03_jprb, &
       8.923315e-03_jprb,8.741803e-03_jprb,8.571472e-03_jprb,8.411860e-03_jprb,8.262543e-03_jprb, &
       8.123136e-03_jprb/)

! For LIQFLAG = 0.
      absliq0 = 0.0903614_jprb

! For LIQFLAG = 1.  In each band, the absorption
! coefficients are listed for a range of effective radii from 2.5
! to 59.5 microns in increments of 1.0 micron.
      absliq1(:, 1) = (/ &
! band  1
       1.64047e-03_jprb, 6.90533e-02_jprb, 7.72017e-02_jprb, 7.78054e-02_jprb, 7.69523e-02_jprb, &
       7.58058e-02_jprb, 7.46400e-02_jprb, 7.35123e-02_jprb, 7.24162e-02_jprb, 7.13225e-02_jprb, &
       6.99145e-02_jprb, 6.66409e-02_jprb, 6.36582e-02_jprb, 6.09425e-02_jprb, 5.84593e-02_jprb, &
       5.61743e-02_jprb, 5.40571e-02_jprb, 5.20812e-02_jprb, 5.02245e-02_jprb, 4.84680e-02_jprb, &
       4.67959e-02_jprb, 4.51944e-02_jprb, 4.36516e-02_jprb, 4.21570e-02_jprb, 4.07015e-02_jprb, &
       3.92766e-02_jprb, 3.78747e-02_jprb, 3.64886e-02_jprb, 3.53632e-02_jprb, 3.41992e-02_jprb, &
       3.31016e-02_jprb, 3.20643e-02_jprb, 3.10817e-02_jprb, 3.01490e-02_jprb, 2.92620e-02_jprb, &
       2.84171e-02_jprb, 2.76108e-02_jprb, 2.68404e-02_jprb, 2.61031e-02_jprb, 2.53966e-02_jprb, &
       2.47189e-02_jprb, 2.40678e-02_jprb, 2.34418e-02_jprb, 2.28392e-02_jprb, 2.22586e-02_jprb, &
       2.16986e-02_jprb, 2.11580e-02_jprb, 2.06356e-02_jprb, 2.01305e-02_jprb, 1.96417e-02_jprb, &
       1.91682e-02_jprb, 1.87094e-02_jprb, 1.82643e-02_jprb, 1.78324e-02_jprb, 1.74129e-02_jprb, &
       1.70052e-02_jprb, 1.66088e-02_jprb, 1.62231e-02_jprb/)
      absliq1(:, 2) = (/ &
! band  2
       2.19486e-01_jprb, 1.80687e-01_jprb, 1.59150e-01_jprb, 1.44731e-01_jprb, 1.33703e-01_jprb, &
       1.24355e-01_jprb, 1.15756e-01_jprb, 1.07318e-01_jprb, 9.86119e-02_jprb, 8.92739e-02_jprb, &
       8.34911e-02_jprb, 7.70773e-02_jprb, 7.15240e-02_jprb, 6.66615e-02_jprb, 6.23641e-02_jprb, &
       5.85359e-02_jprb, 5.51020e-02_jprb, 5.20032e-02_jprb, 4.91916e-02_jprb, 4.66283e-02_jprb, &
       4.42813e-02_jprb, 4.21236e-02_jprb, 4.01330e-02_jprb, 3.82905e-02_jprb, 3.65797e-02_jprb, &
       3.49869e-02_jprb, 3.35002e-02_jprb, 3.21090e-02_jprb, 3.08957e-02_jprb, 2.97601e-02_jprb, &
       2.86966e-02_jprb, 2.76984e-02_jprb, 2.67599e-02_jprb, 2.58758e-02_jprb, 2.50416e-02_jprb, &
       2.42532e-02_jprb, 2.35070e-02_jprb, 2.27997e-02_jprb, 2.21284e-02_jprb, 2.14904e-02_jprb, &
       2.08834e-02_jprb, 2.03051e-02_jprb, 1.97536e-02_jprb, 1.92271e-02_jprb, 1.87239e-02_jprb, &
       1.82425e-02_jprb, 1.77816e-02_jprb, 1.73399e-02_jprb, 1.69162e-02_jprb, 1.65094e-02_jprb, &
       1.61187e-02_jprb, 1.57430e-02_jprb, 1.53815e-02_jprb, 1.50334e-02_jprb, 1.46981e-02_jprb, &
       1.43748e-02_jprb, 1.40628e-02_jprb, 1.37617e-02_jprb/)
      absliq1(:, 3) = (/ &
! band  3
       2.95174e-01_jprb, 2.34765e-01_jprb, 1.98038e-01_jprb, 1.72114e-01_jprb, 1.52083e-01_jprb, &
       1.35654e-01_jprb, 1.21613e-01_jprb, 1.09252e-01_jprb, 9.81263e-02_jprb, 8.79448e-02_jprb, &
       8.12566e-02_jprb, 7.44563e-02_jprb, 6.86374e-02_jprb, 6.36042e-02_jprb, 5.92094e-02_jprb, &
       5.53402e-02_jprb, 5.19087e-02_jprb, 4.88455e-02_jprb, 4.60951e-02_jprb, 4.36124e-02_jprb, &
       4.13607e-02_jprb, 3.93096e-02_jprb, 3.74338e-02_jprb, 3.57119e-02_jprb, 3.41261e-02_jprb, &
       3.26610e-02_jprb, 3.13036e-02_jprb, 3.00425e-02_jprb, 2.88497e-02_jprb, 2.78077e-02_jprb, &
       2.68317e-02_jprb, 2.59158e-02_jprb, 2.50545e-02_jprb, 2.42430e-02_jprb, 2.34772e-02_jprb, &
       2.27533e-02_jprb, 2.20679e-02_jprb, 2.14181e-02_jprb, 2.08011e-02_jprb, 2.02145e-02_jprb, &
       1.96561e-02_jprb, 1.91239e-02_jprb, 1.86161e-02_jprb, 1.81311e-02_jprb, 1.76673e-02_jprb, &
       1.72234e-02_jprb, 1.67981e-02_jprb, 1.63903e-02_jprb, 1.59989e-02_jprb, 1.56230e-02_jprb, &
       1.52615e-02_jprb, 1.49138e-02_jprb, 1.45791e-02_jprb, 1.42565e-02_jprb, 1.39455e-02_jprb, &
       1.36455e-02_jprb, 1.33559e-02_jprb, 1.30761e-02_jprb/)
      absliq1(:, 4) = (/ &
! band  4
       3.00925e-01_jprb, 2.36949e-01_jprb, 1.96947e-01_jprb, 1.68692e-01_jprb, 1.47190e-01_jprb, &
       1.29986e-01_jprb, 1.15719e-01_jprb, 1.03568e-01_jprb, 9.30028e-02_jprb, 8.36658e-02_jprb, &
       7.71075e-02_jprb, 7.07002e-02_jprb, 6.52284e-02_jprb, 6.05024e-02_jprb, 5.63801e-02_jprb, &
       5.27534e-02_jprb, 4.95384e-02_jprb, 4.66690e-02_jprb, 4.40925e-02_jprb, 4.17664e-02_jprb, &
       3.96559e-02_jprb, 3.77326e-02_jprb, 3.59727e-02_jprb, 3.43561e-02_jprb, 3.28662e-02_jprb, &
       3.14885e-02_jprb, 3.02110e-02_jprb, 2.90231e-02_jprb, 2.78948e-02_jprb, 2.69109e-02_jprb, &
       2.59884e-02_jprb, 2.51217e-02_jprb, 2.43058e-02_jprb, 2.35364e-02_jprb, 2.28096e-02_jprb, &
       2.21218e-02_jprb, 2.14700e-02_jprb, 2.08515e-02_jprb, 2.02636e-02_jprb, 1.97041e-02_jprb, &
       1.91711e-02_jprb, 1.86625e-02_jprb, 1.81769e-02_jprb, 1.77126e-02_jprb, 1.72683e-02_jprb, &
       1.68426e-02_jprb, 1.64344e-02_jprb, 1.60427e-02_jprb, 1.56664e-02_jprb, 1.53046e-02_jprb, &
       1.49565e-02_jprb, 1.46214e-02_jprb, 1.42985e-02_jprb, 1.39871e-02_jprb, 1.36866e-02_jprb, &
       1.33965e-02_jprb, 1.31162e-02_jprb, 1.28453e-02_jprb/)
      absliq1(:, 5) = (/ &
! band  5
       2.64691e-01_jprb, 2.12018e-01_jprb, 1.78009e-01_jprb, 1.53539e-01_jprb, 1.34721e-01_jprb, &
       1.19580e-01_jprb, 1.06996e-01_jprb, 9.62772e-02_jprb, 8.69710e-02_jprb, 7.87670e-02_jprb, &
       7.29272e-02_jprb, 6.70920e-02_jprb, 6.20977e-02_jprb, 5.77732e-02_jprb, 5.39910e-02_jprb, &
       5.06538e-02_jprb, 4.76866e-02_jprb, 4.50301e-02_jprb, 4.26374e-02_jprb, 4.04704e-02_jprb, &
       3.84981e-02_jprb, 3.66948e-02_jprb, 3.50394e-02_jprb, 3.35141e-02_jprb, 3.21038e-02_jprb, &
       3.07957e-02_jprb, 2.95788e-02_jprb, 2.84438e-02_jprb, 2.73790e-02_jprb, 2.64390e-02_jprb, &
       2.55565e-02_jprb, 2.47263e-02_jprb, 2.39437e-02_jprb, 2.32047e-02_jprb, 2.25056e-02_jprb, &
       2.18433e-02_jprb, 2.12149e-02_jprb, 2.06177e-02_jprb, 2.00495e-02_jprb, 1.95081e-02_jprb, &
       1.89917e-02_jprb, 1.84984e-02_jprb, 1.80269e-02_jprb, 1.75755e-02_jprb, 1.71431e-02_jprb, &
       1.67283e-02_jprb, 1.63303e-02_jprb, 1.59478e-02_jprb, 1.55801e-02_jprb, 1.52262e-02_jprb, &
       1.48853e-02_jprb, 1.45568e-02_jprb, 1.42400e-02_jprb, 1.39342e-02_jprb, 1.36388e-02_jprb, &
       1.33533e-02_jprb, 1.30773e-02_jprb, 1.28102e-02_jprb/)
      absliq1(:, 6) = (/ &
! band  6
       8.81182e-02_jprb, 1.06745e-01_jprb, 9.79753e-02_jprb, 8.99625e-02_jprb, 8.35200e-02_jprb, &
       7.81899e-02_jprb, 7.35939e-02_jprb, 6.94696e-02_jprb, 6.56266e-02_jprb, 6.19148e-02_jprb, &
       5.83355e-02_jprb, 5.49306e-02_jprb, 5.19642e-02_jprb, 4.93325e-02_jprb, 4.69659e-02_jprb, &
       4.48148e-02_jprb, 4.28431e-02_jprb, 4.10231e-02_jprb, 3.93332e-02_jprb, 3.77563e-02_jprb, &
       3.62785e-02_jprb, 3.48882e-02_jprb, 3.35758e-02_jprb, 3.23333e-02_jprb, 3.11536e-02_jprb, &
       3.00310e-02_jprb, 2.89601e-02_jprb, 2.79365e-02_jprb, 2.70502e-02_jprb, 2.62618e-02_jprb, &
       2.55025e-02_jprb, 2.47728e-02_jprb, 2.40726e-02_jprb, 2.34013e-02_jprb, 2.27583e-02_jprb, &
       2.21422e-02_jprb, 2.15522e-02_jprb, 2.09869e-02_jprb, 2.04453e-02_jprb, 1.99260e-02_jprb, &
       1.94280e-02_jprb, 1.89501e-02_jprb, 1.84913e-02_jprb, 1.80506e-02_jprb, 1.76270e-02_jprb, &
       1.72196e-02_jprb, 1.68276e-02_jprb, 1.64500e-02_jprb, 1.60863e-02_jprb, 1.57357e-02_jprb, &
       1.53975e-02_jprb, 1.50710e-02_jprb, 1.47558e-02_jprb, 1.44511e-02_jprb, 1.41566e-02_jprb, &
       1.38717e-02_jprb, 1.35960e-02_jprb, 1.33290e-02_jprb/)
      absliq1(:, 7) = (/ &
! band  7
       4.32174e-02_jprb, 7.36078e-02_jprb, 6.98340e-02_jprb, 6.65231e-02_jprb, 6.41948e-02_jprb, &
       6.23551e-02_jprb, 6.06638e-02_jprb, 5.88680e-02_jprb, 5.67124e-02_jprb, 5.38629e-02_jprb, &
       4.99579e-02_jprb, 4.86289e-02_jprb, 4.70120e-02_jprb, 4.52854e-02_jprb, 4.35466e-02_jprb, &
       4.18480e-02_jprb, 4.02169e-02_jprb, 3.86658e-02_jprb, 3.71992e-02_jprb, 3.58168e-02_jprb, &
       3.45155e-02_jprb, 3.32912e-02_jprb, 3.21390e-02_jprb, 3.10538e-02_jprb, 3.00307e-02_jprb, &
       2.90651e-02_jprb, 2.81524e-02_jprb, 2.72885e-02_jprb, 2.62821e-02_jprb, 2.55744e-02_jprb, &
       2.48799e-02_jprb, 2.42029e-02_jprb, 2.35460e-02_jprb, 2.29108e-02_jprb, 2.22981e-02_jprb, &
       2.17079e-02_jprb, 2.11402e-02_jprb, 2.05945e-02_jprb, 2.00701e-02_jprb, 1.95663e-02_jprb, &
       1.90824e-02_jprb, 1.86174e-02_jprb, 1.81706e-02_jprb, 1.77411e-02_jprb, 1.73281e-02_jprb, &
       1.69307e-02_jprb, 1.65483e-02_jprb, 1.61801e-02_jprb, 1.58254e-02_jprb, 1.54835e-02_jprb, &
       1.51538e-02_jprb, 1.48358e-02_jprb, 1.45288e-02_jprb, 1.42322e-02_jprb, 1.39457e-02_jprb, &
       1.36687e-02_jprb, 1.34008e-02_jprb, 1.31416e-02_jprb/)
      absliq1(:, 8) = (/ &
! band  8
       1.41881e-01_jprb, 7.15419e-02_jprb, 6.30335e-02_jprb, 6.11132e-02_jprb, 6.01931e-02_jprb, &
       5.92420e-02_jprb, 5.78968e-02_jprb, 5.58876e-02_jprb, 5.28923e-02_jprb, 4.84462e-02_jprb, &
       4.60839e-02_jprb, 4.56013e-02_jprb, 4.45410e-02_jprb, 4.31866e-02_jprb, 4.17026e-02_jprb, &
       4.01850e-02_jprb, 3.86892e-02_jprb, 3.72461e-02_jprb, 3.58722e-02_jprb, 3.45749e-02_jprb, &
       3.33564e-02_jprb, 3.22155e-02_jprb, 3.11494e-02_jprb, 3.01541e-02_jprb, 2.92253e-02_jprb, &
       2.83584e-02_jprb, 2.75488e-02_jprb, 2.67925e-02_jprb, 2.57692e-02_jprb, 2.50704e-02_jprb, &
       2.43918e-02_jprb, 2.37350e-02_jprb, 2.31005e-02_jprb, 2.24888e-02_jprb, 2.18996e-02_jprb, &
       2.13325e-02_jprb, 2.07870e-02_jprb, 2.02623e-02_jprb, 1.97577e-02_jprb, 1.92724e-02_jprb, &
       1.88056e-02_jprb, 1.83564e-02_jprb, 1.79241e-02_jprb, 1.75079e-02_jprb, 1.71070e-02_jprb, &
       1.67207e-02_jprb, 1.63482e-02_jprb, 1.59890e-02_jprb, 1.56424e-02_jprb, 1.53077e-02_jprb, &
       1.49845e-02_jprb, 1.46722e-02_jprb, 1.43702e-02_jprb, 1.40782e-02_jprb, 1.37955e-02_jprb, &
       1.35219e-02_jprb, 1.32569e-02_jprb, 1.30000e-02_jprb/)
      absliq1(:, 9) = (/ &
! band  9
       6.72726e-02_jprb, 6.61013e-02_jprb, 6.47866e-02_jprb, 6.33780e-02_jprb, 6.18985e-02_jprb, &
       6.03335e-02_jprb, 5.86136e-02_jprb, 5.65876e-02_jprb, 5.39839e-02_jprb, 5.03536e-02_jprb, &
       4.71608e-02_jprb, 4.63630e-02_jprb, 4.50313e-02_jprb, 4.34526e-02_jprb, 4.17876e-02_jprb, &
       4.01261e-02_jprb, 3.85171e-02_jprb, 3.69860e-02_jprb, 3.55442e-02_jprb, 3.41954e-02_jprb, &
       3.29384e-02_jprb, 3.17693e-02_jprb, 3.06832e-02_jprb, 2.96745e-02_jprb, 2.87374e-02_jprb, &
       2.78662e-02_jprb, 2.70557e-02_jprb, 2.63008e-02_jprb, 2.52450e-02_jprb, 2.45424e-02_jprb, &
       2.38656e-02_jprb, 2.32144e-02_jprb, 2.25885e-02_jprb, 2.19873e-02_jprb, 2.14099e-02_jprb, &
       2.08554e-02_jprb, 2.03230e-02_jprb, 1.98116e-02_jprb, 1.93203e-02_jprb, 1.88482e-02_jprb, &
       1.83944e-02_jprb, 1.79578e-02_jprb, 1.75378e-02_jprb, 1.71335e-02_jprb, 1.67440e-02_jprb, &
       1.63687e-02_jprb, 1.60069e-02_jprb, 1.56579e-02_jprb, 1.53210e-02_jprb, 1.49958e-02_jprb, &
       1.46815e-02_jprb, 1.43778e-02_jprb, 1.40841e-02_jprb, 1.37999e-02_jprb, 1.35249e-02_jprb, &
       1.32585e-02_jprb, 1.30004e-02_jprb, 1.27502e-02_jprb/)
      absliq1(:,10) = (/ &
! band 10
       7.97040e-02_jprb, 7.63844e-02_jprb, 7.36499e-02_jprb, 7.13525e-02_jprb, 6.93043e-02_jprb, &
       6.72807e-02_jprb, 6.50227e-02_jprb, 6.22395e-02_jprb, 5.86093e-02_jprb, 5.37815e-02_jprb, &
       5.14682e-02_jprb, 4.97214e-02_jprb, 4.77392e-02_jprb, 4.56961e-02_jprb, 4.36858e-02_jprb, &
       4.17569e-02_jprb, 3.99328e-02_jprb, 3.82224e-02_jprb, 3.66265e-02_jprb, 3.51416e-02_jprb, &
       3.37617e-02_jprb, 3.24798e-02_jprb, 3.12887e-02_jprb, 3.01812e-02_jprb, 2.91505e-02_jprb, &
       2.81900e-02_jprb, 2.72939e-02_jprb, 2.64568e-02_jprb, 2.54165e-02_jprb, 2.46832e-02_jprb, &
       2.39783e-02_jprb, 2.33017e-02_jprb, 2.26531e-02_jprb, 2.20314e-02_jprb, 2.14359e-02_jprb, &
       2.08653e-02_jprb, 2.03187e-02_jprb, 1.97947e-02_jprb, 1.92924e-02_jprb, 1.88106e-02_jprb, &
       1.83483e-02_jprb, 1.79043e-02_jprb, 1.74778e-02_jprb, 1.70678e-02_jprb, 1.66735e-02_jprb, &
       1.62941e-02_jprb, 1.59286e-02_jprb, 1.55766e-02_jprb, 1.52371e-02_jprb, 1.49097e-02_jprb, &
       1.45937e-02_jprb, 1.42885e-02_jprb, 1.39936e-02_jprb, 1.37085e-02_jprb, 1.34327e-02_jprb, &
       1.31659e-02_jprb, 1.29075e-02_jprb, 1.26571e-02_jprb/)
      absliq1(:,11) = (/ &
! band 11
       1.49438e-01_jprb, 1.33535e-01_jprb, 1.21542e-01_jprb, 1.11743e-01_jprb, 1.03263e-01_jprb, &
       9.55774e-02_jprb, 8.83382e-02_jprb, 8.12943e-02_jprb, 7.42533e-02_jprb, 6.70609e-02_jprb, &
       6.38761e-02_jprb, 5.97788e-02_jprb, 5.59841e-02_jprb, 5.25318e-02_jprb, 4.94132e-02_jprb, &
       4.66014e-02_jprb, 4.40644e-02_jprb, 4.17706e-02_jprb, 3.96910e-02_jprb, 3.77998e-02_jprb, &
       3.60742e-02_jprb, 3.44947e-02_jprb, 3.30442e-02_jprb, 3.17079e-02_jprb, 3.04730e-02_jprb, &
       2.93283e-02_jprb, 2.82642e-02_jprb, 2.72720e-02_jprb, 2.61789e-02_jprb, 2.53277e-02_jprb, &
       2.45237e-02_jprb, 2.37635e-02_jprb, 2.30438e-02_jprb, 2.23615e-02_jprb, 2.17140e-02_jprb, &
       2.10987e-02_jprb, 2.05133e-02_jprb, 1.99557e-02_jprb, 1.94241e-02_jprb, 1.89166e-02_jprb, &
       1.84317e-02_jprb, 1.79679e-02_jprb, 1.75238e-02_jprb, 1.70983e-02_jprb, 1.66901e-02_jprb, &
       1.62983e-02_jprb, 1.59219e-02_jprb, 1.55599e-02_jprb, 1.52115e-02_jprb, 1.48761e-02_jprb, &
       1.45528e-02_jprb, 1.42411e-02_jprb, 1.39402e-02_jprb, 1.36497e-02_jprb, 1.33690e-02_jprb, &
       1.30976e-02_jprb, 1.28351e-02_jprb, 1.25810e-02_jprb/)
      absliq1(:,12) = (/ &
! band 12
       3.71985e-02_jprb, 3.88586e-02_jprb, 3.99070e-02_jprb, 4.04351e-02_jprb, 4.04610e-02_jprb, &
       3.99834e-02_jprb, 3.89953e-02_jprb, 3.74886e-02_jprb, 3.54551e-02_jprb, 3.28870e-02_jprb, &
       3.32576e-02_jprb, 3.22444e-02_jprb, 3.12384e-02_jprb, 3.02584e-02_jprb, 2.93146e-02_jprb, &
       2.84120e-02_jprb, 2.75525e-02_jprb, 2.67361e-02_jprb, 2.59618e-02_jprb, 2.52280e-02_jprb, &
       2.45327e-02_jprb, 2.38736e-02_jprb, 2.32487e-02_jprb, 2.26558e-02_jprb, 2.20929e-02_jprb, &
       2.15579e-02_jprb, 2.10491e-02_jprb, 2.05648e-02_jprb, 1.99749e-02_jprb, 1.95704e-02_jprb, &
       1.91731e-02_jprb, 1.87839e-02_jprb, 1.84032e-02_jprb, 1.80315e-02_jprb, 1.76689e-02_jprb, &
       1.73155e-02_jprb, 1.69712e-02_jprb, 1.66362e-02_jprb, 1.63101e-02_jprb, 1.59928e-02_jprb, &
       1.56842e-02_jprb, 1.53840e-02_jprb, 1.50920e-02_jprb, 1.48080e-02_jprb, 1.45318e-02_jprb, &
       1.42631e-02_jprb, 1.40016e-02_jprb, 1.37472e-02_jprb, 1.34996e-02_jprb, 1.32586e-02_jprb, &
       1.30239e-02_jprb, 1.27954e-02_jprb, 1.25728e-02_jprb, 1.23559e-02_jprb, 1.21445e-02_jprb, &
       1.19385e-02_jprb, 1.17376e-02_jprb, 1.15417e-02_jprb/)
      absliq1(:,13) = (/ &
! band 13
       3.11868e-02_jprb, 4.48357e-02_jprb, 4.90224e-02_jprb, 4.96406e-02_jprb, 4.86806e-02_jprb, &
       4.69610e-02_jprb, 4.48630e-02_jprb, 4.25795e-02_jprb, 4.02138e-02_jprb, 3.78236e-02_jprb, &
       3.74266e-02_jprb, 3.60384e-02_jprb, 3.47074e-02_jprb, 3.34434e-02_jprb, 3.22499e-02_jprb, &
       3.11264e-02_jprb, 3.00704e-02_jprb, 2.90784e-02_jprb, 2.81463e-02_jprb, 2.72702e-02_jprb, &
       2.64460e-02_jprb, 2.56698e-02_jprb, 2.49381e-02_jprb, 2.42475e-02_jprb, 2.35948e-02_jprb, &
       2.29774e-02_jprb, 2.23925e-02_jprb, 2.18379e-02_jprb, 2.11793e-02_jprb, 2.07076e-02_jprb, &
       2.02470e-02_jprb, 1.97981e-02_jprb, 1.93613e-02_jprb, 1.89367e-02_jprb, 1.85243e-02_jprb, &
       1.81240e-02_jprb, 1.77356e-02_jprb, 1.73588e-02_jprb, 1.69935e-02_jprb, 1.66392e-02_jprb, &
       1.62956e-02_jprb, 1.59624e-02_jprb, 1.56393e-02_jprb, 1.53259e-02_jprb, 1.50219e-02_jprb, &
       1.47268e-02_jprb, 1.44404e-02_jprb, 1.41624e-02_jprb, 1.38925e-02_jprb, 1.36302e-02_jprb, &
       1.33755e-02_jprb, 1.31278e-02_jprb, 1.28871e-02_jprb, 1.26530e-02_jprb, 1.24253e-02_jprb, &
       1.22038e-02_jprb, 1.19881e-02_jprb, 1.17782e-02_jprb/)
      absliq1(:,14) = (/ &
! band 14
       1.58988e-02_jprb, 3.50652e-02_jprb, 4.00851e-02_jprb, 4.07270e-02_jprb, 3.98101e-02_jprb, &
       3.83306e-02_jprb, 3.66829e-02_jprb, 3.50327e-02_jprb, 3.34497e-02_jprb, 3.19609e-02_jprb, &
       3.13712e-02_jprb, 3.03348e-02_jprb, 2.93415e-02_jprb, 2.83973e-02_jprb, 2.75037e-02_jprb, &
       2.66604e-02_jprb, 2.58654e-02_jprb, 2.51161e-02_jprb, 2.44100e-02_jprb, 2.37440e-02_jprb, &
       2.31154e-02_jprb, 2.25215e-02_jprb, 2.19599e-02_jprb, 2.14282e-02_jprb, 2.09242e-02_jprb, &
       2.04459e-02_jprb, 1.99915e-02_jprb, 1.95594e-02_jprb, 1.90254e-02_jprb, 1.86598e-02_jprb, &
       1.82996e-02_jprb, 1.79455e-02_jprb, 1.75983e-02_jprb, 1.72584e-02_jprb, 1.69260e-02_jprb, &
       1.66013e-02_jprb, 1.62843e-02_jprb, 1.59752e-02_jprb, 1.56737e-02_jprb, 1.53799e-02_jprb, &
       1.50936e-02_jprb, 1.48146e-02_jprb, 1.45429e-02_jprb, 1.42782e-02_jprb, 1.40203e-02_jprb, &
       1.37691e-02_jprb, 1.35243e-02_jprb, 1.32858e-02_jprb, 1.30534e-02_jprb, 1.28270e-02_jprb, &
       1.26062e-02_jprb, 1.23909e-02_jprb, 1.21810e-02_jprb, 1.19763e-02_jprb, 1.17766e-02_jprb, &
       1.15817e-02_jprb, 1.13915e-02_jprb, 1.12058e-02_jprb/)
      absliq1(:,15) = (/ &
! band 15
       5.02079e-03_jprb, 2.17615e-02_jprb, 2.55449e-02_jprb, 2.59484e-02_jprb, 2.53650e-02_jprb, &
       2.45281e-02_jprb, 2.36843e-02_jprb, 2.29159e-02_jprb, 2.22451e-02_jprb, 2.16716e-02_jprb, &
       2.11451e-02_jprb, 2.05817e-02_jprb, 2.00454e-02_jprb, 1.95372e-02_jprb, 1.90567e-02_jprb, &
       1.86028e-02_jprb, 1.81742e-02_jprb, 1.77693e-02_jprb, 1.73866e-02_jprb, 1.70244e-02_jprb, &
       1.66815e-02_jprb, 1.63563e-02_jprb, 1.60477e-02_jprb, 1.57544e-02_jprb, 1.54755e-02_jprb, &
       1.52097e-02_jprb, 1.49564e-02_jprb, 1.47146e-02_jprb, 1.43684e-02_jprb, 1.41728e-02_jprb, &
       1.39762e-02_jprb, 1.37797e-02_jprb, 1.35838e-02_jprb, 1.33891e-02_jprb, 1.31961e-02_jprb, &
       1.30051e-02_jprb, 1.28164e-02_jprb, 1.26302e-02_jprb, 1.24466e-02_jprb, 1.22659e-02_jprb, &
       1.20881e-02_jprb, 1.19131e-02_jprb, 1.17412e-02_jprb, 1.15723e-02_jprb, 1.14063e-02_jprb, &
       1.12434e-02_jprb, 1.10834e-02_jprb, 1.09264e-02_jprb, 1.07722e-02_jprb, 1.06210e-02_jprb, &
       1.04725e-02_jprb, 1.03269e-02_jprb, 1.01839e-02_jprb, 1.00436e-02_jprb, 9.90593e-03_jprb, &
       9.77080e-03_jprb, 9.63818e-03_jprb, 9.50800e-03_jprb/)
      absliq1(:,16) = (/ &
! band 16
       5.64971e-02_jprb, 9.04736e-02_jprb, 8.11726e-02_jprb, 7.05450e-02_jprb, 6.20052e-02_jprb, &
       5.54286e-02_jprb, 5.03503e-02_jprb, 4.63791e-02_jprb, 4.32290e-02_jprb, 4.06959e-02_jprb, &
       3.74690e-02_jprb, 3.52964e-02_jprb, 3.33799e-02_jprb, 3.16774e-02_jprb, 3.01550e-02_jprb, &
       2.87856e-02_jprb, 2.75474e-02_jprb, 2.64223e-02_jprb, 2.53953e-02_jprb, 2.44542e-02_jprb, &
       2.35885e-02_jprb, 2.27894e-02_jprb, 2.20494e-02_jprb, 2.13622e-02_jprb, 2.07222e-02_jprb, &
       2.01246e-02_jprb, 1.95654e-02_jprb, 1.90408e-02_jprb, 1.84398e-02_jprb, 1.80021e-02_jprb, &
       1.75816e-02_jprb, 1.71775e-02_jprb, 1.67889e-02_jprb, 1.64152e-02_jprb, 1.60554e-02_jprb, &
       1.57089e-02_jprb, 1.53751e-02_jprb, 1.50531e-02_jprb, 1.47426e-02_jprb, 1.44428e-02_jprb, &
       1.41532e-02_jprb, 1.38734e-02_jprb, 1.36028e-02_jprb, 1.33410e-02_jprb, 1.30875e-02_jprb, &
       1.28420e-02_jprb, 1.26041e-02_jprb, 1.23735e-02_jprb, 1.21497e-02_jprb, 1.19325e-02_jprb, &
       1.17216e-02_jprb, 1.15168e-02_jprb, 1.13177e-02_jprb, 1.11241e-02_jprb, 1.09358e-02_jprb, &
       1.07525e-02_jprb, 1.05741e-02_jprb, 1.04003e-02_jprb/)

      end subroutine lwcldpr

      end module rrtmg_lw_cldprop
