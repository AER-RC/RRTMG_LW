      module rrlw_ref

      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw reference atmosphere 
! Based on standard mid-latitude summer profile
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! pref   :  real   : Reference pressure levels
! preflog:  real   : Reference pressure levels, ln(pref)
! tref   :  real   : Reference temperature levels for MLS profile
! chi_mls:  real   : 
!------------------------------------------------------------------

      real(kind=jprb) , dimension(59) :: pref
      real(kind=jprb) , dimension(59) :: preflog
      real(kind=jprb) , dimension(59) :: tref
      real(kind=jprb) :: chi_mls(7,59)

      end module rrlw_ref
