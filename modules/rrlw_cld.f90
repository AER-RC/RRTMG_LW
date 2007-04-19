      module rrlw_cld

      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw cloud property coefficients

! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! abscld1:  real   : 
! absice0:  real   : 
! absice1:  real   : 
! absice2:  real   : 
! absice3:  real   : 
! absliq0:  real   : 
! absliq1:  real   : 
!------------------------------------------------------------------

      real(kind=jprb) :: abscld1
      real(kind=jprb) , dimension(2) :: absice0
      real(kind=jprb) , dimension(2,5) :: absice1
      real(kind=jprb) , dimension(43,16) :: absice2
      real(kind=jprb) , dimension(46,16) :: absice3
      real(kind=jprb) :: absliq0
      real(kind=jprb) , dimension(58,16) :: absliq1

      end module rrlw_cld

