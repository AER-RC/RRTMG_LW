      module rrlw_kg15

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kao_mn2 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no15 = 16

      real(kind=jprb) :: fracrefao(no15,9)
      real(kind=jprb) :: kao(9,5,13,no15)
      real(kind=jprb) :: kao_mn2(9,19,no15)
      real(kind=jprb) :: selfrefo(10,no15)
      real(kind=jprb) :: forrefo(4,no15)


!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 15
! band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! ka_mn2  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng15 = 2

      real(kind=jprb) :: fracrefa(ng15,9)
      real(kind=jprb) :: ka(9,5,13,ng15) ,absa(585,ng15)
      real(kind=jprb) :: ka_mn2(9,19,ng15)
      real(kind=jprb) :: selfref(10,ng15)
      real(kind=jprb) :: forref(4,ng15)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg15
