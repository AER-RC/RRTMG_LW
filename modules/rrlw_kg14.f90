      module rrlw_kg14

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 14
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
!fracrefbo: real    
! kao     : real     
! kbo     : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no14 = 16

      real(kind=jprb) , dimension(no14) :: fracrefao
      real(kind=jprb) , dimension(no14) :: fracrefbo

      real(kind=jprb) :: kao(5,13,no14)
      real(kind=jprb) :: kbo(5,13:59,no14)
      real(kind=jprb) :: selfrefo(10,no14)
      real(kind=jprb) :: forrefo(4,no14)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 14
! band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
!fracrefb : real    
! ka      : real     
! kb      : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng14 = 2

      real(kind=jprb) , dimension(ng14) :: fracrefa
      real(kind=jprb) , dimension(ng14) :: fracrefb

      real(kind=jprb) :: ka(5,13,ng14)   ,absa(65,ng14)
      real(kind=jprb) :: kb(5,13:59,ng14),absb(235,ng14)
      real(kind=jprb) :: selfref(10,ng14)
      real(kind=jprb) :: forref(4,ng14)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrlw_kg14
