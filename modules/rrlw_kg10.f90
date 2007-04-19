      module rrlw_kg10

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
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

      integer(kind=jpim), parameter :: no10 = 16

      real(kind=jprb) , dimension(no10) :: fracrefao
      real(kind=jprb) , dimension(no10) :: fracrefbo

      real(kind=jprb) :: kao(5,13,no10)
      real(kind=jprb) :: kbo(5,13:59,no10)
      real(kind=jprb) :: selfrefo(10,no10)
      real(kind=jprb) :: forrefo(4,no10)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
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
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng10 = 6

      real(kind=jprb) , dimension(ng10) :: fracrefa
      real(kind=jprb) , dimension(ng10) :: fracrefb

      real(kind=jprb) :: ka(5,13,ng10)   , absa(65,ng10)
      real(kind=jprb) :: kb(5,13:59,ng10), absb(235,ng10)
      real(kind=jprb) :: selfref(10,ng10)
      real(kind=jprb) :: forref(4,ng10)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg10
