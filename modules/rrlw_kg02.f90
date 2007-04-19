      module rrlw_kg02

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
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

      integer(kind=jpim), parameter :: no2  = 16

      real(kind=jprb) :: fracrefao(no2)   , fracrefbo(no2)
      real(kind=jprb) :: kao(5,13,no2)
      real(kind=jprb) :: kbo(5,13:59,no2)
      real(kind=jprb) :: selfrefo(10,no2) , forrefo(4,no2)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 2
! band 2:  250-500 cm-1 (low - h2o; high - h2o)
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
! absa    : real
! absb    : real
! selfref : real     
! forref  : real
!
! refparam: real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng2  = 12

      real(kind=jprb) :: fracrefa(ng2)  , fracrefb(ng2)
      real(kind=jprb) :: ka(5,13,ng2)   , absa(65,ng2)
      real(kind=jprb) :: kb(5,13:59,ng2), absb(235,ng2)
      real(kind=jprb) :: selfref(10,ng2), forref(4,ng2)

      real(kind=jprb) :: refparam(13)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg02


