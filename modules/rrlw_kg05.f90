      module rrlw_kg05

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
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
! kao_mo3 : real     
! selfrefo: real     
! forrefo : real     
! ccl4o   : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no5  = 16

      real(kind=jprb) :: fracrefao(no5,9) ,fracrefbo(no5,5)
      real(kind=jprb) :: kao(9,5,13,no5)
      real(kind=jprb) :: kbo(5,5,13:59,no5)
      real(kind=jprb) :: kao_mo3(9,19,no5)
      real(kind=jprb) :: selfrefo(10,no5)
      real(kind=jprb) :: forrefo(4,no5)
      real(kind=jprb) :: ccl4o(no5)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
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
! ka_mo3  : real     
! selfref : real     
! forref  : real     
! ccl4    : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng5  = 16

      real(kind=jprb) :: fracrefa(ng5,9) ,fracrefb(ng5,5)
      real(kind=jprb) :: ka(9,5,13,ng5)   ,absa(585,ng5)
      real(kind=jprb) :: kb(5,5,13:59,ng5),absb(1175,ng5)
      real(kind=jprb) :: ka_mo3(9,19,ng5)
      real(kind=jprb) :: selfref(10,ng5)
      real(kind=jprb) :: forref(4,ng5)
      real(kind=jprb) :: ccl4(ng5)
      
      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      end module rrlw_kg05

