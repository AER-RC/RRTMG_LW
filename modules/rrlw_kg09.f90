      module rrlw_kg09

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
! kao_mn2o: real     
! kbo_mn2o: real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no9  = 16

      real(kind=jprb) , dimension(no9) :: fracrefbo

      real(kind=jprb) :: fracrefao(no9,9)
      real(kind=jprb) :: kao(9,5,13,no9)
      real(kind=jprb) :: kbo(5,13:59,no9)
      real(kind=jprb) :: kao_mn2o(9,19,no9)
      real(kind=jprb) :: kbo_mn2o(19,no9)
      real(kind=jprb) :: selfrefo(10,no9)
      real(kind=jprb) :: forrefo(4,no9)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
! ka_mn2o : real     
! kb_mn2o : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng9  = 12

      real(kind=jprb) , dimension(ng9) :: fracrefb
      real(kind=jprb) :: fracrefa(ng9,9)
      real(kind=jprb) :: ka(9,5,13,ng9) ,absa(585,ng9)
      real(kind=jprb) :: kb(5,13:59,ng9) ,absb(235,ng9)
      real(kind=jprb) :: ka_mn2o(9,19,ng9)
      real(kind=jprb) :: kb_mn2o(19,ng9)
      real(kind=jprb) :: selfref(10,ng9)
      real(kind=jprb) :: forref(4,ng9)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg09
