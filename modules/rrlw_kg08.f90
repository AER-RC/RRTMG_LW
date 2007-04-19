      module rrlw_kg08

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! kao_mco2: real     
! kbo_mco2: real     
! kao_mn2o: real     
! kbo_mn2o: real     
! kao_mo3 : real     
! selfrefo: real     
! forrefo : real     
! cfc12o  : real     
!cfc22adjo: real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no8  = 16

      real(kind=jprb) , dimension(no8) :: fracrefao
      real(kind=jprb) , dimension(no8) :: fracrefbo
      real(kind=jprb) , dimension(no8) :: cfc12o
      real(kind=jprb) , dimension(no8) :: cfc22adjo

      real(kind=jprb) :: kao(5,13,no8)
      real(kind=jprb) :: kao_mco2(19,no8)
      real(kind=jprb) :: kao_mn2o(19,no8)
      real(kind=jprb) :: kao_mo3(19,no8)
      real(kind=jprb) :: kbo(5,13:59,no8)
      real(kind=jprb) :: kbo_mco2(19,no8)
      real(kind=jprb) :: kbo_mn2o(19,no8)
      real(kind=jprb) :: selfrefo(10,no8)
      real(kind=jprb) :: forrefo(4,no8)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! ka_mco2 : real     
! kb_mco2 : real     
! ka_mn2o : real     
! kb_mn2o : real     
! ka_mo3  : real     
! selfref : real     
! forref  : real     
! cfc12   : real     
! cfc22adj: real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng8  = 8

      real(kind=jprb) , dimension(ng8) :: fracrefa
      real(kind=jprb) , dimension(ng8) :: fracrefb
      real(kind=jprb) , dimension(ng8) :: cfc12
      real(kind=jprb) , dimension(ng8) :: cfc22adj

      real(kind=jprb) :: ka(5,13,ng8)    ,absa(65,ng8)
      real(kind=jprb) :: kb(5,13:59,ng8) ,absb(235,ng8)
      real(kind=jprb) :: ka_mco2(19,ng8)
      real(kind=jprb) :: ka_mn2o(19,ng8)
      real(kind=jprb) :: ka_mo3(19,ng8)
      real(kind=jprb) :: kb_mco2(19,ng8)
      real(kind=jprb) :: kb_mn2o(19,ng8)
      real(kind=jprb) :: selfref(10,ng8)
      real(kind=jprb) :: forref(4,ng8)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg08

