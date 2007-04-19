      module rrlw_kg07

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no7  = 16

      real(kind=jprb) , dimension(no7) :: fracrefbo
      real(kind=jprb) :: fracrefao(no7,9)
      real(kind=jprb) :: kao(9,5,13,no7)
      real(kind=jprb) :: kbo(5,13:59,no7)
      real(kind=jprb) :: kao_mco2(9,19,no7)
      real(kind=jprb) :: kbo_mco2(19,no7)
      real(kind=jprb) :: selfrefo(10,no7)
      real(kind=jprb) :: forrefo(4,no7)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 7
! band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
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
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng7  = 12

      real(kind=jprb) , dimension(ng7) :: fracrefb
      real(kind=jprb) :: fracrefa(ng7,9)
      real(kind=jprb) :: ka(9,5,13,ng7) ,absa(585,ng7)
      real(kind=jprb) :: kb(5,13:59,ng7),absb(235,ng7)
      real(kind=jprb) :: ka_mco2(9,19,ng7)
      real(kind=jprb) :: kb_mco2(19,ng7)
      real(kind=jprb) :: selfref(10,ng7)
      real(kind=jprb) :: forref(4,ng7)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg07
