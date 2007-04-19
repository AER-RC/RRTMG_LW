      module rrlw_kg13

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kao_mco2: real     
! kao_mco : real     
! kbo_mo3 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: no13 = 16

      real(kind=jprb) , dimension(no13) :: fracrefbo

      real(kind=jprb) :: fracrefao(no13,9)
      real(kind=jprb) :: kao(9,5,13,no13)
      real(kind=jprb) :: kao_mco2(9,19,no13)
      real(kind=jprb) :: kao_mco(9,19,no13)
      real(kind=jprb) :: kbo_mo3(19,no13)
      real(kind=jprb) :: selfrefo(10,no13)
      real(kind=jprb) :: forrefo(4,no13)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! ka_mco2 : real     
! ka_mco  : real     
! kb_mo3  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng13 = 4

      real(kind=jprb) , dimension(ng13) :: fracrefb

      real(kind=jprb) :: fracrefa(ng13,9)
      real(kind=jprb) :: ka(9,5,13,ng13) ,absa(585,ng13)
      real(kind=jprb) :: ka_mco2(9,19,ng13)
      real(kind=jprb) :: ka_mco(9,19,ng13)
      real(kind=jprb) :: kb_mo3(19,ng13)
      real(kind=jprb) :: selfref(10,ng13)
      real(kind=jprb) :: forref(4,ng13)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg13
