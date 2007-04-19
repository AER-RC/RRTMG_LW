      module rrlw_kg04

      use parkind ,only : jpim, jprb

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
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

      integer(kind=jpim), parameter :: no4  = 16

      real(kind=jprb) :: fracrefao(no4,9)  ,fracrefbo(no4,6)
      real(kind=jprb) :: kao(9,5,13,no4)
      real(kind=jprb) :: kbo(5,5,13:59,no4)
      real(kind=jprb) :: selfrefo(10,no4)  ,forrefo(4,no4)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 4
! band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! absa    : real
! absb    : real
!fracrefa : real    
!fracrefb : real
! ka      : real     
! kb      : real     
! selfref : real     
! forref  : real     
!-----------------------------------------------------------------

      integer(kind=jpim), parameter :: ng4  = 14

      real(kind=jprb) :: fracrefa(ng4,9)  ,fracrefb(ng4,6)
      real(kind=jprb) :: ka(9,5,13,ng4)   ,absa(585,ng4)
      real(kind=jprb) :: kb(5,5,13:59,ng4),absb(1175,ng4)
      real(kind=jprb) :: selfref(10,ng4)  ,forref(4,ng4)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      end module rrlw_kg04
