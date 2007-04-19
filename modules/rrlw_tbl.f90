      module rrlw_tbl

      use parkind, only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw lookup table arrays

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ntbl   :  integer: Lookup table dimension
! tblint :  real   : Lookup table conversion factor
! tautbl :  real   : Clear-sky optical depth (used in cloudy radiative
!                    transfer)
! trans  :  real   : Transmittance lookup table
! tf     :  real   : Tau transition function; i.e. the transition of
!                    the Planck function from that for the mean layer
!                    temperature to that for the layer boundary
!                    temperature as a function of optical depth.
!                    The "linear in tau" method is used to make 
!                    the table.
! bpade  :  real   : Pade constant   
!------------------------------------------------------------------

      integer(kind=jpim), parameter :: ntbl = 10000

      real(kind=jprb), parameter :: tblint = 10000.0

      real(kind=jprb) , dimension(0:ntbl) :: tautbl
      real(kind=jprb) , dimension(0:ntbl) :: trans
      real(kind=jprb) , dimension(0:ntbl) :: tf
      real(kind=jprb) :: bpade

      end module rrlw_tbl

