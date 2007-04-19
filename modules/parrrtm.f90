
      module parrrtm

      use parkind ,only : jpim, jprb

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw main parameters
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! nlon   :  integer: number of columns or longitudes
! mxlay  :  integer: maximum number of layers
! mg     :  integer: number of original g-intervals per spectral band
! nbands :  integer: number of spectral bands
! maxxsec:  integer: maximum number of cross-section molecules
!                    (e.g. cfcs)
! maxinpx:  integer: 
! ngpt   :  integer: total number of reduced g-intervals for rrtmg_lw
! ngNN   :  integer: number of reduced g-intervals per spectral band
! ngsNN  :  integer: cumulative number of g-intervals per band
!------------------------------------------------------------------

! Settings for single column mode.
! For GCM use, set nlon to number of longitudes, and
! mxlay to number of model layers
      integer(kind=jpim), parameter :: nlon  = 1
      integer(kind=jpim), parameter :: mxlay  = 203
      integer(kind=jpim), parameter :: mg     = 16
      integer(kind=jpim), parameter :: nbands = 16
      integer(kind=jpim), parameter :: maxxsec= 4
      integer(kind=jpim), parameter :: mxmol  = 38
      integer(kind=jpim), parameter :: maxinpx= 38
      integer(kind=jpim), parameter :: nmol   = 7
! Use for 140 g-point model 
      integer(kind=jpim), parameter :: ngpt   = 140
! Use for 256 g-point model 
!      integer(kind=jpim), parameter :: ngpt   = 256

! Use for 140 g-point model
      integer(kind=jpim), parameter :: ng1  = 10
      integer(kind=jpim), parameter :: ng2  = 12
      integer(kind=jpim), parameter :: ng3  = 16
      integer(kind=jpim), parameter :: ng4  = 14
      integer(kind=jpim), parameter :: ng5  = 16
      integer(kind=jpim), parameter :: ng6  = 8
      integer(kind=jpim), parameter :: ng7  = 12
      integer(kind=jpim), parameter :: ng8  = 8
      integer(kind=jpim), parameter :: ng9  = 12
      integer(kind=jpim), parameter :: ng10 = 6
      integer(kind=jpim), parameter :: ng11 = 8
      integer(kind=jpim), parameter :: ng12 = 8
      integer(kind=jpim), parameter :: ng13 = 4
      integer(kind=jpim), parameter :: ng14 = 2
      integer(kind=jpim), parameter :: ng15 = 2
      integer(kind=jpim), parameter :: ng16 = 2

      integer(kind=jpim), parameter :: ngs1  = 10
      integer(kind=jpim), parameter :: ngs2  = 22
      integer(kind=jpim), parameter :: ngs3  = 38
      integer(kind=jpim), parameter :: ngs4  = 52
      integer(kind=jpim), parameter :: ngs5  = 68
      integer(kind=jpim), parameter :: ngs6  = 76
      integer(kind=jpim), parameter :: ngs7  = 88
      integer(kind=jpim), parameter :: ngs8  = 96
      integer(kind=jpim), parameter :: ngs9  = 108
      integer(kind=jpim), parameter :: ngs10 = 114
      integer(kind=jpim), parameter :: ngs11 = 122
      integer(kind=jpim), parameter :: ngs12 = 130
      integer(kind=jpim), parameter :: ngs13 = 134
      integer(kind=jpim), parameter :: ngs14 = 136
      integer(kind=jpim), parameter :: ngs15 = 138

! Use for 256 g-point model
!      integer(kind=jpim), parameter :: ng1  = 16
!      integer(kind=jpim), parameter :: ng2  = 16
!      integer(kind=jpim), parameter :: ng3  = 16
!      integer(kind=jpim), parameter :: ng4  = 16
!      integer(kind=jpim), parameter :: ng5  = 16
!      integer(kind=jpim), parameter :: ng6  = 16
!      integer(kind=jpim), parameter :: ng7  = 16
!      integer(kind=jpim), parameter :: ng8  = 16
!      integer(kind=jpim), parameter :: ng9  = 16
!      integer(kind=jpim), parameter :: ng10 = 16
!      integer(kind=jpim), parameter :: ng11 = 16
!      integer(kind=jpim), parameter :: ng12 = 16
!      integer(kind=jpim), parameter :: ng13 = 16
!      integer(kind=jpim), parameter :: ng14 = 16
!      integer(kind=jpim), parameter :: ng15 = 16
!      integer(kind=jpim), parameter :: ng16 = 16

!      integer(kind=jpim), parameter :: ngs1  = 16
!      integer(kind=jpim), parameter :: ngs2  = 32
!      integer(kind=jpim), parameter :: ngs3  = 48
!      integer(kind=jpim), parameter :: ngs4  = 64
!      integer(kind=jpim), parameter :: ngs5  = 80
!      integer(kind=jpim), parameter :: ngs6  = 96
!      integer(kind=jpim), parameter :: ngs7  = 112
!      integer(kind=jpim), parameter :: ngs8  = 128
!      integer(kind=jpim), parameter :: ngs9  = 144
!      integer(kind=jpim), parameter :: ngs10 = 160
!      integer(kind=jpim), parameter :: ngs11 = 176
!      integer(kind=jpim), parameter :: ngs12 = 192
!      integer(kind=jpim), parameter :: ngs13 = 208
!      integer(kind=jpim), parameter :: ngs14 = 224
!      integer(kind=jpim), parameter :: ngs15 = 240
!      integer(kind=jpim), parameter :: ngs16 = 256

      end module parrrtm
