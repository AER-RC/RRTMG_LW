C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

* Copyright 2002, 2003, 2004, Atmospheric & Environmental Research, Inc. (AER).
* This software may be used, copied, or redistributed as long as it is
* not sold and this copyright notice is reproduced on each copy made.
* This model is provided as is without any express or implied warranties.
*                      (http://www.rtweb.aer.com/)
*
****************************************************************************
*                                                                          *
*                               RRTM                                       *
*                                                                          *
*                                                                          *
*                                                                          *
*                   A RAPID RADIATIVE TRANSFER MODEL                       *
*                       FOR THE LONGWAVE REGION                            * 
*                                                                          *
*                                                                          *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                  *
*                        131 HARTWELL AVENUE                               *
*                        LEXINGTON, MA 02421                               *
*                                                                          *
*                                                                          *
*                         ELI J. MLAWER                                    *
*                         JENNIFER S. DELAMERE                             *
*                         STEVEN J. TAUBMAN~                               *
*                         SHEPARD A. CLOUGH                                *
*                                                                          *
*                                                                          *
*                         ~currently at GFDL                               *
*                                                                          *
*                                                                          *
*                                                                          *
*                       email:  mlawer@aer.com                             *
*                       email:  jdelamer@aer.com                           *
*                                                                          *
*        The authors wish to acknowledge the contributions of the          *
*        following people:  Karen Cady-Pereira, Patrick D. Brown,          *
*        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom. *
*                                                                          *
****************************************************************************

C *** This version of RRTM has been modified to use a reduced set of 
C *** g-points for application of the model to GCMs.  
C *** Michael J. Iacono, AER, Inc.

      SUBROUTINE RRTMG_LW(pplay,paprs,t,tlev,tsol,q,wo,co2mmr,ch4mmr,
     1             n2ommr,cfc11mmr,cfc12mmr,emis,icovlp,cldfra,
     2             inflglw,iceflglw,liqflglw,
     3             cldd1lw,cldd2lw,cldd3lw,cldd4lw,
     4             fulaer,fdlaer,qrlaer,fulcaer,fdlcaer,qrlcaer)

C *** This program is the driver for RRTM, the AER LW radiation model.  
C     This routine:
C     a) calls INATM to read in the atmospheric profile from GCM
C     b) calls SETCOEF to calculate various quantities needed for 
C        the radiative transfer algorithm
C     c) calls TAUGBn to calculate gaseous optical depths for each 
C        of the 16 spectral bands
C     d) calls CLDPROP to set cloud optical depth based on input
C        cloud optical properties
C     e) calls RTRNMR (for both clear and cloudy profiles) to do the
C        radiative transfer calculation with a maximum-random cloud
C        overlap method, or calls RTRN to use random cloud overlap
C     f) passes the necessary fluxes and cooling rates back to GCM

C Parameters
      PARAMETER (MXLAY=203)
      PARAMETER (NGPT=140)
      PARAMETER (NBANDS=16)

C Input
      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/   ONEMINUS
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /IFIL/      IRD,IPR,IPU,IDUM(15)
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)

C Output
      COMMON /CLDFLG/    ICLDLYR(MXLAY)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /OUTCLR/    TOTUCLFL(0:MXLAY), TOTDCLFL(0:MXLAY),
     &                   FNETC(0:MXLAY), HTRC(0:MXLAY)

C RRTM Definitions
C    NGPT                         ! Total number of g-point subintervals
C    MXLAY                        ! Maximum number of model layers
C    NBANDS                       ! Number of longwave spectral bands
C    PI                           ! Geometric constant
C    FLUXFAC                      ! Radiance to flux conversion factor 
C    HEATFAC                      ! Heating rate conversion factor 
C    NG(NBANDS)                   ! Number of g-points per band for input
C                                   absorption coefficient data
C    NSPA(NBANDS),NSPB(NBANDS)    ! Number of reference atmospheres per band
C    WAVENUM1(NBANDS)             ! Longwave band lower limit (wavenumbers)
C    WAVENUM2(NBANDS)             ! Longwave band upper limit (wavenumbers)
C    DELWAVE                      ! Longwave band width (wavenumbers)
C    NUMANGS, IOUT                ! Inactive in this version
C    ISTART,IEND                  ! Beginning and ending bands to include
C    NLAYERS                      ! Number of model layers (klev+1)
C    PAVEL(MXLAY)                 ! Layer pressures (mb)
C    PZ(0:MXLAY)                  ! Level (interface) pressures (mb)
C    TAVEL(MXLAY)                 ! Layer temperatures (K)
C    TZ(0:MXLAY)                  ! Level (interface) temperatures(mb)
C    TBOUND                       ! Surface temperature (K)
C    IREFLECT                     ! Inactive
C    SEMISS(NBANDS)               ! LW surface emissivity
C    CLDFRAC(MXLAY)               ! Layer cloud fraction
C    TAUCLOUD(MXLAY)              ! Layer cloud optical depth
C    ITR(NGPT,MXLAY)              ! Integer look-up table index
C    ICLDLYR(MXLAY)               ! Flag for cloudy layers
C    TOTUFLUX(0:MXLAY)            ! Upward longwave flux (W/m2)
C    TOTDFLUX(0:MXLAY)            ! Downward longwave flux (W/m2)
C    FNET(0:MXLAY)                ! Net longwave flux (W/m2)
C    HTR(0:MXLAY)                 ! Longwave heating rate (K/day)
C    TOTUCLFL(0:MXLAY)            ! Clear sky upward longwave flux (W/m2)
C    TOTDCLFL(0:MXLAY)            ! Clear sky downward longwave flux (W/m2)
C    FNETC(0:MXLAY)               ! Clear sky net longwave flux (W/m2)
C    HTRC(0:MXLAY)                ! Clear sky longwave heating rate (K/day)
C

C Imported commons and dimensions from GCM
C
C------------------------------Modules----------------------------------
C
C------------------------------Parameters-------------------------------
C 
C Add GCM include statements here to pass in the following fields:
C   KLON    Number of longitudes
C   KLEV    Number of model layers
C
C------------------------------Commons----------------------------------
C
C------------------------------Dimensions-------------------------------
C
C Input arrays from GCM
C
      real pplay(klon,klev)           ! Layer pressures (Pa)
      real paprs(klon,klev+1)         ! Interface pressures (Pa)
      real t(klon,klev)               ! Layer temperatures (K)
      real tlev(klon,klev+1)          ! Interface temperatures (K)
      real tsol(klon)                 ! Surface temperature (K)
      real q(klon,klev)               ! H2O specific humidity (mmr)
      real wo(klon,klev)              ! O3 mass mixing ratio
      real co2mmr, ch4mmr, n2ommr     ! CO2, CH4, and N2O mmr
      real cfc11mmr, cfc12mmr         ! CFC11, CFC12 mass mixing ratio
      real emis(klon)                 ! Surface emissivity
      real cldfra(klon,klev)          ! Cloud fraction
      real cldd1lw(klon,klev)         ! Cloud optical depth or water path
      real cldd2lw(klon,klev)         ! Cloud ice fraction
      real cldd3lw(klon,klev)         ! Cloud ice effective radius (microns)
C                                     !   or ice generalized effective size
      real cldd4lw(klon,klev)         ! Cloud water drop effective radius
C See cldprop for further information:
      integer inflglw                 ! Flag for cloud optical properties
      integer iceflglw                ! Flag for ice particle specification
      integer liqflglw                ! Flag for liquid droplet specification

C Output arrays to GCM
      real fulaer(klon,klev+1)        ! Total sky upward flux
      real fdlaer(klon,klev+1)        ! Total sky downward flux
      real qrlaer(klon,klev)          ! Radiative heating rate
      real fulcaer(klon,klev+1)       ! Clear sky upward flux
      real fdlcaer(klon,klev+1)       ! Clear sky downward flux
      real qrlcaer(klon,klev)         ! Clear sky Radiative heating rate


      ONEMINUS = 1. - 1.E-6
      PI = 2.*ASIN(1.)
      FLUXFAC = PI * 2.D4  
      ISTART = 1
      IEND = 16

C  This is the main longitude/column loop within RRTM.
      NUMATMOS = 1
      do 4000 iplon = 1, NUMATMOS

C  Prepare atmospheric profile from GCM for use in RRTM, and define
C  other RRTM input parameters.  Arrays are passed back through the
C  existing RRTM commons and arrays.  Cloud fraction and cloud optical
C  properties are passed to RRTM_LW arrays here.

         call inatm(iplon,paprs,pplay,t,tlev,tsol,q,wo,co2mmr,
     1          ch4mmr,n2ommr,cfc11mmr,cfc12mmr,cldfra,emis,
     2          inflglw,iceflglw,liqflglw,
     3          cldd1lw,cldd2lw,cldd3lw,cldd4lw)

C  For cloudy atmosphere, use CLDPROP to set cloud optical depth based on
C  input cloud optical properties.  Select method based on choices described
C  in CLDPROP.  Cloud fraction, water path, liquid droplet and ice particle
C  effective radius must be passed into CLDFRAC and CLDDAT arrays.  Cloud 
C  fraction and cloud optical depth are transferred to RRTM_LW arrays in 
C  CLDPROP.  By default, cloud fraction and cloud optical depth are passed 
C  directly from GCM to RRTM_LW arrays in INATM.

         IF (ICOVLP.GE.1) THEN
            CALL CLDPROP(ICLDATM)
         ENDIF

C  Calculate information needed by the radiative transfer routine
C  that is specific to this atmosphere, especially some of the 
C  coefficients and indices needed to compute the optical depths
C  by interpolating data from stored reference atmospheres. 

         CALL SETCOEF

C ***    Calculate the gaseous optical depths and Planck fractions for 
C        each longwave spectral band.
         CALL TAUGB1
         CALL TAUGB2
         CALL TAUGB3
         CALL TAUGB4
         CALL TAUGB5
         CALL TAUGB6
         CALL TAUGB7
         CALL TAUGB8
         CALL TAUGB9
         CALL TAUGB10
         CALL TAUGB11
         CALL TAUGB12
         CALL TAUGB13
         CALL TAUGB14
         CALL TAUGB15
         CALL TAUGB16

C  Call the radiative transfer routine.
C  Either routine can be called to do clear sky calculation.  If clouds
C  are present, then select routine based on cloud overlap assumption
C  to be used.  Clear sky calculation is done simultaneously.

        IF (ICOVLP .EQ. 1) THEN
           CALL RTRN 
        ELSE
           CALL RTRNMR
        ENDIF

C  Pass total sky up and down flux, and heating rate profiles to 
C  output arrays.

           do 2400 k = 0, NLAYERS
              fulaer(iplon,k+1) = TOTUFLUX(k)
              fdlaer(iplon,k+1) = TOTDFLUX(k)
              fulcaer(iplon,k+1) = TOTUCLFL(k)
              fdlcaer(iplon,k+1) = TOTDCLFL(k)
 2400      continue
           do 2450 k = 0, NLAYERS-1
              qrlaer(iplon,k+1) = HTR(k)
              qrlcaer(iplon,k+1) = HTRC(k)
 2450      continue

 4000 CONTINUE

      CLOSE(IWR)

      RETURN
      END

C***************************************************************************
      BLOCK DATA RRTMDAT
C***************************************************************************

C Parameters
      PARAMETER (NBANDS=16)
      PARAMETER (MAXINPX=35)

C Commons
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /XSECCTRL/  NXMOL,IXINDX(MAXINPX)

C  Longwave spectral band data
      DATA WAVENUM1(1) /10./, WAVENUM2(1) /350./, DELWAVE(1) /340./
      DATA WAVENUM1(2) /350./, WAVENUM2(2) /500./, DELWAVE(2) /150./
      DATA WAVENUM1(3) /500./, WAVENUM2(3) /630./, DELWAVE(3) /130./
      DATA WAVENUM1(4) /630./, WAVENUM2(4) /700./, DELWAVE(4) /70./
      DATA WAVENUM1(5) /700./, WAVENUM2(5) /820./, DELWAVE(5) /120./
      DATA WAVENUM1(6) /820./, WAVENUM2(6) /980./, DELWAVE(6) /160./
      DATA WAVENUM1(7) /980./, WAVENUM2(7) /1080./, DELWAVE(7) /100./
      DATA WAVENUM1(8) /1080./, WAVENUM2(8) /1180./, DELWAVE(8) /100./
      DATA WAVENUM1(9) /1180./, WAVENUM2(9) /1390./, DELWAVE(9) /210./
      DATA WAVENUM1(10) /1390./,WAVENUM2(10) /1480./,DELWAVE(10) /90./
      DATA WAVENUM1(11) /1480./,WAVENUM2(11) /1800./,DELWAVE(11) /320./
      DATA WAVENUM1(12) /1800./,WAVENUM2(12) /2080./,DELWAVE(12) /280./
      DATA WAVENUM1(13) /2080./,WAVENUM2(13) /2250./,DELWAVE(13) /170./
      DATA WAVENUM1(14) /2250./,WAVENUM2(14) /2380./,DELWAVE(14) /130./
      DATA WAVENUM1(15) /2380./,WAVENUM2(15) /2600./,DELWAVE(15) /220./
      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3250./,DELWAVE(16) /650./

      DATA NG /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/
      DATA NSPB /1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(86400)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/

C     NXMOL     - number of cross-sections input by user
C     IXINDX(I) - index of cross-section molecule corresponding to Ith
C                 cross-section specified by user
C                 = 0 -- not allowed in RRTM
C                 = 1 -- CCL4
C                 = 2 -- CFC11
C                 = 3 -- CFC12
C                 = 4 -- CFC22
      DATA NXMOL  /2/
      DATA IXINDX /0,2,3,0,31*0/
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     units are generally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c          RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07
c          RADCN2 = PLANCK*CLIGHT/BOLTZ

      END

C***************************************************************************
      SUBROUTINE INATM(iplon,paprs,pplay,t,tlev,tsol,q,wo,co2mmr,
     1             ch4mmr,n2ommr,cfc11mmr,cfc12mmr,cldfra,emis,
     2             inflglw,iceflglw,liqflglw,
     3             cldd1lw,cldd2lw,cldd3lw,cldd4lw)
C***************************************************************************
C  RRTM Longwave Radiative Transfer Model
C  Atmospheric and Environmental Research, Inc., Cambridge, MA
C
C  Revision for use with GCMs:  Michael J. Iacono; August, 2003
C
C  Input atmospheric profile from host model, and prepare it for use in RRTM.
C  Set other RRTM input parameters.  Values are passed back through existing
C  RRTM arrays and commons.
C***************************************************************************

C Parameters
      PARAMETER (NBANDS=16)
      PARAMETER (MXLAY=203)
      PARAMETER (MAXINPX=35)
      PARAMETER (MXCBANDS=5)
      PARAMETER (MAXXSEC=4)
      PARAMETER (MAXPRDW = MXLAY*35)
      PARAMETER (MAXPROD = MXLAY*MAXXSEC)
C Input                              
      COMMON /CONTROL/  NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
C Output
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY)
      COMMON /CLOUDDAT/ NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)
      COMMON /SURFACE/  TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /PWV/      PWVCM
C RRTM Definitions
C    MXLAY                        ! Maximum number of model layers
C    MAXXSEC                      ! Maximum number of cross sections
C    NUMANGS, IOUT                ! Inactive in this version
C    ISTART,IEND                  ! Beginning and ending bands to include
C    NLAYERS                      ! Number of model layers (plev+1)
C    PAVEL(MXLAY)                 ! Layer pressures (mb)
C    PZ(0:MXLAY)                  ! Level (interface) pressures (mb)
C    TAVEL(MXLAY)                 ! Layer temperatures (K)
C    TZ(0:MXLAY)                  ! Level (interface) temperatures(mb)
C    TBOUND                       ! Surface temperature (K)
C    COLDRY(MXLAY)                ! Dry air column density (molecules/cm2)
C    WBRODL(MXLAY)                ! Broadening gas column density (molec/cm2)
C    WKL(35,MXLAY)                ! Molecular amounts (molecules/cm2)
C    WX(MAXXSEC)                  ! Cross-section amounts (molecules/cm2)
C    CLDFRAC(MXLAY)               ! Layer cloud fraction
C    TAUCLOUD(MXLAY)              ! Layer cloud optical depth
C    PWVCM                        ! Precipitable water vapor (cm)
C    SBC                          ! Stefan-Boltzmann constant
C    AMD                          ! Atomic weight of dry air
C    AMW                          ! Atomic weight of water
C    AMO                          ! Atomic weight of ozone
C    AMCH4                        ! Atomic weight of methane
C    AMN2O                        ! Atomic weight of nitrous oxide
C    AMC11                        ! Atomic weight of CFC-11
C    AMC12                        ! Atomic weight of CFC-12
C    NXMOL                        ! Number of cross-section molecules
C    IXINDX                       ! Cross-section molecule index (see below)
C    IXSECT                       ! On/off flag for cross-sections (inactive)
C    IXMAX                        ! Maximum number of cross-sections (inactive)
C

C Imported commons and dimensions from GCM
C
C------------------------------Modules----------------------------------
C
C------------------------------Parameters-------------------------------
C 
C Add GCM include statements here to pass in the following fields:
C   KLON    Number of longitudes
C   KLEV    Number of model layers
C
C------------------------------Commons----------------------------------
C
C------------------------------Dimensions-------------------------------
C
C Input arrays from GCM
C
      real pplay(klon,klev)           ! Layer pressures (Pa)
      real paprs(klon,klev+1)         ! Interface pressures (Pa)
      real t(klon,klev)               ! Layer temperatures (K)
      real tlev(klon,klev+1)          ! Interface temperatures (K)
      real tsol(klon)                 ! Surface temperature (K)
      real q(klon,klev)               ! H2O specific humidity (mmr)
      real wo(klon,klev)              ! O3 mass mixing ratio
      real co2mmr, ch4mmr, n2ommr     ! CO2, CH4, and N2O mmr
      real cfc11mmr, cfc12mmr         ! CFC11, CFC12 mass mixing ratio
      real cldfra(klon,klev)          ! Cloud fraction
      real emis(klon)                 ! Surface emissivity
      real cldd1lw(klon,klev)         ! Cloud optical depth or water path
      real cldd2lw(klon,klev)         ! Cloud ice fraction
      real cldd3lw(klon,klev)         ! Cloud ice effective radius (microns)
C                                     !   or ice generalized effective size
      real cldd4lw(klon,klev)         ! Cloud water drop effective radius
      integer inflglw                 ! Flag for cloud optical properties
      integer iceflglw                ! Flag for ice particle specification
      integer liqflglw                ! Flag for liquid droplet specification

C Internal arrays
      real amd                  ! Effective molecular weight of dry air (g/mol)
      real amw                  ! Molecular weight of water vapor (g/mol)
      real amc                  ! Molecular weight of carbon dioxide (g/mol)
      real amo                  ! Molecular weight of ozone (g/mol)
      real amo2                 ! Molecular weight of oxygen (g/mol)
      real amch4                ! Molecular weight of methane (g/mol)
      real amn2o                ! Molecular weight of nitrous oxide (g/mol)
      real amc11                ! Molecular weight of CFC11 (g/mol) - CFCL3
      real amc12                ! Molecular weight of CFC12 (g/mol) - CF2CL2
c      real avgdro               ! Avogadro's number (molecules/mole)
      real sbc                  ! Stefan-Boltzmann constant (W/m2K4)
      real o2mmr                ! O2 mass mixing ratio

C Stefan-Boltzmann constant
      parameter (sbc = 5.67e-8)
C Acceleration of gravity (m/s2)
      parameter (grav = 9.8066)
C Oxygen mass mixing ratio (as defined in ccm3.6)
      data o2mmr /  0.23143   /

C Atomic weights for conversion from mass to volume mixing ratios
      data amd   /  28.9644   /
      data amw   /  18.0154   /
      data amc   /  44.0098   /
      data amo   /  47.9998   /
      data amo2  /  31.9999   /
      data amch4 /  16.0430   /
      data amn2o /  44.0128   /
      data amc11 / 137.3684   / 
      data amc12 / 120.9138   /

C Set molecular weight ratios
      real amdw,                ! Molecular weight of dry air / water vapor
     $     amdc,                ! Molecular weight of dry air / carbon dioxide
     $     amdo,                ! Molecular weight of dry air / ozone
     $     amdm,                ! Molecular weight of dry air / methane
     $     amdn,                ! Molecular weight of dry air / nitrous oxide
     $     amdo2,               ! Molecular weight of dry air / oxygen
     $     amdc1,               ! Molecular weight of dry air / CFC11
     $     amdc2                ! Molecular weight of dry air / CFC12
      data amdw /  1.607758 /
      data amdc /  0.658114 /
      data amdo /  0.603428 /
      data amdm /  1.805423 /
      data amdn /  0.658090 /
      data amdo2/  0.905140 /
      data amdc1/  0.210852 /
      data amdc2/  0.239546 /

C  Activate cross section molecules: (see BLOCK DATA RRTMDAT)
C     NXMOL     - number of cross-sections input by user
C     IXINDX(I) - index of cross-section molecule corresponding to Ith
C                 cross-section specified by user
C                 = 0 -- not allowed in RRTM
C                 = 1 -- CCL4
C                 = 2 -- CFC11
C                 = 3 -- CFC12
C                 = 4 -- CFC22

C  Initialize all molecular amounts to zero here, then pass input amounts
C  into RRTM arrays WKL and WX below.

      DO 1000 ILAY = 1,MXLAY
         DO 1100 ISP = 1,35
 1100       WKL(ISP,ILAY) = 0.0
         DO 1200 ISP = 1,MAXXSEC
 1200       WX(ISP,ILAY) = 0.0
 1000 CONTINUE
      AMTTL = 0
      WVTTL = 0
 
C  Set parameters needed for RRTM execution:
      NUMANGS = 0
      IOUT = -1
      IXSECT = 1
      IXMAX = 4

C  Set surface temperature.
      TBOUND = tsol(iplon)

C  Install input arrays into RRTM arrays for pressure, temperature,
C  and molecular amounts.  Pressures are converted from Pascals
C  (GCM) to mb (RRTM).  Molecular amounts are converted from mass
C  mixing ratio to volume mixing ratio.  CO2 and trace gas mmr are
C  constant at all levels.  The dry air column COLDRY (in molec/cm2) 
C  is calculated from the level pressures PZ (in mb) based on the 
C  hydrostatic equation and includes a correction to account for H2O 
C  in the layer.  The molecular weight of moist air (amm) is calculated 
C  for each layer.  
C  Note: By default, both the incoming pressure levels and RRTM levels 
C  count from bottom to top.

c      write(23,*) 'INATM'
c      write(23,*) 'iplon,klon,klev:',iplon, klon, klev
c      write(23,*) 'paprs:'
c      write(23,9990) (paprs(iplon,k),k=1,klev+1)
c      write(23,*) 'pplay:'
c      write(23,9990) (pplay(iplon,k),k=1,klev)
c      write(23,*) 'tlev:'
c      write(23,9990) (tlev(iplon,k),k=1,klev+1)
c      write(23,*) 't:'
c      write(23,9990) (t(iplon,k),k=1,klev)
c      write(23,*) 'q:'
c      write(23,9990) (q(iplon,k),k=1,klev)
c      write(23,*) 'wo:'
c      write(23,9990) (wo(iplon,k),k=1,klev)
c      write(23,*) 'co2mmr, ch4mmr, n2ommr:', co2mmr, ch4mmr, n2ommr
c      write(23,*) 'cfc11mmr, cfc12mmr:', cfc11mmr, cfc12mmr
c      write(23,*) 'grav,avogad: ', grav, avogad
c 9990 format(1p6e12.5)

      NLAYERS = klev
      NMOL = 7
      PZ(0) = paprs(iplon,1)*1.E-2
      TZ(0) = tlev(iplon,1)
      DO 2000 L = 1, NLAYERS
         PAVEL(L) = pplay(iplon,L)*1.E-2
         TAVEL(L) = t(iplon,L)
         PZ(L) = paprs(iplon,L+1)*1.E-2
         TZ(L) = tlev(iplon,L+1)
         WKL(1,L) = q(iplon,L)*amdw
         WKL(2,L) = co2mmr*amdc
         WKL(3,L) = wo(iplon,L)*amdo
         WKL(4,L) = n2ommr*amdn
         WKL(6,L) = ch4mmr*amdm
         WKL(7,L) = o2mmr*amdo2
         amm = (1-WKL(1,L))*amd + WKL(1,L)*amw            
         COLDRY(L) = (PZ(L-1)-PZ(L))*1.E3*avogad/
     1                      (1.E2*grav*amm*(1+WKL(1,L)))
 2000    CONTINUE

C  Set cross section molecule amounts from input; convert to vmr
      DO 2100 L=1, NLAYERS
         WX(2,L) = cfc11mmr*amdc1
         WX(3,L) = cfc12mmr*amdc2
 2100 CONTINUE

c      write(23,*) ' '
c      write(23,*) 'pz:'
c      write(23,9990) (pz(k),k=0,klev)
c      write(23,*) 'pavel:'
c      write(23,9990) (pavel(k),k=1,klev)
c      write(23,*) 'tz:'
c      write(23,9990) (tz(k),k=0,klev)
c      write(23,*) 'tavel:'
c      write(23,9990) (tavel(k),k=1,klev)
c      write(23,*) 'coldry:'
c      write(23,9990) (coldry(k),k=1,klev)
c      write(23,*) 'wbrodl:'
c      write(23,9990) (wbrodl(k),k=1,klev)
c      write(23,*) 'wkl(1), h2o:'
c      write(23,9990) (wkl(1,k),k=1,klev)
c      write(23,*) 'wkl(2), co2:'
c      write(23,9990) (wkl(2,k),k=1,klev)
c      write(23,*) 'wkl(3), ozone:'
c      write(23,9990) (wkl(3,k),k=1,klev)
c      write(23,*) 'wkl(4), n2o:'
c      write(23,9990) (wkl(4,k),k=1,klev)
c      write(23,*) 'wkl(6), ch4:'
c      write(23,9990) (wkl(6,k),k=1,klev)
c      write(23,*) 'wx(2), cfc11:'
c      write(23,9990) (wx(2,k),k=1,klev)
c      write(23,*) 'wx(3), cfc12:'
c      write(23,9990) (wx(3,k),k=1,klev)
      
C  Following section is commented.  It can be used to set profile for
C  additional layer from model top to 0 mb.
C  Set up values for extra layer at top of the atmosphere.
C  The top layer temperature for all gridpoints is set to the top layer-1
C  temperature plus a constant (17.7 K) that represents a global average
C  temperature increase above 3 mb.  Top layer interface temperatures are
C  linearly interpolated from the layer temperatures.
C  Note: The top layer temperature and ozone amount are based on a 0-3mb
C  top layer and must be modified if the layering is changed.   

C      PAVEL(NLAYERS) = 0.5*PZ(NLAYERS-1)
C      TAVEL(NLAYERS) = TAVEL(NLAYERS-1)+17.7
C      PZ(NLAYERS) = 0.00
C      TZ(NLAYERS-1) = 0.5*(TAVEL(NLAYERS)+TAVEL(NLAYERS-1))
C      TZ(NLAYERS) = TZ(NLAYERS-1)+17.7
C      WKL(1,NLAYERS) = WKL(1,NLAYERS-1)
C      WKL(2,NLAYERS) = rco2*amdc
C      WKL(3,NLAYERS) = 0.6*WKL(3,NLAYERS-1)
C      WKL(4,NLAYERS) = WKL(4,NLAYERS-1)
C      WKL(6,NLAYERS) = WKL(6,NLAYERS-1)
C      WKL(7,NLAYERS) = WKL(7,NLAYERS-1)
C      amm = (1-WKL(1,NLAYERS-1))*amd + WKL(1,NLAYERS-1)*amw            
C      COLDRY(NLAYERS) = (PZ(NLAYERS-1))*1.E3*avgdro/
C     1                      (1.E2*grav*amm*(1+WKL(1,NLAYERS-1)))
C      WX(2,NLAYERS) = WX(2,NLAYERS-1)
C      WX(3,NLAYERS) = WX(3,NLAYERS-1)

C  Here, all molecules in WKL and WX are in volume mixing ratio; convert to
C  molec/cm2 based on COLDRY for use in RRTM.  Also, compute precipitable
C  water vapor for diffusivity angle adjustments in RTRN and RTRNMR.

      DO 5000 L = 1, NLAYERS
         SUMMOL = 0.0
         DO 4100 IMOL = 2, NMOL
            SUMMOL = SUMMOL + WKL(IMOL,L)
 4100    CONTINUE
         WBRODL(L) = COLDRY(L) * (1. - SUMMOL)
         DO 4200 IMOL = 1, NMOL
            WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
 4200    CONTINUE
         AMTTL = AMTTL + COLDRY(L)+WKL(1,L)
         WVTTL = WVTTL + WKL(1,L)
         DO 4400 IX = 1,MAXXSEC
            IF (IXINDX(IX) .NE. 0) THEN
               WX(IXINDX(IX),L) = COLDRY(L) * WX(IX,L) * 1.E-20
            ENDIF
 4400    CONTINUE
 5000 CONTINUE

      WVSH = (amw*WVTTL)/(amd*AMTTL)
      PWVCM = WVSH*(1.E3*PZ(0))/(1.E2*GRAV)

c      write(23,*) 'INATM'
c      write(23,*) 'cldfra:'
c      write(23,9990) (cldfra(iplon,k),k=1,klev)
c      write(23,*) 'emis:'
c      write(23,9990) emis(iplon)
c      write(23,*) 'pwvcm: ', pwvcm

C  Set spectral surface emissivity for each longwave band.  
      DO 5500 N=1,NBANDS
         SEMISS(N) = emis(iplon)
C         SEMISS(N) = 1.0
 5500 CONTINUE

C  Transfer cloud fraction and cloud optical properties to RRTM variables

      INFLAG = INFLGLW
      ICEFLAG = ICEFLGLW
      LIQFLAG = LIQFLGLW
      DO 7000 L = 1, NLAYERS
C  Cloud fraction
         CLDFRAC(L) = cldfra(iplon,L)
C  Cloud optical depth or cloud water path
         CLDDAT1(L) = cldd1lw(iplon,L)
C  Cloud ice fraction
         CLDDAT2(L) = cldd2lw(ilpon,L)
C  Ice particle effective radius or generalized effective size (microns)
         CLDDAT3(L) = cldd3lw(iplon,L)
C  Liquid droplet effective radius (microns)
         CLDDAT4(L) = cldd4lw(iplon,L)
 7000 CONTINUE

c      write(23,*) ' '
c      write(23,*) 'cldfrac:'
c      write(23,9990) (cldfrac(k),k=1,klev)
c      write(23,*) 'clddat1:'
c      write(23,9990) (clddat1(k),k=1,klev)
c      write(23,*) 'clddat2:'
c      write(23,9990) (clddat2(k),k=1,klev)
c      write(23,*) 'clddat3:'
c      write(23,9990) (clddat3(k),k=1,klev)
c      write(23,*) 'clddat4:'
c      write(23,9990) (clddat4(k),k=1,klev)
c      write(23,*) 'semiss:'
c      write(23,9990) (semiss(n),n=1,nbands)

      RETURN
      END 










