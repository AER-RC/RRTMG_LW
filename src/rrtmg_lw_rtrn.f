C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002-2005, Atmospheric & Environmental Research, Inc. (AER).  |
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
C  --------------------------------------------------------------------------

C***************************************************************************
      SUBROUTINE RTRN
C***************************************************************************
C  RRTM Longwave Radiative Transfer Model
C  Atmospheric and Environmental Research, Inc., Cambridge, MA
C
C  Original version:   E. J. Mlawer, et al. RRTM_V3.0
C  Revision for GCMs:  Michael J. Iacono; October, 2002
C
C  This program calculates the upward fluxes, downward fluxes, and
C  heating rates for an arbitrary clear or cloudy atmosphere.  The input
C  to this program is the atmospheric profile, all Planck function
C  information, and the cloud fraction by layer.  A variable diffusivity 
C  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
C  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
C  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
C  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
C  use of the emissivity angle for the flux integration can cause errors of 
C  1 to 4 W/m2 within cloudy layers.  Clouds are treated with a random 
C  cloud overlap method.
C***************************************************************************

C Parameters
      PARAMETER (MG=16)
      PARAMETER (NGPT=140)
      PARAMETER (MXLAY=203)
      PARAMETER (MXANG = 4)
      PARAMETER (NBANDS = 16)
      PARAMETER (NTBL = 10000,TBLINT=10000.0)

C Input
      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /FEATUREC/  NGC(NBANDS), NGS(NBANDS), NGN(NGPT), NGB(NGPT)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SURFACE/   TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)
      COMMON /PLNKDAT/   PLANKLAY(MXLAY,NBANDS),
     &                   PLANKLEV(0:MXLAY,NBANDS),PLANKBND(NBANDS)
      COMMON /PLANKG/    FRACS(MXLAY,NGPT)
      COMMON /TAUGCOM/   TAUG(MXLAY,NGPT)
      COMMON /RTTBL/     BPADE,
     &                   TAUTBL(0:NTBL),TRANS(0:NTBL),TF(0:NTBL)
      COMMON /PWV/       PWVCM

      COMMON /CVRRTR/    HNAMRTR,HVRRTR

      CHARACTER*18       HNAMRTR,HVRRTR

C Output
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /OUTCLR/    TOTUCLFL(0:MXLAY), TOTDCLFL(0:MXLAY),
     &                   FNETC(0:MXLAY), HTRC(0:MXLAY)

C RRTM Definitions
C Input
C    MXLAY                        ! Maximum number of model layers
C    NGPT                         ! Total number of g-point subintervals
C    NBANDS                       ! Number of longwave spectral bands
C    NCBANDS                      ! Number of spectral bands for clouds
C    SECDIFF                      ! Diffusivity angle
C    WTDIFF                       ! Weight for radiance to flux conversion
C    NLAYERS                      ! Number of model layers (plev+1)
C    PAVEL                        ! Layer pressures (mb)
C    PZ                           ! Level (interface) pressures (mb)
C    TAVEL                        ! Layer temperatures (K)
C    TZ                           ! Level (interface) temperatures(mb)
C    TBOUND                       ! Surface temperature (K)
C    CLDFRAC                      ! Layer cloud fraction
C    TAUCLOUD                     ! Layer cloud optical depth
C    ITR                          ! Integer look-up table index
C    ICLDLYR                      ! Flag for cloudy layers
C    ICLDDN                       ! Flag for cloud in column at any layer
C    SEMISS                       ! Surface emissivities for each band
C    REFLECT                      ! Surface reflectance
C    BPADE                        ! 1/(Pade constant)
C    TAUTBL                       ! Clear sky optical depth look-up table
C    TF                           ! Tau transition function look-up table
C    TRANS                        ! Clear sky transmittance look-up table

C Local
C    ATRANS                       ! Gaseous absorptivity
C    ABSCLD                       ! Cloud absorptivity
C    ATOT                         ! Combined gaseous and cloud absorptivity
C    ODCLR                        ! Clear sky (gaseous) optical depth
C    ODCLD                        ! Cloud optical depth
C    ODTOT                        ! Optical depth of gas and cloud
C    TFACGAS                      ! Gas-only Pade factor, used for Planck fn
C    TFACTOT                      ! Gas and cloud Pade factor, used for Planck fn
C    BBDGAS                       ! Gas-only Planck function for downward rt
C    BBUGAS                       ! Gas-only Planck function for upward rt
C    BBDTOT                       ! Gas and cloud Planck function for downward rt
C    BBUTOT                       ! Gas and cloud Planck function for upward calc.
C    GASSRC                       ! Source radiance due to gas only
C    EFCLFRAC                     ! Effective cloud fraction
C    RADLU                        ! Spectrally summed upward radiance 
C    RADCLRU                      ! Spectrally summed clear sky upward radiance 
C    URAD                         ! Upward radiance by layer
C    CLRURAD                      ! Clear sky upward radiance by layer
C    RADLD                        ! Spectrally summed downward radiance 
C    RADCLRD                      ! Spectrally summed clear sky downward radiance 
C    DRAD                         ! Downward radiance by layer
C    CLRDRAD                      ! Clear sky downward radiance by layer

C Output
C    TOTUFLUX(0:MXLAY)            ! Upward longwave flux (W/m2)
C    TOTDFLUX(0:MXLAY)            ! Downward longwave flux (W/m2)
C    FNET(0:MXLAY)                ! Net longwave flux (W/m2)
C    HTR(0:MXLAY)                 ! Longwave heating rate (K/day)
C    TOTUCLFL(0:MXLAY)            ! Clear sky upward longwave flux (W/m2)
C    TOTDCLFL(0:MXLAY)            ! Clear sky downward longwave flux (W/m2)
C    FNETC(0:MXLAY)               ! Clear sky net longwave flux (W/m2)
C    HTRC(0:MXLAY)                ! Clear sky longwave heating rate (K/day)
C

C Dimensions for radiative transfer
      DIMENSION BBUGAS(MXLAY)
      DIMENSION BBUTOT(MXLAY)
      DIMENSION ATRANS(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD(0:MXLAY),URAD(0:MXLAY)
      DIMENSION UCLFL(0:MXLAY),DCLFL(0:MXLAY)
      DIMENSION CLRDRAD(0:MXLAY),CLRURAD(0:MXLAY)
      DIMENSION ODCLD(MXLAY,NBANDS)
      DIMENSION ATOT(MXLAY)
      DIMENSION ABSCLD(MXLAY,NBANDS)
      DIMENSION EFCLFRAC(MXLAY,NBANDS)
      DIMENSION IPAT(16,0:2)
      DIMENSION SECDIFF(NBANDS),A0(NBANDS),A1(NBANDS),A2(NBANDS)

C Dimensions for cloud 
      DIMENSION ICLDLYR(MXLAY)

C     These arrays indicate the spectral 'region' (used in the 
C     calculation of ice cloud optical depths) corresponding
C     to each spectral band.  See cldprop.f for more details.
      DATA IPAT /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1,
     &           1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5,
     &           1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/


C *** This secant and weight corresponds to the standard diffusivity 
C     angle.  This initial value is redefined below for some bands.
C      DATA SECDIFF /16*1.66/
      DATA WTDIFF /0.5/

      DATA REC_6 /0.166667/

C Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
C and 1.80) as a function of total column water vapor.  The function
C has been defined to minimize flux and cooling rate errors in these bands
C over a wide range of precipitable water values.
      DATA A0 / 1.66, 1.55, 1.58, 1.66, 1.54,1.454, 1.89, 1.33, 
     &         1.668, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66 /
      DATA A1 / 0.00, 0.25, 0.22, 0.00, 0.13,0.446,-0.10, 0.40,
     &        -0.006, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /
      DATA A2 / 0.00,-12.0,-11.7, 0.00,-0.72,-0.243,0.19,-0.062,
     &         0.414, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /

      DO 50 IBND = 1,NBANDS
         IF (IBND.EQ.1 .OR. IBND.EQ.4 .OR. IBND.GE.10) THEN
           SECDIFF(IBND) = 1.66
         ELSE
           SECDIFF(IBND) = A0(IBND) + A1(IBND)*EXP(A2(IBND)*PWVCM)
         ENDIF
50    CONTINUE
      IF (PWVCM.LT.1.0) SECDIFF(6) = 1.80
      IF (PWVCM.GT.7.1) SECDIFF(7) = 1.50


      HVRRTR = '$Revision$'

      URAD(0) = 0.0
      DRAD(0) = 0.0
      TOTUFLUX(0) = 0.0
      TOTDFLUX(0) = 0.0
      CLRURAD(0) = 0.0
      CLRDRAD(0) = 0.0
      TOTUCLFL(0) = 0.0
      TOTDCLFL(0) = 0.0

      DO 200 LAY = 1, NLAYERS
         URAD(LAY) = 0.0
         DRAD(LAY) = 0.0
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
         CLRURAD(LAY) = 0.0
         CLRDRAD(LAY) = 0.0
         TOTUCLFL(LAY) = 0.0
         TOTDCLFL(LAY) = 0.0

         DO 100 IB = 1, NCBANDS
            IF (CLDFRAC(LAY) .GE. 1.E-6) THEN
               ODCLD(LAY,IB) = SECDIFF(IB) * TAUCLOUD(LAY,IB)
               TRANSCLD = EXP(-ODCLD(LAY,IB))
               ABSCLD(LAY,IB) = 1. - TRANSCLD
               EFCLFRAC(LAY,IB) = ABSCLD(LAY,IB) * 
     &                 CLDFRAC(LAY)
               ICLDLYR(LAY) = 1
            ELSE
               ODCLD(LAY,IB) = 0.0
               ABSCLD(LAY,IB) = 0.0
               EFCLFRAC(LAY,IB) = 0.0
               ICLDLYR(LAY) = 0
            ENDIF
 100     CONTINUE
 200  CONTINUE

      IGC = 1
C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND

C *** Reinitialize g-point counter for each band if output for each
C     band is requested.
         IF (IOUT.GT.0.AND.IBAND.GE.2) IGC = NGS(IBAND-1)+1
         ICLDDN = 0
         IF (NCBANDS .EQ. 1) THEN
            IB = IPAT(IBAND,0)
         ELSEIF (NCBANDS .EQ.  5) THEN
            IB = IPAT(IBAND,1)
         ELSEIF (NCBANDS .EQ. 16) THEN
            IB = IPAT(IBAND,2)
         ENDIF

C ***    Loop over g-channels.
 1000    CONTINUE

C ***    Radiative transfer starts here.
         RADLD = 0.
         RADCLRD = 0.

C ***    DOWNWARD RADIATIVE TRANSFER loop.  

         DO 2500 LEV = NLAYERS, 1, -1
               PLFRAC = FRACS(LEV,IGC)
               BLAY = PLANKLAY(LEV,IBAND)
               DPLANKUP = PLANKLEV(LEV,IBAND) - BLAY
               DPLANKDN = PLANKLEV(LEV-1,IBAND) - BLAY
               ODEPTH = SECDIFF(IBAND) * TAUG(LEV,IGC)
               IF (ODEPTH .LT. 0.0) ODEPTH = 0.0
C Cloudy layer
               IF (ICLDLYR(LEV).EQ.1) THEN
                  ICLDDN = 1
                  ODTOT = ODEPTH + ODCLD(LEV,IB)
                  IF (ODTOT .LT. 0.06) THEN
                     ATRANS(LEV) = ODEPTH - 0.5*ODEPTH*ODEPTH
                     ODEPTH_REC = REC_6*ODEPTH
                     GASSRC = PLFRAC*(BLAY+DPLANKDN*ODEPTH_REC)
     &                    *ATRANS(LEV)

                     ATOT(LEV) =  ODTOT - 0.5*ODTOT*ODTOT
                     ODTOT_REC = REC_6*ODTOT
                     BBDTOT =  PLFRAC * (BLAY+DPLANKDN*ODTOT_REC)
                     RADLD = RADLD - RADLD * (ATRANS(LEV) +
     &                    EFCLFRAC(LEV,IB) * (1. - ATRANS(LEV))) +
     &                    GASSRC + CLDFRAC(LEV) * 
     &                    (BBDTOT * ATOT(LEV) - GASSRC)
                     DRAD(LEV-1) = DRAD(LEV-1) + RADLD
                  
                     BBUGAS(LEV) =  PLFRAC *
     &                    (BLAY+DPLANKUP*ODEPTH_REC)
                     BBUTOT(LEV) =  PLFRAC * 
     &                    (BLAY+DPLANKUP*ODTOT_REC)

                  ELSEIF (ODEPTH .LE. 0.06) THEN
                     ATRANS(LEV) = ODEPTH - 0.5*ODEPTH*ODEPTH
                     ODEPTH_REC = REC_6*ODEPTH
                     GASSRC = PLFRAC*(BLAY+DPLANKDN*ODEPTH_REC)
     &                    *ATRANS(LEV)

                     ODTOT = ODEPTH + ODCLD(LEV,IB)
                     TBLIND = ODTOT/(BPADE+ODTOT)
                     ITTOT = TBLINT*TBLIND + 0.5
                     TFACTOT = TF(ITTOT)
                     BBDTOT = PLFRAC * (BLAY + TFACTOT*DPLANKDN)
                     ATOT(LEV) = 1. - TRANS(ITTOT)

                     RADLD = RADLD - RADLD * (ATRANS(LEV) +
     &                    EFCLFRAC(LEV,IB) * (1. - ATRANS(LEV))) +
     &                    GASSRC + CLDFRAC(LEV) * 
     &                    (BBDTOT * ATOT(LEV) - GASSRC)
                     DRAD(LEV-1) = DRAD(LEV-1) + RADLD

                     BBUGAS(LEV) = PLFRAC * 
     &                 (BLAY + DPLANKUP*ODEPTH_REC)
                     BBUTOT(LEV) = PLFRAC * 
     &                 (BLAY + TFACTOT * DPLANKUP)

                  ELSE

                     TBLIND = ODEPTH/(BPADE+ODEPTH)
                     ITGAS = TBLINT*TBLIND+0.5
                     ODEPTH = TAUTBL(ITGAS)
                     ATRANS(LEV) = 1. - TRANS(ITGAS)
                     TFACGAS = TF(ITGAS)
                     GASSRC = ATRANS(LEV) * PLFRAC * 
     &                    (BLAY + TFACGAS*DPLANKDN)

                     ODTOT = ODEPTH + ODCLD(LEV,IB)
                     TBLIND = ODTOT/(BPADE+ODTOT)
                     ITTOT = TBLINT*TBLIND + 0.5
                     TFACTOT = TF(ITTOT)
                     BBDTOT = PLFRAC * (BLAY + TFACTOT*DPLANKDN)
                     ATOT(LEV) = 1. - TRANS(ITTOT)

                  RADLD = RADLD - RADLD * (ATRANS(LEV) +
     &              EFCLFRAC(LEV,IB) * (1. - ATRANS(LEV))) +
     &              GASSRC + CLDFRAC(LEV) * 
     &              (BBDTOT * ATOT(LEV) - GASSRC)
                  DRAD(LEV-1) = DRAD(LEV-1) + RADLD
                  BBUGAS(LEV) = PLFRAC * 
     &                 (BLAY + TFACGAS * DPLANKUP)
                  BBUTOT(LEV) = PLFRAC * 
     &                 (BLAY + TFACTOT * DPLANKUP)
                  ENDIF
C Clear layer
               ELSE
                  IF (ODEPTH .LE. 0.06) THEN
                     ATRANS(LEV) = ODEPTH-0.5*ODEPTH*ODEPTH
                     ODEPTH = REC_6*ODEPTH
                     BBD = PLFRAC*(BLAY+DPLANKDN*ODEPTH)
                     BBUGAS(LEV) = PLFRAC*
     &                    (BLAY+DPLANKUP*ODEPTH)
                  ELSE
                     TBLIND = ODEPTH/(BPADE+ODEPTH)
                     ITR = TBLINT*TBLIND+0.5
                     TRANSC = TRANS(ITR)
                     ATRANS(LEV) = 1.-TRANSC
                     TAUSFAC = TF(ITR)
                     BBD = PLFRAC*(BLAY+TAUSFAC*DPLANKDN)
                     BBUGAS(LEV) = PLFRAC * 
     &                    (BLAY + TAUSFAC * DPLANKUP)
                  ENDIF   
                  RADLD = RADLD + (BBD-RADLD)*ATRANS(LEV)
                  DRAD(LEV-1) = DRAD(LEV-1) + RADLD
               ENDIF
C    Set clear sky stream to total sky stream as long as layers
C    remain clear.  Streams diverge when a cloud is reached (ICLDDN=1),
C    and clear sky stream must be computed separately from that point.
                  IF (ICLDDN.EQ.1) THEN
                     RADCLRD = RADCLRD + (BBD-RADCLRD) * ATRANS(LEV) 
                     CLRDRAD(LEV-1) = CLRDRAD(LEV-1) + RADCLRD
                  ELSE
                     RADCLRD = RADLD
                     CLRDRAD(LEV-1) = DRAD(LEV-1)
                  ENDIF
 2500       CONTINUE

C ***    SPECTRAL EMISSIVITY & REFLECTANCE
C    Include the contribution of spectrally varying longwave emissivity
C    and reflection from the surface to the upward radiative transfer.
C    Note: Spectral and Lambertian reflection are identical for the
C    diffusivity angle flux integration used here.

         RAD0 = FRACS(1,IGC) * PLANKBND(IBAND)
C    Add in specular reflection of surface downward radiance.
         REFLECT = 1. - SEMISS(IBAND)
         RADLU = RAD0 + REFLECT * RADLD
         RADCLRU = RAD0 + REFLECT * RADCLRD

C ***    UPWARD RADIATIVE TRANSFER loop.
         URAD(0) = URAD(0) + RADLU
         CLRURAD(0) = CLRURAD(0) + RADCLRU

         DO 2600 LEV = 1, NLAYERS
C Cloudy layer
            IF (ICLDLYR(LEV) .EQ. 1) THEN
               GASSRC = BBUGAS(LEV) * ATRANS(LEV)
               RADLU = RADLU - RADLU * (ATRANS(LEV) +
     &              EFCLFRAC(LEV,IB) * (1. - ATRANS(LEV))) +
     &              GASSRC + CLDFRAC(LEV) * 
     &              (BBUTOT(LEV) * ATOT(LEV) - GASSRC)
               URAD(LEV) = URAD(LEV) + RADLU
C Clear layer
            ELSE
               RADLU = RADLU + (BBUGAS(LEV)-RADLU)*ATRANS(LEV)
               URAD(LEV) = URAD(LEV) + RADLU
            ENDIF
C    Set clear sky stream to total sky stream as long as all layers
C    are clear (ICLDDN=0).  Streams must be calculated separately at 
C    all layers when a cloud is present (ICLDDN=1), because surface 
C    reflectance is different for each stream.
               IF (ICLDDN.EQ.1) THEN
                  RADCLRU = RADCLRU + (BBUGAS(LEV)-RADCLRU)*ATRANS(LEV) 
                  CLRURAD(LEV) = CLRURAD(LEV) + RADCLRU
               ELSE
                  RADCLRU = RADLU
                  CLRURAD(LEV) = URAD(LEV)
               ENDIF
 2600    CONTINUE

         IGC = IGC + 1
         IF (IGC .LE. NGS(IBAND)) GO TO 1000

C ***    Process longwave output from band for total and clear streams.
C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = URAD(LEV)*WTDIFF
            DFLUX(LEV) = DRAD(LEV)*WTDIFF
            URAD(LEV) = 0.0
            DRAD(LEV) = 0.0
            TOTUFLUX(LEV) = TOTUFLUX(LEV) + UFLUX(LEV) * DELWAVE(IBAND)
            TOTDFLUX(LEV) = TOTDFLUX(LEV) + DFLUX(LEV) * DELWAVE(IBAND)
            UCLFL(LEV) = CLRURAD(LEV)*WTDIFF
            DCLFL(LEV) = CLRDRAD(LEV)*WTDIFF
            CLRURAD(LEV) = 0.0
            CLRDRAD(LEV) = 0.0
            TOTUCLFL(LEV) = TOTUCLFL(LEV) + UCLFL(LEV) * DELWAVE(IBAND)
            TOTDCLFL(LEV) = TOTDCLFL(LEV) + DCLFL(LEV) * DELWAVE(IBAND)
 5000    CONTINUE
 6000 CONTINUE

      TOTUFLUX(0) = TOTUFLUX(0) * FLUXFAC
      TOTDFLUX(0) = TOTDFLUX(0) * FLUXFAC
      FNET(0) = TOTUFLUX(0) - TOTDFLUX(0)
      TOTUCLFL(0) = TOTUCLFL(0) * FLUXFAC
      TOTDCLFL(0) = TOTDCLFL(0) * FLUXFAC
      FNETC(0) = TOTUCLFL(0) - TOTDCLFL(0)

      DO 7000 LEV = 1, NLAYERS
         TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
         TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
         FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
         TOTUCLFL(LEV) = TOTUCLFL(LEV) * FLUXFAC
         TOTDCLFL(LEV) = TOTDCLFL(LEV) * FLUXFAC
         FNETC(LEV) = TOTUCLFL(LEV) - TOTDCLFL(LEV)
         L = LEV - 1

C        Calculate Heating Rates.
         HTR(L)=HEATFAC*(FNET(L)-FNET(LEV))/(PZ(L)-PZ(LEV)) 
         HTRC(L)=HEATFAC*(FNETC(L)-FNETC(LEV))/(PZ(L)-PZ(LEV)) 
 7000 CONTINUE
      HTR(NLAYERS) = 0.0
      HTRC(NLAYERS) = 0.0

 9000 CONTINUE

      RETURN
      END   



