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

*******************************************************************************
*                                                                             *
*                  Optical depths developed for the                           *
*                                                                             *
*                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
*                                                                             *
*                                                                             *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
*                        131 HARTWELL AVENUE                                  *
*                        LEXINGTON, MA 02421                                  *
*                                                                             *
*                                                                             *
*                           ELI J. MLAWER                                     * 
*                         JENNIFER DELAMERE                                   * 
*                         STEVEN J. TAUBMAN                                   *
*                         SHEPARD A. CLOUGH                                   *
*                                                                             *
*                                                                             *
*                                                                             *
*                                                                             *
*                       email:  mlawer@aer.com                                *
*                       email:  jdelamer@aer.com                              *
*                                                                             *
*        The authors wish to acknowledge the contributions of the             *
*        following people:  Karen Cady-Pereira, Patrick D. Brown,             *  
*        Michael J. Iacono, Ronald E. Farren, Luke Chen, Robert Bergstrom.    *
*                                                                             *
*******************************************************************************
*                                                                             *
*  Revision for g-point reduction: Michael J. Iacono, AER, Inc.               *
*                                                                             *
*******************************************************************************
*     TAUMOL                                                                  *
*                                                                             *
*     This file contains the subroutines TAUGBn (where n goes from            *
*     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
*     per g-value and layer for band n.                                       *
*                                                                             *
*  Output:  optical depths (unitless)                                         *
*           fractions needed to compute Planck functions at every layer       *
*               and g-value                                                   *
*                                                                             *
*     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
*     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
*                                                                             *
*  Input                                                                      *
*                                                                             *
*     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
*     COMMON /PRECISE/  ONEMINUS                                              *
*     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
*     &                 PZ(0:MXLAY),TZ(0:MXLAY)                               *
*     COMMON /PROFDATA/ LAYTROP,                                              *
*    &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),             *
*    &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),             *
*    &                  COLO2(MXLAY)
*     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
*    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
*     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
*     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
*                                                                             *
*     Description:                                                            *
*     NG(IBAND) - number of g-values in band IBAND                            *
*     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
*                   atmospheres that are stored for band IBAND per            *
*                   pressure level and temperature.  Each of these            *
*                   atmospheres has different relative amounts of the         *
*                   key species for the band (i.e. different binary           *
*                   species parameters).                                      *
*     NSPB(IBAND) - same for upper atmosphere                                 *
*     ONEMINUS - since problems are caused in some cases by interpolation     *
*                parameters equal to or greater than 1, for these cases       *
*                these parameters are set to this value, slightly < 1.        *
*     PAVEL - layer pressures (mb)                                            *
*     TAVEL - layer temperatures (degrees K)                                  *
*     PZ - level pressures (mb)                                               *
*     TZ - level temperatures (degrees K)                                     *
*     LAYTROP - layer at which switch is made from one combination of         *
*               key species to another                                        *
*     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
*               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
*               respectively (molecules/cm**2)                                *

*     FACij(LAY) - for layer LAY, these are factors that are needed to        *
*                  compute the interpolation factors that multiply the        *
*                  appropriate reference k-values.  A value of 0 (1) for      *
*                  i,j indicates that the corresponding factor multiplies     *
*                  reference k-value for the lower (higher) of the two        *
*                  appropriate temperatures, and altitudes, respectively.     *
*     JP - the index of the lower (in altitude) of the two appropriate        *
*          reference pressure levels needed for interpolation                 *
*     JT, JT1 - the indices of the lower of the two appropriate reference     *
*               temperatures needed for interpolation (for pressure           *
*               levels JP and JP+1, respectively)                             *
*     SELFFAC - scale factor needed for water vapor self-continuum, equals    *
*               (water vapor density)/(atmospheric density at 296K and        *
*               1013 mb)                                                      *
*     SELFFRAC - factor needed for temperature interpolation of reference     *
*                water vapor self-continuum data                              *
*     INDSELF - index of the lower of the two appropriate reference           *
*               temperatures needed for the self-continuum interpolation      *
*     FORFAC  - scale factor needed for water vapor foreign-continuum.        *
*     FORFRAC - factor needed for temperature interpolation of reference      *
*                water vapor foreign-continuum data                           *
*     INDFOR  - index of the lower of the two appropriate reference           *
*               temperatures needed for the foreign-continuum interpolation   *
*                                                                             *
*  Data input                                                                 *
*     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG),*
*                 FORREF(4,MG), KA_M'MGAS', KB_M'MGAS'                        *
*        (note:  n is the band number,'MGAS' is the species name of the minor *
*         gas)                                                                *
*                                                                             *
*     Description:                                                            *
*     KA - k-values for low reference atmospheres (key-species only)          *
*          (units: cm**2/molecule)                                            *
*     KB - k-values for high reference atmospheres (key-species only)         *
*          (units: cm**2/molecule)                                            *
*     KA_M'MGAS' - k-values for low reference atmosphere minor species        *
*          (units: cm**2/molecule)                                            *
*     KB_M'MGAS' - k-values for high reference atmosphere minor species       *
*          (units: cm**2/molecule)                                            *
*     SELFREF - k-values for water vapor self-continuum for reference         *
*               atmospheres (used below LAYTROP)                              *
*               (units: cm**2/molecule)                                       *
*     FORREF  - k-values for water vapor foreign-continuum for reference      *
*               atmospheres (used below/above LAYTROP)                        *
*               (units: cm**2/molecule)                                       *
*                                                                             *
*     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
*     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
*                                                                             *
*******************************************************************************

      SUBROUTINE TAUGB1

C     Written by Eli J. Mlawer, Atmospheric & Environmental Research.
C     Revised by Michael J. Iacono, Atmospheric & Environmental Research.

C     BAND 1:  10-350 cm-1 (low key - H2O; low minor - N2)
C                          (high key - H2O; high minor - N2)

C     NOTE: Previous versions of RRTM BAND 1: 
C           10-250 cm-1 (low - H2O; high - H2O)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG1=10)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K1C/      KA(5,13,NG1), KB(5,13:59,NG1),FORREF(4,NG1), 
     &                  SELFREF(10,NG1),KA_MN2(19,NG1),KB_MN2(19,NG1)
      COMMON /PF1C/     FRACREFA(NG1), FRACREFB(NG1)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(65,NG1),ABSB(235,NG1)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA, KB, KA_MN2, KB_MN2, MINORFRAC

C Minor gas mapping levels:
C     LOWER - N2, P = 142.5490 mbar, T = 215.70 K
C     UPPER - N2, P = 142.5490 mbar, T = 215.70 K

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(1) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(1) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         PP = PAVEL(LAY)
         CORRADJ =  1.
         IF (PP .LT. 250.) THEN
            CORRADJ = 1. - 0.15 * (250.-PP) / 154.4
         ENDIF

         SCALEN2 = COLBRD(LAY) * SCALEMINORN2(LAY)
         DO 2000 IG = 1, NG1
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUN2 = SCALEN2*(KA_MN2(INDM,IG) +  
     &           MINORFRAC(LAY) *
     &           (KA_MN2(INDM+1,IG) - KA_MN2(INDM,IG)))
            TAUG(LAY,IG) = CORRADJ * (COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))  
     &           + TAUSELF + TAUFOR
     &           + TAUN2)
             FRACS(LAY,IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(1) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(1) + 1
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         PP = PAVEL(LAY)
         CORRADJ =  1. - 0.15 * (PP / 95.6)

         SCALEN2 = COLBRD(LAY) * SCALEMINORN2(LAY)
         DO 3000 IG = 1, NG1
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUN2 = SCALEN2*(KB_MN2(INDM,IG) +  
     &           MINORFRAC(LAY) *
     &           (KB_MN2(INDM+1,IG) - KB_MN2(INDM,IG)))
            TAUG(LAY,IG) = CORRADJ * (COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))   
     &           + TAUFOR
     &           + TAUN2)
            FRACS(LAY,IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB2

C     BAND 2:  350-500 cm-1 (low key - H2O; high key - H2O)

C     NOTE: Previous version of RRTM BAND 2: 
C           250 - 500 cm-1 (low - H2O; high - H2O)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, NBANDS=16)
      PARAMETER (NGPT=140, NG2=12, NGS1=10)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K2C/      KA(5,13,NG2), KB(5,13:59,NG2), FORREF(4,NG2), 
     &                  SELFREF(10,NG2)
      COMMON /PF2C/     FRACREFA(NG2), FRACREFB(NG2) 

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(65,NG2),ABSB(235,NG2)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.
      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(2) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(2) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         PP = PAVEL(LAY)
         CORRADJ = 1. - .05 * (PP - 100.) / 900.
         DO 2000 IG = 1, NG2
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUG(LAY,NGS1+IG) = CORRADJ * (COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR)
            FRACS(LAY,NGS1+IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(2) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(2) + 1
         INDF = INDFOR(LAY)
         DO 3000 IG = 1, NG2
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUG(LAY,NGS1+IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + TAUFOR
            FRACS(LAY,NGS1+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB3

C     BAND 3:  500-630 cm-1 (low key - H2O,CO2; low minor - n2o)
C                           (high key - H2O,CO2; high minor - n2o)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, NBANDS=16)
      PARAMETER (NGPT=140, NG3=16, NGS2=22)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)  
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K3C/   KA(9,5,13,NG3), KB(5,5,13:59,NG3), FORREF(4,NG3),
     &                  SELFREF(10,NG3), KA_MN2O(9,19,NG3), 
     &                  KB_MN2O(5,19,NG3)
      COMMON /PF3C/     FRACREFA(NG3,9), FRACREFB(NG3,5) 

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      REAL KA,KB
      REAL KA_MN2O, KB_MN2O, MINORFRAC
      REAL N2OM1,N2OM2

      DIMENSION ABSA(585,NG3),ABSB(1175,NG3)

C Minor gas mapping levels:
C     LOWER - N2O, P = 706.272 mbar, T = 278.94 K
C     UPPER - N2O, P = 95.58 mbar, T = 215.7 K

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C     P = 212.725 mb
      REFRAT_PLANCK_A = CHI_MLS(1,9)/CHI_MLS(2,9)

C     P = 95.58 mb
      REFRAT_PLANCK_B = CHI_MLS(1,13)/CHI_MLS(2,13)

C     P = 706.270mb
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(2,3)

C     P = 95.58 mb 
      REFRAT_M_B = CHI_MLS(1,13)/CHI_MLS(2,13)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water vapor 
C     self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)        

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)
         FMN2OMF = MINORFRAC(LAY)*FMN2O
c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/COLDRY(LAY)
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(3) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(3) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG3
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)
               TAUG(LAY,NGS2+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
               FRACS(LAY,NGS2+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000       CONTINUE
         ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG3
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)
               TAUG(LAY,NGS2+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O 
               FRACS(LAY,NGS2+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010          CONTINUE
         ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG3
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG) - 
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG) - 
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY) *
     &              (N2OM2 - N2OM1)
               TAUG(LAY,NGS2+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O   
               FRACS(LAY,NGS2+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
        ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 4.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS1) * FAC01(LAY)
         FAC011 = (1. - FS1) * FAC11(LAY)
         FAC101 = FS1 * FAC01(LAY)
         FAC111 = FS1 * FAC11(LAY)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_B*COLCO2(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 4.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)
         FMN2OMF = MINORFRAC(LAY)*FMN2O
c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/COLDRY(LAY)
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_B*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 4.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(3) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(3) + JS1
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG3
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            N2OM1 = KB_MN2O(JMN2O,INDM,IG) + FMN2O*
     &           (KB_MN2O(JMN2O+1,INDM,IG)-KB_MN2O(JMN2O,INDM,IG))
            N2OM2 = KB_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &           (KB_MN2O(JMN2O+1,INDM+1,IG)-KB_MN2O(JMN2O,INDM+1,IG))
            ABSN2O = N2OM1 + MINORFRAC(LAY) * (N2OM2 - N2OM1)
            TAUG(LAY,NGS2+IG) = SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &          FAC100 * ABSB(IND0+1,IG) +
     &          FAC010 * ABSB(IND0+5,IG) +
     &          FAC110 * ABSB(IND0+6,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSB(IND1,IG) + 
     &          FAC101 * ABSB(IND1+1,IG) +
     &          FAC011 * ABSB(IND1+5,IG) +
     &          FAC111 * ABSB(IND1+6,IG)) 
     &          + TAUFOR
     &          + ADJCOLN2O*ABSN2O            
            FRACS(LAY,NGS2+IG) = FRACREFB(IG,JPL) + FPL *
     &          (FRACREFB(IG,JPL+1)-FRACREFB(IG,JPL))
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB4

C     BAND 4:  630-700 cm-1 (low key - H2O,CO2; high key - O3,CO2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG4=14, NGS3=38)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K4C/   KA(9,5,13,NG4), KB(5,5,13:59,NG4), FORREF(4,NG4),
     &                  SELFREF(10,NG4)
      COMMON /PF4C/     FRACREFA(NG4,9), FRACREFB(NG4,5) 

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG4),ABSB(1175,NG4)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     P =   142.5940 mb
      REFRAT_PLANCK_A = CHI_MLS(1,11)/CHI_MLS(2,11)

C     P = 95.58350 mb
      REFRAT_PLANCK_B = CHI_MLS(3,13)/CHI_MLS(2,13)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water 
C     vapor self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(4) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(4) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG4
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS3+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS3+IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
         ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG4
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS3+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS3+IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
         ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG4
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS3+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS3+IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
        ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         SPECCOMB = COLO3(LAY) + RAT_O3CO2(LAY)*COLCO2(LAY)
         SPECPARM = COLO3(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLO3(LAY) + RAT_O3CO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLO3(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 4.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS1) * FAC01(LAY)
         FAC011 = (1. - FS1) * FAC11(LAY)
         FAC101 = FS1 * FAC01(LAY)
         FAC111 = FS1 * FAC11(LAY)

         SPECCOMB_PLANCK = COLO3(LAY)+REFRAT_PLANCK_B*COLCO2(LAY)
         SPECPARM_PLANCK = COLO3(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 4.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(4) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(4) + JS1

         DO 3000 IG = 1, NG4
            TAUG(LAY,NGS3+IG) =  SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &          FAC100 * ABSB(IND0+1,IG) +
     &          FAC010 * ABSB(IND0+5,IG) +
     &          FAC110 * ABSB(IND0+6,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSB(IND1,IG) + 
     &          FAC101 * ABSB(IND1+1,IG) +
     &          FAC011 * ABSB(IND1+5,IG) +
     &          FAC111 * ABSB(IND1+6,IG)) 
            FRACS(LAY,NGS3+IG) = FRACREFB(IG,JPL) + FPL *
     &          (FRACREFB(IG,JPL+1)-FRACREFB(IG,JPL))
 3000    CONTINUE

C EMPIRICAL MODIFICATION TO CODE TO IMPROVE STRATOSPHERIC COOLING RATES
C FOR CO2.  REVISED TO APPLY WEIGHTING FOR G-POINT REDUCTION IN THIS BAND.

         TAUG(LAY,NGS3+8)=TAUG(LAY,NGS3+8)*0.92
         TAUG(LAY,NGS3+9)=TAUG(LAY,NGS3+9)*0.88
         TAUG(LAY,NGS3+10)=TAUG(LAY,NGS3+10)*1.07
         TAUG(LAY,NGS3+11)=TAUG(LAY,NGS3+11)*1.1
         TAUG(LAY,NGS3+12)=TAUG(LAY,NGS3+12)*0.99
         TAUG(LAY,NGS3+13)=TAUG(LAY,NGS3+13)*0.88
         TAUG(LAY,NGS3+14)=TAUG(LAY,NGS3+14)*0.943

 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB5

C     BAND 5:  700-820 cm-1 (low key - H2O,CO2; low minor - O3, CCL4)
C                           (high key - O3,CO2)

      PARAMETER (MG=16, MXLAY=203, MAXXSEC=4, NBANDS=16)
      PARAMETER (NGPT=140, NG5=16, NGS4=52)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY),INDSELF(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K5C/      KA(9,5,13,NG5), KB(5,5,13:59,NG5),
     &              FORREF(4,NG5), SELFREF(10,NG5), KA_MO3(9,19,NG5),
     &                  CCL4(NG5)
      COMMON /PF5C/     FRACREFA(NG5,9), FRACREFB(NG5,5)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG5),ABSB(1175,NG5)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB
      REAL KA_MO3, MINORFRAC
      REAL O3M1, O3M2

C Minor gas mapping level :
C     LOWER - O3, P = 317.34 mbar, T = 240.77 K
C     LOWER - CCL4

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb
      REFRAT_PLANCK_A = CHI_MLS(1,5)/CHI_MLS(2,5)

C     P = 0.2369 mb
      REFRAT_PLANCK_B = CHI_MLS(3,43)/CHI_MLS(2,43)

C     P = 317.3480
      REFRAT_M_A = CHI_MLS(1,7)/CHI_MLS(2,7)

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the 
C     water vapor self-continuum and foreign continuum is 
C     interpolated (in temperature) separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MO3 = COLH2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MO3 = COLH2O(LAY)/SPECCOMB_MO3
         IF (SPECPARM_MO3 .GE. ONEMINUS) SPECPARM_MO3 = ONEMINUS
         SPECMULT_MO3 = 8.*SPECPARM_MO3
         JMO3 = 1 + INT(SPECMULT_MO3)
         FMO3 = AMOD(SPECMULT_MO3,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(5) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(5) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG5
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))

               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,NGS4+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + WX(1,LAY) * CCL4(IG)
               FRACS(LAY,NGS4+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG5
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,NGS4+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF+ TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + WX(1,LAY) * CCL4(IG)
                FRACS(LAY,NGS4+IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG5
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               O3M1 = KA_MO3(JMO3,INDM,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM,IG)-KA_MO3(JMO3,INDM,IG))
               O3M2 = KA_MO3(JMO3,INDM+1,IG) + FMO3*
     &              (KA_MO3(JMO3+1,INDM+1,IG)-KA_MO3(JMO3,INDM+1,IG))
               ABSO3 = O3M1 + MINORFRAC(LAY)*(O3M2-O3M1)
               TAUG(LAY,NGS4+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ABSO3*COLO3(LAY)
     &              + WX(1,LAY) * CCL4(IG)
            FRACS(LAY,NGS4+IG) = FRACREFA(IG,JPL) + FPL *
     &          (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         SPECCOMB = COLO3(LAY) + RAT_O3CO2(LAY)*COLCO2(LAY)
         SPECPARM = COLO3(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 4.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLO3(LAY) + RAT_O3CO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLO3(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 4.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         FAC000 = (1. - FS) * FAC00(LAY)
         FAC010 = (1. - FS) * FAC10(LAY)
         FAC100 = FS * FAC00(LAY)
         FAC110 = FS * FAC10(LAY)
         FAC001 = (1. - FS1) * FAC01(LAY)
         FAC011 = (1. - FS1) * FAC11(LAY)
         FAC101 = FS1 * FAC01(LAY)
         FAC111 = FS1 * FAC11(LAY)

         SPECCOMB_PLANCK = COLO3(LAY)+REFRAT_PLANCK_B*COLCO2(LAY)
         SPECPARM_PLANCK = COLO3(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 4.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(5) + JS
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(5) + JS1
         
         DO 3000 IG = 1, NG5
            TAUG(LAY,NGS4+IG) =  SPECCOMB * 
     &          (FAC000 * ABSB(IND0,IG) +
     &          FAC100 * ABSB(IND0+1,IG) +
     &          FAC010 * ABSB(IND0+5,IG) +
     &          FAC110 * ABSB(IND0+6,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSB(IND1,IG) + 
     &          FAC101 * ABSB(IND1+1,IG) +
     &          FAC011 * ABSB(IND1+5,IG) +
     &          FAC111 * ABSB(IND1+6,IG)) 
     &          + WX(1,LAY) * CCL4(IG)
            FRACS(LAY,NGS4+IG) = FRACREFB(IG,JPL) + FPL *
     &          (FRACREFB(IG,JPL+1)-FRACREFB(IG,JPL))
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB6

C     BAND 6:  820-980 cm-1 (low key - H2O; low minor - CO2)
C                           (high key - nothing; high minor - CFC11, CFC12)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, MAXXSEC=4, NBANDS=16)
      PARAMETER (NGPT=140, NG6=8, NGS5=68)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K6C/      KA(5,13,NG6), FORREF(4,NG6), SELFREF(10,NG6), 
     &                  KA_MCO2(19,NG6), CFC11ADJ(NG6), CFC12(NG6)
      COMMON /PF6C/     FRACREFA(NG6) 

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(65,NG6)
      REAL KA_MCO2, MINORFRAC

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C Minor gas mapping level:
C     LOWER - CO2, P = 706.2720 mb, T = 294.2 K
C     UPPER - CFC11, CFC12

C     Compute the optical depth by interpolating in ln(pressure) and
C     temperature. The water vapor self-continuum and foreign continuum
C     is interpolated (in temperature) separately.  

      HVRTAU = '$Revision$'
 
      DO 2500 LAY = 1, LAYTROP

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.77
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(6) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(6) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 2000 IG = 1, NG6
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG)))
            ABSCO2 =  (KA_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MCO2(INDM+1,IG) - KA_MCO2(INDM,IG)))
            TAUG(LAY,NGS5+IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR
     &           + ADJCOLCO2 * ABSCO2
     &           + WX(2,LAY) * CFC11ADJ(IG)
     &           + WX(3,LAY) * CFC12(IG)
            FRACS(LAY,NGS5+IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

C     Nothing important goes on above LAYTROP in this band.
      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG6
            TAUG(LAY,NGS5+IG) = 0.0 
     &           + WX(2,LAY) * CFC11ADJ(IG)
     &           + WX(3,LAY) * CFC12(IG)
            FRACS(LAY,NGS5+IG) = FRACREFA(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB7

C     BAND 7:  980-1080 cm-1 (low key - H2O,O3; low minor - CO2)
C                            (high key - O3; high minor - CO2)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, NBANDS=16)
      PARAMETER (NGPT=140, NG7=12, NGS6=76)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY) 
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K7C/      KA(9,5,13,NG7), KB(5,13:59,NG7), FORREF(4,NG7),
     &              SELFREF(10,NG7), KA_MCO2(9,19,NG7), KB_MCO2(19,NG7)
      COMMON /PF7C/     FRACREFA(NG7,9), FRACREFB(NG7)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG7),ABSB(235,NG7)
      REAL KA_MCO2, KB_MCO2, MINORFRAC

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C Minor gas mapping level :
C     LOWER - CO2, P = 706.2620 mbar, T= 278.94 K
C     UPPER - CO2, P = 12.9350 mbar, T = 234.01 K

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.

C     P = 706.2620 mb
      REFRAT_PLANCK_A = CHI_MLS(1,3)/CHI_MLS(3,3)

C     P = 706.2720 mb
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(3,3)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately. 

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OO3(LAY)*COLO3(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OO3_1(LAY)*COLO3(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MCO2 = COLH2O(LAY) + REFRAT_M_A*COLO3(LAY)
         SPECPARM_MCO2 = COLH2O(LAY)/SPECCOMB_MCO2
         IF (SPECPARM_MCO2 .GE. ONEMINUS) SPECPARM_MCO2 = ONEMINUS
         SPECMULT_MCO2 = 8.*SPECPARM_MCO2

         JMCO2 = 1 + INT(SPECMULT_MCO2)
         FMCO2 = AMOD(SPECMULT_MCO2,1.0)

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 3.0+(RATCO2-3.0)**0.79
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLO3(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(7) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(7) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG7

               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY) * (CO2M2 - CO2M1)
               TAUG(LAY,NGS6+IG) = SPECCOMB *
     &          (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC200 * ABSA(IND0+2,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG) +
     &          FAC210 * ABSA(IND0+11,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC201 * ABSA(IND1+2,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG) +
     &          FAC211 * ABSA(IND1+11,IG)) 
     &          + TAUSELF + TAUFOR
     &          + ADJCOLCO2*ABSCO2
               FRACS(LAY,NGS6+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG7
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY) * (CO2M2 - CO2M1)
               TAUG(LAY,NGS6+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
               FRACS(LAY,NGS6+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG7
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + MINORFRAC(LAY) * (CO2M2 - CO2M1)
               TAUG(LAY,NGS6+IG) = SPECCOMB * 
     &         (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG)) 
     &          + TAUSELF + TAUFOR
     &          + ADJCOLCO2*ABSCO2
               FRACS(LAY,NGS6+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020    CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.79
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(7) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(7) + 1
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG7
            ABSCO2 = KB_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MCO2(INDM+1,IG) - KB_MCO2(INDM,IG))
            TAUG(LAY,NGS6+IG) = COLO3(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + ADJCOLCO2 * ABSCO2
            FRACS(LAY,NGS6+IG) = FRACREFB(IG)
 3000    CONTINUE

C EMPIRICAL MODIFICATION TO CODE TO IMPROVE STRATOSPHERIC COOLING RATES
C FOR O3.  REVISED TO APPLY WEIGHTING FOR G-POINT REDUCTION IN THIS BAND.

         TAUG(LAY,NGS6+6)=TAUG(LAY,NGS6+6)*0.92
         TAUG(LAY,NGS6+7)=TAUG(LAY,NGS6+7)*0.88
         TAUG(LAY,NGS6+8)=TAUG(LAY,NGS6+8)*1.07
         TAUG(LAY,NGS6+9)=TAUG(LAY,NGS6+9)*1.1
         TAUG(LAY,NGS6+10)=TAUG(LAY,NGS6+10)*0.99
         TAUG(LAY,NGS6+11)=TAUG(LAY,NGS6+11)*0.855

 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB8

C     BAND 8:  1080-1180 cm-1 (low key - H2O; low minor - CO2,O3,N2O)
C                             (high key - O3; high minor - CO2, N2O)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, MAXXSEC=4, NBANDS=16)
      PARAMETER (NGPT=140, NG8=8, NGS7=88)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY),SELFFRAC(MXLAY),INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K8C/      KA(5,13,NG8), KB(5,13:59,NG8), FORREF(4,NG8),
     &              SELFREF(10,NG8), KA_MCO2(19,NG8), KA_MO3(19,NG8),
     &              KA_MN2O(19,NG8), KB_MCO2(19,NG8), KB_MN2O(19,NG8),
     &                 CFC12(NG8), CFC22ADJ(NG8)
      COMMON /PF8C/    FRACREFA(NG8), FRACREFB(NG8)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      REAL KA,KB,KA_MCO2,KA_MO3,KA_MN2O,KB_MCO2,KB_MN2O, MINORFRAC              

      DIMENSION ABSA(65,NG8),ABSB(235,NG8)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C Minor gas mapping level:
C     LOWER - CO2, P = 1053.63 mb, T = 294.2 K
C     LOWER - O3,  P = 317.348 mb, T = 240.77 K
C     LOWER - N2O, P = 706.2720 mb, T= 278.94 K
C     LOWER - CFC12,CFC11
C     UPPER - CO2, P = 35.1632 mb, T = 223.28 K
C     UPPER - N2O, P = 8.716e-2 mb, T = 226.03 K

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature, and appropriate species.  Below LAYTROP, the water vapor 
C     self-continuum and foreign continuum is interpolated (in temperature) 
C     separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.65
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           *COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(8) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(8) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         DO 2000 IG = 1, NG8
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG)))
            ABSCO2 =  (KA_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MCO2(INDM+1,IG) - KA_MCO2(INDM,IG)))
            ABSO3 =  (KA_MO3(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MO3(INDM+1,IG) - KA_MO3(INDM,IG)))
            ABSN2O =  (KA_MN2O(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MN2O(INDM+1,IG) - KA_MN2O(INDM,IG)))
            TAUG(LAY,NGS7+IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR
     &           + ADJCOLCO2*ABSCO2
     &           + COLO3(LAY) * ABSO3
     &           + COLN2O(LAY) * ABSN2O
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)
            FRACS(LAY,NGS7+IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/COLDRY(LAY)
         RATCO2 = 1.E20*CHI_CO2/CHI_MLS(2,JP(LAY)+1)
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.65
            ADJCOLCO2 = ADJFAC*CHI_MLS(2,JP(LAY)+1)
     &           * COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(8) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(8) + 1
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG8
            ABSCO2 =  (KB_MCO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MCO2(INDM+1,IG) - KB_MCO2(INDM,IG)))
            ABSN2O =  (KB_MN2O(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MN2O(INDM+1,IG) - KB_MN2O(INDM,IG)))
            TAUG(LAY,NGS7+IG) = COLO3(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + ADJCOLCO2*ABSCO2
     &           + COLN2O(LAY)*ABSN2O 
     &           + WX(3,LAY) * CFC12(IG)
     &           + WX(4,LAY) * CFC22ADJ(IG)
            FRACS(LAY,NGS7+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB9

C     BAND 9:  1180-1390 cm-1 (low key - H2O,CH4; low minor - N2O)
C                             (high key - CH4; high minor - N2O)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, NBANDS=16)
      PARAMETER (NGPT=140, NG9=12, NGS8=96)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K9C/      KA(9,5,13,NG9), KB(5,13:59,NG9), FORREF(4,NG9),
     &              SELFREF(10,NG9), KA_MN2O(9,19,NG9), KB_MN2O(19,NG9)
      COMMON /PF9C/     FRACREFA(NG9,9), FRACREFB(NG9)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      REAL KA,KB
      REAL KA_MN2O,KB_MN2O,MINORFRAC,N2OM1,N2OM2

      DIMENSION ABSA(585,NG9),ABSB(235,NG9)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

C Minor gas mapping level :
C     LOWER - N2O, P = 706.272 mbar, T = 278.94 K
C     UPPER - N2O, P = 95.58 mbar, T = 215.7 K

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 212 mb
      REFRAT_PLANCK_A = CHI_MLS(1,9)/CHI_MLS(6,9)

C     P = 706.272 mb 
      REFRAT_M_A = CHI_MLS(1,3)/CHI_MLS(6,3)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCH4(LAY)*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCH4_1(LAY)*COLCH4(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2O = COLH2O(LAY) + REFRAT_M_A*COLCH4(LAY)
         SPECPARM_MN2O = COLH2O(LAY)/SPECCOMB_MN2O
         IF (SPECPARM_MN2O .GE. ONEMINUS) SPECPARM_MN2O = ONEMINUS
         SPECMULT_MN2O = 8.*SPECPARM_MN2O
         JMN2O = 1 + INT(SPECMULT_MN2O)
         FMN2O = AMOD(SPECMULT_MN2O,1.0)

c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/(COLDRY(LAY))
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCH4(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(9) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(9) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG9
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               TAUG(LAY,NGS8+IG) = SPECCOMB *
     &          (FAC000 * ABSA(IND0,IG) +
     &          FAC100 * ABSA(IND0+1,IG) +
     &          FAC200 * ABSA(IND0+2,IG) +
     &          FAC010 * ABSA(IND0+9,IG) +
     &          FAC110 * ABSA(IND0+10,IG) +
     &          FAC210 * ABSA(IND0+11,IG))
     &          + SPECCOMB1 *
     &          (FAC001 * ABSA(IND1,IG) + 
     &          FAC101 * ABSA(IND1+1,IG) +
     &          FAC201 * ABSA(IND1+2,IG) +
     &          FAC011 * ABSA(IND1+9,IG) +
     &          FAC111 * ABSA(IND1+10,IG) +
     &          FAC211 * ABSA(IND1+11,IG)) 
     &          + TAUSELF + TAUFOR
     &          + ADJCOLN2O*ABSN2O            
               FRACS(LAY,NGS8+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG9
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               TAUG(LAY,NGS8+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
               FRACS(LAY,NGS8+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG9
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2OM1 = KA_MN2O(JMN2O,INDM,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM,IG)-
     &              KA_MN2O(JMN2O,INDM,IG))
               N2OM2 = KA_MN2O(JMN2O,INDM+1,IG) + FMN2O*
     &              (KA_MN2O(JMN2O+1,INDM+1,IG)-
     &              KA_MN2O(JMN2O,INDM+1,IG))
               ABSN2O = N2OM1 + MINORFRAC(LAY)
     &              * (N2OM2 - N2OM1)
               TAUG(LAY,NGS8+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLN2O*ABSN2O            
               FRACS(LAY,NGS8+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
c     In atmospheres where the amount of N2O is too great to be considered
c     a minor species, adjust the column amount of N2O by an empirical factor 
c     to obtain the proper contribution.
         CHI_N2O = COLN2O(LAY)/(COLDRY(LAY))
         RATN2O = 1.E20*CHI_N2O/CHI_MLS(4,JP(LAY)+1)
         IF (RATN2O .GT. 1.5) THEN
            ADJFAC = 0.5+(RATN2O-0.5)**0.65
            ADJCOLN2O = ADJFAC*CHI_MLS(4,JP(LAY)+1)*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLN2O = COLN2O(LAY)
         ENDIF

         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(9) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(9) + 1
         INDM = INDMINOR(LAY)

         DO 3000 IG = 1, NG9
            ABSN2O = KB_MN2O(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MN2O(INDM+1,IG) - KB_MN2O(INDM,IG))
            TAUG(LAY,NGS8+IG) = COLCH4(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + ADJCOLN2O*ABSN2O
            FRACS(LAY,NGS8+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB10

C     BAND 10:  1390-1480 cm-1 (low key - H2O; high key - H2O)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG10=6, NGS9=108)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K10C/    KA(5,13,NG10), KB(5,13:59,NG10), FORREF(4,NG10),
     &                  SELFREF(10,NG10)
      COMMON /PF10C/    FRACREFA(NG10), FRACREFB(NG10)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(65,NG10),ABSB(235,NG10)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(10) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(10) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         DO 2000 IG = 1, NG10
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUG(LAY,NGS9+IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG)) 
     &           + TAUSELF + TAUFOR
            FRACS(LAY,NGS9+IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(10) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(10) + 1
         INDF = INDFOR(LAY)
         DO 3000 IG = 1, NG10
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUG(LAY,NGS9+IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
     &           + TAUFOR
            FRACS(LAY,NGS9+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB11

C     BAND 11:  1480-1800 cm-1 (low - H2O; low minor - O2)
C                              (high key - H2O; high minor - O2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG11=8, NGS10=114)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K11C/  KA(5,13,NG11), KB(5,13:59,NG11), FORREF(4,NG11),
     &               SELFREF(10,NG11), KA_MO2(19,NG11), KB_MO2(19,NG11)
      COMMON /PF11C/   FRACREFA(NG11), FRACREFB(NG11)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(65,NG11),ABSB(235,NG11)
      REAL KA_MO2, KB_MO2, MINORFRAC

      EQUIVALENCE (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C Minor gas mapping level :
C     LOWER - O2, P = 706.2720 mbar, T = 278.94 K
C     UPPER - O2, P = 4.758820 mbarm T = 250.85 K

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum and
C     foreign continuum is interpolated (in temperature) separately.

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(11) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(11) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         SCALEO2 = COLO2(LAY)*SCALEMINOR(LAY)
         DO 2000 IG = 1, NG11
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG)))
            TAUO2 =  SCALEO2 *
     &           (KA_MO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KA_MO2(INDM+1,IG) - KA_MO2(INDM,IG)))
            TAUG(LAY,NGS10+IG) = COLH2O(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))
     &           + TAUSELF + TAUFOR
     &           + TAUO2
            FRACS(LAY,NGS10+IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(11) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(11) + 1
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         SCALEO2 = COLO2(LAY)*SCALEMINOR(LAY)
         DO 3000 IG = 1, NG11
            TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) *
     &           (FORREF(INDF+1,IG) - FORREF(INDF,IG))) 
            TAUO2 =  SCALEO2*
     &           (KB_MO2(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MO2(INDM+1,IG) - KB_MO2(INDM,IG)))
            TAUG(LAY,NGS10+IG) = COLH2O(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG)) 
     &           + TAUFOR
     &           + TAUO2
            FRACS(LAY,NGS10+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB12

C     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG12=8, NGS11=122)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)  
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)            
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K12C/     KA(9,5,13,NG12),FORREF(4,NG12),SELFREF(10,NG12)
      COMMON /PF12C/    FRACREFA(NG12,9)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG12)

      EQUIVALENCE (KA,ABSA)
      REAL KA

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P =   174.164 mb 
      REFRAT_PLANCK_A = CHI_MLS(1,10)/CHI_MLS(2,10)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum adn foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(12) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(12) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG12
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS11+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS11+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG12
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
               TAUG(LAY,NGS11+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS11+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG12
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS11+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS11+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG12
            TAUG(LAY,NGS11+IG) = 0.0
            FRACS(LAY,NGS11+IG) = 0.0
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB13

C     BAND 13:  2080-2250 cm-1 (low key - H2O,N2O; high minor - O3 minor)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, NBANDS=16)
      PARAMETER (NGPT=140, NG13=4, NGS12=130)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K13C/    KA(9,5,13,NG13),FORREF(4,NG13),SELFREF(10,NG13),
     &                  KA_MCO2(9,19,NG13), KA_MCO(9,19,NG13),
     &                  KB_MO3(19,NG13)
      COMMON /PF13C/    FRACREFA(NG13,9), FRACREFB(NG13)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG13)
      REAL KA_MCO2,KA_MCO,KB_MO3,MINORFRAC

      EQUIVALENCE (KA,ABSA)
      REAL KA

C Minor gas mapping levels :
C     LOWER - CO2, P = 1053.63 mb, T = 294.2 K
C     LOWER - CO, P = 706 mb, T = 278.94 K
C     UPPER - O3, P = 95.5835 mb, T = 215.7 K

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower/upper atmosphere.

C     P = 473.420 mb (Level 5)
      REFRAT_PLANCK_A = CHI_MLS(1,5)/CHI_MLS(4,5)

C     P = 1053. (Level 1)
      REFRAT_M_A = CHI_MLS(1,1)/CHI_MLS(4,1)

C     P = 706. (Level 3)
      REFRAT_M_A3 = CHI_MLS(1,3)/CHI_MLS(4,3)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2ON2O(LAY)*COLN2O(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2ON2O_1(LAY)*COLN2O(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MCO2 = COLH2O(LAY) + REFRAT_M_A*COLN2O(LAY)
         SPECPARM_MCO2 = COLH2O(LAY)/SPECCOMB_MCO2
         IF (SPECPARM_MCO2 .GE. ONEMINUS) SPECPARM_MCO2 = ONEMINUS
         SPECMULT_MCO2 = 8.*SPECPARM_MCO2
         JMCO2 = 1 + INT(SPECMULT_MCO2)
         FMCO2 = AMOD(SPECMULT_MCO2,1.0)

c     In atmospheres where the amount of CO2 is too great to be considered
c     a minor species, adjust the column amount of CO2 by an empirical factor 
c     to obtain the proper contribution.
         CHI_CO2 = COLCO2(LAY)/(COLDRY(LAY))
         RATCO2 = 1.E20*CHI_CO2/3.55E-4
         IF (RATCO2 .GT. 3.0) THEN
            ADJFAC = 2.0+(RATCO2-2.0)**0.68
            ADJCOLCO2 = ADJFAC*3.55E-4*COLDRY(LAY)*1.E-20
         ELSE
            ADJCOLCO2 = COLCO2(LAY)
         ENDIF

         SPECCOMB_MCO = COLH2O(LAY) + REFRAT_M_A3*COLN2O(LAY)
         SPECPARM_MCO = COLH2O(LAY)/SPECCOMB_MCO
         IF (SPECPARM_MCO .GE. ONEMINUS) SPECPARM_MCO = ONEMINUS
         SPECMULT_MCO = 8.*SPECPARM_MCO
         JMCO = 1 + INT(SPECMULT_MCO)
         FMCO = AMOD(SPECMULT_MCO,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLN2O(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(13) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(13) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG13
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + 
     &              MINORFRAC(LAY) * (CO2M2 - CO2M1)
               COM1 = KA_MCO(JMCO,INDM,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM,IG)-
     &              KA_MCO(JMCO,INDM,IG))
               COM2 = KA_MCO(JMCO,INDM+1,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM+1,IG)-
     &              KA_MCO(JMCO,INDM+1,IG))
               ABSCO = COM1 + 
     &              MINORFRAC(LAY) * (COM2 - COM1)
               TAUG(LAY,NGS12+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))+
     &              SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
     &              + COLCO(LAY)*ABSCO
               FRACS(LAY,NGS12+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG13
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 +
     &              MINORFRAC(LAY) * (CO2M2 - CO2M1)
               COM1 = KA_MCO(JMCO,INDM,IG) + FMCO*
     &           (KA_MCO(JMCO+1,INDM,IG)-
     &              KA_MCO(JMCO,INDM,IG))
               COM2 = KA_MCO(JMCO,INDM+1,IG) + FMCO*
     &           (KA_MCO(JMCO+1,INDM+1,IG)-
     &              KA_MCO(JMCO,INDM+1,IG))
               ABSCO = COM1 +
     &              MINORFRAC(LAY) * (COM2 - COM1)
               TAUG(LAY,NGS12+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
     &              + COLCO(LAY)*ABSCO
                FRACS(LAY,NGS12+IG) = FRACREFA(IG,JPL) + FPL *
     &               (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG13
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               CO2M1 = KA_MCO2(JMCO2,INDM,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM,IG)-
     &              KA_MCO2(JMCO2,INDM,IG))
               CO2M2 = KA_MCO2(JMCO2,INDM+1,IG) + FMCO2*
     &              (KA_MCO2(JMCO2+1,INDM+1,IG)-
     &              KA_MCO2(JMCO2,INDM+1,IG))
               ABSCO2 = CO2M1 + 
     &              MINORFRAC(LAY) * (CO2M2 - CO2M1)
               COM1 = KA_MCO(JMCO,INDM,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM,IG)-
     &              KA_MCO(JMCO,INDM,IG))
               COM2 = KA_MCO(JMCO,INDM+1,IG) + FMCO*
     &              (KA_MCO(JMCO+1,INDM+1,IG)-
     &              KA_MCO(JMCO,INDM+1,IG))
               ABSCO = COM1 +
     &              MINORFRAC(LAY) * (COM2 - COM1)
               TAUG(LAY,NGS12+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + ADJCOLCO2*ABSCO2
     &              + COLCO(LAY)*ABSCO
               FRACS(LAY,NGS12+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))

 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         INDM = INDMINOR(LAY)
         DO 3000 IG = 1, NG13
            ABSO3 = KB_MO3(INDM,IG) + 
     &           MINORFRAC(LAY) *
     &           (KB_MO3(INDM+1,IG) - KB_MO3(INDM,IG))
            TAUG(LAY,NGS12+IG) = COLO3(LAY)*ABSO3
            FRACS(LAY,NGS12+IG) =  FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB14

C     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG14=2, NGS13=134)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /K14C/    KA(5,13,NG14), KB(5,13:59,NG14), FORREF(4,NG14),
     &                  SELFREF(10,NG14)
      COMMON /PF14C/    FRACREFA(NG14), FRACREFB(NG14)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(65,NG14),ABSB(235,NG14)

      Equivalence (KA,ABSA),(KB,ABSB)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     and foreign continuum is interpolated (in temperature) separately.  

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(14) + 1
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(14) + 1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         DO 2000 IG = 1, NG14
            TAUSELF = SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
            TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &           FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &           FORREF(INDF,IG))) 
            TAUG(LAY,NGS13+IG) = COLCO2(LAY) *
     &          (FAC00(LAY) * ABSA(IND0,IG) +
     &           FAC10(LAY) * ABSA(IND0+1,IG) +
     &           FAC01(LAY) * ABSA(IND1,IG) + 
     &           FAC11(LAY) * ABSA(IND1+1,IG))
     &           + TAUSELF + TAUFOR
            FRACS(LAY,NGS13+IG) = FRACREFA(IG)
 2000    CONTINUE
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(14) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(14) + 1
         DO 3000 IG = 1, NG14
            TAUG(LAY,NGS13+IG) = COLCO2(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
            FRACS(LAY,NGS13+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB15

C     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; low minor - N2)
C                              (high - nothing)

      PARAMETER (MG=16, MXLAY=203, MXMOL=38, NBANDS=16)
      PARAMETER (NGPT=140, NG15=2, NGS14=136)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /MINOR/    MINORFRAC(MXLAY), INDMINOR(MXLAY), 
     &                  SCALEMINOR(MXLAY),SCALEMINORN2(MXLAY)
      COMMON /K15C/     KA(9,5,13,NG15), 
     &              FORREF(4,NG15), SELFREF(10,NG15), KA_MN2(9,19,NG15)
      COMMON /PF15C/    FRACREFA(NG15,9)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG15)

      EQUIVALENCE (KA,ABSA)
      REAL KA, KA_MN2, MINORFRAC
      REAL N2M1,N2M2

C Minor gas mapping level : 
C     Lower - Nitrogen Continuum, P = 1053., T = 294.

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.
C     P = 1053. mb (Level 1)
      REFRAT_PLANCK_A = CHI_MLS(4,1)/CHI_MLS(2,1)

C     P = 1053.
      REFRAT_M_A = CHI_MLS(4,1)/CHI_MLS(2,1)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLN2O(LAY) + RAT_N2OCO2(LAY)*COLCO2(LAY)
         SPECPARM = COLN2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLN2O(LAY) + RAT_N2OCO2_1(LAY)*COLCO2(LAY)
         SPECPARM1 = COLN2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_MN2 = COLN2O(LAY) + REFRAT_M_A*COLCO2(LAY)
         SPECPARM_MN2 = COLN2O(LAY)/SPECCOMB_MN2
         IF (SPECPARM_MN2 .GE. ONEMINUS) SPECPARM_MN2 = ONEMINUS
         SPECMULT_MN2 = 8.*SPECPARM_MN2
         JMN2 = 1 + INT(SPECMULT_MN2)
         FMN2 = AMOD(SPECMULT_MN2,1.0)

         SPECCOMB_PLANCK = COLN2O(LAY)+REFRAT_PLANCK_A*COLCO2(LAY)
         SPECPARM_PLANCK = COLN2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(15) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(15) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)
         INDM = INDMINOR(LAY)
         
         SCALEN2 = COLBRD(LAY)*SCALEMINOR(LAY)
         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG15
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG)-
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG)-
     &              KA_MN2(JMN2,INDM+1,IG))
               TAUN2 = SCALEN2*
     &              (N2M1 + MINORFRAC(LAY) * (N2M2 - N2M1))
               TAUG(LAY,NGS14+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
     &              + TAUN2
               FRACS(LAY,NGS14+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2010 IG = 1, NG15
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG)-
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG)-
     &              KA_MN2(JMN2,INDM+1,IG))
               TAUN2 = SCALEN2*
     &              (N2M1 + MINORFRAC(LAY) * (N2M2 - N2M1))
               TAUG(LAY,NGS14+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
     &              + TAUN2
               FRACS(LAY,NGS14+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)

            DO 2020 IG = 1, NG15
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               N2M1 = KA_MN2(JMN2,INDM,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM,IG)-
     &              KA_MN2(JMN2,INDM,IG))
               N2M2 = KA_MN2(JMN2,INDM+1,IG) + FMN2*
     &              (KA_MN2(JMN2+1,INDM+1,IG)-
     &              KA_MN2(JMN2,INDM+1,IG))
               TAUN2 = SCALEN2 *
     &              (N2M1 + MINORFRAC(LAY) * (N2M2 - N2M1))
               TAUG(LAY,NGS14+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
     &              + TAUN2
               FRACS(LAY,NGS14+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF
 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG15
            TAUG(LAY,NGS14+IG) = 0.0
            FRACS(LAY,NGS14+IG) = 0.0
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB16

C     BAND 16:  2600-3250 cm-1 (low key- H2O,CH4; high key - CH4)

      PARAMETER (MG=16, MXLAY=203, NBANDS=16)
      PARAMETER (NGPT=140, NG16=2, NGS15=138)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,NGPT)
      COMMON /PLANKG/   FRACS(MXLAY,NGPT)                                       

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PRECISE/  ONEMINUS
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,                                   
     &                  COLH2O(MXLAY),COLCO2(MXLAY),COLO3(MXLAY),  
     &                  COLN2O(MXLAY),COLCO(MXLAY),COLCH4(MXLAY),  
     &                  COLO2(MXLAY),COLBRD(MXLAY)
      COMMON /MLS_REF/  PREF(59),PREFLOG(59),TREF(59),CHI_MLS(7,59)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            
     &                  FAC10(MXLAY),FAC11(MXLAY)                             
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /REFRAT_ETA/ RAT_H2OCO2(MXLAY),RAT_H2OCO2_1(MXLAY),
     &                  RAT_H2OO3(MXLAY),RAT_H2OO3_1(MXLAY),
     &                  RAT_H2ON2O(MXLAY),RAT_H2ON2O_1(MXLAY),
     &                  RAT_H2OCH4(MXLAY),RAT_H2OCH4_1(MXLAY),
     &                  RAT_N2OCO2(MXLAY),RAT_N2OCO2_1(MXLAY),
     &                  RAT_O3CO2(MXLAY),RAT_O3CO2_1(MXLAY)
      COMMON /FOREIGN/  FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K16C/     KA(9,5,13,NG16), KB(5,13:59,NG16),
     &                  FORREF(4,NG16), SELFREF(10,NG16)
      COMMON /PF16C/    FRACREFA(NG16,9), FRACREFB(NG16)

      COMMON /CVRTAU/    HNAMTAU,HVRTAU

      CHARACTER*18       HNAMTAU,HVRTAU

      DIMENSION ABSA(585,NG16), ABSB(235,NG16)

      EQUIVALENCE (KA,ABSA), (KB,ABSB)
      REAL KA,KB

C     Calculate reference ratio to be used in calculation of Planck
C     fraction in lower atmosphere.

C     P = 387. mb (Level 6)
      REFRAT_PLANCK_A = CHI_MLS(1,6)/CHI_MLS(6,6)

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature,and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum and foreign continuum is interpolated 
C     (in temperature) separately.  

      HVRTAU = '$Revision$'

      DO 2500 LAY = 1, LAYTROP

         SPECCOMB = COLH2O(LAY) + RAT_H2OCH4(LAY)*COLCH4(LAY)
         SPECPARM = COLH2O(LAY)/SPECCOMB
         IF (SPECPARM .GE. ONEMINUS) SPECPARM = ONEMINUS
         SPECMULT = 8.*(SPECPARM)
         JS = 1 + INT(SPECMULT)
         FS = AMOD(SPECMULT,1.0)

         SPECCOMB1 = COLH2O(LAY) + RAT_H2OCH4_1(LAY)*COLCH4(LAY)
         SPECPARM1 = COLH2O(LAY)/SPECCOMB1
         IF (SPECPARM1 .GE. ONEMINUS) SPECPARM1 = ONEMINUS
         SPECMULT1 = 8.*(SPECPARM1)
         JS1 = 1 + INT(SPECMULT1)
         FS1 = AMOD(SPECMULT1,1.0)

         SPECCOMB_PLANCK = COLH2O(LAY)+REFRAT_PLANCK_A*COLCH4(LAY)
         SPECPARM_PLANCK = COLH2O(LAY)/SPECCOMB_PLANCK
         IF (SPECPARM_PLANCK .GE. ONEMINUS) SPECPARM_PLANCK=ONEMINUS
         SPECMULT_PLANCK = 8.*SPECPARM_PLANCK
         JPL= 1 + INT(SPECMULT_PLANCK)
         FPL = AMOD(SPECMULT_PLANCK,1.0)

         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(16) + JS
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(16) + JS1
         INDS = INDSELF(LAY)
         INDF = INDFOR(LAY)

         IF (SPECPARM .LT. 0.125 .AND. SPECPARM1 .LT. 0.125) THEN
            P = FS - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = FS1 - 1
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)

            DO 2000 IG = 1, NG16
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR =  FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS15+IG) = SPECCOMB *
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC200 * ABSA(IND0+2,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG) +
     &              FAC210 * ABSA(IND0+11,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC201 * ABSA(IND1+2,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG) +
     &              FAC211 * ABSA(IND1+11,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS15+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2000          CONTINUE
      ELSE IF (SPECPARM .GT. 0.875 .AND. SPECPARM1 .GT. 0.875) THEN
            P = -FS 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC000 = FK0*FAC00(LAY)
            FAC100 = FK1*FAC00(LAY)
            FAC200 = FK2*FAC00(LAY)
            FAC010 = FK0*FAC10(LAY)
            FAC110 = FK1*FAC10(LAY)
            FAC210 = FK2*FAC10(LAY)

            P = -FS1 
            P4 = P**4
            FK0 = P4
            FK1 = 1 - P - 2.0*P4
            FK2 = P + P4
            FAC001 = FK0*FAC01(LAY)
            FAC101 = FK1*FAC01(LAY)
            FAC201 = FK2*FAC01(LAY)
            FAC011 = FK0*FAC11(LAY)
            FAC111 = FK1*FAC11(LAY)
            FAC211 = FK2*FAC11(LAY)
            DO 2010 IG = 1, NG16
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS15+IG) = SPECCOMB * 
     &              (FAC200 * ABSA(IND0-1,IG) +
     &              FAC100 * ABSA(IND0,IG) +
     &              FAC000 * ABSA(IND0+1,IG) +
     &              FAC210 * ABSA(IND0+8,IG) +
     &              FAC110 * ABSA(IND0+9,IG) +
     &              FAC010 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC201 * ABSA(IND1-1,IG) +
     &              FAC101 * ABSA(IND1,IG) +
     &              FAC001 * ABSA(IND1+1,IG) +
     &              FAC211 * ABSA(IND1+8,IG) +
     &              FAC111 * ABSA(IND1+9,IG) +
     &              FAC011 * ABSA(IND1+10,IG))
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS15+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2010           CONTINUE
       ELSE
            FAC000 = (1. - FS) * FAC00(LAY)
            FAC010 = (1. - FS) * FAC10(LAY)
            FAC100 = FS * FAC00(LAY)
            FAC110 = FS * FAC10(LAY)

            FAC001 = (1. - FS1) * FAC01(LAY)
            FAC011 = (1. - FS1) * FAC11(LAY)
            FAC101 = FS1 * FAC01(LAY)
            FAC111 = FS1 * FAC11(LAY)
            DO 2020 IG = 1, NG16
               TAUSELF = SELFFAC(LAY)* (SELFREF(INDS,IG) +
     &              SELFFRAC(LAY)  *
     &              (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
               TAUFOR = FORFAC(LAY) * (FORREF(INDF,IG) +
     &              FORFRAC(LAY) * (FORREF(INDF+1,IG) - 
     &              FORREF(INDF,IG))) 
               TAUG(LAY,NGS15+IG) = SPECCOMB * 
     &              (FAC000 * ABSA(IND0,IG) +
     &              FAC100 * ABSA(IND0+1,IG) +
     &              FAC010 * ABSA(IND0+9,IG) +
     &              FAC110 * ABSA(IND0+10,IG))
     &              + SPECCOMB1 *
     &              (FAC001 * ABSA(IND1,IG) + 
     &              FAC101 * ABSA(IND1+1,IG) +
     &              FAC011 * ABSA(IND1+9,IG) +
     &              FAC111 * ABSA(IND1+10,IG)) 
     &              + TAUSELF + TAUFOR
               FRACS(LAY,NGS15+IG) = FRACREFA(IG,JPL) + FPL *
     &              (FRACREFA(IG,JPL+1)-FRACREFA(IG,JPL))
 2020          CONTINUE
      ENDIF

 2500 CONTINUE

      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(16) + 1
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(16) + 1
         DO 3000 IG = 1, NG16
            TAUG(LAY,NGS15+IG) = COLCH4(LAY) * 
     &          (FAC00(LAY) * ABSB(IND0,IG) +
     &           FAC10(LAY) * ABSB(IND0+1,IG) +
     &           FAC01(LAY) * ABSB(IND1,IG) + 
     &           FAC11(LAY) * ABSB(IND1+1,IG))
            FRACS(LAY,NGS15+IG) = FRACREFB(IG)
 3000    CONTINUE
 3500 CONTINUE

      RETURN
      END
