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
C *** g-points for application of the model to GCMs.  This stand-alone 
C *** version will run with the standard RRTM input files.
C *** Michael J. Iacono, AER, Inc.; September, 2002

       PROGRAM RRTMG_LW
                    
C *** This program is the driver for RRTM, the AER rapid model.  
C     For each atmosphere the user wishes to analyze, this routine
C     a) calls READPROF to read in the atmospheric profile
C     b) calls SETCOEF to calculate various quantities needed for 
C        the radiative transfer algorithm
C     c) calls RTRNMR to do the radiative transfer for clear sky or
C        for cloudy sky using a maximum-random cloud overlap method,
C        OR calls RTRN to use a random cloud overlap method
C     d) writes out the upward, downward, and net flux for each
C        level and the heating rate for each layer

      PARAMETER (MXLAY=203)
      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)

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
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /OUTCLR/    TOTUCLFL(0:MXLAY), TOTDCLFL(0:MXLAY),
     &                   FNETC(0:MXLAY), HTRC(0:MXLAY)

      COMMON /CVRRTM/    HNAMRTM,HVRRTM
      COMMON /CVRRTR/    HNAMRTR,HVRRTR
      COMMON /CVRATM/    HNAMATM,HVRATM
      COMMON /CVRSET/    HNAMSET,HVRSET
      COMMON /CVRTAU/    HNAMTAU,HVRTAU
      COMMON /CVRINI/    HNAMINI,HVRINI
      COMMON /CVRCLD/    HNAMCLD,HVRCLD
      COMMON /CVRUTL/    HNAMUTL,HVRUTL
      COMMON /CVREXT/    HNAMEXT,HVREXT
      COMMON /CVRRTX/    HNAMRTX,HVRRTX

      COMMON /CVRSN/     HNAMKG,HVRKG

      CHARACTER*18 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *             HVRINI,HVRCLD,HVRUTL,HVREXT,HVRRTX

      CHARACTER*18 HNAMRTM,HNAMRTR,HNAMATM,HNAMSET,HNAMTAU,
     *             HNAMINI,HNAMCLD,HNAMUTL,HNAMEXT,HNAMRTX

      CHARACTER*18 HVRKG

      CHARACTER*18 HNAMKG

      CHARACTER PAGE

C      DATA WAVENUM1(1) /10./, WAVENUM2(1) /350./, DELWAVE(1) /340./
C      DATA WAVENUM1(2) /350./, WAVENUM2(2) /500./, DELWAVE(2) /150./
C      DATA WAVENUM1(3) /500./, WAVENUM2(3) /630./, DELWAVE(3) /130./
C      DATA WAVENUM1(4) /630./, WAVENUM2(4) /700./, DELWAVE(4) /70./
C      DATA WAVENUM1(5) /700./, WAVENUM2(5) /820./, DELWAVE(5) /120./
C      DATA WAVENUM1(6) /820./, WAVENUM2(6) /980./, DELWAVE(6) /160./
C      DATA WAVENUM1(7) /980./, WAVENUM2(7) /1080./, DELWAVE(7) /100./
C      DATA WAVENUM1(8) /1080./, WAVENUM2(8) /1180./, DELWAVE(8) /100./
C      DATA WAVENUM1(9) /1180./, WAVENUM2(9) /1390./, DELWAVE(9) /210./
C      DATA WAVENUM1(10) /1390./,WAVENUM2(10) /1480./,DELWAVE(10) /90./
C      DATA WAVENUM1(11) /1480./,WAVENUM2(11) /1800./,DELWAVE(11) /320./
C      DATA WAVENUM1(12) /1800./,WAVENUM2(12) /2080./,DELWAVE(12) /280./
C      DATA WAVENUM1(13) /2080./,WAVENUM2(13) /2250./,DELWAVE(13) /170./
C      DATA WAVENUM1(14) /2250./,WAVENUM2(14) /2380./,DELWAVE(14) /130./
C      DATA WAVENUM1(15) /2380./,WAVENUM2(15) /2600./,DELWAVE(15) /220./
C      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3250./,DELWAVE(16) /650./

C      DATA NG  /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
C      DATA NSPA /1, 1,10, 9, 9, 1, 9, 1,11, 1, 1, 9, 9, 1, 9, 9/
C      DATA NSPB /1, 1, 5, 6, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
C      DATA HEATFAC /8.4391/

      ONEMINUS = 1. - 1.E-6
      FLUXFAC = PI * 2.D4  

      IWR = 10
      PAGE = CHAR(12)

      HVRRTM = '$Revision$'

C  Run model initialization, compute lookup tables, perform
C  reduction of g-points from 256 to 140 for input absorption
C  coefficient data and other arrays.

      CALL RRTMINIT

C  Open the INPUT set of atmospheres
      IRD = 9
      OPEN (IRD,FILE='INPUT_RRTM',FORM='FORMATTED')
      
c  Multiple atmosphere option implemented      
      NUMATMOS = 1
      DO 4000 IATMOS = 1, NUMATMOS

C ***    Input atmospheric profile from INPUT_RRTM.
         CALL READPROF

         ISTART = 1
         IEND = 16
         IFLAG = IOUT
         NUMANGS = 0

 1000    CONTINUE
         IF (IFLAG .GT. 0 .AND. IFLAG .LE. 40) THEN
            ISTART = IFLAG
            IEND = IFLAG
         ENDIF

C ***    Calculate information needed by the radiative transfer routine
C        that is specific to this atmosphere, especially some of the 
C        coefficients and indices needed to compute the optical depths
C        by interpolating data from stored reference atmospheres. 
         ICLDATM = 0
         IF (ICLD .GE. 1) CALL CLDPROP(ICLDATM)

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

C ***    Call the radiative transfer routine.
          IF (ICLD .EQ. 1) THEN
             CALL RTRN 
          ELSE
             CALL RTRNMR
          ENDIF

         IF (IOUT .LT. 0) GO TO 4000

C ***    Process output for this atmosphere.
         OPEN (IWR,FILE='OUTPUT_RRTM',FORM='FORMATTED')
         WRITE(IWR,9899)WAVENUM1(ISTART),WAVENUM2(IEND),IATMOS
         WRITE(IWR,9900)
         WRITE(IWR,9901)
C
         DO 3000 I = NLAYERS, 0, -1
            IF (PZ(I) .LT. 1.E-2) THEN
               WRITE(IWR,9952) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 1.E-1) THEN
               WRITE(IWR,9953) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 1.) THEN
               WRITE(IWR,9954) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 10.) THEN
               WRITE(IWR,9955) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 100.) THEN
               WRITE(IWR,9956) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSEIF (PZ(I) .LT. 1000.) THEN
               WRITE(IWR,9957) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ELSE
               WRITE(IWR,9958) I,PZ(I),TOTUFLUX(I),TOTDFLUX(I),
     &              FNET(I), HTR(I)
            ENDIF
 3000    CONTINUE
         WRITE(IWR,9903)PAGE
         
         IF (IOUT .LE. 40 .OR. IFLAG .EQ. 16) GO TO 3500
         IF (IFLAG .EQ. 99) THEN
            IFLAG = 1
         ELSEIF (IOUT .EQ. 99) THEN
            IFLAG = IFLAG + 1
         ENDIF
         GO TO 1000

 3500    CONTINUE
C
C ***    Output module version numbers
C
         WRITE(IWR,9910) HNAMRTM,HVRRTM,HNAMINI,HVRINI,HNAMCLD,HVRCLD,
     *     HNAMRTR,HVRRTR,HNAMRTX,HVRRTX,HNAMSET,HVRSET,HNAMTAU,HVRTAU,
     *     HNAMATM,HVRATM,HNAMUTL,HVRUTL,HNAMEXT,HVREXT,HNAMKG,HVRKG

 4000 CONTINUE
 4444 CONTINUE

      CLOSE(IRD)
      CLOSE(IWR)

 9952 FORMAT(1X,I3,9X,F7.6,3X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9953 FORMAT(1X,I3,9X,F6.5,4X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9954 FORMAT(1X,I3,8X,F6.4,5X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9955 FORMAT(1X,I3,7X,F6.3,6X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9956 FORMAT(1X,I3,6X,F6.2,7X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9957 FORMAT(1X,I3,5X,F6.1,8X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9958 FORMAT(1X,I3,5X,F6.1,8X,F8.4,6X,F8.4,6X,F12.7,10X,F9.5)
 9899 FORMAT(1X,'Wavenumbers: ',F6.1,' - ',F6.1,' cm-1, ATM ',i6)
 9900 FORMAT(1X,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET
     &FLUX       HEATING RATE')
 9901 FORMAT(1X,'            mb          W/m2          W/m2           W/
     &m2          degree/day')
 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,2(G12.6,2X),G13.6,3X,G16.9,0P)
 9903 FORMAT(A)
 9910 FORMAT('  Modules and versions used in this calculation:',/,/,
     *         6(5X,a18,2X,a18,10x,a18,2x,a18,/))

      STOP
      END

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF                                                     
                                                                         
C     Read in atmospheric profile.

      IMPLICIT DOUBLE PRECISION (V)                                      
                                                                         
      PARAMETER (MXLAY=203, MXMOL=38)
      PARAMETER (NBANDS = 16)
      PARAMETER (MAXINPX=MXMOL)
      PARAMETER (MAXXSEC=4)
      PARAMETER (MAXPROD = MXLAY*MAXXSEC)

      DIMENSION ALTZ(0:MXLAY),IXTRANS(14),SEMIS(NBANDS)

      COMMON /CONTROL/  NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SURFACE/  TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBRODL(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /IFIL/     IRD,IPR,IPU,IDUM(15)
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /PATHX/    IXMAX,NXMOL0,IXINDX0(MAXINPX),WX0(MAXINPX,MXLAY)    
      COMMON /XRRTATM/  IXSECT
      COMMON /PWV/      PWVCM

      CHARACTER*80 FORM1(0:1),FORM2(0:1),FORM3(0:1)
      CHARACTER*1 CTEST, CDOLLAR, CPRCNT,CDUM

      DATA CDOLLAR /'$'/
      DATA CPRCNT /'%'/
      DATA IXTRANS /0,0,0,1,2,3,0,0,0,0,0,4,0,0/
c      DATA WX /MAXPROD*0.0/

      FORM1(0) = '(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
      FORM2(0) = '(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))'
      FORM3(0) = '(8E10.3)'
      FORM1(1) = '(G15.7,G10.4,G10.4,A3,I2,1X,2(G7.2,G8.3,G7.2))'
      FORM2(1) = '(G15.7,G10.4,G10.4,A3,I2,23X,(G7.2,G8.3,G7.2))'
      FORM3(1) = '(8G15.7)'

      IXMAX = MAXINPX

 1000 CONTINUE
      READ (IRD,9010,END=8800) CTEST
      IF (CTEST .EQ. CPRCNT) GO TO 8900 
      IF (CTEST .NE. CDOLLAR) GO TO 1000

      READ (IRD,9011) IATM, IXSECT, NUMANGS, IOUT, ICLD
c     If numangs set to -1, reset to default rt code for
c     backwards compatibility with original rrtm
      IF (NUMANGS .EQ. -1) NUMANGS = 0

C     If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      IF (ICLD .GE. 1) CALL READCLD

C     Read in surface information.
      READ (IRD,9012) TBOUND,IEMISS,IREFLECT,(SEMIS(I),I=1,16)
      DO 1500 IBAND = 1, NBANDS
         SEMISS(IBAND) = 1.0
         IF (IEMISS .EQ. 1 .AND. SEMIS(1) .NE. 0.) THEN
            SEMISS(IBAND) = SEMIS(1)
         ELSEIF (IEMISS .EQ. 2) THEN
            IF (SEMIS(IBAND) .NE. 0.) THEN
               SEMISS(IBAND) = SEMIS(IBAND)
            ENDIF
         ENDIF
 1500 CONTINUE

      IF (IATM .EQ. 0) THEN
         READ (IRD,9013) IFORM,NLAYERS,NMOL
         IF (NMOL.EQ.0) NMOL = 7                                    
         READ (IRD,FORM1(IFORM)) PAVEL(1),TAVEL(1),SECNTK,CINP,
     &        IPTHAK,ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
         READ (IRD,FORM3(IFORM)) (WKL(M,1),M=1,7), WBRODL(1)
         IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,1),M=8,NMOL)

         DO 2000 L = 2, NLAYERS
            READ (IRD,FORM2(IFORM)) PAVEL(L),TAVEL(L),SECNTK,CINP,
     &           IPTHRK,ALTZ(L),PZ(L),TZ(L)
            READ (IRD,FORM3(IFORM)) (WKL(M,L),M=1,7), WBRODL(L)
            IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,L),M=8,NMOL)
 2000    CONTINUE                                                            
           
         IF (IXSECT .EQ. 1) THEN                                 
            READ (IRD,9300) NXMOL0
            NXMOL = NXMOL0
            CALL XSIDENT(IRD)
            READ (IRD,9301) IFORMX
C     
            DO 3000 L = 1, NLAYERS       
               READ (IRD,9010) CDUM
               READ (IRD, FORM3(IFORMX)) (WX0(M,L),M=1,7),WBRODX    
               IF (NXMOL0 .GT. 7) READ (IRD,FORM3(IFORMX)) 
     &              (WX0(M,L),M=8,NXMOL0)
 3000       CONTINUE
         ENDIF
      ELSE
         IPU = 7
         IPR = 66
         OPEN(UNIT=IPR,FILE='TAPE6',STATUS='UNKNOWN')
         CALL RRTATM
         IF (IXSECT .EQ. 1) THEN
            DO 3300 MX = 1, NXMOL0
               IXINDX(MX) = IXTRANS(IXINDX0(MX))
 3300       CONTINUE
         ENDIF
      ENDIF
      IF (TBOUND .LT. 0) TBOUND = TZ(0)

C     Test for mixing ratio input.
      IMIX = 1
      DO 3500 M = 1, NMOL
         IF (WKL(M,1) .GT. 1.0) THEN
            IMIX = 0
            GO TO 3600
         ENDIF
 3500 CONTINUE
 3600 CONTINUE

      IF (IXSECT .EQ. 1) THEN
         IMIXX = 0
         IF (WX0(1,1) .LE. 1.0) IMIXX = 1
      ENDIF
      DO 5000 L = 1, NLAYERS
         SUMMOL = 0.0
         DO 4100 IMOL = 2, NMOL
            SUMMOL = SUMMOL + WKL(IMOL,L)
 4100    CONTINUE
         IF (IMIX .EQ. 1) THEN
            COLDRY(L) = WBRODL(L) / (1. - SUMMOL)
            DO 4200 IMOL = 1, NMOL
               WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
 4200       CONTINUE
         ELSE
            COLDRY(L) = WBRODL(L) + SUMMOL
         ENDIF
         AMTTL = AMTTL + COLDRY(L)+WKL(1,L)
         WVTTL = WVTTL + WKL(1,L)
         IF (IXSECT .EQ. 1) THEN
            DO 4400 IX = 1, NXMOL0
               IF (IXINDX(IX) .NE. 0) THEN
                  IF (IMIXX .EQ. 1) THEN
                     WX(IXINDX(IX),L) = COLDRY(L) * WX0(IX,L) * 1.E-20
                  ELSE
                     WX(IXINDX(IX),L) = WX0(IX,L) * 1.E-20
                  ENDIF
               ENDIF
 4400       CONTINUE
         ENDIF
 5000 CONTINUE
 
      WVSH = (18.0154*WVTTL)/(28.9644*AMTTL)
      PWVCM = WVSH*(1.E3*PZ(0))/(1.E2*9.8066)
 
      GO TO 9000

 8800 CONTINUE
 8900 IF (CTEST.EQ.'%') STOP 'END OF INPUT FILE'
 9000 CONTINUE

 9010 FORMAT (A1)
 9011 FORMAT (49X,I1,19X,I1,13X,I2,2X,I3,4X,I1)
 9012 FORMAT (E10.3,1X,I1,2X,I1,16E5.3)
 9013 FORMAT (1X,I1,I3,I5)                                     
 9300 FORMAT (I5)
 9301 FORMAT (1X,I1)

      RETURN
      END 

C************************  SUBROUTINE READCLD  *****************************C

      SUBROUTINE READCLD

C     Purpose:  To read in IN_CLD_RRTM, the file that contains input 
C               cloud properties.

      PARAMETER (MXLAY=203)
      PARAMETER (NBANDS = 16)

      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)

      CHARACTER*1 CTEST, CPERCENT

      DATA CPERCENT /'%'/
      IRDCLD = 11

      OPEN(IRDCLD,FILE='IN_CLD_RRTM',FORM='FORMATTED')

C     Read in cloud input option.  
      READ(IRDCLD,9050) INFLAG, ICEFLAG, LIQFLAG
      DO 500 LAY = 1, NLAYERS
         CLDFRAC(LAY) = 0.
 500  CONTINUE

 1000 CONTINUE
C     For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
C     DAT1 are pertinent.  If CTEST = '%', then there are no more 
C     cloudy layers to process.
      READ (IRDCLD,9100,END=9000) CTEST,LAY,FRAC,DAT1,DAT2,DAT3,DAT4
      IF (CTEST .EQ. CPERCENT) GO TO 9000
      CLDFRAC(LAY) = FRAC
      CLDDAT1(LAY) = DAT1
      CLDDAT2(LAY) = DAT2
      CLDDAT3(LAY) = DAT3
      CLDDAT4(LAY) = DAT4
      GO TO 1000

 9000 CONTINUE
      CLOSE(IRDCLD)

 9050 FORMAT (3X,I2,4X,I1,4X,I1)
 9100 FORMAT (A1,1X,I3,5E10.5)

      RETURN
      END

C************************  SUBROUTINE XSIDENT  *****************************C

      SUBROUTINE XSIDENT(IRD)
C                                                                         
C     This subroutine identifies which cross-sections are to be used.

      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)

      IMPLICIT DOUBLE PRECISION (V)                                     ! 
C                                                                         
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
C                                                                         
C     NXMOL     - number of cross-sections input by user
C     IXINDX(I) - index of cross-section molecule corresponding to Ith
C                 cross-section specified by user
C                 = 0 -- not allowed in RRTM
C                 = 1 -- CCL4
C                 = 2 -- CFC11
C                 = 3 -- CFC12
C                 = 4 -- CFC22
C                                                                         
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          
C                                                                         
      CHARACTER*10 XSNAME(MAXINPX),ALIAS(MAXXSEC,4),BLANK               
C                                                                         
      DATA (ALIAS(1,I),I=1,4)/                                           
     *    'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/ 
      DATA (ALIAS(2,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/         
      DATA (ALIAS(3,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/         
      DATA (ALIAS(4,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/        

      DATA BLANK / '          '/                                          
C                                                                         
      DO 10 I = 1, NXMOL                                                 
         XSNAME(I) = BLANK                                                
   10 CONTINUE                                                            
C                                                                         
C     READ IN THE NAMES OF THE MOLECULES                                  
C                                                                         
      IF (NXMOL.GT.7) THEN                                               
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7)                            
         READ (IRD,'(8A10)') (XSNAME(I),I=8,NXMOL)                       
      ELSE                                                                
         READ (IRD,'(7A10)') (XSNAME(I),I=1,NXMOL)                       
      ENDIF                                                               
C                                                                         
C     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS           
C     AND DETERMINE THE INDEX VALUE.  
      IXMAX = 4                                                          
      DO 40 I = 1, NXMOL                                                 
C        Left-justify all inputed names.                                      
         CALL CLJUST (XSNAME(I),10)
         IXINDX(I) = 0
         DO 20 J = 1, IXMAX
            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(2,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(3,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(4,J))) THEN                           
               IXINDX(I) = J                                              
            ENDIF                                                         
   20    CONTINUE
   40 CONTINUE                                                            

      RETURN
      END

      BLOCK DATA

      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)
      PARAMETER (MXLAY=203)
      PARAMETER (MAXXSEC=4)
      PARAMETER (MAXPROD = MXLAY*MAXXSEC)

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)


      COMMON /CVRRTM/    HNAMRTM,HVRRTM
      COMMON /CVRRTR/    HNAMRTR,HVRRTR
      COMMON /CVRATM/    HNAMATM,HVRATM
      COMMON /CVRSET/    HNAMSET,HVRSET
      COMMON /CVRTAU/    HNAMTAU,HVRTAU
      COMMON /CVRINI/    HNAMINI,HVRINI
      COMMON /CVRCLD/    HNAMCLD,HVRCLD
      COMMON /CVRUTL/    HNAMUTL,HVRUTL
      COMMON /CVREXT/    HNAMEXT,HVREXT
      COMMON /CVRRTX/    HNAMRTX,HVRRTX

      COMMON /CVRSN/     HNAMKG,HVRKG

      CHARACTER*18 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *             HVRINI,HVRCLD,HVRUTL,HVREXT,HVRRTX

      CHARACTER*18 HNAMRTM,HNAMRTR,HNAMATM,HNAMSET,HNAMTAU,
     *             HNAMINI,HNAMCLD,HNAMUTL,HNAMEXT,HNAMRTX

      CHARACTER*18 HVRKG

      CHARACTER*18 HNAMKG

      DATA HVRRTM /'NOT USED'/, HVRRTR /'NOT USED'/,
     *     HVRATM /'NOT USED'/, HVRSET /'NOT USED'/,
     *     HVRTAU /'NOT USED'/, HVRINI /'NOT USED'/,
     *     HVRCLD /'NOT USED'/, HVREXT /'NOT USED'/,
     *     HVRUTL /'NOT USED'/, HVRRTX /'NOT USED'/,
     *     HVRKG  /'        '/

      DATA HNAMRTM / '       rrtmg_lw.f:' /,
     *     HNAMINI / '  rrtmg_lw_init.f:' /,
     *     HNAMCLD / 'rrtmglw_cldprop.f:' /,
     *     HNAMRTR / '  rrtmg_lw_rtrn.f:' /,
     *     HNAMRTX / 'rrtmg_lw_rtrnmr.f:' /,
     *     HNAMSET / 'rrtmglw_setcoef.f:' /,
     *     HNAMTAU / 'rrtmg_lw_taumol.f:' /,
     *     HNAMATM / '         rrtatm.f:' /,
     *     HNAMUTL / '       util_xxx.f:' /,
     *     HNAMEXT / '          extra.f:' /,
     *     HNAMKG  / '   rrtmg_lw_k_g.f:' /

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
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/

      DATA WX /MAXPROD*0.0/

      END
c**********************************************************************
      Block Data phys_consts
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
c
      DATA PI / 3.1415926535 /
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     Pi was obtained from   PI = 2.*ASIN(1.)         
c
c     units are genrally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c              RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07
c              RADCN2 = PLANCK*CLIGHT/BOLTZ           
      end
c



