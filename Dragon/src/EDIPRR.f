*DECK EDIPRR
      SUBROUTINE EDIPRR(IPRINT,NL,ITRANC,NGCOND,NMERGE,ILEAKS,NW,NTAUXT,
     >                  B2,VOLMER,NENER,WENERG,RATECM,FLUXCM,SCATTD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print reaction rates.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPRINT  print level;
*         = 0 no print;
*         = 1 print fluxes;
*         = 2 1+print reaction rates;
*         = 3 2+print homogenized cross sections.
* NL      number of Legendre orders.
* ITRANC  type of transport correction.
* NGCOND  number of condensed groups.
* NMERGE  number of merged regions.
* ILEAKS  type of leakage calculation:
*         = 0 no leakage;
*         = 1 homogeneous leakage (Diffon);
*         = 2 isotropic streaming (Ecco);
*         = 3 anisotropic streaming (Tibere);
*         = 10 isotropic diffusion coefficients recovered from input
*           macrolib;
*         = 11 anisotropic diffusion coefficients recovered from input
*           macrolib.
* NW      type of weighting for PN cross section info (=0 P0; =1 P1).
* NTAUXT  number of reaction rate edits.
* B2      square buckling:
*         for ILEAKS=1,2: B2(4) is homogeneous;
*         for ILEAKS=3: B2(1),B2(2),B2(3) are directional heterogeneous
*         and B2(4) is homogeneous.
* VOLMER  volume of region merged.
* NENER   number of energy groups limits.
* WENERG  energy group limits.
* RATECM  averaged region/group cross sections:
*         = RATECM(*,1) = total P0;
*         = RATECM(*,2) = total P1;
*         = RATECM(*,NW+2) = absorption;
*         = RATECM(*,NW+3) = fission;
*         = RATECM(*,NW+4) = fixed sources / productions;
*         = RATECM(*,NW+5) = leakage;
*         = RATECM(*,NW+6) = total out of group scattering;
*         = RATECM(*,NW+7) = diagonal scattering x-s;
*         = RATECM(*,NW+8) = chi;
*         = RATECM(*,NW+9) = wims type transport correction;
*         = RATECM(*,NW+10) = x-directed leakage;
*         = RATECM(*,NW+11) = y-directed leakage;
*         = RATECM(*,NW+12) = z-directed leakage;
*         = RATECM(*,NW+13) = nu-sigf for delayed neutrons;
*         = RATECM(*,NW+13+NDEL) = fission spectra for delayed neutrons.
* FLUXCM  integrated region/group fluxes:
*         = FLUXCM(*,1) = fluxes P0;
*         = FLUXCM(*,2) = fluxes P1.
* SCATTD  scattering rates.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    IPRINT,NL,ITRANC,NGCOND,NMERGE,ILEAKS,NW,NTAUXT,NENER
      REAL       B2(4),VOLMER(NMERGE),WENERG(NGCOND+1),
     >           RATECM(NMERGE,NGCOND,NTAUXT),FLUXCM(NMERGE,NGCOND,NW+1)
      DOUBLE PRECISION SCATTD(NMERGE,NGCOND,NGCOND,NL)
*----
*  LOCAL VARIABLES
*----
      CHARACTER   APG*3
      PARAMETER  (IUNOUT=6,APG=' > ')
      DOUBLE PRECISION SCATWG,SCATTN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLDMC
*----
*  SCRATCH STORAGE ALLOCATION
*   FLDMC   flux merged and condensed.
*----
      ALLOCATE(FLDMC(NMERGE,NGCOND))
*----
*  COMPUTE AVERAGE FLUX
*----
      DO 224 IGRC=1,NGCOND
        DO 225 IKK=1,NMERGE
          FLDMC(IKK,IGRC)=FLUXCM(IKK,IGRC,1)/VOLMER(IKK)
 225    CONTINUE
 224  CONTINUE
*----
*  PRINT REACTION RATES
*----
      WRITE(IUNOUT,6000)
      WRITE(IUNOUT,6001) (JJ,VOLMER(JJ),JJ=1,NMERGE)
      IF( (NENER.GT.0) .AND. (IPRINT.GT.1) ) THEN
        WRITE(IUNOUT,6002) (WENERG(IG),APG,IG,APG,IG=1,NGCOND),
     >                        WENERG(NGCOND+1)
      ENDIF
      WRITE(IUNOUT,6003)
      DO 154 IGR=1,NGCOND
        IF(IPRINT.EQ.1) THEN
          WRITE(IUNOUT,6010) IGR
          WRITE(IUNOUT,6012) (FLUXCM(IKK,IGR,1),IKK=1,NMERGE)
          WRITE(IUNOUT,6011)
          WRITE(IUNOUT,6012) (FLDMC(IKK,IGR),IKK=1,NMERGE)
          GO TO 154
        ENDIF
        IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2).OR.(ILEAKS.EQ.10)) THEN
          WRITE(IUNOUT,6013) IGR
        ELSE
          WRITE(IUNOUT,6014) IGR
        ENDIF
        DO 155 IKK=1,NMERGE
*----
*  UNCOMMENT THE 2 LINES TO PERFORM TRANSPORT CORRECTION
*----
          TOTAL=RATECM(IKK,IGR,1)
          SCATWG=SCATTD(IKK,IGR,IGR,1)
          IF(ITRANC.NE.0) THEN
*           TOTAL=TOTAL-RATECM(IKK,IGR,NW+9)
*           SCATWG=SCATWG-RATECM(IKK,IGR,NW+9)
          ENDIF
*
          SCATTN=0.0D0
          DO 153 JGR=1,NGCOND
            IF(JGR.NE.IGR) SCATTN=SCATTN+SCATTD(IKK,JGR,IGR,1)
 153      CONTINUE
          IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2).OR.(ILEAKS.EQ.10)) THEN
            WRITE(IUNOUT,6020) IKK,FLDMC(IKK,IGR),FLUXCM(IKK,IGR,1),
     >        TOTAL,RATECM(IKK,IGR,NW+2),RATECM(IKK,IGR,NW+3),
     >        RATECM(IKK,IGR,NW+5)*B2(4),RATECM(IKK,IGR,NW+4),SCATWG,
     >        SCATTN
          ELSE
            WRITE(IUNOUT,6021) IKK, FLDMC(IKK,IGR),FLUXCM(IKK,IGR,1),
     >        TOTAL,RATECM(IKK,IGR,NW+2),RATECM(IKK,IGR,NW+3),
     >        RATECM(IKK,IGR,NW+4),SCATWG,SCATTN
          ENDIF
 155    CONTINUE
        IF((ILEAKS.EQ.3).OR.(ILEAKS.EQ.11)) THEN
          WRITE(IUNOUT,6022)
          DO 156 IKK=1,NMERGE
            WRITE(IUNOUT,6023) IKK,RATECM(IKK,IGR,NW+10)*B2(1)+
     >      RATECM(IKK,IGR,NW+11)*B2(2)+RATECM(IKK,IGR,NW+12)*B2(3),
     >      RATECM(IKK,IGR,NW+10)*B2(1),RATECM(IKK,IGR,NW+11)*B2(2),
     >      RATECM(IKK,IGR,NW+12)*B2(3),RATECM(IKK,IGR,NW+5)*B2(4)
 156      CONTINUE
        ENDIF
 154  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLDMC)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(////5(5X,'REGION',6X,'VOLUME  '))
 6001 FORMAT(1P,5(5X,I4,4X,E12.5))
 6002 FORMAT(/' E N E R G Y   L I M I T S   (EV)'/1P,
     >6(E12.4,A3,I3,A3))
 6003 FORMAT(/' F L U X E S   A N D    R E A C T I O N    R A T E S'/
     >1X,51(1H-))
 6010 FORMAT(/' G R O U P   :',I4/' REGION INTEGRATED FLUX')
 6011 FORMAT(' AVERAGED REGIONAL FLUX')
 6012 FORMAT(1P,7(3X,E15.7))
 6013 FORMAT(/14H G R O U P   :,I4/7H REGION,3X,7HAVERAGE,5X,3HINT,
     > 7HEGRATED,5X,9HCOLLISION,4X,10HABSORPTION,4X,10HNU*FISSION,6X,
     > 7HLEAKAGE,5X,10HPRODUCTION,8X,16HSCATTERING RATES/11X,4HFLUX,
     > 10X,4HFLUX,10X,4HRATE,10X,4HRATE,10X,4HRATE,10X,4HRATE,10X,
     > 4HRATE,6X,26HWITHIN GROUP  OUT OF GROUP)
 6014 FORMAT(/14H G R O U P   :,I4/7H REGION,3X,7HAVERAGE,5X,3HINT,
     > 7HEGRATED,5X,9HCOLLISION,4X,10HABSORPTION,4X,10HNU*FISSION,
     > 4X,10HPRODUCTION,8X,16HSCATTERING RATES/11X,4HFLUX,
     > 10X,4HFLUX,10X,4HRATE,10X,4HRATE,10X,4HRATE,10X,
     > 4HRATE,6X,26HWITHIN GROUP  OUT OF GROUP)
 6020 FORMAT(1X,I4,1P,9E14.5)
 6021 FORMAT(1X,I4,1P,8E14.5,3E14.5)
 6022 FORMAT(/' REGION  TOTAL LEAKAGE     X-LEAKAGE',
     >        '     Y-LEAKAGE     Z-LEAKAGE   HOMOGENEOUS'/
     >        '              RATE           RATE   ',
     >        '       RATE          RATE     LEAKAGE RATE')
 6023 FORMAT(1X,I6,1X,1P,5E14.5)
      END
