*DECK SCRSXS
      SUBROUTINE SCRSXS(NGRP,NL,NREA,IREAF,NOMREA,LXS,B2SAP,FACT0,
     1 WEIGHT,SPH,FLUXS,XSB,SIGSB,SS2DB,LPURE,XS,SIGS,SS2D,TAUXFI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Weight microscopic cross section data in an interpolated microlib.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* NGRP    number of energy groups
* NL      maximum Legendre order (NL=1 is for isotropic scattering)
* NREA    number of reactions in the Saphyb object
* IREAF   position of 'NU*FISSION' reaction in NOMREA array
* NOMREA  names of reactions in the Saphyb object
* LXS     existence flag of each reaction
* B2SAP   buckling as recovered from the Saphyb object
* FACT0   number density ratio for the isotope
* WEIGHT  interpolation weight
* SPH     SPH factors
* FLUXS   averaged flux
* XSB     cross sections per reaction for a unique calculation
* SIGSB   scattering cross sections for a unique calculation
* SS2DB   scattering matrix for a unique calculation
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
*
*Parameters: input/output
* XS      interpolated cross sections per reaction
* SIGS    interpolated scattering cross sections
* SS2D    interpolated scattering matrix
* TAUXFI  interpolated fission rate
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,NL,NREA,IREAF
      INTEGER I, IGR, IL, IREA, IRF, J, JGR
      REAL TAUXF, XSECT
      REAL B2SAP,FACT0,WEIGHT,SPH(NGRP),FLUXS(NGRP),XSB(NGRP*NREA),
     1 SIGSB(NGRP*NL),SS2DB(NL*NGRP*NGRP),XS(NGRP*NREA),SIGS(NGRP*NL),
     2 SS2D(NGRP*NGRP*NL),TAUXFI
      CHARACTER NOMREA(NREA)*12
      LOGICAL LXS(NREA),LPURE
*----
*  COMPUTE FISSION RATE FOR AN ELEMENTARY CALCULATION
*----
      TAUXF=0.0
      IF(.NOT.LPURE.AND.(IREAF.GT.0)) THEN
        DO IGR=1,NGRP
          IRF=(IREAF-1)*NGRP+IGR
          TAUXF=TAUXF+XSB(IRF)*FLUXS(IGR)
        ENDDO
        TAUXFI=TAUXFI+WEIGHT*FACT0*TAUXF
      ENDIF
*----
*  MICROLIB INTERPOLATION
*----
      DO IGR=1,NGRP
        DO IREA=1,NREA
          IF(.NOT.LXS(IREA)) CYCLE
          I=(IREA-1)*NGRP+IGR
          IF(LPURE.AND.NOMREA(IREA).EQ.'SPECTRE') THEN
            XS(I)=XS(I)+WEIGHT*XSB(I)
          ELSE IF(NOMREA(IREA).EQ.'SPECTRE') THEN
            XS(I)=XS(I)+WEIGHT*FACT0*TAUXF*XSB(I)
          ELSE IF(NOMREA(IREA).EQ.'FUITES') THEN
            IF(B2SAP.NE.0.0) THEN
              XSECT=XSB(I)/B2SAP
              XS(I)=XS(I)+SPH(IGR)*FACT0*WEIGHT*XSECT
            ENDIF
          ELSE IF(NOMREA(IREA).EQ.'TOTALE P1') THEN
            XS(I)=XS(I)+FACT0*WEIGHT*XSB(I)/SPH(IGR)
          ELSE
            XS(I)=XS(I)+FACT0*SPH(IGR)*WEIGHT*XSB(I)
          ENDIF
        ENDDO
        DO IL=1,NL
          I=(IL-1)*NGRP+IGR
          IF(MOD(IL,2).EQ.1) THEN
            SIGS(I)=SIGS(I)+FACT0*SPH(IGR)*WEIGHT*SIGSB(I)
          ELSE
            DO JGR=1,NGRP
              J=(IL-1)*NGRP*NGRP+(IGR-1)*NGRP+JGR
              SIGS(I)=SIGS(I)+FACT0*WEIGHT*SS2DB(J)/SPH(JGR)
            ENDDO
          ENDIF
        ENDDO
        DO JGR=1,NGRP
          DO IL=1,NL
            I=(IL-1)*NGRP*NGRP+(JGR-1)*NGRP+IGR
            IF(MOD(IL,2).EQ.1) THEN
              SS2D(I)=SS2D(I)+FACT0*SPH(JGR)*WEIGHT*SS2DB(I)
            ELSE
              SS2D(I)=SS2D(I)+FACT0*WEIGHT*SS2DB(I)/SPH(IGR)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
