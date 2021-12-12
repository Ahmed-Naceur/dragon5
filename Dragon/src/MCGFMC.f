*DECK MCGFMC
      SUBROUTINE MCGFMC(KPN,K,NREG,M,NANI,NFUNL,NZON,KEYFLX,KEYCUR,
     1                  PHIOUT,V,S,ST,KEYANI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Addition of the regional source to the flux when the 'MOCC/MCI'
* integration strategy is turned on for the method of characteristics
* integration.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* KPN     total number of unknowns per group in flux vector.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* NZON    index-number of the mixture type assigned to each volume.
* KEYFLX  position of flux elements in flux vector.
* KEYCUR  position of current elements in flux vector.
* V       volumes.
* S       source vector.
* ST      total cross sections array.
* KEYANI  'mode to l' index: l=KEYANI(nu).
*
*Parameters: input/output
* PHIOUT  flux vector. 
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER KPN,K,NREG,M,NANI,NFUNL,NZON(K),KEYFLX(NREG,NFUNL),
     1 KEYCUR(K-NREG),KEYANI(NFUNL)
      REAL V(K),ST(0:M)
      DOUBLE PRECISION PHIOUT(KPN),S(KPN)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,IBM,IL,IND
*
      IF(NANI.LE.0) CALL XABORT('MCGFMC: INVALID VALUE OF NANI.')
      DO I=1,K
         IF(V(I).GT.0.) THEN 
            IBM=NZON(I)
            IF (IBM.LT.0) THEN
               IND=KEYCUR(I-NREG)
               PHIOUT(IND)=PHIOUT(IND)/DBLE(V(I))
            ELSE
               DO IL=1,NFUNL
                  IND=KEYFLX(I,IL)
                  PHIOUT(IND)=PHIOUT(IND)/DBLE(V(I)*ST(IBM))
     1                       +S(IND)/DBLE(2*KEYANI(IL)+1)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
*
      RETURN
      END
