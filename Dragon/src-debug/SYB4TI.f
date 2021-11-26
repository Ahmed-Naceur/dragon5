*DECK SYB4TI
      SUBROUTINE SYB4TI (NHMAX,IXRAYO,IS1,NSECT4,IFAC,NUMREG,NLMAX,
     & IREGI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the interception region indices of a track in a rectangular cell.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NHMAX   number of interceptions.
* IXRAYO  tube indices.
* IS1     index of the first sector.
* NSECT4  number of sectors.
* IFAC    index of the symmetry.
* NUMREG  region indices of the tube sectors.
*
*Parameters: output
* NLMAX   number of interception regions.
* IREGI   indices of the interception regions.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NHMAX,IXRAYO(NHMAX),IS1,NSECT4,IFAC,NUMREG(NSECT4,*),
     &        NLMAX,IREGI(NHMAX)
*----
*  LOCAL VARIABLES
*----
*:IHE    No de l'abcisse Est
*:IREGC  No de la Region Courante
*:IREGS  No de la Region Suivante
*:ISC    No du Secteur Courant
*:IS2    No du Secteur Suivant
*:JRC    No de la Couronne Courante
*:JRS    No de la Couronne Suivante
*
      JRC = IXRAYO(1)
      ISC = IS1
      ISF = ISC
      IF (NSECT4 .GT. 1) THEN
        IF (IFAC .EQ. 1) THEN
          ISF = 1 - ISF
        ELSEIF (IFAC .EQ. 2) THEN
          ISF = ISF - NSECT4 / 2
        ELSEIF (IFAC .EQ. 3) THEN
          ISF = 1 + NSECT4 / 2 - ISF
        ENDIF
        IF (ISF .LE. 0) ISF = ISF + NSECT4
      ENDIF
      IREGC = NUMREG(ISF, JRC) + 3
*
* Dernier Intervalle
      NLMAX = 1
*----
*  Boucle des Volumes internes
*  Debut
*----
      DO IHE = 2, NHMAX
*
* No de Couronne Suivante
        JRS = IXRAYO(IHE)
*
* Soit  : Meme Couronne => Changement de Secteur
* Sinon : Changement de Couronne => No de Rayon
        IF (JRC .EQ. JRS) THEN
          IF (ISC .LT. NSECT4) THEN
            IS2 = ISC + 1
          ELSE
            IS2 = 1
          ENDIF
        ELSE
          IS2 = ISC
        ENDIF
*
* +++
* Debut
* Region Suivante
        ISF = IS2
        IF (NSECT4 .GT. 1) THEN
          IF (IFAC .EQ. 1) THEN
            ISF = 1 - ISF
          ELSEIF (IFAC .EQ. 2) THEN
            ISF = ISF - NSECT4 / 2
          ELSEIF (IFAC .EQ. 3) THEN
            ISF = 1 + NSECT4 / 2 - ISF
          ENDIF
          IF (ISF .LE. 0) ISF = ISF + NSECT4
        ENDIF
        IREGS = NUMREG(ISF, JRS) + 3
        IF (IREGS .NE. IREGC) THEN
          IREGI(NLMAX) = IREGC
          NLMAX = NLMAX + 1
          IREGC = IREGS
        ENDIF
* Fin Region Suivante
* +++
*
* Suivants
        ISC = IS2
        JRC = JRS
*
* - - - - - - - - - - - - -
* Boucle des Volumes internes
* Fin
* - - - - - - - - - - - - -
      ENDDO
*
      IREGI(NLMAX) = IREGC
*
      RETURN
      END
