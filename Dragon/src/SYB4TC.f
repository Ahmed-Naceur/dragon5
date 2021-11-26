*DECK SYB4TC
      SUBROUTINE SYB4TC (DELTAR,DDELTA,ANGLES,NHMAX,IXRAYO,IS1,NSECT4,
     & IFAC,NUMREG,RAYONS,ZZW,ZZE,ZZR,HXRAYO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the track lengths and interception lengths in a rectangular
* cell.
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
* DELTAR  used to compute the interception.
* DDELTA  undefined.
* ANGLES  angular values.
* NHMAX   number of interceptions.
* IXRAYO  tube indices.
* IS1     index of the first sector.
* NSECT4  number of sectors.
* IFAC    undefined.
* NUMREG  region indices of the tube sectors.
* RAYONS  radius.
* ZZW     position of the west interception (left).
* ZZE     position of the east interception (right).
*
*Parameters: input/output
* ZZR     tracking information.
* HXRAYO  preceding/next interception lengths.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NHMAX,IXRAYO(NHMAX),IS1,NSECT4,IFAC,NUMREG(NSECT4,*)
      REAL    DELTAR,DDELTA,ANGLES(NSECT4),RAYONS(*),ZZW,ZZE,ZZR(*),
     &        HXRAYO(NHMAX+1)
*----
*  LOCAL VARIABLES
*----
*:AE     Aire de la courbure Est  (a Droite)
*:AW     Aire de la courbure West (a Gauche)
*:GE     Abcisse du point Est  (Trajectoire precedente)
*:GW     Abcisse du point West (Trajectoire precedente)
*:HE     Abcisse du point Est
*:HW     Abcisse du point West (Intersection)
*:XW     Abcisse du point West (Region)
*:IHE    No de l'abcisse Est
*:IRE    No du Rayon Cote Est
C
*:IL     No de l'Intersection
*:IREGC  No de la Region Courante
*:IREGS  No de la Region Suivante
*:ISC    No du Secteur Courant
*:IS2    No du Secteur Suivant
*:JRC    No de la Couronne Courante
*:JRS    No de la Couronne Suivante
*
      DDELT2 = DDELTA * DDELTA
      DELTA2 = DELTAR * DELTAR
*
      AW  = 0.
      GW  = HXRAYO(1)
      HW  = ZZW
      XW  = HW
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
      IREGC = NUMREG(ISF, JRC)
      IL  = 0
*
      HXRAYO(1) = HW
*----
* Boucle des Volumes internes
* Debut
*----
      DO IHE = 2, NHMAX
*
* No de Couronne Suivante
        JRS = IXRAYO(IHE)
*
* Soit  : Meme Couronne => Changement de Secteur
* Sinon : Changement de Couronne => No de Rayon
        IRE = 0
        IF (JRC .EQ. JRS) THEN
          IF (ISC .LT. NSECT4) THEN
            IS2 = ISC + 1
          ELSE
            IS2 = 1
          ENDIF
        ELSE
          IS2 = ISC
          IRE = MIN(JRC, JRS)
        ENDIF
*
* Soit  Changement de Secteur
        IF (JRC .EQ. JRS) THEN
          HE = DELTAR * ANGLES(ISC)
*
* Sinon Changement de Couronne
        ELSE
          H2 = RAYONS(IRE) * RAYONS(IRE) - DELTA2
          IF (H2 .GT. 0) THEN
            HE = SQRT(H2)
            IF (JRS .EQ. IRE) THEN
              HE = - HE
            ENDIF
          ELSE
            HE = 0.
          ENDIF
*
        ENDIF
*
* Protection contre les longueurs negatives
        IF (HE .LT. HW) THEN
          HE = HW
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
        IREGS = NUMREG(ISF, JRS)
        IF (IREGS .NE. IREGC) THEN
*
* Nouvelles Valeurs
        AE  = 0.
        GE  = HXRAYO(IHE)
*
* Moyenne entre l'intervalle precedent et Actuel
        ZZH = HE - XW
        ZZH = (ZZH + GE - GW) * 0.5
*
* Ajout de la courbure
        IF (JRC .NE. JRS) THEN
          H2 = HE - GE
          H2 = H2 * H2
          XCORDE = H2 + DDELT2
          IF (XCORDE .GT. 0.) THEN
            XCORDE = SQRT(XCORDE)
            XUNITE = XCORDE / RAYONS(IRE)
            XUNITE = XUNITE / 2.
            XALPHA = ASIN(XUNITE)
            XUNITE = XALPHA - COS(XALPHA) * XUNITE
            AE = XUNITE * RAYONS(IRE) * RAYONS(IRE)
            AE  = AE / DDELTA
          ENDIF
*
          IF (JRS .EQ. IRE) THEN
            AE = - AE
          ENDIF
        ENDIF
*
* Longueur Moyenne
* Nouvelle Abcisse
        IL = IL + 1
        ZZR(IL) = ZZH + AE - AW
C
C Valeurs Precedentes (Region)
        AW  = AE
        GW  = GE
        XW  = HE
        ENDIF
* Fin Region Suivante
* +++
*
* Valeurs Precedentes (Intersection)
        HXRAYO(IHE) = HE
        HW  = HE
*
* Suivants
        IREGC = IREGS
        ISC = IS2
        JRC = JRS
*
* - - - - - - - - - - - - -
* Boucle des Volumes internes
* Fin
* - - - - - - - - - - - - -
      ENDDO
*
* - - - - - - - - - - - - -
* Dernier Intervalle
* Debut
* - - - - - - - - - - - - -
*
* Intervalle
        ZZH = ZZE - XW
*
* Protection contre les Volumes Negatifs
        ZZH = MAX(ZZH, 0.)
*
* Moyenne entre l'intervalle precedent et Actuel
        ZZH = (ZZH + HXRAYO(NHMAX+1) - GW) * 0.5
*
* Ajout de la Courbure Eventuelle
* Mise a jour de la Nouvelle abcisse
        IL = IL + 1
        ZZR(IL) = ZZH - AW
        HXRAYO(NHMAX+1) = ZZE
* - - - - - - - - - - - - -
* Dernier Intervalle
* Fin
* - - - - - - - - - - - - -
*
      RETURN
      END
