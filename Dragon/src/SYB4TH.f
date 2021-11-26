*DECK SYB4TH
      SUBROUTINE SYB4TH (NRMAX,RAYONS,PENTES,ORIGIN,ISCMIN,ISCMAX,
     & ANGLES,DXMIN,DELTAC,ISCW,ISCE,ISXW,ISXE,NHMAX,IXRAYO,HXRAYO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the intersection lenghts of a track in a sectorized Cartesian
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
* NRMAX   number of radius.
* RAYONS  radius of each tube.
* PENTES  slope of the track with respect to the axis or to the sides.
* ORIGIN  origin of the track.
* ISCMIN  first sector with a possible intersection.
* ISCMAX  last sector with a possible intersection.
* ANGLES  slope of the sectors.
* DXMIN   accuracy.
* DELTAC  position of the track.
*
*Parameters: output
* ISCW    index of the first sector (west/left).
* ISCE    index of the last sector (east/right).
* ISXW    index of the first side (west/left). Equal to 1 or 2.
* ISXE    index of the last side (east/right). Equal to 3 or 4.
* NHMAX   number of intersections.
* IXRAYO  indices of the intersecting tubes.
* HXRAYO  position of the intersections (boundary).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NRMAX,ISCMIN,ISCMAX,ISCW,ISCE,ISXW,ISXE,NHMAX,
     &        IXRAYO(NHMAX)
      REAL    RAYONS(NRMAX),PENTES(4),ORIGIN(4),ANGLES(ISCMAX),DXMIN,
     &        DELTAC,HXRAYO(NHMAX+1)
*----
*  LOCAL VARIABLES
*----
      CHARACTER * 4 AJOU(2)
      LOGICAL AJRAYO
      LOGICAL AJSECT
*
* Remarque Importante :
*  Pour choisir entre deux lignes Proches
*    Il faut prendre le cas correspondant a
* DELTAC + DXMIN
* ZHC    + DXMIN
*
      DXMAX = 2.
      D2 = DELTAC * DELTAC
      DR = DELTAC + DXMIN
      ZH1 = DR * PENTES(1) + ORIGIN(1)
      ZH2 = DR * PENTES(2) + ORIGIN(2)
      IF (ZH2 .GE. (ZH1 - DXMIN)) THEN
        ZHW = DELTAC * PENTES(2) + ORIGIN(2)
        ISXW = 2
      ELSE
        ZHW = DELTAC * PENTES(1) + ORIGIN(1)
        ISXW = 1
      ENDIF
      ZH3 = DELTAC * PENTES(3) + ORIGIN(3)
      ZH4 = DELTAC * PENTES(4) + ORIGIN(4)
      IF (ZH3 .LE. (ZH4 + DXMIN)) THEN
        ZHE  = ZH3
        ISXE = 3
      ELSE
        ZHE  = ZH4
        ISXE = 4
      ENDIF
*
* Recherche du Premier Secteur
* L'Angle ISC delimite les Secteurs ISC,ISC+1
      DO 101 ISCW = ISCMIN, ISCMAX
        ZHS  = DELTAC * ANGLES(ISCW)
        IF (ZHS .GT. (ZHW + DXMIN)) THEN
          GOTO 102
        ENDIF
  101 CONTINUE
        ISCW = ISCMAX + 1
        ZHS = ZHE + DXMAX
  102 CONTINUE
*
* Couronne Interne
      DO 201 NRMIN = 1, NRMAX
        IF (RAYONS(NRMIN) .GT. DR) THEN
          R2 = RAYONS(NRMIN) * RAYONS(NRMIN)
          H2 = R2 - D2
          IF (H2 .GT. 0.) THEN
            H0 = SQRT(H2)
            GOTO 202
          ENDIF
        ENDIF
  201 CONTINUE
      H0 = 0.
      NRMIN = NRMAX + 1
  202 CONTINUE
* NRMIN = Couronne la Plus Interne
* NRMIN = Le Plus Petit Rayon a Prendre en Compte
* DXMIN - H0 = Limite entre Couronnes West et Est
*
* Recherche du Premier Rayon
      H2 = ZHW * ZHW
      R2 = D2 + H2
      RW = SQRT(R2)
*
* Pente pour le Rayon West
* Pente pour le Cote
*     PR = - DELTAC / ZHW
*     PC = PENTES(ISXW)
*
* -------------
* Debut
* Couronne West
* -------------
*
* Aucune Couronne
      IF (NRMIN .GT. NRMAX) THEN
        ZHR = ZHE + DXMAX
        ICC = NRMIN
*
* Couronne a Gauche (ZHW <0)
      ELSEIF (ZHW .LE. (- 0.5 * H0)) THEN
        IF ((- ZHW * PENTES(ISXW)) .GT. DR) THEN
          RW = RW - DXMIN
        ELSE
          RW = RW + DXMIN
        ENDIF
        DO 221 ICC = NRMAX, NRMIN, -1
          ZRR = RAYONS(ICC)
          IF (ZRR .LE. RW) THEN
            GOTO 222
          ENDIF
  221   CONTINUE
        ICC = NRMIN - 1
  222   CONTINUE
        IF (ICC .LT. NRMIN) THEN
          ZHR = H0
        ELSE
          R2 = ZRR * ZRR
          H2 = R2 - D2
          ZHR = - SQRT(H2)
        ENDIF
        ICC = ICC + 1
*
* Couronne Centrale
      ELSEIF (ZHW .LE. (0.5 * H0)) THEN
        ICC = NRMIN
        ZHR = H0
*
* Couronne a Droite (ZHW >0)
      ELSE
        IF ((ZHW * PENTES(ISXW)) .GT. - DR) THEN
          RW = RW + DXMIN
        ELSE
          RW = RW - DXMIN
        ENDIF
        DO 231 ICC = NRMIN, NRMAX
          ZRR = RAYONS(ICC)
          IF (ZRR .GE. RW) THEN
            GOTO 232
          ENDIF
  231   CONTINUE
        ICC = NRMAX + 1
  232   CONTINUE
        IF (ICC .GT. NRMAX) THEN
          ZHR = ZHE + DXMAX
        ELSE
          R2 = ZRR * ZRR
          H2 = R2 - D2
          ZHR = SQRT(H2)
        ENDIF
*
      ENDIF
* -------------
* Couronne West
* Fin
* -------------
*
* Premiere Position Courante
      IHC = 1
      ISC = ISCW
      ZHC = ZHW
*
      IXRAYO(IHC) = ICC
      HXRAYO(IHC) = ZHC
*
* Cote West Intersecte Rayon
      IF (ZHR .LE. (ZHC + DXMIN)) THEN
        IF (ZHR .GT. 0.) THEN
          ICC = ICC + 1
        ELSE
          ICC = ICC - 1
        ENDIF
        ZHC = MAX (ZHR, ZHC)
*
        IHC = IHC + 1
        IXRAYO(IHC) = ICC
        HXRAYO(IHC) = ZHC
*
        IF (ICC .GT. NRMAX) THEN
          ZHR = ZHE + DXMAX
        ELSEIF (ICC .EQ. NRMIN) THEN
          ZHR = H0
        ELSEIF (ZHR .LT. 0.) THEN
          R2 = RAYONS(ICC-1) * RAYONS(ICC-1)
          H2 = R2 - D2
          ZHR = - SQRT(H2)
        ELSE
          R2 = RAYONS(ICC) * RAYONS(ICC)
          H2 = R2 - D2
          ZHR =   SQRT(H2)
        ENDIF
      ENDIF
*
* Suivant
  555 CONTINUE
      ZHN = MIN (ZHE, ZHR, ZHS)
*
* +++
* Debut
* point Courant (<Est)
      IF ((ZHN + DXMIN) .LT. ZHE) THEN
*
* Increment Rayon
* Increment Secteur
        AJRAYO = (ZHN + DXMIN) .GE. ZHR
        AJSECT = (ZHN + DXMIN) .GE. ZHS
*
* soit Intersection Rayon et Secteur
* soit Intersection SEULEMENT Rayon
* soit Intersection SEULEMENT Secteur (non Rayon, non Cote)
*
        IF (AJRAYO) THEN
        IF (AJSECT) THEN
          NBAJ = 2
* Choix pentes
*     PR = - DELTAC / ZHN
*     PS = ANGLES(ISC)
          IF (ZHN .GT. 0.) THEN
            ZZ = ZHN * ANGLES(ISC)
            IF (ZZ .LT. - DR) THEN
              AJOU(1) = 'Sect'
              AJOU(2) = 'Rayo'
            ELSE
              AJOU(1) = 'Rayo'
              AJOU(2) = 'Sect'
            ENDIF
          ELSE
            ZZ = - ZHN * ANGLES(ISC)
            IF (ZZ .LT. DR) THEN
              AJOU(1) = 'Sect'
              AJOU(2) = 'Rayo'
            ELSE
              AJOU(1) = 'Rayo'
              AJOU(2) = 'Sect'
            ENDIF
          ENDIF
* Autre Cas
* =soit Intersection SEULEMENT Rayon   (non Secteur)
* =soit Intersection SEULEMENT Secteur (non Rayon, non Cote)
        ELSE
          NBAJ = 1
          AJOU(1) = 'Rayo'
        ENDIF
        ELSE
          NBAJ = 1
          AJOU(1) = 'Sect'
        ENDIF
*
        DO 321 IAJ = 1, NBAJ
*
        IHC = IHC + 1
*
* Avance Rayon
        IF (AJOU(IAJ) .EQ. 'Rayo') THEN
*
          ZZ = MAX (ZHR,  ZHC)
          HXRAYO(IHC) = ZZ
          ZHC = ZZ
*
          IF (ZHR .GT. 0.) THEN
            ICC = ICC + 1
            IF (ICC .GT. NRMAX) THEN
              ZHR = ZHE + DXMAX
            ELSE
              R2 = RAYONS(ICC) * RAYONS(ICC)
              H2 = R2 - D2
              ZHR =   SQRT(H2)
            ENDIF
          ELSE
            ICC = ICC - 1
            IF (ICC .EQ. NRMIN) THEN
              ZHR = H0
            ELSE
              R2 = RAYONS(ICC-1) * RAYONS(ICC-1)
              H2 = R2 - D2
              ZHR = - SQRT(H2)
            ENDIF
          ENDIF
          IXRAYO(IHC) = ICC
*
* Avance Secteur
        ELSE
*
          ZZ = MAX (ZHS,  ZHC)
          HXRAYO(IHC) = ZZ
          ZHC = ZZ
*
          IXRAYO(IHC) = ICC
          ISC = ISC + 1
          IF (ISC .GT. ISCMAX) THEN
            ZHS = ZHE + DXMAX
          ELSE
            ZHS  = DELTAC * ANGLES(ISC)
          ENDIF
        ENDIF
  321   CONTINUE
*
        GOTO 555
*
      ENDIF
* point Courant (<Est)
* Fin
* +++
*
* Debut
* point Final (Est)
*
* Intersection Couronne
      IF ((ZHN + DXMIN) .GE. ZHR) THEN
* ? Pentes
*     PR = - DELTAC / ZHR
*     PC = PENTES(ISXE)
        IF (ZHR .GT. 0.) THEN
          ZZ = ZHR * PENTES(ISXE)
          AJRAYO = ZZ .GT. - DR
        ELSE
          ZZ = - ZHR * PENTES(ISXE)
          AJRAYO = ZZ .GT. DR
        ENDIF
        IF (AJRAYO) THEN
          IHC = IHC + 1
          ZZ = MIN(ZHR,  ZHE)
          ZZ = MAX(ZZ ,  ZHC)
          HXRAYO(IHC) = ZZ
          ZHC = ZZ
          IF (ZHR .GT. 0.) THEN
            ICC = ICC + 1
          ELSE
            ICC = ICC - 1
          ENDIF
          IXRAYO(IHC) = ICC
*
        ENDIF
      ENDIF
*
      HXRAYO(IHC+1) = MAX(ZHC, ZHE)
      NHMAX = IHC
      ISCE  = ISC
*
      RETURN
      END
