*DECK SYB7T0
      SUBROUTINE SYB7T0(MNA,NRD,COTE,RAYONS,JMINR,XCOTE,LFAIRE,DELR,
     1 IQW,PWA2,ZWA2,NXMIN,NXMAX,NZR,ZZR,NZI,ZZI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the tracking information related to an hexagonal sectorized
* heterogeneous cell.
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
* MNA     number of angles in (0,$\\pi$/6).
* NRD     one plus the number of tubes in the cell.
* COTE    length of the hexagon side.
* RAYONS  radius of each cylinder.
* JMINR   first interception with side.
* XCOTE   interceptions with side.
* LFAIRE  tracking calculation flag (=.FALSE. only compute the number
*         of real and integer tracking elements).
* DELR    half distance between the tracks.
* IQW     equal weight quadrature flag (=1 to use equal weight
*         quadratures in angle and space).
* PWA2    weights of the angular quadrature set.
* ZWA2    base points of the angular quadrature set.
*
*Parameters: output
* NXMIN   minimum number of tracks per region.
* NXMAX   maximum number of tracks per region.
* NZR     number of real tracking elements.
* ZZR     real tracking information.
* NZI     number of integer tracking elements.
* ZZI     integer tracking information.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER        MNA,NRD,JMINR,IQW,NXMIN,NXMAX,NZR,NZI,ZZI(*)
      REAL           COTE,RAYONS(NRD-1),XCOTE(NRD),DELR,PWA2(64),
     &               ZWA2(64),ZZR(*)
      LOGICAL        LFAIRE
*----
*  LOCAL VARIABLES
*----
      PARAMETER    (PI314  = 3.141592653589793)
      PARAMETER    (PI6    = PI314 /  6)
      PARAMETER    (PI12   = PI314 / 12)
      PARAMETER    (SQRT3  = 1.732050807568877)
      PARAMETER    (SQRT32 = SQRT3 / 2)
      REAL           ANGLES(5)
      REAL           COSECT(3)
      LOGICAL        LGTRAE
      LOGICAL        LGTRAW
      LOGICAL        LVERIF
      LOGICAL        NTRIOR
      REAL           ORIGIN(4)
      REAL           PENTES(4)
      CHARACTER *  4 TYSUIT
      REAL           WX(64)
      REAL           ZX(64)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IXRAYO
      REAL, ALLOCATABLE, DIMENSION(:) :: HXRAYO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IXRAYO(NRD*2+3),HXRAYO(NRD*2+3))
*----
*  Interpolation des Trajectoires
*----
      IZI = 0
      IZR = 0
      LVERIF = .FALSE.
      DELR2  = DELR * 2.
      HAUTEU = COTE * SQRT32
*----
*  /Debut/ Boucle sur les Angles
*----
      DO 350 IA = 1, MNA
      MNT  = 0
      IZI  = IZI + 1
      IZIT = IZI
*
      WANGIA = PWA2(IA)
      ZAIA   = ZWA2(IA) + 1.
      YAIA   = 1. - ZWA2(IA)
      WANGIA = WANGIA / 6.
*
      PHI  = PI12 * ZAIA
      PHI6 = PI12 * YAIA
      PHI3 = PHI + PI6
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
      COSPH3 = COS(PHI3)
      SINPH3 = SIN(PHI3)
      COSPH6 = COS(PHI6)
      SINPH6 = SIN(PHI6)
      TANPH6 = TAN(PHI6)
      IF(LFAIRE) THEN
        ZZR(IZR+1) = COSPH6
        ZZR(IZR+2) = SINPH6
        ZZR(IZR+3) = SINPHI
        ZZR(IZR+4) = COSPHI
        ZZR(IZR+5) = COSPH3
        ZZR(IZR+6) = -SINPH3
        ZZR(IZR+7) = WANGIA
      ENDIF
      IZR = IZR + 7
      ANGLES(1) = COSPHI
      ANGLES(2) = SINPHI
      ANGLES(3) = - TAN(PHI3)
      ANGLES(4) =   TANPH6
      ANGLES(5) = COSPHI / SINPHI
*
      PENTES(1) = TANPH6
      PENTES(2) = ANGLES(5)
      PENTES(3) = - TAN(PHI3)
      PENTES(4) = PENTES(1)
      ORIGIN(1) = - HAUTEU / COS(PHI6)
      ORIGIN(2) = - HAUTEU / SINPHI
      ORIGIN(3) =   HAUTEU / COS(PHI3)
      ORIGIN(4) = - ORIGIN(1)
*
      COSECT(1) = COS(PHI3)
      COSECT(2) = COS(PHI6)
      COSECT(3) = SINPHI
*----
*  Rayon Max Pour cet angle(Rext)
*----
      XIZERO = HAUTEU * TANPH6
      DO MRA = JMINR, NRD-1
        IF(XCOTE(MRA) .GE. XIZERO)  GOTO 140
      ENDDO
      MRA = NRD
  140 MRAE = MRA
      MRAW = MRA
*----
*  Les Lignes d'integrations sont limites
*  et par les secteurs
*  et par les couronnes
*
*  Demarrage au centre de l'hexagone
*----
      IHMAX  = NRD - MRA
      IHMIN  = IHMAX + 1
      HXRAYO(IHMIN) = ORIGIN(1)
      DO IR = MRAW, 2, -1
        IHMAX = IHMAX + 1
        IXRAYO(IHMAX) = IR
        HXRAYO(IHMAX+1) = - RAYONS(IR-1)
      ENDDO
      DO IS = 1, 3
        IXRAYO(IHMAX+IS) = 1
        HXRAYO(IHMAX+IS+1) = 0.
      ENDDO
      IHMAX = IHMAX + 3
      HXRAYO(IHMAX+1) = 0.
      DO IR = 1, MRAE-1
        IHMAX = IHMAX + 1
        IXRAYO(IHMAX) = IR
        HXRAYO(IHMAX+1) = RAYONS(IR)
      ENDDO
      IHMAX = IHMAX + 1
      IXRAYO(IHMAX) = MRAE
      HXRAYO(IHMAX+1) = ORIGIN(4)
*
      DELTCW = COTE * COSECT(1)
      DELTAS = COTE * COSECT(2)
      DELTCE = COTE * COSECT(3)
      ORIPHI = HAUTEU * COSPHI
      ORIPH6 = HAUTEU * SINPH6
      ORIPH3 = HAUTEU * SINPH3
*
      ISW2 = 2
      CALL SYB7TW(NRD,    JMINR,    XCOTE, IXRAYO(IHMIN),
     &             DELTCW, DELTAS, ORIPH6, ORIPHI,
     &             COSPH6, SINPHI,
     &             ISW2,   LGTRAW, DELTAW)
      ISW = ISW2 / 2
*
      ISE2   = 8
      CALL SYB7TE(NRD,    JMINR,    XCOTE, IXRAYO(IHMAX),
     &             DELTCE, DELTAS, ORIPH6, ORIPH3,
     &             COSPH6, COSPH3,
     &             ISE2,   LGTRAE, DELTAE)
*
      DELTAC = 0.
      NTRIOR = .TRUE.
      IRSUIT = 0
*----
*  /Debut/ Boucle sur les Trajectoires
*----
      DO WHILE(NTRIOR)
        NHMAX = IHMAX + 1 - IHMIN
        DELTAP = DELTAC
*
        TYSUIT = 'Sud'
        IHSUIT = IHMIN
        DELTAC = DELTAS
        IF(DELTAW .LT. DELTAC) THEN
          TYSUIT = 'West'
          IHSUIT = IHMIN
          DELTAC = DELTAW
        ENDIF
        IF(DELTAE .LT. DELTAC) THEN
          TYSUIT = 'Est'
          IHSUIT = IHMAX
          DELTAC = DELTAE
        ENDIF
        CALL SYB7TN(IHMIN,  IHMAX,  IXRAYO, ISW,
     &               COSECT, NRD-1, RAYONS,
     &               TYSUIT, IHSUIT, DELTAC, IRSUIT)
*
        DELTAX = DELTAC - DELTAP
        IF(DELTAX .LE. DELR2) THEN
          NX = 1
          ZX(1) = 0.0
          WX(1) = 2.0
        ELSE
          NX  = INT(DELTAX / DELR2 + 1)
          IF(IQW.EQ.0) THEN
*           GAUSS-LEGENDRE INTEGRATION POINTS. ZX(I) IS NOT USED.
            IF(NX.GT.20) THEN
              IF(NX.LT.24) THEN
                 NX=24
              ELSE IF(NX.LT.28) THEN
                 NX=28
              ELSE IF(NX.LT.32) THEN
                 NX=32
              ELSE IF(NX.LT.64) THEN
                 NX=64
              ELSE IF(NX.GT.64) THEN
                 CALL XABORT('SYB7T0: GAUSS OVERFLOW.')
              ENDIF
            ENDIF
            CALL ALGPT(NX,-1.0,1.0,ZX,WX)
          ELSE
*           EQUAL WEIGHT INTEGRATION POINTS.
            DO 30 I=1,NX
            ZX(I)=(2.0*REAL(I)-1.0)/REAL(NX)-1.0
            WX(I)=2.0/REAL(NX)
   30       CONTINUE
          ENDIF
        ENDIF
        NXMIN = MIN(NX, NXMIN)
        NXMAX = MAX(NX, NXMAX)
*
        MNT = MNT + 1
        IF(LFAIRE) THEN
          ZZI(IZI+1) = NHMAX
          ZZI(IZI+2) = NX
          DO I=0,NHMAX-1
            ZZI(IZI+3+I)=IXRAYO(IHMIN+I)
          ENDDO
*
          DELTAR = DELTAP
          DO 250 IX = 1, NX
*
          WW = 0.5 * WX(IX)
*
          IZR = IZR + 1
          ZZR(IZR) = WW * WANGIA * DELTAX
          DDELTA = DELTAX * WW
          DELTAR = DELTAR + DDELTA
*
* Position de L'intersection Gauche(West)
          ZZW = DELTAR * PENTES(ISW) + ORIGIN(ISW)
*
* Position de L'intersection Droite(Est)
          ISE = ISE2 / 2
          ZZE = DELTAR * PENTES(ISE) + ORIGIN(ISE)
*----
*  Longueur des intersections
*----
          CALL SYB7TC(DELTAR, DDELTA, ANGLES(ISW+2), NHMAX, ZZI(IZI+3),
     &               NRD-1,  RAYONS, ZZW, ZZE,
     &               ZZR(IZR+1), HXRAYO(IHMIN))
          IZR = IZR + NHMAX
*
  250     CONTINUE
        ELSE
*----
* Pour le comptage il faudra compter
* pour plus tard (les douzes symetries)
* deux indices de plus pour les surfaces Entrante et Sortante
*----
          IZR = IZR +(NHMAX + 1) * NX
          IZI = IZI + 2
        ENDIF
*
        IZI = IZI + NHMAX + 2
*
        IF(TYSUIT .EQ. 'Est') THEN
          IRC   = IXRAYO(IHMAX)
          IF(LGTRAE) THEN
            IHMAX = IHMAX - 1
            IF(IRC .EQ. IXRAYO(IHMAX)) ISE2 = ISE2 - 1
          ELSEIF(IRC .GE. NRD) THEN
            NTRIOR = .FALSE.
          ELSE
            IHMAX = IHMAX + 1
            IXRAYO(IHMAX) = IRC + 1
            HXRAYO(IHMAX+1) = HXRAYO(IHMAX)
          ENDIF
          IF(NTRIOR) THEN
            CALL SYB7TE(NRD,    JMINR,    XCOTE, IXRAYO(IHMAX),
     &               DELTCE, DELTAS, ORIPH6, ORIPH3,
     &               COSPH6, COSPH3,
     &               ISE2,   LGTRAE, DELTAE)
          ENDIF
        ELSEIF(TYSUIT .EQ. 'West') THEN
          IRC   = IXRAYO(IHMIN)
          IF(LGTRAW) THEN
            IHMIN = IHMIN + 1
            IF(IRC .EQ. IXRAYO(IHMIN)) ISW2 = ISW2 + 1
          ELSEIF(IRC .GE. NRD) THEN
            NTRIOR = .FALSE.
          ELSE
            IHMIN = IHMIN - 1
            IXRAYO(IHMIN) = IRC + 1
            HXRAYO(IHMIN) = HXRAYO(IHMIN+1)
          ENDIF
          NTRIOR = ISW .LE. 2
          IF(NTRIOR) THEN
            CALL SYB7TW(NRD,    JMINR,    XCOTE, IXRAYO(IHMIN),
     &                 DELTCW, DELTAS, ORIPH6, ORIPHI,
     &                 COSPH6, SINPHI,
     &                 ISW2,   LGTRAW, DELTAW)
            IF((ISW2 / 2) .NE. ISW) THEN
              IF(LFAIRE) THEN
                ZZI(IZI+1) = 0
                ZZI(IZI+2) = 0
                MNT = MNT + 1
              ENDIF
              IZI = IZI+2
            ENDIF
            ISW = ISW2 / 2
            IF(ISW2 .EQ. 5) THEN
              IF(LGTRAW) THEN
                HC = DELTAW * PENTES(ISW) + ORIGIN(ISW)
                LGTRAW = HC .GE. 0.
              ENDIF
            ENDIF
          ENDIF
        ELSEIF(TYSUIT .EQ. 'Coin') THEN
          IXRAYO(IHSUIT) = IRSUIT
        ELSEIF(TYSUIT .EQ. 'Tang') THEN
*
**        Decalage du tableau entier de 2 Cases
          DO I = 2, IHMAX - IHSUIT
            IXRAYO(IHSUIT+I-2) = IXRAYO(IHSUIT+I)
          ENDDO
*
**        Decalage du tableau reel de 2 Cases
          DO I = 2, IHMAX + 1 - IHSUIT
            HXRAYO(IHSUIT+I-2) = HXRAYO(IHSUIT+I)
          ENDDO
*
          IHMAX = IHMAX - 2
        ELSE
          NTRIOR = .FALSE.
        ENDIF
*
* /Fin/ Boucle sur les Trajectoires
      ENDDO
*
* /Fin/ Boucle sur les Angles
      IF(LFAIRE)  ZZI(IZIT) = MNT
  350 CONTINUE
      NZI = IZI
      NZR = IZR
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(HXRAYO,IXRAYO)
*
      RETURN
      END
