*DECK SYB4TR
      SUBROUTINE SYB4TR(MNA,NRD,NSECT4,COTEA,COTEB,RAYONS,IFAC,NUMREG,
     1 JMINRA,XCOTEA,JMINRB,XCOTEB,LFAIRE,DXMIN,DELR,IQW,PWA2,ZWA2,
     2 NXMIN,NXMAX,NZR,ZZR,NZI,ZZI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the tracking information related to a square or rectangular
* sectorized heterogeneous cell (called 4 times).
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
* MNA     number of angles in (0,$\\pi$/2).
* NRD     one plus the number of tubes in the cell.
* NSECT4  number of sectors.
* COTEA   X-axis Cartesian dimensions of the cell.
* COTEB   Y-axis Cartesian dimensions of the cell.
* RAYONS  radius of each cylinder.
* IFAC    starting side (=0, 1, 2 or 3).
* NUMREG  merged volume number in each sector.
* JMINRA  first interception with side a.
* XCOTEA  interceptions with side a.
* JMINRB  first interception with side b.
* XCOTEB  interceptions with side b.
* LFAIRE  tracking calculation flag (=.FALSE. only compute the number
*         of real and integer tracking elements).
* DXMIN   geometrical epsilon.
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
      INTEGER  MNA,NRD,NSECT4,IFAC,NUMREG(NSECT4,NRD),JMINRA,JMINRB,
     1         IQW,NXMIN,NXMAX,NZR,NZI,ZZI(*)
      REAL     COTEA,COTEB,RAYONS(NRD-1),XCOTEA(JMINRA:NRD-1),
     1         XCOTEB(JMINRB:NRD-1),DXMIN,DELR,PWA2(MNA),ZWA2(MNA),
     2         ZZR(*)
      LOGICAL  LFAIRE
*----
*  LOCAL VARIABLES
*----
*ISCW   No du Premier Secteur A L'Ouest (a gauche)
*ISXW   No de la Surface Externe A L'Ouest (a gauche)
*ISXE   No de la Surface Externe A L'Est (a droite)
      PARAMETER     (PI314  = 3.141592653589793)
      PARAMETER     (PIS2   = PI314 /  2)
      PARAMETER     (PIS4   = PI314 /  4)
      LOGICAL        NTRIOR
      REAL           ORIGIN(4)
      REAL           PENTES(4)
      REAL           WX(64)
      REAL           ZX(64)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IXRAYO
      REAL, ALLOCATABLE, DIMENSION(:) :: HXRAYO,ANGLES,COSECT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IXRAYO(NRD*2+NSECT4+3))
      ALLOCATE(HXRAYO(NRD*2+NSECT4+3),ANGLES(NSECT4),COSECT(NSECT4))
*
* Interpolation des Trajectoires
      DCOTEA = COTEA * 0.5
      DCOTEB = COTEB * 0.5
      IZI = 0
      IZR = 0
*
      DELR2  = DELR * 2.
*
* /Debut/ Boucle sur les Angles (0 a Pi/2)
      DO 350 IA = 1, MNA
      MNT  = 0
      IZI  = IZI + 2
      IZIT = IZI
*
      WANGIA = PWA2(IA) * 0.5
      ZAIA   = ZWA2(IA) + 1.
      YAIA   = 1. - ZWA2(IA)
*
      PHI  = PIS4 * ZAIA
      COSPHI = COS(PHI)
      SINPHI = SIN(PHI)
      PENTES(1) = - SINPHI / COSPHI
      ORIGIN(1) = - DCOTEA / COSPHI
      PENTES(2) =   COSPHI / SINPHI
      ORIGIN(2) = - DCOTEB / SINPHI
      PENTES(3) = - SINPHI / COSPHI
      ORIGIN(3) =   DCOTEA / COSPHI
      PENTES(4) =   COSPHI / SINPHI
      ORIGIN(4) =   DCOTEB / SINPHI
      IF (LFAIRE) THEN
        ZZR(IZR+1) = SINPHI
        ZZR(IZR+3) = COSPHI
        IF ((IFAC .EQ. 0) .OR. (IFAC .EQ. 2)) THEN
          ZZR(IZR+2) = COSPHI
          ZZR(IZR+4) = - SINPHI
        ELSE
          ZZR(IZR+2) = - COSPHI
          ZZR(IZR+4) = SINPHI
        ENDIF
        DO 10 I=1,4
        ZZR(IZR+4+I) = ZZR(IZR+I)
   10   CONTINUE
        ZZR(IZR+9) = WANGIA
      ENDIF
        IZR = IZR + 9
      PHIJ = PHI + PI314 / 2
      PHIJ = - PHIJ
      PHIK = PI314 * 2 / NSECT4
      ISCMIN = 1
      ISCMAX = NSECT4
      DO I = 1, NSECT4
        PHIJ = PHIJ + PHIK
        IF     (PHIJ .LE. (DXMIN - PIS2)) THEN
          ISCMIN = I + 1
        ELSEIF (PHIJ .GE. (PIS2 - DXMIN)) THEN
          ISCMAX = I - 1
          GOTO 20
        ELSE
        ANGLES(I) = TAN(PHIJ)
        COSECT(I) = COS(PHIJ)
        ENDIF
      ENDDO
*
   20 DELTAH = (DCOTEB * COSPHI) - (DCOTEA * SINPHI)
*
* Rayon Max Pour cet angle (Rext)
      DELTAS = DCOTEA * SINPHI + DCOTEB * COSPHI
*
* Suivant West (Gauche,Bas=Sud)
      DELTAC = 0.
*----
*  /Debut/ Boucle sur les Trajectoires
*----
      NTRIOR = .TRUE.
      DO WHILE (NTRIOR)
        MNT = MNT + 1
        DELTAP = DELTAC
        CALL SYB4TH (NRD-1,  RAYONS,
     &               PENTES, ORIGIN,
     &               ISCMIN, ISCMAX, ANGLES,
     &               DXMIN,  DELTAC,
     &               ISCW,   ISCE,   ISXW,   ISXE,   NHMAX,
     &               IXRAYO, HXRAYO)
*
        DELTAC = DELTAS
*
        CALL SYB4TN (NHMAX,  IXRAYO, ISCW,
     &               COSECT, NRD-1, RAYONS,
     &               DELTAC)
*
        IF (ISXW .EQ.1) THEN
          CALL SYB4T1 (NRD-1, JMINRB, XCOTEB,
     &                 COSPHI, - DCOTEA * SINPHI,
     &                 ORIGIN(1), PENTES(1), ANGLES(ISCW),
     &                 IXRAYO, DELTAH,
     &                 DXMIN,  DELTAC)
        ELSE
          CALL SYB4T2 (NRD-1, JMINRA, XCOTEA,
     &                 SINPHI, DCOTEB, COSPHI,
     &                 ORIGIN(2), PENTES(2), ANGLES(ISCW),
     &                 IXRAYO, HXRAYO,
     &                 DXMIN,  DELTAC)
        ENDIF
*
        IF (ISXE .EQ. 4) THEN
          CALL SYB4T4 (NRD-1, JMINRA, XCOTEA,
     &                 SINPHI, DCOTEB, COSPHI,
     &                 ORIGIN(4), PENTES(4), ANGLES(ISCE-1),
     &                 IXRAYO(NHMAX), - DELTAH,
     &                 DXMIN,  DELTAC)
        ELSE
          CALL SYB4T3 (NRD-1, JMINRB, XCOTEB,
     &                 SINPHI, DCOTEA, COSPHI,
     &                 ORIGIN(3), PENTES(3), ANGLES(ISCE-1),
     &                 IXRAYO(NHMAX), HXRAYO(NHMAX),
     &                 DXMIN,  DELTAC)
        ENDIF
*
* Intervalle entre 2 Intersections, Decoupage
      DELTAX = DELTAC - DELTAP
      IF (DELTAX .LE. DELR2) THEN
        NX = 1
        ZX(1) = 0.0
        WX(1) = 2.0
      ELSE
        NX  = INT(DELTAX / DELR2 + 1)
        IF(IQW.EQ.0) THEN
*          GAUSS-LEGENDRE INTEGRATION POINTS. ZX(I) IS NOT USED.
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
                 CALL XABORT('SYB4TR: GAUSS OVERFLOW.')
              ENDIF
           ENDIF
           CALL ALGPT(NX,-1.0,1.0,ZX,WX)
        ELSE
*          EQUAL WEIGHT INTEGRATION POINTS.
           DO 30 I=1,NX
           ZX(I)=(2.0*REAL(I)-1.0)/REAL(NX)-1.0
           WX(I)=2.0/REAL(NX)
   30      CONTINUE
        ENDIF
      ENDIF
      NXMIN = MIN(NX, NXMIN)
      NXMAX = MAX(NX, NXMAX)
* Intervalle entre 2 Intersections, Fin du Decoupage
*
      IF (LFAIRE) THEN
*
* Debut : Nombre et Numeros des Regions
      CALL SYB4TI  (NHMAX, IXRAYO, ISCW,
     &              NSECT4, IFAC, NUMREG,
     &              NLMAX, ZZI(IZI+4))
      ZZI(IZI+1) = NLMAX
      ZZI(IZI+2) = NX
      ZZI(IZI+3) = ISXW
      IF (IFAC .EQ. 1) THEN
        ZZI(IZI+3) = 6 - ZZI(IZI+3)
      ELSEIF (IFAC .EQ. 2) THEN
        ZZI(IZI+3) = ZZI(IZI+3) + 2
      ELSEIF (IFAC .EQ. 3) THEN
        ZZI(IZI+3) = 4 - ZZI(IZI+3)
      ENDIF
      IF (ZZI(IZI+3) .GE. 4) ZZI(IZI+3) = ZZI(IZI+3) - 4
*
      IZI = IZI + 4 + NLMAX
      ZZI(IZI) = ISXE
      IF (IFAC .EQ. 1) THEN
        ZZI(IZI) = 6 - ZZI(IZI)
      ELSEIF (IFAC .EQ. 2) THEN
        ZZI(IZI) = ZZI(IZI) + 2
      ELSEIF (IFAC .EQ. 3) THEN
        ZZI(IZI) = 4 - ZZI(IZI)
      ENDIF
      IF (ZZI(IZI) .GE. 4) ZZI(IZI) = ZZI(IZI) - 4
* Fin : Nombre et Numeros des Regions
*
*----
* Boucle sur les trajectoires
*----
      DELTAR = DELTAP
      DO 250 IX = 1, NX
*
        WW = 0.5 *  WX(IX)
*
        IZR = IZR + 1
        ZZR(IZR) = WW * WANGIA * DELTAX
        DDELTA = DELTAX * WW
        DELTAR = DELTAR + DDELTA
*
* Position de L'intersection Gauche (West)
        ZZW = DELTAR * PENTES(ISXW) + ORIGIN(ISXW)
*
* Position de L'intersection Droite (Est)
        ZZE = DELTAR * PENTES(ISXE) + ORIGIN(ISXE)
*
* Longueur des intersections
        CALL SYB4TC (DELTAR, DDELTA, ANGLES,
     &               NHMAX, IXRAYO,
     &               ISCW,  NSECT4, IFAC, NUMREG,
     &               RAYONS, ZZW, ZZE,
     &               ZZR(IZR+1), HXRAYO)
        IZR = IZR + NLMAX
*
  250 CONTINUE
      ELSE
*
* Pour le comptage des Numeros de Regions
* deux indices de plus pour les surfaces Entrante et Sortante
        IZR = IZR + (NHMAX + 1) * NX
        IZI = IZI + 5 + NHMAX
      ENDIF
*
* Position des limites Gauche et Droite (West et Est)
      DR = DELTAC + DXMIN
      ZZW = DR * PENTES(ISXW) + ORIGIN(ISXW)
      ZZE = DR * PENTES(ISXE) + ORIGIN(ISXE)
      NTRIOR = ZZW .LT. ZZE
      IF (NTRIOR) THEN
        IF (DELTAC .LE. DELTAP) THEN
          CALL XABORT ('SYB4TR: INFINITE LOOP.')
        ENDIF
      ENDIF
*
* /Fin/ Boucle sur les Trajectoires
      ENDDO
*
* /Fin/ Boucle sur les Angles
      IF (LFAIRE) THEN
        ZZI(IZIT-1) = MNT
        ZZI(IZIT) = IFAC
      ENDIF
  350 CONTINUE
      NZI = IZI
      NZR = IZR
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(COSECT,ANGLES,HXRAYO)
      DEALLOCATE(IXRAYO)
*
      RETURN
      END
