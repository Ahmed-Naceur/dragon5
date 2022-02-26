*DECK SYB4T2
      SUBROUTINE SYB4T2 (NIR,JMINRA,XCOTEA,SINPHI,DCOTEB,COSPHI,ORIGIN,
     & PENTES,ANGLES,IXRAYO,HXRAYO,DXMIN,DELTAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the position of the next track intersecting with side 2 (west
* side) in a sectorized Cartesian cell.
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
* NIR     number of radius.
* JMINRA  index of the first radius intercepting the side.
* XCOTEA  radius interceptions.
* SINPHI  sinus
* DCOTEB  half-side of the rectangle (center-side position).
* COSPHI  cosinus.
* ORIGIN  origin of the track.
* PENTES  slope of the track.
* ANGLES  slope of the sector.
* IXRAYO  index of the radius.
* HXRAYO  central radius (generally a negative number).
* DXMIN   accuracy.
*
*Parameters: output
* DELTAC  position of the next trajectory.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NIR,JMINRA,IXRAYO
      REAL    XCOTEA(JMINRA:NIR),SINPHI,DCOTEB,COSPHI,ORIGIN,PENTES,
     &        ANGLES,HXRAYO,DXMIN,DELTAC
*
* Rayon Exterieur Suivant
      IF (JMINRA .LE. NIR) THEN
      H0 = - DCOTEB * SINPHI
      IF ((IXRAYO .EQ. JMINRA) .OR. (HXRAYO .GE. H0)) THEN
        IR = IXRAYO
        S  = 1.
      ELSE
        IR = IXRAYO - 1
        S  = - 1.
      ENDIF
      IF (IR .LE. NIR) THEN
        X = S * XCOTEA(IR)
        D = DCOTEB * COSPHI + X * SINPHI
        DELTAC = MIN(DELTAC, D)
      ENDIF
      ENDIF
*
* Secteur Suivant
      DA = ANGLES - PENTES
      D  = ABS(ANGLES) + ABS(PENTES)
      IF (ABS(DA) .GT. DXMIN*D) THEN
        D = ORIGIN / DA
        DELTAC = MIN(DELTAC, D)
      ENDIF
*
      RETURN
      END
