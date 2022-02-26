*DECK SYB4T4
      SUBROUTINE SYB4T4 (NIR,JMINRA,XCOTEA,SINPHI,DCOTEB,COSPHI,ORIGIN,
     & PENTES,ANGLES,IXRAYO,DELTAH,DXMIN,DELTAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the position of the next track intersecting with side 4 in a
* sectorized Cartesian cell.
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
* SINPHI  sinus.
* DCOTEB  half-side of the rectangle (center-side position).
* COSPHI  cosinus.
* ORIGIN  origin of the track.
* PENTES  slope of the track.
* ANGLES  slope of the sector.
* IXRAYO  index of the radius.
* DELTAH  position of the corner.
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
     &        ANGLES,DELTAH,DXMIN,DELTAC
*
* Rayon Exterieur Suivant
      IF (JMINRA .LE. NIR) THEN
      IF (IXRAYO .LE. NIR) THEN
        Y = XCOTEA(IXRAYO)
        D = Y * SINPHI - DCOTEB * COSPHI
        DELTAC = MIN(DELTAC, D)
      ENDIF
      ENDIF
*
* Secteur Suivant
      DA = ANGLES - PENTES
      IF (ABS(DA) .GT. DXMIN) THEN
        D = ORIGIN / DA
        DELTAC = MIN(DELTAC, D)
      ENDIF
*
* Coin Sud-West
      DELTAC = MIN(DELTAC, DELTAH)
*
      RETURN
      END
