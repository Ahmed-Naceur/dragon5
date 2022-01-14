*DECK SYB4T3
      SUBROUTINE SYB4T3 (NIR,JMINRB,XCOTEB,SINPHI,DCOTEA,COSPHI,ORIGIN,
     & PENTES,ANGLES,IXRAYO,HXRAYO,DXMIN,DELTAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the position of the next track intersecting with side 3 in a
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
* JMINRB  index of the first radius intercepting the side.
* XCOTEB  radius interceptions.
* SINPHI  sinus.
* DCOTEA  half-side of the rectangle (center-side position).
* COSPHI  cosinus.
* ORIGIN  origin of the track.
* PENTES  slope of the track.
* ANGLES  slope of the sector.
* IXRAYO  index of the radius.
* HXRAYO  central radius.
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
      INTEGER NIR,JMINRB,IXRAYO
      REAL    XCOTEB(JMINRB:NIR),SINPHI,DCOTEA,COSPHI,ORIGIN,PENTES,
     &        ANGLES,HXRAYO,DXMIN,DELTAC
*
* Rayon Exterieur Suivant
      IF (JMINRB .LE. NIR) THEN
      H0 = DCOTEA * COSPHI
      IF ((IXRAYO .EQ. JMINRB) .OR. (HXRAYO .LE. H0)) THEN
        IR = IXRAYO
        S  = 1.
      ELSE
        IR = IXRAYO - 1
        S  = - 1.
      ENDIF
      IF (IR .LE. NIR) THEN
        X = S * XCOTEB (IR)
        D = DCOTEA * SINPHI + X * COSPHI
        DELTAC = MIN (DELTAC, D)
      ENDIF
      ENDIF
*
* Secteur Suivant
      DA = ANGLES - PENTES
      IF (ABS (DA) .GT. DXMIN) THEN
        D = ORIGIN / DA
        DELTAC = MIN (DELTAC, D)
      ENDIF
*
      RETURN
      END
