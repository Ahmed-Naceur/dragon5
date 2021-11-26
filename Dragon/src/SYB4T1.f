*DECK SYB4T1
      SUBROUTINE SYB4T1 (NIR,JMINRB,XCOTEB,COSPHI,Y0,ORIGIN,PENTES,
     & ANGLES,IXRAYO,DELTAH,DXMIN,DELTAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the position of the next track intersecting with side 1 in a
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
* COSPHI  cosinus.
* Y0      position of the tube center.
* ORIGIN  origin of the track.
* PENTES  slope of the track.
* ANGLES  slope of the sector.
* IXRAYO  index of the radius.
* DELTAH  position of the south-west corner.
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
      REAL    XCOTEB(JMINRB:NIR),COSPHI,Y0,ORIGIN,PENTES,ANGLES,DELTAH,
     &        DXMIN,DELTAC
*
* Rayon Exterieur Suivant
      IF (JMINRB .LE. NIR) THEN
      IF (IXRAYO .LE. NIR) THEN
        Y = XCOTEB(IXRAYO)
        D = Y * COSPHI - Y0
        DELTAC = MIN(DELTAC, D)
      ENDIF
      ENDIF
*
* Secteur Suivant
      DA = ANGLES - PENTES
      IF (ABS (DA) .GT. DXMIN) THEN
        D = ORIGIN / DA
        DELTAC = MIN(DELTAC, D)
      ENDIF
*
* Coin Sud-West
      DELTAC = MIN(DELTAC, DELTAH)
*
      RETURN
      END
