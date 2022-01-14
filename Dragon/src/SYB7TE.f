*DECK SYB7TE
      SUBROUTINE SYB7TE (NRD,MRE,XICOTE,IRMAX,DELTCE,DELTCS,ORIPH6,
     & ORIPH3,COSPH6,COSPH3,ISE2,LGTRAE,DELTAE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the intersection lenghts of a track in a sectorized hexagonal
* cell, following a singular point at east.
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
* NRD     one plus the number of tube in the cell.
* MRE     first radius intersection the side or first region in contact
*         with the external side.
* XICOTE  intersections with the side.
* IRMAX   region index in contact with the external side.
* DELTCE  distance from the east corner.
* DELTCS  distance from the south-east corner.
* ORIPH6  distance at the middle of north-east side.
* ORIPH3  distance at the middle of south-east side.
* COSPH6  cosinus ($\\pi$/6-Phi).
* COSPH3  cosinus ($\\pi$/6+Phi).
*
*Parameters: input/output
* ISE2    half side on east (8 followed by 7-6 or NE followed by
*         SE (north and south)).
*
*Parameters: output
* LGTRAE  removing/addition flag of the next intersection. A removing
*         occurs if and only if the track cross a sector and the
*         south-east north (upper) zone.
* DELTAE  distance of the next track.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NRD,MRE,IRMAX,ISE2
      REAL    XICOTE(NRD-1),DELTCE,DELTCS,ORIPH6,ORIPH3,COSPH6,COSPH3,
     &        DELTAE
      LOGICAL LGTRAE
*
        IF (ISE2 .EQ. 8) THEN
          LGTRAE =  (IRMAX .EQ. NRD)
          IF (LGTRAE) THEN
      DELTAE = DELTCE
          ELSE
      DELTAE = XICOTE(IRMAX) * COSPH6 - ORIPH6
          ENDIF
        ELSE
          LGTRAE = ISE2 .EQ. 7
          IF (LGTRAE) THEN
            LGTRAE = IRMAX .GT. MRE
          ENDIF
          IF (LGTRAE) THEN
      DELTAE = ORIPH3 - XICOTE(IRMAX-1) * COSPH3
          ELSE
            ISE2   = 6
            IF (IRMAX .LT. NRD) THEN
      DELTAE = ORIPH3 + XICOTE(IRMAX) * COSPH3
            ELSE
      DELTAE = DELTCS
            ENDIF
          ENDIF
        ENDIF
*
      RETURN
      END
