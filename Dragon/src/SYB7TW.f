*DECK SYB7TW
      SUBROUTINE SYB7TW (NRD,MRE,XICOTE,IRMIN,DELTCW,DELTCS,ORIPH6,
     & ORIPHI,COSPH6,SINPHI,ISW2,LGTRAW,DELTAW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the intersection lenghts of a track in a sectorized hexagonal
* cell, following a singular point at west.
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
* NRD     one plus the number of tubes in the cell.
* MRE     first radius intersection the side or first region in contact
*         with the external side.
* XICOTE  intersections with the side.
* IRMIN   region index in contact with the external side.
* DELTCW  distance from the south-west corner.
* DELTCS  distance from the south-east corner.
* ORIPH6  distance at the middle of south-west side.
* ORIPHI  distance at the middle of south side.
* COSPH6  cosinus ($\\pi$/6-Phi).
* SINPHI  sinus (Phi).
*
*Parameters: input/output
* ISW2    half side on west (2-3 followed by 4-5 or SW (N-S) followed by
*         south (W-E)).
*
*Parameters: output
* LGTRAW  removing/addition flag of the next intersection. The formula
*         corresponding to case ISW2=5 is complex. In this case, the
*         tangent can pass on both sides.
* DELTAW  distance of the next track.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NRD,MRE,IRMIN,ISW2
      REAL    XICOTE(NRD-1),DELTCW,DELTCS,ORIPH6,ORIPHI,COSPH6,SINPHI,
     &        DELTAW
      LOGICAL LGTRAW
*
** Le cas le plus courant est : Secteur Suivant (Sud)
** Le cas general est nettement plus complexe :
*
* NE PAS CONFONDRE
*  Intervalle et Limite
*
* Le Suivant depasse le Milieu du Cote Sud-West
      IF (ISW2 .EQ. 2) THEN
        IF (IRMIN .LE. MRE) THEN
          ISW2 = 3
        ENDIF
      ENDIF
*
* Le Suivant depasse le Milieu du Cote Sud
      IF (ISW2 .EQ. 4) THEN
        IF (IRMIN .LE. MRE) THEN
          ISW2 = 5
        ENDIF
      ENDIF
*
* Limite suivante Sud-West Nord (Haut)
      IF (ISW2 .EQ. 2) THEN
        LGTRAW = .TRUE.
      DELTAW = ORIPH6 - XICOTE(IRMIN-1) * COSPH6
*
* Limite suivante de Sud-West Sud (Bas)
      ELSEIF (ISW2 .EQ. 3) THEN
*
          LGTRAW = IRMIN .EQ. NRD
*
* Limite suivante (secteur) coin SSW
          IF (LGTRAW) THEN
      DELTAW = DELTCW
*
* Limite suivante egalement Sud-West Sud (Bas)
          ELSE
      DELTAW = ORIPH6 + XICOTE(IRMIN) * COSPH6
          ENDIF
*
* Limite suivante de Sud West (Horizontal Gauche)
      ELSEIF (ISW2 .EQ. 4) THEN
        LGTRAW = .TRUE.
      DELTAW = ORIPHI - XICOTE(IRMIN-1) * SINPHI
*
* Limite suivante de Sud West (Horizontal Gauche)
      ELSE
*
* Limite suivante (secteur) coin Extreme SSE
        IF (IRMIN .GE. NRD) THEN
          LGTRAW = .FALSE.
      DELTAW = DELTCS
        ELSE
          LGTRAW = .TRUE.
      DELTAW = ORIPHI + XICOTE(IRMIN) * SINPHI
        ENDIF
      ENDIF
*
      RETURN
      END
