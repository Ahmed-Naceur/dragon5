*DECK XELPSC
      FUNCTION XELPSC(RANN,PLANE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute annular surface below Cartesian plane.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):G. Marleau
*
*Parameters: input
* RANN    annular radius.
* PLANE   Cartesian plane location.
*
*Parameters: output
* XELPSC  annular surface below plane.
*
*-----------------------------------------------------------------------
*
      IMPLICIT          NONE
      DOUBLE PRECISION  XELPSC,RANN,PLANE
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  RANN2,PLANE2,ALPHA
*----
*  METHOD
*  1) FIND HALF-ANGLE COVERED BY THE TWO INTERSECTION POINTS
*     BETWEEN PLANE AND ANNULAR REGION
*     -- ALPHA=ACOS(-PLANE/RANN)
*  2) COMPUTED ANNULAR SURFACE COVERED BY THIS HALF-ANGLE
*     -- 0.5*RANN2*ALPHA
*  3) ADD SURFACE COVERED BY INTERNAL RECTANGLE IN THIS HALF-ANGLE
*     -- 0.5*PLANE*SQRT(RANN2-PLANE2)
*  4) DOULBLE SURFACE FOR FULL ANGLE
*----
      RANN2=RANN*RANN
      PLANE2=PLANE*PLANE
      ALPHA=ACOS(-PLANE/RANN)
      XELPSC=RANN2*ALPHA+PLANE*SQRT(RANN2-PLANE2)
      RETURN
      END
