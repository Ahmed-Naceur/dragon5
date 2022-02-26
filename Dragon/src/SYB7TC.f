*DECK SYB7TC
      SUBROUTINE SYB7TC (DELTAR, DDELTA, ANGLES, NHMAX, IXRAYO,
     &                   NRI,    RAYONS, ZZW, ZZE, ZZR, HXRAYO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the intersection lenghts of a track in a sectorized hexagonal
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
* DELTAR  used to compute the intersection.
* DDELTA  used to compute the intersection.
* ANGLES  angular values (begin at 1 or 2):
*         1=   $\\tan$(-$\\pi$/6-PHI);
*         2=   $\\tan$(+$\\pi$/6-PHI);
*         3=   $\\tan$(3$\\pi$/6-PHI).
* NHMAX   number of intersections.
* IXRAYO  tube indices.
* NRI     number of radii (= NRD-1)
* RAYONS  radius of each cylinder.
* ZZW     position of the west intersection (left).
* ZZE     position of the east intersection (right).
*
*Parameters: output
* ZZR     intersection lenghts.
*
*Parameters: input/output
* HXRAYO  preceding/next intersection lenghts.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NHMAX,IXRAYO(NHMAX),NRI
      REAL    DELTAR,DDELTA,ANGLES(3),RAYONS(NRI),ZZW,ZZE,ZZR(NHMAX),
     &        HXRAYO(NHMAX+1)
*
      DDELT2 = DDELTA * DDELTA
      DELTA2 = DELTAR * DELTAR
      GC  = 0.
      HC  = ZZW
      IRC = IXRAYO(1)
      ISC = 0
      HPRAYO = HXRAYO(1)
      HXRAYO(1) = HC
*
      DO IH = 1, NHMAX-1
        GP  = GC
        GC  = 0.
        HP  = HC
        IRP = IRC
        IRR = 0
*
        IRC = IXRAYO(IH+1)
        IF (IRP .EQ. IRC) THEN
          ISC = ISC + 1
          HC = DELTAR * ANGLES(ISC)
        ELSE
          IRR = MIN(IRP, IRC)
*
* Distance
          H2 = RAYONS(IRR) * RAYONS(IRR) - DELTA2
          IF (H2 .GT. 0) THEN
            HC = SQRT(H2)
            IF (IRC .EQ. IRR) THEN
              HC = - HC
            ENDIF
          ELSE
            HC = 0.
          ENDIF
*
        ENDIF
*
* Protection contre les longueurs negatives
        IF (HC .LT. HP) THEN
          HC = HP
          ZZH = 0.
        ELSE
          ZZH = HC - HP
        ENDIF
        ZZH = (ZZH + HXRAYO(IH+1) - HPRAYO) * 0.5
        HPRAYO = HXRAYO(IH+1)
        HXRAYO(IH+1) = HC
*
* Ajout de la courbure
        IF (IRP .NE. IRC) THEN
          H2 = HC - HPRAYO
          H2 = H2 * H2
          XCORDE = H2 + DDELT2
          IF (XCORDE .GT. 0.) THEN
* Surface entre la corde et l'arc
            XUNITE = SQRT(XCORDE) / RAYONS(IRR) / 2.
            XALPHA = ASIN(XUNITE)
            XUNITE = XALPHA - COS(XALPHA) * XUNITE
            GC = XUNITE * RAYONS(IRR) * RAYONS(IRR) / DDELTA
          ELSE
            GC = 0.
          ENDIF
*
          IF (IRC .EQ. IRR) THEN
            GC = - GC
          ENDIF
        ENDIF
*
* Longueur Moyenne
        ZZR(IH) = ZZH + GC - GP
*
      ENDDO
*
* Dernier
        IF (ZZE .LT. HC) THEN
          ZZH = 0.
        ELSE
          ZZH = ZZE - HC
        ENDIF
        ZZH = (ZZH + HXRAYO(NHMAX+1) - HPRAYO) * 0.5
        ZZR(NHMAX) = ZZH - GC
        HXRAYO(NHMAX+1) = ZZE
*
      RETURN
      END
