*DECK SYB7TN
      SUBROUTINE SYB7TN (IHMIN,IHMAX,IXRAYO,ISDEBU,COSECT,NRI,RAYONS,
     & TYSUIT,IHSUIT,DELTAC,IRSUIT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Search of the next singular point in an hexagonal cell.
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
* IHMIN   index of the first tube.
* IHMAX   index of the last tube.
* IXRAYO  tube indices from west to east.
* ISDEBU  index of the first sector.
* COSECT  sector angle cosinus.
* NRI     number of radius.
* RAYONS  radius of the tubes.
*
*Parameters: input/output
* DELTAC  next distance.
*
*Parameters: output
* TYSUIT  type of the next singular point.
* IHSUIT  index of the next singular point.
* IRSUIT  index of the preceding tube.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER        IHMIN,IHMAX,IXRAYO(IHMAX),ISDEBU,NRI,IHSUIT,IRSUIT
      REAL           COSECT(3),RAYONS(NRI),DELTAC
      CHARACTER      TYSUIT*4
*----
*  LOCAL VARIABLES
*----
      LOGICAL        LGSEC1
      LOGICAL        LGSEC2
*
      IR1 = IXRAYO(IHMIN)
      IR2 = IXRAYO(IHMIN+1)
      LGSEC2 = IR1 .EQ. IR2
      IF (LGSEC2) THEN
        ISC = ISDEBU
      ELSE
        ISC = ISDEBU - 1
      ENDIF
*
      DO IHC = IHMIN + 1, IHMAX - 1
        LGSEC1 = LGSEC2
        IR0 = IR1
        IR1 = IR2
*
        IR2 = IXRAYO(IHC+1)
        LGSEC2 = IR1 .EQ. IR2
        IF (LGSEC2) THEN
          ISC = ISC + 1
        ENDIF
*
* Tangente = Intersection Couronne
        IF (IR2 .EQ. IR0) THEN
        IF (IR2 .EQ. IR1+1) THEN
          IF (RAYONS(IR1) .LT. DELTAC) THEN
            IHSUIT = IHC
            DELTAC = RAYONS(IR1)
            TYSUIT = 'Tang'
          ENDIF
        ENDIF
*
* Coin du 1er Secteur
        ELSEIF (ISC .EQ. 1) THEN
         IF (LGSEC2) THEN
          RR  = RAYONS(IR1)
          DD  = RR * COSECT(ISC)
          IF (DD .LT. DELTAC) THEN
            IHSUIT = IHC
            DELTAC = DD
            TYSUIT = 'Coin'
            IRSUIT = IR1 + 1
          ENDIF
         ENDIF
*
* Coin d'un Secteur Est
        ELSEIF (LGSEC1) THEN
          RR  = RAYONS(IR1)
          DD  = RR * COSECT(ISC)
          IF (DD .LT. DELTAC) THEN
            IHSUIT = IHC
            DELTAC = DD
            TYSUIT = 'Coin'
            IRSUIT = IR1 + 1
          ENDIF
        ENDIF
*
      ENDDO
*
      RETURN
      END
