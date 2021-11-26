*DECK SYB4TN
      SUBROUTINE SYB4TN (NHMAX,IXRAYO,ISDEBU,COSECT,NRI,RAYONS,DELTAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Search of the next singular point in a Cartesian cell.
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
* NHMAX   number of intervals (the number of interceptions is NHMAX+1).
* IXRAYO  tubes indices from west to east.
* ISDEBU  index of the first sector.
* COSECT  cosinus of the sector angles.
* NRI     number of radius.
* RAYONS  radius of the tubes.
*
*Parameters: input/output
* DELTAC  next distance.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER        NHMAX,IXRAYO(NHMAX),ISDEBU,NRI
      REAL           COSECT(*),RAYONS(NRI),DELTAC
*----
*  LOCAL VARIABLES
*----
      LOGICAL        LGSEC1,LGSEC2
*
      IF (NHMAX .LT. 2) RETURN
      IHX = 0
      IR1 = IXRAYO(1)
      IR2 = IXRAYO(2)
      LGSEC2 = IR1 .EQ. IR2
      IF (LGSEC2) THEN
        ISC = ISDEBU + 1
      ELSE
        ISC = ISDEBU
      ENDIF
*
      DO IHC = 2, NHMAX - 1
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
            DELTAC = RAYONS(IR1)
            IHX = IHC
          ENDIF
        ENDIF
*
* Coin Secteur/Rayon
        ELSEIF (LGSEC1) THEN
        IF (IR2 .EQ. IR1+1) THEN
          RR  = RAYONS(IR1)
          DD  = RR * COSECT(ISC-1)
          IF (DD .LT. DELTAC) THEN
            DELTAC = DD
            IHX = IHC
          ENDIF
        ENDIF
*
* Coin Secteur/Rayon
        ELSEIF (LGSEC2) THEN
        IF (IR0 .EQ. IR1+1) THEN
          RR  = RAYONS(IR1)
          DD  = RR * COSECT(ISC-1)
          IF (DD .LT. DELTAC) THEN
            DELTAC = DD
            IHX = IHC
          ENDIF
        ENDIF
*
        ENDIF
*
      ENDDO
*
      RETURN
      END
