*DECK MCGDS6
      SUBROUTINE MCGDS6(NGEFF,NPJJM,NREG,PJJD,VOL,PJJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Volume normalization and conversion from double to single precision
* of the PJJ.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* NGEFF   number of groups to process.
* NPJJM   number of pjj modes to store for LPJJAN option.
* NREG    number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* PJJD    PJJ in double precision to be normalized.
* VOL     region volumes.
*
*Parameters: output
* PJJ     PJJ in sigle precision, normalized.
* 
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER NGEFF,NPJJM,NREG
      DOUBLE PRECISION PJJD(NREG,NPJJM,NGEFF)
      REAL VOL(NREG),PJJ(NREG,NPJJM,NGEFF)
*---
* LOCAL VARIABLES
*---
      INTEGER II,I,IMOD
      DOUBLE PRECISION VID
*
      DO 30 I=1,NREG
         VID=2.D0/DBLE(VOL(I))
         DO 20 II=1,NGEFF
         DO 10 IMOD=1,NPJJM
         PJJ(I,IMOD,II)=REAL(PJJD(I,IMOD,II)*VID)
 10      CONTINUE
 20      CONTINUE
 30   CONTINUE
*
      RETURN
      END
