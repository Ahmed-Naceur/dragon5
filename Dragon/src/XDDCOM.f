*DECK XDDCOM
      FUNCTION XDDCOM(DBLE1,DBLE2,DEPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To compare two double precision values.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Le Tellier.
*
*Parameters: input
* DBLE1   first double precision value.
* DBLE2   second double precision value.
* DEPS    comparison criterion.
*
*Parameters: output
* XDDCOM  comparison flag.
*
*-----------------------------------------------------------------------
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      DOUBLE PRECISION DBLE1,DBLE2,DEPS
      LOGICAL XDDCOM
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION DLIM,DIFF
*
      DLIM=DEPS*1.D-3
      IF (ABS(DBLE1).GT.DLIM) THEN
         DIFF=(DBLE2/DBLE1-1.D0)
         XDDCOM=(ABS(DIFF).LT.DEPS)
      ELSE
         XDDCOM=(ABS(DBLE2).LT.DLIM)
      ENDIF
*
      END
