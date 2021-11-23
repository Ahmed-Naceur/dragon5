*DECK ADDBLE
      SUBROUTINE ADDBLE(NSIZE,IPOS,DSUM,DINC)
*-----------------------------------------------------------------------
*
*Purpose:
* add a double to a vector of doubles.
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
* NSIZE   size of double array.
* IPOS    position in double array.
* DINC    double to be added.
*
*Parameters: input/output
* DSUM    double array.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NSIZE,IPOS
      DOUBLE PRECISION DSUM(NSIZE),DINC
*
      DSUM(IPOS) = DSUM(IPOS) + DINC
*
      RETURN
      END
