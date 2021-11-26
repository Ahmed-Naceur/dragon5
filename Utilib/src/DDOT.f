*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT(N,SX,INCX,SY,INCY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* forms the dot product of two vectors. Uses unrolled loops for
* increments equal to one.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): Jack Dongarra, linpack, 3/11/78.
*
*Parameters: input
* N       number of components in the vectors.
* SX      first vector.
* INCX    increment in first vector.
* SY      second vector.
* INCY    increment in second vector.
*
*Parameters: output
* DDOT    dot product.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,INCX,INCY
      DOUBLE PRECISION SX(N*INCX),SY(N*INCY)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,IX,IY,M,MP1
      DOUBLE PRECISION DTEMP
*
      DTEMP = 0.0D0
      DDOT = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
*----
*  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1.
*----
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*----
*   CODE FOR BOTH INCREMENTS EQUAL TO 1. CLEAN-UP LOOP.
*----
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     *   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
