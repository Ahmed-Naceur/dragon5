*DECK SAXPY
       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* constant times a vector plus a vector. Uses unrolled loop for
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
* N       number of components in the vector.
* SA      constant.
* SX      RHS vector.
* INCX    increment in RHS vector.
*
*Parameters: output
* SY      LHS vector.
* INCY    increment in LHS vector.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL SX(N*INCX),SY(N*INCY),SA
      INTEGER N,INCX,INCY
*----
*  LOCAL VARIABLES
*----
      INTEGER I,IX,IY,M,MP1
*
      IF(N.LE.0)RETURN
      IF (SA .EQ. 0.0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
*----
*  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1.
*----
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
*----
*  CODE FOR BOTH INCREMENTS EQUAL TO 1. CLEAN-UP LOOP.
*----
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
      END
