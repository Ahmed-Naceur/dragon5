*DECK VALPL
      FUNCTION VALPL(L,U)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Return the Legendre function coefficients for the nodal collocation
* method.
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
* L     order of the Legendre polynomial.
* U     indemendent variable.
*                                                                      
*----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      IF (L.EQ.0) P=1.0
      IF (L.EQ.1) P=2.0*U
      IF (L.EQ.2) P=(6.0*U*U-0.5)
      IF (L.EQ.3) P=(20.0*U*U-3.0)*U
      IF (L.EQ.4) P=(U*U*(70.0*U*U-15.0)+0.375)
      VALPL=SQRT(REAL(2*L+1))*P
      RETURN
      END
