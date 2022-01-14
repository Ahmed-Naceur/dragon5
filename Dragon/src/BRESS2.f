      SUBROUTINE BRESS2(ITRIAL,XX,DIFF,SIGR,SIGT,B11)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of non-leakage system matrices for the nodal expansion
* method.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ITRIAL  type of base (=1: polynomial; =2: hyperbolic)
* XX      mesh width.
* DIFF    diffusion coefficient.
* SIGR    macroscopic removal cross section.
* SIGT    macroscopic cross section.
*
*Parameters: output
* B11     nodal matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITRIAL
      REAL XX,DIFF,SIGR,SIGT,B11(4,4)
*----
*  WEIGHT RESIDUAL EQUATIONS
*----
      B11(:4,:4)=0.0
      ETA=XX*SQRT(SIGR/DIFF)
      B11(1,1)=SIGT/12.0
      B11(2,2)=SIGT/180.0
      IF (ITRIAL == 1) THEN
        B11(1,3)=-SIGT/120.0
        B11(2,4)=-SIGT/2100.0
      ELSE IF (ITRIAL == 2) THEN
        ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
        ALP2=((12.0+ETA**2)*SINH(ETA/2)-6.0*ETA*COSH(ETA/2))/(3.0*ETA)
        B11(1,3)=SIGT*ALP1/(ETA**2)
        B11(2,4)=SIGT*ALP2/(ETA**2)
      ENDIF
      RETURN
      END
