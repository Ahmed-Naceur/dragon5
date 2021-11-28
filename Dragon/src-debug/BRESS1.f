      SUBROUTINE BRESS1(ITRIAL,XX,DIFF,SIGR,A11)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of leakage system matrices for the nodal expansion method.
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
* SIGR    macroscopic removal cross section.
* SIGT    macroscopic cross section.
*
*Parameters: output
* A11     nodal matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITRIAL
      REAL XX,DIFF,SIGR,A11(4,4)
*----
*  WEIGHT RESIDUAL EQUATIONS
*----
      A11(:4,:4)=0.0
      DX2=XX**2
      ETA=XX*SQRT(SIGR/DIFF)
      A11(1,1)=SIGR/12.0
      A11(2,2)=SIGR/180.0
      IF (ITRIAL == 1) THEN
        A11(1,3)=-SIGR/120.0-DIFF/(2.0*DX2)
        A11(2,4)=-SIGR/2100.0-DIFF/(15.0*DX2)
      ELSE IF (ITRIAL == 2) THEN
        ALP1=ETA*COSH(ETA/2)-2.0*SINH(ETA/2)
        ALP2=((12.0+ETA**2)*SINH(ETA/2)-6.0*ETA*COSH(ETA/2))/(3.0*ETA)
        A11(1,3)=(SIGR/(ETA**2)-DIFF/DX2)*ALP1
        A11(2,4)=(SIGR/(ETA**2)-DIFF/DX2)*ALP2
      ENDIF
      ! LEFT CURRENT CONDITION
      A11(3,1)=-DIFF/XX
      A11(3,2)=DIFF/XX
      IF (ITRIAL == 1) THEN
        A11(3,3)=-DIFF/(2.0*XX)
        A11(3,4)=DIFF/(5.0*XX)
      ELSE IF (ITRIAL == 2) THEN
        A11(3,3)=-(DIFF/XX)*ETA*COSH(ETA/2)
        A11(3,4)=(DIFF/XX)*ETA*SINH(ETA/2)
      ENDIF
      ! RIGHT CURRENT CONDITION
      A11(4,1)=-DIFF/XX
      A11(4,2)=-DIFF/XX
      IF (ITRIAL == 1) THEN
        A11(4,3)=-DIFF/(2.0*XX)
        A11(4,4)=-DIFF/(5.0*XX)
      ELSE IF (ITRIAL == 2) THEN
        A11(4,3)=-(DIFF/XX)*ETA*COSH(ETA/2)
        A11(4,4)=-(DIFF/XX)*ETA*SINH(ETA/2)
      ENDIF
      RETURN
      END
