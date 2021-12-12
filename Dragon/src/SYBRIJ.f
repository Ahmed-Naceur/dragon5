*DECK SYBRIJ
      SUBROUTINE SYBRIJ(PP,SG,TAU0,TAUI,TAUJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $R_{ij}$ function in 1D slab geometry.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* SG     specular parity sign.
* TAU0   side to side optical path.
* TAUI   optical path in volume i (or volume j).
* TAUJ   optical path in volume j (or volume i).
*
*Parameters: output
* PP   value of the expression.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL PP,SG,TAU0,TAUI,TAUJ
*----
*  SYMMETRIC FORMULA IN (I,J). WE SET POPI <= POPJ
*----
      IF(TAUJ.LT.TAUI) THEN
         POPI = TAUJ
         POPJ = TAUI
      ELSE
         POPI = TAUI
         POPJ = TAUJ
      ENDIF
      IF((POPI.GT.1.0E-10).AND.(POPJ.GT.1.0E-10)) THEN
         PP=0.5*(TABEN(3,TAU0)-TABEN(3,TAU0+SG*POPI)-
     1           TABEN(3,TAU0+SG*POPJ)+TABEN(3,TAU0+SG*POPI+SG*POPJ))
      ELSE IF(POPJ.GT.1.0E-10) THEN
         PP=SG*0.5*POPI*(TABEN(2,TAU0)-TABEN(2,TAU0+SG*POPJ))
      ELSE IF(POPI.GT.1.0E-10) THEN
         PP=SG*0.5*POPJ*(TABEN(2,TAU0)-TABEN(2,TAU0+SG*POPI))
      ELSE
         PP=0.5*POPI*POPJ*TABEN(1,TAU0)
      ENDIF
      RETURN
      END
