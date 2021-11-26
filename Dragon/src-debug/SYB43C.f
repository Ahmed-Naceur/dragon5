*DECK SYB43C
      SUBROUTINE SYB43C(PP,TAUI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $D_i$ function in 1D spherical geometry.
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
* TAUI   optical path in volume i.
*
*Parameters: output
* PP     value of the expression.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL PP,TAUI
*
      IF(TAUI.GT.0.004) THEN
         PP=TAUI-1.0+EXP(-TAUI)
      ELSE IF(TAUI.GT.0.002) THEN
         PP=TAUI-1.0+EXP(-TAUI)
         PQ=0.5*(TAUI**2)
         FACT=500.0*(TAUI-0.002)
         PP=PP*FACT+PQ*(1.0-FACT)
      ELSE
         PP=0.5*(TAUI**2)
      ENDIF
      RETURN
      END
