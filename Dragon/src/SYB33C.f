*DECK SYB33C
      SUBROUTINE SYB33C (PPLUS,TAUP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $D_i$ function in 1D cylindrical and 2D geometry.
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
* TAUP   optical path.
*
*Parameters: output
* PPLUS  value of the expression.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL PPLUS,TAUP
*
      IF (TAUP .EQ. 0.0) THEN
         PPLUS=0.0
      ELSE IF (TAUP .GT. 0.004) THEN
         PPLUS=TAUP+TABKI(3,TAUP)-TABKI(3,0.0)
      ELSE IF (TAUP .GT. 0.002) THEN
         PPLUS=TAUP+TABKI(3,TAUP)-TABKI(3,0.0)
         PQLUS= 0.5*TAUP*TAUP*TABKI(1,TAUP*0.5)
         FACT=500.0*(TAUP-0.002)
         PPLUS=PPLUS*FACT+PQLUS*(1.0-FACT)
      ELSE
         PPLUS=0.5*TAUP*TAUP*TABKI(1,TAUP*0.5)
      ENDIF
      RETURN
      END
