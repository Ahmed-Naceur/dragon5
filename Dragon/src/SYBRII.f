*DECK SYBRII
      SUBROUTINE SYBRII(PP,SG,TAU0,TAUI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $R_{ii}$ function in 1D slab geometry.
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
* TAU0   side-to-side optical path.
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
      REAL PP,SG,TAU0,TAUI
*
      IF(TAUI.GT.1.0E-10) THEN
         PP=SG*TABEN(2,TAU0)*TAUI-(TABEN(3,TAU0)-TABEN(3,TAU0+SG*TAUI))
      ELSE
         PP=0.5*(TAUI**2)*TABEN(1,TAU0)
      ENDIF
      RETURN
      END
