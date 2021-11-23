*DECK RENINS
      SUBROUTINE RENINS(SIZE,LEV,DEG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Sort a level by increasing degree using the "Insertion method".
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
* SIZE    number of nodes in the level.
* LEV     level to sort.
* DEG     degrees of the level.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER SIZE,LEV(SIZE),DEG(SIZE)
*---
* LOCAL VARIABLES
*---
      INTEGER I,INDD,INDL,J
*
      DO I=2,SIZE
         INDD=DEG(I)
         INDL=LEV(I)
         J=I
         DO WHILE ((J.GT.1).AND.(DEG(J-1).GT.INDD))
            DEG(J)=DEG(J-1)
            LEV(J)=LEV(J-1)
            J=J-1
         ENDDO
         DEG(J)=INDD
         LEV(J)=INDL
      ENDDO
*
      RETURN
      END
