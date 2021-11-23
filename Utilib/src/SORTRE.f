*DECK SORTRE
      SUBROUTINE SORTRE(N,ARRAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Sort an array of reals by increasing order using the "Insertion sort"
* method. Based on the C routine available at
* http://linux.wku.edu/~lamonml/algor/sort/insertion.html
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
* N       size of the array.
*    
*Parameters: input/output
* ARRAY   array of reals.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N
      REAL ARRAY(N)
*---
* LOCAL VARIABLES
*---
      INTEGER I,J
      REAL WORK
*
      DO I = 2, N
         WORK = ARRAY(I)
         J = I
         DO WHILE ((J.GT.1).AND.(ARRAY(J-1).GT.WORK))
            ARRAY(J) = ARRAY(J-1)
            J = J - 1
         ENDDO
         ARRAY(J) = WORK
      ENDDO
*
      RETURN
      END
