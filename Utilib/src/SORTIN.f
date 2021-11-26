*DECK SORTIN
      SUBROUTINE SORTIN(N,ARRAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Sort an array of integers by increasing order using the "Insertion sort"
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
* Parameters: input/output
* ARRAY   array of integers.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N,ARRAY(N)
*---
* LOCAL VARIABLES
*---
      INTEGER I,J,INDEX
*
      DO I = 2, N
         INDEX = ARRAY(I)
         J = I
         DO WHILE ((J.GT.1).AND.(ARRAY(J-1).GT.INDEX))
            ARRAY(J) = ARRAY(J-1)
            J = J - 1
         ENDDO
         ARRAY(J) = INDEX
      ENDDO
*
      RETURN
      END
