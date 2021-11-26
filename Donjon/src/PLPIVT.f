*DECK PLPIVT
      SUBROUTINE PLPIVT(N,M,L,K,P,IROW,ICOL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Performs a pivot (L,K) on matrix P.
* PLPIVT = Linear Programmation PIVoT
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert and R. Chambon
*
*Parameters: input
* N       number of control variables.
* M       number of constraints.
* L       row index.
* K       column index.
* P       coefficient matrix.
* IROW    permutation vector for row elements.
* ICOL    permutation vector for column elements.
*
*Parameters: ouput
* P       coefficient matrix.
* IROW    permutation vector for row elements.
* ICOL    permutation vector for column elements.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     M,N,L,K,IROW(N),ICOL(N+1)
      DOUBLE PRECISION  P(N,M)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  PVT
      INTEGER     I,J
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PLJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PLJ(M))
*
      PVT = 1.0D0/P(L,K)
      DO 10 J=1,M
         PLJ(J) = -PVT*P(L,J)
         P(L,J) =  PLJ(J)
   10 CONTINUE
*
      DO 30 I=1,N
         IF (I.EQ.L) GO TO 30
         DO 20 J=1,M
            IF (J.EQ.K) GO TO 20
            P(I,J) = P(I,J) + PLJ(J)*P(I,K)
   20    CONTINUE
         P(I,K) = PVT*P(I,K)
   30 CONTINUE
*
      P(L,K) = PVT
      I = IROW(L)
      IROW(L) = ICOL(K)
      ICOL(K) = I
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PLJ)
      RETURN
      END
