*DECK ALINVC
      SUBROUTINE ALINVC(N,A,MAX,IER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* in-place inversion of a non singular matrix using gaussian elimination
* with partial pivoting. COMPLEX*16 version.
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* N       order of the coefficient matrix.
* A       coefficient matrix to be inverted.
* MAX     first dimention of matrix A.
*
*Parameters: output
* A       inverted matrix.
* IER     error flag (execution failure if IER.ne.0).
*
*-----------------------------------------------------------------------
*
      IMPLICIT COMPLEX*16 (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,MAX,IER
      COMPLEX*16 A(MAX,N)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TEST
*----
*  ALLOCATABLE ARRAY
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: IND
*
      ALLOCATE(IND(N))
      IN=0
      IER=0
      DO 1 I=1,N
      IND(I)=I
1     CONTINUE
      DO 12 J=1,N
      TEST=0.0D0
      DO 2 I=J,N
      IF (ABS(A(I,J)).LE.TEST) GO TO 2
      TEST=ABS(A(I,J))
      IN=I
2     CONTINUE
      IF (TEST.NE.0.0D0) GO TO 3
      IER=1
      DEALLOCATE(IND)
      RETURN
3     PMX=A(IN,J)
      A(IN,J)=1.0D0
      DO 4 I=1,N
      PER=A(IN,I)/PMX
      A(IN,I)=A(J,I)
      A(J,I)=PER
4     CONTINUE
      IPER=IND(IN)
      IND(IN)=IND(J)
      IND(J)=IPER
      DO 11 I=1,N
      IF (I.EQ.J) GO TO 11
      PMX=A(I,J)
      A(I,J)=0.0D0
      DO 9 K=1,N
      A(I,K)=A(I,K)-PMX*A(J,K)
9     CONTINUE
11    CONTINUE
12    CONTINUE
      DO 16 J=1,N
      DO 13 K=J,N
      IF (IND(K).NE.J) GO TO 13
      IN=K
      GO TO 14
13    CONTINUE
14    DO 15 I=1,N
      PER=A(I,J)
      A(I,J)=A(I,IN)
      A(I,IN)=PER
      IND(IN)=IND(J)
15    CONTINUE
16    CONTINUE
      DEALLOCATE(IND)
      RETURN
      END
