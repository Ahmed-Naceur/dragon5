*DECK ALPINV
      SUBROUTINE ALPINV(M,N,A,AINV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* pseudo inversion of a non singular matrix using Gaussian elimination
* with partial pivoting. Real version.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* M       first dimension of matrix A.
* N       second dimension of matrix A.
* A       coefficient matrix to be inverted.
*
*Parameters: output
* AINV    pseudo inverted matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER M,N
      REAL A(M,N),AINV(N,M)
*----
*  LOCAL VARIABLES
*----
      REAL, DIMENSION(:,:), ALLOCATABLE :: B
*
      ALLOCATE(B(N,N))
      B=MATMUL(TRANSPOSE(A),A)
      CALL ALINV(N, B, N, IER)
      IF(IER.NE.0) CALL XABORT('ALPINV: PSEUDO INVERSION FAILURE.')
      AINV=MATMUL(B, TRANSPOSE(A))
      DEALLOCATE(B)
      RETURN
      END
