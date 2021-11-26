*DECK ALST2F
      SUBROUTINE ALST2F(MDIM,M,N,A,TAU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* to obtain the QR factorization of the matrix a using Householder
* transformations. Use LAPACK's DGEQRF routine storage. Douple precision
* routine.
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
*Reference:
* P.A. BUSINGER, Num. Math. 7, 269-276 (1965).
*
*Parameters: input
* MDIM    dimensioned column length of A.
* M       number of rows of A
* N       number of columns of A. N.le.M is assumed.
* A       matrix A.
*
*Parameters: output
* A       decomposed matrix. On exit, the elements on and above the
*         diagonal of the array contain the m by n upper trapezoidal
*         matrix R (R is upper triangular if m >= n); the elements
*         below the diagonal, with the array TAU, represent the
*         orthogonal matrix Q as a product of elementary reflectors.
* TAU     scalar factors of the elementary reflectors.
*
*-----------------------------------------------------------------------
*
      IMPLICIT REAL(KIND=8)(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MDIM,M,N
      REAL(KIND=8) A(MDIM,N),TAU(N)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: W
*----
*  CHECK THE INPUT
*----
      IF(MDIM.LT.M) CALL XABORT('ALST2F: MDIM.LT.M')
      IF(N.LT.1) CALL XABORT('ALST2F: N.LT.1')
      IF(N.GT.M) THEN
         WRITE(HSMG,'(18HALST2F: N.GT.M (N=,I3,3H M=,I3,2H).)') N,M
         CALL XABORT(HSMG)
      ENDIF
*----
*  PERFORM QR FACTORIZATION.
*----
      ALLOCATE(W(M,1))
      DO J=1,N
        M1 = M-J+1; W(:M1,1) = A(J:M,J); X1 = W(1,1);
        AX = SQRT(DOT_PRODUCT(W(:M1,1),W(:M1,1)))
        A1 = ABS(X1); S = SIGN(1.0D0,W(1,1));
        SSSS = -AX*S; A1 = A1+AX;
        W(1,1) = A1*S
        DD2 = A1*AX
        IF(DD2 == 0.0D0) CALL XABORT('ALST2F: SINGULAR REFLECTION')
        W(:M1,1) = W(:M1,1)/SQRT(DD2)
        A(J:M,J) = W(:M1,1)
        IF(J < N) THEN
          A(J:M,J+1:N) = A(J:M,J+1:N)
     1      -MATMUL(W(:M1,:),(MATMUL(TRANSPOSE(W(:M1,:)),A(J:M,J+1:N))))
        ENDIF
        DIAG = A(J,J)
        A(J:M,J) = A(J:M,J)/DIAG
        A(J,J) = SSSS
        TAU(J) = -DIAG*DIAG
      ENDDO
      DEALLOCATE(W)
      RETURN
      END
