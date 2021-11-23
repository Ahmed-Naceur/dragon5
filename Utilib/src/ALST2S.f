*DECK ALST2S
      SUBROUTINE ALST2S(MDIM,M,N,A,TAU,B,X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* to solve the least squares problem A*X=B when the matrix a has
* already been decomposed by ALST2F.
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
* MDIM    dimensioned column length of A.
* M       number of rows of A
* N       number of columns of A. N.le.M is assumed.
* A       decomposed matrix.
* TAU     scalar factors of the elementary reflectors.
* B       right-hand side.
*
*Parameters: output
* B       B has been clobbered.
*         SQRT(SUM(I=N+1,M)(B(I)**2)) is the L2 norm of the residual
*         in the solution of the equations.
* X       solution vectors. X=B IS OK.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MDIM,M,N
      DOUBLE PRECISION A(MDIM,N),TAU(N),B(M),X(N)
*----
*  CHECK THE INPUT.
*----
      IF(MDIM.LT.M) CALL XABORT('ALST2S: MDIM.LT.M')
      IF(N.LT.1) CALL XABORT('ALST2S: N.LT.1')
      IF(N.GT.M) CALL XABORT('ALST2S: N.GT.M')
*----
*  APPLY Q-TRANSPOSE TO B.
*----
      DO J=1,N
        IF((TAU(J).EQ.0.0D0).OR.(A(J,J).EQ.0.0D0)) THEN
          CALL XABORT('ALST2S: TAU(J)=0 OR A(J,J)=0')
        ENDIF
        S=B(J)
        DO I=J+1,M
          S=S+A(I,J)*B(I)
        ENDDO
        S=S*TAU(J)
        B(J)=B(J)+S
        DO I=J+1,M
          B(I)=B(I)+S*A(I,J)
        ENDDO
      ENDDO
*----
*  BACK-SOLVE THE TRIANGULAR SYSTEM U*X=(Q-TRANSPOSE)*B.
*----
      X(N)=B(N)/A(N,N)
      DO II=2,N
        I=N+1-II
        S=B(I)
        DO J=I+1,N
          S=S-A(I,J)*X(J)
        ENDDO
        X(I)=S/A(I,I)
      ENDDO
      RETURN
      END
