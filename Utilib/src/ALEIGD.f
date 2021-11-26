*DECK ALEIGD
      SUBROUTINE ALEIGD (A,B,N,EVAL,EVECT,EPS,ITER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* find the fondamental eigenvalue and corresponding eigenvector of
* equation (A-EVAL*B)*EVECT=0 using the inverse power method.
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
* A      first coefficient matrix
* B      second coefficient matrix
* N      number of unknowns
* EVECT  initial estimate
* EPS2   stopping criterion
*
*Parameters: output
* EVAL   fondamental eigenvalue
* EVECT  corresponding eigenvector
* ITER   number of iterations
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,ITER
      DOUBLE PRECISION A(N,N),B(N,N),EVAL,EVECT(N),EPS
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MMAX=1000)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GAR
*----
*  COMPUTE THE ITERATIVE MATRIX
*----
      CALL ALINVD (N,A,N,IER)
      IF(IER.EQ.1) CALL XABORT('ALEIGD: SINGULAR MATRIX.')
      ALLOCATE(GAR(N))
      DO 30 I=1,N
      DO 10 J=1,N
      GAR(J)=A(I,J)
10    CONTINUE
      DO 25 J=1,N
      A(I,J)=0.0D0
      DO 20 K=1,N
      A(I,J)=A(I,J)+GAR(K)*B(K,J)
20    CONTINUE
25    CONTINUE
30    CONTINUE
*----
*  PERFORM POWER ITERATIONS
*----
      TEST=0.0D0
      ITER=0
      EVAL=0.0D0
40    ITER=ITER+1
      IF(ITER.GT.MMAX) CALL XABORT('ALEIGD: UNABLE TO CONVERGE(1).')
      S1=0.0D0
      S2=0.0D0
      DO 60 I=1,N
      GAR(I)=0.0D0
      DO 50 J=1,N
      GAR(I)=GAR(I)+A(I,J)*EVECT(J)
50    CONTINUE
      S1=S1+GAR(I)*EVECT(I)
      S2=S2+GAR(I)**2
60    CONTINUE
      IF(S2.EQ.0.0D0) CALL XABORT('ALEIGD: DIVIDE CHECK.')
      ZZ=ABS(EVAL-S1/S2)
      EVAL=S1/S2
      ERR1=0.0D0
      ERR2=0.0D0
      DO 70 I=1,N
      ERR1=MAX(ERR1,ABS(GAR(I)*EVAL))
      ERR2=MAX(ERR2,ABS(GAR(I)*EVAL-EVECT(I)))
      EVECT(I)=GAR(I)*EVAL
70    CONTINUE
      IF((ZZ.LE.EPS).AND.(ERR2.LE.ERR1*EPS)) THEN
         DEALLOCATE(GAR)
         RETURN
      ENDIF
      IF(ITER.EQ.1) TEST=ZZ
      IF((ITER.GE.10).AND.(ZZ.GT.TEST)) CALL XABORT('ALEIGD: UNABLE TO'
     1 //' CONVERGE(2).')
      GO TO 40
      END
