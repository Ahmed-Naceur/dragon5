*DECK AL1EIG
      SUBROUTINE AL1EIG(N,A,EPSOUT,MAXOUT,ITER,EVECT,EVAL,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the fundamental eigenvalue and corresponding eigenvector of
* equation (A-EVAL)*EVECT=0 using the power method.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* N       number of unknowns
* A       coefficient matrix
* EPSOUT  convergence epsilon for the power method
* MAXOUT  maximum number of iterations for the power method
* EVECT   initial estimate
* IPRINT  print parameter
*
*Parameters: output
* ITER    number of iterations
* EVECT   corresponding eigenvector
* EVAL    fondamental eigenvalue
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,MAXOUT,ITER,IPRINT
      REAL A(N,N),EPSOUT,EVECT(N),EVAL
*----
*  LOCAL VARIABLES
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR
*----
*  POWER METHOD
*----
      EVECT(:N)=1.0
      ITER=0;
      ALLOCATE(GAR(N))
      DO
        ITER=ITER+1
        IF (ITER > MAXOUT) CALL XABORT('AL1EIG: UNABLE TO CONVERGE.')
        GAR(:)=EVECT(:)
        EVECT(:)=MATMUL(A(:,:),EVECT(:))
        EVAL=SQRT(DOT_PRODUCT(EVECT(:),EVECT(:)))
        EVECT(:)=EVECT(:)/EVAL
        ERR1=MAXVAL(ABS(EVECT))
        ERR2=MAXVAL(ABS(GAR(:)-EVECT(:)))
        IF(IPRINT.GT.1) THEN
          IF (MOD(ITER,5) == 1) WRITE(6,10) ITER,EVAL,ERR2
        ENDIF
        IF(ERR2 <= ERR1*EPSOUT) EXIT
      ENDDO
      DEALLOCATE(GAR)
      RETURN
   10 FORMAT(14H AL1EIG: ITER=,I6,6H EVAL=,1P,E12.5,7H ERROR=,E11.4)
      END
