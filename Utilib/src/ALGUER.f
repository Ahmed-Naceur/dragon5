*DECK ALGUER
      SUBROUTINE ALGUER(A,M,X,ITS,LFAIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* find one root of a polynomial.
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
* A       polynomial coefficients. DIMENSION A(M+1)
* M       polynomial order.
*
*Parameters: output
* X       complex single root.
* ITS     number of iterations.
* LFAIL   set to .true. in case of failure.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER M,ITS
      COMPLEX*16 A(M+1),X
      LOGICAL LFAIL
*----
*  LOCAL VARIABLES
*----
      INTEGER MAXIT,MR,MT
      DOUBLE PRECISION EPSS
      PARAMETER (EPSS=2.D-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER ITER,J
      DOUBLE PRECISION ABX,ABP,ABM,ERR,FRAC(MR)
      COMPLEX*16 DX,X1,B,D,F,G,H,SQ,GP,GM,G2,TMP
      SAVE FRAC
      DATA FRAC /.5D0,.25D0,.75D0,.13D0,.38D0,.62D0,.88D0,1.D0/
*
      LFAIL=.FALSE.
      DO 12 ITER=1,MAXIT
        ITS=ITER
        B=A(M+1)
        ERR=ABS(B)
        D=CMPLX(0.D0,0.D0,KIND=KIND(D))
        F=CMPLX(0.D0,0.D0,KIND=KIND(F))
        ABX=ABS(X)
        DO 11 J=M,1,-1
          F=X*F+D
          D=X*D+B
          B=X*B+A(J)
          ERR=ABS(B)+ABX*ERR
11      CONTINUE
        ERR=EPSS*ERR
        IF(ABS(B).LE.ERR) THEN
          RETURN
        ELSE
          G=D/B
          G2=G*G
          H=G2-2.D0*F/B
          SQ=SQRT((M-1)*(M*H-G2))
          GP=G+SQ
          GM=G-SQ
          ABP=ABS(GP)
          ABM=ABS(GM)
          IF(ABP.LT.ABM) GP=GM
          IF (MAX(ABP,ABM).GT.0.D0) THEN
            DX=M/GP
          ELSE
            TMP=CMPLX(LOG(1.D0+ABX),DBLE(ITER),KIND=KIND(TMP))
            DX=EXP(TMP)
          ENDIF
        ENDIF
        X1=X-DX
        IF(X.EQ.X1)RETURN
        IF (MOD(ITER,MT).NE.0) THEN
          X=X1
        ELSE
          X=X-DX*FRAC(ITER/MT)
        ENDIF
12    CONTINUE
      LFAIL=.TRUE.
      RETURN
      END
