*DECK ALSVDS
      SUBROUTINE ALSVDS(U,W,V,M,N,MP,NP,B,X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* linear system solution after singular value decomposition.
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
* U       first decomposed matrix. DIMENSION U(MP,NP)
* W       singular values. DIMENSION W(NP)
* V       second decomposed matrix. DIMENSION V(NP,NP)
* M,N     first/second mathematical dimension of matrix A
* MP,NP   first/second physical dimension of matrix A
* B       RHS vector. DIMENSION B(MP)
*
*Parameters: output
* X       solution vector. DIMENSION X(NP)
*
*-----------------------------------------------------------------------
*
      INTEGER M,MP,N,NP
      DOUBLE PRECISION B(MP),U(MP,NP),V(NP,NP),W(NP),X(NP),S
      INTEGER I,J,JJ
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TMP
*
      ALLOCATE(TMP(NP))
      DO 12 J=1,N
        S=0.0D0
        IF(W(J).NE.0.0D0)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.0D0
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      DEALLOCATE(TMP)
      RETURN
      END
