*DECK ALLDLS
      SUBROUTINE ALLDLS (L4,MU1,ASS,F)
*
*-----------------------------------------------------------------------
*
*Purpose:
* solution of a symmetric linear system where the coefficient matrix
* have been factorized by a preceding call to ALLDLF.
*
*Copyright:
* Copyright (C) 1989 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* L4      order of the coefficient matrix.
* MU1     position of each diagonal element in vector ASS.
* ASS     LDL(T) factors of the coefficient matrix in compressed
*         diagonal storage mode. DIMENSION ASS(MU1(L4)-MU1(1)+1)
* F       right-hand side of the linear system.
*
*Parameters: output
* F       solution of the linear system.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER L4,MU1(L4)
      REAL ASS(*),F(L4)
*
      IF (L4.EQ.1) GO TO 60
      K1=MU1(1)+1
      DO 20 I=2,L4
      K2=MU1(I)
      KJ=I-K2+K1
      T=-F(I)
      DO 10 K=K1,K2-1
      T=T+F(KJ)*ASS(K)
      KJ=KJ+1
   10 CONTINUE
      K1=K2+1
      F(I)=-T
   20 CONTINUE
C
      DO 30 I=1,L4
      F(I)=F(I)*ASS(MU1(I))
   30 CONTINUE
C
      K2=MU1(L4)
      DO 50 I=L4,2,-1
      T=-F(I)
      K1=MU1(I-1)+1
      KJ=I-K2+K1
      DO 40 K=K1,K2-1
      F(KJ)=F(KJ)+ASS(K)*T
      KJ=KJ+1
   40 CONTINUE
      K2=K1-1
   50 CONTINUE
      RETURN
C
   60 F(1)=F(1)*ASS(MU1(1))
      RETURN
      END
