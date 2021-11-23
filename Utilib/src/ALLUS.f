*DECK ALLUS
      SUBROUTINE ALLUS(L4,MU1,IMA,ASS,F)
*
*-----------------------------------------------------------------------
*
*Purpose:
* solution of a linear system where the coefficient matrix have been
* factorized by a preceding call to ALLUF.
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
* IMA     position of the first non-zero column element in vector ASS.
* ASS     LU factors of the coefficient matrix in compressed diagonal
*         storage mode. DIMENSION ASS(IMA(L4))
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
      INTEGER L4,MU1(L4),IMA(L4)
      REAL ASS(*),F(L4)
*
      F(1)=F(1)/ASS(MU1(1))
      DO 20 I=2,L4
      K1=IMA(I-1)+1
      K2=MU1(I)
      KJ=I-K2+K1
      T=-F(I)
      DO 10 K=K1,K2-1
      T=T+F(KJ)*ASS(K)
      KJ=KJ+1
   10 CONTINUE
      F(I)=-T/ASS(MU1(I))
   20 CONTINUE
*
      DO 40 I=L4,2,-1
      K1=IMA(I)
      K2=MU1(I)
      KJ=I-K1+K2
      T=-F(I)
      DO 30 K=K1,K2+1,-1
      F(KJ)=F(KJ)+ASS(K)*T
      KJ=KJ+1
   30 CONTINUE
   40 CONTINUE
      RETURN
      END
