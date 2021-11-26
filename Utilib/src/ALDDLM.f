*DECK ALDDLM
      SUBROUTINE ALDDLM (L4,ASS,VEC,Z,MU1,ITY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* multiplication of a symmetric matrix in compressed diagonal storage
* mode by a vector.
* Double precision version.
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
* ASS     coefficient matrix in compressed diagonal storage mode.
*         DIMENSION ASS(MU1(L4))
* VEC     vector to multiply.
* Z       vector that will be added to the result if ITY=2.
* MU1     position of each diagonal element in vector ASS.
* ITY     type of multiplication (ITY=1: Z=ASS*VEC;
*         ITY=2: Z=Z+(ASS-DIAG(ASS))*VEC).
*
*Parameters: output
* Z       solution of the multiplication.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER L4,MU1(L4),ITY
      DOUBLE PRECISION ASS(*),VEC(L4),Z(L4)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION ZK
*
      GO TO (10,60),ITY
*
* CALCULATION OF Z=ASS*VEC.
   10 Z(1)=ASS(MU1(1))*VEC(1)
      I1=MU1(1)+1
      DO 50 K=2,L4
      I2=MU1(K)
      KEY1=I2-K
      ZK=0.0D0
      DO 30 L=I1-I2+K,K-1
      ZK=ZK+ASS(KEY1+L)*VEC(L)
      Z(L)=Z(L)+ASS(KEY1+L)*VEC(K)
   30 CONTINUE
      Z(K)=ZK+ASS(KEY1+K)*VEC(K)
      I1=I2+1
   50 CONTINUE
      RETURN
*
* CALCULATION OF Z=Z+(ASS-DIAG(ASS))*VEC.
   60 I1=MU1(1)+1
      DO 80 K=2,L4
      I2=MU1(K)
      KEY1=I2-K
      DO 70 L=I1-I2+K,K-1
      Z(K)=Z(K)+ASS(KEY1+L)*VEC(L)
      Z(L)=Z(L)+ASS(KEY1+L)*VEC(K)
   70 CONTINUE
      I1=I2+1
   80 CONTINUE
      RETURN
      END
