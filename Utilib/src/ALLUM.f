*DECK ALLUM
      SUBROUTINE ALLUM(L4,ASS,VEC,Z,MU1,IMA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* multiplication of a general matrix in compressed diagonal storage
* mode by a vector. Z=ASS*VEC
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
*         DIMENSION ASS(IMA(L4))
* VEC     vector to multiply.
* MU1     position of each diagonal element in vector ASS.
* IMA     position of the first non-zero column element in vector ASS.
*
*Parameters: output
* Z       solution of the multiplication.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER L4,MU1(L4),IMA(L4)
      REAL ASS(*),VEC(L4),Z(L4)
*
      KEY1=MU1(1)
      Z(1)=ASS(KEY1)*VEC(1)
      DO 50 I=2,L4
      ZK=0.0
      DO 30 J=IMA(I-1)-MU1(I)+I+1,I
      KEY1=KEY1+1
      ZK=ZK+ASS(KEY1)*VEC(J)
   30 CONTINUE
      Z(I)=ZK
      ZK=VEC(I)
      DO 40 J=I-1,MU1(I)+I-IMA(I),-1
      KEY1=KEY1+1
      Z(J)=Z(J)+ASS(KEY1)*ZK
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
