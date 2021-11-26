*DECK ALDDLF
      SUBROUTINE ALDDLF (L4,ASS,MU1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* in-place L-D-L(T) factorization of a symmetric positive definite
* matrix in compressed diagonal storage mode.
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
*         A(I,J)=ASS(MU1(I)-I+J) if J.le.I and J.gt.I+MU1(I-1)-MU1(I)
*               =A(J,I)          if I.lt.J
*               =0.0             else
*         DIMENSION ASS(MU1(L4)-MU1(1)+1)
* MU1     position of each diagonal element in vector ASS.
*
*Parameters: output
* ASS     LDL(T) factors.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER L4,MU1(L4)
      DOUBLE PRECISION ASS(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION S,R
*
      ASS(MU1(1))=1.0D0/ASS(MU1(1))
      IF (L4.EQ.1) RETURN
      DO 10 K=2,L4
      K1=MU1(K)-K
      KM=MU1(K-1)+1-K1
      IF(KM+1.GT.K) GO TO 7
      DO 2 I=KM+1,K-1
      R=ASS(K1+I)
      ASS(K1+I)=0.0D0
      S=0.0D0
      I1=MU1(I)-I
      IM=MU1(I-1)+1-I1
      DO 5 J=MAX0(IM,KM),I
      S=S+ASS(K1+J)*ASS(I1+J)
 5    CONTINUE
      ASS(K1+I)=R-S
 2    CONTINUE
      S=0.0D0
      DO 6 I=KM,K-1
      R=ASS(K1+I)
      ASS(K1+I)=R*ASS(MU1(I))
      S=S+R*ASS(K1+I)
 6    CONTINUE
      ASS(MU1(K))=ASS(MU1(K))-S
 7    ASS(MU1(K))=1.0D0/ASS(MU1(K))
 10   CONTINUE
      RETURN
      END
