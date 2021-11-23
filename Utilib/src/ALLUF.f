*DECK ALLUF
      SUBROUTINE ALLUF(L4,ASS,MU1,IMA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LU factorization of a general positive definite matrix in compressed
* diagonal storage mode.
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
*         A(I,J)=ASS(MU1(I)-I+J) if J.le.I and J.gt.I+IMA(I-1)-MU1(I)
*               =ASS(MU1(J)+J-I) if I.le.J and I.ge.J-IMA(J)+MU1(J)
*               =0.0             else
*         DIMENSION ASS(IMA(L4))
* MU1     position of each diagonal element in vector ASS.
* IMA     position of the first non-zero column element in vector ASS.
*
*Parameters: output
* ASS     LU factors.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER L4,MU1(L4),IMA(L4)
      REAL ASS(*)
*
      DO 120 I=2,L4
      MU1IMI=MU1(I)-I
      MU1IPI=MU1(I)+I
      IND1=IMA(I-1)-MU1IMI+1
      IND2=MU1IPI-IMA(I)
      IF (IND1.LE.IND2) THEN
         DO 20 J=IND1,IND2-1
         MU1JPJ=MU1(J)+J
         SUM=0.0
         DO 10 K=MAX(IND1,MU1JPJ-IMA(J)),J-1
         SUM=SUM+ASS(MU1IMI+K)*ASS(MU1JPJ-K)
   10    CONTINUE
         ASS(MU1IMI+J)=ASS(MU1IMI+J)-SUM
   20    CONTINUE
         DO 50 J=IND2,I-1
         MU1JMJ=MU1(J)-J
         MU1JPJ=MU1(J)+J
         SUM=0.0
         DO 30 K=MAX(IND1,MU1JPJ-IMA(J)),J-1
         SUM=SUM+ASS(MU1IMI+K)*ASS(MU1JPJ-K)
   30    CONTINUE
         ASS(MU1IMI+J)=ASS(MU1IMI+J)-SUM
         SUM=0.0
         IF(J.GT.1) THEN
            DO 40 K=MAX(IND2,IMA(J-1)-MU1JMJ+1),J-1
            SUM=SUM+ASS(MU1JMJ+K)*ASS(MU1IPI-K)
   40       CONTINUE
         ENDIF
         ASS(MU1IPI-J)=(ASS(MU1IPI-J)-SUM)/ASS(MU1JMJ+J)
   50    CONTINUE
      ELSE
         DO 70 J=IND2,IND1-1
         MU1JMJ=MU1(J)-J
         SUM=0.0
         IF(J.GT.1) THEN
            DO 60 K=MAX(IND2,IMA(J-1)-MU1JMJ+1),J-1
            SUM=SUM+ASS(MU1JMJ+K)*ASS(MU1IPI-K)
   60       CONTINUE
         ENDIF
         ASS(MU1IPI-J)=(ASS(MU1IPI-J)-SUM)/ASS(MU1JMJ+J)
   70    CONTINUE
         DO 100 J=IND1,I-1
         MU1JMJ=MU1(J)-J
         MU1JPJ=MU1(J)+J
         SUM=0.0
         DO 80 K=MAX(IND1,MU1JPJ-IMA(J)),J-1
         SUM=SUM+ASS(MU1IMI+K)*ASS(MU1JPJ-K)
   80    CONTINUE
         ASS(MU1IMI+J)=ASS(MU1IMI+J)-SUM
         SUM=0.0
         DO 90 K=MAX(IND2,IMA(J-1)-MU1JMJ+1),J-1
         SUM=SUM+ASS(MU1JMJ+K)*ASS(MU1IPI-K)
   90    CONTINUE
         ASS(MU1IPI-J)=(ASS(MU1IPI-J)-SUM)/ASS(MU1JMJ+J)
  100    CONTINUE
      ENDIF
      SUM=0.0
      DO 110 K=MAX(IND1,IND2),I-1
      SUM=SUM+ASS(MU1IMI+K)*ASS(MU1IPI-K)
  110 CONTINUE
      ASS(MU1IMI+I)=ASS(MU1IMI+I)-SUM
  120 CONTINUE
      RETURN
      END
