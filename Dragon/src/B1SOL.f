*DECK B1SOL
      SUBROUTINE B1SOL(NGRO,B,IER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of a system of linear equations that appear in a B1 method.
* Use ALSBD.f for solution in thermal groups.
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
* NGRO    order of the coefficient matrix.
*
*Parameters: input/output
* B       coefficient matrix augmented with the right hand vector on
*         input and solution vector, starting at B(1,NGRO+1) at output.
*
*Parameters: output
* IER     error flag (execution failure if IER.ne.0).
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRO,IER
      DOUBLE PRECISION B(NGRO,NGRO+1)
*
      IER=0
      NGROIN=0
      DO 30 I=1,NGRO
      DO 10 J=I+1,NGRO
      IF(B(I,J).NE.0.0D0) GO TO 40
   10 CONTINUE
      NGROIN=I
      IF(B(I,I).EQ.0.0D0) THEN
        IER=-1
        RETURN
      ENDIF
      ZGAR=B(I,NGRO+1)
      DO 20 J=1,I-1
      ZGAR=ZGAR-B(I,J)*B(J,NGRO+1)
   20 CONTINUE
      B(I,NGRO+1)=ZGAR/B(I,I)
   30 CONTINUE
   40 IF(NGROIN.EQ.NGRO) RETURN
      DO 60 I=NGROIN+1,NGRO
      ZGAR=B(I,NGRO+1)
      DO 50 J=1,NGROIN
      ZGAR=ZGAR-B(I,J)*B(J,NGRO+1)
   50 CONTINUE
      B(I,NGRO+1)=ZGAR
   60 CONTINUE
      CALL ALSBD(NGRO-NGROIN,1,B(NGROIN+1,NGROIN+1),IER,NGRO)
      RETURN
      END
