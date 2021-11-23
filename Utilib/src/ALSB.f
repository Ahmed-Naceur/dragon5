*DECK ALSB
      SUBROUTINE ALSB (N,IS,B,IER,MAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* solution of a system of linear equations using gaussian elimination
* with partial pivoting. Simple precision version.
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
* N       order of the coefficient matrix.
* IS      number of right hand vectors.
* B       coefficient matrix augmented with the right hand vectors.
*         DIMENSION B(MAX,N+IS)
* MAX     first dimention of matrix B.
*
*Parameters: output
* B       solution vectors, starting at B(1,N+1).
* IER     error flag. Execution failure if IER.ne.0.
*
*-----------------------------------------------------------------------
*
      DIMENSION B(MAX,*)
      IN=0
      M=N+IS
      IER=0
      IF (N.EQ.1) GO TO 100
*
* SEARCH FOR MAXIMUM PIVOT ON COLUMN JCOL.
      NM1=N-1
      NP1=N+1
      DO 60 JCOL=1,NM1
      TEST=0.0
      DO 10 I=JCOL,N
      IF (ABS(B(I,JCOL)).LE.TEST) GO TO 10
      TEST=ABS(B(I,JCOL))
      IN=I
10    CONTINUE
      IF (TEST.EQ.0.0) GO TO 120
*
* TRIANGULARIZATION.
      PMX=B(IN,JCOL)
      B(IN,JCOL)=B(JCOL,JCOL)
      IP1=JCOL+1
      DO 50 J=IP1,M
      PER=B(IN,J)/PMX
      B(IN,J)=B(JCOL,J)
      B(JCOL,J)=PER
      DO 40 I=IP1,N
      B(I,J)=B(I,J)-B(I,JCOL)*PER
40    CONTINUE
50    CONTINUE
60    CONTINUE
      PER=B(N,N)
      IF (PER.EQ.0.0) GO TO 120
      DO 70 J=NP1,M
      B(N,J)=B(N,J)/PER
70    CONTINUE
*
* BACK SUBSTITUTION.
      DO 95 IN=2,N
      I=N-IN+1
      IP1=I+1
      DO 90 J=NP1,M
      PER=B(I,J)
      DO 80 K=IP1,N
      PER=PER-B(I,K)*B(K,J)
80    CONTINUE
      B(I,J)=PER
90    CONTINUE
95    CONTINUE
      RETURN
*
100   PER=B(1,1)
      IF (PER.EQ.0.0) GO TO 120
      DO 110 J=2,M
      B(1,J)=B(1,J)/PER
110   CONTINUE
      RETURN
120   IER=1
      RETURN
      END
