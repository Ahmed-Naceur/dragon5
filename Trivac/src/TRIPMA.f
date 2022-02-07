*DECK TRIPMA
      SUBROUTINE TRIPMA(LC,T,TS,Q,QS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the unit matrices for a primal finite element method in 3-D.
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
* LC      order of the unit matrices.
* T       cartesian linear product vector.
* TS      cylindrical linear product vector.
* Q       cartesian stiffness matrix.
* QS      cylindrical stiffness matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE PARAMETERS
*----
      INTEGER LC,IJ1,IJ2,IJ3,ISR
      REAL T(LC),TS(LC),Q(LC,LC),QS(LC,LC)
*----
*  LOCAL VARIABLES
*----
      COMMON /ELEM2/LL,LCC,IJ1(125),IJ2(125),IJ3(125),ISR(6,25),
     1 Q3DP1(125,125),Q3DP2(125,125),Q3DP3(125,125),R3DP(125),
     2 Q3DC1(125,125),Q3DC2(125,125),Q3DC3(125,125),R3DC(125),
     3 R2DP(25),R2DC(25)
      SAVE /ELEM2/
*----
*  CALCULATION OF COMMON /ELEM2/
*----
      LCC=LC*LC
      LL=LC*LC*LC
      DO 5 L=1,LL
      L1=1+MOD(L-1,LC)
      L2=1+(L-L1)/LC
      L3=1+MOD(L2-1,LC)
      IJ1(L)=L1
      IJ2(L)=L3
      IJ3(L)=1+(L2-L3)/LC
    5 CONTINUE
*----
*  CALCULATION OF MATRIX ISR.
*----
      K1=0
      K2=0
      J1=0
      J2=0
      I1=0
      I2=0
      L=0
      DO 8 I=1,LC
      DO 7 J=1,LC
      DO 6 K=1,LC
      L=L+1
      IF(K.EQ.1) THEN
         K1=K1+1
         ISR(1,K1)=L
      ELSE IF(K.EQ.LC) THEN
         K2=K2+1
         ISR(2,K2)=L
      ENDIF
      IF(J.EQ.1) THEN
         J1=J1+1
         ISR(3,J1)=L
      ELSE IF(J.EQ.LC) THEN
         J2=J2+1
         ISR(4,J2)=L
      ENDIF
      IF(I.EQ.1) THEN
         I1=I1+1
         ISR(5,I1)=L
      ELSE IF(I.EQ.LC) THEN
         I2=I2+1
         ISR(6,I2)=L
      ENDIF
    6 CONTINUE
    7 CONTINUE
    8 CONTINUE
*----
*  CALCULATION OF 3-D MASS AND STIFFNESS MATRICES FROM TENSORIAL PRODUCT
*  OF 1-D MATRICES.
*----
      DO 20 I=1,LL
      I1=IJ1(I)
      I2=IJ2(I)
      I3=IJ3(I)
      DO 10 J=1,LL
      J1=IJ1(J)
      J2=IJ2(J)
      J3=IJ3(J)
      IF((I2.EQ.J2).AND.(I3.EQ.J3)) THEN
         Q3DP1(I,J)=Q(I1,J1)*T(I2)*T(I3)
         Q3DC1(I,J)=QS(I1,J1)*T(I2)*T(I3)
      ELSE
         Q3DP1(I,J)=0.0
         Q3DC1(I,J)=0.0
      ENDIF
      IF((I1.EQ.J1).AND.(I3.EQ.J3)) THEN
         Q3DP2(I,J)=T(I1)*Q(I2,J2)*T(I3)
         Q3DC2(I,J)=TS(I1)*Q(I2,J2)*T(I3)
      ELSE
         Q3DP2(I,J)=0.0
         Q3DC2(I,J)=0.0
      ENDIF
      IF((I1.EQ.J1).AND.(I2.EQ.J2)) THEN
         Q3DP3(I,J)=T(I1)*T(I2)*Q(I3,J3)
         Q3DC3(I,J)=TS(I1)*T(I2)*Q(I3,J3)
      ELSE
         Q3DP3(I,J)=0.0
         Q3DC3(I,J)=0.0
      ENDIF
   10 CONTINUE
      R3DP(I)=T(I1)*T(I2)*T(I3)
      R3DC(I)=TS(I1)*T(I2)*T(I3)
   20 CONTINUE
*----
*  CALCULATION OF 2-D MASS MATRICES FROM TENSORIAL PRODUCT OF 1-D
*  MATRICES.
*----
      DO 30 I=1,LC*LC
      I1=IJ1(I)
      I2=IJ2(I)
      R2DP(I)=T(I1)*T(I2)
      R2DC(I)=TS(I1)*T(I2)
   30 CONTINUE
      RETURN
      END
