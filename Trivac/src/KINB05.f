*DECK KINB05
      SUBROUTINE KINB05(SGD,IELEM,NBLOS,LL4,NBMIX,SIDE,MAT,IPERT,
     1 KN,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in Thomas-Raviart-Schneider
* (dual) finite element diffusion approximation (hexagonal geometry).
* Special version for Bivac.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
* SGD     mixture-ordered cross sections.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* LL4     number of unknowns per group in Bivac.
* NBMIX   number of macro-mixtures.
* SIDE    side of the hexagons.
* MAT     mixture index per region.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* F2      vector to multiply.
*
*Parameters: output
* F3      result of the multiplication.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,NBLOS,LL4,NBMIX,MAT(3,NBLOS),IPERT(NBLOS),
     1 KN(NBLOS,4+6*IELEM*(IELEM+1))
      REAL SGD(NBMIX),SIDE,F2(LL4),F3(LL4)
*----
*  ASSEMBLY OF A SYSTEM MATRIX.
*----
      TTTT=0.5*SQRT(3.0)*SIDE*SIDE
      NUM=0
      DO 20 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 20
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 20
      SIG=SGD(IBM)
      DO 15 K2=0,IELEM-1
      DO 10 K1=0,IELEM-1
      JND1=KN(NUM,1)+K2*IELEM+K1
      JND2=KN(NUM,2)+K2*IELEM+K1
      JND3=KN(NUM,3)+K2*IELEM+K1
      F3(JND1)=F3(JND1)+TTTT*SIG*F2(JND1)
      F3(JND2)=F3(JND2)+TTTT*SIG*F2(JND2)
      F3(JND3)=F3(JND3)+TTTT*SIG*F2(JND3)
   10 CONTINUE
   15 CONTINUE
   20 CONTINUE
      RETURN
      END
