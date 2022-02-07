*DECK KINB02
      SUBROUTINE KINB02(SGD,IELEM,NREG,LL4,NBMIX,MAT,KN,VOL,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mixed-dual finite element
* diffusion approximation (Cartesian geometry). Special version for
* Bivac.
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
*Parameters: input
* SGD     mixture-ordered cross sections.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NREG    number of elements in Bivac.
* LL4     number of unknowns per group in Bivac.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* VOL     volume of regions.
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
      INTEGER IELEM,NREG,LL4,NBMIX,MAT(NREG),KN(5*NREG)
      REAL SGD(NBMIX),VOL(NREG),F2(LL4),F3(LL4)
*
      NUM1=0
      DO 30 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 30
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 20
      DO 15 I0=1,IELEM
      DO 10 J0=1,IELEM
      JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
      F3(JND1)=F3(JND1)+VOL0*SGD(L)*F2(JND1)
   10 CONTINUE
   15 CONTINUE
   20 NUM1=NUM1+5
   30 CONTINUE
      RETURN
      END
