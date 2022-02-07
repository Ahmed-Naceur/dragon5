*DECK KINT05
      SUBROUTINE KINT05(NBMIX,NEL,LL4,SGD,VOL,MAT,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mesh centered finite
* difference approximation (hexagonal geometry). Special version for
* Trivac.
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
* NBMIX   maximum number of material mixtures.
* NEL     total number of finite elements.
* LL4     number of unknowns (order of the system matrices).
* SGD     cross section per material mixture.
* VOL     volumes.
* MAT     index-number of the mixture type assigned to each volume.
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
      INTEGER NBMIX,NEL,LL4,MAT(NEL)
      REAL SGD(NBMIX),VOL(NEL),F2(LL4),F3(LL4)
*
      KEL=0
      DO 10 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 10
      KEL=KEL+1
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 10
      F3(KEL)=F3(KEL)+SGD(L)*VOL0*F2(KEL)
   10 CONTINUE
      RETURN
      END
