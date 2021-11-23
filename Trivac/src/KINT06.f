*DECK KINT06
      SUBROUTINE KINT06(ISPLH,NBMIX,NEL,LL4,VOL,MAT,SGD,KN,IPW,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mesh centered finite
* difference approximation (hexagonal geometry with triangular
* submeshes). Special version for Trivac.
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
* ISPLH   related to the triangular submesh. The number of triangles is
*         6*(ISPLH-1)**2.
* NBMIX   maximum number of material mixtures.
* NEL     total number of finite elements.
* LL4     order of the system matrices.
* VOL     volume of each element.
* MAT     mixture index assigned to each hexagon.
* SGD     nuclear properties per material mixtures.
* KN      element-ordered unknown list.
* IPW     permutation matrices.
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
      INTEGER ISPLH,NBMIX,NEL,LL4,MAT(NEL),KN(NEL*(18*(ISPLH-1)**2+8)),
     1 IPW(LL4)
      REAL VOL(NEL),SGD(NBMIX),F2(LL4),F3(LL4)
*----
*  MULTIPLICATION
*----
      NUM1 = 0
      NTPH = 6 * (ISPLH-1)**2
      NTPL = 1 + 2 * (ISPLH-1)
      NVT1 = NTPL + 2 * (ISPLH-2) + NTPH / 2
      NVT2 = NTPH - NTPL - (ISPLH-4) * (NTPL+2)
      NVT3 = NTPH - (ISPLH-4) * NTPL
      IVAL = 3*NTPH+8
      IF(ISPLH.EQ.3) NVT2 = NTPH
      IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
      IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
      ICR = ISAU*(1+2*(ISPLH-2))
      DO 40 K=1,NEL
         L = MAT(K)
         IF(L.EQ.0) GO TO 40
         VOL0 = VOL(K)/NTPH
         IF(VOL0.EQ.0.0) GO TO 30
         DO 20 I = 1,NTPH
*
            CALL TRINEI (3,1,1,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,
     >                   NTPH,NTPL,NVT1,NVT2,NVT3,IVAL,KN)
*
            IND1=IPW(KEL)
            F3(IND1)=F3(IND1)+SGD(L)*VOL0*F2(IND1)
   20    CONTINUE
   30    NUM1=NUM1+IVAL
   40 CONTINUE
      RETURN
      END
