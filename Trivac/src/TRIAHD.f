*DECK TRIAHD
      SUBROUTINE TRIAHD (IR,NEL,LL4,SGD,VOL,MAT,VEC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a diagonal system matrix corresponding to a single cross
* section type. Mesh centered finite difference case in hexagonal
* geometry. Note: vector VEC should be initialized by the calling
* program.
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
* IR      maximum number of material mixtures.
* NEL     total number of finite elements.
* LL4     number of unknowns (order of the system matrices).
* SGD     cross section per material mixture.
* VOL     volumes.
* MAT     index-number of the mixture type assigned to each volume.
*
*Parameters: output
* VEC     diagonal system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL)
      REAL SGD(IR),VOL(NEL),VEC(LL4)
*
      KEL=0
      DO 10 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 10
      KEL=KEL+1
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 10
      VEC(KEL)=VEC(KEL)+SGD(L)*VOL0
   10 CONTINUE
      RETURN
      END
