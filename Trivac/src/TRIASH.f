*DECK TRIASH
      SUBROUTINE TRIASH(IELEM,NBMIX,LL4F,NBLOS,MAT,SIDE,ZZ,FRZ,SGD,KN,
     > IPERT,VEC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a diagonal system matrix corresponding to a single cross
* section type (Thomas-Raviart-Schneider dual cases). Note: vector VEC
* should be initialized by the calling program.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IELEM   degree of the Lagrangian finite elements.
* NBMIX   maximum number of material mixtures.
* LL4F    total number of flux unknowns per group.
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* MAT     mixture index assigned to each element.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* FRZ     volume fractions for the axial SYME boundary condition.
* SGD     cross section per material mixture.
* KN      ADI permutation indices for the volumes.
* IPERT   mixture permutation index.
*
*Parameters: output
* VEC     diagonal system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER(MAXIEL=3)
      INTEGER IELEM,NBMIX,LL4F,NBLOS,MAT(3,NBLOS),KN(NBLOS,3),
     1 IPERT(NBLOS)
      REAL SIDE,ZZ(3,NBLOS),FRZ(NBLOS),SGD(NBMIX),VEC(LL4F)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT,VOL0,SIG
*
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 20 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 20
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 20
      NUM=NUM+1
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      SIG=SGD(IBM)
      DO 12 K3=0,IELEM-1
      DO 11 K2=0,IELEM-1
      DO 10 K1=0,IELEM-1
      JND1=(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      JND2=(KN(NUM,1)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      JND3=(KN(NUM,2)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      VEC(JND1)=VEC(JND1)+REAL(VOL0*SIG)
      VEC(JND2)=VEC(JND2)+REAL(VOL0*SIG)
      VEC(JND3)=VEC(JND3)+REAL(VOL0*SIG)
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
   20 CONTINUE
      RETURN
      END
