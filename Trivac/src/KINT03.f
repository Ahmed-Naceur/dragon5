*DECK KINT03
      SUBROUTINE KINT03(MAXKN,ISPLH,NBMIX,NEL,LL4,SGD,SIDE,ZZ,VOL,MAT,
     1 KN,R,RH,RT,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mesh-corner finite
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
* ISPLH   type of mesh-splitting: =1 for complete hexagons; .gt.1 for
*         triangle mesh-splitting.
* NBMIX   number of material mixtures.
* NEL     total number of finite elements.
* LL4     order of system matrices.
* SGD     cross section per material mixture.
* SIDE    dide of an hexagon.
* ZZ      height of each hexagon.
* VOL     volume of each element.
* MAT     mixture index assigned to each element.
* KN      element-ordered unknown list.
* R       unit matrix.
* RH      unit matrix.
* RT      unit matrix.
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
      INTEGER MAXKN,ISPLH,NBMIX,NEL,LL4,MAT(NEL),KN(MAXKN)
      REAL SGD(NBMIX),SIDE,ZZ(NEL),VOL(NEL),R(2,2),RH(6,6),RT(3,3),
     1 F2(LL4),F3(LL4)
*----
*  LOCAL VARIABLES
*----
      INTEGER ILIEN(6,3),IJ17(14),IJ27(14),IJ16(12),IJ26(12),IJ1(14),
     1 IJ2(14)
      REAL RH2(7,7)
      DOUBLE PRECISION RR,VOL1,RTHG(14,14)
      DATA ILIEN/6*4,2,1,5,6,7,3,1,5,6,7,3,2/
      DATA IJ16,IJ26 /1,2,3,4,5,6,1,2,3,4,5,6,6*1,6*2/
      DATA IJ17,IJ27 /1,2,3,4,5,6,7,1,2,3,4,5,6,7,7*1,7*2/
*----
*  COMPUTE THE HEXAGONAL MASS (RH2).
*----
      IF(ISPLH.EQ.1) THEN
         LC=6
         DO 20 I=1,2*LC
         IJ1(I)=IJ16(I)
         IJ2(I)=IJ26(I)
   20    CONTINUE
         DO 41 I=1,LC
         DO 40 J=1,LC
         RH2(I,J)=RH(I,J)
   40    CONTINUE
   41    CONTINUE
      ELSE
         LC=7
         DO 60 I=1,2*LC
         IJ1(I)=IJ17(I)
         IJ2(I)=IJ27(I)
   60    CONTINUE
         DO 76 I=1,LC
         DO 75 J=1,LC
         RH2(I,J)=0.0
   75    CONTINUE
   76    CONTINUE
         DO 82 K=1,6
         DO 81 I=1,3
         NUMI=ILIEN(K,I)
         DO 80 J=1,3
         NUMJ=ILIEN(K,J)
         RH2(NUMI,NUMJ)=RH2(NUMI,NUMJ)+RT(I,J)
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
      ENDIF
      LL=2*LC
*----
*  CALCULATION OF 3-D MASS AND STIFFNESS MATRICES FROM TENSORIAL PRODUCT
*  OF 1-D AND 2-D MATRICES.
*----
      DO 91 I=1,LL
      I1=IJ1(I)
      I2=IJ2(I)
      DO 90 J=1,LL
      J1=IJ1(J)
      J2=IJ2(J)
      RTHG(I,J)=RH2(I1,J1)*R(I2,J2)
   90 CONTINUE
   91 CONTINUE
*
      NUM1=0
      VOL1=SIDE*SIDE
      DO 160 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 160
      IF(VOL(K).EQ.0.0) GO TO 150
      DO 110 I=1,LL
      INW1=KN(NUM1+I)
      IF(INW1.EQ.0) GO TO 110
      RR=RTHG(I,I)*VOL1*ZZ(K)
      F3(INW1)=F3(INW1)+REAL(RR)*SGD(L)*F2(INW1)
  110 CONTINUE
  150 NUM1=NUM1+LL
  160 CONTINUE
      RETURN
      END
