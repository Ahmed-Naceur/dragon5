*DECK TRIAHP
      SUBROUTINE TRIAHP (MAXKN,ISPLH,IR,NEL,LL4,SGD,SIDE,ZZ,VOL,MAT,KN,
     1 R,RH,RT,VEC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a diagonal system matrix corresponding to a single cross
* section type (primal formulation) in hexagonal geometry.
* Note: vector VEC should be initialized by the calling program.
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
* ISPLH   type of mesh-splitting: =1 for complete hexagons; .gt.1 for
*         triangle mesh-splitting.
* IR      number of material mixtures.
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
*
*Parameters: output
* VEC     diagonal matrix corresponding to the cross section term.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXKN,ISPLH,IR,NEL,LL4,MAT(NEL),KN(MAXKN)
      REAL SGD(IR),SIDE,ZZ(NEL),VOL(NEL),R(2,2),RH(6,6),RT(3,3),VEC(LL4)
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
*
* COMPUTE THE HEXAGONAL MASS (RH2).
      IF(ISPLH.EQ.1) THEN
         LC=6
         DO 10 I=1,2*LC
         IJ1(I)=IJ16(I)
         IJ2(I)=IJ26(I)
   10    CONTINUE
         DO 30 I=1,LC
         DO 20 J=1,LC
         RH2(I,J)=RH(I,J)
   20    CONTINUE
   30    CONTINUE
      ELSE
         LC=7
         DO 40 I=1,2*LC
         IJ1(I)=IJ17(I)
         IJ2(I)=IJ27(I)
   40    CONTINUE
         DO 60 I=1,LC
         DO 50 J=1,LC
         RH2(I,J)=0.0
   50    CONTINUE
   60    CONTINUE
         DO 85 K=1,6
         DO 80 I=1,3
         NUMI=ILIEN(K,I)
         DO 70 J=1,3
         NUMJ=ILIEN(K,J)
         RH2(NUMI,NUMJ)=RH2(NUMI,NUMJ)+RT(I,J)
   70    CONTINUE
   80    CONTINUE
   85    CONTINUE
      ENDIF
      LL=2*LC
*
* CALCULATION OF 3-D MASS AND STIFFNESS MATRICES FROM TENSORIAL PRODUCT
* OF 1-D AND 2-D MATRICES.
      DO 100 I=1,LL
      I1=IJ1(I)
      I2=IJ2(I)
      DO 90 J=1,LL
      J1=IJ1(J)
      J2=IJ2(J)
      RTHG(I,J)=RH2(I1,J1)*R(I2,J2)
   90 CONTINUE
  100 CONTINUE
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
      VEC(INW1)=VEC(INW1)+REAL(RR*SGD(L))
  110 CONTINUE
  150 NUM1=NUM1+LL
  160 CONTINUE
      RETURN
      END
