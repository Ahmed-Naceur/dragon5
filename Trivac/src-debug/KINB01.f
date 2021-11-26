*DECK KINB01
      SUBROUTINE KINB01(MAXKN,SGD,CYLIND,NREG,LL4,NBMIX,XX,DD,MAT,KN,
     1 VOL,LC,R,RS,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in primal finite element
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
* MAXKN   dimension of array KN.
* SGD     mixture-ordered cross sections.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry).
* NREG    number of elements in Bivac.
* LL4     order of matrix SYS.
* NBMIX   number of macro-mixtures.
* XX      X-directed mesh spacings.
* DD      value used with a cylindrical geometry.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* VOL     volume of regions.
* LC      number of polynomials in a complete 1-D basis.
* R       Cartesian mass matrix.
* RS      cylindrical mass matrix.
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
      INTEGER MAXKN,NREG,LL4,NBMIX,MAT(NREG),KN(MAXKN),LC
      REAL SGD(NBMIX),XX(NREG),DD(NREG),VOL(NREG),R(LC,LC),RS(LC,LC),
     1 F2(LL4),F3(LL4)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(25),IJ2(25)
      REAL R2DP(25,25),R2DC(25,25)
*----
*  COMPUTE VECTORS IJ1 AND IJ2.
*----
      LL=LC*LC
      DO 10 I=1,LL
      IJ1(I)=1+MOD(I-1,LC)
      IJ2(I)=1+(I-IJ1(I))/LC
   10 CONTINUE
*----
*  COMPUTE THE CARTESIAN 2-D MASS MATRICES FROM TENSORIAL PRODUCTS OF
*  1-D MATRICES.
*----
      DO 25 I=1,LL
      I1=IJ1(I)
      I2=IJ2(I)
      DO 20 J=1,LL
      J1=IJ1(J)
      J2=IJ2(J)
      R2DP(I,J)=R(I1,J1)*R(I2,J2)
      R2DC(I,J)=RS(I1,J1)*R(I2,J2)
   20 CONTINUE
   25 CONTINUE
*----
*  MULTIPLICATION.
*----
      NUM1=0
      DO 60 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 60
      IF(VOL(K).EQ.0.0) GO TO 50
      DX=XX(K)
      DO 40 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 40
      DO 30 J=1,LL
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 30
      IF(CYLIND) THEN
        RR=R2DP(I,J)+R2DC(I,J)*DX/DD(K)
      ELSE
        RR=R2DP(I,J)
      ENDIF
      IF(RR.EQ.0.0) GO TO 30
      F3(IND1)=F3(IND1)+RR*SGD(L)*VOL(K)*F2(IND2)
   30 CONTINUE
   40 CONTINUE
   50 NUM1=NUM1+LL
   60 CONTINUE
      RETURN
      END
