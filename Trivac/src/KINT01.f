*DECK KINT01
      SUBROUTINE KINT01(MAXKN,SGD,CYLIND,NREG,LL4,NBMIX,XX,DD,MAT,KN,
     1 VOL,LC,T,TS,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in primal finite element
* diffusion approximation (Cartesian geometry). Special version for
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
* MAXKN   dimension of array KN.
* SGD     mixture-ordered cross sections.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry).
* NREG    number of elements in TRIVAC.
* LL4     order of matrix SYS.
* NBMIX   number of macro-mixtures.
* XX      X-directed mesh spacings.
* DD      value used with a cylindrical geometry.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* VOL     volume of regions.
* LC      number of polynomials in a complete 1-D basis.
* T       Cartesian linear product vector.
* TS      cylindrical linear product vector.
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
      REAL SGD(NBMIX),XX(NREG),DD(NREG),VOL(NREG),T(LC),TS(LC),F2(LL4),
     1 F3(LL4)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      REAL R3DP(125),R3DC(125)
*----
*  CALCULATION OF 3-D MASS MATRICES FROM TENSORIAL PRODUCT OF 1-D
*  MATRICES
*----
      LL=LC*LC*LC
      DO 20 L=1,LL
      L1=1+MOD(L-1,LC)
      L2=1+(L-L1)/LC
      L3=1+MOD(L2-1,LC)
      I1=L1
      I2=L3
      I3=1+(L2-L3)/LC
      R3DP(L)=T(I1)*T(I2)*T(I3)
      R3DC(L)=TS(I1)*T(I2)*T(I3)
   20 CONTINUE
*----
*  MULTIPLICATION.
*----
      NUM1=0
      DO 90 K=1,NREG
      L=MAT(K)
      IF(L.EQ.0) GO TO 90
      IF(VOL(K).EQ.0.0) GO TO 80
      DX=XX(K)
      DO 50 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 50
      IF(CYLIND) THEN
         RR=(R3DP(I)+R3DC(I)*DX/DD(K))
      ELSE
         RR=R3DP(I)
      ENDIF
      F3(IND1)=F3(IND1)+RR*SGD(L)*VOL(K)*F2(IND1)
   50 CONTINUE
   80 NUM1=NUM1+LL
   90 CONTINUE
      RETURN
      END
