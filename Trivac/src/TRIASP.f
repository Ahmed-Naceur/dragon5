*DECK TRIASP
      SUBROUTINE TRIASP (IELEM,IR,NEL,LL4,CYLIND,SGD,XX,DD,VOL,MAT,KN,
     1 LC,T,TS,VEC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a diagonal system matrix corresponding to a single cross
* section type (primal formulation). Note: vector VEC should be
* initialized by the calling program.
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
* IELEM   degree of the Lagrangian finite elements.
* IR      number of material mixtures.
* NEL     total number of finite elements.
* ll4     order of system matrices.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry).
* SGD     cross section per material mixture.
* XX      X-directed mesh spacings.
* DD      used with cylindrical geometry.
* VOL     volume of each element.
* MAT     mixture index assigned to each element.
* KN      element-ordered unknown list.
* LC      order of the unit matrices.
* T       Cartesian linear product vector.
* TS      cylindrical linear product vector.
*
*Parameters: output
* VEC     diagonal matrix corresponding to the cross section term.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,IR,NEL,LL4,MAT(NEL),KN(NEL*(IELEM+1)**3),LC
      REAL SGD(IR),XX(NEL),DD(NEL),VOL(NEL),T(LC),TS(LC),VEC(LL4)
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
*
      NUM1=0
      DO 90 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 90
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 80
      DX=XX(K)
      DO 50 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 50
      IF(CYLIND) THEN
         RR=(R3DP(I)+R3DC(I)*DX/DD(K))*VOL0
      ELSE
         RR=R3DP(I)*VOL0
      ENDIF
      VEC(IND1)=VEC(IND1)+RR*SGD(L)
   50 CONTINUE
   80 NUM1=NUM1+LL
   90 CONTINUE
      RETURN
      END
