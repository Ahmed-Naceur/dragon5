*DECK FLDTN2
      SUBROUTINE FLDTN2 (NEL,LL4,IELEM,CYLIND,EVECT,XX,DD,MAT,VOL,IDL,
     1 KN,LC,T,TS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the integrated flux in each finite element for a Lagrangian
* finite element discretization in Cartesian or cylindrical geometry.
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
* NEL     number of finite elements.
* LL4     order of system matrices.
* IELEM   degree of the finite elements.
* CYLIND  cylindrical geometry flag (set with CYLIND=.true.).
* EVECT   unknown vector containing the variational coefficients in
*         locations 1 to LL4.
* XX      X-directed mesh spacings.
* DD      used with cylindrical geometry.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* IDL     indices pointing to integrated fluxes in EVECT array.
* KN      element-ordered unknown list.
* LC      order of the finite element basis.
* T       linear product vector.
* TS      linear product vector.
*
*Parameters: output
* EVECT   unknown vector containing the integrated fluxes in location
*         IDL(I).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NEL,LL4,IELEM,MAT(NEL),IDL(NEL),KN(NEL*(IELEM+1)**3),LC
      REAL EVECT(LL4+NEL),XX(NEL),DD(NEL),VOL(NEL),T(LC),TS(LC)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(125),IJ2(125),IJ3(125)
*
      LL=LC*LC*LC
      DO 100 L=1,LL
      L1=1+MOD(L-1,LC)
      L2=1+(L-L1)/LC
      L3=1+MOD(L2-1,LC)
      IJ1(L)=L1
      IJ2(L)=L3
      IJ3(L)=1+(L2-L3)/LC
  100 CONTINUE
*
      NUM1=0
      DO 130 K=1,NEL
      IF(MAT(K).EQ.0) GO TO 130
      EVECT(IDL(K))=0.0
      IF(VOL(K).EQ.0.0) GO TO 120
      DO 110 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 110
      I1=IJ1(I)
      I2=IJ2(I)
      I3=IJ3(I)
      IF(CYLIND) THEN
         SS=(T(I1)+TS(I1)*XX(K)/DD(K))*T(I2)*T(I3)
      ELSE
         SS=T(I1)*T(I2)*T(I3)
      ENDIF
      EVECT(IDL(K))=EVECT(IDL(K))+SS*EVECT(IND1)
  110 CONTINUE
  120 NUM1=NUM1+LL
  130 CONTINUE
      RETURN
      END
