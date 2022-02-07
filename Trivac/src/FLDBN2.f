*DECK FLDBN2
      SUBROUTINE FLDBN2 (NEL,LL4,IELEM,CYLIND,EVECT,XX,DD,MAT,VOL,IDL,
     1 KN,LC,T,TS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the averaged flux with a primal finite element method.
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
* NEL     total number of finite elements.
* LL4     order of the system matrices.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* CYLIND  cylindrical geometry flag (set with CYLIND=.true.).
* EVECT   variational coefficients of the flux. The information is
*         contained in position EVECT(1) to EVECT(LL4).
* XX      X-side of each element.
* DD      used with cylindrical geometry.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* IDL     position of the average flux component associated with each
*         volume.
* KN      element-ordered unknown list.
* LC      order of the finite element basis.
* T       linear product vector.
* TS      linear product vector.
*
*Parameters: output
* EVECT   averaged fluxes. The information is contained in positions
*         EVECT(IDL(I)).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NEL,LL4,IELEM,MAT(NEL),IDL(NEL),KN(NEL*IELEM*IELEM),LC
      REAL EVECT(LL4+NEL),XX(NEL),DD(NEL),VOL(NEL),T(LC),TS(LC)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(25),IJ2(25)
*----
*  COMPUTE VECTORS IJ1 AND IJ2
*----
      LL=LC*LC
      DO 10 I=1,LL
      IJ1(I)=1+MOD(I-1,LC)
      IJ2(I)=1+(I-IJ1(I))/LC
   10 CONTINUE
*
      NUM1=0
      DO 40 K=1,NEL
      IF(MAT(K).EQ.0) GO TO 40
      EVECT(IDL(K))=0.0
      IF(VOL(K).EQ.0.0) GO TO 30
      DO 20 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 20
      I1=IJ1(I)
      I2=IJ2(I)
      IF(CYLIND) THEN
         SS=(T(I1)+TS(I1)*XX(K)/DD(K))*T(I2)
      ELSE
         SS=T(I1)*T(I2)
      ENDIF
      EVECT(IDL(K))=EVECT(IDL(K))+SS*EVECT(IND1)
   20 CONTINUE
   30 NUM1=NUM1+LL
   40 CONTINUE
      RETURN
      END
