*DECK FLDTH2
      SUBROUTINE FLDTH2 (ISPLH,NEL,NUN,EVECT,MAT,VOL,IDL,KN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the averaged flux with a linear Lagrangian finite
* element or mesh corner finite difference method in hexagonal geometry.
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
* ISPLH   type of hexagonal mesh-splitting: =1 for complete hexagons;
*         =2 for triangular mesh-splitting.
* NEL     total number of finite elements.
* NUN     total number of unknowns per energy group.
* EVECT   variational coefficients of the flux. The information is
*         contained in position EVECT(1) to EVECT(LL4) where LL4 is
*         the order of the system matrices.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* IDL     position of the average flux component associated with each
*         volume.
* KN      element-ordered unknown list.
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
      INTEGER ISPLH,NEL,NUN,MAT(NEL),IDL(NEL),KN(14*NEL)
      REAL EVECT(NUN),VOL(NEL)
*----
*  LOCAL VARIABLES
*----
      REAL TH(14)
      SAVE TH
      DATA TH/3*0.055555555556,0.166666666667,6*0.055555555556,
     1 0.166666666667,3*0.055555555556/
*
      NUM1=0
      IF(ISPLH.EQ.1) THEN
         SS=1.0/12.0
         DO 40 K=1,NEL
         IF(MAT(K).EQ.0) GO TO 40
         EVECT(IDL(K))=0.0
         IF(VOL(K).EQ.0.0) GO TO 30
         DO 20 I=1,12
         IND1=KN(NUM1+I)
         IF(IND1.EQ.0) GO TO 20
         EVECT(IDL(K))=EVECT(IDL(K))+SS*EVECT(IND1)
   20    CONTINUE
   30    NUM1=NUM1+12
   40    CONTINUE
      ELSE IF(ISPLH.EQ.2) THEN
         DO 70 K=1,NEL
         IF(MAT(K).EQ.0) GO TO 70
         EVECT(IDL(K))=0.0
         IF(VOL(K).EQ.0.0) GO TO 60
         DO 50 I=1,14
         IND1=KN(NUM1+I)
         IF(IND1.EQ.0) GO TO 50
         EVECT(IDL(K))=EVECT(IDL(K))+TH(I)*EVECT(IND1)
   50    CONTINUE
   60    NUM1=NUM1+14
   70    CONTINUE
      ENDIF
      RETURN
      END
