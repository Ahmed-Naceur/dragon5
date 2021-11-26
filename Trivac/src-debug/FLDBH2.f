*DECK FLDBH2
      SUBROUTINE FLDBH2 (ISPLH,NEL,NUN,NELEM,EVECT,VOL,IDL,KN,QFR,RH,RT)
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
*         >1 for triangular mesh-splitting.
* NEL     number of hexagons.
* NUN     number of unknowns per energy group.
* NELEM   number of finite elements (hexagons or triangles) excluding
*         the virtual elements.
* EVECT   variational coefficients of the flux. The information is
*         contained in position EVECT(1) to EVECT(LL4) where LL4 is the
*         order of the system matrices.
* VOL     volume of each hexagon.
* IDL     position of the average flux component associated with each
*         hexagon.
* KN      element-ordered unknown list. The dimension of KN is equal
*         to (LC+1)*NELEM where LC=6 (hexagons) or 3 (triangles).
* QFR     element-ordered albedo information. The dimension of QFR is
*         equal to (LC+1)*NELEM.
* RH      unit matrix
* RT      unit matrix
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
      INTEGER ISPLH,NEL,NUN,NELEM,IDL(NEL),KN(*)
      REAL EVECT(NUN),VOL(NEL),QFR(*),RH(6,6),RT(3,3)
*----
*  LOCAL VARIABLES
*----
      REAL T(6)
*----
*  COMPUTE THE LINEAR PRODUCT VECTOR T
*----
      IF(ISPLH.EQ.1) THEN
*        HEXAGONAL BASIS.
         LC=6
         DO 15 I=1,6
         T(I)=0.0
         DO 10 J=1,6
         T(I)=T(I)+RH(I,J)
   10    CONTINUE
   15    CONTINUE
         CONST=1.5*SQRT(3.0)
      ELSE
*        TRIANGULAR BASIS.
         LC=3
         DO 25 I=1,3
         T(I)=0.0
         DO 20 J=1,3
         T(I)=T(I)+RT(I,J)
   20    CONTINUE
   25    CONTINUE
         CONST=0.25*SQRT(3.0)
      ENDIF
*
      DO 30 KHEX=1,NEL
      IF(IDL(KHEX).NE.0) EVECT(IDL(KHEX))=0.0
   30 CONTINUE
      NUM1=0
      DO 60 K=1,NELEM
      KHEX=KN(NUM1+LC+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 50
      DO 40 I=1,LC
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 40
      SS=T(I)*QFR(NUM1+LC+1)/(CONST*VOL(KHEX))
      EVECT(IDL(KHEX))=EVECT(IDL(KHEX))+SS*EVECT(IND1)
   40 CONTINUE
   50 NUM1=NUM1+LC+1
   60 CONTINUE
      RETURN
      END
