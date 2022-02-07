*DECK FLDTH1
      SUBROUTINE FLDTH1 (ISPLH,NEL,LL4,EVECT,MAT,VOL,IDL,KN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the averaged flux for a linear primal formulation of
* the diffusion equation in hexagonal geometry with triangular
* mesh-splitting.
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
* ISPLH   type of triangular mesh-splitting (ISPLH.GT.1).
* NEL     total number of finite elements.
* LL4     order of the system matrices.
* EVECT   variational coefficients of the flux. The information is
*         contained in position EVECT(1) to EVECT(LL4).
* MAT     mixture index assigned to each element.
* VOL     volume of each element
* IDL     position of the average flux component associated with each
*         volume.
* KN      element-ordered unknown list.
*
*Parameters: output
* EVECT   averaged fluxes. The information is contained in positions
*         EVECT(LL4+1) to EVECT(LL4+NEL).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    ISPLH,NEL,LL4,MAT(NEL),IDL(NEL),
     1           KN(NEL*(18*(ISPLH-1)**2+8))
      REAL       EVECT(LL4+NEL),VOL(NEL)
*
      IVAL=18*(ISPLH-1)**2+8
      NUM1=0
      SS=1.0/REAL(6*(ISPLH-1)**2)
      DO 70 K=1,NEL
      IF(MAT(K).EQ.0) GO TO 70
      EVECT(IDL(K))=0.0
      IF(VOL(K).EQ.0.0) GO TO 60
      DO 50 I=1,6*(ISPLH-1)**2
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 50
      EVECT(IDL(K))=EVECT(IDL(K))+SS*EVECT(IND1)
   50 CONTINUE
   60 NUM1=NUM1+IVAL
   70 CONTINUE
      RETURN
      END
