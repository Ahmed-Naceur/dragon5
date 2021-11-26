*DECK FLDBH1
      SUBROUTINE FLDBH1 (NEL,NUN,LL4,EVECT,VOL,IDL,KN,QFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the averaged flux with a mesh centered finite
* difference method in hexagonal geometry with triangular
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
* NEL     number of hexagons.
* NUN     number of unknowns per energy group.
* LL4     order of the system matrices.
* EVECT   variational coefficients of the flux. The information is
*         contained in position EVECT(1) to EVECT(LL4).
* VOL     volume of each element.
* IDL     position of the average flux component associated with each
*         volume.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
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
      INTEGER NEL,NUN,LL4,IDL(NEL),KN(4*LL4)
      REAL EVECT(NUN),VOL(NEL),QFR(4*LL4)
*
      NSURF=3
      DO 10 K=1,NEL
      IF(IDL(K).NE.0) EVECT(IDL(K))=0.0
   10 CONTINUE
      NUM1=0
      DO 20 IND1=1,LL4
      K=KN(NUM1+NSURF+1)
      EVECT(IDL(K))=EVECT(IDL(K))+QFR(NUM1+NSURF+1)*EVECT(IND1)/VOL(K)
      NUM1=NUM1+NSURF+1
   20 CONTINUE
      RETURN
      END
