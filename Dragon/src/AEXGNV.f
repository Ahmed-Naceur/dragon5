*DECK AEXGNV
      SUBROUTINE AEXGNV(ICHAMP,SEGMEN,TCHDIM,TCHTYP,TCHDKL,IDK,NV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Localization of a specific segment component in a SAPHYR archive.
* Component of a FORTRAN-77 emulator of the SAPHYR archive system.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ICHAMP  component index in segment (low integer value).
* SEGMEN  integer segment.
* TCHDIM  dimension of segment components.
* TCHTYP  character type of segment components.
* TCHDKL  position of segment components.
*
*Parameters: output
* IDK     position of segment component with index ichamp.
* NV      length of segment component with index ichamp.
*
*-----------------------------------------------------------------------
*
      INTEGER SEGMEN(*),TCHDIM(*),TCHTYP(*),TCHDKL(*)
      CHARACTER TEXT4*4
*
      LONMOT=4
      CALL LCMCAR(TEXT4,.FALSE.,TCHTYP(ICHAMP))
      ND=TCHDIM(ICHAMP)
      IDK=TCHDKL(ICHAMP+1)
      NV=1
      DO 100 ID=1,ND
      NV=NV*SEGMEN(IDK+ID)
  100 CONTINUE
      IDKE=TCHDKL(ICHAMP)
      IF(IDKE.LT.0) THEN
        IDK=SEGMEN(1-IDKE)+1
      ELSE
        IDK=IDKE+1
      ENDIF
      IF(TEXT4(1:1).EQ.'C') IDK=(IDK-1)/LONMOT+1
      RETURN
      END
