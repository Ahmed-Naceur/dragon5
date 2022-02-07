*DECK TRINDX
      SUBROUTINE TRINDX(I,IP,MAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the perdue storage indices.
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
* I     index of the new information element.
* IP    array of perdue storage indices.
* MAX   size of array IP.
*
*Parameters: output
* IP    array of perdue storage indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER I,MAX,IP(MAX)
*
      DO 10 I0=1,MAX
      IF(IP(I0).EQ.I) THEN
         RETURN
      ELSE IF(IP(I0).EQ.0) THEN
         IP(I0)=I
         RETURN
      ENDIF
   10 CONTINUE
      CALL XABORT('TRINDX: INDEX SEARCH FAILURE.')
      RETURN
      END
