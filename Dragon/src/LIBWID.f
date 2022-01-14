*DECK LIBWID
      FUNCTION LIBWID(NEL,IWISO,ISOID)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find isotope number associated with isotope id on WIMS-NEA and
* WIMS-AECL library.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* NEL     number of isotopes on library.
* IWISO   id of isotope on library.
* ISOID   isotope id requested.
*
*Parameters: output
* LIBWID  isotope number associated with id  (= 0 when isotope id is
*         not found).
*
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    NEL,IWISO(NEL),ISOID,LIBWID
*----
*  LOCAL VARIABLES
*----
      INTEGER    IEL
*----
      LIBWID=0
      IF(ISOID.GT.0) THEN
        DO 100 IEL=1,NEL
          IF(IWISO(IEL).EQ.ISOID) THEN
            LIBWID=IEL
            GO TO 105
          ENDIF
 100    CONTINUE
 105    CONTINUE
      ENDIF
      RETURN
      END
