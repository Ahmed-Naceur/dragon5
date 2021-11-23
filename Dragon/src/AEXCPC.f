*DECK AEXCPC
      SUBROUTINE AEXCPC(IDKSGT,LNGTAB,TABSGT,TABLUE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Extraction of character information from an integer segment. Component
* of a FORTRAN-77 emulator of the SAPHYR archive system.
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
* IDKSGT  position (in byte) of information in segment.
* LNGTAB  number of characters in the information.
* TABSGT  integer segment.
*
*Parameters: output
* TABLUE  CHARACTER*(LNGTAB) information.
*
*-----------------------------------------------------------------------
*
      INTEGER IDKSGT,LNGTAB,TABSGT(*)
      CHARACTER TABLUE*(*),TEXT4*4
*
      IDK=(IDKSGT+3)/4
      IOF=1
      DO 100 I=1,(LNGTAB+3)/4
        WRITE(TEXT4,'(A4)') TABSGT(IDK+I)
        TABLUE(IOF:IOF+3)=TEXT4
        IOF=IOF+4
  100 CONTINUE
      RETURN
      END
