*DECK XDSSET
      SUBROUTINE XDSSET(VECTOR,NWORDS,VALUE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To initialize the \verb|NWORDS| first elements of the character
* string array \verb|VECTOR| to the character string variable 
* \verb|VALUE|.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert.
*
*Parameters: input
* NWORDS  number of words to process.
* VALUE   variable to store in the array.
*
*Parameters: input/output
* VECTOR  array to initialize.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NWORDS
      CHARACTER        VECTOR(NWORDS)*(*),VALUE*(*)
*----
*  Local variables
*----
      INTEGER          IWORD,ICHAR,ILONVE,ILONVA
      ILONVE=LEN(VECTOR(1))
      ILONVA=LEN(VALUE) 
      DO 100 IWORD=1,NWORDS
        DO 101 ICHAR=1,MIN(ILONVA,ILONVE)
          VECTOR(IWORD)(ICHAR:ICHAR)=VALUE(ICHAR:ICHAR)
 101    CONTINUE
        DO 102 ICHAR=ILONVA+1,ILONVE
          VECTOR(IWORD)(ICHAR:ICHAR)=' '
 102    CONTINUE
 100  CONTINUE
      RETURN
      END
