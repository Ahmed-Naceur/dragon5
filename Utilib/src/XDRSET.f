*DECK XDRSET
      SUBROUTINE XDRSET(VECTOR,NWORDS,VALUE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To initialize the \verb|NWORDS| first elements of the real
* array \verb|VECTOR| to the real variable \verb|VALUE|.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
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
      REAL             VECTOR(NWORDS),VALUE
*----
*  Local variables
*----
      INTEGER          IWORD
      DO 100 IWORD=1,NWORDS
        VECTOR(IWORD)=VALUE
 100  CONTINUE
      RETURN
      END
