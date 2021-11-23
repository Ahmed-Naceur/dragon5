*DECK LIBCOV
      SUBROUTINE LIBCOV(TEXT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Convert a lower-case character variable to upper case.
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
*Parameters: input/output
* TEXT  variable to be converted.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*(*) TEXT
*
      DO 10 I=1,LEN(TEXT)
      IF(TEXT(I:I).EQ.'a') TEXT(I:I)='A'
      IF(TEXT(I:I).EQ.'b') TEXT(I:I)='B'
      IF(TEXT(I:I).EQ.'c') TEXT(I:I)='C'
      IF(TEXT(I:I).EQ.'d') TEXT(I:I)='D'
      IF(TEXT(I:I).EQ.'e') TEXT(I:I)='E'
      IF(TEXT(I:I).EQ.'f') TEXT(I:I)='F'
      IF(TEXT(I:I).EQ.'g') TEXT(I:I)='G'
      IF(TEXT(I:I).EQ.'h') TEXT(I:I)='H'
      IF(TEXT(I:I).EQ.'i') TEXT(I:I)='I'
      IF(TEXT(I:I).EQ.'j') TEXT(I:I)='J'
      IF(TEXT(I:I).EQ.'k') TEXT(I:I)='K'
      IF(TEXT(I:I).EQ.'l') TEXT(I:I)='L'
      IF(TEXT(I:I).EQ.'m') TEXT(I:I)='M'
      IF(TEXT(I:I).EQ.'n') TEXT(I:I)='N'
      IF(TEXT(I:I).EQ.'o') TEXT(I:I)='O'
      IF(TEXT(I:I).EQ.'p') TEXT(I:I)='P'
      IF(TEXT(I:I).EQ.'q') TEXT(I:I)='Q'
      IF(TEXT(I:I).EQ.'r') TEXT(I:I)='R'
      IF(TEXT(I:I).EQ.'s') TEXT(I:I)='S'
      IF(TEXT(I:I).EQ.'t') TEXT(I:I)='T'
      IF(TEXT(I:I).EQ.'u') TEXT(I:I)='U'
      IF(TEXT(I:I).EQ.'v') TEXT(I:I)='V'
      IF(TEXT(I:I).EQ.'w') TEXT(I:I)='W'
      IF(TEXT(I:I).EQ.'x') TEXT(I:I)='X'
      IF(TEXT(I:I).EQ.'y') TEXT(I:I)='Y'
      IF(TEXT(I:I).EQ.'z') TEXT(I:I)='Z'
   10 CONTINUE
      RETURN
      END
