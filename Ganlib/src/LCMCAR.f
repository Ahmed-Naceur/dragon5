*DECK LCMCAR
      SUBROUTINE LCMCAR(TEXT,LACTIO,NITMA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transform a character variable into integer vector back and forth.
* This routine is portable and based on the *ascii* collating sequence,
* equivalence between: text='    ' <=> nitma=0, is imposed.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* TEXT    character variable.
* LACT    logical conversion flag: .true.   character to integer;
*         .false.  integer   to character.
* NITMA   integer (32 bits) vector.
*
*Limitations:
*           it is assumed that:  0 <= ichar() <= 255,
*           otherwise a character would not stand in one byte.
*
*Internal parameters:
* ALPHAB  limited alphabet used for variable names (character*96).
* TASCII  table to convert ichar() values into *ascii* codes.
* IASCII  inversion of tascii().
* IBASE1  integral basis defined as maximum value of ichar()+1;
*         to optimize calculations, it is a power of 2 (128 or 256).
*
*-----------------------------------------------------------------------
*
      IMPLICIT    NONE
      CHARACTER   TEXT*(*)
      LOGICAL     LACTIO
      CHARACTER   ALPHAB*96
      INTEGER     NITMA(*)
      INTEGER     IBASE1,IBASE2,TASCII(0:255),IASCII(0:127)
      INTEGER     I0,I1,I2,I3,J01,J23,K,LMAX,L1,NBDIM,IBDIM
      INTEGER     IWRITE
      PARAMETER ( IWRITE= 6 )
      SAVE IBASE1,IBASE2,TASCII,IASCII
      DATA IBASE1/0/
*
      IF(IBASE1.EQ.0) THEN
*        PREPARE TABLES TASCII() AND IASCII() AND SET INTEGERS IBASE1
*        + IBASE2 FOR CHARACTER/INTEGER CONVERSIONS.
*            0         1         2         3
*            0123456789012345678901234567890123456789
         ALPHAB=' !..$%&.()*+,-./0123456789:;<=>?.ABCDEF'//
     >      'GHIJKLMNOPQRSTUVWXYZ...._.abcdefghijklmn'//
     >      'opqrstuvwxyz.....'
*
         LMAX= 0
         DO 30 K=0,95
            L1= ICHAR(ALPHAB(K+1:K+1))
            LMAX= MAX(LMAX,L1)
            TASCII(L1)= K
            IASCII(K)= L1
   30    CONTINUE
         IF( LMAX.LT.128 )THEN
            IBASE1= 128
         ELSE
            IBASE1= 256
         ENDIF
         IBASE2= IBASE1*IBASE1
      ENDIF
*
      NBDIM= LEN(TEXT)
      IF( MOD(NBDIM,4).NE.0 )THEN
         WRITE(IWRITE,*) 'LCMCAR: LEN(TEXT)=',NBDIM,' NOT / BY 4'
         CALL XABORT('LCMCAR: INVALID CHARACTER <-> INTEGER CONVERSION')
      ELSE
         NBDIM= NBDIM/4
      ENDIF
      IF( LACTIO )THEN
*
*        CONVERT EACH CHARACTER*4 TO INTEGER
         DO 10 IBDIM= 1, NBDIM
            I0= TASCII(ICHAR(TEXT(IBDIM*4-3:IBDIM*4-3)))
            I1= TASCII(ICHAR(TEXT(IBDIM*4-2:IBDIM*4-2)))
            I2= TASCII(ICHAR(TEXT(IBDIM*4-1:IBDIM*4-1)))
            I3= TASCII(ICHAR(TEXT(IBDIM*4  :IBDIM*4  )))
            NITMA(IBDIM)= (I0+IBASE1*I1) + IBASE2*(I2+IBASE1*I3)
   10    CONTINUE
      ELSE
*
*        CONVERT INTEGER TO CHARACTER*4
         DO 20 IBDIM= 1, NBDIM
            J23=   NITMA(IBDIM)/IBASE2
            I3 =     J23/IBASE1
            I2 =     J23-IBASE1*I3
            J01=   NITMA(IBDIM)-J23*IBASE2
            I1 =     J01/IBASE1
            I0 =     J01-IBASE1*I1
            TEXT(IBDIM*4-3:IBDIM*4-3)= CHAR(IASCII(I0))
            TEXT(IBDIM*4-2:IBDIM*4-2)= CHAR(IASCII(I1))
            TEXT(IBDIM*4-1:IBDIM*4-1)= CHAR(IASCII(I2))
            TEXT(IBDIM*4  :IBDIM*4  )= CHAR(IASCII(I3))
   20    CONTINUE
      ENDIF
      RETURN
      END
