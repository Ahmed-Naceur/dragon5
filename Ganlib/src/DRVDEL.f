*DECK DRVDEL
      SUBROUTINE DRVDEL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* multigroup assembly operator for system matrices.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =-2 the LCM object or file have been destroyed;
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,TEXT4*4
      DOUBLE PRECISION DFLOTT
      LOGICAL LDEL
*
      IMPX=1
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 20
      IF(INDIC.NE.3) CALL XABORT('DRVDEL: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DRVDEL: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 20
      ELSE
         WRITE(HSMG,'(8HDRVDEL: ,A4,22H IS AN INVALID ACTION.)') TEXT4
         CALL XABORT(HSMG)
      ENDIF
      GO TO 10
*
   20 LDEL=.TRUE.
      DO 30 I=1,NENTRY
      IF(JENTRY(I).EQ.1) THEN
         LDEL=.FALSE.
         JENTRY(I)=-2
         IF((IMPX.GT.0).AND.(IENTRY(I).LE.2)) THEN
            WRITE(6,40) HENTRY(I)
         ELSE IF(IMPX.GT.0) THEN
            WRITE(6,50) FILUNIT(KENTRY(I)),HENTRY(I)
         ENDIF
      ENDIF
   30 CONTINUE
      IF(LDEL) CALL XABORT('DRVDEL: NO ENTRY DELETED.')
      RETURN
*
   40 FORMAT(/33H DRVDEL: DELETE LCM OBJECT NAMED ,A)
   50 FORMAT(/25H DRVDEL: DELETE FILE UNIT,I5,7H NAMED ,A)
      END
