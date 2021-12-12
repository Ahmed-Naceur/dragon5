*DECK FMTGET
      SUBROUTINE FMTGET(IPRINT,NOPT,IOPT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read formatting command for the FMT module.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input/output
* IPRINT  print level.
* NOPT    number of options.
* IOPT    processing option.
*
*Comments:
* Input data is of the form:
*    [ EDIT iprint ]
*    {
*      SUS3D { SN | CP } |
*      DIRFLX
*      }
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NOPT,IOPT(NOPT)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMTGET')
*----
*  Variables for input via REDGET
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Initialize default values for IPRINT
*----
      IPRINT=1
*----
*  Get data from input file
*----
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 110  CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected.')
      IF(CARLIR .EQ. ';') THEN
        GO TO 105
      ELSE IF(CARLIR .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- print level expected after EDIT.')
        IPRINT=INTLIR
      ELSE IF(CARLIR .EQ. 'SUS3D') THEN
        IF(IOPT(1) .NE. 0) CALL XABORT(NAMSBR//
     >  ': Only one formatting option permitted.')
        IOPT(1)=1
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Character variable for SUS3D option.')
        IOPT(2)=0
        IF(CARLIR .EQ. 'CP') THEN
          IOPT(2)=1
        ELSE IF(CARLIR .EQ. 'SN') THEN
          IOPT(2)=0
        ELSE
          GO TO 110
        ENDIF
      ELSE IF(CARLIR .EQ. 'DIRFLX') THEN
        IF(IOPT(1) .NE. 0) CALL XABORT(NAMSBR//
     >  ': Only one formatting option permitted.')
        IOPT(1)=2
      ELSE IF(CARLIR .EQ. 'BURNUP') THEN
        IF(IOPT(1) .NE. 0) CALL XABORT(NAMSBR//
     >  ': Only one formatting option permitted.')
        IOPT(1)=3
*----
*  Note : input will be completed in FMTBRN
*----
        GO TO 105
      ELSE
        CALL XABORT(NAMSBR//': Keyword '//CARLIR//' is invalid.')
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  Processing finished, return
*----
      RETURN
      END
