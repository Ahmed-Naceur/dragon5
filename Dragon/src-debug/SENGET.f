*DECK SENGET
      SUBROUTINE SENGET(IPRINT,NL,NANIS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read from the input instructions for the SENS: module.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* G. Marleau
*
*Parameters: output
* IPRINT  print level.
* NL      Legendre contribution in scattering on library.
* NANIS   Legendre contribution in scattering for SENS:
*         = 1  only isotropic contribution; 
*         = ilana contribution up to ilana=NL considered
*         (default is ilana=1).
*
*Comments:
* Input data is of the form:
*    [ EDIT iprint ]
*    [ ANIS ilana ]
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NL,NANIS
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='SENGET')
*----
*  Variables for input via REDGET
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*72
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
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      IF(CARLIR(1:4) .EQ. ';') THEN
        GO TO 105 
      ELSE IF(CARLIR(1:4) .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- print level expected after EDIT.')
        IPRINT=INTLIR
      ELSE IF(CARLIR(1:5) .EQ. 'ANIS') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- scattering order expected')
        IF(INTLIR .LE. 1) THEN
          NANIS=1 
        ELSE 
          NANIS=MIN(NL,INTLIR)
        ENDIF
      ELSE
        CALL XABORT(NAMSBR//': Keyword '//CARLIR(1:5)//' is invalid.')
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  Processing finished, return
*----
      RETURN 
      END 
