*DECK KDRVER
      SUBROUTINE KDRVER(REV,DATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract CVS or SVN version and production date.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: output
* REV     revision character identification
* DATE    revision character date 
*
*-----------------------------------------------------------------------
*
      CHARACTER REV*48,DATE*64
*
      REV='Version 5.0.7 ($Revision: 2060 $)'
      DATE='$Date: 2021-02-02 16:03:29 -0500 (Tue, 02 Feb 2021) $'
      IF(REV(22:).EQ.'ion$)') THEN
*        CVS or SVN keyword expansion not performed
         REV='Version 5.0.7'
         DATE='February 2, 2021'
      ENDIF
      RETURN
      END
