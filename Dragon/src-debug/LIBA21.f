*DECK LIBA21
      SUBROUTINE LIBA21(TYPSEG,RETCAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Description of the archive segments used in the APOLIB-2 file
* used as external subroutine by subroutine AEXTRT.
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
*Parameters: input
* TYPSEG  type of segment.
*
*Parameters: output
* RETCAR  encoded generic character type of segment.
*
*-----------------------------------------------------------------------
*
      CHARACTER*(*) TYPSEG
      CHARACTER*(*) RETCAR
*
      IF (TYPSEG .EQ. 'APOLIB') THEN
        RETCAR = 'LIIIIIIIIIII'
*
      ELSEIF (TYPSEG .EQ. 'PCOM') THEN
        RETCAR = 'C1C1'
*
      ELSEIF (TYPSEG .EQ. 'PCONST') THEN
        RETCAR = 'I1I1I1I1I1I1R1'
*
      ELSEIF (TYPSEG .EQ. 'PFIX') THEN
        RETCAR = 'IIILLLIIIRRI1I1I1L1R1R1R1C1C1'
*
      ELSEIF (TYPSEG .EQ. 'PFLUXC') THEN
        RETCAR = 'R1'
*
      ELSEIF (TYPSEG .EQ. 'PHEAD ') THEN
        RETCAR = 'CC1C1R1'
*
      ELSEIF (TYPSEG .EQ. 'PMAIL ') THEN
        RETCAR = 'RIR1R1I1I1'
*
      ELSEIF (TYPSEG .EQ. 'PNUMF') THEN
        RETCAR = 'C1C1R3'
*
      ELSEIF (TYPSEG .EQ. 'PPPSN') THEN
        RETCAR = 'IIIIIII1I1I1'
*
      ELSEIF (TYPSEG .EQ. 'PPSN') THEN
        RETCAR = 'R1'
*
      ELSEIF (TYPSEG .EQ. 'PSECT') THEN
        RETCAR = 'R1'
*
      ELSEIF (TYPSEG .EQ. 'PTHOM1') THEN
        RETCAR = 'RCI1I1I1I2I2R1R1R1R1'
*
      ELSEIF (TYPSEG .EQ. 'PTHOM2') THEN
        RETCAR = 'R3R3R3'
*
      ELSEIF (TYPSEG .EQ. 'PTHOM3') THEN
        RETCAR = 'R1R1R1R1'
*
      ELSEIF (TYPSEG .EQ. 'PTHOM4') THEN
        RETCAR = 'R1R1R1R1R1R1R1R1R1R1R1'
*
      ELSEIF (TYPSEG .EQ. 'PTHOM5') THEN
        RETCAR = 'R1R1R1'
*
      ELSEIF (TYPSEG .EQ. 'QFIX') THEN
        RETCAR = 'I1'
*
      ELSEIF (TYPSEG .EQ. 'QFIXS') THEN
        RETCAR = 'I1'
*
      ELSEIF (TYPSEG .EQ. 'QFLUXC') THEN
        RETCAR = 'I1'
      ELSE
        RETCAR = ' '
      ENDIF
*
      RETURN
      END
