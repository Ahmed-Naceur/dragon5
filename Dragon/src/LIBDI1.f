*DECK LIBDI1
      SUBROUTINE LIBDI1 (MAXDIL,IPDRL,HNISOR,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in Draglib format.
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
* MAXDIL  maximum number of dilutions.
* IPDRL   pointer to the Draglib (L_DRAGLIB signature).
* HNISOR  library name of the resonant isotope.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDRL
      INTEGER MAXDIL,NDIL
      CHARACTER HNISOR*12
      REAL DILUT(MAXDIL)
*
      CALL LCMLEN (IPDRL,HNISOR,ILEN,ITYLCM)
      IF(ILEN.EQ.0) CALL XABORT('LIBDI1: ISOTOPE '//HNISOR//' NOT AVAI'
     1 //'LABLE IN THE DRAGLIB.')
      CALL LCMSIX (IPDRL,HNISOR,1)
      CALL LCMLEN (IPDRL,'TEMPERATURE',NTMP,ITYLCM)
      IF(NTMP.GT.0) CALL LCMSIX (IPDRL,'SUBTMP0001',1)
      CALL LCMLEN (IPDRL,'DILUTION',NDIL,ITYLCM)
      IF(NDIL+1.GT.MAXDIL) CALL XABORT('LIBDI1: MAXDIL IS TOO SMALL.')
      IF(NDIL.GT.0) CALL LCMGET (IPDRL,'DILUTION',DILUT)
      DILUT(NDIL+1)=1.0E10
      IF(NTMP.GT.0) CALL LCMSIX (IPDRL,' ',2)
      CALL LCMSIX (IPDRL,' ',2)
      RETURN
      END
