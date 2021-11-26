*DECK LIBXS6
      SUBROUTINE LIBXS6 (MAXDIL,NAMFIL,HSHI,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in Apolib-XSM format.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXDIL  maximum number of dilutions.
* NAMFIL  name of the APOLIB-XSM file.
* HSHI    library name of the self-shielding data.
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
      INTEGER MAXDIL,NDIL
      CHARACTER HSHI*12
      CHARACTER NAMFIL*(*)
      REAL DILUT(MAXDIL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPAP
      CHARACTER TEXT20*20,TEXT12*12,HSMG*131
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NOMS
*----
*  OPEN AND PROBE THE APOLIB-XSM FILE.
*----
      CALL LCMOP(IPAP,NAMFIL,2,2,0)
      CALL LCMSIX(IPAP,'PHEAD',1)
      CALL LCMLEN(IPAP,'NOMS',NV,ITYLCM)
      NISOTS=NV/5
      ALLOCATE(NOMS(5*NISOTS))
      CALL LCMGET(IPAP,'NOMS',NOMS)
      KISEG=0
      DO ISO=1,NISOTS
        WRITE(TEXT20,'(5A4)') (NOMS((ISO-1)*5+II),II=1,5)
        IF(TEXT20(:12).EQ.HSHI) THEN
          KISEG=ISO
          EXIT
        ENDIF
      ENDDO
      DEALLOCATE(NOMS)
      IF(KISEG.EQ.0) THEN
        WRITE(HSMG,'(45HLIBXS6: UNABLE TO FIND SELF-SHIELDED ISOTOPE ,
     1  A12,1H.)') HSHI
        CALL XABORT(HSMG)
      ENDIF
      CALL LCMSIX(IPAP,' ',2)
*----
*  RECOVER DILUTIONS
*----
      CALL LCMSIX(IPAP,'QFIXS',1)
      WRITE(TEXT12,'(4HISOT,I8.8)') KISEG
      CALL LCMSIX(IPAP,TEXT12,1)
      CALL LCMSIX(IPAP,'SSDATA',1)
      CALL LCMLEN(IPAP,'SEQHOM',NDIL,ITYLCM)
      IF(NDIL.EQ.0) THEN
        WRITE(HSMG,'(47HLIBXS6: NO DILUTIONS FOR SELF-SHIELDED ISOTOPE ,
     1  A12,1H.)') HSHI
        CALL XABORT(HSMG)
      ELSE IF(NDIL.GT.MAXDIL) THEN
        WRITE(HSMG,'(46HLIBXS6: MAXDIL OVERFLOW SELF-SHIELDED ISOTOPE ,
     1  A12,1H.)') HSHI
        CALL XABORT(HSMG)
      ENDIF
      NDIL=NDIL-1
      CALL LCMGET(IPAP,'SEQHOM',DILUT)
      CALL LCMSIX(IPAP,' ',2)
      CALL LCMSIX(IPAP,' ',2)
      CALL LCMSIX(IPAP,' ',2)
      CALL LCMCL(IPAP,1)
      RETURN
      END
