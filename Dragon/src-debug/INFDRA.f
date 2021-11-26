*DECK INFDRA
      SUBROUTINE INFDRA(CFILNA,IPRINT,NBISO,HNAMIS,AWRISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover mass for isotopes of DRAGON libraries.
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
*Parameters: input
* CFILNA  DRAGLIB file name.
* IPRINT  print flag.
* NBISO   number of isotopes.
* HNAMIS  isotope names.
*
*Parameters: output
* AWRISO  isotope weights.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT     NONE
      INTEGER      IPRINT,NBISO
      CHARACTER    CFILNA*64,HNAMIS(NBISO)*8
      REAL         AWRISO(NBISO)
C----
C FUNCTIONS
C----
      DOUBLE PRECISION XDRCST
C----
C  DRAGON LIBRARY PARAMETERS
C----
      TYPE(C_PTR)  IPDRL
      INTEGER      IOUT,ISO,LENGT,ITYLCM
      PARAMETER   (IOUT=6)
      CHARACTER    NAMLOC*12,HSMG*131
      REAL         CONVM
*----
*  For INFDRA, file name is limited to 12 characters
*  because of the requirements for compatibility with 
*  LINKED_LIST
*----
      NAMLOC=CFILNA(1:12)
C----
C  TEST IF FILE NAME EXISTS
C----
      CONVM=REAL(XDRCST('Neutron mass','amu'))
      IF(NAMLOC.EQ.' ' )THEN
        CALL XABORT('INFDRA: DRAGON LIBRARY HAS NOT BEEN SET')
      ENDIF
C----
C  OPEN FILE AND READ INFORMATION DATA RECORDS
C----
      CALL LCMOP(IPDRL,NAMLOC,2,2,0)
      DO 100 ISO=1,NBISO
        CALL LCMLEN(IPDRL,HNAMIS(ISO),LENGT,ITYLCM)
        IF(LENGT.EQ.0) THEN
          CALL LCMLIB(IPDRL)
          WRITE(HSMG,9000) HNAMIS(ISO),CFILNA
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMSIX(IPDRL,HNAMIS(ISO),1)
        CALL LCMGET(IPDRL,'AWR',AWRISO(ISO))
        AWRISO(ISO)=AWRISO(ISO)*CONVM
        IF(IPRINT.GE.100) THEN
          WRITE(IOUT,6000) HNAMIS(ISO),AWRISO(ISO)
        ENDIF
        CALL LCMSIX(IPDRL,' ',2)
 100  CONTINUE
C----
C  CLOSE FILE
C----
      CALL LCMCL(IPDRL,1)
C----
C  RETURN
C----
      RETURN
C----
C  PRINT FORMAT
C----
 6000 FORMAT(' DRAGON ISOTOPE =',A8,
     >       ' HAS ATOMIC WEIGHT RATIO = ',F12.5)
C----
C  ABORT FORMAT
C----
 9000 FORMAT('INFDRA: MATERIAL/ISOTOPE ',A8,
     >       ' IS MISING ON DRAGON LIBRARY FILE ',A64)
      END
