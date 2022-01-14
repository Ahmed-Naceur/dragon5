*DECK EPCRPD
      SUBROUTINE EPCRPD(NENTRY,KENTRY,IPRINT,NOPT,IOPT,CARRET)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract library parameters with normal distribution
* around the average.
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
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* KENTRY  data structure pointer.
* IPRINT  print level.
* NOPT    number of options.
* IOPT    processing option with:
*         IOPT(1)  type of processing (=0 for current option);
*         IOPT(2)  entry position for L_EPC structure;
*         IOPT(3)  number of parameters;
*         IOPT(4)  entry for normal distribution file;
*         IOPT(5)  number of records on normal distribution file.
* CARRET  last input option read.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NENTRY
      TYPE(C_PTR)      KENTRY(NENTRY)
      INTEGER          IPRINT,NOPT,IOPT(NOPT)
      CHARACTER*12     CARRET
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRPD')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Input and output parameters
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Local variables
*----
      TYPE(C_PTR)      IPNDI
      INTEGER          IKNDI
      INTEGER          ISTATE(NSTATE)
      INTEGER          NOREC,NTREC,NFREC,MAXD
*----
*  Output structure
*----
      IKNDI=1
      IF(IOPT(3) .EQ. 0) IKNDI=-1
      IPNDI=KENTRY(IOPT(2))
      CALL LCMGET(IPNDI,'STATE-VECTOR',ISTATE)
*----
*  Input structure
*----
      NOREC=ISTATE(1)
      MAXD=ISTATE(2)
      NTREC=IOPT(5)
*----
*  Recover parameter names from output and input structures
*----
      NFREC=NTREC+NOREC
      ISTATE(1)=NFREC
      ISTATE(2)=MAXD
      CALL LCMPUT(IPNDI,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  INPUT/OUTPUT VARIABLES
*  Input data is of the form
*    [ EDIT iprint ]
*    [ SET (dataSET) ]
*    [ GET (dataGET) ]
*  where
*  EDIT               = keyword for print level
*  iprint             = integer print level
*  SET                = keyword to set reference value
*                       and extract number of parameters
*  GET                = keyword to set next value
*                       and to extract parameter
*  (dataSET)          = SET data processed by NDISET
*  (dataGET)          = GET data processed by NDIGET
*----
      IPRINT=1
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      IF(CARLIR .EQ. ';') THEN
        GO TO 105
      ELSE IF(CARLIR .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer value for EDIT expected.')
        IPRINT=INTLIR
      ELSE IF(CARLIR .EQ. 'SET') THEN
        CARRET=CARLIR
        CARLIR=CARRET
        GO TO 101
      ELSE IF(CARLIR .EQ. 'GET') THEN
        CARRET=CARLIR
        CARLIR=CARRET
        GO TO 101
      ELSE
        CALL XABORT(NAMSBR//
     >  ': '//CARLIR//' is an invalid keyword.')
      ENDIF
      GO TO 100
 105  CONTINUE
      RETURN
      END
