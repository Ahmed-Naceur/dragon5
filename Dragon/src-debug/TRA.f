*DECK TRA
      SUBROUTINE TRA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transposition of a macrolib.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) creation type(L_MACROLIB);
*         HENTRY(2) read-only type(L_MACROLIB).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
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
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) IPMAC1,IPMAC2
      CHARACTER HSIGN*12,TEXT12*12
      INTEGER ISTATE(NSTATE)
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.2) CALL XABORT('T: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('T: LI'
     1 //'NKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('T: ENTRY IN CREATE OR MODE E'
     1 //'XPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('T: LINKED LIST OR XSM FILE IN READ-ONLY MODE E'
     2 //'XPECTED AT RHS.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IPMAC1=KENTRY(1)
      IF(HSIGN.EQ.'L_MACROLIB') THEN
         IPMAC2=KENTRY(2)
      ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
         IPMAC2=LCMGID(KENTRY(2),'MACROLIB')
      ELSE
         TEXT12=HENTRY(2)
         CALL XABORT('T: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MACROLIB OR L_LIBRARY EXPECTED.')
      ENDIF
*----
*  TRANSPOSITION
*----
      CALL LCMGET(IPMAC2,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      NMIL=ISTATE(2)
      NL=ISTATE(3)
      NF=ISTATE(4)
      NDEL=ISTATE(7)
      ISTEP=ISTATE(11)
      CALL TRAXS(IPMAC1,IPMAC2,NG,NMIL,NL,NF,NDEL,ISTEP)
*----
*  SAVE THE SIGNATURE AND STATE VECTOR
*----
      HSIGN='L_MACROLIB'
      CALL LCMPTC(IPMAC1,'SIGNATURE',12,1,HSIGN)
      IF(ISTATE(13).EQ.0) THEN
        ISTATE(13)=1
      ELSE IF(ISTATE(13).EQ.1) THEN
        ISTATE(13)=0
      ENDIF
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
*
      RETURN
      END
