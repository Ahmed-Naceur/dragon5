*DECK DUO
      SUBROUTINE DUO(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perturbative analysis of two systems and determination of the origins
* of Keff discrepancies.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
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
*         HENTRY(1): read-only type(L_LIBRARY) first system;
*         HENTRY(2): read-only type(L_LIBRARY) second system.
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
      TYPE(C_PTR) IPLIB1,IPLIB2
      CHARACTER HSIGN*12,CARLIR*12
      INTEGER ISTATE(NSTATE)
      REAL REALIR
      DOUBLE PRECISION DBLLIR
      LOGICAL LENER,LISOT,LMIXT,LREAC
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.2) CALL XABORT('DUO: TWO PARAMETERS EXPECTED.')
      DO IEN=1,2
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2)) CALL XABORT('DUO'
     1  //': LCM OBJECT EXPECTED AT LHS.')
        IF(JENTRY(IEN).NE.2) CALL XABORT('DUO: ENTRY IN READ-ONLY MODE'
     1  //' EXPECTED.')
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_LIBRARY') THEN
          CARLIR=HENTRY(IEN)
          CALL XABORT('DUO: SIGNATURE OF '//CARLIR//' IS '//HSIGN//
     1    '. L_LIBRARY EXPECTED.')
        ENDIF
      ENDDO
      IPLIB1=KENTRY(1)
      IPLIB2=KENTRY(2)
      CALL LCMGET(IPLIB1,'STATE-VECTOR',ISTATE)
      NG=ISTATE(3)
      CALL LCMGET(IPLIB2,'STATE-VECTOR',ISTATE)
      IF(ISTATE(3).NE.NG) CALL XABORT('DUO: INVALID NUMBER OF ENERGY G'
     1 //'ROUPS.')
*---
*  READ DATA
*---
      IPRINT=1
      LENER=.FALSE.
      LISOT=.FALSE.
      LMIXT=.FALSE.
      LREAC=.FALSE.
   10 CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.EQ.10) GO TO 100
      IF(ITYPLU.NE.3) CALL XABORT('DUO: READ ERROR - CHARACTER VARI'
     > //'ABLE EXPECTED')
      IF(CARLIR.EQ.';') THEN
        GO TO 100
      ELSE IF(CARLIR.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IPRINT,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('DUO: READ ERROR - INTEGER VARI'
     >  //'ABLE EXPECTED')
      ELSE IF(CARLIR.EQ.'ENERGY') THEN
        LENER=.TRUE.
      ELSE IF(CARLIR.EQ.'ISOTOPE') THEN
        LISOT=.TRUE.
      ELSE IF(CARLIR.EQ.'MIXTURE') THEN
        LMIXT=.TRUE.
      ELSE IF(CARLIR.EQ.'REAC') THEN
        LREAC=.TRUE.
        GO TO 100
      ELSE
        CALL XABORT('DUO: ILLEGAL KEYWORD '//CARLIR)
      ENDIF
      GO TO 10
  100 CALL LCMGET(IPLIB1,'STATE-VECTOR',ISTATE)
      NMIX=ISTATE(1)
      NISOT=ISTATE(2)
      NGRP=ISTATE(3)
      CALL LCMGET(IPLIB2,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NMIX) CALL XABORT('DUO: THE TWO MICROLIBS HAVE A'
     > //' DIFFERENT NUMBER OF MIXTURES.')
      IF(ISTATE(2).NE.NISOT) CALL XABORT('DUO: THE TWO MICROLIBS HAVE '
     > //'A DIFFERENT NUMBER OF ISOTOPES.')
      IF(ISTATE(3).NE.NGRP) CALL XABORT('DUO: THE TWO MICROLIBS HAVE A'
     > //' DIFFERENT NUMBER OF GROUPS.')
*
      CALL DUODRV(IPLIB1,IPLIB2,IPRINT,LENER,LISOT,LMIXT,LREAC,NMIX,
     > NISOT,NGRP)
      RETURN
      END
