*DECK GEO
      SUBROUTINE GEO(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Geometry definition operator.
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
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_GEOM).
*         HENTRY(2): optional read-only type(L_GEOM).
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
      CHARACTER    TEXT12*12,TEXT13*12
      TYPE(C_PTR)  IPLIST
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.EQ.0) CALL XABORT('GEO: PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('GEO: LCM O'
     1 //'BJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('GEO: CREAT'
     1 //'E OR MODIFICATION MODE EXPECTED.')
      ITYPE=JENTRY(1)
      IPLIST=KENTRY(1)
*
      IMPX=1
      IF((ITYPE.EQ.0).AND.(NENTRY.GT.1)) THEN
*        CREATE A NEW GEOMETRY BASED ON AN EXISTING ONE.
         IF(JENTRY(2).NE.2) CALL XABORT('GEO: RHS GEOMETRY EXPECTED OPE'
     1   //'N IN READ-ONLY MODE.')
         IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)) CALL XABORT('GEO: LC'
     1   //'M OBJECT EXPECTED AT RHS.')
         CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,TEXT12)
         IF(TEXT12.NE.'L_GEOM') THEN
            TEXT13=HENTRY(2)
            CALL XABORT('GEO: SIGNATURE OF '//TEXT13//' IS '//TEXT12//
     1      '. L_GEOM EXPECTED(1).')
         ENDIF
         CALL LCMEQU(KENTRY(2),IPLIST)
      ELSE IF(ITYPE.EQ.1) THEN
*        MODIFY AN EXISTING GEOMETRY USING THE SAME NAME.
         CALL LCMGTC(IPLIST,'SIGNATURE',12,1,TEXT12)
         IF(TEXT12.NE.'L_GEOM') THEN
            TEXT13=HENTRY(1)
            CALL XABORT('GEO: SIGNATURE OF '//TEXT13//' IS '//TEXT12//
     1      '. L_GEOM EXPECTED(2).')
         ENDIF
      ENDIF
*
      TEXT12='/'
      CALL GEOIN1 (TEXT12,IPLIST,1,IMPX,MAXMIX)
      RETURN
      END
