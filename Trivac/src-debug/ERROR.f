*DECK ERROR
      SUBROUTINE ERROR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reaction rate comparison operator.
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
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): read-only reference macrolib type(L_MACROLIB);
*         HENTRY(2): read-only macrolib type(L_MACROLIB);
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Comments:
* The ERROR: calling specifications are:
* ERROR: MACRO1 MACRO2 :: [ HREA hname ] [ NREG nreg ] ;
* where
*   MACRO1 : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     extended \emph{macrolib} used to compute the reference reaction rates.
*   MACRO2 : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     extended \emph{macrolib} used to compute the approximate reaction rates.
*   HREA   : keyword used to set the character name hname.
*   hname  : name of the nuclear reaction used to compute the power map. By 
*     default, reaction H-FACTOR is used.
*   NREG   : keyword used to set the nreg number.
*   nreg   : integer number set to the number of regions used in statistics. By 
*     default, all available regions are used.
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
      CHARACTER TITLE*72,TEXT12*12,HSIGN*12,TEXT4*4,TEXT6*6,HREAC*8
      INTEGER IDATA(NSTATE)
      DOUBLE PRECISION DFLOTT
      TYPE(C_PTR) IPMAC1,IPMAC2
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LE.1) CALL XABORT('ERROR: TWO PARAMETERS EXPECTED.')
      IF((JENTRY(1).NE.2).OR.(IENTRY(1).LT.1).OR.(IENTRY(1).GT.4))
     1 CALL XABORT('ERROR: LINKED LIST OR FILE IN READ-ONLY MODE EXPE'
     2 //'CTED AT FIRST RHS.')
      IF((JENTRY(2).NE.2).OR.(IENTRY(2).LT.1).OR.(IENTRY(2).GT.4))
     1 CALL XABORT('ERROR: LINKED LIST OR FILE IN READ-ONLY MODE EXPE'
     2 //'CTED AT SECOND RHS.')
*----
*  PROCESS FIRST AND SECOND RHS.
*----
      IF(IENTRY(1).GE.3) THEN
         IFTRAK=FILUNIT(KENTRY(1))
         CALL LCMOP(IPMAC1,'COPY1',0,1,0)
         CALL LCMEXP(IPMAC1,0,IFTRAK,IENTRY(1)-2,2)
      ELSE
         IPMAC1=KENTRY(1)
      ENDIF
      CALL LCMGTC(IPMAC1,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB') THEN
         TEXT12=HENTRY(1)
         CALL XABORT('ERROR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MACROLIB EXPECTED.')
      ENDIF
      IF(IENTRY(2).GE.3) THEN
         IFTRAK=FILUNIT(KENTRY(2))
         CALL LCMOP(IPMAC2,'COPY2',0,1,0)
         CALL LCMEXP(IPMAC2,0,IFTRAK,IENTRY(2)-2,2)
      ELSE
         IPMAC2=KENTRY(2)
      ENDIF
      CALL LCMGTC(IPMAC2,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('ERROR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MACROLIB EXPECTED.')
      ENDIF
      CALL LCMLEN(IPMAC2,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPMAC2,'TITLE',72,1,TITLE)
      ELSE
         TITLE='*** NO TITLE PROVIDED ***'
      ENDIF
      WRITE(6,'(/1X,A72)') TITLE
*
      CALL LCMGET(IPMAC1,'STATE-VECTOR',IDATA)
      NGRP=IDATA(1)
      NREG=IDATA(2)
      CALL LCMGET(IPMAC2,'STATE-VECTOR',IDATA)
      IF((NREG.NE.IDATA(2)).OR.(NGRP.NE.IDATA(1))) THEN
         WRITE (6,'(/16H REFERENCE NREG=,I7,6H NGRP=,I7)') NREG,NGRP
         WRITE (6,'(/16H   APPROX. NREG=,I7,6H NGRP=,I7)') IDATA(2),
     1   IDATA(1)
         CALL XABORT('ERROR: REFERENCE AND APPROX. DATA ARE NOT COMPA'
     1   //'TIBLE.')
      ENDIF
*----
*  READ THE MAC: MODULE OPTIONS.
*----
      NREG2=NREG
      IMPX=1
      IPICK=0
      HREAC='H-FACTOR'
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 20
      IF(INDIC.NE.3) CALL XABORT('ERROR: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
*        SET EDITION
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('ERROR: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'NREG') THEN
*        SET NUMBER OF REGIONS
         CALL REDGET(INDIC,NREG2,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('ERROR: INTEGER DATA EXPECTED.')
         IF((NREG2.LE.0).OR.(NREG2.GT.NREG)) THEN
            CALL XABORT('ERROR: INVALID NUMBER OF REGIONS AFTER NREG.')
         ENDIF
      ELSE IF(TEXT4.EQ.'HREA') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,HREAC,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('ERROR: CHARACTER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 20
      ELSE IF(TEXT4.EQ.'PICK') THEN
         IPICK=1
         GO TO 20
      ELSE
         CALL XABORT('ERROR: '//TEXT4//' IS AN INVALID KEY-WORD.')
      ENDIF
      GO TO 10
*----
*  COMPUTE STATISTICS
*----
   20 CALL ERRDRV(IPMAC1,IPMAC2,NREG,NREG2,NGRP,HREAC,ERAMAX,ERASUM,
     1 ERQMAX,ERQSUM)
      IF(IENTRY(1).GE.3) CALL LCMCL(IPMAC1,2)
      IF(IENTRY(2).GE.3) CALL LCMCL(IPMAC2,2)
*----
*  PICK STATISTICS AS CLE200 VARIABLES
*----
      IF(IPICK.EQ.1) THEN
   30    CALL REDGET(INDIC,NITMA,FLOTT,TEXT6,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('ERROR: CHARACTER DATA EXPECTED.')
         IF(TEXT6.EQ.';') RETURN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.-2) CALL XABORT('ERROR: OUTPUT REAL EXPECTED.')
         INDIC=2
         IF(TEXT6.EQ.'ERAMAX') THEN
           FLOTT=ERAMAX
         ELSE IF(TEXT6.EQ.'ERASUM') THEN
           FLOTT=ERASUM
         ELSE IF(TEXT6.EQ.'ERQMAX') THEN
           FLOTT=ERQMAX
         ELSE IF(TEXT6.EQ.'ERQSUM') THEN
           FLOTT=ERQSUM
         ELSE
           CALL XABORT('ERROR: INVALID KEYWORD: '//TEXT6//'.')
         ENDIF
         IF(IMPX.GT.0) WRITE(6,40) TEXT6,FLOTT
         CALL REDPUT(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         GO TO 30
      ENDIF
      RETURN
   40 FORMAT(/13H ERROR: PICK ,A,1H=,1P,E12.4,2H %)
      END
