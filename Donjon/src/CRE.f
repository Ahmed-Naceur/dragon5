*DECK CRE
      SUBROUTINE CRE(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover and interpolate a macrolib from one or many compo objects;
* generate a fuel-map macrolib.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input/output
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         IENTRY=1 for LCM memory object;
*         IENTRY=2 for XSM file;
*         IENTRY=3 for sequential binary file;
*         IENTRY=4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         JENTRY=0 for a data structure in creation mode;
*         JENTRY=1 for a data structure in modifications mode;
*         JENTRY=2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
*
*Comments:
* The CRE: module specifications are:
* Option 1: 
* MACRO := CRE: [ MACRO ] [[ CPO ]]  :: (desccre1) ;
* Option 2
* MACFL := CRE: [[ CPO ]] FMAP :: (desccre2) ;
* where
*   MACRO : name of the \emph{macrolib}
*     object to be created or updated for the few reactor material properties.
*     Note that if MACRO appears on the RHS, the information previously
*     stored in MACRO is kept.
*   CPO   : name of the \emph{compo} 
*     object containing the mono-parameter database from transport calculations.
*   MACFL : name of the fuel-map \emph{macrolib}
*     that will be created only for the fuel properties over the fuel lattice.
*   FMAP  : name of the \emph{fmap} 
*     object containing the fuel-map specification and burnup informations.
*   (desccre1) : structure describing the input data to the CRE:
*     module when the \emph{fmap} object is not specified.
*   (desccre2) : structure describing the input data to the CRE:
*     module for the fuel-map \emph{macrolib} construction.
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
      PARAMETER(NSTATE=40)
      TYPE(C_PTR) IPMAC,IPMAP
      CHARACTER TEXT*12,HSMG*131,HSIGN*12
      INTEGER ISTATE(NSTATE)
      LOGICAL LMAC
      DOUBLE PRECISION DFLOT
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LE.1)CALL XABORT('@CRE: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@CRE'
     1 //': LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).EQ.0)THEN
        HSIGN='L_MACROLIB'
        CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        LMAC=.FALSE.
      ELSEIF(JENTRY(1).EQ.1)THEN
        CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_MACROLIB')THEN
          TEXT=HENTRY(1)
          CALL XABORT('@CRE: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_MACROLIB EXPECTED.')
        ENDIF
        LMAC=.TRUE.
      ELSE
        CALL XABORT('@CRE: MACROLIB IN CREATE OR MODIFICATION MOD'
     1   //'E EXPECTED.')
      ENDIF
      IPMAC=KENTRY(1)
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@CRE:'
     1 //' LCM OBJECT EXPECTED AT RHS.')
      IF(JENTRY(2).NE.2)CALL XABORT('@CRE: COMPO IN READ-ONLY MOD'
     1 //'E EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_COMPO')THEN
        TEXT=HENTRY(2)
        CALL XABORT('@CRE: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_COMPO EXPECTED.')
      ENDIF
      IPMAP=C_NULL_PTR
      IF(NENTRY.EQ.2)GOTO 10
      DO 5 IEN=3,NENTRY
      IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))CALL XABORT('@C'
     1 //'RE: LCM OBJECT EXPECTED AT RHS.')
      IF(JENTRY(IEN).NE.2)CALL XABORT('@CRE: LCM OBJECT IN READ-ON'
     1 //'LY MODE EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_COMPO')THEN
        IF(HSIGN.EQ.'L_MAP')THEN
          IF(LMAC)CALL XABORT('@CRE: MACROLIB IN CREATE MODE EXPEC'
     1     //'TED WITH FUEL-MAP OBJECT.')
          IF(IEN.EQ.NENTRY)THEN
            IPMAP=KENTRY(IEN)
          ELSE
            CALL XABORT('@CRE: FUEL-MAP OBJECT EXPECTED TO BE THE '
     1       //'LAST PARAMETER.')
          ENDIF
        ELSE
          TEXT=HENTRY(IEN)
          CALL XABORT('@CRE: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1    '. L_COMPO EXPECTED.')
        ENDIF
      ENDIF
    5 CONTINUE
*----
*  RECOVER INFORMATION
*----
   10 CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(2),'STATE-VECTOR',ISTATE)
      NGRP=ISTATE(2)
      NL=ISTATE(4)
      NMIXT=0
      IF(C_ASSOCIATED(IPMAP)) CALL LCMLEN(IPMAP,'FLMIX',NMIXT,ITYP)
      IF(LMAC)THEN
        CALL XDISET(ISTATE,NSTATE,0)
        CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NGRP)THEN
          WRITE(HSMG,'(40HCRE: INCONSISTENT NB OF GROUPS. IN MACRO,
     1    5HLIB =,I5,11H IN COMPO =,I5)') ISTATE(1),NGRP
          CALL XABORT(HSMG)
        ENDIF
        IF(ISTATE(3).NE.NL)THEN
          WRITE(HSMG,'(40HCRE: INCONSISTENT NB OF LEGENDRE ORDERS.,
     1    14H IN MACROLIB =,I5,11H IN COMPO =,I5)') ISTATE(3),NL
          CALL XABORT(HSMG)
        ENDIF
        NMIXT=ISTATE(2)
      ENDIF
*----
* READ THE INPUT DATA
*----
      IMPX=0
   20 CALL REDGET(INDIC,NITMA,FLOT,TEXT,DFLOT)
      IF(INDIC.NE.3)CALL XABORT('@CRE: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'EDIT')THEN
*       READ THE PRINT INDEX.
        CALL REDGET(INDIC,IMPX,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1)CALL XABORT('@CRE: INTEGER DATA EXPECTED(1).')
      ELSEIF(TEXT.EQ.'NMIX')THEN
*       READ THE MAXIMUM NUMBER OF MATERIAL MIXTURES.
        IF(NMIXT.NE.0)CALL XABORT('@CRE: NMIX IS ALREADY DEFINED.')
        CALL REDGET(INDIC,NMIXT,FLOT,TEXT,DFLOT)
        IF(INDIC.NE.1)CALL XABORT('@CRE: INTEGER DATA EXPECTED(2).')
      ELSEIF(TEXT.EQ.'READ')THEN
        IF(NMIXT.EQ.0)CALL XABORT('@CRE: ZERO NUMBER OF MIXTURES.')
        IF(NGRP.EQ.0)CALL XABORT('@CRE: ZERO NUMBER OF GROUPS.')
        CALL CREDRV(IPMAC,IPMAP,NENTRY,HENTRY,KENTRY,LMAC,NMIXT,NGRP,
     1  NL,ILEAK,IMPX)
        GOTO 30
      ELSE
        CALL XABORT('@CRE: '//TEXT//' IS AN INVALID KEYWORD.')
      ENDIF
      GOTO 20
   30 CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGRP
      ISTATE(2)=NMIXT
      ISTATE(3)=NL
      ISTATE(4)=1
      ISTATE(9)=ILEAK
      CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
      IF(IMPX.GT.1)CALL LCMLIB(IPMAC)
      RETURN
      END
