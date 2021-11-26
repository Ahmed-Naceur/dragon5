*DECK AFM
      SUBROUTINE AFM(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Generate a macrolib using the AFM feedback model
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):  
* M.T. Sissaoui
*
*Update(s):
*  E. Varin, B. Dionne
*
*Parameters: input
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
*Reference:
* M. T. Sissaoui, G. Marleau and D. Rozon, "CANDU Reactor Simulations
* Using the Feedback Model with Actinide Burnup History," Nucl.
* Technology, 125, 197 (1999).
*
*Comments:
* The AFM: calling specifications are:
* MACRO := AFM: [ MACRO ] DBASE [ MAPFL ] :: (descafm) ;
* where
*   MACRO : name of the extended \emph{macrolib}
*   DBASE : name of the \emph{database} object containing fuel properties with 
*     respect to local parameters.
*   MAPFL : name of the \emph{map} object containing fuel regions description 
*     and burnupinformations. This file is only required when a \emph{MACRO is 
*     created for fuel area.
*   (descafm) : structure containing the data to module AFM:.
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
      CHARACTER TEXT12*12,HSMG*131,HSIGN*12,TEXT*5,CTITRE*72,TINFO*72
      LOGICAL LMCR,LMAP
      DOUBLE PRECISION DFLOTT
      INTEGER  IPAR(NSTATE),IDATA(NSTATE)
      TYPE(C_PTR) IPLIST
*
      LMCR=.FALSE.
      LMAP=.FALSE.
      MSFT=0
*
* PARAMETER VALIDATION.
      IF(NENTRY.LE.1) CALL XABORT('AFM: 2 PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('AFM:'
     1 //'  MACROLIB LINKED LIST OR XSM FILE EXPECTED AT LHS.')
*
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)) CALL XABORT('AFM:'
     1 //' DATABASE LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(2).NE.2) CALL XABORT('AFM: DATABASE IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
*
      CALL REDGET (INDIC,NITMA,FLOTT,TEXT,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('AFM: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'MCR') THEN
        LMCR=.TRUE.
        WRITE(6,'(A37)') 'AFM: GENERATION OF A SINGLE MACROLIB'
        CALL REDGET (INDIC,MXSH,FLOTT,TEXT,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('AFM: INTEGER DATA EXPECTED.')
      ELSEIF(TEXT.EQ.'MAP') THEN
        LMAP=.TRUE.
        IF((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2)) CALL XABORT('AFM:'
     1    //' FUEL MAP LINKED LIST OR XSM FILE EXPECTED AT RHS.')
        IF(JENTRY(3).NE.2) CALL XABORT('AFM: COMPO IN READ-ONLY '
     1    //'MODE EXPECTED AT RHS.')
      ELSE
         CALL XABORT('AFM: MAP OR MCR KEY WORD EXPECTED')
      ENDIF
*
      ITYPE=JENTRY(1)
      IPLIST=KENTRY(1)
      MMIX=1
*---------------------------------------------------------------*
* CHECK THE SIGNTURE OF THE LINKED LIST OR XSM FILE
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'REACTOR_XSDB') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('AFM: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. REACTOR_XSDB EXPECTED.')
      ENDIF
*---------------------------------------------------------------*
* IF L_MAP IS NOT AVAILABLE AFM GENERATE ONLY A TABLE
      IF(LMAP) THEN
        CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_MAP') THEN
          TEXT12=HENTRY(3)
          CALL XABORT('AFM: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MAP EXPECTED.')
        ENDIF
        WRITE(6,'(A42)') 'AFM: GENARATION OF A MACROLIB USING L_MAP'
      ENDIF
*---------------------------------------------------------------*
*        READ THE INFORMATION TITLE.
      CALL REDGET (INDIC,NITMA,FLOTT,TEXT,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('AFM: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'INFOR') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TINFO,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('AFM: CHARACTER DATA EXPECTED.')
      ELSE
          CALL XABORT('AFM:KEY WORD INFOR EXPECTED.')
      ENDIF
*
* RECOVER SOME INFORMATIONS FROM THE DATABASE.
      TEXT12='INFORMATION'
      CALL LCMLEN(KENTRY(2),TEXT12,LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
        CALL LCMGTC(KENTRY(2),'INFORMATION',72,1,CTITRE)
        WRITE(6,*)'INFORMATION TITLE ',TINFO
*
        IF(CTITRE.NE.TINFO) THEN
            CALL XABORT('AFM: INCONSISTENT TITLES  '//CTITRE//
     1      '  EXPECTED. INSTEAD OF ' //TINFO// ' ')
        ENDIF
      ELSE
        CALL XABORT('AFM: DATA BASE TITLE IS NOT PROVIDED ')
      ENDIF
*---------------------------------------------------------------*
*        CHECK THE NAMES OF THE DIFFERENTS DIRECTORIES
      CALL REDGET (INDIC,NITMA,FLOTT,TEXT,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('AFM: CHARACTER DATA EXPECTED.')
      IF(TEXT.NE.'DNAME') CALL XABORT('AFM:KEY WORD DNAME EXPECTED.')
      CALL REDGET(INDIC,NUT,FLOTT,TEXT,DFLOTT)
      IF(INDIC.NE.1) CALL XABORT('AFM: INTEGER DATA EXPECTED.')
      IF(NUT.GT.1.AND.LMCR) CALL XABORT('AFM: INVALID NUMBER.')
      DO 100 IJ=1,NUT
        CALL REDGET (INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.3) CALL XABORT('AFM: CHARACTER DATA EXPECTED.')
        CALL LCMLEN(KENTRY(2),TEXT12,LENGT,ITYLCM)
        IF(LENGT.EQ.0) THEN
           CALL XABORT('AFM: DATA NAME  '//TEXT12//
     1        '  DO NOT EXIST ')
         ENDIF
 100  CONTINUE
*---------------------------------------------------------------*
*     RECOVER SOME INFORMATIONS FROM DATABASE.
*     TEXT12='SIGNATURE'
      CALL LCMSIX(KENTRY(2),TEXT12,1)
*     RECOVER THE TITLE.
      CALL LCMLEN(KENTRY(2),'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(KENTRY(2),'TITLE',72,1,CTITRE)
      ELSE
         CTITRE='*** NO TITLE PROVIDED ***'
      ENDIF
*     READ PARAMETERS DANS L_FBM 
      CALL LCMGET(KENTRY(2),'PARAM',IPAR)
      CALL LCMSIX(KENTRY(2),' ',2)
      NGRP =IPAR(1)
      NISO =IPAR(2)
      NL   =IPAR(3)
      NBURN=IPAR(4)
      IXYZ =IPAR(5)
*---------------------------------------------------------------*
      IF(LMCR) THEN
        NBCH=MXSH
        NCCO=1
        NCZO=1
        ISC=4
        MMIX=NBCH*NCCO
        MSFT =0
      ELSEIF (LMAP) THEN
* RECOVER INFORMATIONS FROM L_MAP.
*     READ PARAMETERS
        CALL LCMGET(KENTRY(3),'STATE-VECTOR',IPAR)
        NBCH =IPAR(1)
        NCCO =IPAR(2)
        NCZO =IPAR(3)
        ISC  =IPAR(5)
        MSFT =IPAR(6)
        NPARM =IPAR(8)
        MMIX=NBCH*NCCO
C HISTORY PARAMETER
        IF(IPAR(4).NE.NGRP) THEN
            WRITE(HSMG,'(A40,I5,A18,I5)') 'AFM: INCONSISTENT NB OF '
     1        //'GROUPS. IN MAP =',IPAR(4),' IN REACTOR_XSDB =',NGRP
            CALL XABORT(HSMG)
        ENDIF
      ENDIF
*  MSFT IS THE TOTAL NUMBER OF SHIFT
      MNPS=MSFT+2
*  READ THE INPUT DATA.
*  TO USE THE SAME VECTOR TO GET THE REFERENCE LOCAL PARAMETER
*      IF(NISO.LT.8) THEN
*        NISM=8
      IF(NISO.LT.7) THEN
        NISM=7
      ELSE
        NISM=NISO
      ENDIF
      IF(ITYPE.NE.0) THEN
         CALL LCMGET(IPLIST,'STATE-VECTOR',IDATA)
         IF(NGRP.NE.IDATA(1)) CALL XABORT('WRONG NUMBER OF ENERGY'
     1   //' GROUPS IN UPDATED MACROLIB')
         IF(MMIX.NE.IDATA(2)) CALL XABORT('WRONG NUMBER OF MATER'
     1   //'IAL MIXTURES IN UPDATED MACROLIB')
         IF(NL.NE.IDATA(3)) CALL XABORT('WRONG ORDER OF ANISOTROPY'
     1   //'IN UPDATED MACROLIB')
      ENDIF
*---------------------------------------------------------------*
*        NTYP TYPE OF CROSS-SECTIONS CONSIDERED
         NTYP=5+NL+IXYZ*2
*---------------------------------------------------------------*
*        DRIVER TO COMPUTE THE FEEDBACK COEFFICIENTS.
*---------------------------------------------------------------*
      CALL AFMDRV(KENTRY,NENTRY,NPARM,ITYPE,NBURN,NGRP,NISO,ISC,MNPS,
     1 NL,ILEAK,NTYP,NBCH,NCCO,NCZO,NUT,CTITRE,LMCR,IXYZ,MMIX,MSFT,
     2 NISM)
*---------------------------------------------------------------*
      IF(JENTRY(1).EQ.0) THEN
        CALL XDISET(IDATA,NSTATE,0)
        HSIGN='L_MACROLIB'
        CALL LCMPTC(IPLIST,'SIGNATURE',12,1,HSIGN)
        IDATA(1)=NGRP
        IDATA(2)=MMIX
        IDATA(3)=NL
        IDATA(4)=1
        IDATA(9)=ILEAK
        CALL LCMPUT(IPLIST,'STATE-VECTOR',NSTATE,1,IDATA)
      ELSE
        IDATA(1)=NGRP
        IDATA(2)=MMIX
        IDATA(3)=1
        IDATA(4)=1
        IDATA(9)=ILEAK
        CALL LCMPUT(IPLIST,'STATE-VECTOR',NSTATE,1,IDATA)
      ENDIF
      RETURN
      END
