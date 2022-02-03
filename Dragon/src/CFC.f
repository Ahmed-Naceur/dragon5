*DECK CFC
      SUBROUTINE CFC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Module to compute the feedback coefficients (CFC) for both the cell
* and the reflector. A CANDU-6 database is generated using the feedback
* model parametrization of the cross-sections. Recover and interpolate
* Macrolib information from one or many Compo objects.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): M. T. Sissaoui
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) : Create or modification type(REACTOR_XSDB).
*         HENTRY(2) : Read-only type(L_COMPO) for reference(nominal).
*         HENTRY(3) : Read-only type(L_COMPO) for fuel-up temp.
*         HENTRY(4) : Read-only type(L_COMPO) for fuel-down temp.
*         HENTRY(5) : Read-only type(L_COMPO) for cool-up temp.
*         HENTRY(6) : Read-only type(L_COMPO) for cool-down temp.
*         HENTRY(7) : Read-only type(L_COMPO) for mod-up temp
*         HENTRY(8) : Read-only type(L_COMPO) for mod-down temp.
*         HENTRY(9) : Read-only type(L_COMPO) for cool-up den.
*         HENTRY(10) : Read-only type(L_COMPO) for cool-down den.
*         HENTRY(11) : Read-only type(L_COMPO) for mod-up den.
*         HENTRY(12) : Read-only type(L_COMPO) for mod-down den.
*         HENTRY(13) : Read-only type(L_COMPO) for boron.
*         HENTRY(14) : Read-only type(L_COMPO) for purity.
*         HENTRY(15) : Read-only type(L_COMPO) for XE135.
*         HENTRY(16) : Read-only type(L_COMPO) for SM149.
*         HENTRY(17) : Read-only type(L_COMPO) for NP239.
*         HENTRY(18) : Read-only type(L_COMPO) for (t-fuel,den).
*         HENTRY(19) : Read-only type(L_COMPO) for (t-cool,den).
*         HENTRY(20) : Read-only type(L_COMPO) for power-up.
*         HENTRY(21) : Read-only type(L_COMPO) for power-in.
*         HENTRY(22) : Read-only type(L_COMPO) for power-d.
*         HENTRY(23) : Read-only type(L_COMPO) for reflec. ref.
*         HENTRY(24) : Read-only type(L_COMPO) for reflec. mod-up temp.
*         HENTRY(25) : Read-only type(L_COMPO) for reflec. mod-dw temp.
*         HENTRY(26) : Read-only type(L_COMPO) for reflec. mod-up den.
*         HENTRY(27) : Read-only type(L_COMPO) for reflec. mod-dw den.
*         HENTRY(28) : Read-only type(L_COMPO) for reflec. boron.
*         HENTRY(29) : Read-only type(L_COMPO) for reflec. purity.
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
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR)   IPLIST
      INTEGER       IOUT,NSTATE
      CHARACTER     NAMSBR*6
      PARAMETER    (IOUT=6,NSTATE=40,NAMSBR='CFC   ')
      CHARACTER     TEXT12*12,HSMG*131,TINFO*72,
     >              HSIGN*12,DBNAME*9,CTITRE*72
      INTEGER       ISTATE(NSTATE)
      INTEGER       IPRINT
      INTEGER       NBPARA
      PARAMETER    (NBPARA=18)
      REAL          DBPARA(NBPARA)
*-----
*  PARAMETER VALIDATION.
*-----
      IF(NENTRY.LE.28) CALL XABORT('CFC: 29 PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('CFC:'
     1 //'  MACROLIB LINKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('CFC: '
     1 //' MACROLIB IN CREATE OR MODIFICATION MODE EXPECTED.')
*-----
*     INDIVIDUAL LOCAL PARAMETER.
*-----
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)) CALL XABORT('CFC:'
     1 //' COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(2).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(3).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(4).NE.1).AND.(IENTRY(4).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(4).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(5).NE.1).AND.(IENTRY(5).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(5).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(6).NE.1).AND.(IENTRY(6).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(6).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(7).NE.1).AND.(IENTRY(7).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(7).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY'
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(8).NE.1).AND.(IENTRY(8).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(8).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(9).NE.1).AND.(IENTRY(9).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(9).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(10).NE.1).AND.(IENTRY(10).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(10).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(11).NE.1).AND.(IENTRY(11).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(11).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(12).NE.1).AND.(IENTRY(12).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(12).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY'
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(13).NE.1).AND.(IENTRY(13).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(13).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(14).NE.1).AND.(IENTRY(14).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(14).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
*-----
*     STURATING ISOTOPES
*-----
      IF((IENTRY(15).NE.1).AND.(IENTRY(15).NE.2)) CALL XABORT('CFC:'
     1 //' COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(15).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
      IF((IENTRY(16).NE.1).AND.(IENTRY(16).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(16).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
*-----
*      MIXED PARAMETERS
*-----
      IF((IENTRY(17).NE.1).AND.(IENTRY(17).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(17).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(18).NE.1).AND.(IENTRY(18).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(18).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(19).NE.1).AND.(IENTRY(19).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(19).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
*-----
*     HISTORY PARAMETERS (FISSILE ISOTOPES)
*-----
      IF((IENTRY(20).NE.1).AND.(IENTRY(20).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(20).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(21).NE.1).AND.(IENTRY(21).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(21).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(22).NE.1).AND.(IENTRY(22).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(22).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
*-----
*     MODERATOR PROPERTIES
*-----
      IF((IENTRY(23).NE.1).AND.(IENTRY(23).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(23).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(24).NE.1).AND.(IENTRY(24).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(24).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(25).NE.1).AND.(IENTRY(25).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(25).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(26).NE.1).AND.(IENTRY(26).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(26).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
       IF((IENTRY(27).NE.1).AND.(IENTRY(27).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(27).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(28).NE.1).AND.(IENTRY(28).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(28).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
      IF((IENTRY(29).NE.1).AND.(IENTRY(29).NE.2)) CALL XABORT('CFC: '
     1 //'COMPO LINKED LIST OR XSM FILE EXPECTED AT RHS.')
      IF(JENTRY(29).NE.2) CALL XABORT('CFC: COMPO IN READ-ONLY '
     1 //' MODE EXPECTED AT RHS.')
*-----
*  END OF L_COMPO FILES REQUIRED
*-----
      ITYPE=JENTRY(1)
      IPLIST=KENTRY(1)
      IPRINT=1
*-----
*  CHECK THE SIGNTURE OF THE LCM OBJECT
*-----
      DO 200 I=2,29
      CALL LCMGTC(KENTRY(I),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_COMPO') THEN
         TEXT12=HENTRY(I)
         CALL XABORT('CFC: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_COMPO EXPECTED.')
      ENDIF
 200  CONTINUE
*----
*  READ THE INPUT FILE
*----
      CALL CFCGET(TINFO,DBNAME,IPRINT,NBPARA,DBPARA)
*-----
*        CREATED OR MODIFIED LCM OBJECT TYPE(L_COMPO)
*-----
         IF(JENTRY(1).EQ.0) THEN
         HSIGN='REACTOR_XSDB'
         CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
*-----
*        STORE THE INFORMATION TITLE
*-----
         CALL LCMPTC(KENTRY(1),'INFORMATION',72,1,TINFO)
         ELSE
         CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN.NE.'REACTOR_XSDB') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('CFC: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. REACTOR_XSDB EXPECTED.')
            ENDIF
*-----
*        RECOVER THE INFORMATION TITLE
*-----
         CALL LCMGTC(KENTRY(1),'INFORMATION',72,1,CTITRE)
            IF(TINFO.NE.CTITRE) THEN
            CALL XABORT('CFC: '//CTITRE//' INFOR TITLE EXPECTED.')
            ENDIF
      ENDIF
*-----
*  RECOVER SOME INFORMATIONS FROM THE FIRST COMPO (NOMINAL).
*-----
      TEXT12='SIGNATURE'
      CALL LCMNXT(KENTRY(2),TEXT12)
      IF(TEXT12.EQ.'SIGNATURE') CALL XABORT('CFC: INVALID INPUT COMPO.')
      CALL LCMSIX(KENTRY(2),TEXT12,1)
*-----
*  RECOVER THE CPO TITLE.
*-----
      CALL LCMLEN(KENTRY(2),'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(KENTRY(2),'TITLE',72,1,CTITRE)
      ELSE
         CTITRE='*** NO TITLE PROVIDED ***'
      ENDIF
      CALL LCMSIX(KENTRY(2),' ',2)
*-----
*  READ PARAMETERS
*-----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(2),'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).GT.1) THEN
      WRITE(HSMG,'(28HNUMBER OF MIXUTRE SAVED IS =,I5,
     1 28H ONLY ONE REGION IS REQUIRED)') ISTATE(1)
      CALL XABORT(HSMG)
      ENDIF
      NGRP=ISTATE(2)
      NISO=ISTATE(3)
      NL=ISTATE(4)
      NBURN=ISTATE(5)
      NXS=21+NL
*-----
*  CHECK OTHERS CPO INFORMATIONS.
*-----
      DO 100 I=3,29
      TEXT12='SIGNATURE'
      CALL LCMNXT(KENTRY(I),TEXT12)
      IF(TEXT12.EQ.'SIGNATURE') CALL XABORT('CFC: INVALID '
     1 //'INPUT COMPO.')
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(I),'STATE-VECTOR',ISTATE)
*
*     ENERGY GROUP AND ORDER OF SCATTERING ANISOTROPY.
*
       IF(ISTATE(2).NE.NGRP) THEN
       WRITE(HSMG,'(7HNGRP = ,I5,13H IN REF NGRP=,I5)')
     1 ISTATE(2),NGRP
       CALL XABORT('CFC:  INCONSISTENT NB OF GROUPS '
     1 //TEXT12//' IS '//HSMG//' ')
       ELSE IF(ISTATE(4).LT.NL) THEN
       WRITE(HSMG,'(5HNL = ,I5,11H IN REF NL=,I5)')
     1 ISTATE(4),NL
       CALL XABORT('CFC:  INCONSISTENT NB OF LEGENDRE ORDERS.'
     1 //TEXT12//'IS '//HSMG//' ')
       ENDIF
 100   CONTINUE
*-----
*  CALL CFC DRIVER.
*-----
      CALL CFCDRV(IPRINT,NENTRY,KENTRY,HENTRY,NBURN,NGRP,NISO,NL,
     1 CTITRE,DBNAME,NXS,NBPARA,DBPARA)
      RETURN
      END
