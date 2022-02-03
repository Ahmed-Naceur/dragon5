*DECK MACD
      SUBROUTINE MACD(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Macroscopic cross sections and diffusion coefficients input module.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
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
*         HENTRY(1) : create or modification type(L_MACROLIB).
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
      CHARACTER TEXT4*4,TEXT12*12,HSMG*131,HSIGN*12
      DOUBLE PRECISION DFLOTT
      INTEGER IPAR(NSTATE)
      TYPE(C_PTR) IPLIST
      REAL, DIMENSION(:,:), ALLOCATABLE :: ALBP
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.EQ.0) CALL XABORT('MACD: PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('MACD: LCM'
     1 //' OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('MACD: ENT'
     1 //'RY IN CREATE OR MODIFICATION MODE EXPECTED.')
      ITYPE=JENTRY(1)
      IPLIST=KENTRY(1)
*----
*  READ THE INPUT DATA.
*----
*     DEFAULT OPTIONS:
      IND=1
      IMPX=1
      ISTEP=0
      IF(ITYPE.EQ.0) THEN
         NL=1
         NGRP=0
         NMIXT=0
         NIFISS=1
         NDG=0
         NALBP=0
         NSTEP=0
         IF(NENTRY.EQ.2) THEN
            IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)) CALL XABORT('MACD'
     1      //': LCM OBJECT EXPECTED AT RHS.')
            IF(JENTRY(2).NE.2) CALL XABORT('MACD: RHS ENTRY IN READ-ONL'
     1      //'Y MODE EXPECTED.')
            CALL LCMEQU(KENTRY(2),IPLIST)
            IND=2
         ENDIF
      ELSE IF(ITYPE.EQ.1) THEN
         IND=2
      ENDIF
      IF(IND.EQ.2) THEN
         CALL LCMGTC(IPLIST,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_MACROLIB') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('MACD: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_MACROLIB EXPECTED.')
         ENDIF
         IND=2
         CALL LCMGET(IPLIST,'STATE-VECTOR',IPAR)
         NGRP=IPAR(1)
         NMIXT=IPAR(2)
         NL=IPAR(3)
         NIFISS=IPAR(4)
         NDG=IPAR(7)
         NALBP=IPAR(8)
         NSTEP=IPAR(11)
      ENDIF
*----
*  READ THE MAC: MODULE OPTIONS.
*----
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('MACD: CHARACTER DATA EXPECTED(1).')
   20 IF(TEXT4.EQ.'EDIT') THEN
*        READ THE PRINT INDEX.
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'NGRO') THEN
*        READ THE NUMBER OF ENERGY GROUPS.
         IF(IND.EQ.2) CALL XABORT('MACD: NGRO IS ALREADY DEFINED.')
         CALL REDGET(INDIC,NGRP,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(2).')
      ELSE IF(TEXT4.EQ.'NMIX') THEN
*        READ THE MAXIMUM NUMBER OF MATERIAL MIXTURES.
         IF(IND.EQ.2) CALL XABORT('MACD: NMIX IS ALREADY DEFINED.')
         CALL REDGET(INDIC,NMIXT,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(3).')
      ELSE IF(TEXT4.EQ.'DELP') THEN
*        READ THE MAXIMUM NUMBER OF PRECURSORS.
         IF(IND.EQ.2) CALL XABORT('MACD: DELP IS ALREADY DEFINED.')
         CALL REDGET(INDIC,NDG,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(3).')
      ELSE IF(TEXT4.EQ.'ANIS') THEN
*        READ THE SCATTERING ANISOTROPY FOR TRANSPORT THEORY CASES.
         IF(IND.EQ.2) CALL XABORT('MACD: NMIX IS ALREADY DEFINED.')
         CALL REDGET(INDIC,NL,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(4).')
      ELSE IF(TEXT4.EQ.'NIFI') THEN
*        READ THE NUMBER OF FISSILE ISOTOPES
         IF(IND.EQ.2) CALL XABORT('MACD: NIFISS IS ALREADY DEFINED.')
         CALL REDGET(INDIC,NIFISS,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(5).')
      ELSE IF(TEXT4.EQ.'ALBP') THEN
*        READ GROUP-INDEPENDENT PHYSICAL ALBEDOS
         CALL REDGET(INDIC,NALBP,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(6).')
         IF(NALBP.GT.0) THEN
           ALLOCATE(ALBP(NALBP,NGRP))
           DO IAL=1,NALBP
             DO IGR=1,NGRP
               CALL REDGET(INDIC,NITMA,ALBP(IAL,IGR),TEXT4,DFLOTT)
               IF(INDIC.NE.2) CALL XABORT('MACD: REAL DATA EXPECTED.')
             ENDDO
           ENDDO
           CALL LCMPUT(IPLIST,'ALBEDO',NALBP*NGRP,2,ALBP)
           DEALLOCATE(ALBP)
         ELSE
           CALL XABORT('MACD: INVALID NUMBER OF ALBEDOS.')
         ENDIF
         IF(ITYPE.EQ.1) THEN
           CALL LCMGET(IPLIST,'STATE-VECTOR',IPAR)
           IPAR(8)=NALBP
           CALL LCMPUT(IPLIST,'STATE-VECTOR',NSTATE,1,IPAR)
         ENDIF
      ELSE IF(TEXT4.EQ.'STEP') THEN
*        STEP TO A SON DIRECTORY AND WRITE PERTURBATION VALUES.
         CALL REDGET(INDIC,ISTEP,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('MACD: INTEGER DATA EXPECTED(7).')
         WRITE(TEXT12,'(4HSTEP,I8)') ISTEP
         IF(IND.EQ.1) THEN
            CALL LCMLEN(IPLIST,TEXT12,ILENG,ITYLCM)
            IF(ILENG.GT.0) THEN
               WRITE(HSMG,'(30HMACD: PERTURBATION DIRECTORY '',A12,
     1         21H'' ALREADY EXISTS IN '',A12,2H''.)') TEXT12,HENTRY(1)
               CALL XABORT(HSMG)
            ENDIF
         ENDIF
         NSTEP=MAX(NSTEP,ISTEP)
         CALL LCMSIX(IPLIST,TEXT12,1)
         IF(IMPX.GT.0) WRITE(6,'(/34H MACD: WRITE PERTURBATION VALUES O,
     1   13HN DIRECTORY '',A12,6H'' OF '',A12,2H''.)') TEXT12,HENTRY(1)
      ELSE IF(TEXT4.EQ.'READ') THEN
*        INPUT NON-PERTURBED OR PERTURBED DIFFUSION COEFFICIENTS AND
*        CROSS SECTIONS PER MIXTURE.
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF((INDIC.NE.3).OR.(TEXT4.NE.'INPU')) CALL XABORT('MACD: INPU'
     1   //'T KEYWORD EXPECTED.')
         CALL MACXSI(IPLIST,IND,NMIXT,NGRP,NDG,NL,IMPX,NBMIX,JND)
         IF(ISTEP.GT.0) THEN
            IF(IMPX.GT.1) CALL LCMLIB(IPLIST)
            CALL LCMSIX(IPLIST,' ',2)
         ENDIF
         IF(JND.EQ.1) THEN
            GO TO 40
         ELSE IF(JND.EQ.2) THEN
            TEXT4='STEP'
            GO TO 20
         ENDIF
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 40
      ELSE
         CALL XABORT('MACD: '//TEXT4//' IS AN INVALID KEY-WORD.')
      ENDIF
      GO TO 10
*
   40 IF(ITYPE.EQ.0) THEN
         HSIGN='L_MACROLIB'
         CALL LCMPTC(IPLIST,'SIGNATURE',12,1,HSIGN)
         CALL XDISET(IPAR,NSTATE,0)
         IPAR(1)=NGRP
         IPAR(2)=NMIXT
         IPAR(3)=NL
         IPAR(4)=NIFISS
         IPAR(5)=0
         IPAR(6)=0
         IPAR(7)=NDG
         IPAR(8)=NALBP
         IPAR(11)=NSTEP
         CALL LCMPUT(IPLIST,'STATE-VECTOR',NSTATE,1,IPAR)
      ENDIF
      IF(IMPX.GT.1) CALL LCMLIB(IPLIST)
      RETURN
      END
