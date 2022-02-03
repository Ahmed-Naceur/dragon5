*DECK M2T
      SUBROUTINE M2T(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
*  Recover information from a macrolib and translate the requested data
*  towards an Apotrim interface file.
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
*         HENTRY(1) create or modification ascii file containing
*         Apotrim data;
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
      TYPE(C_PTR) IPMAC
      PARAMETER (NSTATE=40,IOUT=6)
      CHARACTER TEXT12*12,TEXT20*20,HSIGN*12
      DOUBLE PRECISION DFLOTT
      INTEGER ISTATE(NSTATE)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBM,HBM
      REAL, ALLOCATABLE, DIMENSION(:) :: BUP,TEMP
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LE.1) CALL XABORT('M2T: MINIMUM OF 2 OBJECTS EXPECTED.')
      TEXT12=HENTRY(1)
      IF(IENTRY(1).NE.4) CALL XABORT('M2T: ASCII FILE NAMED '//TEXT12
     1 //' EXPECTED AT LHS.')
      IF(JENTRY(1).EQ.2) CALL XABORT('M2T: ASCII FILE IN CREATE OR MOD'
     1 //'IFICATION MODE EXPECTED.')
      LOUT=FILUNIT(KENTRY(1))
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)) CALL XABORT('M2T: LCM '
     1 //'OBJECT EXPECTED AT RHS.')
      IF(JENTRY(2).NE.2) CALL XABORT('M2T: LCM OBJECTS IN READ-ONLY MO'
     1 //'DE EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB') CALL XABORT('M2T: MACROLIB OBJECT EXPE'
     1 //'CTED AT RHS.')
      IPMAC=KENTRY(2)
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NGRP=ISTATE(1)
      MAXMIX=ISTATE(2)
      MAXNL=ISTATE(3)
      NBFIS=ISTATE(4)
      IF(NBFIS.GT.1) CALL XABORT('M2T: THE CAPABILITY TO MERGE MANY FI'
     1 //'SSION SPECTRA IS NOT IMPLEMENTED.')
*----
*  ALLOCATE MEMORY
*----
      ALLOCATE(NBM(MAXMIX),HBM(5*MAXMIX),BUP(MAXMIX),TEMP(MAXMIX))
      CALL XDISET(NBM,MAXMIX,1)
      CALL XDRSET(BUP,MAXMIX,0.0)
      CALL XDRSET(TEMP,MAXMIX,0.0)
*----
*  READ THE INPUT DATA
*----
      NL=1
      NBMIX=0
      ICTR=0
      IGMAIL=0
      IMPX=1
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('M2T: CHARACTER DATA EXPECTED(1).')
      IF(TEXT12.EQ.'EDIT') THEN
*        READ THE PRINT INDEX.
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('M2T: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT12.EQ.'MIX') THEN
*        READ A MATERIAL MIXTURE.
         TEXT20=' '
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT20,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('M2T: CHARACTER DATA EXPECTED(2).')
         NBMIX=NBMIX+1
         IF(NBMIX.GT.MAXMIX) CALL XABORT('M2T: MAXMIX OVERFLOW.')
         READ(TEXT20,'(5A4)') (HBM(5*(NBMIX-1)+I0),I0=1,5)
   30    CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('M2T: CHARACTER DATA EXPECTED(3).')
         IF(TEXT12.EQ.'FROM') THEN
            CALL REDGET(INDIC,NBM(NBMIX),FLOTT,TEXT12,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('M2T: INTEGER DATA EXPECTED(2).')
            GO TO 30
         ELSE IF(TEXT12.EQ.'BURN') THEN
*           READ A BURNUP.
            CALL REDGET(INDIC,NITMA,BUP(NBMIX),TEXT12,DFLOTT)
            IF(INDIC.NE.2) CALL XABORT('M2T: REAL DATA EXPECTED(1).')
            GO TO 30
         ELSE IF(TEXT12.EQ.'TEMP') THEN
*           READ A TEMPERATURE.
            CALL REDGET(INDIC,NITMA,TEMP(NBMIX),TEXT12,DFLOTT)
            IF(INDIC.NE.2) CALL XABORT('M2T: REAL DATA EXPECTED(2).')
            GO TO 30
         ELSE IF(TEXT12.NE.'ENDMIX') THEN
            CALL XABORT('M2T: FROM, BURN, TEMP OR ENDMIX EXPECTED.')
         ENDIF
      ELSE IF(TEXT12.EQ.'PN') THEN
*        READ THE ANISOTROPY ORDER
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT20,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('M2T: INTEGER DATA EXPECTED(3).')
         NL=MIN(NITMA+1,MAXNL)
      ELSE IF(TEXT12.EQ.'TRAN') THEN
*        PERFORM TRANSPORT CORRECTION
         ICTR=1
      ELSE IF(TEXT12.EQ.'NOMA') THEN
*        DO NOT WRITE ENERGY MESH ON APOTRIM FILE
         IGMAIL=1
      ELSE IF(TEXT12.EQ.';') THEN
         GO TO 40
      ELSE
         CALL XABORT('M2T: '//TEXT12//' IS AN INVALID KEYWORD.')
      ENDIF
      GO TO 20
*----
*  RECOVER INFORMATION
*----
   40 CALL M2TDRV(IMPX,LOUT,IPMAC,NGRP,NBMIX,MAXMIX,NL,NBFIS,ICTR,
     1 IGMAIL,BUP,TEMP,HBM,NBM)
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(TEMP,BUP,HBM,NBM)
      RETURN
      END
