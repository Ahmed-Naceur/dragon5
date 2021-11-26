*DECK NSSF
      SUBROUTINE NSSF(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux solution for the nodal expansion method (NEM).
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
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
*         HENTRY(1): create type(L_FLUX) nodal flux;
*         HENTRY(2): read-only type(L_TRACK) nodal tracking;
*         HENTRY(3): read-only type(L_MACROLIB) nodal macrolib.
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
      TYPE(C_PTR) IPFLX,IPTRK,IPMAC
      CHARACTER HSIGN*12,TEXT4*4,TEXT12*12,HSMG*131
      LOGICAL LNODF
      INTEGER ISTATE(NSTATE)
      REAL REALIR
      DOUBLE PRECISION DBLLIR
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ITRIAL
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.NE.3) CALL XABORT('NSSF: 3 PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('NSSF: LCM'
     1 //' OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('NSSF: ENTRY IN CREATE MODE EXPEC'
     1 //'TED.')
      HSIGN='L_FLUX'
      IPFLX=KENTRY(1)
      CALL LCMPTC(IPFLX,'SIGNATURE',12,1,HSIGN)
      DO IEN=2,3
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2)) CALL XABORT('NSS'
     1  //'F: LCM OBJECT EXPECTED AT LHS.')
        IF(JENTRY(IEN).NE.2) CALL XABORT('NSSF: ENTRY IN READ-ONLY MOD'
     1  //'E EXPECTED.')
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
        TEXT12=HENTRY(IEN)
        IF(IEN.EQ.2) THEN
          IF(HSIGN.NE.'L_TRACK') THEN
            CALL XABORT('NSSF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_TRACK EXPECTED.')
          ENDIF
          IPTRK=KENTRY(2)
        ELSE IF(IEN.EQ.3) THEN
          IF(HSIGN.NE.'L_MACROLIB') THEN
            CALL XABORT('NSSF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_MACROLIB EXPECTED.')
          ENDIF
          IPMAC=KENTRY(3)
        ENDIF
      ENDDO
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,TEXT12)
      IF(TEXT12.NE.'TRIVAC') CALL XABORT('NSSF: TRIVAC TRACKING EXPECT'
     1 //'ED.')
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      LX1=ISTATE(1)
      NMIX=ISTATE(4)
      IGMAX=ISTATE(39)
      IF(ISTATE(6).NE.2) CALL XABORT('NSSF: 1D SLAB GEOMETRY EXPECTED.')
      IF(ISTATE(12).NE.4) CALL XABORT('NSSF: NEM DISCRETIZATION EXPEC'
     1 //'TED.')
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      IF(ISTATE(2).NE.NMIX) THEN
        WRITE(HSMG,'(39HNSSF: INVALID NUMBER OF MIXTURES (GEOM=,I5,
     1  10H MACROLIB=,I5,2H).)') NMIX,ISTATE(2)
        CALL XABORT(HSMG)
      ENDIF
*---
*  READ DATA
*---
      ALLOCATE(ITRIAL(NMIX,NG))
      IPRINT=1
      MAXOUT=1000
      EPSOUT=1.0E-6
      LNODF=.FALSE.
      BB2=0.0
      ITRIAL(:,:)=1
   10 CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
      IF(ITYPLU.EQ.10) GO TO 100
   20 IF(ITYPLU.NE.3) CALL XABORT('NSSF: READ ERROR - CHARACTER VARIAB'
     > //'LE EXPECTED')
      IF(TEXT4.EQ.';') THEN
        GO TO 100
      ELSE IF(TEXT4.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IPRINT,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('NSSF: READ ERROR - INTEGER VARIAB'
     >  //'LE EXPECTED')
      ELSE IF(TEXT4.EQ.'EXTE') THEN
   30   CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.EQ.1) THEN
          MAXOUT=INTLIR
        ELSE IF(ITYPLU.EQ.2) THEN
          EPSOUT=REALIR
        ELSE
          GO TO 20
        ENDIF
        GO TO 30
      ELSE IF(TEXT4.EQ.'NODF') THEN
        LNODF=.TRUE.
      ELSE IF(TEXT4.EQ.'BUCK') THEN
        CALL REDGET(ITYPLU,INTLIR,BB2,TEXT4,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('NSSF: READ ERROR - REAL VARIABLE '
     >  //'EXPECTED')
      ELSE
        CALL XABORT('NSSF: ILLEGAL KEYWORD '//TEXT4)
      ENDIF
      GO TO 10
  100 IF(IGMAX.GT.NG) CALL XABORT('NSSF: IGMAX>NG.')
      IF(IGMAX.GT.0) ITRIAL(:NMIX,IGMAX:NG)=2
      CALL NSSDRV(IPTRK,IPMAC,IPFLX,NG,LX1,NMIX,ITRIAL,EPSOUT,MAXOUT,
     1 LNODF,BB2,IPRINT)
      DEALLOCATE(ITRIAL)
      RETURN
      END
