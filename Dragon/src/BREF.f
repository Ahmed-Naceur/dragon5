*DECK BREF
      SUBROUTINE BREF(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the discontinuity factors in a 1D reflector model.
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
*         HENTRY(1): create type(L_GEOM) nodal geometry;
*         HENTRY(2): create type(L_MACROLIB) nodal macrolib;
*         HENTRY(3): read-only type(L_GEOM) sn geometry;
*         HENTRY(4): read-only type(L_EDITION) sn edition.
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
      TYPE(C_PTR) IPGEO1,IPMAC1,IPGEO2
      CHARACTER HSIGN*12,TEXT4*4,TEXT12*12,HSMG*131,HMREFL*12
      INTEGER ISTATE(NSTATE)
      REAL REALIR
      DOUBLE PRECISION DBLLIR
      LOGICAL LALB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITRIAL,IMIX1,IGAP
      REAL, ALLOCATABLE, DIMENSION(:) :: ADFREF
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPEDI2
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LT.4) CALL XABORT('BREF: >=4 PARAMETERS EXPECTED.')
      NC=NENTRY-3
      ALLOCATE(IPEDI2(NC))
      DO IEN=1,2
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2)) CALL XABORT('BREF'
     1  //': LCM OBJECT EXPECTED AT LHS.')
        IF(JENTRY(IEN).NE.0) CALL XABORT('BREF: ENTRY IN CREATE MODE EX'
     1  //'PECTED.')
        IF(IEN.EQ.1) THEN
          HSIGN='L_GEOM'
          IPGEO1=KENTRY(1)
        ELSE IF(IEN.EQ.2) THEN
          HSIGN='L_MACROLIB'
          IPMAC1=KENTRY(2)
        ENDIF
        CALL LCMPTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
      ENDDO
      DO IEN=3,NENTRY
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2)) CALL XABORT('BREF'
     1  //': LCM OBJECT EXPECTED AT LHS.')
        IF(JENTRY(IEN).NE.2) CALL XABORT('BREF: ENTRY IN READ-ONLY MODE'
     1  //' EXPECTED.')
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
        TEXT12=HENTRY(IEN)
        IF(IEN.EQ.3) THEN
          IF(HSIGN.NE.'L_GEOM') THEN
            CALL XABORT('BREF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_GEOM EXPECTED.')
          ENDIF
          IPGEO2=KENTRY(3)
        ELSE IF(IEN.GE.4) THEN
          IF(HSIGN.NE.'L_EDIT') THEN
            CALL XABORT('BREF: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_EDIT EXPECTED.')
          ENDIF
          IPEDI2(IEN-3)=KENTRY(IEN)
        ENDIF
      ENDDO
      CALL LCMGET(IPEDI2(1),'STATE-VECTOR',ISTATE)
      NMIX2=ISTATE(1)
      NG=ISTATE(2)
*---
*  READ DATA
*---
      ALLOCATE(ITRIAL(NG),ADFREF(NG))
      IPRINT=1
      ITRIAL(:)=1
      HMREFL=' '
      ISPH=0
      LX1=0
      LALB=.TRUE.
      NGET=0
   10 CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
      IF(ITYPLU.EQ.10) CALL XABORT('BREF: MISSING USER DATA.')
      IF(ITYPLU.NE.3) CALL XABORT('BREF: READ ERROR - CHARACTER VARIAB'
     > //'LE EXPECTED')
   20 IF(TEXT4.EQ.';') THEN
        GO TO 100
      ELSE IF(TEXT4.EQ.'EDIT') THEN
        CALL REDGET(ITYPLU,IPRINT,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('BREF: READ ERROR - INTEGER VARIAB'
     >  //'LE EXPECTED')
      ELSE IF(TEXT4.EQ.'HYPE') THEN
        CALL REDGET(ITYPLU,IGMAX,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT('BREF: READ ERROR - INTEGER VARIAB'
     >  //'LE EXPECTED')
        IF(IGMAX.LE.0) CALL XABORT('BREF: IGMAX<=0.')
        IF(IGMAX.GT.NG) CALL XABORT('BREF: IGMAX>NG.')
        ITRIAL(IGMAX:NG)=2
      ELSE IF(TEXT4.EQ.'MIX') THEN
        ALLOCATE(IMIX1(NMIX2))
   30   CALL REDGET(ITYPLU,INTLIR,REALIR,TEXT4,DBLLIR)
        IF(ITYPLU.EQ.1) THEN
          LX1=LX1+1
          IF(LX1.GT.NMIX2) CALL XABORT('BREF: LX1 OVERFLOW.')
          IF(INTLIR.GT.NMIX2) THEN
            WRITE(HSMG,'(12HBREF: IMIX1=,I5,10H > NMIX2=,I5,1H.)')
     >      INTLIR,NMIX2
            CALL XABORT(HSMG)
          ENDIF
          IMIX1(LX1)=INTLIR
          GO TO 30
        ELSE IF(ITYPLU.EQ.3) THEN
          GO TO 20
        ELSE
          CALL XABORT('BREF: READ ERROR - INTEGER OR CHARACTER VARIABL'
     >    //'E EXPECTED')
        ENDIF
      ELSE IF(TEXT4.EQ.'GAP') THEN
        ALLOCATE(IGAP(LX1))
        DO IBM1=1,LX1
          CALL REDGET(ITYPLU,IGAP(IBM1),REALIR,TEXT4,DBLLIR)
          IF(ITYPLU.NE.1) CALL XABORT('BREF: READ ERROR - INTEGER VARI'
     >    //'ABLE EXPECTED')
          IF(IGAP(IBM1).GT.NMIX2) THEN
            WRITE(HSMG,'(11HBREF: IGAP=,I5,10H > NMIX2=,I5,1H.)')
     >      IGAP(IBM1),NMIX2
            CALL XABORT(HSMG)
          ENDIF
        ENDDO
      ELSE IF(TEXT4.EQ.'MODE') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,HMREFL,DBLLIR)
        IF(ITYPLU.NE.3) CALL XABORT('BREF: READ ERROR - CHARACTER VARI'
     >  //'ABLE EXPECTED')
      ELSE IF(TEXT4.EQ.'SPH') THEN
        ISPH=1
      ELSE IF(TEXT4.EQ.'NOAL') THEN
        LALB=.FALSE.
      ELSE IF(TEXT4.EQ.'NGET') THEN
        NGET=1
        DO IGR=1,NG
          CALL REDGET(ITYPLU,INTLIR,ADFREF(IGR),TEXT4,DBLLIR)
          IF(ITYPLU.EQ.2) THEN
            CYCLE
          ELSE IF(ITYPLU.EQ.3) THEN
            NGET=2
            GO TO 20
          ELSE
            CALL XABORT('BREF: READ ERROR - REAL OR CHARACTER VARIABLE'
     >      //' EXPECTED')
          ENDIF
        ENDDO
      ELSE
        CALL XABORT('BREF: ILLEGAL KEYWORD '//TEXT4)
      ENDIF
      GO TO 10
  100 NMIX1=0
      DO IBM1=1,LX1
        IF(IMIX1(IBM1).NE.0) NMIX1=NMIX1+1
      ENDDO
      CALL BREDRV(NC,IPGEO1,IPMAC1,IPGEO2,IPEDI2,NG,LX1,NMIX1,NMIX2,
     > ITRIAL,IMIX1,IGAP,HMREFL,ISPH,LALB,NGET,ADFREF,IPRINT)
      DEALLOCATE(IMIX1,IGAP,ADFREF,ITRIAL,IPEDI2)
      IF(IPRINT.GT.0) THEN
        CALL LCMGET(IPMAC1,'STATE-VECTOR',ISTATE)
        WRITE(6,110) IPRINT,(ISTATE(I),I=1,9),ISTATE(12)
      ENDIF
      RETURN
*
  110 FORMAT(/17H MACROLIB OPTIONS/17H ----------------/
     1 7H IMPX  ,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H NGROUP,I6,28H   (NUMBER OF ENERGY GROUPS)/
     3 7H NBMIX ,I6,39H   (NUMBER OF MIXTURES IN THE MACROLIB)/
     4 7H NANISO,I6,34H   (MAXIMUM SCATTERING ANISOTROPY)/
     5 7H NIFISS,I6,45H   (MAXIMUM NUMBER OF FISSILE ISOTOPES IN A M,
     6 7HIXTURE)/
     7 7H NEDMAC,I6,34H   (NUMBER OF CROSS SECTION EDITS)/
     8 7H ITRANC,I6,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     9 43H=RECOVER FROM LIBRARY/4=LEAKAGE CORRECTION)/
     1 7H NLG   ,I6,39H   (NUMBER OF DELAYED PRECURSOR GROUPS)/
     2 7H NALB  ,I6,31H   (NUMBER OF PHYSICAL ALBEDOS)/
     3 7H ILEAK ,I6,40H   (1=DIFF AVAILABLE; 2=DIFFX AVAILABLE)/
     4 7H IDF   ,I6,42H   (=0/2/3 ADF INFORMATION ABSENT/PRESENT))
      END
