*DECK SOUR
      SUBROUTINE SOUR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set the fixed external sources.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Bienvenue
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): creation type(L_SOURCE);
*         HENTRY(2): read-onlu type(L_MACROLIB);
*         HENTRY(3): read-only type(L_TRACKING);
*         HENTRY(4): read-only type(L_GEOM). 
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
      PARAMETER (IOUT=6,NSTATE=40,MAXPAR=10)
      TYPE(C_PTR) IPSOUR,IPMAC,IPTRK,IPGEOM
      CHARACTER HSIGN*12,TEXT12*12
      INTEGER ISTATE(NSTATE)
      REAL ISOUR
      DOUBLE PRECISION DFLOTT
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,YYY,ZZZ
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.4) CALL XABORT('SOUR: FOUR PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('SOUR: LI'
     1 //'NKED LIST OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('SOUR: ENTRY IN CREATE MODE EXPE'
     1 //'CTED.')
      IPSOUR=KENTRY(1)
      DO I=2,NENTRY
        IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)) CALL XABORT('SOUR: '
     1  //'LINKED LIST OR XSM FILE EXPECTED AT RHS.')
        IF(JENTRY(I).NE.2) CALL XABORT('SOUR: ENTRY IN READ-ONLY MODE'
     1  //' EXPECTED.')
      ENDDO
*----
*  RECOVER MACROLIB INFORMATION
*----
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_MACROLIB') THEN
        IPMAC=KENTRY(2)
      ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
        IPMAC=LCMGID(KENTRY(2),'MACROLIB')
      ELSE
        TEXT12=HENTRY(2)
        CALL XABORT('SOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_MACROLIB OR L_LIBRARY EXPECTED')
      ENDIF
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
      NMAT=ISTATE(2)
*----
*  RECOVER TRACKING INFORMATION
*----
      IPTRK=KENTRY(3)
      CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
        TEXT12=HENTRY(3)
        CALL XABORT('SOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NDIM=ISTATE(9)
      LX=ISTATE(12)
      LY=ISTATE(13)
      LZ=ISTATE(14)
      NLF=ISTATE(15)
      IF(ISTATE(4).NE.NMAT) CALL XABORT('SOUR: INVALID NUMBER OF MI'
     1 //'XTURES.')
*----
*  RECOVER GEOMETRY INFORMATION
*----
      IPGEOM=KENTRY(4)
      CALL LCMGTC(IPGEOM,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_GEOM') THEN
        TEXT12=HENTRY(4)
        CALL XABORT('SOUR: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_GEOM EXPECTED.')
      ENDIF

      ALLOCATE(XXX(LX+1),YYY(LY+1),ZZZ(LZ+1))

      CALL LCMGET(IPGEOM,'MESHX',XXX)
      IF(NDIM.GE.2) CALL LCMGET(IPGEOM,'MESHY',YYY)
      IF(NDIM.GE.3) CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
*----
*  READ THE INPUT DATA
*----

! PARAMETERS VALUES
      STYPE=-1    
      ISOUR=1.0   
      GN=1         
      XMIN=0.0
      XMAX=0.0
      ISXDEF=0
      YMIN=0.0
      YMAX=0.0
      ISYDEF=0
      ZMIN=0.0
      ZMAX=0.0
      ISZDEF=0
      XB=0
      YB=0
      ZB=0
      ! direction parameter (mu,nu,eta? 2 angles?)

  10  CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('SOUR: CHARACTER DATA EXPECTED.')
! ISOTROPIC SOURCE KEYWORD
      IF(TEXT12.EQ.'ISO') THEN
        STYPE=0
! MONODIRECTIONNAL BOUNDARY SOURCE KEYWORD
      ELSE IF(TEXT12.EQ.'MONO') THEN
        STYPE=1
! SOURCE INTENSITY KEYWORD
      ELSE IF(TEXT12.EQ.'IN') THEN
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        ISOUR=FLOTT
! ENERGY GROUP NUMBER KEYWORD
      ELSE IF(TEXT12.EQ.'GN') THEN
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('SOUR: INTEGER DATA EXPECTED.')
        GN=NITMA
        IF(GN.GT.NG.OR.GN.LT.1) CALL XABORT('SOUR: INVALID GROUP NUMBER'
     1  //'.')
! (X,Y,Z) LIMITS KEYWORDS
      ELSE IF(TEXT12.EQ.'XLIM') THEN
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        XMIN=FLOTT
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        XMAX=FLOTT
        ISXDEF=1
      ELSE IF(TEXT12.EQ.'YLIM') THEN       
        IF(NDIM.LT.2) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 2 ')
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        YMIN=FLOTT
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        YMAX=FLOTT
        ISYDEF=1
      ELSE IF(TEXT12.EQ.'ZLIM') THEN
        IF(NDIM.LT.3) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 3 ')
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        ZMIN=FLOTT
        CALL REDGET(INDIC,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('SOUR: REAL DATA EXPECTED.')
        ZMAX=FLOTT
        ISZDEF=1
! BOUNDARY SOURCE LOCATION KEYWORD
      ELSE IF(TEXT12.EQ.'X-') THEN
        XB=-1
      ELSE IF(TEXT12.EQ.'X+') THEN
        XB=1
      ELSE IF(TEXT12.EQ.'Y-') THEN
        IF(NDIM.LT.2) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 2 ')
        YB=-1
      ELSE IF(TEXT12.EQ.'Y+') THEN
        IF(NDIM.LT.2) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 2 ')
        YB=1
      ELSE IF(TEXT12.EQ.'Z-') THEN
        IF(NDIM.LT.3) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 3 ')
        ZB=-1
      ELSE IF(TEXT12.EQ.'Z+') THEN
        IF(NDIM.LT.3) CALL XABORT('SOUR: INVALID USE OF YLIM, DIM < 3 ')
        ZB=1
! END OF INPUT DATA
      ELSE IF(TEXT12.EQ.';') THEN
        GO TO 20
! ERROR IF KEYWORD IS UNRECONIZED
      ELSE
        CALL XABORT('SOUR: '//TEXT12//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
  20  CONTINUE

! OTHER VERIFICATIONS

      IF(ISXDEF.EQ.1.AND.(XMIN.LT.XXX(1).OR.XMAX.GT.XXX(LX+1))) THEN
        CALL XABORT('SOUR: SOURCE X-LIMITS OUTSIDE THE GEOMETRY.')
      ENDIF
      IF(NDIM.GE.2) THEN
      IF(ISYDEF.EQ.1.AND.(YMIN.LT.YYY(1).OR.YMAX.GT.YYY(LX+1))) THEN
        CALL XABORT('SOUR: SOURCE Y-LIMITS OUTSIDE THE GEOMETRY.')
      ENDIF
      IF(NDIM.GE.3) THEN
      IF(ISZDEF.EQ.1.AND.(ZMIN.LT.ZZZ(1).OR.ZMAX.GT.ZZZ(LX+1))) THEN
        CALL XABORT('SOUR: SOURCE Z-LIMITS OUTSIDE THE GEOMETRY.')
      ENDIF
      ENDIF
      ENDIF


      IF(STYPE.EQ.0) THEN
      IF(XMIN.GE.XMAX.OR.ISXDEF.EQ.0) CALL XABORT('SOUR: SOURCE VOLUME'
     1 //' X-BOUNDARY NOT DEFINED OR INVALID.')
      IF(NDIM.GE.2) THEN
      IF(YMIN.GE.YMAX.OR.ISYDEF.EQ.0) CALL XABORT('SOUR: SOURCE VOLUME'
     1 //' Y-BOUNDARY NOT DEFINED OR INVALID.')
      IF(NDIM.GE.3) THEN
      IF(ZMIN.GE.ZMAX.OR.ISZDEF.EQ.0) CALL XABORT('SOUR: SOURCE VOLUME'
     1 //' Z-BOUNDARY NOT DEFINED OR INVALID.')
      ENDIF
      ENDIF
      ENDIF

      IF(STYPE.EQ.1) THEN
      IF(XB.EQ.0.AND.YB.EQ.0.AND.ZB.EQ.0) CALL XABORT('SOUR: BOUNDARY'
     1 //' SOURCE LOCATION IS NEEDED.')
      IF(XB.NE.0.AND.ISXDEF.EQ.1) CALL XABORT('SOUR: CONFLICTING INPUT'
     1 //': NO X-ORIENTED THICKNESS SHOULD BE DEFINE FOR X-BOUNDARY'
     2 //' SOURCES.')
      IF(YB.NE.0.AND.ISYDEF.EQ.1) CALL XABORT('SOUR: CONFLICTING INPUT'
     1 //': NO Y-ORIENTED THICKNESS SHOULD BE DEFINE FOR Y-BOUNDARY'
     2 //' SOURCES.')
      IF(ZB.NE.0.AND.ISZDEF.EQ.1) CALL XABORT('SOUR: CONFLICTING INPUT'
     1 //': NO Z-ORIENTED THICKNESS SHOULD BE DEFINE FOR Z-BOUNDARY'
     2 //' SOURCES.')
      IF(XB.EQ.0.AND.ISXDEF.EQ.0) THEN
        XMIN=XXX(1)
        XMAX=XXX(LX+1)
      ENDIF
      IF(NDIM.GE.2) THEN
      IF(YB.EQ.0.AND.ISYDEF.EQ.0) THEN
        YMIN=YYY(1)
        YMAX=YYY(LY+1)
      ENDIF
      IF(NDIM.GE.3) THEN
      IF(ZB.EQ.0.AND.ISZDEF.EQ.0) THEN
        ZMIN=ZZZ(1)
        ZMAX=ZZZ(LZ+1)
      ENDIF
      ENDIF
      ENDIF
      ENDIF

*----
*  COMPUTE THE FIXED EXTERNAL SOURCES
*----

      IF(STYPE.EQ.0) THEN
*        CALL SOURISO()
      ELSE IF(STYPE.EQ.1) THEN
*        CALL SOURMONO()
      ELSE
        CALL XABORT('SOUR: SOURCE TYPE EXPECTED.')
      ENDIF

*----
*  SAVE THE FIXED EXTERNAL SOURCES ON LCM
*----


*----
*  SAVE THE SIGNATURE AND STATE VECTOR
*----






      END
