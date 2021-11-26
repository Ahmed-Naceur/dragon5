*DECK DSET
      SUBROUTINE DSET(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Set new parameters for the user-selected devices and/or for the
* groups of devices.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
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
      PARAMETER(NSTATE=40,IOUT=6)
      INTEGER ISTATE(NSTATE),RGRP
      DOUBLE PRECISION DFLOT
      CHARACTER TEXT*12,HSIGN*12
      LOGICAL LROD
      TYPE(C_PTR) IPDEV
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.GT.1) CALL XABORT('@DSET: ONE PARAMETER ALLOWED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('@DSET:'
     1 //' LCM OBJECT EXPECTED AT LHS.')
      CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_DEVICE') THEN
        TEXT=HENTRY(1)
        CALL XABORT('@DSET: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_DEVICE EXPECTED.')
      ENDIF
      IF(JENTRY(1).NE.1) CALL XABORT('@DSET: MODIFICATION MODE EX'
     1 //'PECTED FOR L_DEVICE.')
      IPDEV=KENTRY(1)
*----
*  RECOVER INFORMATION
*----
      CALL LCMGET(IPDEV,'STATE-VECTOR',ISTATE)
      IGEO=ISTATE(1)
      IF(IGEO.NE.7) CALL XABORT('@DSET: ONLY 3D-CARTESIAN GEOMETRY ALL'
     1 //'OWED.')
      NROD=ISTATE(2)
      RGRP=ISTATE(3)
      NLZC=ISTATE(4)
      LGRP=ISTATE(5)
      IMODE=ISTATE(6)
      IF((IMODE.EQ.0).AND.(NROD.GT.0)) CALL XABORT('@DSET: IMODE NOT S'
     1 //'ET.')
*     READ PRINTING INDEX
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3) CALL XABORT('@DSET: CHARACTER DATA EXPECTED.')
      IF(TEXT.NE.'EDIT') CALL XABORT('@DSET: KEYWORD EDIT EXPECTED.')
      CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1) CALL XABORT('@DSET: INTEGER FOR EDIT EXPECTED.')
      NDEV=0
      NGRP=0
   10 NDEV=NDEV+1
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
*----
*  ROD OPTION
*----
      IF(TEXT.EQ.'ROD') THEN
        CALL REDGET(ITYP,ID,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1) CALL XABORT('@DSET: INTEGER ROD-ID NUMB'
     1   //'ER EXPECTED.')
        IF((ID.GT.NROD).OR.(ID.EQ.0)) THEN
          WRITE(IOUT,*)'@DSET: READ CURRENT ROD-ID #',ID
          CALL XABORT('@DSET: WRONG ROD-ID NUMBER.')
        ENDIF
        IF(IMPX.GT.0) WRITE(IOUT,1000)ID
        LROD=.TRUE.
        CALL DSET1D(IPDEV,IMODE,ID,LROD,IMPX)
*----
*  LZC OPTION
*----
      ELSEIF(TEXT.EQ.'LZC') THEN
        CALL REDGET(ITYP,ID,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@DSET: INTEGER LZC-ID NUMB'
     1   //'ER EXPECTED.')
        IF((ID.GT.NLZC).OR.(ID.EQ.0)) THEN
          WRITE(IOUT,*)'@DSET: READ CURRENT LZC-ID #',ID
          CALL XABORT('@DSET: WRONG LZC-ID NUMBER.')
        ENDIF
        IF(IMPX.GT.0) WRITE(IOUT,1001)ID
        LROD=.FALSE.
        CALL DSET1D(IPDEV,IMODE,ID,LROD,IMPX)
*----
*  ROD-GROUP OPTION
*----
      ELSEIF(TEXT.EQ.'ROD-GROUP') THEN
        CALL REDGET(ITYP,IGRP,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@DSET: INTEGER GROUP-ID NUM'
     1   //'BER EXPECTED.')
        IF((IGRP.GT.RGRP).OR.(IGRP.LE.0)) THEN
          WRITE(IOUT,*)'@DSET: READ CURRENT GROUP-ID #',IGRP
          CALL XABORT('@DSET: WRONG ROD GROUP-ID NUMBER.')
        ENDIF
        IF(IMPX.GT.0) WRITE(IOUT,1002)IGRP
        LROD=.TRUE.
        CALL DSETGR(IPDEV,IMODE,IGRP,NDGR,LROD,IMPX)
        NDEV=NDEV+NDGR-1
        NGRP=NGRP+1
*----
*  LZC-GROUP OPTION
*----
      ELSEIF(TEXT.EQ.'LZC-GROUP') THEN
        CALL REDGET(ITYP,IGRP,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@DSET: INTEGER GROUP-ID NUM'
     1   //'BER EXPECTED.')
        IF((IGRP.GT.LGRP).OR.(IGRP.LE.0)) THEN
          WRITE(IOUT,*)'@DSET: READ CURRENT GROUP-ID #',IGRP
          CALL XABORT('@DSET: WRONG LZC GROUP-ID NUMBER.')
        ENDIF
        IF(IMPX.GT.0) WRITE(IOUT,1003)IGRP
        LROD=.FALSE.
        CALL DSETGR(IPDEV,IMODE,IGRP,NDGR,LROD,IMPX)
        NDEV=NDEV+NDGR-1
        NGRP=NGRP+1
*
      ELSEIF(TEXT.EQ.';') THEN
        GOTO 20
      ELSE
        CALL XABORT('@DSET: WRONG KEYWORD '//TEXT)
      ENDIF
      GOTO 10
   20 IF(IMPX.GT.0) WRITE(IOUT,1004)NGRP,NDEV-1
      IF(IMPX.GT.4) CALL LCMLIB(IPDEV)
      RETURN
*
 1000 FORMAT(/5X,'DSET: ** SETING PARAMETERS FOR ROD #',I3.3)
 1001 FORMAT(/5X,'DSET: ** SETING PARAMETERS FOR LZC #',I2.2)
 1002 FORMAT(/5X,'DSET: ** SETING PARAMETERS FOR ROD-GROUP #',I2.2)
 1003 FORMAT(/5X,'DSET: ** SETING PARAMETERS FOR LZC-GROUP #',I2.2)
 1004 FORMAT(/5X,'--------------------------------------'/
     1        5X,'TOTAL NUMBER OF UPDATED GROUPS  :',I4/
     2        5X,'TOTAL NUMBER OF UPDATED DEVICES :',I4/)
      END
