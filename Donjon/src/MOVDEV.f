*DECK MOVDEV
      SUBROUTINE MOVDEV(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Simulate the time-dependent displacement of individual devices
* and/or of groups of devices in the reactor core.
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
*Comments:
* The MOVDEV: module specification is:
* DEVICE := MOVDEV: DEVICE :: (descmove) ;
* where
*   DEVICE : name of the \emph{device} object that will be modified by the 
*     module. The rods positions are updated according to the current time step
*     of movement.
*   (descmove) : structure describing the input data to the MOVDEV: module.
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
      CHARACTER TEXT*12,HSIGN*12
      INTEGER ISTATE(NSTATE),DGRP
      DOUBLE PRECISION DFLOT
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.GT.1)CALL XABORT('@MOVDEV: ONE PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@MOV'
     1 //'DEV: LCM OBJECT EXPECTED AT LHS.')
      CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_DEVICE')THEN
        TEXT=HENTRY(1)
        CALL XABORT('@MOVDEV: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     1  '. L_DEVICE EXPECTED.')
      ENDIF
      IF(JENTRY(1).NE.1)CALL XABORT('@MOVDEV: MODIFICATION MODE EX'
     1 //'PECTED FOR L_DEVICE.')
*
      CALL LCMGET(KENTRY(1),'STATE-VECTOR',ISTATE)
      IGEO=ISTATE(1)
      IF(IGEO.NE.7)CALL XABORT('@MOVDEV: ONLY 3D-CARTESIAN GEOMETR'
     1 //'Y  ALLOWED.')
      NDEV=ISTATE(2)
      DGRP=ISTATE(3)
      IMODE=ISTATE(6)
      IF(IMODE.EQ.0)CALL XABORT('@MOVDEV: IMODE NOT SET.')
*----
*  RECOVER INFORMATION
*----
      IMPX=1
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@MOVDEV: CHARACTER DATA EXPECTED.')
      IF(TEXT.NE.'EDIT')GOTO 10
*     PRINTING INDEX
      CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1)CALL XABORT('@MOVDEV: INTEGER FOR EDIT EXPECTED.')
*     TIME STEP INCREMENT
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
   10 IF(TEXT.NE.'DELT')CALL XABORT('@MOVDEV: KEYWORD DELT EXPECTED.')
      CALL REDGET(ITYP,NITMA,DELT,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@MOVDEV: REAL FOR DELT EXPECTED.')
      IF(DELT.LE.0.)CALL XABORT('@MOVDEV: VALUE OF DELT SHOULD B'
     1 //'E POSITIVE.')
      ND=0
      NG=0
   20 ND=ND+1
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(TEXT.EQ.'ROD')THEN
*----
*  ROD OPTION
*----
        CALL REDGET(ITYP,ID,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@MOVDEV: INTEGER ROD-ID NUMB'
     1   //'ER EXPECTED.')
        IF((ID.GT.NDEV).OR.(ID.EQ.0))THEN
          WRITE(IOUT,*)'@MOVDEV: READ CURRENT ROD-ID #',ID
          CALL XABORT('@MOVDEV: WRONG ROD-ID NUMBER.')
        ENDIF
        IF(IMPX.GT.0)WRITE(IOUT,1000)ID
        CALL MOVPOS(KENTRY(1),IMODE,ID,DELT,IMPX)
      ELSEIF(TEXT.EQ.'GROUP')THEN
*----
*  GROUP OPTION
*----
        CALL REDGET(ITYP,IGRP,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@MOVDEV: INTEGER GROUP-ID NUM'
     1   //'BER EXPECTED.')
        IF((IGRP.GT.DGRP).OR.(IGRP.EQ.0))THEN
          WRITE(IOUT,*)'@MOVDEV: READ CURRENT GROUP-ID #',IGRP
          CALL XABORT('@MOVDEV: WRONG GROUP-ID NUMBER.')
        ENDIF
        IF(IMPX.GT.0)WRITE(IOUT,1001)IGRP
        CALL MOVGRP(KENTRY(1),IMODE,IGRP,NDGR,DELT,IMPX)
        ND=ND+NDGR-1
        NG=NG+1
*
      ELSEIF(TEXT.EQ.';')THEN
        GOTO 30
      ELSE
        WRITE(IOUT,*)'@MOVDEV: WRONG KEYWORD : ',TEXT
        CALL XABORT('@MOVDEV: KEYWORD ROD OR GROUP EXPECTED.')
      ENDIF
      GOTO 20
   30 IF(IMPX.GT.0)WRITE(IOUT,1002)NG,ND-1
      IF(IMPX.GT.4)CALL LCMLIB(KENTRY(1))
      RETURN
*
 1000 FORMAT(/5X,'MOVING ROD #',I3.3)
 1001 FORMAT(/5X,'MOVING GROUP #',I2.2)
 1002 FORMAT(
     1 /5X,'-------------------------------------'/
     2  5X,'TOTAL NUMBER OF DISPLACED GROUPS : ',I2/
     3  5X,'TOTAL NUMBER OF DISPLACED RODS : ',I3/)
      END
