*DECK TAVG
      SUBROUTINE TAVG(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform computations according to the time-average model.
* 
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
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
* The TAVG: module specification is: 
* FMAP := TAVG: FMAP POWER :: (desctavg) ;
* where
*   FMAP : name of a \emph{fmap} object, that will be updated by the TAVG: 
*     module. The FMAP object must contain the average exit burnups and 
*     refuelling schemes of channels.
*   POWER name of a \emph{power} object containing the channel and bundle 
*     powers, previously computed by the FLPOW: module. The channel and bundle 
*     powers are used by the TAVG: module to compute the normalized axial 
*     power-shape over each channel. 
*   (desctavg) : structure describing the input data to the TAVG: module.
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
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION DFLOT
      LOGICAL LEXIT,LSHAP
      TYPE(C_PTR) IPMAP,IPPOW,JPMAP
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.2)CALL XABORT('@TAVG: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT(' '
     1 //'@TAVG: LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.1)CALL XABORT('@TAVG: MODIFICATION MODE '
     1 //'FOR L_MAP EXPECTED.')
      CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      TEXT=HENTRY(1)
      IF(HSIGN.NE.'L_MAP')CALL XABORT('@TAVG: SIGNATURE '
     1 //' OF '//TEXT//' IS '//HSIGN//'. L_MAP EXPECTED.')
      IPMAP=KENTRY(1)
      IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT(' '
     1 //'@TAVG: LCM OBJECT EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      TEXT=HENTRY(2)
      IF(HSIGN.NE.'L_POWER')CALL XABORT('@TAVG: SIGNATURE '
     1 //' OF '//TEXT//' IS '//HSIGN//'. L_POWER EXPECTED.')
      IF(JENTRY(2).NE.2)CALL XABORT('@TAVG: READ-ONLY MODE '
     1 //'FOR L_POWER EXPECTED.')
      IPPOW=KENTRY(2)
*----
*  READ INPUT DATA
*----
      IMPX=1
      ARP=0.5
      LEXIT=.FALSE.
      LSHAP=.FALSE.
*     PRINTING INDEX
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@TAVG: CHARACTER DATA EXPECTED(1).')
      IF(TEXT.NE.'EDIT')GOTO 10
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1)CALL XABORT('@TAVG: INTEGER DATA EXPECTED.')
      IMPX=MAX(0,NITMA)
*     AX-SHAPE OPTION
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@TAVG: CHARACTER DATA EXPECTED(2).')
   10 IF(TEXT.NE.'AX-SHAPE')GOTO 20
      LSHAP=.TRUE.
*     RELAXATION PARAMETER
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@TAVG: CHARACTER DATA EXPECTED(3).')
   20 IF(TEXT.NE.'RELAX')GOTO 30
      CALL REDGET(ITYP,NITMA,ARP,TEXT,DFLOT)
      IF(ITYP.NE.2)CALL XABORT('@TAVG: REAL DATA EXPECTED.')
      IF(ARP.LE.0.)CALL XABORT('@TAVG: POSITIVE AND NON-ZERO RELAX'
     1 //'ATION PARAMETER EXPECTED.')
*     B-EXIT OPTION
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@TAVG: CHARACTER DATA EXPECTED(4).')
   30 IF(TEXT.NE.'B-EXIT')GOTO 40
      LEXIT=.TRUE.
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@TAVG: CHARACTER DATA EXPECTED(5).')
   40 IF(TEXT.NE.';')CALL XABORT('@TAVG: END TO MODULE ; EXPECTED.')
      IF((.NOT.LSHAP).AND.(.NOT.LEXIT))CALL XABORT('@TAVG: MODULE'
     1 //' OPTION WAS NOT SPECIFIED.')
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NB=ISTATE(1)
      NCH=ISTATE(2)
      NCOMB=ISTATE(3)
*     FUEL-MAP GEOMETRY
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(JPMAP,'STATE-VECTOR',ISTATE)
      NX=ISTATE(3)
      NY=ISTATE(4)
      NZ=ISTATE(5)
*     CHECK EXISTING DATA
      CALL LCMLEN(IPMAP,'BURN-AVG',LENGT,ITYP)
      IF(LENGT.EQ.0)CALL XABORT('@TAVG: MISSING BURNUP DATA IN FUEL'
     1 //'-MAP OBJECT.')
      CALL LCMLEN(IPMAP,'REF-SCHEME',LENGT,ITYP)
      IF(LENGT.EQ.0)CALL XABORT('@TAVG: MISSING REF-SCHEME DATA IN '
     1 //'FUEL-MAP OBJECT.')
      CALL LCMLEN(IPPOW,'POWER-CHAN',LENGT,ITYP)
      IF(LENGT.EQ.0)CALL XABORT('@TAVGCL: MISSING POWER-CHAN DATA I'
     1 //'N L_POWER OBJECT.')
*----
*  PERFORM CALCULATION
*----
      IF(LSHAP)CALL TAVGCL(IPMAP,IPPOW,NCH,NB,NCOMB,NX,NY,NZ,ARP,IMPX)
      IF(LEXIT)CALL TAVGEX(IPMAP,IPPOW,NCH,NCOMB,NX,NY,NZ,IMPX)
      IF(IMPX.GT.2)CALL LCMLIB(IPMAP)
      RETURN
      END
