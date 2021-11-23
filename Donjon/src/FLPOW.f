*DECK FLPOW
      SUBROUTINE FLPOW(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute and print power and flux distributions over the reactor core.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Update(s):
*  M. Guyot 15/07/10 : Creation of L_FLUX object to be used by
*  module DETECT:, 
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
* The FLPOW: module specifications are:
* Option 1
* POWER [ NRMFLUX ] [ FMAP ] := FLPOW:  [ POWOLD ] FMAP
*  { FLUX | KINET } TRACK MATEX [ MACRO ] :: (descflpow} ;
* Option 2
* POWER := FLPOW:  [ POWOLD ] { FLUX | KINET }  TRACK MACRO :: (descflpow) ;
* where
*   POWER   : name of the \emph{power} object that will be created by the 
*     module. It will contain the information related to the reactor fluxes 
*     and powers.
*   NRMFLUX : name of the \emph{flux} object, in creation mode. According to 
*     the chosen option, this object contains either the fluxes normalized to 
*     the given total reactor power or the fluxes per bundle. Is it useful if 
*     you want to compute the detectors readings with the DETECT: module.
*   POWOLD  : name of the read-only \emph{power} object. It must contain the 
*     previously computed flux normalization factor, which corresponds to the 
*     reactor nominal or equilibrium conditions.
*   FMAP    : name of the \emph{fmap} object containing the fuel lattice 
*     specification. When FMAP is specified on the RHS, the fluxes and powers 
*     calculations are performed over the fuel lattice as well as over the 
*     whole reactor geometry. If FMAP is specified on the LHS, its records 
*     'BUND-PW' and 'FLUX-AV' will be set according to the information present 
*     in POWER.
*   FLUX    : name of the \emph{flux} object, previously created by the 
*     FLUD: module. The numerical flux solution contained in FLUX is 
*     recovered and all flux are normalized to the given total reactor power.
*   KINET   : name of the \emph{kinet} object, previously created by the 
*     KINSOL: module. The numerical flux solution contained in KINET is 
*     recovered.
*   TRACK   : name of the \emph{track} object, created by the TRIVAT: module. 
*     The information stored in TRACK is recovered and used for the average 
*     flux calculation.
*   MATEX   : name of the \emph{matex} object, containing the reactor material 
*     index and the h-factors that will be recovered and used for the power 
*     calculation.
*   MACRO name of the \emph{macrolib} object, containing the h-factors that 
*     will be recovered and used for the power calculation.
*   (descflpow) : structure describing the input data to the FLPOW: module .
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
      CHARACTER HSIGN*12,TEXT*12
      LOGICAL LNEW,LMAP,LFLX,LRAT,LPOW,LFSTH,LFLU,LNRM,LBUN
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) IPPOW,IPFLX,IPKIN,IPTRK,IPMTX,IPMAP,IPMAC,IPNFX
*----
*  PARAMETER VALIDATION
*----
      LFLU=.FALSE.
      IF(NENTRY.LT.4)CALL XABORT('@FLPOW: PARAMETER EXPECTED.')
      TEXT=HENTRY(1)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@FLPOW'
     1 //': LCM OBJECT FOR L_POWER EXPECTED AT LHS ('//TEXT//').')
      IF(JENTRY(1).NE.0)CALL XABORT('@FLPOW: CREATE MODE FOR L_POW'
     1 //'ER EXPECTED AT LHS ('//TEXT//').')
      IPPOW=KENTRY(1)
      IF(JENTRY(2).EQ.0)THEN
        LFLU=.TRUE.
        IPNFX=KENTRY(2)
      ENDIF
      IPFLX=C_NULL_PTR
      IPKIN=C_NULL_PTR
      IPTRK=C_NULL_PTR
      IPMTX=C_NULL_PTR
      IPMAP=C_NULL_PTR
      IPMAC=C_NULL_PTR
      LNEW=.FALSE.
      JMOD=0
      IF(LFLU)THEN
        NRHS=3
      ELSE
        NRHS=2
        IPNFX=C_NULL_PTR
      ENDIF
      DO 10 IEN=NRHS,NENTRY
      IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))CALL XABORT('@F'
     1 //'LPOW: LCM OBJECT EXPECTED AT THE RHS.')
      CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_POWER')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_POWER OBJECT AT RHS.')
        IF(LNEW)CALL XABORT('@FLPOW: L_POWER ALREADY DEFINED AT RHS.')
        CALL LCMEQU(KENTRY(IEN),IPPOW)
        LNEW=.TRUE.
      ELSEIF(HSIGN.EQ.'L_MATEX')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_MATEX OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPMTX))THEN
          IPMTX=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_MATEX ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_FLUX')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_FLUX OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPFLX))THEN
          IPFLX=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_FLUX ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_KINET')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_KINET OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPKIN))THEN
          IPKIN=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_KINET ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_TRACK')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_TRACK OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPTRK))THEN
          IPTRK=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_TRACK ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_MACROLIB')THEN
        IF(JENTRY(IEN).NE.2)CALL XABORT('@FLPOW: READ-ONLY MODE EXPE'
     1  //'CTED FOR THE L_MACROLIB OBJECT AT RHS.')
        IF(.NOT.C_ASSOCIATED(IPMAC))THEN
          IPMAC=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_MACROLIB ALREADY DEFINED.')
        ENDIF
      ELSEIF(HSIGN.EQ.'L_MAP')THEN
        IF(JENTRY(IEN).EQ.1) JMOD=1
        IF(.NOT.C_ASSOCIATED(IPMAP))THEN
          IPMAP=KENTRY(IEN)
        ELSE
          CALL XABORT('@FLPOW: L_MAP ALREADY DEFINED.')
        ENDIF
      ENDIF
   10 CONTINUE
      IF((.NOT.C_ASSOCIATED(IPFLX)).AND.(.NOT.C_ASSOCIATED(IPKIN))) THEN
         CALL XABORT('@FLPOW: MISSING L_FLUX OR L_KINET OBJECT.')
      ELSE IF((C_ASSOCIATED(IPFLX)).AND.(C_ASSOCIATED(IPKIN))) THEN
         CALL XABORT('@FLPOW: L_FLUX AND L_KINET OBJECTS BOTH DEFINED.')
      ELSE IF(.NOT.C_ASSOCIATED(IPTRK)) THEN
         CALL XABORT('@FLPOW: MISSING L_TRACK OBJECT.')
      ELSE IF((C_ASSOCIATED(IPMAP)).AND.(.NOT.C_ASSOCIATED(IPMTX))) THEN
         CALL XABORT('@FLPOW: MISSING L_MATEX OBJECT.')
      ELSE IF((.NOT.C_ASSOCIATED(IPMAP)).AND.(C_ASSOCIATED(IPMTX))) THEN
         CALL XABORT('@FLPOW: MISSING L_MAP OBJECT.')
      ELSE IF((.NOT.C_ASSOCIATED(IPMTX)).AND.
     1        (.NOT.C_ASSOCIATED(IPMAC))) THEN
         CALL XABORT('@FLPOW: MISSING L_MATEX OR L_MACROLIB OBJECT.')
      ELSE IF((.NOT.C_ASSOCIATED(IPMAP)).AND.
     1        (.NOT.C_ASSOCIATED(IPMAC))) THEN
         CALL XABORT('@FLPOW: MISSING L_MAP OR L_MACROLIB OBJECT.')
      ENDIF
*----
*  READ KEYWORD
*----
      IMPX=1
      PTOT=0.0
      LFSTH=.FALSE.
      LNRM=.FALSE.
      LBUN=.FALSE.
      FSTH=0.0
      LFLX=.FALSE.
      LPOW=.FALSE.
      LMAP=.FALSE.
      LRAT=.FALSE.
   20 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.EQ.10) GO TO 40
   30 IF(ITYP.NE.3)CALL XABORT('@FLPOW: CHARACTER DATA EXPECTED.')
      IF(TEXT.EQ.'EDIT') THEN
*       PRINTING INDEX
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@FLPOW: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT.EQ.'P-NEW') THEN
        IF(.NOT.LNEW)CALL XABORT('@FLPOW: MISSING READ-ONLY L_POWER'
     1  //' OBJECT AT RHS.')
      ELSE IF(TEXT.EQ.'PTOT') THEN
        IF(LNEW)CALL XABORT('@FLPOW: ONLY ONE L_POWER OBJECT IN CRE'
     1  //'ATE MODE EXPECTED WITH PTOT OPTION.')
        CALL REDGET(ITYP,NITMA,PTOT,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@FLPOW: REAL FOR PTOT EXPECTED.')
        IF(PTOT.LE.0.)CALL XABORT('@FLPOW: INVALID VALUE PTOT < 0.')
      ELSE IF(TEXT.EQ.'FSTH') THEN
        CALL REDGET(ITYP,NITMA,FSTH,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@FLPOW: REAL DATA EXPECTED FOR FSTH.')
        IF((FSTH.GT.1.0).OR.(FSTH.LE.0.0)) CALL XABORT('@FLPOW: FSTH '
     1  //'SHOULD BE BETWEEN 0.0 AND 1.0.')
        LFSTH=.TRUE.
      ELSE IF(TEXT.EQ.'NORM') THEN
        LNRM=.TRUE.
      ELSE IF(TEXT.EQ.'BUND') THEN
        IF(.NOT.C_ASSOCIATED(IPMAP)) CALL XABORT('@FLPOW: NO RHS FUELM'
     1  //'AP DEFINED.')
        LBUN=.TRUE.
      ELSE IF(TEXT.EQ.'PRINT') THEN
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.EQ.'MAP')THEN
          IF(.NOT.C_ASSOCIATED(IPMAP))CALL XABORT('@FLPOW: INVALID KEY'
     1    //'WORD MAP. MISSING L_MAP OBJECT FOR PRINT.')
          LMAP=.TRUE.
        ELSEIF(TEXT.EQ.'ALL')THEN
          LFLX=.TRUE.
          LPOW=.TRUE.
          IF(C_ASSOCIATED(IPMAP))LMAP=.TRUE.
          LRAT=.TRUE.
        ELSEIF(TEXT.EQ.'DISTR')THEN
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.NE.3)CALL XABORT('@FLPOW: CHARACTER DATA EXPECTED AF'
     1    //'TER DISTR.')
          IF(TEXT.EQ.'FLUX')THEN
            IF(LFLX)CALL XABORT('@FLPOW: KEYWORD FLUX ALREADY READ.')
            LFLX=.TRUE.
          ELSEIF(TEXT.EQ.'POWER')THEN
            IF(LPOW)CALL XABORT('@FLPOW: KEYWORD POWER ALREADY READ.')
            LPOW=.TRUE.
          ELSEIF(TEXT.EQ.'RATIO')THEN
            IF(LRAT)CALL XABORT('@FLPOW: KEYWORD RATIO ALREADY READ.')
            LRAT=.TRUE.
          ELSE
            GO TO 30
          ENDIF
        ELSE
          CALL XABORT('@FLPOW: KEYWORD MAP/DISTR/ALL EXPECTED.')
        ENDIF
      ELSE IF(TEXT.EQ.'INIT') THEN
        IF(JENTRY(IEN).EQ.1) JMOD=2
      ELSE IF(TEXT.EQ.';') THEN
        GO TO 40
      ELSE
        CALL XABORT('@FLPOW: INVALID KEYWORD '//TEXT//'.')
      ENDIF
      GO TO 20
*----
*  CHECK CONSISTENCY
*----
   40 IF(LMAP) THEN
        IF(.NOT.C_ASSOCIATED(IPMAP)) THEN
           CALL XABORT('@FLPOW: MISSING L_MAP OBJECT.')
        ELSE IF(.NOT.C_ASSOCIATED(IPMTX)) THEN
           CALL XABORT('@FLPOW: MISSING L_MATEX OBJECT.')
        ELSE IF(.NOT.C_ASSOCIATED(IPMAC)) THEN
           CALL XABORT('@FLPOW: MISSING L_MACROLIB OBJECT.')
        ENDIF
      ENDIF
*----
*  PERFORM CALCULATION
*----
      CALL FLPDRV(IPPOW,IPNFX,IPFLX,IPKIN,IPTRK,IPMTX,IPMAP,IPMAC,PTOT,
     1 LNEW,LMAP,JMOD,LFLX,LPOW,LRAT,IMPX,FSTH,LFSTH,LFLU,LBUN,LNRM)
      RETURN
      END
