*DECK CREXSR
      SUBROUTINE CREXSR(IPCPO,LTAB,HCOMPO,NMIXT,IMPX,NISO,IBM,DERIV,UPS,
     1 NAMDIR,NISOR,HISO,ITY,CONC,NBURN,KBURN,IVARTY,IBTYP,BURN0,BURN1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read input data stream for MIX record and recover the information
* from l_compo linked list.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input/output
* IPCPO   pointer to l_compo information.
* HCOMPO  name of l_compo linked list.
* LTAB    flag: =.true. table option; =.false. compo option.
* NMIXT   maximum number of material mixtures.
* IMPX    printing index (=0 for no print).
* NISO    1+maximum number of extracted isotopes.
* IBM     mixture number to be treated.
* NAMDIR  character*12 name of directory in l_compo object.
* DERIV   flag: =.true. derivative of the macrolib is computed with
*         respect to burn1.
* UPS     flag: =.true. no up-scatering cross section will be stored.
* NISOR   1+number of extracted isotopes.
* HISO    hollerith name information for extracted isotopes.
* ITY     =0: do not process the isotope; =1: use number density
*         stored in conc(i); =2: use number density stored in compo.
* CONC    user defined number density.
* NBURN   number of burnup steps in compo linked list.
* BURN0   user defined initial burnup.
* BURN1   user defined final burnup:  if burn0=burn1 => a simple 
*         interpolation is performed; if burn0<burn1 => a time-average
*         calculation is performed.
* KBURN   =0: no burnup parameters; =1: use mw day/tonne of initial
*         heavy elements.
* IVARTY  index of the exit burnup used to compute derivatives. Set to
*         zero to avoid taking the derivative.
* IBTYP   type of interpolation: =1 time-average; =2 instantaneous;
*         derivative with respect to a single exit burnup.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER NMIXT,IMPX,NISO,IBM,NISOR,NBURN,KBURN,IVARTY,IBTYP,
     1        HISO(3*NISO),ITY(NISO)
      LOGICAL DERIV,UPS,LTAB
      CHARACTER NAMDIR*12,HCOMPO*12
      REAL CONC(NISO),BURN0,BURN1
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6)
      TYPE(C_PTR) JPCPO
      INTEGER IPAR(NSTATE)
      CHARACTER TEXT12*12,TEXT*12,CGRPNM*12
      DOUBLE PRECISION DFLOT
*
      KBURN=0
      CALL XDISET(ITY,NISO,0)
      ITY(1)=2
*----
*  RECOVER INFORMATION
*----
      IBM=0
      TEXT12='MIX'
   10 IF(TEXT12.EQ.'MIX')THEN
        IVARTY=0
        IBTYP=0
        IF(IBM.NE.0)CALL XABORT('@CREXSR: MIX ALREADY SELECTED.')
        CALL REDGET(ITYP,IBM,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@CREXSR: INTEGER DATA EXPECTED.')
        IF(IBM.GT.NMIXT)CALL XABORT('@CREXSR: INVALID MIX INDEX.')
        CALL REDGET(ITYP,NITMA,FLOT,NAMDIR,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@CREXSR: CHARACTER DATA EXPECTED.')
        IF(IMPX.GT.0)WRITE(6,'(/27H CREXSR: ACCESS DIRECTORY '',A12,
     1  17H'' IN COMPO FILE '',A12,2H''.)') NAMDIR,HCOMPO
        JPCPO=LCMGID(IPCPO,NAMDIR)
        CALL LCMGET(JPCPO,'PARAM',IPAR)
        NISOR=IPAR(2)
        NBURN=IPAR(4)
        IF(NISOR.GT.1)CALL LCMGET(JPCPO,'ISOTOPESNAME',HISO)
      ELSEIF(TEXT12.EQ.'I-BURNUP')THEN
        IF(LTAB )CALL XABORT('@CREXSR: INVALID OPTION I-BURNUP WITH'
     1   //' FUEL MAP OBJECT.')
        IF(NBURN.LE.1)CALL XABORT('@CREXSR: NO BURNUP INFORMATION.')
        KBURN=1
        CALL REDGET(ITYP,NITMA,BURN0,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@CREXSR: REAL DATA EXPECTED(1).')
        BURN1=BURN0
      ELSEIF(TEXT12.EQ.'T-BURNUP')THEN
        IF(LTAB )CALL XABORT('@CREXSR: INVALID OPTION T-BURNUP WITH'
     1   //' FUEL MAP OBJECT.')
        IF(NBURN.LE.1)CALL XABORT('@CREXSR: NO BURNUP INFORMATION.')
        KBURN=1
        CALL REDGET(ITYP,NITMA,BURN0,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@CREXSR: REAL DATA EXPECTED(2).')
        CALL REDGET(ITYP,NITMA,BURN1,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@CREXSR: REAL DATA EXPECTED(3).')
        IF(BURN1.LE.BURN0)CALL XABORT('@CREXSR: INVALID BURN1.')
      ELSEIF(TEXT12.EQ.'MICRO')THEN
        IF(NISO.LE.1)CALL XABORT('NO EXTRACTED ISOTOPES IN L_COMPO.')
   20   CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@CREXSR: CHARACTER DATA EXPECTED.')
        IF(TEXT.EQ.'ALL')THEN
          DO 30 I=2,NISO
          ITY(I)=2
   30     CONTINUE
        ELSEIF(TEXT.EQ.'ENDMIX')THEN
          TEXT12=TEXT
          GOTO 10
        ELSEIF(TEXT.EQ.'UPS')THEN
          TEXT12=TEXT
          GOTO 10
        ELSE
          DO 50 I=1,NISO
          WRITE(CGRPNM,'(3A4)') (HISO(3*(I-1)+J),J=1,3)
          IF(CGRPNM.EQ.TEXT)THEN
            CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
            IF(ITYP.EQ.2)THEN
              CONC(I)=FLOT
              ITY(I)=1
            ELSEIF((ITYP.EQ.3).AND.(TEXT.EQ.'*'))THEN
              ITY(I)=2
            ELSE
              CALL XABORT('@CREXSR: REAL NUMBER OR * EXPECTED.')
            ENDIF
            GOTO 20
          ENDIF
   50     CONTINUE
          CALL XABORT('@CREXSR: UNABLE TO MATCH ISOTOPE'//TEXT//'.')
        ENDIF
      ELSEIF(TEXT12.EQ.'DERIV')THEN
        DERIV=.TRUE.
      ELSEIF(TEXT12.EQ.'TIMAV-BURN')THEN
        IBTYP=1
      ELSEIF(TEXT12.EQ.'INST-BURN')THEN
        IBTYP=2
      ELSEIF(TEXT12.EQ.'AVG-EX-BURN') THEN
        IBTYP=3
        CALL REDGET(INDIC,IVARTY,FLOTT,TEXT12,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('CREXSR: INTEGER DATA EXPECTED'
     1  //'(AVG-EX-BURN).')
      ELSEIF(TEXT12.EQ.'UPS')THEN
        UPS=.TRUE.
      ELSEIF(TEXT12.EQ.'ENDMIX')THEN
        IF(LTAB)THEN
          IF(NBURN.LE.0)CALL XABORT('@CREXSR: NO BURNUP INFORMATION '
     1     //'FOR THIS MIXTURE.')
        ELSE
          IF((KBURN.EQ.0).AND.(NBURN.GT.1))CALL XABORT('@CREXSR: BUR'
     1     //'NUP INTEGRATION OPTION REQUIRED.')
        ENDIF
        RETURN
      ELSE
        WRITE(IOUT,'(A40)')'@CREXSR: MIX SHOULD FINISH WITH ENDMIX.'
        CALL XABORT('@CREXSR: WRONG KEYWORD '//TEXT12//'.')
      ENDIF
      CALL REDGET(ITYP,NITMA,FLOT,TEXT12,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@CREXSR: CHARACTER DATA EXPECTED.')
      GOTO 10
      END
