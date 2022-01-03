*DECK ROD
      SUBROUTINE ROD(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Control rod management module based on SAPHYB or MULTICOMPO
* interpolation.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Tixier
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
* The ROD: module specifications are:
* FLMAP := ROD: FLMAP :: (descrod1) ;
* where
*   FLMAP  :name of the \emph{MAP} object that will contain the 3-D rod file. 
*     The FLMAP has to be modified for the module and must appear on both LHS 
*     and RHS.
*   (descrod1) : structure describing the main input data to the ROD: module. 
*     Note that this input data is mandatory and must be specified.
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
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),NPARAM,NREB,INSM
      INTEGER MAXMIX,NGRP,RODSIZE,NASS,RODINFO,NCALL
      REAL INI,INSS
      LOGICAL :: EXISTENCE=.FALSE.
      CHARACTER HSIGN*12,TEXT*40,PAR1*12,PNAME*12
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) IPMAP,JPMAP,KPMAP,MPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: INTROD,HMIX,RMIX
      CHARACTER(LEN=3), ALLOCATABLE, DIMENSION(:) :: RNAME
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INS,NUMMIX
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.NE.1)CALL XABORT('@ROD: 1 OBJECT EXPECTED.')
      IPMAP=KENTRY(1)
      IF(IENTRY(1).NE.1)CALL XABORT('@ROD:'
     > //' LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.1)CALL XABORT('@ROD: FLMAP MUST BE IN'
     > //' MODIFICATION MODE AND NOT IN CREATION MODE.')
      CALL LCMGTC(IPMAP,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MAP')THEN
        TEXT=HENTRY(1)
        CALL XABORT('@ROD: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     >  '. L_MAP EXPECTED.')
      ENDIF
*----
*  RECOVER L_MAP STATE-VECTOR
*----
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NB=ISTATE(1)
      NCH=ISTATE(2)
      NPARAM=ISTATE(8)
*----
*  READ INPUT DATA
*----
      NCALL=0
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@ROD: CHARACTER DATA EXPECTED.')
* Read printing index
      IF(TEXT.EQ.'EDIT') THEN
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@ROD: INTEGER FOR EDIT EXPECTED.')
* Name of the parameter record that is to be created
      ELSE IF(TEXT.EQ.'PARA') THEN
        NCALL=1
        CALL REDGET(ITYP,NITMA,FLOT,PAR1,DFLOT)
        IF(ITYP.NE.3) CALL XABORT('@ROD: CHARACTER'
     1   //' DATA FOR PARAMETER NAME EXPECTED.')
* Checking of the record existence
        JPMAP=LCMGID(IPMAP,'PARAM')
        EXISTENCE=.FALSE.
        DO IPAR=1,NPARAM,1
          KPMAP=LCMGIL(JPMAP,IPAR)
          CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
          IF(PNAME.EQ.PAR1) THEN
            EXISTENCE=.TRUE.
            EXIT  
          ENDIF
        ENDDO
        IF(.NOT.EXISTENCE) THEN
* If PARAM doesn't exist, it is created
          NPARAM=NPARAM+1
          JPMAP=LCMLID(IPMAP,'PARAM',NPARAM)
          KPMAP=LCMDIL(JPMAP,NPARAM)
          CALL LCMPTC(KPMAP,'P-NAME',12,1,PAR1)
          CALL LCMPTC(KPMAP,'PARKEY',12,1,PAR1)
          IPTYP=2 
          CALL LCMPUT(KPMAP,'P-TYPE',1,1,IPTYP)
          RODINFO=4
          MPMAP=LCMDID(IPMAP,'ROD-INFO')
          CALL REDGET(ITYP,NITMA,INI,TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@ROD: REAL DATA FOR'
     1     //' STEP EXPECTED.')
          CALL LCMPUT(MPMAP,'ROD-INIT',1,2,INI)
        ENDIF
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
* Check if LINS is defined
        IF(TEXT.NE.'LINS')CALL XABORT('@ROD: KEYWORD'
     1   //' LINS EXPECTED.')
        CALL REDGET(ITYP,INSM,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@ROD: INTEGER DATA FOR'
     1   //' LINS EXPECTED.')
        IF(INSM.LT.0)CALL XABORT('@ROD: LINS MUST BE POSITIVE.')
        CALL LCMPUT(MPMAP,'INS-MAX',1,1,INSM)
* Check if STEP is defined
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(TEXT.NE.'STEP')CALL XABORT('@ROD: KEYWORD'
     1   //' STEP EXPECTED.')
        CALL REDGET(ITYP,NITMA,INSS,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@ROD: REAL DATA FOR'
     1   //' STEP EXPECTED.')
        IF(INSS.LT.0.0)CALL XABORT('@ROD: STEP MUST BE POSITIVE.')
        CALL LCMPUT(MPMAP,'STEP-CM',1,2,INSS)
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
* Check if NRFB is defined
        IF(TEXT.NE.'NRFB')CALL XABORT('@ROD: KEYWORD '
     1   //'NRFB EXPECTED.')
        CALL REDGET(ITYP,NREB,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@ROD: INTEGER DATA FOR'
     1   //' NRFB EXPECTED.')
        IF(NREB.LT.0)CALL XABORT('@ROD: NRFB MUST BE POSITIVE.')
        CALL LCMPUT(MPMAP,'REFL-BOTTOM',1,1,NREB)
* Definition of rod groups
      ELSE IF(TEXT.EQ.'RGRP') THEN
        JPMAP=LCMGID(IPMAP,'PARAM')
        KPMAP=LCMGIL(JPMAP,NPARAM)
        IF(NCALL.EQ.1) THEN
* Creation of records with the number of rod groups and the maximum of
* rod zones
        CALL REDGET(ITYP,NGRP,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@ROD: INTEGER'
     1   //' DATA FOR GROUP NUMBER EXPECTED.')
        CALL LCMPUT(MPMAP,'NB-GROUP',1,1,NGRP)
        CALL REDGET(ITYP,MAXMIX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@ROD: INTEGER'
     1   //' DATA FOR MAXIMUM MIX NUMBER EXPECTED.')
        CALL LCMPUT(MPMAP,'MAX-MIX',1,1,MAXMIX)
        ALLOCATE(RNAME(NGRP),INS(NGRP),NUMMIX(NGRP))
        ALLOCATE(HMIX(NGRP*MAXMIX),RMIX(NGRP*MAXMIX))
        RODSIZE=NCH*NB
        ALLOCATE(INTROD(RODSIZE))
        CALL XDRSET(HMIX,NGRP*MAXMIX,0.0)
        CALL XDRSET(RMIX,NGRP*MAXMIX,-999.0)
        CALL XDISET(INS,NGRP,-1)
        CALL RODTYP(IPMAP,NGRP,MAXMIX,RNAME,INS,HMIX,RMIX,NUMMIX)
        ELSE
* Recovering rod parameters in order to modify only groups defined
* in the second call of the module.
        MPMAP=LCMGID(IPMAP,'ROD-INFO')
        CALL LCMGET(MPMAP,'NB-GROUP',NGRP)
        CALL LCMGET(MPMAP,'MAX-MIX',MAXMIX)
        ALLOCATE(RNAME(NGRP),INS(NGRP),NUMMIX(NGRP))
        ALLOCATE(HMIX(NGRP*MAXMIX),RMIX(NGRP*MAXMIX))
        RODSIZE=NCH*NB
        ALLOCATE(INTROD(RODSIZE))
* Store rod insertion modification in the fuel map
        CALL RODMOV(IPMAP,NGRP,RNAME,INS)
        CALL XDRSET(INTROD,RODSIZE,INI)
        CALL RODMOD(IPMAP,NGRP,MAXMIX,NCH,NB,RNAME,INS,INSS,HMIX,
     >  RMIX,NREB,RODSIZE,INTROD,INI,NUMMIX,NCALL)        
        ENDIF
* Definition of the rod map
      ELSE IF(TEXT.EQ.'RMAP') THEN
        CALL XDRSET(INTROD,RODSIZE,INI)
        CALL REDGET(ITYP,NASS,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1) CALL XABORT('@ROD: INTEGER'
     1   //' DATA FOR ASSEMBLY NUMBER EXPECTED.')
        IF(NASS.NE.NCH) CALL XABORT('@ROD: NUMBER OF ASSEMBLIES'
     1   //' MUST BE EQUAL TO NCH.')
        CALL RODMOD(IPMAP,NGRP,MAXMIX,NCH,NB,RNAME,INS,INSS,HMIX,
     >  RMIX,NREB,RODSIZE,INTROD,INI,NUMMIX,NCALL)
      ELSE IF(TEXT.EQ.';') THEN
*----
*  SAVE ROD INSERTION INFORMATION ON LCM OBJECT L_MAP
*----
        CALL LCMPUT(KPMAP,'P-VALUE',RODSIZE,2,INTROD)
        ISTATE(8)=NPARAM
        CALL LCMPUT(IPMAP,'STATE-VECTOR',NSTATE,1,ISTATE)
        DEALLOCATE(RNAME,INS,HMIX,RMIX,INTROD)
        GO TO 20
      ELSE
        CALL XABORT('@ROD: INVALID KEYWORD: '//TEXT//'.')
      ENDIF
      GO TO 10
*
   20 IF(IMPX.GT.2) CALL LCMLIB(IPMAP)
      RETURN
      END
