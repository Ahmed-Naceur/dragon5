*DECK MCC
      SUBROUTINE MCC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Fuel map modification module.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M. Cordiez
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
* The MCC: module specifications are:
* [FLMAP1] := MCC: FLMAP1 [FLMAP2] :: (descmcc1) ; 
* where
*   FLMAP1 : name of the \emph{MAP} object that will contain the updated 
*     fuel-lattice information. If FLMAP1 appears on both LHS and RHS, it will 
*     be updated; if it only appears on RHS, it will only be read to display 
*     its contents.
*   FLMAP2 : name of the \emph{MAP} object that contains information to be 
*     recovered to update FLMAP1. If FLMAP2 exists, data to update FLMAP1 will 
*     be taken in it. If not, data to update FLMAP1 will be taken in FLMAP1.
*   (descmcc1) : structure describing the main input data to the MCC: module. 
*     Note that this input data is mandatory and must be specified either if 
*     FLMAP1 is updated or only read.
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
      INTEGER ISTATE(NSTATE),NPARAM
      INTEGER PTYPETCOOL,PTYPEDCOOL,VALSIZE
      REAL TSAT
      CHARACTER HSIGN*12,TEXT*40,REC1*40,REC2*40,PNAME*12
      DOUBLE PRECISION DFLOT
      LOGICAL :: EXISTENCE=.FALSE.,EXISTENCE2=.FALSE.
      LOGICAL :: PRESTCOOL=.FALSE.,PRESDCOOL=.FALSE.
      TYPE(C_PTR) IPMAP,JPMAP,KPMAP,IPMAP2
      REAL, ALLOCATABLE, DIMENSION(:) :: VALTCOOL,VALDCOOL
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LT.1)CALL XABORT('@MCC: MINIMUM OF 1 OBJECT EXPECTED.')
      IPMAP=KENTRY(1)
      IF(NENTRY.EQ.2) IPMAP2=KENTRY(2)
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2))CALL XABORT('@MCC:'
     > //' LCM OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.1)CALL XABORT('@MCC: FLMAP1 MUST BE IN'
     > //' MODIFICATION MODE AND NOT IN CREATION MODE.')
      CALL LCMGTC(IPMAP,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MAP')THEN
        TEXT=HENTRY(2)
        CALL XABORT('@MCC: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     >  '. L_MAP EXPECTED.')
      ENDIF
      IF(NENTRY.EQ.2) THEN
        IPMAP2=KENTRY(2)
        IF((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2))CALL XABORT('@MCC:'
     > //' LCM OBJECT EXPECTED FOR FLMAP2.')
        IF(JENTRY(2).NE.2)CALL XABORT('@MCC: FLMAP2 MUST BE IN READ-'
     >    //'ONLY MODE AND NOT IN CREATION MODE.')
        CALL LCMGTC(IPMAP,'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_MAP')THEN
          TEXT=HENTRY(2)
          CALL XABORT('@MCC: SIGNATURE OF '//TEXT//' IS '//HSIGN//
     >    '. L_MAP EXPECTED.')
        ENDIF
      ENDIF
*----
*  RECOVER L_MAP STATE-VECTOR
*----
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NB=ISTATE(1)
      NCH=ISTATE(2)
      NPARAM=ISTATE(8)
      IMPX=1
*----
*  READ INPUT DATA
*----
   10 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@MCC: CHARACTER DATA EXPECTED.')
* Read printing index
      IF(TEXT.EQ.'EDIT') THEN
        CALL REDGET(ITYP,IMPX,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@MCC: INTEGER FOR EDIT EXPECTED.')
* Name of the record that is to be modified
      ELSE IF(TEXT.EQ.'REC') THEN
        CALL REDGET(ITYP,NITMA,FLOT,REC1,DFLOT)
        IF(ITYP.NE.3) CALL XABORT('@MCC: STRING PARAMETER EXPECTED '
     >               //'FOR THE NAME OF THE RECORD THAT IS '
     >               //'TO BE MODIFIED.')
* Checking of the record existence
        JPMAP=LCMGID(IPMAP,'PARAM')
        EXISTENCE=.FALSE.
        DO IPAR=1,NPARAM,1
          KPMAP=LCMGIL(JPMAP,IPAR)
          CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
          IF(PNAME.EQ.REC1) THEN
            EXISTENCE=.TRUE.
            EXIT  
          ENDIF
        ENDDO
        IF(.NOT.EXISTENCE) CALL XABORT('@MCC: LOCAL PARAMETER: '
     >                 //REC1//' DOES NOT EXIST IN THE FUEL MAP.')
*
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
*********
* Case of a uniform edition
*********
        IF(TEXT.EQ.'UNI') THEN
          CALL REDGET(ITYP,NITMA,VAL1,TEXT,DFLOT)
          IF(ITYP.NE.2) CALL XABORT('@MCC: REAL VALUE EXPECTED FOR '
     >                            //'value1.')
*         Fuel map modification: every value set to VAL1
          CALL MCCMOD(IMPX,IPMAP,NPARAM,NCH,NB,REC1,VAL1,0)
        ELSEIF(TEXT.EQ.'ADD') THEN
          CALL REDGET(ITYP,NITMA,VAL2,TEXT,DFLOT)
          IF(ITYP.NE.2) CALL XABORT('@MCC: REAL VALUE EXPECTED FOR '
     >                            //'value2.')
*         Fuel map modification: every value incremented of VAL2
          CALL MCCMOD(IMPX,IPMAP,NPARAM,NCH,NB,REC1,VAL2,1)
*********
* Case of a copy from a different directory or fuel map
*********
* Same fuel map
        ELSEIF(TEXT.EQ.'SAME') THEN
          CALL REDGET(ITYP,NITMA,FLOT,REC2,DFLOT)
          IF(ITYP.NE.3) CALL XABORT('@MCC: STRING PARAMETER EXPECTED '
     >                            //'FOR THE NAME OF THE RECORD rec2 ')
          IF((REC1.EQ.REC2).AND.(IMPX.GT.0)) WRITE(6,'(A)') 'WARNING: '
     >           //'rec1 AND rec2 ARE IDENTICAL! THIS CALL IS USELESS.'
* Checking of the record existence
          JPMAP=LCMGID(IPMAP,'PARAM')
          EXISTENCE2=.FALSE.
          DO IPAR=1,NPARAM,1
            KPMAP=LCMGIL(JPMAP,IPAR)
            CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
            IF(PNAME.EQ.REC2) THEN
              EXISTENCE2=.TRUE.
              EXIT  
            ENDIF
          ENDDO
          IF(.NOT.EXISTENCE) CALL XABORT('@MCC: LOCAL PARAMETER: '
     >                 //REC1//' DOES NOT EXIST IN THE FUEL MAP.')
          CALL MCCCPY(IMPX,IPMAP,IPMAP,NPARAM,NCH,NB,REC1,REC2)
*
* Different fuel map
        ELSEIF(TEXT.EQ.'READ') THEN
          CALL REDGET(ITYP,NITMA,FLOT,REC2,DFLOT)
          IF(ITYP.NE.3) CALL XABORT('@MCC: STRING PARAMETER EXPECTED '
     >                            //'FOR THE NAME OF THE RECORD rec2 ')
* Checking of the record existence
          JPMAP=LCMGID(IPMAP2,'PARAM')
          EXISTENCE2=.FALSE.
          DO IPAR=1,NPARAM,1
            KPMAP=LCMGIL(JPMAP,IPAR)
            CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
            IF(PNAME.EQ.REC2) THEN
              EXISTENCE2=.TRUE.
              EXIT  
            ENDIF
          ENDDO
          IF(.NOT.EXISTENCE) CALL XABORT('@MCC: LOCAL PARAMETER: '
     >                 //REC1//' DOES NOT EXIST IN THE FUEL MAP.')
          CALL MCCCPY(IMPX,IPMAP,IPMAP2,NPARAM,NCH,NB,REC1,REC2)
        ELSE
          CALL XABORT('@MCC: WRONG KEYWORD.')
        ENDIF
*********
* Calculation of D-COOL from T-COOL
*********
      ELSE IF(TEXT.EQ.'TTD') THEN
        CALL REDGET(ITYP,NITMA,PINLET,TEXT,DFLOT)
        IF(ITYP.NE.2) CALL XABORT('@MCC: REAL PARAMETER EXPECTED '
     >                          //'FOR THE CORE PRESSURE.')
* Checking of the existence of the T-COOL and D-COOL directories
* Recovery of T-COOL data
        JPMAP=LCMGID(IPMAP,'PARAM')
        DO IPAR=1,NPARAM,1
          KPMAP=LCMGIL(JPMAP,IPAR)
          CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
          IF(PNAME.EQ.'T-COOL') THEN
            PRESTCOOL=.TRUE.
            CALL LCMGET(KPMAP,'P-TYPE',PTYPETCOOL)
            IF(PTYPETCOOL.EQ.1) THEN
              VALSIZE=1
              ALLOCATE(VALTCOOL(VALSIZE))
              CALL LCMGET(KPMAP,'P-VALUE',VALTCOOL)
            ELSE
              VALSIZE=NCH*NB
              ALLOCATE(VALTCOOL(VALSIZE))
              CALL LCMGET(KPMAP,'P-VALUE',VALTCOOL)
            ENDIF
          ENDIF
          IF(PNAME.EQ.'D-COOL') THEN
            PRESDCOOL=.TRUE.
            CALL LCMGET(KPMAP,'P-TYPE',PTYPEDCOOL)
            IF(PTYPEDCOOL.EQ.1) THEN
              VALSIZE=1
              ALLOCATE(VALDCOOL(VALSIZE))
            ELSE
              VALSIZE=NCH*NB
              ALLOCATE(VALDCOOL(VALSIZE))
            ENDIF
          ENDIF
        ENDDO
        IF(.NOT.PRESTCOOL) CALL XABORT('@MCC: LOCAL PARAMETER:'
     >                //' T-COOL DOES NOT EXIST IN THE FUEL MAP AND'
     >                //' IS REQUIRED TO COMPUTE D-COOL.')
        IF(.NOT.PRESDCOOL) CALL XABORT('@MCC: LOCAL PARAMETER:'
     >                //' D-COOL DOES NOT EXIST IN THE FUEL MAP.')
        IF(PTYPETCOOL.NE.PTYPEDCOOL) CALL XABORT('@MCC: T-COOL AND'
     >                //' D-COOL HAVE DIFFERENT TYPES (ONE IS GLOBAL'
     >                //' AND THE OTHER IS LOCAL...).')
* Definition of the pressure table size (the same as T-COOL table)
        DO IVAL=1,VALSIZE,1
          CALL THMSAT(PINLET,TSAT)
          IF(VALTCOOL(IVAL).GT.TSAT) CALL XABORT('@MCC: WATER TEMPERA'
     >           //'TURE IS GREATER THAN SATURATION TEMPERATURE (COO'
     >           //'LANT IS BOILING).')
          IF(VALTCOOL(IVAL).LT.273) CALL XABORT('@MCC: WATER TEMPERA'
     >           //'TURE IS LOWER THAN 273K (FROZEN) IN SOME REGIONS.')
          CALL THMPT(PINLET,VALTCOOL(IVAL),VALDCOOL(IVAL),R1,R2,R3,R4)
          VALDCOOL(IVAL)=VALDCOOL(IVAL)/1000
        ENDDO
* Replacement of the old D-COOL values by the new ones
        JPMAP=LCMGID(IPMAP,'PARAM')
        DO IPAR=1,NPARAM,1
          KPMAP=LCMGIL(JPMAP,IPAR)
          CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
          IF(PNAME.EQ.'D-COOL') THEN
            CALL LCMPUT(KPMAP,'P-VALUE',VALSIZE,2,VALDCOOL)
            EXIT
          ENDIF
        ENDDO
        IF(IMPX.GE.1) WRITE(6,'(1X,A/)') 'PARAMETER D-COOL HAS BEEN CO'
     >                //'MPUTED FROM T-COOL USING THE WATER TABLES.'
*
      ELSE IF(TEXT.EQ.';') THEN
        GO TO 20
      ELSE
        CALL XABORT('@MCC: INVALID KEYWORD: '//TEXT//'.')
      ENDIF
      GO TO 10
*
   20 RETURN
      END
