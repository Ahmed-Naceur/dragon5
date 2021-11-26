*DECK TINST
      SUBROUTINE TINST(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform computations according to the time-linear model.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
*
*Author(s): 
* B. Toueg, M. Guyot
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
* The TINST: module specification is:
* Option 1
* FMAP := TINST: FMAP [ POWER ] :: (desctinst) ;
* Option 2  
* MICLIB3 FMAP := TINST: FMAP MICLIB2 MICLIB :: (desctinst) ;
* where
*   FMAP : name of a \emph{fmap} object, that will be updated by the TINST: 
*     module. The FMAP object must contain the instantaneous burnups for each 
*     fuel bundle and the weight of each fuel mixture.
*   POWER : name of a \emph{power} object containing the channel and bundle 
*     powers, previously computed by the FLPOW: module. The channel and bundle 
*     powers are used by the TINST: module to compute the new burn-up of each 
*     bundle. If bundle-powers are previously specified with the module RESINI:,
*     you can refuel your core without a POWER object.
*   MICLIB3 : name of a \emph{library} object, that will be created by the 
*     TINST: module. This \emph{MICROLIB} contains the fuel properties after 
*     refueling when keyword MICRO is used in (desctinst).
* 
*   MICLIB2 : name of a \emph{library} object, that will be read by the TINST: 
*     module. This must be a fuel-map LIBRARY created either created by the 
*     NCR: or the EVO: module.
*   MICLIB : name of a \emph{library} object, that will be read by the TINST: 
*     module. This \emph{MICROLIB} contains the new fuel properties, that 
*     should be used for the refueling.
*   (desctinst) : structure describing the input data to the TINST: module.
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
      CHARACTER TEXT*12,HSIGN*12,NAMCHA*8,NAMCHA2*8,TEXT12*12
      INTEGER IMPX,MSHT,NB,NCH,IMOD,NCOMB,NF,NX,NY,NZ,MAXS,ITYP,
     +        NREG,KREF,I,ISTATE(NSTATE),LENGT,NS,NSS,NITMA
      DOUBLE PRECISION DFLOT
      LOGICAL LNOTHING,LMIC
      REAL TIME,BURNSTEP,FLOT
      TYPE(C_PTR) IPMAP,IPPOW,IPMIC,IPMIC2,IPMIC3,JPMAP,KPMAP,LPMAP,
     +            MPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDX,NSSV,IXN,IYN,MIX,IVS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHN
      REAL, ALLOCATABLE, DIMENSION(:) :: BUNDPOW,WINT,BS,PS,POW,RFCHAN,
     + BURNINST
*----
*  PARAMETER VALIDATION
*----
      IF((NENTRY.LT.1).AND.(NENTRY.GT.5)) CALL XABORT('@TINST: WRONG '
     +  //'NUMBER OF PARAMETERS')
      IPMAP=C_NULL_PTR
      IPPOW=C_NULL_PTR
      IPMIC=C_NULL_PTR
      IPMIC2=C_NULL_PTR
      IPMIC3=C_NULL_PTR
      IF(JENTRY(1).EQ.0) THEN
        IPMIC=KENTRY(1)
        I=2
        TEXT12=HENTRY(1)
        IF(IENTRY(1).GT.2) CALL XABORT('@TINST: LCM OR XSM OBJECT TYPE'
     +  //' FOR ENTRY='//TEXT12//'.')
      ELSE
        I=1
      ENDIF
      DO IEN=I,NENTRY
        TEXT12=HENTRY(IEN)
        IF(IENTRY(IEN).GT.2) CALL XABORT('@TINST: LCM OR XSM OBJECT TY'
     +  //'PE FOR ENTRY='//TEXT12//'.')
        CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.EQ.'L_MAP') THEN
          IPMAP=KENTRY(IEN)
          IF(JENTRY(IEN).NE.1) CALL XABORT('@TINST: MODIFICATION MODE '
     +    //'FOR L_MAP EXPECTED')
        ELSEIF(HSIGN.EQ.'L_POWER') THEN
          IPPOW=KENTRY(IEN)
          IF(JENTRY(IEN).NE.2) CALL XABORT('@TINST: READ-ONLY MODE '
     +    //'FOR L_POWER EXPECTED')
        ELSEIF(HSIGN.EQ.'L_LIBRARY') THEN
          IF(.NOT.C_ASSOCIATED(IPMIC2)) THEN
            IPMIC2=KENTRY(IEN)
            CALL LCMEQU(IPMIC2,IPMIC)
            IF(JENTRY(IEN).NE.2) CALL XABORT('@TINST: READ-ONLY MODE'
     +      //' FOR SECOND L_LIBRARY EXPECTED')
          ELSE
            IPMIC3=KENTRY(IEN)
            IF(JENTRY(IEN).NE.2) CALL XABORT('@TINST: READ-ONLY MODE '
     +      //'FOR THIRD L_LIBRARY EXPECTED')
          ENDIF
        ENDIF          
      ENDDO
*----
*  RECOVER INFORMATION
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      NB=ISTATE(1)
      NCH=ISTATE(2)
      NCOMB=ISTATE(3)
      IMOD=ISTATE(5)
      MAXS=ISTATE(6)
      MSHT=MAXS+1
      NF=ISTATE(7)
      NPARM=ISTATE(8)
      IF(NF.EQ.0) CALL XABORT('@TINST: NO FUEL IN MAP OBJECT.')
*----
*  ONLY TIME INSTANTANEOUS CALCULATIONS IN TINST:
*----
      IF(IMOD.NE.2)
     +    CALL XABORT('@TINST: INST-BURN OPTION '
     + //'SHOULD BE USED IN RESINI.')
      JPMAP=LCMGID(IPMAP,'GEOMAP')
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(JPMAP,'STATE-VECTOR',ISTATE)
      IGEO=ISTATE(1)
      NX=ISTATE(3)
      NY=ISTATE(4)
      NZ=ISTATE(5)
      NREG = ISTATE(6)
*     CHECK EXISTING DATA
      IF(.NOT.C_ASSOCIATED(IPPOW)) THEN
        CALL LCMLEN(IPMAP,'BUND-PW',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@TINST: MISSING BUND-PW DATA IN '
     +   //'L_MAP OBJECT.')
      ELSE
        CALL LCMLEN(IPPOW,'POWER-CHAN',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@TINST: MISSING POWER-CHAN DATA I'
     +   //'N L_POWER OBJECT.')
      ENDIF
*----
*  READ INPUT DATA
*----
      IMPX=0
      LNOTHING=.TRUE.
      LMIC=.FALSE.
      TTIME=0.0
      ALLOCATE(RFCHAN(NCH))
      CALL XDRSET(RFCHAN,NCH,0.0)
   2  TIME=0.0
      BURNSTEP=0.0
*     READ KEYWORD
   1  CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@TINST: CHARACTER DATA EXPECTED(1).')
      IF(TEXT.EQ.'EDIT')THEN
*       PRINTING INDEX
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@TINST: INTEGER DATA EXPECTED.')
        IMPX=MAX(0,NITMA)
        GOTO 1
      ELSEIF(TEXT.EQ.'TIME')THEN
*       TIME VALUE
        IF(TIME.NE.0.0)CALL XABORT('@TINST: TIME ALREADY SPECIFIED(1).')
        IF(BURNSTEP.NE.0.0)CALL XABORT('@TINST: BURNSTEP ALREADY //
     +    //SPECIFIED(1).')
        CALL REDGET(ITYP,NITMA,TIME,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@TINST: REAL DATA EXPECTED(1).')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@TINST: CHARACTER DATA EXPECTED(2).')
        IF(TIME.LT.0.)CALL XABORT('@TINST: EXPECTING REAL > 0 (1).')
        IF(TEXT.EQ.'DAY')THEN
          TIME=TIME
        ELSEIF(TEXT.EQ.'HOUR')THEN
          TIME=TIME/24.
        ELSEIF(TEXT.EQ.'MINUTE')THEN
          TIME=TIME/(24.*60.)
        ELSEIF(TEXT.EQ.'SECOND')THEN
          TIME=TIME/(24.*60.*60.)
        ELSE
          CALL XABORT('@TINST: EXPECTING DAY|HOUR|MINUTE|SECOND.')
        ENDIF
        LNOTHING=.FALSE.
        GOTO 10
      ELSEIF(TEXT.EQ.'BURN-STEP')THEN
*       BURN-STEP
        IF(TIME.NE.0.)CALL XABORT('@TINST: TIME ALREADY SPECIFIED(2).')
        IF(BURNSTEP.NE.0.)CALL XABORT('@TINST: BURNSTEP ALREADY '
     +    //'SPECIFIED(2).')
        CALL REDGET(ITYP,NITMA,BURNSTEP,TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@TINST: REAL DATA EXPECTED(2).')
        IF(BURNSTEP.LE.0.)CALL XABORT('@TINST: EXPECTING REAL > 0 (2).')
        LNOTHING=.FALSE.
        GOTO 10
      ELSEIF(TEXT.EQ.'REFUEL')THEN
*      REFUEL
           KREF=1
           LNOTHING=.FALSE.
           CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
           IF(ITYP.NE.3)
     +      CALL XABORT('@TINST: CHARACTER DATA EXPECTED(3).')
           IF(TEXT.EQ.'MICRO') THEN
             LMIC=.TRUE.
             CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
           IF(ITYP.NE.3)
     +      CALL XABORT('@TINST: CHARACTER DATA EXPECTED(4).')
           ENDIF
           IF(TEXT.EQ.'CHAN') THEN
             CALL REDGET(ITYP,NITMA,FLOT,NAMCHA,DFLOT)
             IF(ITYP.NE.3)
     +        CALL XABORT('@TINST: CHARACTER DATA EXPECTED(5).')
             CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
             IF(ITYP.NE.1)
     +         CALL XABORT('@TINST: INTEGER DATA EXPECTED(2).')
             NS = NITMA
             CALL TINCHA(IPMAP,NCH,IMPX,NAMCHA,TTIME,RFCHAN)
           ELSE
            CALL XABORT('@TINST: INVALID KEYWORD '//TEXT)
           ENDIF
       GOTO 20
      ELSEIF(TEXT.EQ.'NEWFUEL')THEN
*      NEWFUEL
           KREF=2
           LNOTHING=.FALSE.
           CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
           IF(ITYP.NE.3)
     +      CALL XABORT('@TINST: CHARACTER DATA EXPECTED(4).')
           IF(TEXT.EQ.'CHAN') THEN
             CALL REDGET(ITYP,NITMA,FLOT,NAMCHA,DFLOT)
             IF(ITYP.NE.3)
     +        CALL XABORT('@TINST: CHARACTER DATA EXPECTED(5).')
             CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
             IF(ITYP.NE.1)
     +         CALL XABORT('@TINST: INTEGER DATA EXPECTED(3).')
             NS = NITMA
             NSS=ABS(NS)
             ALLOCATE(IDX(NSS))
             CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
             IF(ITYP.NE.3)
     +         CALL XABORT('@TINST: CHARACTER DATA EXPECTED(6).')
             IF(TEXT.EQ.'SOME')THEN
                DO 11 I=1,NSS
                 CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
                 IF(ITYP.NE.1)
     +             CALL XABORT('@TINST: INTEGER DATA EXPECTED(4).')
                IF (NITMA.GT.NF)
     +             CALL XABORT('@TINST: WRONG NUMBER OF FUEL TYPE. ')
                 IDX(I) = NITMA
 11             CONTINUE
             ELSEIF(TEXT.EQ.'ALL')THEN
                CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
                IF(ITYP.NE.1)
     +             CALL XABORT('@TINST: INTEGER DATA EXPECTED(5).')
                IF (NITMA.GT.NF)
     +             CALL XABORT('@TINST: WRONG NUMBER OF FUEL TYPE. ')
                DO 12 I=1,NSS
                 IDX(I) = NITMA
 12             CONTINUE
             ELSE
               CALL XABORT('@TINST: INVALID KEYWORD '//TEXT)
             ENDIF
           ELSE
            CALL XABORT('@TINST: INVALID KEYWORD '//TEXT)
           ENDIF
        GOTO 20
*     SHUFFL
      ELSEIF (TEXT.EQ.'SHUFF') THEN
           KREF = 3
           LNOTHING=.FALSE.
           CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
           IF(ITYP.NE.3)
     +      CALL XABORT('@TINST: CHARACTER DATA EXPECTED(7).')
           IF(TEXT.EQ.'CHAN') THEN
             CALL REDGET(ITYP,NITMA,FLOT,NAMCHA,DFLOT)
             IF(ITYP.NE.3)
     +        CALL XABORT('@TINST: CHARACTER DATA EXPECTED(8).')
             CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
             IF(ITYP.NE.3)
     +        CALL XABORT('@TINST: CHARACTER DATA EXPECTED(8).')
             IF(TEXT.EQ.'TO') THEN
                CALL REDGET(ITYP,NITMA,FLOT,NAMCHA2,DFLOT)
                IF(ITYP.NE.3)
     +            CALL XABORT('@TINST: CHARACTER DATA EXPECTED(9).')
                IF(IMPX.GT.2)
     +            WRITE(6,*) 'TINST : ACTION ',NAMCHA,' TO ',NAMCHA2
             ELSE
               CALL XABORT('@TINST: INVALID KEYWORD '//TEXT)
             ENDIF
           ELSE
             CALL XABORT('@TINST: INVALID KEYWORD '//TEXT)
           ENDIF
        GOTO 20
      ELSEIF(TEXT.EQ.'PICK') THEN
*       RECOVER THE BURNUP AND SAVE IT IN A CLE-2000 VARIABLE
        IF(IMPX.GT.2) WRITE(IOUT,40) BURNAVG
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.-2) CALL XABORT('TINST: OUTPUT REAL EXPECTED.')
        ITYP=2
        CALL REDPUT(ITYP,NITMA,BURNAVG,TEXT,DFLOT)
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF((ITYP.NE.3).OR.(TEXT.NE.';')) THEN
          CALL XABORT('TINST: ; CHARACTER EXPECTED.')
        ENDIF      
        GOTO 30
      ELSEIF(TEXT.EQ.';')THEN
        GOTO 30
      ELSE
*       KEYWORD DOES NOT MATCH
        CALL XABORT('@TINST: WRONG KEYWORD: '//TEXT//'.')
      ENDIF
*----
*  PERFORM CALCULATION
*----
   10 ALLOCATE(BUNDPOW(NCH*NB))
      IF(.NOT.C_ASSOCIATED(IPPOW)) THEN
        CALL LCMGET(IPMAP,'BUND-PW',BUNDPOW)
      ELSE 
        CALL LCMGET(IPPOW,'POWER-BUND',BUNDPOW)
      ENDIF
      IF(LMIC) CALL XABORT('@TINST: NO MICRO-DEPLETION ')
      TTIME  = TTIME + TIME
      ALLOCATE(BURNINST(NCH*NB))
      CALL TINSTB(IPMAP,TIME,BURNSTEP,NCH,NB,NF,BUNDPOW,BURNAVG,
     1 BURNINST,IMPX)
*----
*  SAVE LOCAL PARAMETERS FOR HISTORICAL FOLLOW-UP
*----
      CALL LCMLEN(IPMAP,'_TINST',ILONG,ITYLCM)
      IF(IMPX.GT.0) WRITE(6,50) ILONG+1,BURNAVG
      JPMAP=LCMLID(IPMAP,'_TINST',ILONG+1)
      KPMAP=LCMDIL(JPMAP,ILONG+1)
      CALL LCMPUT(KPMAP,'TIME',1,2,TTIME)
      CALL LCMPUT(KPMAP,'BURNAVG',1,2,BURNAVG)
      CALL LCMPUT(KPMAP,'BURN-INST',NCH*NB,2,BURNINST)
      CALL LCMPUT(KPMAP,'POWER-BUND',NCH*NB,2,BUNDPOW)
      IF(NPARM.GT.0) THEN
        LPMAP=LCMGID(IPMAP,'PARAM')
        MPMAP=LCMLID(KPMAP,'PARAM',NPARM)
        CALL LCMEQU(LPMAP,MPMAP)
        ISTATE(19)=1
        CALL LCMPTC(IPMAP,'CYCLE-NAMES',12,1,'_TINST')
      ENDIF
      DEALLOCATE(BURNINST,BUNDPOW)
      GOTO 1
*
   20 CALL LCMSIX(IPMAP,' ',0)
      ALLOCATE(NSSV(NCH))
      CALL LCMLEN(IPMAP,'REF-SCHEME',ILONG,ITYP)
      IF(ILONG.EQ.0) THEN
         DO 25 I=1,NCH
            NSSV(I) = 0
   25    CONTINUE
      ELSEIF(ILONG.NE.NCH) THEN
         CALL XABORT('@TINST: REF-SCHEME HAS NOT THE CORRECT LENGHT')
      ELSE
         CALL LCMGET(IPMAP,'REF-SCHEME',NSSV)
      ENDIF
      CALL LCMSIX(IPMAP,' ',0)
      IF(IGEO.EQ.7) THEN
*       Cartesian geometry.
        ALLOCATE(IXN(NX),IYN(NY))
        CALL LCMGET(IPMAP,'XNAME',IXN)
        CALL LCMGET(IPMAP,'YNAME',IYN)
        ALLOCATE(WINT(NCH*NB),MIX(NREG),BS(NCH*NB*MSHT),PS(NCH*NB*MSHT),
     1  IVS(NCH*NB*MSHT))
        IF(KREF.EQ.1.OR.KREF.EQ.2) THEN
          IF(KREF.EQ.1) ALLOCATE(IDX(ABS(NS)))
          ALLOCATE(POW(NCH*NB))
          IF(.NOT.C_ASSOCIATED(IPPOW)) THEN
            CALL LCMGET(IPMAP,'BUND-PW',POW)
          ELSE 
            CALL LCMGET(IPPOW,'POWER-BUND',POW)
          ENDIF
          CALL TINREF(IPMAP,IPMIC,IPMIC2,IPMIC3,NCH,NB,NX,NY,NZ,NREG,
     +               NAMCHA,NS,MSHT,WINT,MIX,IXN,IYN,BS,PS,IVS,POW,
     +               MAXS,NSSV,IDX,IMPX,KREF,LMIC)
          DEALLOCATE(POW,IDX)
        ELSE
          CALL TINSHU(IPMAP,NCH,NB,NX,NY,NZ,NREG,MSHT,NAMCHA,NAMCHA2,
     +       WINT,MIX,BS,PS,IVS,IXN,IYN,IMPX)
        ENDIF
        DEALLOCATE(IXN,IYN)
      ELSE IF(IGEO.EQ.9) THEN
*       Hexagonal geometry.
        ALLOCATE(IHN(2,NX))
        CALL LCMGET(IPMAP,'HNAME',IHN)
        ALLOCATE(WINT(NCH*NB),MIX(NREG),BS(NCH*NB*MSHT),PS(NCH*NB*MSHT),
     1  IVS(NCH*NB*MSHT))
        IF(KREF.EQ.1.OR.KREF.EQ.2) THEN
          IF(KREF.EQ.1) ALLOCATE(IDX(ABS(NS)))
          ALLOCATE(POW(NCH*NB))
          IF(.NOT.C_ASSOCIATED(IPPOW)) THEN
            CALL LCMGET(IPMAP,'BUND-PW',POW)
          ELSE 
            CALL LCMGET(IPPOW,'POWER-BUND',POW)
          ENDIF
          CALL TINREH(IPMAP,IPMIC,IPMIC2,IPMIC3,NCH,NB,NX,NZ,NREG,
     +               NAMCHA,NS,MSHT,WINT,MIX,IHN,BS,PS,IVS,POW,MAXS,
     +               NSSV,IDX,IMPX,KREF,LMIC)
          DEALLOCATE(POW,IDX)
        ELSE
          CALL TINSHH(IPMAP,NCH,NB,NX,NZ,NREG,MSHT,NAMCHA,NAMCHA2,
     +       WINT,MIX,BS,PS,IVS,IHN,IMPX)
        ENDIF
        DEALLOCATE(IHN)
      ELSE
        CALL XABORT('TINST: GEOMETRY TYPE NOT SUPPORTED')
      ENDIF
      DEALLOCATE(BS,PS,IVS,MIX)
      DEALLOCATE(WINT)
      DEALLOCATE(NSSV)
      MSHT=MAXS+1
      KREF=0
      GOTO 2
*
  30  IF(LNOTHING)CALL XABORT('@TINST: NO OPTION SPECIFIED.')
      CALL LCMSIX(IPMAP,' ',0)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      ISTATE(6)=MAXS
      CALL LCMPUT(IPMAP,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAP,'DEPL-TIME',1,2,TTIME)
      CALL LCMPUT(IPMAP,'REF-CHAN',NCH,2,RFCHAN)
      DEALLOCATE(RFCHAN)
      RETURN
   40 FORMAT(/20H TINST: PICK BURNUP=,1P,E12.4,10H MWd/tonne)
   50 FORMAT(/38H TINST: STORE INFORMATION IN LIST ITEM,I3,9H OF TINST,
     + 20H DIRECTORY AT BURNUP,1P,E12.4,8H MW-D/T./)
      END
