*DECK HSTGET
      SUBROUTINE HSTGET(IPHST,  IPRINT, MAXG,   MAXL,   NCHA,   NBUN,
     >                  ITYPRO, ITYRED, CARRED, IUPDC,  IUPDB,
     >                  NAMG,   PARAMG, NAML,   PARAML, IDCELL, IDFUEL)
*
*----------
*
*Purpose:
* To read from the input file or send to CLE-2000 variables the 
* local and burnup parameters associated with a fuel cell.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IPHST   address of the \dds{history} data structure.
* IPRINT  print level.
* MAXG    maximum number of global parameters.                   
* MAXL    maximum number of local parameters.                   
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.
* ITYPRO  type of processing where:
*         ITYPRO > 0 if history is in creation or update mode; 
*         ITYPRO < 0 if history is in read-only mode. 
* ITYRED  type of the last variable read.                
* CARRED  last character string read.
*
*Parameters: input/output
* NMAG    global parameter names.
* PARAMG  values of the global parameters.
* NMAL    local parameter names.
* PARAML  values of the local parameters.
* IDCELL  cell identifier for each fuel bundle in each channel.
* IDFUEL  fuel type identifier for each fuel bundle in each channel.
*
*Parameters: output
* IUPDC   number of the channel to analyze.                   
* IUPDB   number of the bundle to analyze.
*
*----------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST
      INTEGER          IPRINT,MAXG,MAXL,NCHA,NBUN,ITYPRO
      INTEGER          ITYRED,IUPDC,IUPDB
      CHARACTER        CARRED*12
      INTEGER          NAMG(3,0:MAXG),NAML(3,0:MAXL)
      REAL             PARAMG(0:MAXG),PARAML(0:MAXL,2)
      INTEGER          IDCELL(NBUN,NCHA),IDFUEL(NBUN,NCHA)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,NTC,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NTC=3,ILCMUP=1,ILCMDN=2,
     >                 NAMSBR='HSTGET')
*----
*  INPUT/OUTPUT VARIABLES 
*  Input data is of the form
*  [ GET (hstpar) ] [ PUT (hstpar) ]
*  [ CELLID icha ibun  [ idfuel ] 
*      [ GET (hstpar) ] 
*      [ PUT { BREFL  (hsrbrn) (hstpar) 
*              AREFL   (hsrbrn) (hstpar)   |
*            [ AREFL ] (hsrbrn) (hstpar)   } ] ]
*
*  HERE:
*  (hstpar)           = NAMPAR valpar
*                       where NAMPAR is the name of a local or global
*                       parameter and valpar its value.
*  (hstbrn)           = BURN period power
*                       where period is the burnup time step 
*                       and power the burnup power density in kW/kg.
*  For global parameter:
*  GET                = implies that (hstpar) is transfered to the 
*                       HISTORY file, 
*  PUT                = implies that (hstpar) is transfered to
*                       CLE-2000 variables. 
*  For local parameters:
*  GET                = implies that (hstpar) is transfered to the 
*                       HISTORY file for the case before and
*                       after refueling.
*  PUT                = implies that (hstbrn) and (hstpar)  
*                       are transfered to CLE-2000 variables. 
*  BREFL              = Indicates that the information before
*                       refueling is considered.
*  AREFL              = Indicates that the information after
*                       refueling is considered.
*                       This is the default option is neither
*                       BREFL nor AREFL is defined.
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
      INTEGER          ITYPUT,INTPUT
      CHARACTER        CARPUT*12
      REAL             REAPUT
      DOUBLE PRECISION DBLPUT
*----
*  LOCAL VARIABLES
*----                  
      INTEGER          ICONTR,IGP,IFTN,ISREF,IUPDL,IUPDG,IUPDF
      INTEGER          ITC,INEXT,IB,IC,IPL,IP
      INTEGER          ICT,IOK
      CHARACTER        NAMP*12
      REAL             TIMPOW(2,2)
*----
*  Initialize input vectors
*----
      CALL XDRSET(PARAML,(MAXL+1)*2,0.0)
      CALL XDRSET(TIMPOW,       2*2,0.0)
*----
*  Initialize variables
*  IUPDC  -> channel number to process or update.
*  IUPDB  -> bundle number to process or update.
*  ICONTR -> indicates processing of ITYRED and CARRED
*            = 0 processing required.
*            = 1 processing has been performed.
*  IGP    -> indicate if a GET or PUT command is in effect.
*            =-1 PUT command in effect
*            = 0 no GET or PUT command in effect
*            = 1 GET command in effect 
*  IFTN      = new fuel type 
*  ISREF  -> indicate the REFUEL state
*            is to be processed
*            = 0 no processing
*            = 1 processing before refuel
*            = 2 processing after refuel 
*  IUPDL  -> indicates local parameters update 
*            = 0 no update
*            > 0 updated
*  IUPDG  -> indicates global parameters update 
*            = 0 no update
*            > 0 updated 
*  IUPDF  -> Fuel type update
*            = 0 no update
*            > 0 updated
*----  
      IUPDC=0
      IUPDB=0
      ICONTR=0
      IGP   =0
      IFTN  =0
      ISREF =0
      IUPDL =0 
      IUPDG =0
      IUPDF =0 
 100  CONTINUE
      IF(ICONTR .EQ. 0) THEN
        ITYPLU=ITYRED
        CARLIR=CARRED
        ICONTR=1
      ELSE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      ENDIF
 101  CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      IF(CARLIR .EQ. ';') THEN
        GO TO 105 
      ELSE IF(CARLIR .EQ. 'CELLID') THEN
        IGP=0
*----
*  Channel number
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer value for channel number expected.')
        IF(INTLIR .LT. 0 ) CALL XABORT(NAMSBR//
     >  ': Read error -- value for channel number must be > 0.')
        IF(IUPDC .NE. 0) CALL XABORT(NAMSBR//
     >  ': Only one channel can be updated for each call to HST.')
        IUPDC=INTLIR
*----
*  Bundle number
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': Read error -- integer value for bundle number expected.')
        IF(INTLIR .LT. 0 ) CALL XABORT(NAMSBR//
     >  ': Read error -- value for bundle number must be > 0')
        IF(IUPDB .NE. 0) CALL XABORT(NAMSBR//
     >  ': Only one bundle can be updated for each call to HST.')
        IUPDB=INTLIR 
*----
*  Fuel type (optional)
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR) 
        IFTN=-1
        IF(ITYPLU .EQ. 1) THEN
          IFTN=INTLIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        ENDIF
*----
*  IF CELL IS NOT IDENTIFIED ASSOCIATE TO CELL NEXT 
*  CELL NUMBER AVAILABLE AND TO FUEL TYPE
*  VALUE PROVIDED IN IFTN
*----   
        IF(IDCELL(IUPDB,IUPDC) .LE. 0) THEN
          DO 110 INEXT=1,NBUN*NCHA
            DO 111 IB=1,NBUN
              DO 112 IC=1,NCHA
                IF(IDCELL(IB,IC) .EQ. INEXT) GO TO 115
 112          CONTINUE
 111        CONTINUE
            IDCELL(IUPDB,IUPDC)=INEXT
            GO TO 116
 115        CONTINUE
 110      CONTINUE
          CALL XABORT(NAMSBR//': No cell id available')
 116      CONTINUE
          IDFUEL(IUPDB,IUPDC)=ABS(IFTN) 
        ELSE
*----
*  CELL EXIST, READ IF POSSIBLE EXISTING LOCAL 
*  PARAMETERS VALUES
*----     
          ICT=IDCELL(IUPDB,IUPDC)
          WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
          CALL LCMSIX(IPHST,NAMP,ILCMUP)
*----
*  Get local parameters from cell before refueling
*----
          IOK=-1              
          CALL HSTGSL(IPHST ,MAXL  ,IOK   ,
     >                TIMPOW(1,1)  ,PARAML(0,1)) 
          IF((IPRINT.GT.0).AND.(IOK.NE.0)) THEN
            WRITE(IOUT,7000) NAMSBR
            WRITE(IOUT,7010) IUPDC,IUPDB,'BEFORE'
          ENDIF
*----
*  Get local parameters from cell after refueling
*----
          IOK=-2              
          CALL HSTGSL(IPHST ,MAXL  ,IOK   ,
     >                TIMPOW(1,2)  ,PARAML(0,2)) 
          IF((IPRINT.GT.0).AND.(IOK.NE.0)) THEN
            WRITE(IOUT,7000) NAMSBR
            WRITE(IOUT,7010) IUPDC,IUPDB,'AFTER '
          ENDIF
          CALL LCMSIX(IPHST,NAMP,ILCMDN)
        ENDIF
        GO TO 101   
      ELSE IF(CARLIR .EQ. 'GET') THEN
        IF(ITYPRO .LT. 0) CALL XABORT(NAMSBR//
     >': Option GET not permitted for history in read only mode')
        IGP=1
        ISREF=2
      ELSE IF(CARLIR .EQ. 'PUT') THEN
        IGP=-1
        ISREF=2
      ELSE IF(CARLIR .EQ. 'BREFL') THEN
        IF(IGP .NE. -1) CALL XABORT(NAMSBR//
     >': Option BREFL permitted for PUT only')
        ISREF=1
      ELSE IF(CARLIR .EQ. 'AREFL') THEN
        IF(IGP .NE. -1) CALL XABORT(NAMSBR//
     >': Option AREFL permitted for PUT only')
        ISREF=2
      ELSE
        IF(IGP .EQ. 0) CALL XABORT(NAMSBR//
     >  ': GET or PUT must be specified ')
        IF(IUPDC*IUPDB .GT. 0) THEN
*----
*  CARLIR contains a local parameter
*----
          IF(CARLIR .EQ. 'BURN') THEN 
            IF(IGP .EQ. 1) CALL XABORT(NAMSBR//
     >': Option GET not permitted for BURN keyword')
            IF(ITYPRO .GT. 0) CALL XABORT(NAMSBR//
     >': Option BURN permitted only for history in read only mode')
            REAPUT=TIMPOW(1,ISREF)
            ITYPUT=2
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU .NE. -ITYPUT) CALL XABORT(NAMSBR//
     >': Real output variable for burnup period expected')
            CALL REDPUT(ITYPUT,INTPUT,REAPUT,CARPUT,DBLPUT)
*----
*  The power density expected is in kW/kg.
*----
            REAPUT=TIMPOW(2,ISREF)
            ITYPUT=2
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU .NE. -ITYPUT) CALL XABORT(NAMSBR//
     >': Real output variable for burnup power expected')
            CALL REDPUT(ITYPUT,INTPUT,REAPUT,CARPUT,DBLPUT)
          ELSE
*----
*  Scan local parameters to see is CARLIR is one of them
*----       
            IP=0
            DO 120 IPL=1,MAXL
              WRITE(NAMP,'(3A4)') (NAML(ITC,IPL),ITC=1,NTC)
              IF(NAMP .EQ. CARLIR) THEN
                IP=IPL
                GO TO 125
              ELSE IF(NAMP .EQ. '            ') THEN 
                IP=IPL
                READ(CARLIR,'(3A4)') (NAML(ITC,IP),ITC=1,NTC)
                GO TO 125
              ENDIF
 120        CONTINUE
            CALL XABORT(NAMSBR//': Number of local parameters '//
     >      'provided larger than number permitted.') 
 125        CONTINUE
            IF(IGP .EQ. -1) THEN
              REAPUT=PARAML(IP,ISREF)
              ITYPUT=2
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
              IF(ITYPLU .NE. -ITYPUT) CALL XABORT(NAMSBR//
     >':  Real output variable for local parameter expected')
              CALL REDPUT(ITYPUT,INTPUT,REAPUT,CARPUT,DBLPUT)
            ELSE IF(IGP .EQ. 1) THEN
              IUPDL=IUPDL+1
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
              IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >        ': Real value for local parameter missing.')
              PARAML(IP,ISREF)=REALIR 
            ENDIF
          ENDIF
        ELSE
*----
*  CARLIR contains a global parameter
*----
          IF(CARLIR .EQ. 'POWER') THEN
            CALL XABORT(NAMSBR//
     >      ': POWER is a local not global parameter') 
          ELSE 
            IP=0
            DO 130 IPL=1,MAXG
              WRITE(NAMP,'(3A4)') (NAMG(ITC,IPL),ITC=1,NTC)
              IF(NAMP .EQ. CARLIR) THEN
                IP=IPL
                GO TO 135
              ELSE IF(NAMP .EQ. '            ') THEN 
                IP=IPL
                READ(CARLIR,'(3A4)') (NAMG(ITC,IP),ITC=1,NTC)
                GO TO 135
              ENDIF
 130        CONTINUE
            CALL XABORT(NAMSBR//': Number of global parameters '//
     >      'provided larger than number permitted.') 
 135        CONTINUE
            IF(IGP .EQ. -1) THEN
              REAPUT=PARAMG(IP)
              ITYPUT=2
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
              IF(ITYPLU .NE. -ITYPUT) CALL XABORT(NAMSBR//
     >':  Real output variable for global parameter expected')
              CALL REDPUT(ITYPUT,INTPUT,REAPUT,CARPUT,DBLPUT)
            ELSE IF(IGP .EQ. 1) THEN
              IUPDG=IUPDG+1
              CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
              IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >        ': Real value for global parameter missing.')
              PARAMG(IP)=REALIR 
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  Save global parameters if some are updated
*----
      IF(IUPDG .GT. 0) THEN 
        CALL LCMPUT(IPHST,'NAMEGLOBAL  ',3*MAXG,3,NAMG(1,1))
        CALL LCMPUT(IPHST,'PARAMGLOBAL ',MAXG,2,PARAMG(1))
      ENDIF       
      IF(IUPDL .GT. 0) THEN 
        CALL LCMPUT(IPHST,'NAMELOCAL   ',3*MAXL,3,NAML(1,1))
        ICT=IDCELL(IUPDB,IUPDC)
        WRITE(NAMP,'(A4,I8.8)') 'CELL',ICT
        CALL LCMSIX(IPHST,NAMP,ILCMUP)
        IOK=2              
        CALL HSTGSL(IPHST ,MAXL  ,IOK   ,
     >              TIMPOW(1,2)  ,PARAML(0,2)) 
        CALL LCMSIX(IPHST,NAMP,ILCMDN)
      ENDIF       
      RETURN 
*----
*  Formats
*  WARNING
*----
 7000 FORMAT(' ***** WARNING IN ',A6,' *****')
 7010 FORMAT(' Local parameters for channel ',I5,' bundle ',I5,
     >       ' not available for ',A6,' state'/
     >       ' Initialize to 0.0')
      END 
