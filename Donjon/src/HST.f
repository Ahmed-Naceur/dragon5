*DECK HST
      SUBROUTINE HST(NENTRY, HENTRY, IENTRY, JENTRY, KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract from or save to a \dds{history} data structure 
* the information related to various cells in a reactor.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau, E. Varin
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
* For HST:, the possible calling specifications  are:
* Option 1: Updating an \emph{history} structure using a \emph{map} structure
* history := HST: [ history ]  map [ :: [ (hstdim) ]  [ GET (hstpar) ] ] ;
* Option 2: Updating an \emph{history} structure using a \emph{burnup} structure
* history := HST: [ history ]  [ burnup ] [ :: [ (hstdim) ]
*   [ GET (hstpar) ] [ CELLID icha ibun [ idfuel ] [ GET (hstpar) ] ] ] ; 
* Option 3: Updating a \emph{burnup} structure using an \emph{history} structure
* burnup := HST: history [ :: [ (hstdim ] 
*    [ PUT (hstpar) ]
*    CELLID icha ibun
*    [ PUT { BREFL  (hstbrn) (hstpar) AREFL (hstbrn) (hstpar) 
*            | [ AREFL ] (hstbrn) (hstpar) } ] ] ;
* Option 4: Updating a \emph{map} data structure from the information available 
*   on an \emph{history} data structure:
* map := HST:  map  history ;
* where
*   history : name of an \emph{history} data structure. 
*   burnup  : name of a \emph{burnup} data structure. 
*   map     : name of a \emph{map} data structure. 
*   (hstdim) : structure containing the dimensions for the \emph{history} 
*     data structure.
*   CELLID  : keyword to identify the cell for which history information is 
*     to be processed.
*   icha    : channel number for which history information is to be processed. 
*   ibun    : bundle number for which history information is to be processed. 
*   idfuel  : fuel type number associated with this cell. One can associate to 
*     each fuel cell a different fuel type. By default a single fuel type is
*     defined and it fills every fuel cell. Only the initial properties of each
*     fuel type are saved. These properties are used for refueling. 
*   GET     : keyword to specify that the values of the parameters selected in 
*     (brnpar will be read from the input stream or CLE-2000 local variables
*     and stored on the \emph{history data structure.
*   PUT     : keyword to specify that the values of the parameters selected in
*     (brnpar will be read from the \emph{history data structure and
*      transferred to local CLE-2000 variables.
*   BREFL   : to specify that the information to extract from the \emph{history}
*     data structure is related to the properties of the cell before refueling 
*     takes place.
*   AREFL   : to specify that the information to extract from the \emph{history}
*     data base is related to the properties of the cell after refueling took
*     place.
*   (hstbrn) : structure containing the burnup options.
*   (hstpar) : structure containing the local parameters options.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,ILCMUP,ILCMDN,NSTATE,NTC,MAXENT
      CHARACTER        NAMSBR*6,TEXT12*12
      PARAMETER       (IOUT=6,ILCMUP=1,ILCMDN=2,NSTATE=40,
     >                 NTC=3,MAXENT=2,NAMSBR='HST   ')
      INTEGER          NSTOLD
      PARAMETER       (NSTOLD=20)
*----
*  Debug print flag
*  IDEB = 0  -> no print debug
*       > 0  -> print debug
*----
      INTEGER          IDEB 
      PARAMETER       (IDEB=0) 
*----
*  LOCAL VARIABLES
*----
      CHARACTER        CBLANK*4,SIGENT(MAXENT)*12
      INTEGER          IBLANK,NAMTMP(NTC)
      INTEGER          ISTATB(NSTATE),ISTATH(NSTATE),ISTATM(NSTATE)
      INTEGER          ILCMLN,ILCMTY 
      INTEGER          IEN,ITC,ITYPRO
      INTEGER          IKHST,IKEVO,IKMAP
      INTEGER          NCELL,IUPDC,IUPDB
      TYPE(C_PTR)      IPHST,IPEVO,IPMAP
*----
*  HISTORY Parameters
*----
      INTEGER          MAXG,MAXL,NBUNH,NCHAH,
     >                 ITSOLH,ITBURH,MAXIH,NREGH
      REAL             BUNLEN 
*----
*  BURNUP Parameters
*----
      INTEGER          ITSOLB,ITBURB,NBBTS,MAXIB
      REAL             REVOL(5)
*----
*  MAP Parameters
*----
      INTEGER          NBUNM,NCHAM,NBFUEL 
*----
*  Variables from HSTGDM
*----
      INTEGER          IPRINT,NGLO,NLOC,NBUN,NCHA,ITYRED
      CHARACTER*12     CARRED
      INTEGER          II
*----
*  MEMORY ALLOCATION
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NAMG,NAML,IDCELL,IDFUEL,
     > IREFUS
      REAL, ALLOCATABLE, DIMENSION(:) :: PARAG,PARAL,REFUT,DENI,POWR,
     > BURN
*----
*  initialize blank signatures
*----
      DO 100 IEN=1,MAXENT
        SIGENT(IEN)='            '
 100  CONTINUE
      CBLANK='    '
      READ(CBLANK,'(A4)') IBLANK
      CALL XDISET(ISTATB,NSTATE,0)
      CALL XDISET(ISTATH,NSTATE,0)
      CALL XDISET(ISTATM,NSTATE,0)
*----
*  PARAMETER VALIDATION.
*  1 or 2 data structures permitted
*  If one data structure it must be an 
*  HISTORY structure, 
*  If two data structure, one of them must be and history
*  while the second one can be a BURNUP or MAP structure
*  Options:
*  2)   [History] := HST: [History] [Burnup] :: ... ;
*  3)    History  := HST: [History] Map      :: ... ;
*  3)    Burnup   := HST: History            :: ... ;
*---- 
      IF(NENTRY .EQ. 0) THEN
        CALL XABORT(NAMSBR//
     >': At least one data structure expected.')
      ELSE IF(NENTRY .GT. MAXENT) THEN
        CALL XABORT(NAMSBR//
     >': Maximum number of structures exceeded.')
      ENDIF
      DO 110 IEN=1,NENTRY
        TEXT12=HENTRY(IEN)
        IF(IENTRY(IEN) .NE. 1 .AND. IENTRY(IEN) .NE. 2)
     >  CALL XABORT(NAMSBR//
     >': Data structure '//TEXT12//' must be of type LCM or XSM.')
 110  CONTINUE
      IEN = 1 
      IF(JENTRY(IEN) .EQ. 2 ) THEN
        IF(NENTRY .EQ. 2) CALL XABORT(NAMSBR//
     >    ': First data structure must be in create or update mode.')
        CALL LCMLEN(KENTRY(IEN),'SIGNATURE',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0) THEN
          CALL LCMGET(KENTRY(IEN),'SIGNATURE',NAMTMP)
          WRITE(SIGENT(IEN),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
        ENDIF
      ELSE IF(JENTRY(IEN) .EQ. 1 ) THEN
        CALL LCMLEN(KENTRY(IEN),'SIGNATURE',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0) THEN
          CALL LCMGET(KENTRY(IEN),'SIGNATURE',NAMTMP)
          WRITE(SIGENT(IEN),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
        ENDIF
      ENDIF
      IF(NENTRY .EQ. 2) THEN 
        IEN = 2 
        IF(JENTRY(IEN) .NE. 2 ) CALL XABORT(NAMSBR//
     >    ': Second data structure must be in read-only mode.')
        CALL LCMLEN(KENTRY(IEN),'SIGNATURE',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0) CALL XABORT(NAMSBR//
     >': No signature found on second data structure')
        CALL LCMGET(KENTRY(IEN),'SIGNATURE',NAMTMP) 
        WRITE(SIGENT(IEN),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
      ENDIF
      IKHST=0
      IKEVO=0
      IKMAP=0
      DO 111 IEN=1,NENTRY
        IF     (SIGENT(IEN) .EQ. 'L_HISTORY   ') THEN
          IF(IKHST .NE. 0) CALL XABORT(NAMSBR//
     >      ': Two history structure forbidden.')
          IKHST=IEN
        ELSE IF(SIGENT(IEN) .EQ. 'L_BURNUP    ') THEN 
          IF(IKEVO .NE. 0) CALL XABORT(NAMSBR//
     >      ': Two burnup structure forbidden.')
          IKEVO=IEN
        ELSE IF(SIGENT(IEN) .EQ. 'L_MAP       ') THEN 
          IF(IKMAP .NE. 0) CALL XABORT(NAMSBR//
     >      ': Two map structure forbidden.')
          IKMAP=IEN
        ELSE IF(SIGENT(IEN) .NE. '            ') THEN
          CALL XABORT(NAMSBR//
     >    ': At least on structure type is invalid.')
        ENDIF
 111  CONTINUE
      BUNLEN=1.0
*----
*  For structures with SIGNATURE read STATE-VECTOR
*----
      IF(IKHST .GT. 0) THEN
        CALL LCMGET(KENTRY(IKHST),'STATE-VECTOR',ISTATH)
        CALL LCMGET(KENTRY(IKHST),'BUNDLELENGTH',BUNLEN)
      ENDIF 
      IF(IKEVO .GT. 0) THEN
        CALL LCMGET(KENTRY(IKEVO),'STATE-VECTOR',ISTATB)
        CALL LCMGET(KENTRY(IKEVO),'EVOLUTION-R ',REVOL)
      ENDIF 
      IF(IKMAP .GT. 0) THEN
        CALL LCMGET(KENTRY(IKMAP),'STATE-VECTOR',ISTATM)
      ENDIF 
*----
*  Select type of processing depending
*  on order of structures
*  ITYPRO =  1 : History := HST:                  ::
*  ITYPRO =  2 : History := HST: History          ::
*  ITYPRO =  3 : History := HST: Burnup           ::
*  ITYPRO =  4 : History := HST: History Burnup   ::
*  ITYPRO =  5 : History := HST: Map              ::
*  ITYPRO =  6 : History := HST: History Map      ::
*  ITYPRO = -1 :         := HST: History          ::
*  ITYPRO = -3 : Burnup  := HST: History          ::
*  ITYPRO = -4 : Burnup  := HST: Burnup History   ::
*  ITYPRO = -5 : Map     := HST: Map  History     ::
*----
      IF(NENTRY .EQ. 1) THEN
        IF(IKEVO .NE. 0 .OR. IKMAP .NE. 0) CALL XABORT(NAMSBR//
     >  ': A single burnup or map structure forbidden.')
        IF(IKHST .EQ. 1) THEN
          ITYPRO=2
          IF(JENTRY(1) .EQ. 2) THEN
            ITYPRO=-1
          ENDIF
        ELSE
          IKHST=1
          ITYPRO=1
          SIGENT(IKHST)='L_HISTORY   '
          READ(SIGENT(IKHST),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
          CALL LCMPUT(KENTRY(IKHST),'SIGNATURE',NTC,3,NAMTMP)
        ENDIF
      ELSE
        IF(IKHST .EQ. 2) THEN
          IF(IKMAP .EQ. 1) THEN
            ITYPRO = -5
          ELSE IF(IKEVO .EQ. 1) THEN
            ITYPRO=-4
          ELSE 
            ITYPRO=-3
            IKEVO=1 
            SIGENT(IKEVO)='L_BURNUP    '
            READ(SIGENT(IKEVO),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
            CALL LCMPUT(KENTRY(IKEVO),'SIGNATURE',NTC,3,NAMTMP)
          ENDIF
        ELSE IF(IKEVO.EQ.2) THEN
          IF(IKHST .EQ. 1) THEN
            ITYPRO=4
          ELSE
            ITYPRO=3
            IKHST=1 
            SIGENT(IKHST)='L_HISTORY   '
            READ(SIGENT(IKHST),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
            CALL LCMPUT(KENTRY(IKHST),'SIGNATURE',NTC,3,NAMTMP)
          ENDIF
        ELSE IF(IKMAP.EQ.2) THEN
          IF(IKHST .EQ. 1) THEN
            ITYPRO=6
          ELSE
            ITYPRO=5
            IKHST=1 
            SIGENT(IKHST)='L_HISTORY   '
            READ(SIGENT(IKHST),'(3A4)') (NAMTMP(ITC),ITC=1,NTC)
            CALL LCMPUT(KENTRY(IKHST),'SIGNATURE',NTC,3,NAMTMP)
          ENDIF
        ELSE
          CALL XABORT(NAMSBR//
     >    ': A read-only burnup or map structure required.')
        ENDIF
      ENDIF
      IF(IKHST .NE. 0) IPHST=KENTRY(IKHST)
      IF(IKEVO .NE. 0) IPEVO=KENTRY(IKEVO)
      IF(IKMAP .NE. 0) IPMAP=KENTRY(IKMAP)
*----
*  Get elements of HISTORY STATE-VECTOR
*----
      MAXG  =ISTATH( 1)
      MAXL  =ISTATH( 2)
      NBUNH =ISTATH( 3)
      NCHAH =ISTATH( 4)
      ITSOLH=ISTATH( 6)
      ITBURH=ISTATH( 7)
      MAXIH =ISTATH( 8)
      NREGH =ISTATH(10)
      IF(IDEB .EQ. 1) THEN
        WRITE(IOUT,7000) (ISTATH(II),II=1,8),ISTATH(10)
      ENDIF
*----
*  Get elements of BURNUP STATE-VECTOR
*----
      ITSOLB=ISTATB(1)
      ITBURB=ISTATB(2)
      NBBTS =ISTATB(3)
      MAXIB =ISTATB(4)
      IF(IDEB .EQ. 1) THEN
        WRITE(IOUT,7001) (ISTATB(II),II=1,6)
      ENDIF
      IF(ITYPRO .EQ. 3 .OR. ITYPRO .EQ. 4) THEN
        ITSOLH=ITSOLB
        ITBURH=ITBURB
        IF(MAXIH .NE. 0 .AND. MAXIH .NE. MAXIB) CALL XABORT(NAMSBR//
     >  ': Different number of isotopes in history and burnup')
        MAXIH=MAXIB
      ELSE IF(ITYPRO .EQ. -3 .OR. ITYPRO .EQ. -4) THEN
        ITSOLB=ITSOLH
        ITBURB=ITBURH
        IF(MAXIB .NE. 0 .AND. MAXIB .NE. MAXIH) CALL XABORT(NAMSBR//
     >  ': Different number of isotopes in history and burnup')
        MAXIB=MAXIH
      ENDIF 
*----
*  Get elements of MAP STATE-VECTOR
*  and verify consistency with history information
*----
      NBUNM =ISTATM(1)
      NCHAM =ISTATM(2)
      IF(NBUNM .NE. 0) THEN
        IF(NBUNH .EQ. 0) THEN
          NBUNH=NBUNM
        ELSE IF(NBUNH .NE. NBUNM) THEN 
          CALL XABORT(NAMSBR//': Different number of bundles in'//
     >                        ' MAP and HISTORY structures')
        ENDIF
      ENDIF
      IF(NCHAM .NE. 0) THEN
        IF(NCHAH .EQ. 0) THEN
          NCHAH=NCHAM
        ELSE IF(NCHAH .NE. NCHAM) THEN 
          CALL XABORT(NAMSBR//': Different number of channels in'//
     >                        ' MAP and HISTORY structures')
        ENDIF
      ENDIF
*----
*  Test compatibility of HISTORY, BURNUP and MAP data structures.
*---- 
      IF(ITYPRO .EQ.  4 .OR. ITYPRO .EQ. -4) THEN
        IF(ITSOLB .NE. ITSOLH .OR. 
     >     ITBURB .NE. ITBURH .OR.
     >     MAXIB  .NE. MAXIH  ) CALL XABORT(NAMSBR//
     >    ': HISTORY and BURNUP parameters incompatible')
      ELSE IF(ITYPRO .EQ.  6) THEN
        IF(NBUNM  .NE. NBUNH .OR. 
     >     NCHAM  .NE. NCHAH      ) CALL XABORT(NAMSBR//
     >    ': HISTORY and MAP parameters incompatible')
      ENDIF
*----
*  Get EDIT level and dimensioning parameters for history structure
*  and test their validity
*----
      IPRINT=1
      NGLO  =MAXG
      NLOC  =MAXL
      NBUN  =NBUNH
      NCHA  =NCHAH
      CALL HSTGDM(IPRINT,NGLO  ,NLOC  ,NCHA  ,NBUN  ,
     >            BUNLEN,ITYRED,CARRED)
*----
*  Test dimensioning parameters for coherence
*  with already defined parameters
*----
      MAXG=MAX(MAXG,NGLO)
      MAXL=MAX(MAXL,NLOC)
      IF(NBUN .LE. 0 ) CALL XABORT(NAMSBR// 
     >': Number of bundles must be larger than 0')
      IF(NCHA .LE. 0 ) CALL XABORT(NAMSBR// 
     >': Number of channels must be larger than 0')
      IF(NBUNH .GT. 0 .AND. NBUN .NE. NBUNH) CALL XABORT(NAMSBR// 
     >': Number of bundles on input'//
     >' different from HISTORY, MAP or BURNUP structures')
      NBUNH=MAX(NBUN,NBUNH)
      IF(NCHAH .GT. 0 .AND. NCHA .NE. NCHAH) CALL XABORT(NAMSBR// 
     >': Number of channels on input'//
     >' different from HISTORY, MAP or BURNUP structures')
      NCHAH=MAX(NCHA,NCHAH)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6010) NGLO,NLOC,NCHA,NBUN
      ENDIF
*----
* Allocate memory for global and local parameters
*---- 
      ALLOCATE(NAMG(3*(MAXG+1)),PARAG(MAXG+1),NAML(3*(MAXL+1)),
     > PARAL((MAXL+1)*2))
      CALL XDISET(NAMG,3*(MAXG+1),IBLANK)
      CALL XDRSET(PARAG,(MAXG+1),0.0)
      CALL XDISET(NAML,3*(MAXL+1),IBLANK)
      IF(ISTATH(1) .GT. 0) THEN
        CALL LCMGET(IPHST,'NAMEGLOBAL  ',NAMG(4))
        CALL LCMGET(IPHST,'PARAMGLOBAL ',PARAG(2))
        IF(IDEB .GE. 1) THEN 
          WRITE(IOUT,'(A18,2I10)') 'Initial NAMEGLOBAL',MAXG,ISTATH(1) 
          WRITE(IOUT,'(6(3A4,2X))') (NAMG(3+II),II=1,3*MAXG) 
        ENDIF
      ENDIF
      IF(ISTATH(2) .GT. 0) THEN
        CALL LCMGET(IPHST,'NAMELOCAL   ',NAML(4))
        IF(IDEB .GE. 1) THEN 
          WRITE(IOUT,'(A18,2I10)') 'Initial NAMELOCAL ',MAXL,ISTATH(2) 
          WRITE(IOUT,'(6(3A4,2X))') (NAML(3+II),II=1,3*MAXL) 
        ENDIF
      ENDIF
      IF(NCHAH .LT. 1 .OR. NBUNH .LT. 1 ) CALL XABORT(NAMSBR//
     >': Both the number of channels and bundles must be > 0')
*----
*  Allocate memory for core description
*----
      NCELL=NCHAH*NBUNH
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6011) NCELL,NCHAH,MAXIH
      ENDIF
      ALLOCATE(IDCELL(NCELL),IDFUEL(NCELL),IREFUS(NCHAH),REFUT(NCHAH))
      CALL XDISET(IDCELL,NCELL,0)
      CALL XDISET(IDFUEL,NCELL,0)
      IF(ISTATH( 3)*ISTATH( 4) .GT. 0) THEN
        CALL LCMGET(IPHST,'CELLID      ',IDCELL)
        CALL LCMGET(IPHST,'FUELID      ',IDFUEL)
      ENDIF
      CALL XDISET(IREFUS,NCHAH,0)
      CALL XDRSET(REFUT,NCHAH,0.0)
      ALLOCATE(DENI(MAXIH+1))
      NBFUEL=0
*----
*  Allocate memory for MAP power
*----
      ALLOCATE(POWR(NCELL),BURN(NCELL))
      CALL XDRSET(POWR,NCELL,0.0)
      CALL XDRSET(BURN,NCELL,0.0)
      IF(ITYPRO .EQ. 5 .OR. ITYPRO .EQ. 6) THEN
*----
*  Read information from MAP data structure
*  and update history using this information
*---- 
        CALL HSTUHM(IPHST, IPMAP, IPRINT, MAXL, NCHAH ,NBUNH, MAXIH,
     >              POWR,BURN,IREFUS,
     >              REFUT,BUNLEN,IDCELL,IDFUEL,PARAL,DENI)        
*----
* Update Map with History
*---- 
      ELSE IF(ITYPRO .EQ. -5) THEN
        CALL HSTUMH(IPMAP, IPHST, IPRINT,NCHAH ,NBUNH, IDCELL, BURN)
      ENDIF
*----
*  Release memory for MAP power
*----
      DEALLOCATE(BURN,POWR)
*----
*  Read or write remaining information on input
*  Also extract information from history if required
*----
      CALL HSTGET(IPHST ,IPRINT,MAXG  ,MAXL  ,NCHAH ,NBUNH ,
     >            ITYPRO,ITYRED,CARRED,IUPDC ,IUPDB ,
     >            NAMG ,PARAG,NAML ,
     >            PARAL,IDCELL,IDFUEL)
      IF(ITYPRO .GT. 0) THEN
        IF(MAXG .GT. 0) THEN
          CALL LCMPUT(IPHST,'NAMEGLOBAL  ',3*MAXG,3,NAMG(4))
          CALL LCMPUT(IPHST,'PARAMGLOBAL ',  MAXG,2,PARAG(2))
          IF(IDEB .GE. 1) THEN 
            WRITE(IOUT,'(A18,2I10)') 'Final NAMEGLOBAL  ',MAXG,ISTATH(1)
            WRITE(IOUT,'(6(3A4,2X))') (NAMG(3+II),II=1,3*MAXG) 
          ENDIF
        ENDIF
        IF(MAXL .GT. 0) THEN
          CALL LCMPUT(IPHST,'NAMELOCAL   ',3*MAXL,3,NAML(4))
          IF(IDEB .GE. 1) THEN 
            WRITE(IOUT,'(A18,2I10)') 'Final NAMELOCAL   ',MAXL,ISTATH(2)
            WRITE(IOUT,'(6(3A4,2X))') (NAML(3+II),II=1,3*MAXL) 
          ENDIF
        ENDIF
        IF(NCELL .GT. 0) THEN
          CALL LCMPUT(IPHST,'CELLID      ',NCELL,1,IDCELL)
          CALL LCMPUT(IPHST,'FUELID      ',NCELL,1,IDFUEL)
        ENDIF
      ENDIF
*----
*  If channel and bundle specified
*  Update information on HISTORY or BURNUP structures
*----
      IF(IUPDC .GT. 0 .AND. IUPDB .GT. 0) THEN
*----
*  Allocate memory for isotopes and burnup
*----
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6000) NAMSBR,IUPDC,IUPDB 
        ENDIF
        IF(ITYPRO .EQ. 3 .OR. ITYPRO .EQ. 4) THEN
*----
*  Update HISTORY information from BURNUP data for 
*  channel IUPDC, bundle IUPDB.
*---- 
          IF(IPRINT .GE. 10) THEN
            WRITE(IOUT,6001) 
          ENDIF
          CALL HSTUHB(IPHST ,IPEVO ,IPRINT,MAXIH ,NBBTS ,
     >                NCHAH ,NBUNH ,IUPDC ,IUPDB ,
     >                IDCELL,IDFUEL,
     >                DENI ,MAXL, PARAL) 
        ELSE IF(ITYPRO .EQ. -3 .OR. ITYPRO .EQ. -4) THEN
*----
*  Update BURNUP information from HISTORY data for 
*  channel IUPDC, bundle IUPDB.
*---- 
          IF(IPRINT .GE. 10) THEN
            WRITE(IOUT,6002) 
          ENDIF
          CALL HSTUBH(IPEVO ,IPHST ,IPRINT,MAXIH ,NBBTS ,
     >                NCHAH ,NBUNH ,IUPDC ,IUPDB ,
     >                IDCELL,IDFUEL,DENI)
        ENDIF 
      ENDIF
      DEALLOCATE(DENI,REFUT,IREFUS,IDFUEL,IDCELL,PARAL,NAML,PARAG,NAMG)
      IF(ITYPRO .GT. 0) THEN
*----
*  Saving updated HISTORY state vector
*----
        CALL LCMPUT(IPHST,'BUNDLELENGTH',1,2,BUNLEN)
        CALL XDISET(ISTATH,NSTATE,0)
        ISTATH( 1) = MAXG  
        ISTATH( 2) = MAXL  
        ISTATH( 3) = NBUNH 
        ISTATH( 4) = NCHAH 
        ISTATH( 5) = 0 
        ISTATH( 6) = ITSOLH
        ISTATH( 7) = ITBURH
        ISTATH( 8) = MAXIH 
        ISTATH(10) = NREGH
        IF(IPRINT .EQ. 10) THEN
          WRITE(IOUT,7010) (ISTATH(II),II=1,8),ISTATH(10)
        ENDIF
        CALL LCMPUT(IPHST,'STATE-VECTOR',NSTATE,1,ISTATH)
      ELSE IF(ITYPRO .EQ. -3 .OR. ITYPRO .EQ. -4 )  THEN
*----
*  Set burnup parameters to default values
*  See subroutine EVO.f
*----
        REVOL(1)=1.0E-5 
        REVOL(2)=1.0E-4 
        REVOL(3)=80.0   
        REVOL(4)=1.0E-4 
        REVOL(5)=0.0 
        CALL LCMPUT(IPEVO,'EVOLUTION-R ',5,2,REVOL)
*----
*  Saving updated BURNUP state vector
*----
        CALL XDISET(ISTATB,NSTATE,0)
        ISTATB( 1) = ITSOLB 
        ISTATB( 2) = ITBURB 
        IF(ISTATB( 1) .EQ. 0) ISTATB( 1) = 2
        IF(ISTATB( 2) .EQ. 0) ISTATB( 2) = 2
        ISTATB( 3) = 1 
        ISTATB( 4) = MAXIH 
        ISTATB( 8) = NCHA*NBUN
        IF(IPRINT .GT. 1) THEN
          WRITE(IOUT,7011) (ISTATB(II),II=1,8)
        ENDIF
        CALL LCMPUT(IPEVO,'STATE-VECTOR',NSTOLD,1,ISTATB)
      ENDIF
*----
*  Module execution completed
*----
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(' ***** OUTPUT FROM ',A6/
     >' Processing: Channel ',I10,5X,'Bundle ',I10) 
 6001 FORMAT(' Updating HISTORY from BURNUP')
 6002 FORMAT(' Updating BURNUP from HISTORY')
 6010 FORMAT(' ***** General dimensioning '/
     >       10X,'NGLO  =',I10,5X,'NLOC   =',I5/
     >       10X,'NCHA  =',I10,5X,'NBUN   =',I5)
 6011 FORMAT(10X,'NCELL =',I10,5X,'NCHAH  =',I5/
     >       10X,'MAXIH =',I10)
 7000 FORMAT(' Initial contents of HISTORY state vector'/
     >5X,'MAXG  = ',I5,5X,'MAXL  = ',I5,5X,'NBUNH = ',I5,/
     >5X,'NCHAH = ',I5,5X,'      = ',I5,5X,'ITSOLH= ',I5,/
     >5X,'ITBURH= ',I5,5X,'MAXIH = ',I5,5X,'NREGH = ',I5) 
 7001 FORMAT(' Initial contents of BURNUP state vector'/
     >5X,'ITSOL = ',I5,5X,'ITBUR = ',I5,5X,'NBBTS = ',I5,/
     >5X,'MAXI  = ',I5,5X,'NGRP  = ',I5,5X,'NREG  = ',I5)
 7010 FORMAT(' Final contents of HISTORY state vector'/
     >5X,'MAXG  = ',I5,5X,'MAXL  = ',I5,5X,'NBUNH = ',I5,/
     >5X,'NCHAH = ',I5,5X,'      = ',I5,5X,'ITSOLH= ',I5,/
     >5X,'ITBURH= ',I5,5X,'MAXIH = ',I5,5X,'NREGH = ',I5) 
 7011 FORMAT(' Final contents of BURNUP state vector'/
     >5X,'ITSOL = ',I5,5X,'ITBUR = ',I5,5X,'NBBTS = ',I5,/
     >5X,'MAXI  = ',I5,5X,'      = ',I5,5X,'      = ',I5,/
     >5X,'      = ',I5,5X,'NBMIX = ',I5) 
      END  
