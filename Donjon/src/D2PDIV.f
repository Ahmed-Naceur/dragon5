*DECK D2PDIV
      SUBROUTINE D2PDIV(  IPDAT, IPSAP , IPRINT,    NGP,   NBU,  NVAR,
     >                     GRID,  NPAR ,   NREA,   NISO,  NMAC,  NMIL,
     >                     NANI, NADRX , STAIDX,  STATE, STAVAR,  NSF,
     >                     LABS,   SCAT,   LADF                      )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the DIVERS directory of an elementary calculation and store
* additional XS recovered directly from IPSAP
* WARNING: the GET_DIVERS_INFO subroutine cannot recover DIVERS
* information in the case where cross sections are interpolated by
* the SCR: module
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of the INFO data block
* IPSAP   address of the saphyb object
* IPRINT  control the printing on screen
* NGP     number of energy groups
* NBU     number of burnup point in IPSAP
* NVAR    number of state parameters in INFO data block
* GRID    type of gridding for branches (0 = default, 1 = Saphyb
*         branching etc )
* NPAR    number of state parameters in saphyb (including FLUE and
*         TIME)
* NREA    number of reactions in IPSAP
* NISO    number of isotopoes in IPSAP
* NMAC    number of macros in IPSAP
* NMIL    number of mixtrures  in IPSAP
* NANI    number of anisotropy
* STAIDX  index of state variables
* STATE   state variables of current branch calculation
* STAVAR  state variables in INFO data block
* NSF     nummber of surface in IPSAP
* LABS    information for absorption reconstruction
* SCAT    information for scattering XS reconstruction
* LADF    flag for ADF reconstrcution
*
*Parameters: 
* NADRX   
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPSAP
      INTEGER NPAR,NMIL,GRID,NVAR,NBU,NSF,NREA,NISO,NADRX
      INTEGER NGP,IPRINT,NMAC,NANI,STAIDX (NVAR)
      REAL STATE(NVAR)
      CHARACTER(LEN=12) STAVAR(NVAR)
      LOGICAL LABS(3),SCAT,LADF
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPTH,KPTH
      ! LOOP INDEX
      INTEGER i, It, Ib,PK
      ! LOOP INDEX OF : PARAMETERS (ISV=1..NPAR), STATES (INP=1..NVAR)
      INTEGER ISV,INP
      ! DIMENSION OF ARBVAL
      INTEGER DIMARB
      ! NUMBER OF ELEMENTARY CALCULATIONS
      INTEGER NCALS
      ! TYPE OF DATA RECOVERED FROM GANLIB SUBROUTINES
      INTEGER ITYLCM
      ! NUMBER OF VALUES IN IDVAL ET VALDIV
      INTEGER NVDIV
      ! ORDER NUMBERS OF FLUE PARAMETERS IN SAPHYB
      INTEGER :: FLUE_ID = 0
      ! ORDER NUMBERS OF TIME PARAMETERS IN SAPHYB
      INTEGER :: TIME_ID = 0
      ! CF : APOLLO2 : NOTICE INFORMATIQUE DE LA VERSION 2.8-1
      INTEGER MUPLET(NPAR)
      ! VECTOR OF : RANK ORDER OF STATE PARAMETERS, NUMBER OF VALUES
      ! FOR EACH STATE PARAMETERS
      INTEGER  RANK_ORDER(NPAR), NVALUE(NPAR)
      REAL B2
      CHARACTER*3 :: ADF_T = 'DRA'
      ! NAME OF DIRECTORIES IN SAPHYB : ELEMENTARY CALCULATION,
      ! CONTROL ROD
      CHARACTER(LEN=12) CALDIR,BARRDIR
      ! NAME OF STATE VARIABLES IN SAPHYB
      CHARACTER(LEN=12) PKNAM(6)
      ! STATE VARIABLES IN SAPHYB
      CHARACTER(LEN=12) PKEY(NPAR)
      LOGICAL LFLAG(6)

      ! CF : APOLLO2 : NOTICE INFORMATIQUE DE LA VERSION 2.8-1
      ! VALUES OF : VALDIV = (KEFF, KINF,B2),  CONTROL ROD KEFF, KINF,B
      INTEGER, ALLOCATABLE, DIMENSION(:) ::  DEBARB,ARBVAL
      REAL, ALLOCATABLE, DIMENSION(:) :: VALDIV,BARR_VAL
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: IDVAL

      ! RECOVER INFOMATION FROM INFO DATA BLOCK AND SAPHYB OBJECT

      ! MOVING INTO INFO DATA BLOCK
      CALL LCMSIX (IPSAP,' ',0)

      CALL LCMSIX (IPSAP,'paramdescrip',1)
      CALL LCMGTC(IPSAP,'PARKEY',4,NPAR,PKEY)
      CALL LCMGET (IPSAP,'NVALUE',NVALUE)

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      IF (LADF) CALL LCMGTC(IPDAT,'ADF_TYPE',3,1,ADF_T)
      PKEY (1:NPAR) (5:12) = "        "
      DO PK=1, 6
        IPTH=LCMGID(IPDAT,'PKEY_INFO')
        KPTH=LCMDIL(IPTH,PK)
        CALL LCMGET(KPTH,'LFLAG',LFLAG(PK))
        IF (PK == 1 .OR. PK==6)THEN
         CALL LCMGTC(KPTH,'NAME',12,1,PKNAM(PK))
        ELSE
         IF(LFLAG(PK)) CALL LCMGTC(KPTH,'NAME',12,1,PKNAM(PK))
        ENDIF
      ENDDO
      ! LOOP TO STORE THE INDEX OF FLUE AND
      ! LINK THE FLUE AND TIME VARIABLES INDEX TO BURN VARIABLE INDEX
      DO It=1, NPAR
        IF(PKEY(It)=="TIME") TIME_ID=It
        IF(PKEY(It)=="FLUE") FLUE_ID=It
      ENDDO
      ! LOOP OVER NUMBER OF STATE PARAMETERS IN SAPHYB
      DO ISV=1, NPAR
        ! LOOP OVER NUMBER OF STATE PARAMETERS IN INFO DATA BLOCK
        DO INP=1, NVAR
         ! IF NAME OF STATE VARIABLE IN INFO AND SAPHYB ARE EQUAL
         IF(PKEY(ISV)==STAVAR(INP)) THEN
          ! SPECIAL CASE FOR BARR parameters
          IF(PKEY(ISV)==PKNAM(1)) THEN
           !SPECIAL CASE FOR CONTROL ROD
           ALLOCATE (BARR_VAL(NVALUE(ISV)))
           WRITE(BARRDIR,'("pval", I8)') ISV
           ! NAME OF DIRECTORY IN SAPHYB  CONTAINING CONTROL ROD VALUES
           IF(LFLAG(1)) THEN
           ! RECOVER CONTROL ROD VALUES
            CALL LCMSIX (IPSAP,' ',0)
            CALL LCMSIX (IPSAP,'paramvaleurs',1)
            CALL LCMGET(IPSAP,BARRDIR,BARR_VAL)

           ! LOOP OVER POSSIBLE VALUES OF CONTROL ROD IN SAPHYB
            DO Ib=1, NVALUE(ISV)
             IF(STATE(INP)==BARR_VAL(Ib)) THEN
              ! STORE THE ORDER NUMBERS OF CURRENT CONTROL VALUES
              ! CORRESPONDING TO THE BRANCH CALCULATED
              RANK_ORDER(ISV)=Ib
             ENDIF
            ENDDO
           ENDIF
           DEALLOCATE (BARR_VAL)

          ! SPECIAL CASE WITH DEFAULT VALUES FOR STATE VARIABLES
          ! (OTHER THAN BARR)
          ELSE IF(GRID==0) THEN
           ! TREATEMENT OF THE MID VALUE OF THE GRID
           IF(STAIDX(INP)==2) THEN
            ! ONLY DMOD,TCOM AND CBOR ARE  AFFECTED BY THE DEFAULT
            ! GRIDDING
            IF((PKEY(ISV)==PKNAM(2))) THEN
             RANK_ORDER(ISV)=NINT (NVALUE(ISV)/2.0)
            ELSE IF((PKEY(ISV)==PKNAM(4)))THEN
             RANK_ORDER(ISV)=NINT (NVALUE(ISV)/2.0)
            ELSE IF((PKEY(ISV)==PKNAM(3)))THEN
             RANK_ORDER(ISV)=NINT (NVALUE(ISV)/2.0)
            ELSE
             RANK_ORDER(ISV)=STAIDX(INP)
            ENDIF
           ! TREATEMENT OF THE LAST VALUE OF THE GRID
           ELSE  IF(STAIDX(INP)==3) THEN
            ! ONLY DMOD,TCOM AND CBOR ARE  AFFECTED BY THE DEFAULT
            ! GRIDDING
            IF((PKEY(ISV)==PKNAM(2))) THEN
             RANK_ORDER(ISV)=NVALUE(ISV)
            ELSE IF((PKEY(ISV)==PKNAM(4)))THEN
             RANK_ORDER(ISV)=NVALUE(ISV)
            ELSE IF((PKEY(ISV)==PKNAM(3)))THEN
             RANK_ORDER(ISV)=NVALUE(ISV)
            ELSE
             RANK_ORDER(ISV)=STAIDX(INP)
            ENDIF
             ! ONLY DMOD,TCOM AND CBOR ARE  AFFECTED BY THE DEFAULT
             ! GRIDDING
           ELSE ! THE FIRST VALUE IS UNCHANGED BY SET_DEFAULT_VALUE
             RANK_ORDER(ISV)=STAIDX(INP)
           ENDIF
         ! IF WE KEEP THE INITIAL STATE VARIABLE GRID OF SAPHYB
          ELSE
           RANK_ORDER(ISV)=STAIDX(INP)
          ENDIF
          !TREATMENT OF FLUE AND TIME VARIABLES
          IF(PKEY(ISV)==PKNAM(6)) THEN
           IF(FLUE_ID>0) RANK_ORDER(FLUE_ID)=RANK_ORDER(ISV)
           IF(TIME_ID>0) RANK_ORDER(TIME_ID)=RANK_ORDER(ISV)
          ENDIF
         ENDIF
        ENDDO
      ENDDO

      ! RECOVER INFORMATION FROM SAPHYB
      CALL LCMSIX (IPSAP,' ',0)
      CALL LCMSIX (IPSAP,'paramarbre',1)
      CALL LCMLEN (IPSAP,'ARBVAL',DIMARB,ITYLCM)
      ALLOCATE (ARBVAL(DIMARB),DEBARB(DIMARB+1))
      CALL LCMGET (IPSAP,'NCALS',NCALS)
      CALL LCMGET (IPSAP,'ARBVAL',ARBVAL)
      CALL LCMGET (IPSAP,'DEBARB',DEBARB)
      ! PROCEDURE TO RECOVER THE NUMBER OF THE ELEMENTARY CALCULATION
      ! CORREPSONDING TO THE CURRENT BRANCH
      ! CF APOLLO2 : NOTICE INFORMATIQUE DE LA VERSION 2.8-1
      II=1
      DO 30 IPAR=1,NPAR
        MUPLET(IPAR) =RANK_ORDER(IPAR)
        DO 10 I=DEBARB(II),DEBARB(II+1)-1
         IF(MUPLET(IPAR).LE.ARBVAL(I))THEN
          IF(MUPLET(IPAR).EQ.ARBVAL(I))THEN
           II=I
           GO TO 30
          ELSE
           GO TO 20
          ENDIF
         ENDIF
10      CONTINUE
20      ICAL=0
       WRITE(6,*) " MUPLET : ", MUPLET
       CALL XABORT ("@D2PDIV: ELEMENTARY CALCULATION UNKNOWN")
       RETURN
30    CONTINUE
      ! END OF APPOLO2 PROCEDURE

      ICAL=DEBARB(II+1) ! number of the elementary calculation

      ! MOVING IN THE ELEMENTARY CALCULATION AND RECONVER THE B2, KEFF
      ! AND KINF DATA
      WRITE(CALDIR,'("calc", I8)') ICAL
      CALL LCMSIX (IPSAP,' ',0)
      CALL LCMSIX (IPSAP,CALDIR,1)
      CALL LCMSIX(IPSAP,'divers',1)
      CALL LCMGET(IPSAP,'NVDIV',NVDIV)

      ALLOCATE(IDVAL(NVDIV),VALDIV(NVDIV))
      CALL LCMGTC(IPSAP,'IDVAL',4,NVDIV,IDVAL)
      CALL LCMGET(IPSAP,'VALDIV',VALDIV)


      ! STORE RESULTS (IF CORRESPONDING DATA IS AVAILABLE) INTO INFO
      ! data block at :
      ! INFO/BRANCH_INFO/KEFF
      ! INFO/BRANCH_INFO/B2
      ! INFO/BRANCH_INFO/KINF

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      IF(STAIDX(NVAR)==1) THEN
        IPTH=LCMLID(IPDAT,'DIVERS',NBU)
      ELSE
        IPTH=LCMGID(IPDAT,'DIVERS')
      ENDIF
      KPTH=LCMDIL(IPTH,STAIDX(NVAR))

      IF(IPRINT>1) THEN
        WRITE(6,*)
        WRITE(6,*) "****          DIVERS INFORMATION           ****"
      ENDIF
      DO Idiv=1, NVDIV
        IF(IDVAL(Idiv)=="KEFF") THEN
         CALL LCMPUT(KPTH,'KEFF',1,2,VALDIV(Idiv))
         IF(IPRINT>1) WRITE(6,*)"KEFF                    :",VALDIV(Idiv)
        ENDIF
        IF(IDVAL(Idiv)=="KINF") THEN
         CALL LCMPUT(KPTH,'KINF',1,2,VALDIV(Idiv))
         IF(IPRINT>1) WRITE(6,*)"KINF                    :",VALDIV(Idiv)
        ENDIF
        IF(IDVAL(Idiv)=="B2") THEN
         CALL LCMPUT(KPTH,'B2',1,2,VALDIV(Idiv))
         B2=VALDIV(Idiv)
         IF(IPRINT>1) WRITE(6,*)"B2                      :",VALDIV(Idiv)
        ENDIF
      ENDDO
      ! TEMPORARY SUBROUTINE WAITING FOR FURTHER DEVELOMENTS TO RECOVER
      ! ADDITIONAL INFORMATION
      CALL D2PXSA(IPDAT,IPSAP,ICAL,IPRINT,NGP,NREA,NISO,NMAC,NMIL,
     1   NANI,NVAR,NADRX,STAIDX,B2,ADF_T,NSF,LABS,SCAT,LADF)
      DEALLOCATE (ARBVAL,DEBARB,VALDIV,IDVAL)            ! FREE MEMORY
      END
