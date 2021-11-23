      SUBROUTINE D2PPRC ( IPDAT,IPPRC, HEQUI, HMASL, ISOTVAL, ISOTOPT,
     >                     LMEM,IPRINT,MIXDIR,JOBOPT                 )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build a procedure file for the interpolation of cross sections
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   adress of info data block
* HEQUI   name of the equivalence record in the saphyb|MCO object
* HMASL   name of heavy metal density record in the saphyb|MCO object
* ISOTVAL concentration of particularized isotopes
* ISOTOPT otpion for paticularised isotopes
*
*Parameters: 
* IPPRC    
* LMEM     
* IPRINT   
* MIXDIR   
* JOBOPT   
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPTH,KPTH
      INTEGER IPPRC,PK,IPRINT
      CHARACTER*4 HEQUI,HMASL
      CHARACTER*1 ISOTOPT,JOBOPT(14)
      CHARACTER*12,ISOTOPES(8)
      REAL ISOTVAL
      LOGICAL LMEM,LFLAG(6)
*----
*  LOCAL VARIABLES
*----
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PKEY
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: OTHPK
      INTEGER, ALLOCATABLE, DIMENSION(:) :: OTHTYP
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: OTHVAC
      REAL, ALLOCATABLE, DIMENSION(:) :: OTHVAR
      CHARACTER*12 PKNAM(6),MIXDIR
      INTEGER STAVEC(40),NVAR,ITYP,NOTH
      INTEGER :: NTOT = 0
      INTEGER :: NPKEY = 0
      INTEGER :: ORDER(6) = -1
      CHARACTER*6 :: NAMSAP='XSLIB'
      CHARACTER*4,DIMENSION(6) :: REFNAM
      DATA REFNAM/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/

      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)'!    Auto Generation of input file for D2P      *'
      WRITE(IPPRC,*)'! - Recovering of information from D2P PHASE 1  *'
      WRITE(IPPRC,*)'! - call to the interpolation module(SCR|NCR)   *'
      WRITE(IPPRC,*)'! - call of D2P for PHASE 2 and 3               *'
      WRITE(IPPRC,*)'! Author(s) : J. TAFOREAU (2016)                *'
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)

      WRITE(IPPRC,*)" SEQ_ASCII GENPMAXS :: FILE 'GENPMAXS.inp' ; "
      WRITE(IPPRC,*)" SEQ_ASCII HELIOS :: FILE 'HELIOS.dra' ; "
      WRITE(IPPRC,*)" XSM_FILE XSLIB :: FILE 'XSLIB' ;  "
      WRITE(IPPRC,*)" XSM_FILE D2PINFO :: FILE 'Info.xsm' ; "
      WRITE(IPPRC,*)" LINKED_LIST INFO ; "
      IF (LMEM) THEN
       WRITE(IPPRC,*)'LINKED_LIST XSL ; '
       NAMSAP='XSL'
      ENDIF

      WRITE(IPPRC,*)'LINKED_LIST Micro ; '
      WRITE(IPPRC,*)'MODULE END: D2P: SCR: NCR: GREP: DELETE: UTL: ;'
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)'!  STEP 0 :   Initializing state parameters     *'
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)

      CALL LCMGET(IPDAT,'STATE-VECTOR',STAVEC)
      NVAR=STAVEC(2)
      ITYP=STAVEC(18)
      NOTH=STAVEC(20)
      ALLOCATE(PKEY(NVAR))
      ALLOCATE(OTHPK(NOTH),OTHTYP(NOTH),OTHVAC(NOTH),OTHVAR(NOTH))
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMGTC(IPDAT,'ISOTOPES',12,4,ISOTOPES)
      CALL LCMGTC(IPDAT,'STATE_VAR',12,NVAR,PKEY)
      CALL LCMGET(IPDAT,'NTOT',NTOT)
      IF (NOTH>0) THEN
       CALL LCMGTC(IPDAT,'OTHPK',12,NOTH,OTHPK)
       CALL LCMGET(IPDAT,'OTHTYP',OTHTYP)
       CALL LCMGTC(IPDAT,'OTHVAC',12,NOTH,OTHVAC)
       CALL LCMGET(IPDAT,'OTHVAR',OTHVAR)
      ENDIF

      DO PK=1, 6
        IPTH=LCMGID(IPDAT,'PKEY_INFO')
        KPTH=LCMDIL(IPTH,PK)
        CALL LCMGET(KPTH,'LFLAG',LFLAG(PK))
        IF(LFLAG(PK)) THEN
         NPKEY=NPKEY+1
         CALL LCMGTC(KPTH,'NAME',12,1,PKNAM(PK))
         CALL LCMGTC(KPTH,'NAME',12,1,PKNAM(PK))
         WRITE(IPPRC,*) 'STRING  ',
     >    REFNAM(PK),' := "',TRIM(PKNAM(PK)),'" ; '
         WRITE(IPPRC,*) 'REAL    ',REFNAM(PK),'_VAL     ; '
         DO I=1,NVAR
          IF (PKNAM(PK).EQ.PKEY(I)) THEN
          ORDER(PK)=I
          ENDIF
         ENDDO
        ENDIF
      ENDDO

      IF (NTOT.NE.(NOTH+NPKEY)) THEN
       WRITE(6,*) "@D2PPROC: INCONSISTENT D2P INPUT DATA WITH",
     > "XS LIBRARY"
       WRITE(6,*) "D2P INPUT DATA     : "
       WRITE(6,*) "    STATE VARIABLE : ", NPKEY
       WRITE(6,*) "    OTHER VARIABLE : ", NOTH
       WRITE(6,*) "D2P  TOTAL         = ", NPKEY+NOTH
       WRITE(6,*) "XS LIBRARY CONTENT = ", NTOT
       CALL XABORT ("=>PLEASE USE THE D2P CARD 'PKEY'AND/OR 'OTHER'")
      ENDIF

      IF (NPKEY .EQ. 0) THEN
       WRITE(6,*) "@D2PPROC : NUMBER OF STATE VARIABLES IS ZERO"
       CALL XABORT ("=> PLEASE CHECK THE D2P DATA INPUT ")
      ENDIF
      WRITE(IPPRC,*)'INFO := D2PINFO ; '
      IF (LMEM) WRITE(IPPRC,*)'XSL := XSLIB ;'
      WRITE(IPPRC,*)'INTEGER NVAR := ',NPKEY,' ; '
      WRITE(IPPRC,*)'INTEGER STOP REWIND ITER := 0 0 0 ; '

      WRITE(IPPRC,*)
      WRITE(IPPRC,*)'WHILE STOP 1 <> DO'

      WRITE(IPPRC,*)
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)'!  STEP 1 :   recovering state parameters       *'
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)

      DO PK=1, 6
      IF (LFLAG(PK)) THEN
       WRITE(IPPRC,*) "GREP: INFO :: STEP UP 'BRANCH_INFO'"
       WRITE(IPPRC,*) "GETVAL  STATE ",ORDER(PK)," NVAL 1 >>",
     >  REFNAM(PK),"_VAL<< ;"
      ENDIF
      ENDDO

      WRITE(IPPRC,*)
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)'!  STEP 2 :   interpolation of cross sections   *'
      WRITE(IPPRC,*)'!  warning => check the isotopes names          *'
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)

      WRITE(IPPRC,*)'EVALUATE ITER := ITER 1 + ;'
      IF (ITYP.EQ.0) WRITE(IPPRC,*)'  Micro := SCR: ',NAMSAP,' ::'
      IF (ITYP.EQ.1) WRITE(IPPRC,*)'  Micro := NCR: ',NAMSAP,' ::'
      WRITE(IPPRC,*)'   EDIT ',IPRINT
      IF (ITYP.EQ.0) THEN
       IF (HEQUI.NE.'NONE') WRITE(IPPRC,*)'   EQUI ',HEQUI
       IF (HMASL.NE.'NONE') WRITE(IPPRC,*)'   MASL ',HMASL
      ENDIF

      WRITE(IPPRC,*)'   MICRO LINEAR NMIX 1'
      IF (ITYP.EQ.0)WRITE(IPPRC,*)'   SAPHYB ',NAMSAP
      IF (ITYP.EQ.1)WRITE(IPPRC,*)'   COMPO ',NAMSAP,' ',
     >              TRIM(MIXDIR)

      WRITE(IPPRC,*)'   MIX 1'
      DO IOTH=1,NOTH
      WRITE(IPPRC,'(A,A)',advance='no')'     SET LINEAR ',
     > TRIM(OTHPK(IOTH))
      SELECT CASE (OTHTYP(IOTH))
       CASE (1)
        WRITE(IPPRC,*) ' ',INT(OTHVAR(IOTH))
       CASE (2)
        WRITE(IPPRC,*) ' ',OTHVAR(IOTH)
       CASE (3)
        WRITE(IPPRC,*) " '", TRIM(OTHVAC(IOTH)),"'"
      END SELECT
      ENDDO
      DO PK=1,6
       IF (LFLAG(PK)) THEN
      WRITE(IPPRC,*)'     SET LINEAR <<',REFNAM(PK),'>> <<',
     >               REFNAM(PK),'_VAL>>'
       ENDIF
      ENDDO

      IF (JOBOPT(2).EQ.'T') THEN
       CALL LCMSIX(IPDAT,'',0)
       CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
       CALL LCMSIX(IPDAT,'ISOTOPES',1)
       CALL LCMGTC(IPDAT,'XE135',12,1,ISOTOPES(1))
       CALL LCMGTC(IPDAT,'I135',12,1,ISOTOPES(2))
       CALL LCMGTC(IPDAT,'SM149',12,1,ISOTOPES(3))
       CALL LCMGTC(IPDAT,'PM149',12,1,ISOTOPES(4))
       CALL LCMGTC(IPDAT,'PM148',12,1,ISOTOPES(5))
       CALL LCMGTC(IPDAT,'PM148M',12,1,ISOTOPES(6))
       CALL LCMGTC(IPDAT,'ND147',12,1,ISOTOPES(7))
       CALL LCMGTC(IPDAT,'PM147',12,1,ISOTOPES(8))
       WRITE(IPPRC,*)'     MICRO ALL'

      DO I=1,8
      SELECT CASE (ISOTOPT)
      CASE ('*')
       WRITE(IPPRC,*)"     '",TRIM(ISOTOPES(I)),"'  *"
      CASE DEFAULT
       IF ((I.EQ.1).OR.(I.EQ.3).OR.(I.EQ.8)) THEN
         WRITE(IPPRC,*)"     '",TRIM(ISOTOPES(I)),"'  *"
       ELSE
        WRITE(IPPRC,*)"     '",TRIM(ISOTOPES(I)),"'",ISOTVAL
       ENDIF
      END SELECT
      ENDDO
      ENDIF
      WRITE(IPPRC,*)'    ENDMIX'

      IF ((JOBOPT(9).EQ.'T').AND.(ITYP.EQ.0) ) THEN
       WRITE(IPPRC,*)"    CHAIN"
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(2)),"  NG 0.0"
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(7)),"  NG 0.0"
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(1)),
     >  "  NG 0.0 FROM DECAY 1.0E+00 ",TRIM(ISOTOPES(2))
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(8)),
     >  "  NG 0.0 FROM DECAY 1.0E+00 ",TRIM(ISOTOPES(7))
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(5)),
     >  "  NG 0.0 FROM NG    5.3E-01 ",TRIM(ISOTOPES(8))
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(6)),
     >  " NG 0.0 FROM NG    4.7E-01 ",TRIM(ISOTOPES(8))
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(4)),
     > "  NG 0.0 FROM NG    1.0E+00 ",TRIM(ISOTOPES(5)),
     > " NG 1.0E+00 ",TRIM(ISOTOPES(6))
       WRITE(IPPRC,*)"     ",TRIM(ISOTOPES(3)),
     >  "  NG 0.0 FROM DECAY 1.0E+00 ",TRIM(ISOTOPES(4))
       WRITE(IPPRC,*)"      MACR     NFTOT 0.0"
       WRITE(IPPRC,*)"    ENDCHAIN"
      ENDIF
      WRITE(IPPRC,*)'    ;'
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)'!  STEP 3 :   branching calculation             *'
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)"IF ITER 1 = THEN "
      WRITE(IPPRC,*)"HELIOS GENPMAXS INFO Micro := D2P: "
      WRITE(IPPRC,*)"Micro INFO ",
     >  NAMSAP," ::"
      WRITE(IPPRC,*)"PHASE 2 EDIT",IPRINT,";"
      WRITE(IPPRC,*)"ELSE"
      WRITE(IPPRC,*)"HELIOS GENPMAXS INFO Micro := D2P: "
      WRITE(IPPRC,*)"Micro INFO GENPMAXS ",
     >  NAMSAP," HELIOS ::"
      WRITE(IPPRC,*)"PHASE 2 EDIT",IPRINT,";"
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)"ENDIF ;"
      WRITE(IPPRC,*)"Micro := DELETE: Micro ;"
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)"GREP: INFO :: STEP UP 'BRANCH_INFO'"
      WRITE(IPPRC,*)"GETVAL REWIND 1 NVAL 1 >>REWIND<< ;"

      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)'!  STEP 4 :   storing the current branch        *'
      WRITE(IPPRC,*)'!************************************************'
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)"IF REWIND 1 = THEN"
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)" HELIOS GENPMAXS INFO := D2P: INFO "
      WRITE(IPPRC,*)" GENPMAXS HELIOS ::"
      WRITE(IPPRC,*)" PHASE 3 EDIT",IPRINT," ;"

      WRITE(IPPRC,*)" GREP: INFO :: STEP UP 'BRANCH_INFO'"
      WRITE(IPPRC,*)" GETVAL STOP 1 NVAL 1 >>STOP<< ;"

      WRITE(IPPRC,*)"ENDIF ;"
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)"ENDWHILE ;"
      WRITE(IPPRC,*)
      WRITE(IPPRC,*)"END: ;"
      WRITE(IPPRC,*)"QUIT ."
      DEALLOCATE(PKEY)
      DEALLOCATE(OTHPK,OTHTYP,OTHVAC,OTHVAR)
      END
