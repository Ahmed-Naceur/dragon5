*DECK PMAXS
      SUBROUTINE D2P(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* PMAXS interface file generation.
*
*Copyright:
* Copyright (C) 2015 IRSN
*
*Author(s): 
* J. Taforeau
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
* None
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
      INTEGER :: STAVEC(40) = 0
      CHARACTER TEXT*72
      INTEGER UV
      INTEGER :: NFC1 = 0
      INTEGER :: NFC2 = 0
      INTEGER :: NFC3 = 0
      INTEGER :: NFC4 = 0
      INTEGER :: NXC  = 0
      INTEGER :: JOBT  = 0
      INTEGER :: NGP = 2
      INTEGER :: NCRD = 0
      INTEGER :: MIX = 1
      INTEGER :: FA_K = -1
      INTEGER :: IUPS = 0
      INTEGER :: USRSTA = 0
      INTEGER :: XESM = 3
      INTEGER :: ITEMP = 0
      INTEGER :: NOTHPK = 0
      INTEGER :: IOTHPK = 0
      REAL :: VERS = 3.0
      REAL :: SFAC = 1.0
      REAL :: BFAC = 1.0
      REAL :: THCK = -1.
      REAL FLOTT
      INTEGER PHASE, ITYPLIR, NITMA
      DOUBLE PRECISION DFLOT
      INTEGER :: IPRINT = -1
      INTEGER,DIMENSION(20) :: CRDINF = -1
      INTEGER,DIMENSION(12) :: USRVAL = 0
      INTEGER ,DIMENSION(12) :: OTHTYP = 2
      REAL,DIMENSION(5) :: LOCYLD = (/0.,-1.,-1.,-1.,-1. /)
      REAL,DIMENSION(5)::FC1=(/17.0,17.0,3.0,0.0,0.73659 /)
      REAL,DIMENSION(8)::FC2
      DATA FC2/6.2506E-01,1E-04,6*0.0/
      REAL,DIMENSION(7)::FC3
      DATA FC3/2.4921E+02, 2.4921E+02, 2.4921E+02, 2.3020E+01,
     1    1.4407E+02, 4.5099E+01, 4.5099E+01/
      REAL,DIMENSION(3)::FC4
      DATA FC4/1.44270E+00, 7.21350E-01, 7.21350E-01/
      REAL,DIMENSION(3)::XSC
      DATA XSC/ 1.0, 1.0, 5.32151E-01/
      REAL,DIMENSION(3)::YLD
      DATA YLD/ 0.06386, 0.00228, 0.0113/
      CHARACTER*16 :: JOBTIT = 'D2P.PMAXS'
      CHARACTER*40 :: COM = 'PWR CASE : UOX/MOX CORE FUEL'
      CHARACTER*12 :: FILNAM = 'HELIOS.dra'
      CHARACTER*12 :: MIXDIR = 'default     '
      CHARACTER*12 :: HDET = 'NULL        '
      CHARACTER*4 :: DER = 'T'
      CHARACTER*1 :: JOBOPT(16)
      CHARACTER*5 :: MESH = 'SAP'
      CHARACTER*12 :: USRPAR(12) = '            '
      CHARACTER*12 :: OTHPK(12) = '            '
      CHARACTER*8 :: HCUR(2)= 'NUL'
      CHARACTER*8 :: HFLX(2)= 'NUL'

      CHARACTER*12,DIMENSION(12) :: OTHVAL = '            '
      REAL :: OTHVAR(12)

      REAL USRVAPK(12,10)
      CHARACTER*4 :: CRDMOD = ' '
      CHARACTER*3 :: ADF = 'NUL'
      CHARACTER*3 :: CDF = 'NUL'
      CHARACTER*8,DIMENSION(4) :: ADFD = 'FD_B    '
      CHARACTER*8,DIMENSION(8) :: CDFD = 'FD_C    '
      CHARACTER*3 :: GFF = 'NUL'
      CHARACTER*12,DIMENSION(6) :: PKEY
      DATA PKEY/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/
      CHARACTER*12,DIMENSION(6) :: REFNAM
      DATA REFNAM/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/
      CHARACTER*12,DIMENSION(8) :: ISOT
      DATA ISOT/ "XE135PF","SM149PF","I135PF","PM149PF","PM148PF",
     >  "PM148MPF","ND147PF","PM147PF"/
      DATA JOBOPT/14*'F',"",""/
      CHARACTER*3 :: YLDOPT = 'REF'
      CHARACTER*4 :: OPT = 'NONE'
      CHARACTER*4 :: HEQUI = 'NONE'
      CHARACTER*4 :: HMASL = 'NONE'
      CHARACTER*1 :: ISOTOPT = '*'
      REAL :: ISOTVAL = 0.
      LOGICAL :: SAP=.FALSE.
      LOGICAL :: MIC=.TRUE.
      LOGICAL :: EXCESS=.FALSE.
      LOGICAL :: SCAT=.FALSE.
      LOGICAL :: LADD=.FALSE.
      LOGICAL :: LNEW=.FALSE.
      LOGICAL :: LPRC=.FALSE.
      LOGICAL :: LMEM=.FALSE.
      LOGICAL :: LCOR=.FALSE.
      OTHVAR(:) = -1
*----
*  parameters VALIDATION
*----
*----
*  RECOVER iPHASE AND iPRINT INDICES
*----
      WRITE(6,*) "****************************************************"
      WRITE(6,*) "*          RECOVERING D2P: DATA INPUT              *"
      WRITE(6,*) "****************************************************"

 100  CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
      IF (ITYPLIR.NE.3) THEN
       CALL XABORT ('@D2P: KEYWORD EXPECTED AS INPUT OF D2P: MODULE')
      ELSE IF (ITYPLIR.EQ.3) THEN
       IF (TEXT.EQ.'PHASE' ) THEN
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.EQ.1) THEN
         PHASE=NITMA
         GO TO 100
        ELSE
         CALL XABORT('@D2P: INTEGER EXPECTED AFTER PHASE KEYWORD')
        ENDIF
       ELSE IF (TEXT.EQ.'EDIT' ) THEN
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.EQ.1) THEN
         IPRINT=NITMA
         GO TO 100
        ELSE
         CALL XABORT('@D2P: INTEGER EXPECTED AFTER EDIT KEYWORD')
        ENDIF
       ELSE IF (TEXT.EQ.'MIX') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF (ITYPLIR.EQ.1) THEN
          MIX = NITMA
          GO TO 100
         ELSE
           CALL XABORT('@D2P: INTEGER EXPECTED AFTER MIX KEYWORD')
         ENDIF
       ELSE IF (TEXT.EQ.'NAMDIR' ) THEN
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.EQ.3) THEN
         MIXDIR=TEXT(1:12)
         IF(NITMA.GT.12) CALL XABORT('@D2P: C*12 EXPECTED FOR NAMDIR')
         GO TO 100
        ELSE
         CALL XABORT('@D2P: CHARACTER EXPECTED AFTER NAMDIR')
        ENDIF
       ELSE IF (TEXT.EQ. 'TEMP') THEN
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.EQ.3) THEN
         IF (NITMA.NE.1) THEN
          CALL XABORT('@D2P: "C" or "K" EXPECTED AFTER TEMP KEYWORD')
         ELSE
          IF (TEXT.EQ. 'C') THEN
           ITEMP=0
           GO TO 100
          ELSE IF (TEXT.EQ. 'K') THEN
           ITEMP=1
           GO TO 100
          ELSE
           CALL XABORT('@D2P: "C" or "K" EXPECTED AFTER TEMP KEYWORD')
          ENDIF
         ENDIF
        ELSE
          CALL XABORT('@D2P: CHARACTER EXPECTED AFTER TEMP KEYWORD')
        ENDIF
       ELSE IF (TEXT .EQ. 'PKEY') THEN
  15    CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.NE.3) THEN
         CALL XABORT('@D2P: CHARACTER EXPECTED AFTER PKEY KEYWORD')
        ELSE
         IF (TEXT.EQ.REFNAM(1)) THEN
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          PKEY(1)=TEXT(:12)
          GO TO 15
         ELSE IF (TEXT.EQ.REFNAM(2)) THEN
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          PKEY(2)=TEXT(:12)
          GO TO 15
         ELSE IF (TEXT.EQ.REFNAM(3)) THEN
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          PKEY(3)=TEXT(:12)
          GO TO 15
         ELSE IF (TEXT.EQ.REFNAM(4)) THEN
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          PKEY(4)=TEXT(:12)
          GO TO 15
         ELSE IF (TEXT.EQ.REFNAM(5)) THEN
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          PKEY(5)=TEXT(:12)
          GO TO 15
         ELSE IF (TEXT.EQ.REFNAM(6)) THEN
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          PKEY(6)=TEXT(:12)
          GO TO 15
         ELSE IF (TEXT.EQ.'ENDPKEY') THEN
         GO TO 100
         ELSE
          CALL XABORT('@D2P: UNKNOWN PKEY NAME : '//TEXT//'.')
         ENDIF
        ENDIF
      ELSE IF (TEXT .EQ. 'OTHER') THEN
       IOTHPK=0
       CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF(ITYPLIR.NE.1)THEN
         CALL XABORT('@D2P: INTEGER EXPECTED AFTER OTHPK CARD')
        ENDIF
       NOTHPK=NITMA
       STAVEC(20)=NOTHPK
       DO WHILE (IOTHPK.LT.NOTHPK)
       CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF(ITYPLIR.NE.3)THEN
         CALL XABORT('@D2P: C*12 (othnam) EXPECTED AFTER OTHPK CARD')
        ELSE
         IF(NITMA.GT.12)THEN
          CALL XABORT('@D2P: C*12 EXPECTED AFTER OTHPK CARD')
         ELSE
          IOTHPK=IOTHPK+1
          OTHPK(IOTHPK)=TEXT(:12)
         ENDIF
        ENDIF
       CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF(ITYPLIR.NE.3)THEN
         CALL XABORT('@D2P: C*1 ((othtyp) EXPECTED AFTER OTHPK CARD')
        ELSE
         IF(NITMA.GT.1)THEN
          CALL XABORT('@D2P: C*1 (othtyp) EXPECTED AFTER OTHPK CARD')
         ELSE
          IF (TEXT.EQ.'R') THEN
           OTHTYP(IOTHPK)=2
          ELSE IF (TEXT.EQ.'I') THEN
           OTHTYP(IOTHPK)=1
          ELSE IF (TEXT.EQ.'S') THEN
           OTHTYP(IOTHPK)=3
          ELSE
           WRITE(6,*) '@D2P: UNKNOWN TYPE (',TEXT(:1),') FOR (',
     >     OTHPK(IOTHPK),') PKEY.'
           CALL XABORT('@D2P: PLEASE USE I/R or S')
          ENDIF
         ENDIF
        ENDIF
       CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          IF(ITYPLIR.EQ.OTHTYP(IOTHPK))THEN
           IF (ITYPLIR.EQ.1) THEN
           WRITE(OTHVAL(IOTHPK),*)NITMA
           OTHVAR(IOTHPK)=NITMA
           ENDIF
           IF (ITYPLIR.EQ.2) THEN
            WRITE(OTHVAL(IOTHPK),'(f12.5)')FLOTT
            OTHVAR(IOTHPK)=FLOTT
           ENDIF
           IF (ITYPLIR.EQ.3) OTHVAL(IOTHPK)=TEXT(:12)
          ELSE
           CALL XABORT('@D2P: INCONSISTENT VALUE (othval)')
          ENDIF
       ENDDO
       GO TO 100
      ELSE IF (TEXT .EQ. 'ADF') THEN
  17    CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF(ITYPLIR.NE.3)THEN
         CALL XABORT('@D2P: C*3 EXPECTED AFTER ADF CARD')
        ENDIF
        IF(NITMA.GT.5) THEN
         CALL XABORT('@D2P: C*3 OR C*5 EXPECTED AFTER ADF CARD')
        ENDIF
        ADF=TEXT(:3)
        IF (TEXT(:5).EQ.'MERGE') THEN
         STAVEC(21)=1
         GO TO 17
        ELSE IF ((ADF.NE.'SEL') .AND. (ADF .NE.'GET')
     >       .AND. (ADF .NE.'DRA').AND. (ADF .NE.'GEN')) THEN
         WRITE(6,*) "@D2P: UNKNOWN KEYWORD :", ADF
         CALL XABORT('@D2P: DRA, SEL OR GET EXPECTED AFTER ADF CARD')
        ENDIF
        IF (ADF.EQ.'DRA') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF(ITYPLIR.NE.1) THEN
          CALL XABORT('@D2P: INTEGER EXPECTED AFTER ADF DRA CARD')
         ENDIF
         STAVEC(13)=NITMA !NADF
         IF((NITMA.NE.1).AND.(NITMA.NE.4)) THEN
          CALL XABORT('@D2P: 1 or 4 EXPECTED AFTER ADF DRA CARD')
         ENDIF
         DO I=1,STAVEC(13)
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          IF(ITYPLIR.NE.3) THEN
           CALL XABORT('@D2P: NADF STRING EXPECTED AFTER ADF DRA'
     >     //' CARD')
          ENDIF
          ADFD(I)=TEXT(:8)
         ENDDO
         GO TO 100
        ELSEIF (ADF .EQ. 'GEN') THEN
        STAVEC(13)= 1 !NADF
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF((ITYPLIR.NE.3).AND.(TEXT(:5).NE.'THICK')) THEN
          CALL XABORT('@D2P: REFLECTOR THICKNESS (THICK) EXPECTED')
        ELSE
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF (ITYPLIR.NE.2) THEN
          CALL XABORT('@D2P: REAL EXPECTED FOR REFLECTOR THICKNESS')
         ELSE
          THCK=FLOTT
         ENDIF
        ENDIF
        DO J=1,2
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF(ITYPLIR.NE.3) THEN
          CALL XABORT('@D2P: CHARACTER EXPECTED AFTER ADF GEN CARD')
        ENDIF
        IF ((TEXT(:4).NE.'FLUX').AND.(TEXT.NE.'CURR'))THEN
         CALL XABORT('@D2P: FLUX OR CURR KEYWORD EXPECTED AFTER GEN')
        ELSE
         DO I=1,2
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          IF(ITYPLIR.NE.3) THEN
           CALL XABORT('@D2P: CHARACTER EXPECTED CURR OR FLUX')
          ENDIF
          IF (J.EQ.1)HFLX(I)=TEXT(:8)
          IF (J.EQ.2)HCUR(I)=TEXT(:8)
         ENDDO
        ENDIF
        ENDDO
        GO TO 100
        ENDIF
       ELSE IF (TEXT .EQ. 'CDF') THEN
       CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF(ITYPLIR.NE.3)THEN
         CALL XABORT('@D2P: C*3 EXPECTED AFTER CDF CARD')
        ENDIF
        IF(NITMA.NE.3) THEN
         CALL XABORT('@D2P: C*3 EXPECTED AFTER CDF CARD')
        ENDIF
        CDF=TEXT(:3)
        IF ((CDF .NE.'DRA')) THEN
         WRITE(6,*) "@D2P: UNKNOWN KEYWORD :", CDF
         CALL XABORT('@D2P: DRA EXPECTED AFTER CDF CARD')
        ENDIF
        IF (CDF.EQ.'DRA') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF(ITYPLIR.NE.1) THEN
          CALL XABORT('@D2P: Integer EXPECTED AFTER CDF DRA CARD')
         ENDIF
         STAVEC(15)=NITMA !NCDF
         IF((NITMA.NE.1).AND.(NITMA.NE.2).AND.(NITMA.NE.3).AND.
     >      (NITMA.NE.4).AND.(NITMA.NE.5).AND.(NITMA.NE.8)) THEN
          CALL XABORT('@D2P: 1 to 5 or 8 EXPECTED AFTER CDF DRA'
     >    //' CARD')
         ENDIF
         DO I=1,STAVEC(15)
          CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
          IF(ITYPLIR.NE.3) THEN
           CALL XABORT('@D2P: NCDF String EXPECTED AFTER CDF DRA'
     >     //' CARD')
          ENDIF
          CDFD(I)=TEXT(:8)
         ENDDO
        ENDIF
        GO TO 100
       ELSE IF (TEXT .EQ. 'GFF') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF(ITYPLIR.NE.3)THEN
          CALL XABORT('@D2P: C*3 EXPECTED AFTER GFF CARD')
         ENDIF
         IF(NITMA.NE.3) THEN
          CALL XABORT('@D2P: C*3 EXPECTED AFTER GFF CARD')
         ENDIF
         GFF=TEXT(:3)
         IF ((GFF .NE.'DRA')) THEN
          WRITE(6,*) "@D2P: UNKNOWN KEYWORD :", GFF
          CALL XABORT('@D2P: DRA EXPECTED AFTER GFF CARD')
         ENDIF
         GO TO 100
       ELSE IF (TEXT.EQ.'FUEL' ) THEN
        FA_K=1
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.NE.3) THEN
         CALL XABORT ('@D2P: KEYWORD BARR EXPECTED AFTER FUEL CARD')
        ELSE IF (TEXT.EQ.'BARR') THEN
  10     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF (ITYPLIR.EQ.1 .AND. TEXT .NE.'ENDBARR') THEN
          IF (TEXT .EQ. 'DEF' .OR. TEXT .EQ. 'USER') THEN
           NCRD = NCRD + 1
           IF (NCRD <= 20) THEN
            CRDINF(NCRD)=NITMA
            GO TO 10
           ELSE
            CALL XABORT('@D2P: NUMBER OF BARR COMPOSITIONS EXCEED 20')
           ENDIF
          ELSE
            CALL XABORT('@D2P: DEF OR USER KEYWORD EXPECTED AFTER BARR'
     1      //' KEYWORD')
          ENDIF
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (TEXT .NE. ';' ) THEN
  11       IF (TEXT .EQ. 'GRID') THEN
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.NE.3) THEN
             CALL XABORT('@D2P: CHARACTER EXPECTED AFTER GRID KEYWORD')
            ELSE
             IF((TEXT.EQ.'SAP').OR.(TEXT.EQ.'DEF').OR.
     1         (TEXT.EQ.'USER')) THEN
              IF (TEXT.EQ.'SAP') THEN
              MESH=TEXT(:5)
               GO TO 10
              ELSE IF (TEXT.EQ.'DEF') THEN
              MESH=TEXT(:5)
               GO TO 10
              ELSE IF (TEXT.EQ.'USER') THEN
  12           CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
               IF (ITYPLIR.NE.3) THEN
                CALL XABORT('@D2P: KEYWORD EXPECTED AFTER USER KEYWORD')
               ELSE IF (TEXT.EQ.'NEW') THEN
                LNEW=.TRUE.
                GO TO 12
               ELSE IF (TEXT.EQ.'GLOBAL') THEN
                IF (LNEW) THEN
                 CALL XABORT('@D2P: INCOMPATIBLE OPT GLOBAL WITH NEW')
                ENDIF
                MESH='GLOB'
  90            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
                IF  (ITYPLIR.NE.3) THEN
                 CALL XABORT('@D2P: PKEY NAME EXPECTED IN GLOBAL OPT')
                ELSE
                 IF (TEXT.EQ.'ENDGLOBAL') GO TO 10
                 IF(NITMA > 12) THEN
                  CALL XABORT('@D2P: PKEY NAME IN GLOBAL MUST BE C*12')
                 ELSE
                  CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
                  USRSTA = USRSTA+1
                  USRPAR(USRSTA)=TEXT(:12)
                  IF  (ITYPLIR.NE.1) THEN
                   CALL XABORT ('@D2P: NB OF VALUES FOR STATE '//TEXT//
     1              ' EXPECTED')
                  ELSE
                   USRVAL(USRSTA)=NITMA
                   GO TO 90
                  ENDIF
                 ENDIF
                ENDIF
               ELSE IF (TEXT.EQ.'ADD') THEN
                MESH='ADD'
                LADD=.TRUE.
  95            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
                IF  (ITYPLIR.NE.3) THEN
                 CALL XABORT('@D2P: PKEY NAME EXPECTED IN USER ADD OPT')
                ELSE
                 IF (TEXT.EQ.'ENDADD') GO TO 10
                 IF(NITMA.GE.12) THEN
                  CALL XABORT('@D2P: STATE NAME IN ADD MUST BE C*12')
                 ELSE
                  USRSTA = USRSTA+1
                  USRPAR(USRSTA)=TEXT(:12)
                  CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
                  IF  (ITYPLIR.NE.1) THEN
                   CALL XABORT('@D2P: NB OF VALUES FOR STATE '//TEXT//
     1              'EXPECTED')
                  ELSE
                   USRVAL(USRSTA)=NITMA
                   DO UV=1,USRVAL(USRSTA)
                   CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
                    IF (ITYPLIR .NE. 2) THEN
                     CALL XABORT ('@D2P: REAL EXPECTED IN USER ADD OPT')
                    ELSE
                     USRVAPK(USRSTA,UV)=FLOTT
                    ENDIF
                   ENDDO
                   UV=1
                   GO TO 95
                  ENDIF
                 ENDIF
                ENDIF
               ELSE
                CALL XABORT('@D2P: UNKNOWN OPTION '//TEXT//
     1                'FOR USER OPT')
               ENDIF
              ENDIF
             ELSE
              CALL XABORT('@D2P: UNKNOWN OPTION FOR GRID KEYWORD')
             ENDIF
            ENDIF
           ELSE IF (TEXT.EQ.'DEF' .OR. TEXT.EQ.'USER') THEN
            CRDMOD=TEXT(:4)
            GO TO 10
           ELSE IF (TEXT.EQ.'ENDBARR') THEN
            IF (CRDMOD=='DEF') THEN
              CALL XABORT('@D2P: ENDBARR KEYWORD IS EXPECTED ONLY FOR'
     1        //' USER BARR COMPOSITION')
            ELSE
              GO TO 10
            ENDIF
           ELSE IF (TEXT .EQ. 'SCATTERING') THEN
             SCAT=.TRUE.
             GO TO 10
           ELSE IF (TEXT .EQ. 'DET') THEN
             CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
             IF(ITYPLIR.NE.3)THEN
              CALL XABORT('@D2P: C*12 EXPECTED AFTER DET CARD')
             ENDIF
             IF(NITMA.GT.12) THEN
              CALL XABORT('@D2P: C*12 EXPECTED AFTER GFF CARD')
             ENDIF
             HDET=TEXT(:12)
             GO TO 10
           ELSE IF (TEXT .EQ. 'ABSORPTION') THEN
   5        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.NE.3) THEN
             CALL XABORT('@D2P: C*6 EXPECTED AFTER ABSORPTION KEYWORD')
            ELSE IF (TEXT .EQ. 'MIC') THEN
             MIC =.TRUE.
             GO TO 10
            ELSE IF (TEXT .EQ. 'SAP') THEN
             SAP =.TRUE.
             MIC =.FALSE.
             GO TO 5
            ELSE IF (TEXT .EQ. 'EXCESS') THEN
             IF (SAP .EQV. .FALSE.) THEN
              CALL XABORT('@D2P: SAP KEYWORD EXPECTED BEFORE EXCESS')
             ELSE
              EXCESS = .TRUE.
              GO TO 10
             ENDIF
            ELSE
             GO TO 11
            ENDIF
           ELSE IF (TEXT .EQ. 'ISOTOPES') THEN
  25        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.NE.3) THEN
             CALL XABORT('@D2P: CHARACTER EXPECTED AFTER PKEY KEYWORD')
            ELSE
             IF (TEXT.EQ.'XE135') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(1)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'SM149') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(2)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'I135') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(3)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'PM149') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(4)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'PM148') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(5)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'PM148M') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(6)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'PM147') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(8)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'ND147') THEN
              CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
              ISOT(7)=TEXT(:12)
              GO TO 25
             ELSE IF (TEXT.EQ.'ENDISOTOPES') THEN
             GO TO 10
             ELSE
              CALL XABORT('@D2P: UNKNOWN NAME OF ISOTOPE: '//TEXT//'.')
             ENDIF
            ENDIF
           ELSE IF (TEXT .EQ. 'YLD') THEN
  37        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.NE.3) THEN
             CALL XABORT('@D2P: CHARACTER EXPECTED AFTER YLD KEYWORD')
            ELSE
             IF (TEXT(:3).EQ.'COR') THEN
              LCOR=.TRUE.
              GO TO 37
             ELSE
              IF (.NOT.LCOR) STAVEC(22)=0
              YLDOPT=TEXT(:3)
              IF (YLDOPT.EQ.'REF') THEN
               IF (LCOR) STAVEC(22)=1
               GO TO 10
              ELSE IF (YLDOPT.EQ.'FIX') THEN
               IF (LCOR) THEN
                WRITE (6,*) '@D2P : NO CORRECTION POSSIBLE OF FISSION'
                CALL XABORT ('YIELDS WITH THE FIX OPTION')
               ENDIF
               DO I=1, 3
                CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
                IF (ITYPLIR==2) THEN
                 YLD(I)=FLOTT
                ELSE
                 CALL XABORT('REAL EXPECTED FOR YIELD VALUES')
                ENDIF
               ENDDO
               GO TO 10
              ELSE IF (YLDOPT.EQ.'MAN') THEN
               IF (LCOR) STAVEC(22)=2
  35           CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
               I=1
               DO WHILE (TEXT.NE.REFNAM(I).AND.(I.LE.5))
                I=I+1
               ENDDO
               IF (I.GT.5) THEN
                IF (TEXT.EQ.'ENDMAN') GO TO 10
                CALL XABORT('@D2P: PKEY NAME ('//TEXT(:12)//') NOT '
     >          //'ALLOWED IN YIELD CARD')
               ELSE
                CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

                IF (ITYPLIR.EQ.2) THEN
                 LOCYLD(I)=FLOTT
                 GO TO 35
                ELSE
                CALL XABORT('@D2P: SOMETHING WRONG OCCURS IN YLD CARD')
                ENDIF
               ENDIF
              ELSE
               CALL XABORT('@D2P: UNKNOWN OPTION FOR YLD: '//TEXT//'.')
              ENDIF
             ENDIF
            ENDIF
           ELSE IF (TEXT.EQ.'GENPMAXS') THEN
            GO TO 120
           ELSE IF (TEXT.EQ.'HELIOS') THEN
            GO TO 21
           ELSE IF (TEXT.EQ.'PROC') THEN
            GO TO 220
           ELSE
            CALL XABORT('@D2P: SOMETHING WRONG OCCURS IN INPUT DATA')
           ENDIF
          ELSE
           GOTO 200
          ENDIF
         ENDIF
        ELSE
          CALL XABORT('@D2P: UNKNOWN KEYWORD '//TEXT//', BARR EXPECTED')
        ENDIF
       ELSE IF (TEXT .EQ. 'REFLECTOR') THEN
        FA_K=0
        NCRD=1
        CRDINF(1)=1
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
        IF (ITYPLIR.NE.3) THEN
         CALL XABORT('@D2P: SOMETHING WRONG OCCURS IN REFLECTOR DATA')
        ELSE IF (ITYPLIR.EQ.3) THEN
         IF ((TEXT.NE.';'))THEN
          IF (TEXT.EQ.'GENPMAXS') THEN
           GO TO 120
          ELSE IF (TEXT.EQ.'HELIOS') THEN
           GO TO 21
          ELSE IF (TEXT.EQ.'PROC') THEN
           GO TO 220
          ELSE
           CALL XABORT('@D2P: SOMETHING WRONG OCCURS IN INPUT DATA')
          ENDIF
         ELSE
           GOTO 200
         ENDIF
        ELSE
         CALL XABORT('@D2P: UNKNOWN KEYWORD'//TEXT//', IN INPUT DATA')
        ENDIF
       ENDIF

  21   IF (TEXT .EQ. 'HELIOS') THEN
        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
  30    IF (ITYPLIR.NE.3) THEN
         CALL XABORT ('@D2P: KEYWORD EXPECTED AFTER HELIOS CARD')
        ELSE IF (TEXT.EQ.'FILE_CONT_1') THEN
  40     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.EQ.1).OR.(ITYPLIR.EQ.2)) THEN
          NFC1 = NFC1 +1
          IF (NFC1 <= 5) THEN
           IF (ITYPLIR.EQ.1) FC1(NFC1) = NITMA
           IF (ITYPLIR.EQ.2) FC1(NFC1) = FLOTT
           GO TO 40
          ELSE
           CALL XABORT('@D2P: FIVE VALUES FOR FILE_CONT_1 ARE EXPECTED')
          ENDIF
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (NFC1.NE.5) THEN
           CALL XABORT('@D2P: FIVE VALUES FOR FILE_CONT_1 ARE EXPECTED')
          ENDIF
          IF (TEXT .NE. ';' ) THEN
           GO TO 30
          ELSE
           GOTO 200
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'FILE_CONT_2') THEN

  50     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

         IF (ITYPLIR.EQ.2) THEN
          NFC2 = NFC2 +1
          IF (NFC2 <= 8) THEN
           FC2(NFC2) = FLOTT
           GO TO 50
          ELSE
           CALL XABORT('@D2P: 8 VALUES AT MOST IN FILE_CONT_2')
          ENDIF
         ELSE IF (ITYPLIR.EQ.1) THEN
           CALL XABORT('@D2P: REAL VALUES EXPECTED IN FILE_CONT_2')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (NFC2<2) THEN
           CALL XABORT('@D2P: 2 VALUES AT LEAST IN FILE_CONT_2')
          ENDIF
          IF (TEXT .NE. ';' ) THEN
           GO TO 30
          ELSE
           GOTO 200
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'FILE_CONT_3') THEN

  60     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

         IF (ITYPLIR.EQ.2) THEN
          NFC3 = NFC3 +1
          IF (NFC3 <= 7) THEN
           FC3(NFC3) = FLOTT
           GO TO 60
          ELSE
           CALL XABORT('@D2P: 7 VALUES IN FILE_CONT_3 EXPECTED')
          ENDIF
         ELSE IF (ITYPLIR.EQ.1) THEN
           CALL XABORT('@D2P: REAL VALUES EXPECTED IN FILE_CONT_3  ')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (NFC3<7) THEN
           CALL XABORT('@D2P: 7 VALUES FOR FILE_CONT_3 EXPECTED')
          ENDIF
          IF (TEXT .NE. ';' ) THEN
           GO TO 30
          ELSE
           GOTO 200
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'FILE_CONT_4') THEN

  70     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

         IF (ITYPLIR.EQ.2) THEN
          NFC4 = NFC4 +1
          IF (NFC4 <= 3) THEN
           FC4(NFC4) = FLOTT
           GO TO 70
          ELSE
           CALL XABORT('@D2P: 3 VALUES IN FILE_CONT_4 EXPECTED')
          ENDIF
         ELSE IF (ITYPLIR.EQ.1) THEN
           CALL XABORT('@D2P: REAL VALUES EXPECTED IN FILE_CONT_4')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (NFC4<3) THEN
           CALL XABORT('@D2P: 3 VALUES IN FILE_CONT_4 EXPECTED')
          ENDIF
          IF (TEXT .NE. ';' ) THEN
           GO TO 30
          ELSE
           GOTO 200
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'XS_CONT') THEN

  80     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

         IF ((ITYPLIR.EQ.1).OR.(ITYPLIR.EQ.2)) THEN
          NXC = NXC +1
          IF (NXC <= 3) THEN
           IF (ITYPLIR.EQ.1) XSC(NXC) = NITMA
           IF (ITYPLIR.EQ.2) XSC(NXC) = FLOTT
           GO TO 80
          ELSE
           CALL XABORT('@D2P: 3 VALUES IN XS_CONT ARE EXPECTED')
          ENDIF
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (NXC<3) THEN
           CALL XABORT('@D2P: 3 VALUES FOR XS_CONT ARE EXPECTED')
          ENDIF
          IF (TEXT .NE. ';' ) THEN
           GO TO 30
          ELSE
           GOTO 200
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'GENPMAXS')THEN
         GO TO 120
        ELSE IF (TEXT.EQ.'PROC')THEN
         GO TO 220
        ELSE
         CALL XABORT ('@D2P: UNKNOWN KEYWORD: '//TEXT//'.')
        ENDIF
       ENDIF



  120  IF (TEXT .EQ. 'GENPMAXS') THEN

        CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
  130   IF (ITYPLIR.NE.3) THEN
         CALL XABORT ('@D2P: KEYWORD EXPECTED AFTER GENPMAXS CARD')
        ELSE IF (TEXT.EQ.'JOB_TIT') THEN

         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

         IF ((ITYPLIR.NE.3)) THEN
          CALL XABORT('@D2P: CHARACTER EXPECTED AFTER JOB_TIT CARD')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (TEXT .NE. ';' ) THEN
           JOBTIT=TEXT(:16)
           IF (NITMA>16) THEN
            CALL XABORT('@D2P: JOB_TIT NAME TOO LONG (>C*16)')
           ENDIF
           CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
           GO TO 130
          ELSE
           CALL XABORT('@D2P: CHARACTER EXPECTED AFTER JOB_TIT CARD')
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'FILE_NAME') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.3)) THEN
          CALL XABORT('CHARACTER EXPECTED AFTER JOB_TIT CARD')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (TEXT .NE. ';' ) THEN
           FILNAM=TEXT(:12)
           IF (NITMA>12) THEN
            CALL XABORT('FILE_NAME NAME TOO LONG (>C*12)')
           ENDIF
           CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
           GO TO 130
          ELSE
           CALL XABORT('@D2P: CHARACTER EXPECTED AFTER FILE_NAME CARD')
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'DERIVATIVE') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.3)) THEN
          CALL XABORT('@D2P: CHARACTER EXPECTED AFTER DERIVATIVE CARD')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (TEXT .NE. ';' ) THEN
           IF ((TEXT.EQ.'T').OR.(TEXT.EQ.'F')) THEN
            DER=TEXT(:4)
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            GO TO 130
           ELSE
            CALL XABORT('@D2P: (T/F) EXECTED AFTER DERIVATIVE CARD')
           ENDIF
          ELSE
           CALL XABORT('@D2P: (T/F) EXPECTED AFTER DERIVATIVE CARD')
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'COMMENT') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.3)) THEN
          CALL XABORT('@D2P: CHARACTER EXPECTED AFTER COMMENT CARD')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (TEXT .NE. ';' ) THEN
           COM=TEXT(:40)
           IF (NITMA>40)THEN
            CALL XABORT('@D2P: COMMENT NAME TOO LONG (>C*40)')
           ENDIF
           CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
           GO TO 130
          ELSE
           CALL XABORT('@D2P: CHARACTER EXPECTED AFTER COMMENT CARD')
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'JOB_OPT') THEN
 140     CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

         IF ((ITYPLIR.NE.3)) THEN
          CALL XABORT('@D2P: CHARACTER EXPECTED AFTER JOB_OPT CARD')
         ELSE IF (ITYPLIR.EQ.3) THEN
          IF (TEXT .NE. ';' ) THEN
           JOBT=JOBT+1

           IF ((TEXT.EQ.'T').OR.(TEXT.EQ.'F')) THEN
            IF (JOBT<=14) THEN
             JOBOPT(JOBT)=TEXT(:1)
             GO TO 140
            ELSE
             WRITE (6,*) '@D2P: LAST JOB_OPT VALUE :', TEXT
             CALL XABORT('@D2P: 14 VALUES EXPECTED FOR JOB_OPT CARD')
            ENDIF
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            GO TO 130
           ELSE
             IF (JOBT<=14) THEN
             WRITE (6,*)'@D2P: ',JOBT,'th JOB_OPT VALUE :', TEXT
             CALL XABORT('@D2P: (T/F) VALUES EXPECTED FOR JOB_OPT CARD')
             ELSE
             GO TO 130
             ENDIF
           ENDIF
          ELSE IF (JOBT==14 .and. TEXT==';') THEN
           GO TO 190
          ELSE IF (JOBT==15) THEN
           GO TO 130
          ELSE
           CALL XABORT('@D2P: CHARACTER EXPECTED AFTER JOB_OPT CARD')
          ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'IUPS') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.1)) THEN
          CALL XABORT('@D2P: INTEGER EXPECTED AFTER IUPS CARD')
         ELSE IF (ITYPLIR.EQ.1) THEN
           IF ((NITMA>2).OR.(NITMA<0)) THEN
            CALL XABORT ('@D2P: IUPS INTEGER MUST BE 0,1 or 2')
           ELSE
            IUPS=NITMA
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.EQ.3) THEN
             GO TO 130
            ELSE
             CALL XABORT('@D2P: ONLY 1 VALUE IS EXPECTED FOR IUPS')
            ENDIF
           ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'XESM') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.1)) THEN
          CALL XABORT('@D2P: INTEGER EXPECTED AFTER XESM CARD')
         ELSE IF (ITYPLIR.EQ.1) THEN
           IF ((NITMA>3).OR.(NITMA<1)) THEN
            CALL XABORT ('@D2P: XESM CARD INTEGER MUST BE 1,2 or 3')
           ELSE
            XESM=NITMA
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.EQ.3) THEN
             GO TO 130
            ELSE
             CALL XABORT('@D2P: ONLY 1 VALUE IS EXPECTED FOR XESM')
            ENDIF
           ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'VERSION') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.2)) THEN
          CALL XABORT('@D2P: REAL EXPECTED AFTER VERSION CARD')
         ELSE IF (ITYPLIR.EQ.2) THEN
           IF ((FLOTT<0)) THEN
            CALL XABORT ('@D2P: VERSION NUMBER MUST BE POSITIVE')
           ELSE
            VERS=FLOTT
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.EQ.3) THEN
             GO TO 130
            ELSE
             CALL XABORT('@D2P: ONLY ONE VALUE IS EXPECTED FOR VERSION')
            ENDIF
           ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'SFAC') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.2)) THEN
          CALL XABORT('@D2P: REAL EXPECTED AFTER SFAC CARD')
         ELSE IF (ITYPLIR.EQ.2) THEN
           IF ((FLOTT<0)) THEN
            CALL XABORT ('@D2P: SFAC FACTOR MUST BE POSITIVE')
           ELSE
            SFAC=FLOTT
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.EQ.3) THEN
             GO TO 130
            ELSE
             CALL XABORT('@D2P: ONLY ONE VALUE IS EXPECTED FOR SFAC ')
            ENDIF
           ENDIF
         ENDIF
        ELSE IF (TEXT.EQ.'BFAC') THEN
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         IF ((ITYPLIR.NE.2)) THEN
          CALL XABORT('@D2P: REAL EXPECTED AFTER BFAC CARD')
         ELSE IF (ITYPLIR.EQ.2) THEN
           IF ((FLOTT<0)) THEN
            CALL XABORT ('@D2P: BFAC FACTOR MUST BE POSITIVE')
           ELSE
            BFAC=FLOTT
            CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
            IF (ITYPLIR.EQ.3) THEN
             GO TO 130
            ELSE
             CALL XABORT('@D2P: ONLY 1 VALUE IS EXPETCTED FOR BFAC')
            ENDIF
           ENDIF
         ENDIF
        ELSE IF (TEXT .EQ. ';' ) THEN
         GO TO 190
        ELSE IF (TEXT .EQ. 'PROC') THEN
         GO TO 220
        ELSE
         CALL XABORT ('@D2P: UNKNOWN KEYWORD: '//TEXT//'.')
        ENDIF
       ENDIF

  220  IF (TEXT .EQ. 'PROC') THEN
        LPRC=.TRUE.
  221   CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)

        IF (ITYPLIR.NE.3) THEN
         CALL XABORT ('@D2P: C*4 EXPECTED AFTER PROC CARD')
        ELSE IF (TEXT.EQ.';') THEN
         GO TO 190
        ELSE IF (TEXT.EQ.'MEMO')THEN
          LMEM=.TRUE.
          GO TO 221
        ELSE IF ((NITMA .NE. 4)) THEN
         CALL XABORT ('@D2P: C*4 EXPECTED AFTER PROC CARD')
        ELSE
         OPT=TEXT(:4)
         CALL REDGET(ITYPLIR,NITMA,FLOTT,TEXT,DFLOT)
         SELECT CASE (OPT)
         CASE ('ISOT')
          IF (TEXT.EQ."*") THEN
           ISOTOPT=TEXT(:1)
          ELSE IF(ITYPLIR.NE.2) THEN
           CALL XABORT('@D2P: * OR REAL EXPECTED AFTER ISOT CARD')
          ELSE
           ISOTOPT='R'
           ISOTVAL=FLOTT
          ENDIF
         CASE ('EQUI')
          IF(ITYPLIR.NE.3) THEN
           CALL XABORT('@D2P: HEQUI (C*4) EXPECTED AFTER EQUI CARD')
          ELSE
           IF (NITMA.NE.4) THEN
            CALL XABORT('@D2P: HEQUI (C*4) EXPECTED AFTER EQUI CARD')
           ELSE
            HEQUI=TEXT(:4)
           ENDIF
          ENDIF
         CASE ('MASL')
          IF(ITYPLIR.NE.3) THEN
           CALL XABORT('@D2P: HMASL (C*4) EXPECTED AFTER MASL CARD')
          ELSE
           IF (NITMA.NE.4) THEN
            CALL XABORT('@D2P: HMASL (C*4) EXPECTED AFTER MASL CARD')
           ELSE
            HMASL=TEXT(:4)
           ENDIF
          ENDIF

         CASE DEFAULT
          CALL XABORT('@D2P: UNKNOWN OPTION ('//OPT//') IN PROC CARD')
         END SELECT
         GO TO 221
        ENDIF


       ENDIF

 190   IF (TEXT(:1) .EQ. ';' ) THEN
        IF (NFC2.NE.0) NGP = NFC2
        IF (JOBT.NE.14) THEN
         IF (JOBT.NE.15) THEN
          IF (JOBT.NE.0)THEN
           CALL XABORT('@D2P: 14 VALUES EXPECTED FOR JOB_OPT')
          ENDIF
         ENDIF
        ENDIF
        GO TO 200
       ELSE
        CALL XABORT('@D2P: UNKNOWN KEYWORD:'//TEXT//'.')
       ENDIF
      ENDIF



 200  IF (PHASE.EQ.1) THEN
      IF ((ADF.EQ.'NUL') .and. (JOBOPT(1).EQ.'T')) THEN
       WRITE(6,*)"@D2P: ADF CALCULATION REQUIRED, PLEASE USE THE 'ADF'",
     >    " CARD."
       CALL XABORT("")
      ELSE IF ((ADF.NE.'NUL') .and. (JOBOPT(1).EQ.'F')) THEN
       WRITE(6,*)"@D2P: ADF CALCULATION REQUIRED, PLEASE TURN ON THE ",
     >    "'ladf' FLAG IN JOB_OPT CARD."
       CALL XABORT("")
      ENDIF
      IF ((CDF.EQ.'NUL') .and. (JOBOPT(10).EQ.'T')) THEN
       WRITE(6,*)"@D2P: CDF CALCULATION REQUIRED, PLEASE USE THE 'CDF'",
     >    " CARD."
       CALL XABORT("")
      ELSE IF ((CDF.NE.'NUL') .and. (JOBOPT(10).EQ.'F')) THEN
       WRITE(6,*)"@D2P: CDF CALCULATION REQUIRED, PLEASE TURN ON THE ",
     >    "'lcdf' FLAG IN JOB_OPT CARD."
       CALL XABORT("")
      ENDIF
            IF ((GFF.EQ.'NUL') .and. (JOBOPT(11).EQ.'T')) THEN
       WRITE(6,*)"@D2P: GFF CALCULATION REQUIRED, PLEASE USE THE 'GFF'",
     >    " CARD."
       CALL XABORT("")
      ELSE IF ((CDF.NE.'NUL') .and. (JOBOPT(10).EQ.'F')) THEN
       WRITE(6,*)"@D2P: GFF CALCULATION REQUIRED, PLEASE TURN ON THE ",
     >    "'lgff' FLAG IN JOB_OPT CARD."
       CALL XABORT("")
      ENDIF
      IF (FA_K==0) THEN
       IF ((ADF.EQ.'SEL').OR.(ADF.EQ.'GET')) THEN
        CALL XABORT('@D2P: ADF OF TYPE DRA EXPECTED FOR REFLECTOR CASE')
       ENDIF
       DO I=2, 16
        IF (JOBOPT(I).EQ.'T') THEN
         JOBOPT(I)='F'
         WRITE(6,*)"@D2P: JOB_OPT(",I,") SET TO 'F' FOR RELFECTOR CASE"
        ENDIF
       ENDDO
      ENDIF
      ENDIF
      CALL  D2PDRV(      NENTRY, HENTRY, IENTRY, JENTRY, KENTRY,    NGP,
     >                     NCRD,    MIX,   FA_K,   IUPS, USRSTA,  PHASE,
     >                   IPRINT, STAVEC, CRDINF, USRVAL,   VERS,   SFAC,
     >                     BFAC,   FC1,    FC2,    FC3,     FC4,    XSC,
     >                  USRVAPK,   ADF,    DER,  JOBOPT, USRPAR,   MESH,
     >                    PKEY, FILNAM,   ISOT,  JOBTIT,    COM,    SAP,
     >                     MIC,    EXC,   SCAT,    LADD,   LNEW, MIXDIR,
     >                     CDF,    GFF,   ADFD,    CDFD,    YLD, YLDOPT,
     >                  LOCYLD,   XESM,  ITEMP,   OTHPK, OTHTYP, OTHVAL,
     >                    HDET,   LPRC,  HEQUI,  HMASL ,ISOTOPT,ISOTVAL,
     >                    LMEM, OTHVAR,   THCK,    HFLX,   HCUR        )

      END
