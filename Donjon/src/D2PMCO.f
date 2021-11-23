*DECK D2PMCO
      SUBROUTINE D2PMCO(  IPSAP,  IPDAT, STAVEC, CRDINF,   NCRD,  PKNAM,
     >                    ISOT ,   MESH, USRPAR, USRVAL, USRSTA,USRVAPK,
     >                    SAP  ,    MIC,    EXC,   SCAT,    ADF,  LADD ,
     >                    LNEW ,   LADF, IPRINT, MIXDIR,    MIX,   LCDF,
     >                    LGFF ,    CDF,    GFF,   ADFD,   CDFD,   LYLD,
     >                      YLD, YLDOPT, LOCYLD,  OTHPK, OTHTYP, OTHVAL,
     >                   OTHREA,   THCK,   HFLX,   HCUR                )
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the global stated variable data contained in the SAPHYB object
*
*Author(s): 
* J. Taforeau
*
*Parameters: input/output
* IPDAT   address of the INFO data block
* IPSAP   address of the saphyb object
* NCRD    number of control rod composition recovered from D2P input
*         user
* MIX     index of mixture on which XS are to be extracted (only for
*         reflector cases)
* USRSTA  state variable names recovered from GLOBAL record in D2P:
* USRVAL  number of value for state variables  recovered from GLOBAL
*         record in D2P:
* IPRINT  control the printing on screen
* STAVEC  various parameters associated with the IPDAT structure
* CRDINF  meaning of control rods in the IPSAP object
* USRVAPK value of state prameter set by the user and recoverd from
*         USER ADD option in D2P:
* ADF     type of ADF to be selected
* DER     partials derivative (T) or row cross section (F) to be stored
*         in PMAXS
* USRPAR  name of state variables (sapnam) in IPSAP associated to
*         DMOD TCOM etc. recovered from PKEY card in D2P:
* MESH    type of meshing to be applied for the branching calculation
* PKNAM   name of state variable (refnam) recovered from PKEY card in
*         D2P:
* ISOT    name of isotopes in IPSAP for xenon samarium and spomethium
* SAP     flag to indicate that absorption cross section must be
*         directly recovered from IPSAP
* MIC     flag to indicate that absorption cross section must be
*         directly recovered from IPMIC
* EXC     flag to indicate that excess cross section is to be extracted
*         from absoption xs (only if SAP)
* SCAT    flag to indicate that scattering cross section must be
*         directly reconstructed from IPSAP
* LADD    flag to indicate that new points must be added to the IPSAP
*         original meshing
* LNEW    flag to indicate that only new points must be used during the
*         branching calculation
* LADF    Assembly Discontinuity Factors must be recovered
* MIXDIR  directory that contains homogeneous mixture information
* MIX     Index of mixture that contains homogeneous cross sections
* LCDF    Corner Discontinuity Factors must be recovered
* LGFF    Group Form Factors must be recovered
* CDF     type of CDF to be selected
* GFF     type of GFF to be selected
* ADFD    name of record for 'DRA' type of ADF
* CDFD    name of record for 'DRA' type of CDF
* LYLD    Fission Yield must be recovered
* YLD     user defined values for fission yields (1:I, 2:XE, 3:PM)
* LOCYLD  value for state parameter on which fission yield will be
*         calculated
* YLDOPT  option for fission yields calculation (DEF, MAN, FIX)
* OTHREA  real (or integer) value for OTHER parameter
* LMER    ADF are merged in the cross sections
* THCK    Thickness of reflector
* HFLX    Name of the record for the flux
* HCUR    Name of the record for the current
*
*Parameters: 
* OTHPK   
* OTHTYP  
* OTHVAL  
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPDAT
      INTEGER NCRD,USRSTA,MIX
      INTEGER IPRINT
      REAL THCK
      INTEGER STAVEC(40),CRDINF(20),USRVAL(12),OTHTYP(12)
      REAL USRVAPK(12,10),YLD(3),LOCYLD(5),OTHREA(12)
      CHARACTER*3 ADF,CDF,GFF,YLDOPT
      CHARACTER*8 ADFD(4),CDFD(8),HFLX(2),HCUR(2)
      CHARACTER*12 USRPAR(12),OTHVALC
      CHARACTER*5 MESH
      CHARACTER*12 PKNAM(6),OTHPK(12), OTHVAL(12)
      CHARACTER*12 ISOT(8), MIXDIR
      LOGICAL SAP, MIC, EXC,SCAT,LADD,LNEW,LADF,LCDF,LGFF,LYLD
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR)  IPROOT,IPDIR,IPTH,KPTH,IPMCO,JPMCO

      PARAMETER(NSTATE=40)
      INTEGER :: N_XS = 8
      INTEGER :: NTOT = 0
      INTEGER,DIMENSION(6)  :: ORDER_VAL = 0
      INTEGER DIMMCO(NSTATE),DIMCAL(NSTATE),DIMGEO(NSTATE)
      INTEGER NPAR,NCALS,NSVAR,NOTH
      INTEGER NCRD_SAP,NVALTMP(10)
      INTEGER RKOTH(STAVEC(20))
      INTEGER :: NOTHTH = 0
      REAL OTHR(20,20)
      REAL :: OTHVAR(20) = -1
      INTEGER i, j, k, l , n, UV
      REAL FIRST_VAL,LAST_VAL,PITCH
      LOGICAL LABS(3)
      LOGICAL :: LBARR = .FALSE.
      LOGICAL :: LDMOD = .FALSE.
      LOGICAL :: LCBOR = .FALSE.
      LOGICAL :: LTCOM = .FALSE.
      LOGICAL :: LTMOD = .FALSE.
      LOGICAL :: LBURN = .FALSE.
      LOGICAL :: LOTH(12) =.FALSE.
      CHARACTER(LEN=12)   PKEY_BARR(6), OTHC(20,20)
      CHARACTER(LEN=12)   :: OTHVAC(20) = 'NULL        '
      CHARACTER*12,DIMENSION(6) :: PKREF
      DATA PKREF/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NVAL, RANK,RANK_INDEX,PKIDX
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PKEY,PKEY_TMP
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PVALDIR
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT
!      REAL, ALLOCATABLE, DIMENSION(:) :: SV_VAL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VALPAR

      IPROOT=IPSAP
      NOTH=STAVEC(20)
      LABS(1)=MIC
      LABS(2)=SAP
      LABS(3)=EXC

      CALL LCMPUT(IPDAT,'BARR_INFO',NCRD,1,CRDINF)
      ! RECOVER DIMMCO INFORMATION FROM SAPHYB
      CALL XDISET(DIMMCO,NSTATE,0)
      CALL LCMSIX(IPSAP,MIXDIR,1)
      IPDIR=IPSAP
      CALL LCMGET(IPDIR,'STATE-VECTOR',DIMMCO)
      NGFF  = DIMMCO(14)
      IPMCO=LCMGID(IPDIR,'MIXTURES')
      JPMCO=LCMGIL(IPMCO,MIX)
      IPMCO=LCMGID(JPMCO,'CALCULATIONS')
      JPMCO=LCMGIL(IPMCO,1)
      CALL LCMGET (JPMCO,'STATE-VECTOR',DIMCAL)
      NPAR  = DIMMCO(5)
      NMIL  = DIMMCO(1)
      NCALS = DIMMCO(4)
      NDEL  = DIMCAL(19)
      ! RECOVER NPIN FOR GFF
      NPIN=1
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMPUT(IPDAT,'NPAR',1,1,NPAR)
      CALL LCMSIX(IPDAT,' ',0)
      IF(LGFF) THEN
        CALL LCMSIX(JPMCO,'MACROLIB',1)
        CALL LCMSIX(JPMCO,'GFF',1)
        CALL LCMSIX(JPMCO,'GFF-GEOM',1)
        CALL LCMGET(JPMCO,'STATE-VECTOR',DIMGEO)
        NXG=DIMGEO(3)
        NYG=DIMGEO(4)
        IF(NXG.NE. NYG) THEN
          WRITE(6,*) "@D2PMAC:",
     1     " NPIN NOT THE SAME X AND Y AXES IN MCO"
          CALL XABORT("=> NXG .NE. NYG")
        ENDIF
        NPIN=NXG
        CALL LCMSIX(JPMCO,' ',2)
        CALL LCMSIX(JPMCO,' ',2)
        CALL LCMSIX(JPMCO,' ',2)
      ENDIF

      ! INITIALIZATION OF PARAMETERS
      NSVAR = 0
      k = 1
      ! MEMORY ALLOCATION
      ALLOCATE (PKEY(NPAR))
      ALLOCATE (NVAL(NPAR))
      ALLOCATE (PKEY_TMP(NPAR))
      ALLOCATE (RANK(NPAR))
      ALLOCATE (RANK_INDEX(NPAR+1))
      ALLOCATE (PARFMT(NPAR))

      CALL LCMSIX(IPDIR,'GLOBAL',1)
      CALL LCMGTC(IPDIR,'PARKEY',12,NPAR,PKEY)
      CALL LCMGTC(IPDIR,'PARFMT',8,NPAR,PARFMT)

      CALL LCMGET(IPDIR,'NVALUE',NVAL)

      NVALTMP=NVAL

      ! LOOP OVER STATE VARIABLES OF SAPHYB
      ! CHECK OF EXISTENCE OF STATE PARAMETER

      DO i=1, NPAR
        IF ((PKEY(i).NE.'FLUE').AND.(PKEY(i).NE.'TIME')) NTOT=NTOT+1
        IF(PKEY(i)==PKNAM(1)) THEN      ! BARR
          LBARR=.TRUE.
        ELSE IF(PKEY(i)==PKNAM(2)) THEN ! DMOD
          LDMOD=.TRUE.
        ELSE IF(PKEY(i)==PKNAM(4)) THEN ! TCOM
          LTCOM=.TRUE.
        ELSE IF(PKEY(i)==PKNAM(5)) THEN ! TMOD
          LTMOD=.TRUE.
        ELSE IF(PKEY(i)==PKNAM(3)) THEN ! CBOR
          LCBOR=.TRUE.
        ELSE IF(PKEY(i)==PKNAM(6)) THEN ! BURN
          LBURN =.TRUE.
        ELSE
         DO j=1,NOTH
           IF (PKEY(i)==OTHPK(j)) THEN
            LOTH(j) = .TRUE.
            SELECT CASE (PARFMT(i))
            CASE ('REAL')
             IF (OTHTYP(j) .EQ. 2) GO TO 100
            CASE ('STRING')
             IF (OTHTYP(j) .EQ. 3) GO TO 100
            CASE ('INTEGER')
             IF (OTHTYP(j) .EQ. 1) GO TO 100
            CASE DEFAULT
             WRITE(6,*) '@D2PMCO : UNKNOWN TYPE (',PARFMT(i),') FOR',
     >       ' PKEY (',PKEY(i),').'
            CALL XABORT('')
            END SELECT
            WRITE(6,*) '@D2PMCO : INCONSITENT TYPE FOR',
     >       ' PKEY (',PKEY(i),'), TYPE (',PARFMT(i),') EXPECTED.'
            CALL XABORT ('')
  100      RKOTH(j)=i
           EXIT
           ENDIF
          ENDDO
        ENDIF
          RANK_INDEX(i)=0
      ENDDO
      RANK_INDEX(NPAR+1)=0

      ! DETERMINE ODER_VAL ARRAY
      IF(LBARR) THEN
        ORDER_VAL(1)=1
      ELSE
        NCRD_SAP=1
        IF(NCRD>1) THEN
          WRITE(6,*) "@D2PMCO:",
     1     " CONTROL ROD STATE VARIABLE IS MISSING IN SAPHYB"
          CALL XABORT("=> NUMBER OF CTRL ROD VALUE MUST BE SET TO 1")
        ELSE IF(CRDINF(1).NE. 1) THEN
          WRITE(6,*) "@D2PMCO:",
     1     " CONTROL ROD STATE VARIABLE IS MISSING IN SAPHYB"
          CALL XABORT("=> CTRL ROD UNRODDED INDEX MUST BE SET TO 1")
        ENDIF
      ENDIF
      IF(LDMOD) THEN
         ORDER_VAL(2)=1
         IF(LBARR) ORDER_VAL(2)=2
      ENDIF
      IF(LCBOR) THEN
         IF(LDMOD) THEN
          ORDER_VAL(3)=ORDER_VAL(2)+1
         ELSE IF(LBARR) THEN
          ORDER_VAL(3)=2
         ELSE
          ORDER_VAL(3)=1
         ENDIF
      ENDIF
      IF(LTCOM) THEN
         IF(LCBOR) THEN
          ORDER_VAL(4)=ORDER_VAL(3)+1
         ELSE IF(LDMOD) THEN
          ORDER_VAL(4)=ORDER_VAL(2)+1
         ELSE IF(LBARR) THEN
          ORDER_VAL(4)=2
         ELSE
          ORDER_VAL(4)=1
         ENDIF
      ENDIF
      IF(LTMOD) THEN
         IF(LTCOM) THEN
          ORDER_VAL(5)=ORDER_VAL(4)+1
         ELSE IF(LCBOR) THEN
          ORDER_VAL(5)=ORDER_VAL(3)+1
         ELSE IF(LDMOD) THEN
          ORDER_VAL(5)=ORDER_VAL(2)+1
         ELSE IF(LBARR) THEN
          ORDER_VAL(5)=2
         ELSE
          ORDER_VAL(5)=1
         ENDIF
      ENDIF
      ! STORE THE NAME OF CURENT PKEY IN PKEY_TMP
      DO i=1, NPAR
        PKEY_TMP(i)=PKEY(i)
      ENDDO

      IF(.NOT.LBURN) THEN
        WRITE(6,*)
        WRITE(6,*)('WARNING: BURN VARIABLE IS MISSING IN MCO')
        WRITE(6,*)('=> 0 MWJ/T SINGLE EXPOSURE ASSUMED')
        WRITE(6,*)
        DEALLOCATE (PKEY,NVAL)
        NPAR=NPAR+1
        ALLOCATE (PKEY(NPAR),NVAL(NPAR))
        DO i=1, NPAR-1
          PKEY(i)=PKEY_TMP(i)
          NVAL(i)=NVALTMP(i)
        ENDDO
        PKEY(NPAR)="BURN"
        NVAL(NPAR)=1
        DEALLOCATE (PKEY_TMP)
        ALLOCATE(PKEY_TMP (NPAR))
        PKEY_TMP=PKEY
      ENDIF
      IF(LTMOD) THEN
        ORDER_VAL(6)=ORDER_VAL(5)+1
      ELSE IF(LTCOM) THEN
        ORDER_VAL(6)=ORDER_VAL(4)+1
      ELSE IF(LCBOR) THEN
        ORDER_VAL(6)=ORDER_VAL(3)+1
      ELSE IF(LDMOD) THEN
        ORDER_VAL(6)=ORDER_VAL(2)+1
      ELSE IF(LBARR) THEN
        ORDER_VAL(6)=2
      ELSE
        ORDER_VAL(6)=1
      ENDIF

      ALLOCATE (PVALDIR(NPAR))
      ALLOCATE(VALPAR(NPAR,100))

      OTHR(:,:)=0.
      OTHC(:,:)=''

      DO i=1, NPAR
        ! NAME OF DIRECTORY IN SAPHYB CONTAINING VALUES OF STATE
        ! VARIABLES : PKEY(I)
        IF ((PARFMT(i).NE.'STRING')) THEN
         IF ((PKEY(i).NE.PKNAM(6))) THEN

          WRITE(PVALDIR(i),'("pval", I8.8)') i
         ! STORE VALUES IN VALPAR
          CALL LCMGET(IPDIR,PVALDIR(i),VALPAR(i,1:NVAL(i)))

         ELSE IF(LBURN) THEN

          WRITE(PVALDIR(i),'("pval", I8.8)') i
         ! STORE VALUES IN VALPAR LBURN
          CALL LCMGET(IPDIR,PVALDIR(i),VALPAR(i,1:NVAL(i)))
         ELSE
         ! STORE VALUES IN VALPAR
          VALPAR(i,1:NVAL(i))=0.0
         ENDIF
        ENDIF

        DO j=1,NOTH
         IF (LOTH(j).EQV. .FALSE.) THEN
          WRITE(6,*) '@D2PMCO: UNKNOWN PKEY (',OTHPK(j),') IN MCO'
          CALL XABORT ('=> PLEASE CHECK MCO CONTENT')
         ELSE IF (PKEY(i).EQ.OTHPK(j)) THEN
          WRITE(PVALDIR(i),'("pval", I8.8)') i
          IF (OTHTYP(j).EQ.3) THEN
           CALL LCMGTC(IPDIR,PVALDIR(i),12,NVAL(i),OTHC(i,1:NVAL(i)))
           DO k=1, NVAL(i)
            IF (OTHC(i,k).EQ.OTHVAL(j)) THEN
             OTHVAC(j)=OTHVAL(j)
             EXIT
            ENDIF
            IF (k.EQ.NVAL(i)) THEN
             WRITE (6,*) '@D2PMCO: VALUE (',OTHVAL(j),') FOR PKEY(',
     >       PKEY(i),') IS OUT OF RANGE'
             WRITE (6,*) '=> POSSIBLE VALUES ARE :'
             WRITE (6,'(A12,1X)') OTHC(i,1:NVAL(i))
             CALL XABORT ("")
            ENDIF
           ENDDO
          ELSE
           CALL LCMGET(IPDIR,PVALDIR(i),OTHR(i,1:NVAL(i)))
           DO k=1, NVAL(i)
            WRITE(OTHVALC,'(f12.5)')OTHR(i,k)
            IF (OTHVALC.EQ.OTHVAL(j)) THEN
             OTHVAR(j)=OTHR(i,k)
             EXIT
            ENDIF
            IF (k.EQ.NVAL(i)) THEN
             OTHVAR(j)=OTHREA(j)
             WRITE (6,*) 'WARNING : VALUE (',OTHVAL(j),') FOR PKEY(',
     >       PKEY(i),') IS OUT OF RANGE'
             WRITE (6,*) '=> POSSIBLE VALUES ARE :'
             WRITE (6,'(e12.5,1X)') OTHR(i,1:NVAL(i))
             WRITE (6,*) '=>INTERPOLATION WILL BE NEEDED'
            ENDIF
           ENDDO

          ENDIF
         ENDIF
        ENDDO


        ! CASE OF CONTROL ROD
        IF(PKEY(i)==PKNAM(1)) THEN
           RANK(i)=ORDER_VAL(1);
           RANK_INDEX(ORDER_VAL(1))=i

           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(1)) THEN
              WRITE(6,*)('@D2PMCO: IMPOSSIBLE TO ADD A CONTROL ')
              CALL XABORT ('ROD VALUE IN THE PMAXS TREE')
             ENDIF
            ENDDO
           ENDIF
        ! CASE OF MODERATOR DENSITY
        ELSE IF(PKEY(i)==PKNAM(2)) THEN
           RANK(i)=ORDER_VAL(2)
           RANK_INDEX(ORDER_VAL(2))=i

           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(2)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              NVAL(i)=0
              ENDIF
              VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >         USRVAPK(UV,1:USRVAL(UV))
              NVAL(i)=NVAL(i)+USRVAL(UV)
              CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
             ENDIF
            ENDDO
           ENDIF
        ! CASE OF BORON CONCENTRATION
        ELSE IF(PKEY(i)==PKNAM(3)) THEN
           RANK(i)=ORDER_VAL(3)
           RANK_INDEX(ORDER_VAL(3))=i
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(3)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              NVAL(i)=0
              ENDIF
              VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >         USRVAPK(UV,1:USRVAL(UV))
              NVAL(i)=NVAL(i)+USRVAL(UV)
              CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
             ENDIF
            ENDDO
           ENDIF
        ! CASE OF FUEL TEMPERATURE
        ELSE IF(PKEY(i)==PKNAM(4)) THEN
           RANK(i)=ORDER_VAL(4)
           RANK_INDEX(ORDER_VAL(4))=i

           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(4)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              NVAL(i)=0
              ENDIF
              VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >         USRVAPK(UV,1:USRVAL(UV))
              NVAL(i)=NVAL(i)+USRVAL(UV)
              CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
             ENDIF
            ENDDO
           ENDIF
        ! CASE OF MODERATOR DENSITY
        ELSE IF(PKEY(i)==PKNAM(5)) THEN
           RANK(i)=ORDER_VAL(5)
           RANK_INDEX(ORDER_VAL(5))=i
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(5)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              NVAL(i)=0
              ENDIF
              VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >         USRVAPK(UV,1:USRVAL(UV))
              NVAL(i)=NVAL(i)+USRVAL(UV)
              CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
             ENDIF
            ENDDO
           ENDIF
        ! CASE OF BURN UP
        ELSE IF(PKEY(i)==PKNAM(6))  THEN
           RANK(i)=NPAR
           RANK_INDEX(NPAR)=i
           STAVEC(4)=NVAL(i)
           IF(LADD) THEN
             DO UV=1,USRSTA
              IF(USRPAR(UV)==PKNAM(6)) THEN
               IF(LNEW) THEN
               VALPAR(i,1:NVAL(i))=0.0
               NVAL(i)=0
               ENDIF
               VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >          USRVAPK(UV,1:USRVAL(UV))
               NVAL(i)=NVAL(i)+USRVAL(UV)
               STAVEC(4)=NVAL(i)
               CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
              ENDIF
             ENDDO
           ENDIF
        ELSE
           NOTHTH=NPAR-MAXVAL(ORDER_VAL)
           IF((PKEY(i)=='FLUE').OR.(PKEY(i)=='TIME')) NOTHTH=NOTHTH-1
           RANK(i) = NPAR+i
           RANK_INDEX(NPAR+1)=NPAR+1
        END IF
      ENDDO

      ! D2PSOR STATE VARIABLE INPUT TO MATCH GENPMAXS ORDER
      CALL D2PSOI(RANK,NPAR)

      ! LOOP OVER STATES VARIABLES IN SAPHYB
      DO i=1, NPAR
        ! WE KEEP ONLY "REAL" STATES VARIABLE (IE EXEPT FLUE, TIME ETC.
        IF(RANK(i)<=NPAR) THEN
          ! RESTORE THE NAME OK PKEY AFTER THE CALL TO D2PSOR SUBROUTINE
          PKEY(i)=PKEY_TMP(RANK_INDEX(RANK(i)))
          NSVAR = NSVAR + 1
        ENDIF
      ENDDO

      ! CREATION OF THE SAPHYB_INFO DIRECTORY INTO THE INFO DATA BLOCK
      STAVEC(2) = NSVAR  ! NVAR
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      ! CREATION OF : INFO/SAPHYB_INFO/STATE_VAR
      CALL LCMPUT(IPDAT,'NTOT',1,1,NTOT)
      CALL LCMPTC(IPDAT,'NAMDIR',12,1,MIXDIR)
      CALL LCMPTC(IPDAT,'STATE_VAR',12,NSVAR,PKEY)
      IF(.NOT.(LBARR)) THEN
        PKEY_BARR(1)="BARR        "
        DO j=1, NSVAR
          PKEY_BARR(j+1)=PKEY(j)
        ENDDO
        ! CREATION OF : INFO/SAPHYB_INFO/STATE_VAR
        CALL LCMPTC(IPDAT,'STATE_VAR',12,NSVAR+1,PKEY_BARR)
        ! CREATION OF : INFO/SAPHYB_INFO/BARR
        CALL LCMPUT(IPDAT,'BARR',1,2,1.0)
        STAVEC(2) = NSVAR + 1 ! NVAR
!        NSVAR=NSVAR+1
      ENDIF

       ALLOCATE (PKIDX(STAVEC(2)))
       CALL XDISET(PKIDX,STAVEC(2),0)
       IF (.NOT. LBARR) PKIDX(STAVEC(2))= -1
      DO i=1, NSVAR
        DO j=2,6
           IF(PKEY(i)==PKNAM(j)) THEN
           PKIDX(i)=j
           ENDIF
        ENDDO
        IF(PKEY(i)==PKNAM(1)) THEN
           IF (LBARR) PKIDX(i)=1
           NCRD_SAP=NVAL(RANK_INDEX(RANK(i)))
           ! REORGANIZATION OF BARR PARAMETERS TO MATCH GENPMAXS
           ! FORMALISM. SPECIAL TREATMENT FOR BARR PARAMETERS TO TAKE
           ! INTO ACCOUNT THE MEANING OF BARR VALUES
           IF(NCRD.NE.NCRD_SAP) THEN
            WRITE(6,*) "@D2PMCO: ERROR IN CONTROL ROD COMPOSITION "
            WRITE(6,*) "THE NUMBER OF CONTROL ROD COMPOSITIONS IN ",
     1      "SAP (",NCRD_SAP,") IS DIFFERENT FROM D2P INPUT (",NCRD,")"
            WRITE(6,*) "SAP :",VALPAR(RANK_INDEX(RANK(i)),1:NCRD_SAP)
            WRITE(6,*) "D2P INPUT   :",CRDINF(1:5)
            CALL XABORT('D2PMCO: INPUT ERROR')
           ENDIF
           CALL D2PREO(IPDAT,VALPAR,RANK_INDEX(RANK(i)),NPAR,
     1     NVAL(RANK_INDEX(RANK(i))),IPRINT)
        ENDIF

        IF(MESH.EQ.'GLOB') THEN
            CALL LCMPUT(IPDAT,PKREF(PKIDX(i)),
     1       NVAL(RANK_INDEX(RANK(i))),2,VALPAR(RANK_INDEX(RANK(i)),
     2       1:NVAL(RANK_INDEX(RANK(i)))))
          DO l=1,USRSTA
           IF(USRPAR(l)==PKEY(i)) THEN
            IF(PKEY(i) =='BARR') THEN
             CALL XABORT('@D2PMCO: THE CR STATE CANNOT BE SET BY USER')
            ENDIF
            IF((USRVAL(l)>1).and.NVAL(RANK_INDEX(RANK(i)))==1) THEN
             WRITE(6,*)"@D2PMCO: IMPOSSIBLE TO DEFINE USER MESHING",
     1        " FOR ",PKEY(i)
             CALL XABORT ('ONLY ONE VALUE IS CONTAINED IN THE MCO')
            ENDIF

            FIRST_VAL=VALPAR(RANK_INDEX(RANK(i)),1)
            LAST_VAL=NVAL(RANK_INDEX(RANK(i)))
            LAST_VAL=VALPAR(RANK_INDEX(RANK(i)),INT(LAST_VAL))
            NVAL(RANK_INDEX(RANK(i))) = USRVAL(l)
            IF(USRVAL(l)>1) THEN
             PITCH = (LAST_VAL-FIRST_VAL)/(USRVAL(l)-1)

             DO n=1,USRVAL(l)
              VALPAR(RANK_INDEX(RANK(i)),n)=FIRST_VAL+PITCH*(n-1)
             ENDDO
            ELSE
             VALPAR(RANK_INDEX(RANK(i)),1)=(FIRST_VAL+LAST_VAL)/2.0
            ENDIF

            CALL  LCMPUT(IPDAT,PKREF(PKIDX(i)),USRVAL(l),2,
     1      VALPAR(RANK_INDEX(RANK(i)),1:USRVAL(l)))
           ENDIF
          ENDDO
        ELSE
          ! CREATION OF: INFO/SAPHYB_INFO/SVNAME

          CALL  LCMPUT(IPDAT,PKREF(PKIDX(i)),
     1     NVAL(RANK_INDEX(RANK(i))),2,VALPAR(RANK_INDEX(RANK(i)),
     2     1:NVAL(RANK_INDEX(RANK(i)))) )
        ENDIF
       ENDDO

      CALL LCMPUT(IPDAT,'PKIDX',STAVEC(2),1,PKIDX)
      IF (NOTH>0) THEN
       CALL LCMPTC(IPDAT,'OTHPK',12,NOTH,OTHPK)
       CALL LCMPUT(IPDAT,'OTHTYP',NOTH,1,OTHTYP)
       CALL LCMPTC(IPDAT,'OTHVAC',12,NOTH,OTHVAC)
       CALL LCMPUT(IPDAT,'OTHVAR',NOTH,2,OTHVAR)
      ENDIF
      IF(MESH=='DEF') THEN
        STAVEC(5) = 0
      ELSE IF(MESH=='SAP') THEN
        STAVEC(5) = 1
      ELSE IF(MESH=='GLOB') THEN
        STAVEC(5) = 2
      ELSE IF(MESH=='ADD') THEN
        STAVEC(5) = 3
        IF(LNEW) STAVEC(5) = 4
      ENDIF
      IF (LADF) THEN
       CALL LCMPTC(IPDAT,'ADF_TYPE',3,1,ADF)
       IF (ADF.EQ.'DRA') THEN
         CALL LCMPTC(IPDAT,'HADF',8,STAVEC(13),ADFD)
       ELSE IF (ADF.EQ.'GEN') THEN
         CALL LCMPTC(IPDAT,'HFLX',8,2,HFLX)
         CALL LCMPTC(IPDAT,'HCUR',8,2,HCUR)
         CALL LCMPUT(IPDAT,'THCK',2,1,THCK)
       ENDIF
      ENDIF

      IF (LCDF) THEN
       CALL LCMPTC(IPDAT,'CDF_TYPE',3,1,CDF)
       CALL LCMPTC(IPDAT,'HCDF',8,STAVEC(15),CDFD)
      ENDIF

      IF (LGFF) CALL LCMPTC(IPDAT,'GFF_TYPE',3,1,GFF)

      IF (LYLD) THEN
       CALL LCMPTC(IPDAT,'YLD_OPT',3,1,YLDOPT)
       CALL LCMPUT(IPDAT,'YLD_FIX',3,2,YLD)
       CALL LCMPUT(IPDAT,'YLD_LOC',5,2,LOCYLD)
      ENDIF
      CALL LCMPUT(IPDAT,'LABS', 3,5,LABS)
      CALL LCMPUT(IPDAT,'SCAT', 1,5,SCAT)
      CALL LCMSIX(IPDAT,'ISOTOPES',1)
      CALL LCMPTC(IPDAT,'XE135',12,1,ISOT(1))
      CALL LCMPTC(IPDAT,'SM149',12,1,ISOT(2))
      CALL LCMPTC(IPDAT,'I135',12,1,ISOT(3))
      CALL LCMPTC(IPDAT,'PM149',12,1,ISOT(4))
      CALL LCMPTC(IPDAT,'PM148',12,1,ISOT(5))
      CALL LCMPTC(IPDAT,'PM148M',12,1,ISOT(6))
      CALL LCMPTC(IPDAT,'ND147',12,1,ISOT(7))
      CALL LCMPTC(IPDAT,'PM147',12,1,ISOT(8))

      ! SET THE IPDAT/STAVEC
      STAVEC(1) = DIMMCO(2)                    ! NGROUP
      STAVEC(3) = N_XS                         ! N_XS
      STAVEC(6) = NCRD                         ! NCOMPO

      STAVEC(7) = NDEL                         ! NDLAY

      STAVEC(16)= NGFF                         ! GFF(NGFF,NG)
      STAVEC(17)= NPIN                         ! GFFP(NPIN,NPIN,NG)

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMPUT(IPDAT,'STATE-VECTOR',40,1,STAVEC)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMPUT(IPDAT,'MIX',1,1,MIX)

      IPTH=LCMLID(IPDAT,'PKEY_INFO',6)
      DO J=1, 6
         KPTH=LCMDIL(IPTH,J)
         IF(J==1) THEN
           CALL LCMPTC(KPTH,'NAME',12,1,PKNAM(1))
           CALL LCMPUT(KPTH,'LFLAG',1,5,LBARR)
         ELSE IF(J==2)THEN
           IF(LDMOD) CALL LCMPTC(KPTH,'NAME',12,1,PKNAM(2))
           CALL LCMPUT(KPTH,'LFLAG',1,5,LDMOD)
         ELSE IF(J==3) THEN
           IF(LCBOR) CALL LCMPTC(KPTH,'NAME',12,1,PKNAM(3))
           CALL LCMPUT(KPTH,'LFLAG',1,5,LCBOR)
         ELSE IF(J==4)THEN
           IF(LTCOM) CALL LCMPTC(KPTH,'NAME',12,1,PKNAM(4))
           CALL LCMPUT(KPTH,'LFLAG',1,5,LTCOM)
         ELSE IF(J==5)THEN
           IF(LTMOD) CALL LCMPTC(KPTH,'NAME',12,1,PKNAM(5))
           CALL LCMPUT(KPTH,'LFLAG',1,5,LTMOD)
         ELSE IF(J==6) THEN
           CALL LCMPTC(KPTH,'NAME',12,1,PKNAM(6))
           CALL LCMPUT(KPTH,'LFLAG',1,5,LBURN)
         ENDIF
      ENDDO
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'HELIOS_HEAD',1)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMSIX(IPDAT,' ',0)

      ! EDIT THE LISTING FILE
      IF(IPRINT > 0)  THEN
         WRITE(6,*) "******* CONTENT OF MULTICOMPO RECOVERED **********"
         WRITE(6,*) " DIRECTORY NAME   : ", MIXDIR
         WRITE(6,*) " INDEX OF MIXTURE : ", MIX
         WRITE(6,*)
         WRITE(6,*) "******* CONTENT OF MULTICOMPO RECOVERED **********"
         WRITE(6,*)
         WRITE(6,*)"NB OF STATE VARIBALE IN MCO      :", NPAR
         WRITE(6,*)"NB OF STATE VARIABLES RECOGNIZED :", NSVAR
         WRITE(6,*)"NAME OF STATE VARIABLES IN MCO   :", PKEY_TMP
         WRITE(6,*)"RECOGNIZED STATE VARIABLES       :",PKEY(1:NSVAR)
         IF (NOTH.GE.1) THEN
         WRITE(6,*)"OTHER STATE VARIABLES            :",OTHPK(1:NOTH)
         WRITE(6,*)"OTHER STATE  VALUES              :",OTHVAL(1:NOTH)
         ENDIF
         IF(NOTHTH.NE.NOTH) THEN
         WRITE(6,*) "=> WARNING: UNRECOGNIZED VARIABLES !"
         WRITE(6,*) "==>PLEASE USE THE PKEY CARD OF D2P: MODULE"
         CALL XABORT("")
         ENDIF
         WRITE(6,*) "FLAG FOR STATE VARIABLES : "
         WRITE(6,*) "   CONTROL ROD           : ", LBARR
         WRITE(6,*) "   MODERATOR DENSITY     : ", LDMOD
         WRITE(6,*) "   BORON CONCENTRATION   : ", LCBOR
         WRITE(6,*) "   FUEL TEMPERATURE      : ", LTCOM
         WRITE(6,*) "   MODERATOR TEMPERATURE : ", LTMOD
         WRITE(6,*) "   BURNUP                : ", LBURN
         WRITE(6,*) "ASSEMBLY DISCONTINUITY FACTORS           : ", LADF
         IF(LADF) THEN
          IF(ADF .EQ. 'DRA')  WRITE(6,*) "TYPE OF ADF : DRAGON"
          IF(ADF .EQ. 'GET')  WRITE(6,*) "TYPE OF ADF : GET"
          IF(ADF .EQ. 'SEL')  WRITE(6,*) "TYPE OF ADF : SELENGUT"
          IF(ADF .EQ. 'GEN')  WRITE(6,*) "TYPE OF ADF : GENPMAXS"
         ENDIF
         IF (STAVEC(21).EQ.1) THEN
          WRITE(6,*)'WARNING => ADF ARE INTEGRATED IN CROSS SECTIONS'
         ENDIF
         WRITE(6,*) "CORNER DISCONTINUITY FACTORS             : ", LCDF
         IF(LCDF) THEN
          IF(CDF .EQ. 'DRA')  WRITE(6,*) "TYPE OF CDF : DRAGON"
         ENDIF
         WRITE(6,*) "GROUP FORM FACTORS                       : ", LGFF
         IF(LGFF) THEN
          IF(GFF .EQ. 'DRA')  WRITE(6,*) "TYPE OF GFF : DRAGON"
         ENDIF
         WRITE(6,*) "ABSORPTION TYPE          : "
         WRITE(6,*) "   SAP                   : ", SAP
         WRITE(6,*) "   MIC                   : ", MIC
         WRITE(6,*) "   EXC                   : ", EXC

         WRITE(6,*)
         DO i=1, NSVAR
          WRITE(6,*) "NUMBER OF VALUES FOR ",PKEY(i)," PARAMETER :",
     1    NVAL(RANK_INDEX(RANK(i)))
          WRITE(6,*) "VALUES FOR ",PKEY(i)," PARAMETER :",
     1    VALPAR(RANK_INDEX(RANK(i)),1:NVAL(RANK_INDEX(RANK(i))))
          WRITE(6,*)
         ENDDO
        WRITE(6,*)
        WRITE(6,*) "NAME OF FISSION PRODUCTS FOR FISSION YIELD :"
        WRITE(6,*) "XE135  : ",ISOT(1)
        WRITE(6,*) "SM149  : ",ISOT(2)
        WRITE(6,*) "I135   : ",ISOT(3)
        WRITE(6,*) "PM149  : ",ISOT(4)
        WRITE(6,*) "PM148  : ",ISOT(5)
        WRITE(6,*) "PM148M : ",ISOT(6)
        WRITE(6,*) "ND147  : ",ISOT(7)
        WRITE(6,*) "PM147  : ",ISOT(8)
        WRITE(6,*)
        IF (LYLD) THEN
         WRITE(6,*) "OPTION FOR FISSION YIELD RECOVERY: ",YLDOPT
         IF (STAVEC(22)>0) THEN
         WRITE(6,*)"CORRECTION FOR SAMARIUM PRODUCTION IS APPLIED"
         ENDIF
          IF (YLDOPT.EQ.'MAN')THEN
           WRITE(6,*)"LOCAL CONDITIONS SET BY THE USER :"
           DO I=1,5
            IF (LOCYLD(I).NE.-1) THEN
             WRITE(6,*) PKNAM(I)," = ",LOCYLD(I)
            ENDIF
           ENDDO
          ENDIF
        ENDIF

        WRITE(6,*)
      ENDIF
      ! free memory
      DEALLOCATE (PKEY)
      DEALLOCATE (PKIDX)
      DEALLOCATE (NVAL)
      DEALLOCATE (PVALDIR)
      DEALLOCATE (PKEY_TMP)
      DEALLOCATE (RANK)
      DEALLOCATE (RANK_INDEX)
      DEALLOCATE (VALPAR)
      DEALLOCATE (PARFMT)
      RETURN
      END
