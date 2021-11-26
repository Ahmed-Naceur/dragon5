*DECK D2PSAP
      SUBROUTINE D2PSAP(  IPSAP,  IPDAT, STAVEC, CRDINF,   NCRD,  PKNAM,
     >                    ISOT ,   MESH, USRPAR, USRVAL, USRSTA,USRVAPK,
     >                    SAP  ,    MIC,    EXC,   SCAT,    ADF,  LADD ,
     >                    LNEW ,   LADF, IPRINT,   LYLD,    YLD, YLDOPT,
     >                   LOCYLD,   HDET                                )
*
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
* LYLD    Fission Yield must be recovered
* YLD     user defined values for fission yields (1:I, 2:XE, 3:PM)
* LOCYLD  value for state parameter on which fission yield will be
*         calculated
* YLDOPT  option for fission yields calculation (DEF, MAN, FIX)
* HDET    name of isotope for the detector cross sections
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPDAT
      INTEGER NCRD,USRSTA
      INTEGER IPRINT
      INTEGER STAVEC(40),CRDINF(20),USRVAL(12)
      REAL USRVAPK(12,10),YLD(3),LOCYLD(5)
      CHARACTER*3 ADF,YLDOPT
      CHARACTER*12 USRPAR(12)
      CHARACTER*5 MESH
      CHARACTER*12 PKNAM(6)
      CHARACTER*12 ISOT(8)
      LOGICAL SAP, MIC, EXC,SCAT,LADD,LNEW,LADF,LYLD
      CHARACTER*12 HDET
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR)  IPROOT,IPTH,KPTH

      PARAMETER(NDIMSAP=50)
      INTEGER :: N_XS = 8
      INTEGER,DIMENSION(6)  :: ORDER_VAL = 0
      INTEGER DIMSAP(NDIMSAP)
      INTEGER NPAR,NCALS,NSVAR,NBREA,ITYLCM,VALTOT
      INTEGER NCRD_SAP,NVALTMP(10)
      INTEGER i, j, k, l , n, UV,ILONG
      INTEGER :: NTOT = 0
      REAL FIRST_VAL,LAST_VAL,PITCH
      LOGICAL LABS(3)
      LOGICAL :: LBARR = .FALSE.
      LOGICAL :: LDMOD = .FALSE.
      LOGICAL :: LCBOR = .FALSE.
      LOGICAL :: LTCOM = .FALSE.
      LOGICAL :: LTMOD = .FALSE.
      LOGICAL :: LBURN = .FALSE.
      CHARACTER(LEN=12)   PKEY_BARR(6)
      CHARACTER*12,DIMENSION(6) :: PKREF
      DATA PKREF/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/

      INTEGER, ALLOCATABLE, DIMENSION(:) :: NVAL,RANK,RANK_INDEX,PKIDX
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PKEY,PKEY_TMP
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PVALDIR, NOMREA
      REAL, ALLOCATABLE, DIMENSION(:) :: SV_VAL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VALPAR

      IPROOT=IPSAP
      LABS(1)=MIC
      LABS(2)=SAP
      LABS(3)=EXC
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMPUT(IPDAT,'BARR_INFO',NCRD,1,CRDINF)
      ! RECOVER DIMSAP INFORMATION FROM SAPHYB
      CALL XDISET(DIMSAP,NDIMSAP,0)
      CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)

      NPAR =DIMSAP(8)
      NMIL =DIMSAP(7)
      NREA =DIMSAP(4)
      NISO =DIMSAP(5)
      NMAC =DIMSAP(6)
      ! INITIALIZATION OF PARAMETERS
      VALTOT = 0
      NSVAR = 0
      k = 1

      ! MEMORY ALLOCATION
      ALLOCATE (PKEY(NPAR),NVAL(NPAR),RANK(NPAR))
      ALLOCATE (PKEY_TMP(NPAR),RANK_INDEX(NPAR+1))
      CALL LCMSIX (IPSAP,' ',0)
      CALL LCMSIX (IPSAP,'paramarbre',1)
      CALL LCMGET (IPSAP,'NCALS',NCALS)
      CALL LCMSIX (IPSAP,' ',0)
      CALL LCMSIX (IPSAP,'contenu',1)
      CALL LCMLEN(IPSAP,'NOMREA',NBREA,ITYLCM)
      ALLOCATE (NOMREA(NBREA))
      CALL LCMGTC(IPSAP,'NOMREA',12,NBREA,NOMREA)
      CALL LCMSIX(IPSAP,' ',0)
      CALL LCMSIX(IPSAP,'paramdescrip',1)
      CALL LCMLEN(IPSAP,'PARKEY',ILONG,ITYLCM)
      CALL LCMGTC(IPSAP,'PARKEY',4,NPAR,PKEY)
      CALL LCMGET(IPSAP,'NVALUE',NVAL)
      NVALTMP=NVAL
      CALL LCMSIX(IPSAP,' ',0)
      CALL LCMSIX(IPSAP,'paramvaleurs',1)
      ! LOOP OVER STATE VARIABLES OF SAPHYB
      ! CHECK OF EXISTENCE OF STATE PARAMETER
      PKEY (1:NPAR) (5:12) = "        "

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
          WRITE(6,*) "@D2PSAP:",
     1     " CONTROL ROD STATE VARIABLE IS MISSING IN SAPHYB"
          CALL XABORT("=> NUMBER OF CTRL ROD VALUE MUST BE SET TO 1")
        ELSE IF(CRDINF(1).NE. 1) THEN
          WRITE(6,*) "@D2PSAP:",
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
      DO i=1, NPAR
        PKEY_TMP(i)=PKEY(i)
      ENDDO

      IF(.NOT.LBURN) THEN
        WRITE(6,*)
        WRITE(6,*)('WARNING: BURN VARIABLE IS MISSING IN MCO')
        WRITE(6,*)('=> 0 MWJ/T SINGLE EXPOSURE IS ASSUMED')
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

      ALLOCATE (PVALDIR(NPAR),VALPAR(NPAR,100))

      DO i=1, NPAR
        ! NAME OF DIRECTORY IN SAPHYB CONTAINING VALUES OF STATE
        IF ((PKEY(i).NE.PKNAM(6))) THEN
         WRITE(PVALDIR(i),'("pval", I8)') i
         CALL LCMGET(IPSAP,PVALDIR(i),VALPAR(i,1:NVAL(i)))
        ELSE IF(LBURN) THEN
         WRITE(PVALDIR(i),'("pval", I8)') i
         CALL LCMGET(IPSAP,PVALDIR(i),VALPAR(i,1:NVAL(i)))
        ELSE
         VALPAR(i,1:NVAL(i))=0.0
        ENDIF
        ! CASE OF CONTROL ROD
        IF(PKEY(i)==PKNAM(1)) THEN
           RANK(i)=ORDER_VAL(1);
           RANK_INDEX(ORDER_VAL(1))=i
           VALTOT=VALTOT+NVAL(i);
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(1)) THEN
              WRITE(6,*)('@D2PSAP: IMPOSSIBLE TO ADD A CONTROL ')
              CALL XABORT ('ROD VALUE IN THE PMAXS TREE')
             ENDIF
            ENDDO
           ENDIF
        ! CASE OF MODERATOR DENSITY
        ELSE IF(PKEY(i)==PKNAM(2)) THEN
           RANK(i)=ORDER_VAL(2)
           RANK_INDEX(ORDER_VAL(2))=i
           VALTOT=VALTOT+NVAL(i);
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(2)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              VALTOT=VALTOT-NVAL(i);
              NVAL(i)=0
              ENDIF
              VALTOT=VALTOT+USRVAL(UV)
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
           VALTOT=VALTOT+NVAL(i);
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(3)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              VALTOT=VALTOT-NVAL(i);
              NVAL(i)=0
              ENDIF
              VALTOT=VALTOT+USRVAL(UV)
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
           VALTOT=VALTOT+NVAL(i);
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(4)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              VALTOT=VALTOT-NVAL(i);
              NVAL(i)=0
              ENDIF
              VALTOT=VALTOT+USRVAL(UV)
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
           VALTOT=VALTOT+NVAL(i);
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(5)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              VALTOT=VALTOT-NVAL(i);
              NVAL(i)=0
              ENDIF
              VALTOT=VALTOT+USRVAL(UV)
              VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >         USRVAPK(UV,1:USRVAL(UV))
              NVAL(i)=NVAL(i)+USRVAL(UV)
              CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
             ENDIF
            ENDDO
           ENDIF
        ELSE IF(PKEY(i)==PKNAM(6))  THEN
           RANK(i)=NPAR
           RANK_INDEX(NPAR)=i
           VALTOT=VALTOT+NVAL(i)
           STAVEC(4)=NVAL(i)
           IF(LADD) THEN
            DO UV=1,USRSTA
             IF(USRPAR(UV)==PKNAM(6)) THEN
              IF(LNEW) THEN
              VALPAR(i,1:NVAL(i))=0.0
              VALTOT=VALTOT-NVAL(i);
              NVAL(i)=0
               ENDIF
               VALTOT=VALTOT+USRVAL(UV)
               VALPAR(i,NVAL(i)+1:NVAL(i)+1+USRVAL(UV))=
     >          USRVAPK(UV,1:USRVAL(UV))
               NVAL(i)=NVAL(i)+USRVAL(UV)
               STAVEC(4)=NVAL(i)
               CALL D2PSOR(VALPAR(i,1:NVAL(i)),NVAL(i))
              ENDIF
             ENDDO
           ENDIF
        ELSE
           RANK(i) = NPAR+i
           RANK_INDEX(NPAR+1)=NPAR+1
        END IF
      ENDDO

      ALLOCATE (SV_VAL(VALTOT))
      ! D2PSOR STATE VARIABLE INPUT TO MATCH GENPMAXS ORDER
      CALL D2PSOI(RANK,NPAR)

      ! LOOP OVER STATES VARIABLES IN SAPHYB
      DO i=1, NPAR
        ! WE KEEP ONLY "REAL" STATES VARIABLE (IE EXEPT FLUE, TIME ETC.
        IF(RANK(i)<=NPAR) THEN
          ! RESTORE THE NAME OK PKEY AFTER THE CALL TO D2PSOR SUBROUTINE
          PKEY(i)=PKEY_TMP(RANK_INDEX(RANK(i)))
          NSVAR = NSVAR + 1
          DO j=1, NVAL(RANK_INDEX(RANK(i)))
           SV_VAL(k)=VALPAR(RANK_INDEX(RANK(i)),j)
           k=k+1
          ENDDO
        ENDIF
      ENDDO

      ! CREATION OF THE SAPHYB_INFO DIRECTORY INTO THE INFO DATA BLOCK
      STAVEC(2) = NSVAR  ! NVAR
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      ! CREATION OF : INFO/SAPHYB_INFO/STATE_VAR
      CALL LCMPUT(IPDAT,'NTOT',1,1,NTOT)
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
            WRITE(6,*) "@D2PSAP: ERROR IN CONTROL ROD COMPOSITION "
            WRITE(6,*) "THE NUMBER OF CONTROL ROD COMPOSITIONS IN ",
     1      "SAP (",NCRD_SAP,") IS DIFFERENT FROM D2P INPUT (",NCRD,")"
            WRITE(6,*) "SAP :",VALPAR(RANK_INDEX(RANK(i)),1:NCRD_SAP)
            WRITE(6,*) "D2P INPUT   :",CRDINF(1:5)
            CALL XABORT('')
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
             CALL XABORT('@D2PSAP: THE CR STATE CANNOT BE SET BY USER')
            ENDIF
            IF((USRVAL(l)>1).and.NVAL(RANK_INDEX(RANK(i)))==1) THEN
             WRITE(6,*)"@D2PSAP: IMPOSSIBLE TO DEFINE USER MESHING",
     1        " FOR ",PKEY(i)
             CALL XABORT ('ONLY ONE VALUE IS CONTAINED IN THE SAPHYB')
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
      IF (LYLD) THEN
       CALL LCMPTC(IPDAT,'YLD_OPT',3,1,YLDOPT)
       CALL LCMPUT(IPDAT,'YLD_FIX',3,2,YLD)
       CALL LCMPUT(IPDAT,'YLD_LOC',5,2,LOCYLD)
      ENDIF

      CALL LCMPTC(IPDAT,'ADF',3,1,ADF)
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
      CALL LCMPTC(IPDAT,'DET',12,1,HDET)
      ! SET THE IPDAT/STAVEC
      STAVEC(1) = DIMSAP(20)                     ! NGROUP
      STAVEC(3) = N_XS                           ! N_XS
      STAVEC(6) = NCRD                           ! NCOMPO
      STAVEC(7) = DIMSAP(31)                     ! NDLAY

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMPUT(IPDAT,'STATE-VECTOR',40,1,STAVEC)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
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
         WRITE(6,*) "********** CONTENT OF SAPHYB RECOVERED ***********"
         WRITE(6,*)
         WRITE(6,*) "NUMBER OF STATE VARIBALE IN PARAMDESCRIP : ", NPAR
         WRITE(6,*) "NUMBER OF STATE VARIABLES : ", NSVAR
         WRITE(6,*) "NAME OF STATE VARIABLES IN SAPHYB : ", PKEY
         WRITE(6,*) "STATE VARIABLES RECOGNIZED : ",PKEY(1:NSVAR)
         IF(NSVAR<NPAR-1) THEN
         WRITE(6,*) "WARNING:"
         WRITE(6,*) "STATE VARIABLES UNRECOGNIZED:",PKEY(NSVAR+1:NPAR-1)
         WRITE(6,*) "==>PLEASE USE THE PKEY CARD OF D2P: MODULE"
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
         ENDIF
         IF (STAVEC(21).EQ.1) THEN
         WRITE(6,*)'WARNING => ADF ARE INTEGRATED IN CROSS SECTIONS'
         CALL XABORT('STOP')
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
      DEALLOCATE (PKIDX)
      DEALLOCATE (SV_VAL)
      DEALLOCATE (VALPAR,PVALDIR)
      DEALLOCATE (NOMREA)
      DEALLOCATE (RANK_INDEX,PKEY_TMP)
      DEALLOCATE (RANK,NVAL,PKEY)
      RETURN
      END
