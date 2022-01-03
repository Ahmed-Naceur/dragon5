*DECK D2PINP
      SUBROUTINE  D2PINP( IPSAP, IPDAT , IPRINT, STAVEC, CRDINF,   NCRD,
     >                     PKEY,   ISOT,   MESH, USRPAR, USRVAL, USRSTA,
     >                  USRVAPK,    SAP,    MIC, EXC   ,   SCAT, ADF   ,
     >                      DEB,   FA_K,   LADD,   LNEW,    MIX,    XSC,
     >                   JOBOPT, SIGNAT, MIXDIR,    CDF,    GFF,   ADFD,
     >                     CDFD,    YLD, YLDOPT, LOCYLD,  OTHPK, OTHTYP,
     >                   OTHVAL,   HDET, OTHVAR,   THCK,   HFLX,   HCUR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* 1) Recover data from saphyb or multicompo object.
* 2) Build headers of GenPMAXS and Helios like file
*
*Copyright:
* Copyright (C) 2015 IRSN
*
*Author(s): 
* J. Taforeau
*
*Parameters: input/output
* IPSAP   address of saphyb or multicompo object
* IPDAT   address of data structure INFO
* NCRD    number of control rod composition recovered from D2P
*         input user
* MIX     index of mixture on which XS are to be extracted (only
*         for reflector cases)
* FA_K    assembly type
*         =0 reflector
*         =1 assembly
* USRSTA  state variable names recovered from GLOBAL record in D2P:
* USRVAL  number of value for state variables  recovered from GLOBAL
*         record in D2P:
* IPRINT  control the printing on screen
* STAVEC  various parameters associated with the IPDAT structure
* CRDINF  meaning of control rods in the IPSAP object
* XSC     XS_CONT recovered from D2P: input
* DEB     FLAG to indicate the first call to the D2PGEN subroutine
* USRVAPK value of state prameter set by the user and recoverd from
*         USER ADD option in D2P:
* ADF     type of ADF to be selected
* JOBOPT  flag for JOB_OPT record in IPINP object
* USRPAR  name of state variables (sapnam) in IPSAP associated to
*         DMOD TCOM etc. recovered from PKEY card in D2P:
* MESH    type of meshing to be applied for the branching calculation
* PKEY    name of state variable (refnam) recovered from PKEY card in
*         D2P:
* ISOT    name of isotopes in IPSAP for xenon samarium and promethium
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
* SIGNAT  signature of the object containing cross sections
* MIXDIR  directory that contains homogeneous mixture information
* CDF     type of CDF to be selected
* GFF     type of GFF to be selected
* ADFD    name of record for 'DRA' type of ADF
* CDFD    name of record for 'DRA' type of CDF
* YLD     user defined values for fission yields (1:I, 2:XE, 3:PM)
* LOCYLD  value for state parameter on which fission yield will be
*         calculated
* YLDOPT  option for fission yields calculation (DEF, MAN, FIX)
* HDET    name of isotope for the detector cross sections
* THCK    Thickness of reflector
* HFLX    Name of the record for the flux
* HCUR    Name of the record for the current
*
*Parameters: 
* OTHPK   
* OTHTYP  
* OTHVAL  
* OTHVAR  
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPDAT
      INTEGER NCRD,MIX,FA_K,USRSTA
      INTEGER IPRINT,DEB
      REAL THCK
      INTEGER STAVEC(40),CRDINF(20),USRVAL(12)
      REAL YLD(3),LOCYLD(5)
      REAL XSC(3)
      REAL USRVAPK(12,10),OTHVAR(12)
      CHARACTER JOBOPT(16)
      CHARACTER*3 ADF,CDF,GFF,YLDOPT
      CHARACTER*8 ADFD(4),CDFD(8)
      CHARACTER*5 MESH
      CHARACTER*8 PKEY(6),HFLX(2),HCUR(2)
      CHARACTER*12 ISOT(8), SIGNAT,MIXDIR,USRPAR(12)
      CHARACTER*12 OTHPK(12), OTHTYP(12), OTHVAL(12),HDET
      LOGICAL SAP, MIC, EXC,SCAT,LADD,LNEW
*----
*  LOCAL VARIABLES
*----
      LOGICAL :: LADF=.FALSE.
      LOGICAL :: LCDF=.FALSE.
      LOGICAL :: LGFF=.FALSE.
      LOGICAL :: LYLD=.FALSE.
      INTEGER NADF,NCDF

      IF (JOBOPT(1)=='T') THEN
        NADF=STAVEC(13)
        IF (NADF.NE.XSC(1)) THEN
         WRITE(6,*)'@D2PINP: INCOHERENT NUMBER OF ADF (',NADF,
     >   ')','AND NUMBER OF SIDES IN ASSEMBLY (',XSC(1),').'
         CALL XABORT ("=> CHECK CARD 'ADF' AND 'XS_CONT'")
        ENDIF
        IF ((SIGNAT.EQ.'L_SAPHYB').and.(ADF.EQ.'DRA')) THEN
          WRITE(6,*) "@D2PINP: ADF OF TYPE (",ADF,
     1     ") NOT YET IMPLEMENTED WITH SAPHYB OBJECT"
          WRITE(6,*)"=> WARNING :  ADF CALUCLATION IGNORED"
          LADF = .FALSE.
          JOBOPT(1)='F'
        ELSE IF ((SIGNAT.EQ.'L_MULTICOMPO').and.
     >     ((ADF.EQ.'SEL').OR.(ADF.EQ.'SEL'))) THEN
          WRITE(6,*) "@D2PINP: ADF OF TYPE (",ADF,
     1     " NOT YET IMPLEMENTED WITH MULTICOMPO OBJECT"
          WRITE(6,*)"=> WARNING :  ADF CALUCLATION IGNORED"
          LADF = .FALSE.
          JOBOPT(1)='F'
        ELSE
          LADF = .TRUE.
        ENDIF
      ELSE
        LADF = .FALSE.
      ENDIF
      IF (JOBOPT(10)=='T') THEN
        NCDF=STAVEC(15)
        IF (NCDF.NE.XSC(2)) THEN
         WRITE(6,*)'@D2PINP: INCOHERENT NUMBER OF CDF (',NCDF,
     >   ')','AND NUMBER OF CORNERS IN ASSEMBLY (',XSC(2),').'
         CALL XABORT ("=> CHECK CARD 'CDF' AND 'XS_CONT'")
        ENDIF
        IF (SIGNAT.EQ.'L_SAPHYB') THEN
          WRITE(6,*) "@D2PINP: CDF CALCULATION",
     1     " NOT YET IMPLEMENTED WITH SAPHYB OBJECT"
          WRITE(6,*)"=> WARNING :  CDF CALUCLATION IGNORED"
          LCDF = .FALSE.
          JOBOPT(10)='F'
        ENDIF
        IF (CDF.NE. 'DRA') THEN
         CALL XABORT ("@D2PINP UNKNOW CDF TYPE : "//CDF//'.')
        ENDIF
        LCDF = .TRUE.
      ELSE
        LCDF = .FALSE.
      ENDIF
      IF (JOBOPT(11)=='T') THEN
        IF (SIGNAT.EQ.'L_SAPHYB') THEN
          WRITE(6,*) "@D2PINP: GFF CALCULATION",
     1     " NOT YET IMPLEMENTED WITH SAPHYB OBJECT"
          WRITE(6,*)"=> WARNING :  GFF CALUCLATION IGNORED"
          LGFF = .FALSE.
          JOBOPT(11)='F'
        ENDIF
        IF (GFF.NE. 'DRA') THEN
         CALL XABORT ("@D2PINP UNKNOW GFF TYPE : '"//GFF//"'.")
        ENDIF
        LGFF = .TRUE.
      ELSE
        LGFF = .FALSE.
      ENDIF

      IF (JOBOPT(9)=='T') LYLD = .TRUE.
      IF ((JOBOPT(2)=='T').and.(JOBOPT(9)=='F')) THEN
          WRITE(6,*) "@D2PINP: JOB_OPT : XE/SM ARE REQUESTED (lxes=T) ",
     1     "BUT FISSION YIELDS ARE NOT RECOVERED (lyld=F) "
         WRITE(6,*) "=> THE lyld FLAG IS FORCED TO TRUE"
         JOBOPT(9)='T'
         LYLD = .TRUE.
      ENDIF

      IF((FA_K.EQ.1).OR.(FA_K.EQ.0)) THEN
*       CASE FOR FUEL PMAXS
        IF (SIGNAT.EQ.'L_SAPHYB') THEN
        STAVEC(18)=0
        WRITE(6,*) "*******    EXTRACTION OF DATA FROM SAPHYB      ****"
        CALL    D2PSAP (  IPSAP,  IPDAT, STAVEC, CRDINF,   NCRD,   PKEY,
     >                    ISOT ,   MESH, USRPAR, USRVAL, USRSTA,USRVAPK,
     >                    SAP  ,    MIC,    EXC,   SCAT,    ADF,   LADD,
     >                    LNEW ,   LADF, IPRINT,   LYLD,    YLD, YLDOPT,
     >                    LOCYLD,  HDET                                )

        ELSE  IF (SIGNAT.EQ.'L_MULTICOMPO') THEN
        STAVEC(18)=1
        WRITE(6,*) "*******    EXTRACTION OF DATA FROM MULTICOMPO  ****"
        WRITE(6,*)
        WRITE(6,*) "DIRECTORY:'",MIXDIR,"' AT MIXUTRE INDEX ",MIX,"."
        WRITE(6,*) "=> WARNING: CHECK EXISTENCE OF ",MIXDIR,"DIRECTORY."
        CALL LCMLIB(IPSAP)
         IF (LADF) THEN
          WRITE(6,*) "ADF CALCULATION REQUESTED:"
          WRITE(6,*)"=> WARNING: CHECK EXISTENCE OF ADF RECORDS"
         ENDIF
         IF (LCDF) THEN
          WRITE(6,*) "CDF CALCULATION REQUESTED:"
          WRITE(6,*)"=> WARNING: CHECK EXISTENCE OF '",CDFD(1:NCDF),
     >     "' RECORDS"
         ENDIF

        CALL  D2PMCO   (  IPSAP,  IPDAT, STAVEC, CRDINF,   NCRD,   PKEY,
     >                    ISOT ,   MESH, USRPAR, USRVAL, USRSTA,USRVAPK,
     >                    SAP  ,    MIC,    EXC,   SCAT,    ADF,   LADD,
     >                    LNEW ,   LADF, IPRINT, MIXDIR,    MIX,   LCDF,
     >                    LGFF ,    CDF,    GFF,   ADFD,   CDFD,  LYLD ,
     >                      YLD, YLDOPT, LOCYLD,  OTHPK, OTHTYP, OTHVAL,
     >                   OTHVAR,   THCK,   HFLX,   HCUR                )
        ELSE
         CALL XABORT ('@D2PINP: UNKNOWN SIGNATURE')
        ENDIF
      ELSE
        CALL XABORT('@D2PINP: PHASE 1: FUEL OR REFLECTOR CARD EXPECTED')
      ENDIF

      IF (YLDOPT.EQ.'MAN') THEN
       DEB = -1
      ELSE
       DEB = 0
      ENDIF

      IF (STAVEC(19).EQ.1) THEN
       WRITE(6,*)"=> WARNING: THE TEMPERATURE ARE INDIACTED IN K"
      ENDIF

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMPUT(IPDAT,'STATE-VECTOR',40,1,STAVEC)

      END
