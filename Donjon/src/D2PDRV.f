*DECK D2PDRV
      SUBROUTINE D2PDRV( NENTRY, HENTRY, IENTRY, JENTRY, KENTRY,    NGP,
     >                     NCRD,    MIX,   FA_K,   IUPS, USRSTA,  PHASE,
     >                   IPRINT, STAVEC, CRDINF, USRVAL,   VERS,   SFAC,
     >                     BFAC,   FC1,    FC2,    FC3,     FC4,    XSC,
     >                  USRVAPK,    ADF,    DER, JOBOPT, USRPAR,   MESH,
     >                     PKEY, FILNAM,   ISOT, JOBTIT,    COM,    SAP,
     >                      MIC,    EXC,   SCAT,   LADD,   LNEW, MIXDIR,
     >                      CDF,    GFF,   ADFD,   CDFD,    YLD, YLDOPT,
     >                   LOCYLD,   XESM,  ITEMP,  OTHPK, OTHTYP, OTHVAL,
     >                     HDET,  LPRC,  HEQUI,   HMASL,ISOTOPT,ISOTVAL,
     >                     LMEM,OTHVAR,   THCK,    HFLX,   HCUR        )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store an isotopic data recovered from a Saphyb into a Microlib.
*
*Copyright:
* Copyright (C) 2015 IRSN
*
*Author(s): 
* J. Taforeau
*
*Parameters: input/output
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
* NGP     number of energy groups recovered from D2P input user
* NCRD    number of control rod composition recovered from D2P input
*         user
* MIX     index of mixture on which XS are to be extracted (only for
*         reflector cases)
* FA_K    assembly type
*         =0 reflector
*         =1 assembly
* IUPS    treatment for upscattering
*         =0 keep up scatter XS
*         =1 remove up scatter XS, modify  down scatter with DRAGON
*         spectrum (not available in this version)
*         =2 remove up scatter XS, modify  down scatter with infinite
*         medium spectrum
* USRSTA  state variable names recovered from GLOBAL record in D2P:
* USRVAL  number of value for state variables  recovered from GLOBAL
*         record in D2P:
* PHASE   the current phase of D2P:
* IPRINT  control the printing on screen
* STAVEC  various parameters associated with the IPDAT structure
* CRDINF  meaning of control rods in the IPSAP object
* VERS    version of PARCS to be used
* SFAC    the scattering cross section factor
* BFAC    the multiplier for betas
* FC1     FILE_CONT_1 recovered from D2P: input
* FC2     FILE_CONT_2 recovered from D2P: input
* FC3     FILE_CONT_3 recovered from D2P: input
* FC4     FILE_CONT_4 recovered from D2P: input
* XSC     XS_CONT     recovered from D2P: input
* USRVAPK value of state prameter set by the user and recoverd from
*         USER ADD option in D2P:
* ADF     type of ADF to be selected
* DER     partials derivative (T) or row cross section (F) to be stored
*         in PMAXS
* JOBOPT  flag for JOB_OPT record in IPINP object
* USRPAR  name of state variables (sapnam) in IPSAP associated to
*         DMOD TCOM etc. recovered from PKEY card in D2P:
* MESH    type of meshing to be applied for the branching calculation
* PKEY    name of state variable (refnam) recovered from PKEY card in
*         D2P:
* FILNAM  name of IPINP
* ISOT    name of isotopes in IPSAP for xenon samarium and promethium
* JOBTIT  title of in header of PMAXS file
* COM     comment to be printed in PMAXS file
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
* MIXDIR  directory that contains homogeneous mixture information
* CDF     type of CDF to be selected
* GFF     type of GFF to be selected
* ADFD    name of record for 'DRA' type of ADF
* CDFD    name of record for 'DRA' type of CDF
* YLD     user defined values for fission yields (1:I, 2:XE, 3:PM)
* LOCYLD  value for state parameter on which fission yield will be calcu
* YLDOPT  option for fission yield calculation (DEF, MAN, FIX)
* XESM    option for comparing k-inf in GenPMAX (1: using Pm/Sm data;
*         2: using I/Xe data; 3: using I/Xe/Pm/Sm data)
* ITEMP   indicate if temperature is in C or in K in the SAP/MCO objec
* HDET    name of isotope for the detector cross sections
* LMER    ADF are merged in the cross sections
* THCK    Thickness of reflector
* HFLX    Name of the record for the flux
* HCUR    Name of the record for the current
*
*Parameters: 
* OTHPK   
* OTHTYP  
* OTHVAL  
* LPRC    
* HEQUI   
* HMASL   
* ISOTOPT 
* ISOTVAL 
* LMEM    
* OTHVAR  
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
      INTEGER NGP,NCRD,MIX,FA_K,IUPS,USRSTA,XESM
      INTEGER PHASE,IPRINT,ITEMP
      REAL THCK
      INTEGER STAVEC(40),CRDINF(20),USRVAL(12)
      REAL VERS,SFAC,BFAC,YLD(3),LOCYLD(5)
      REAL FC1(5),FC2(8),FC3(7),FC4(3),XSC(3)
      REAL USRVAPK(12,10),ISOTVAL,OTHVAR(12)
      CHARACTER JOBOPT(16)
      CHARACTER*12 OTHTYP(12),OTHPK(12),OTHVAL(12),HDET
      CHARACTER*3 ADF,CDF,GFF,YLDOPT
      CHARACTER*8 ADFD(4),CDFD(8)
      CHARACTER*4 DER,HEQUI,HMASL,JOB(4)
      CHARACTER*1 ISOTOPT
      CHARACTER*5 MESH
      CHARACTER*8 PKEY(6)
      CHARACTER*12 FILNAM,ISOT(6),SIGNAT,MIXDIR,USRPAR(12)
      CHARACTER*16 JOBTIT
      CHARACTER*8 HCUR(2),HFLX(2)
      CHARACTER*40 COM
      LOGICAL SAP, MIC, EXC,SCAT,LADD,LNEW,LPRC,LMEM
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPSAP,IPDAT,IPMIC
      INTEGER IPHEL,IPINP,IPPRC
      INTEGER DEB,REW
      CHARACTER TEXT12*12,HSIGN*12,HSMG*131

      IF (IPRINT.EQ.-1) IPRINT = 0
*----
*  PHASE 1 : SET HEADER OF GENPMAXS INPUT FILE (.inp) AND HELIOS LIKE FI
*----
      IF (PHASE.EQ.1) THEN
       IF(NENTRY.NE.3) CALL XABORT('@D2PDRV: 3 PARAMETERS EXPECTED.')
       IF((IENTRY(1).NE.4)) THEN
        WRITE(HSMG,'(12H@D2P: ENTRY ,A12,24H IS NOT OF SEQUENTIAL AS,
     >  9HNCII TYPE)') HENTRY(2)
        CALL XABORT(HSMG)
       ELSE IF(JENTRY(1).EQ.2) THEN
        WRITE(HSMG,'(12H@D2P: ENTRY ,A12,24H IS NOT IN CREATION OR I,
     >  19HN MODIFICATION MODE)')
     >  HENTRY(1)
        CALL XABORT(HSMG)
       ENDIF
       IF(IENTRY(2).NE.1) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,19H IS NOT OF LCM TYPE)')
     >  HENTRY(2)
        CALL XABORT(HSMG)
       ELSE IF(JENTRY(2).EQ.2) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT IN CREATION OR I,
     >  19HN MODIFICATION MODE)')
     >  HENTRY(2)
        CALL XABORT(HSMG)
       ENDIF
       CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
       IF(HSIGN.NE.'L_INFO') THEN
        TEXT12=HENTRY(2)
        CALL XABORT('@D2P: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_INFO EXPECTED.')
       ENDIF
      IF(IENTRY(3).GT.2) THEN
       WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,19H IS NOT OF LCM TYPE)')
     >  HENTRY(3)
       CALL XABORT(HSMG)
      ELSE IF(JENTRY(3).NE.2) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,20H IS NOT IN READ ONLY,
     >  5H MODE)')
     >  HENTRY(3)
        CALL XABORT(HSMG)
       ENDIF
       CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
       SIGNAT=HSIGN
       IF((HSIGN.NE.'L_SAPHYB')) THEN
        IF(HSIGN.NE.'L_MULTICOMPO') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('@D2PDRV: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_SAPHYB OR L_MULTICOMPO EXPECTED.')
        ENDIF
       ENDIF

       IPPRC=FILUNIT(KENTRY(1))
       IPDAT=KENTRY(2) ! output INFO address
       IPSAP=KENTRY(3) ! input saphyb address
*
       STAVEC(8)=INT(FC1(1))
       STAVEC(9)=INT(FC1(2))
       STAVEC(10)=INT(FC1(3))
       STAVEC(11)=INT(XSC(1))
       STAVEC(12)=INT(XSC(2))
       STAVEC(19)=ITEMP
       IF ((XSC(1).GT.4).OR.(XSC(2).GT.8)) THEN
        CALL XABORT ('@D2PDRV: CARD XS_CONT : NSIDE AND NCORNER'
     1   //' CANNOT EXCEED 4 AND 8 RESPECTIVELY.')
       ENDIF
       IF (MESH.EQ.'GLOB'.OR.MESH.EQ.'ADD') THEN
        IF ((JOBOPT(1).EQ.'T').AND. (ADF .NE. 'DRA')) THEN
         CALL XABORT('@D2PDRV: ADF OF TYPE (SEL/GET) CANNOT BE EXTRACT'
     1   //'ED WITH USER DEFINED BRANCHING CALCULATION')
        ENDIF
       ENDIF
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*        DRAG2PARCS INPUT DATA RECOVERED           *"
       WRITE(6,*) "****************************************************"
      IF(IPRINT > 0)  THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "* PHASE 1 : RECOVER DATA AND CREATE INPUT FILES    *"
       WRITE(6,*) "****************************************************"
       WRITE(6,*)
      ENDIF
*----
      CALL        D2PINP( IPSAP, IPDAT , IPRINT, STAVEC, CRDINF,   NCRD,
     >                     PKEY,   ISOT,   MESH, USRPAR, USRVAL, USRSTA,
     >                  USRVAPK,    SAP,    MIC, EXC   ,   SCAT, ADF   ,
     >                      DEB,   FA_K,   LADD,   LNEW,    MIX,    XSC,
     >                   JOBOPT, SIGNAT, MIXDIR,    CDF,    GFF,   ADFD,
     >                     CDFD,    YLD, YLDOPT, LOCYLD,   OTHPK,OTHTYP,
     >                   OTHVAL,   HDET, OTHVAR,   THCK,    HFLX,  HCUR)

      CALL      D2PGEN ( IPINP, IPDAT,   STAVEC, JOBTIT, FILNAM,    DER,
     >                    VERS,   COM,   JOBOPT,   IUPS,   FA_K,   SFAC,
     >                    BFAC,   DEB,     XESM,  FC1  ,    FC2,    FC3,
     >                     FC4,   XSC,   IPRINT                        )

      IF (LPRC) THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*                BUILDING PROCEDURE                *"
       WRITE(6,*) "****************************************************"
       CALL  D2PPRC(   IPDAT, IPPRC,HEQUI, HMASL, ISOTVAL, ISOTOPT,LMEM,
     >                IPRINT,MIXDIR,JOBOPT                             )
      ENDIF
      IF(IPRINT > 0)  THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*                END OF PHASE 1                    *"
       WRITE(6,*) "****************************************************"
      ENDIF
*----
*  PHASE 2 : BRANCHING CALCULATION
*----
      ELSE IF (PHASE.EQ.2) THEN
       IF(NENTRY.NE.5) CALL XABORT('@D2PDRV: 5 PARAMETERS EXPECTED.')
       IF(IENTRY(1).NE.4) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT OF SEQUENTIAL AS,
     >  9HNCII TYPE)') HENTRY(1)
        CALL XABORT(HSMG)
       ELSE IF(JENTRY(1).EQ.2) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT IN CREATION OR I,
     >  19HN MODIFICATION MODE)')
     >  HENTRY(1)
        CALL XABORT(HSMG)
       ENDIF
       IF((IENTRY(5).GT.2)) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,19H IS NOT OF XSM TYPE)')
     >  HENTRY(5)
        CALL XABORT(HSMG)
       ELSE IF(JENTRY(5).NE.2) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT IN READ-ONLY MOD,
     >  1HE)')
     >  HENTRY(5)
        CALL XABORT(HSMG)
       ENDIF
       CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
       IF(HSIGN.NE.'L_INFO') THEN
        TEXT12=HENTRY(3)
        CALL XABORT('@D2PDRV: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1  '. L_INFO EXPECTED.')
       ENDIF
       CALL LCMGTC(KENTRY(5),'SIGNATURE',12,1,SIGNAT)
       IF((SIGNAT.NE.'L_SAPHYB')) THEN
        IF(SIGNAT.NE.'L_MULTICOMPO') THEN
         TEXT12=HENTRY(5)
         CALL XABORT('@D2PDRV: SIGNATURE OF '//TEXT12//' IS '//SIGNAT//
     1   '. L_SAPHYB OR L_MULTICOMPO EXPECTED.')
        ENDIF
       ENDIF
      IPHEL=FILUNIT(KENTRY(1))
      IPINP=FILUNIT(KENTRY(2)) ! input GENPMAXS file unit
      IPDAT=KENTRY(3) ! input DATA vector address
      IPMIC=KENTRY(4) ! input Microlib vector address
      IPSAP=KENTRY(5) ! input SAPHYB OBJECT
      IF(IPRINT > 0)  THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*     PHASE 2 : RECOVER CROSS SECTIONS OF BRANCH   *"
       WRITE(6,*) "****************************************************"
      ENDIF
      CALL LCMSIX(IPDAT,' ',0)

      CALL LCMGET(IPDAT,'STATE-VECTOR',STAVEC)
      IF (STAVEC(18).EQ.1) THEN
       CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
       CALL LCMGTC(IPDAT,'NAMDIR',12,1,MIXDIR)
       CALL LCMSIX(IPDAT,' ',0)
      ENDIF
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMGTC(IPDAT,'JOB_OPT',4,4,JOB)
      NGP = STAVEC(1)
      i=1
      DO j=1,4
       DO k=1,4
        JOBOPT(i)= JOB(j)(k:k)
        i=i+1
       ENDDO
      ENDDO
      CALL LCMGET(IPDAT,'FLAG',DEB)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMGET(IPDAT,'REWIND',REW)
      IF ((DEB.LE.0).AND.(REW.EQ.1))THEN

       WRITE(6,*)"*******   CREATION OF HELIOS ANS GENPMAXS FILES *****"
       CALL      D2PHEL( IPHEL,  IPDAT,  IPMIC ,  IPINP,       STAVEC,
     >                  JOBOPT, IPRINT                               )

      ENDIF

      CALL D2PXS(IPDAT,IPMIC,IPSAP,STAVEC,SIGNAT,MIXDIR,JOBOPT,IPRINT)

      CALL LCMSIX(IPDAT,' ',0)

      CALL LCMPUT(IPDAT,'STATE-VECTOR',40,1,STAVEC)

      IF(IPRINT > 0)  THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*                END OF PHASE 2                    *"
       WRITE(6,*) "****************************************************"
      ENDIF
*----
*  PHASE 3 : STORE BRANCHES IN HELIOS FILE
*----
      ELSE IF (PHASE.EQ.3) THEN
       IF(NENTRY.GT.4) CALL XABORT('@D2PDRV: 3 PARAMETERS EXPECTED.')
       IF((IENTRY(1).NE.4)) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT OF SEQUENTIAL AS,
     >  9HNCII TYPE)') HENTRY(1)
        CALL XABORT(HSMG)
       ELSE IF(JENTRY(1).EQ.2) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT IN CREATION OR I,
     >  19HN MODIFICATION MODE)')
     >  HENTRY(1)
        CALL XABORT(HSMG)
       ENDIF
       IF(IENTRY(2).NE.4) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT OF SEQUENTIAL AS,
     >  9HNCII TYPE)') HENTRY(2)
        CALL XABORT(HSMG)
       ELSE IF(JENTRY(2).EQ.2) THEN
        WRITE(HSMG,'(15H@D2PDRV: ENTRY ,A12,24H IS NOT IN CREATION OR I,
     >  19HN MODIFICATION MODE)')
     >  HENTRY(2)
        CALL XABORT(HSMG)
       ENDIF
       CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
       IF(HSIGN.NE.'L_INFO') THEN
        TEXT12=HENTRY(2)
        CALL XABORT('@D2PDRV: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     >  '. L_INFO EXPECTED.')
       ENDIF
       IPHEL=FILUNIT(KENTRY(1)) !  dragon file unit
       IPINP=FILUNIT(KENTRY(2)) !  GENPMAXS file unit
       IPDAT=KENTRY(3) !  DATA vector address
       CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
       CALL LCMGET(IPDAT,'FLAG',DEB)
      IF (DEB<0) THEN
       DEB=1
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*            END OF FISSION YIELD BRANCH           *"
       WRITE(6,*) "****************************************************"
       CALL LCMPUT(IPDAT,'FLAG',1,1,DEB)
      ELSE
      IF(IPRINT > 0)  THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*     PHASE 3 : STORE BRANCHES IN HELIOS FILE      *"
       WRITE(6,*) "****************************************************"
      ENDIF
       WRITE(6,*) "***** STORE CURRENT BRANCH IN HELIOS LIKE FILE *****"
       DEB=1

       CALL LCMSIX(IPDAT,' ',0)
       CALL LCMGET(IPDAT,"STATE-VECTOR",STAVEC)
       CALL D2PBRA( IPDAT,IPINP,IPHEL,STAVEC,DEB,SIGNAT,IPRINT)

      IF(IPRINT > 0)  THEN
       WRITE(6,*) "****************************************************"
       WRITE(6,*) "*                END OF PHASE 3                    *"
       WRITE(6,*) "****************************************************"
      ENDIF
      ENDIF
      ENDIF
      END
