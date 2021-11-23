*DECK D2PBRA
      SUBROUTINE D2PBRA( IPDAT,IPINP,IPHEL,STAVEC,DEB,SIGNAT,IPRINT    )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover information from the INFO data block for a complete branch
* and write it in the IPHEL file . The format of this file is described
* in the DRAG2PARCS: manual. This routine write sequentially the IPHEL
* file, branch after branch
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of info data block
* IPINP   file unit of the input file GENPMAXS.inp
* IPHEL   file unit of the HELIOS.dra file
* STAVEC  various parameters associated with the IPDAT structure
* DEB     flag for D2PGEN
* SIGNAT  signature of the object containing cross sections
* IPRINT  control the printing on screen
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT
      INTEGER IPINP,IPHEL,STAVEC(40),DEB,IPRINT

      CHARACTER*16 SIGNAT
*----
*  LOCAL VARIABLES
*----
      INTEGER GRID,ITBRAN,NVAR,i,j,k
      INTEGER NSF,NGP,NXS,NBU,FLPRIN,LMER
      INTEGER STAIDX(STAVEC(2)),PKIDX(STAVEC(2))
      INTEGER IUPS,FA_K,NADF,NCDF,NPIN,NCOLA,NROWA,XESM
      REAL XS(STAVEC(1),STAVEC(3),STAVEC(4)) ! TABLE FOR XS
      REAL ADF(STAVEC(13),STAVEC(1),STAVEC(4))
      REAL FLXL(STAVEC(1),STAVEC(4))
      REAL FLXR(STAVEC(1),STAVEC(4))
      REAL CURL(STAVEC(1),STAVEC(4))
      REAL CURR(STAVEC(1),STAVEC(4))
      REAL CDF(STAVEC(15),STAVEC(1),STAVEC(4))
      REAL GFF(STAVEC(8),STAVEC(9),STAVEC(1),STAVEC(4))
      REAL SCAT(STAVEC(1)*STAVEC(1),STAVEC(4))
      REAL BURN(STAVEC(4)),XSC(3),DATSRC(5)
      REAL DIV(3,STAVEC(4))
      REAL ND(2,STAVEC(4))
      CHARACTER(len=4) BRANCH,JOB(4)
      CHARACTER*12 FILNAM
      CHARACTER COM
      CHARACTER*16 JOBTIT
      CHARACTER JOBOPT(16)
      CHARACTER*3 ADF_T
      CHARACTER*1 DER
      REAL FC1(5)
      REAL FC2(8)
      REAL FC3(7)
      REAL FC4(3)
      REAL VERS,SFAC,BFAC
      LOGICAL :: LTH = .FALSE.
      LOGICAL :: LADF = .FALSE.
      LOGICAL :: LXES = .FALSE.
      LOGICAL :: LCDF = .FALSE.
      LOGICAL :: LGFF = .FALSE.
      LOGICAL :: LDET = .FALSE.


      ! INITIALIZATION OF VARIABLES
      NGP=STAVEC(1)
      NVAR=STAVEC(2)
      NXS=STAVEC(3)
      NBU=STAVEC(4)
      GRID=STAVEC(5)
      NCOLA=STAVEC(8)
      NROWA=STAVEC(9)
      NPART=STAVEC(10)
      NSF=STAVEC(11)
      NCF=STAVEC(12)
      NADF=STAVEC(13)
      NCDF=STAVEC(15)
      NGFF=STAVEC(16)
      NPIN=STAVEC(17)
      LMER=STAVEC(21)


      IF(IPRINT > 0)  THEN
        WRITE(6,*)
        WRITE(6,*) "**** WRITING CURRENT BRANCH IN HELIOS FILE ****"

      ENDIF
      ! RECOVER INFORMATION FROM INFO DATA BLOCK
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'HELIOS_HEAD',1)
      CALL LCMGET(IPDAT,'FILE_CONT_1',FC1)
      CALL LCMGET(IPDAT,'FILE_CONT_2',FC2)
      CALL LCMGET(IPDAT,'FILE_CONT_3',FC3)
      CALL LCMGET(IPDAT,'FILE_CONT_4',FC4)
      CALL LCMGET(IPDAT,'XS_CONT',XSC)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMGET(IPDAT,'DAT_SRC',DATSRC)
      CALL LCMGTC(IPDAT,'JOB_OPT',4,4,JOB)
      CALL LCMGET(IPDAT,'VERSION',VERS)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMGET(IPDAT,'PRINT',FLPRIN)
      CALL LCMGET(IPDAT,'STATE_INDEX',STAIDX)
      CALL LCMGET(IPDAT,'BRANCH_IT',ITBRAN)
      CALL LCMGTC(IPDAT,'BRANCH',4,1,BRANCH)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMGET(IPDAT,'PKIDX',PKIDX)
      CALL LCMGET(IPDAT,'BURN',BURN)

      i=1
      DO j=1,4
       DO k=1,4
        JOBOPT(i)= JOB(j)(k:k)
        i=i+1
       ENDDO
      ENDDO

      IF(JOBOPT(1)=='T') THEN
      LADF = .TRUE.
      CALL LCMGTC(IPDAT,'ADF_TYPE',3,1,ADF_T)
      ENDIF
      IF(JOBOPT(2)=='T') LXES = .TRUE.
      IF(JOBOPT(8)=='T') LDET = .TRUE.
      IF((JOBOPT(5)=='T').OR.(JOBOPT(7)=='T').OR.
     > (JOBOPT(9)=='T').OR.(JOBOPT(13)=='T').OR.(JOBOPT(12)=='T'))THEN
        LTH =.TRUE.
      ENDIF
      IF(JOBOPT(10)=='T') LCDF = .TRUE.
      IF(JOBOPT(11)=='T') LGFF = .TRUE.

      ! WRITE THE CURRENT BRANCH IN THE HELIOS.DRA FILE
      IF(FLPRIN==1) THEN
         ! RECOVER CROSS SECTIONS FROM THE TEMPORARY FILE
         CALL READXS (IPDAT,     XS,   SCAT,    ND,     DIV,    NGP,
     >                  NXS,    ADF,    CDF,   GFF,     NBU,   NADF,
     >               DATSRC,   GRID,   NCDF, NCOLA,   NROWA,   LADF,
     >                 LCDF,   LGFF,   LXES,  LDET,  SIGNAT,   LMER,
     >               IPRINT,  ADF_T,   FLXL,  FLXR,    CURL,   CURR)
         ! WRITE IN HELIOS.DRA THE SET OF BURNUP POINTS
         CALL SETBU (IPHEL,BRANCH,ITBRAN,XSC,BURN,NBU,          IPRINT)

         ! WRITE IN HELIOS.DRA THE SET OF CROSS SECTIONS
         CALL SETXS   (  IPHEL, BRANCH, ITBRAN,    XS,     NGP,    NXS,
     >                     NBU,   BURN, DATSRC,    LXES,   LDET,IPRINT)

         ! WRITE IN HELIOS.DRA THE ELEMENT OF THE SCATTERING MATRIX
         CALL SETSCT(IPHEL,BRANCH,ITBRAN,SCAT,NGP,NBU,BURN,     IPRINT)

         IF(LADF.AND.(LMER.EQ.0)) THEN
          CALL SETADF(  IPHEL, BRANCH, ITBRAN,   ADF,    NADF,    NGP,
     >                    NBU,   BURN, IPRINT, ADF_T,    FLXR,   FLXL,
     >                   CURL,   CURR)
         ENDIF

         IF(DATSRC(3)==1.0) THEN
           IF(LXES) THEN
            ! WRITE IN HELIOS.DRA THE NUMBRE DENSITIES FOR XENON AND
            ! SAMARIUM
            CALL SETND (IPHEL,BRANCH,ITBRAN,  ND,NBU,BURN,       IPRINT)
           ENDIF
           IF((GRID<2).AND.(SIGNAT.EQ.'L_SAPHYB'))THEN
            ! WRITE IN HELIOS.DRA THE DIVERS INFORMATION
            CALL SETDIV(IPHEL,BRANCH,ITBRAN,DIV,NBU,BURN,IPRINT)
           ENDIF
           IF(LTH) THEN
           ! WRITE IN HELIOS.DRA THE T:H INVARIANT DATA BLOCK
           CALL SETTH  (  IPHEL, BRANCH, ITBRAN,  BURN,     NBU, JOBOPT,
     >                      NGP,  IPDAT, IPRINT                        )
           ENDIF

           IF(LCDF) THEN
            CALL SETCDF(  IPHEL, BRANCH, ITBRAN,   CDF,    NCDF,    NGP,
     >                      NBU,   BURN, IPRINT                        )
           ENDIF
           IF(LGFF) THEN
            IF ((NCOLA .NE. NPIN) .OR. (NROWA .NE.NPIN)) THEN
             WRITE (6,*) "@D2PBRA: NUMBER OF PIN IN MCO (NPIN= ",NPIN,
     >     ") INCOHERENT WITH ncols AND nrows (",NCOLA,') IN D2P: INPUT'
            CALL XABORT ('')
            ENDIF
            CALL SETGFF(  IPHEL, BRANCH, ITBRAN,   GFF,   NCOLA,  NROWA,
     >                    NPART,    NGP,    NBU,  BURN,    NGFF, IPRINT,
     >                     VERS)
           ENDIF
         ENDIF
         ! SIGNATURE OF THE END OF A BRANCH (MANDATORY FOR GENPMAXS
         ! CODE)
         WRITE(IPHEL,*)
         WRITE(IPHEL,30)'*********************************************'
         WRITE(IPHEL,30)'* Normal End,    No warning messages issued *'
         WRITE(IPHEL,30)'*                                           *'
         WRITE(IPHEL,30)'*   Total CPU time used =                   *'
         WRITE(IPHEL,30)'*********************************************'
  30     FORMAT(25X,A)
      ENDIF

      ! UPDATE OF THE INFO DATA BLOCK
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMPUT(IPDAT,'FLAG',1,1,1)

      IF(IPRINT > 0)  THEN
        WRITE(6,*) "******** UPDATING the GENPMAXS.INP FILE *********"
      ENDIF
      ! UPDATE OF THE GENPMAXS.INP FILE (MANY ARGUMENTS IN THIS CALL
      ! ARE NOT USED IN D2PGEN)
      CALL      D2PGEN( IPINP, IPDAT,   STAVEC, JOBTIT, FILNAM,    DER,
     >                   VERS,   COM,   JOBOPT,   IUPS,   FA_K,   SFAC,
     >                   BFAC,   DEB,     XESM,  FC1  ,    FC2,    FC3,
     >                    FC4,   XSC,   IPRINT                        )
      IF(IPRINT > 0)  THEN
        WRITE(6,*)"********* SELECTING A NEW BRANCH CALCULATION *****"
      ENDIF

      CALL D2PSEL   (  IPDAT,  IPINP, STAVEC,BRANCH,  ITBRAN, STAIDX,
     >                  NVAR, JOBOPT,    DEB,  FC1  ,    FC2,    FC3,
     >                   FC4,   XSC,   IPRINT                       )


      WRITE(6,*) "*********       BRANCH SELECTED              *****"

      END

      SUBROUTINE SETBU(IPHEL,BRANCH,ITBRAN,XSC,BURN,NBU,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the inforamtion about burnup points and
* XSC card (sides in assembly (NSIDES),
* corners in assembly (NCORNERS), VFCM).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* XSC        content of the XS_CONT card
* BURN       set of burnup points
* NBU        number of bunup points
* IPRINT     control the printing on screen
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBU,ITBRAN,IPHEL,IPRINT
      REAL XSC(3),BURN (NBU)
      CHARACTER BRANCH*4
*----
*  LOCAL VARIABLES
*----
      ! number of sides and corners in assembly
      INTEGER NSIDE, NCORNER



      NSIDE = NINT(XSC(1))
      NCORNER = NINT(XSC(2))

      ! XS_CONT CARD (Cf DRAG2PARCS Manual for details on HELIOS format)
      IF (IPRINT>5) WRITE(6,*) 'SETBU: WRITE BURNUP INFO'
      ! HEADER OF XS_CONT card
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s)  1) ==========================='
      WRITE(IPHEL,*) '               2) %STAT_xxxx'
      WRITE(IPHEL,*) '               3) ==========================='
      WRITE(IPHEL,*) '               4)%XS_CONT'
      WRITE(IPHEL,*) '               5)Meaning : NBN,NSIDE,NCORNER,'
     1 //'VFCM'


      ! RIEGO block of HELIOS.dra file
      CALL SET_RIEGO(IPHEL)


      ! Set the content of XS_CONT in HELIOS.dra
      WRITE(IPHEL,'(25X,4A14)') '           NBN',
     1   '         NSIDE','       NCORNER','          VFCM'
      WRITE(IPHEL,200) 'Label E','.-.-E-.-.','.-.-E-.-.','.-.-E-.-.',
     1     '1-.-E-.-.'
      WRITE(IPHEL,'(I4,1X,A,I4,A,A,I4,A,I5,3I12,ES12.5E2)')
     1 1,BRANCH(1:2),ITBRAN,' ',BRANCH(1:2),ITBRAN,':',
     2 0,NBU,NSIDE,NCORNER,
     3 XSC(3)

      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'

      ! BURNUP INFORMATION


      ! HEADER OF Burnup card
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s)  1) ==========================='
      WRITE(IPHEL,*) '               2) %XS_STAT'
      WRITE(IPHEL,*) '               3) ==========================='
      WRITE(IPHEL,*) '               4)Meaning : Bunrup'


      ! RIEGO block of HELIOS.dra file
      CALL SET_RIEGO(IPHEL)


      WRITE(IPHEL,'(30X,A6)') 'BURNUP'
      WRITE(IPHEL,210) 'Label E','.-.-E-.-.'
      ! LOOP over burnup points
      DO IT=1, NBU

          WRITE(IPHEL,220) IT,BRANCH(1:2),ITBRAN,' ',
     1    BRANCH(1:2),ITBRAN,':',NINT(BURN(IT)),BURN(IT)/1000.0

      ENDDO

      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'


      ! format of HELIOS.dra file
  200 FORMAT(6X,A,12X,A,3X,A,3X,A,3X,A)
  210 FORMAT(6X,A,12X,A)
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5,6X,F6.3)
      END

      SUBROUTINE READXS (IPDAT,     XS,   SCAT,    ND,     DIV,    NGP,
     >                     NXS,    ADF,    CDF,   GFF,     NBU,   NADF,
     >                    DATSRC,  GRID,  NCDF, NCOLA,   NROWA,   LADF,
     >                      LCDF,  LGFF,  LXES,  LDET,  SIGNAT,   LMER,
     >                    IPRINT, ADF_T,  FLXL,  FLXR,    CURL,   CURR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover cross section from the INFO data block.
*
*parameters: input
* IPDAT      address of info data block
* XS         table of cross sections
* SCAT       scattering matrix
* ND         number densities for xenon and samarium
* DIV        divers info directory
* NGP        number of energy groups
* NXS        number of cross sections
* NBU        number of burnup points
* ADF        assembly dicontinuity factor
* NADF       number of surfaces in assembly
* NCDF       number of corners in assembly
* NCOLA      number of pin in assembly along x-axis
* NROWA      number of pin in assembly along y-axis
* GRID       type of gridding for branching calculation
* LADF       flag for assembly discontinuity factors
* LCDF       flag for corner discontinuity factors
* LGFF       flag for group form factors
* LXES       flag for microscopic cross sections
* DAT SRC    array containing the DATA source (reflector of fuel)
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT
      INTEGER NGP,NBU,NXS,NADF,GRID,NCDF,LMER
      REAL DATSRC(5)
      REAL XS(NGP,NXS,NBU)
      REAL SCAT(NGP*NGP,NBU)
      REAL ND(2,NBU)
      REAL DIV(3,NBU)
      REAL ADF(NADF,NGP,NBU)
      REAL FLXL(NGP,NBU)
      REAL FLXR(NGP,NBU)
      REAL CURL(NGP,NBU)
      REAL CURR(NGP,NBU)
      REAL CDF(NCDF,NGP,NBU)
      REAL GFF(NCOLA,NROWA,NGP,NBU)
      REAL ADFMOY(NGP,NBU)
      LOGICAL LADF,LXES,LCDF,LGFF,LDET
      CHARACTER*16 SIGNAT
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPTH,KPTH
      INTEGER BU
      CHARACTER*3 ADF_T

      IF(IPRINT>5) WRITE(6,*) 'READXS: RECOVER INFO DATA BLOCK'

      ! LOOP over burnup points
      DO BU=1, NBU
         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
         IPTH=LCMGID(IPDAT,'CROSS_SECT')
         KPTH=LCMDIL(IPTH,BU)
         CALL LCMSIX(KPTH,'MACROLIB_XS',1)
         CALL LCMGET(KPTH,'XTR',XS(1:NGP,1,BU))
         CALL LCMGET(KPTH,'ABSORPTION',XS(1:NGP,2,BU))
         CALL LCMGET(KPTH,'X_NU_FI',XS(1:NGP,3,BU))
         CALL LCMGET(KPTH,'KAPPA_FI',XS(1:NGP,4,BU))
         IF(LXES)CALL LCMGET(KPTH,'SFI',XS(1:NGP,7,BU))
         IF(LADF) THEN
         IF (ADF_T.EQ.'DRA')THEN
          CALL LCMGET(KPTH,'ADF',ADF(:,:,BU))
         ELSE IF(ADF_T.EQ.'GEN')THEN
          CALL LCMGET(KPTH,'FLXL',FLXL(:,BU))
          CALL LCMGET(KPTH,'FLXR',FLXR(:,BU))
          CALL LCMGET(KPTH,'CURR',CURR(:,BU))
          CALL LCMGET(KPTH,'CURL',CURL(:,BU))
         ENDIF
         ENDIF
         IF(LCDF)CALL LCMGET(KPTH,'CDF',CDF(:,:,BU))
         IF(LGFF)CALL LCMGET(KPTH,'GFF',GFF(:,:,:,BU))


         CALL LCMGET(KPTH,'SCAT',SCAT(1:NGP*NGP,BU))
         IF(DATSRC(3)==1) THEN

          IF((LXES).OR.(LDET)) THEN
            CALL LCMSIX(KPTH,' ',2)
            CALL LCMSIX(KPTH,'MICROLIB_XS',1)

            IF(LDET) CALL LCMGET(KPTH,'DET',XS(1:NGP,8,BU))
            IF (LXES) THEN
             CALL LCMGET(KPTH,'XENG',XS(1:NGP,5,BU))
             CALL LCMGET(KPTH,'SMNG',XS(1:NGP,6,BU))
             CALL LCMGET(KPTH,'XEND',ND(1,BU))
             CALL LCMGET(KPTH,'SMND',ND(2,BU))
            ENDIF
          ENDIF
          IF((GRID<2).and. (SIGNAT.EQ.'L_SAPHYB')) THEN
            CALL LCMSIX(IPDAT,' ',0)
            CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
            IPTH=LCMGID(IPDAT,'DIVERS')
            KPTH=LCMDIL(IPTH,BU)
            CALL LCMGET(KPTH,'KEFF',DIV(1,BU))
            CALL LCMGET(KPTH,'KINF',DIV(2,BU))
            CALL LCMGET(KPTH,'B2',DIV(3,BU))
          ENDIF
         ENDIF
      ENDDO
      IF (LMER.EQ.1) THEN
       DO I=1,NGP
        DO BU=1,NBU
         ADFMOY(I,BU)=SUM(ADF(1:NADF,I,BU))/NADF
        ENDDO
       ENDDO

       DO I=1,NGP
        DO BU=1,NBU
         SCAT(I,BU)=SCAT(I,BU)/ADFMOY(NGP-1+1,BU)
         SCAT(I+NGP,BU)=SCAT(I+NGP,BU)/ADFMOY(NGP-I+1,BU)
         XS(I,1,BU)=XS(I,1,BU)*ADFMOY(I,BU)
         XS(I,2:NXS,BU)=XS(I,2:NXS,BU)/ADFMOY(I,BU)
        ENDDO
       ENDDO
      ENDIF
      CALL LCMSIX(IPDAT,' ',0)
      END

      SUBROUTINE SETXS(  IPHEL, BRANCH, ITBRAN,    XS,     NGP,    NXS,
     >                     NBU,   BURN, DATSRC,    LXES,  LDET, IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the cross sections for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* XS         table of cross sections
* NGP        number of energy groups
* NXS        number of cross sections
* NBU        number of burnup points
* BURN       set of burnup points
* IPRINT     control the printing on screen
* DATSRC     array containing the DATA source (reflector of fuel)
* LXES       flag for presence of micoscopic cross sections
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,NXS,NGP,ITBRAN
!      REAL XS(NGP,NXS,NBU),BURN (NBU),DATSRC(3)
      REAL XS(NGP,NXS,NBU),BURN (NBU),DATSRC(5)
      CHARACTER(len=4) BRANCH,XS_name
      LOGICAL LXES,LDET
*----
*  LOCAL VARIABLES
*----
      INTEGER XST  ! INDEX OF CROSS SECTIONS
      REAL FA_KIND
      LOGICAL :: LXS = .TRUE.

      IF(IPRINT>5) WRITE(6,*) 'SETXS: WRITE INFO FOR A BANCH'

      FA_KIND=DATSRC(3)

      ! LOOP OVER CROSS SECTIONS TYPE
      DO XST=1, NXS
         LXS = .TRUE.
         SELECT CASE (XST)
         CASE (1)
            XS_name = 'STR'               ! TRANSPORT XS
         CASE (2)
            XS_name = 'SAB'               ! ABSORPTION XS
         CASE (3)
            XS_name = 'SNF'               ! NU SIGMA FISSION XS
         CASE (4)
            XS_name = 'SKF'               ! KAPPA FISSION XS
         CASE (5)
            IF(.NOT. LXES) LXS=.FALSE.
            XS_name = 'XENG'              ! XE MICROSCOPIC ABSORPTION XS
         CASE (6)
            IF(.NOT. LXES) LXS=.FALSE.
            XS_name = 'SMNG'              ! SM MICROSCOPIC ABSORPTION XS
         CASE (7)
            IF(.NOT. LXES) LXS=.FALSE.
            XS_name = 'SFI'               ! FISSION XS
         CASE (8)
            IF(.NOT. LDET) LXS=.FALSE.
            XS_name = 'DET'               ! DETECTOR XS
         END SELECT
         IF(LXS) THEN
          ! LABEL FOR XS TYPE
          WRITE(IPHEL,'()')
          WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
          WRITE(IPHEL,*) 'Labels Array    : PrintVector'
          WRITE(IPHEL,110) XS_name
          WRITE(IPHEL,120) XS_name

  110     FORMAT(29H List Title(s)  1) %XS_PRIN %,A)
  120     FORMAT(34H Meaning : (.-.-E-G-.) G-th Group ,A,
     1    15H cross sections)

          CALL SET_RIEGO(IPHEL)

         ! LOOP OVER ENERGY GROUPS
         ! CREATION OF LABEL FOR CROSS SECTIONS
          DO IT=1, NGP
           IF(IT==1) THEN
            WRITE(IPHEL,'(27X,A4,A2)',advance='no') XS_name,'Xs'
           ELSE IF(IT==NGP .OR. IT==8 ) THEN
            WRITE(IPHEL,'(5X,A4,A2)')XS_name,'Xs'
           ELSE
            WRITE(IPHEL,'(5X,A4,A2)',advance='no')XS_name,'Xs'
           ENDIF
          ENDDO
          DO IT=1, NGP
           IF(IT==1) THEN
            WRITE(IPHEL,'(6X,A,12X,A,I1,A)',advance='no')
     1      'Label E','.-.-E-',IT,'-.'
           ELSE IF(IT==NGP .OR. IT==8 ) THEN
            WRITE(IPHEL,'(3X,A,I1,A)')
     1      '.-.-E-',IT,'-.'
           ELSE
            WRITE(IPHEL,'(3X,A,I1,A)',advance='no')
     1      '.-.-E-',IT,'-.'
           ENDIF
          ENDDO

          ! STORE XS DATA IN HELIOS.DRA FILE
          ! LOOP OVER BURNUP POINTS
          DO NB=1, NBU
            WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1      ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
            DO NG=1, NGP
              IF(NG == 1) THEN
               WRITE(IPHEL,'(ES12.5E2)',advance='no') XS(NG,XST,NB)
              ELSE IF(NG.NE.NGP) THEN
               WRITE(IPHEL,'(ES12.5E2)',advance='no') XS(NG,XST,NB)
              ELSE
                WRITE(IPHEL,'(ES12.5E2)') XS(NG,XST,NB)
              ENDIF
            ENDDO ! NG
          ENDDO ! NB

  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
          WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
         ENDIF
      ENDDO ! XST
      END

      SUBROUTINE SETADF(  IPHEL, BRANCH, ITBRAN,   ADF,    NADF,   NGP,
     >                      NBU,   BURN, IPRINT,ADF_T,    FLXR,   FLXL,
     >                     CURL,                                 CURR )
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the cross sections for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* ADF        Assembly discontinuity factor
* NADF       number of Assembly discontinuity factor
* NGP        number of energy groups
* NBU        number of burnup points
* BURN       set of burnup points
* IPRINT  control the printing on screen
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,NGP,NADF,ITBRAN,NIT,IPRINT
      REAL ADF(NADF,NGP,NBU),BURN (NBU)
      REAL FLXR(NGP,NBU),FLXL(NGP,NBU)
      REAL CURR(NGP,NBU),CURL(NGP,NBU), BNO(NGP,NBU)
      CHARACTER*3 ADF_T
      CHARACTER BRANCH*4
*----
*  LOCAL ARGUMENTS
*----
      INTEGER IT,ITA
      REAL ADF_TMP(NADF,NGP,NBU)
      CHARACTER*4 BOUND
      CHARACTER*12 LABEL,XSPRIN

      IF(IPRINT>5) WRITE(6,*) 'SETADF: RECOVER ADF INFO'
      IF (ADF_T.EQ.'DRA') THEN
      NIT=0
      IF((NADF.NE.1) .AND. (NADF.NE.4)) THEN
          WRITE(6,*) "NUMBER OF ADF : ",NADF
          CALL XABORT (" NUMBER OF ADF MUST BE 4 (SEL/GET/DRA) OR 1 "
     >  //"(DRA)")
      ELSE IF(NADF == 4) THEN
         ! CASE FOR SEL OR GET ADF
         ! REARRANGEMENT OF ADF ORDER TO MATCH HELIOS iN CASE OD SEL OR
         ! GET ADF
         ! SAPHYB SURF => SIDE
         !        1        N
         !        2        E
         !        3        S
         !        4        W
         ! HELIOS SURF => SIDE
         !        1        W
         !        2        S
         !        3        E
         !        4        N

          ADF_TMP(:,:,:)=ADF(:,:,:)
          ADF(1,:,:)=ADF_TMP(4,:,:)
          ADF(2,:,:)=ADF_TMP(3,:,:)
          ADF(3,:,:)=ADF_TMP(2,:,:)
          ADF(4,:,:)=ADF_TMP(1,:,:)
      ENDIF
      NIT =   NGP*NADF
      ! LABEL FOR XS TYPE : ADF
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_PRIN %SDF 2'
      WRITE(IPHEL,*)'2) Meaning : (F-.-E-G-.) F-Face, G-Group'
      IF(NADF==4) THEN
         WRITE(IPHEL,*)'3)          F=1/2/3/4 denotes W/S/E/N Side'
      ELSE
         WRITE(IPHEL,*)'3)          F=1 denotes average ADF'
      ENDIF

      CALL SET_RIEGO(IPHEL)

      !  LOOP OVER ENERGY GROUPS
      !  CREATION OF LABEL FOR CROSS SECTIONS
      ngrp=1
      nsurf=0
      DO ITA=1,NIT,7  ! ITA
       NITTMP=MIN(NIT-ITA+1,7)
       ngrpb=ngrp
       nsurfb=nsurf
       DO IT=1,NITTMP
          IF((IT==1).AND.(IT.LT.NITTMP)) THEN
            WRITE(IPHEL,'(30X,A6)',advance='no') 'SideDF'
          ELSEIF((IT==1).AND.(IT.EQ.NITTMP)) THEN
            WRITE(IPHEL,'(30X,A6)') 'SideDF'
          ELSE IF(IT.EQ.NITTMP) THEN
            WRITE(IPHEL,'(6X,A6)')'SideDF'
          ELSE
            WRITE(IPHEL,'(6X,A6)',advance='no')'SideDF'
          ENDIF
       ENDDO

       DO IT=1,NITTMP
         nsurf=nsurf+1
         IF(nsurf.GT.NADF) THEN
           nsurf=1
           ngrp=ngrp+1
         ENDIF
         IF((IT==1).AND.(IT.LT.NITTMP)) THEN
           WRITE(IPHEL,'(6X,A,14X,I1,A,I1,A)',advance='no')
     >        'Label E',nsurf,'-.-E-',ngrp,'-.'
         ELSEIF((IT==1).AND.(IT.EQ.NITTMP)) THEN
           WRITE(IPHEL,'(6X,A,14X,I1,A,I1,A)')
     >        'Label E',nsurf,'-.-E-',ngrp,'-.'
         ELSE IF(IT.EQ.NITTMP) THEN
           WRITE(IPHEL,'(3X,I1,A,I1,A)') nsurf,'-.-E-',ngrp,'-.'
         ELSE
           WRITE(IPHEL,'(3X,I1,A,I1,A)',advance='no')
     >         nsurf,'-.-E-',ngrp,'-.'
         ENDIF
       ENDDO

      ! STORE XS DATA IN HELIOS.DRA FILE
      ! LOOP OVER BURNUP POINTS

       DO NB=1,NBU
         ngrp=ngrpb
         nsurf=nsurfb
         WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1   ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         DO IT=1,NITTMP
           nsurf=nsurf+1
           IF(nsurf.GT.NADF) THEN
             nsurf=1
             ngrp=ngrp+1
           ENDIF
! in xs_helios_read.f90
! l1015    READ(XS_set_unit,hfnF5) rvector(1:RIEGO%how_many_data)
! in xs_heliosM.f90
! l104         hfnF5='(  X,8F13.5)        '
           IF(IT.EQ.NITTMP) THEN
             WRITE(IPHEL,'(5X,F7.5)') ADF(nsurf,ngrp,NB)
           ELSE
             WRITE(IPHEL,'(5X,F7.5)',advance='no') ADF(nsurf,ngrp,NB)
           ENDIF
         ENDDO
       ENDDO

       WRITE(IPHEL,*)
      ENDDO


      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
      ELSE IF (ADF_T.EQ.'GEN') THEN
       DO I=1,4
       SELECT CASE (I)
       CASE(1)
        XSPRIN='%PHW 1'
        BOUND='West'
        LABEL='FluxWest'
        BNO=FLXL
       CASE(2)
        XSPRIN='%PHE 1'
        BOUND='East'
        LABEL='FluxEast'
        BNO=FLXR
       CASE(3)
        XSPRIN='%JNW 1'
        BOUND='West'
        LABEL='JnetWest'
        BNO=CURL
       CASE(4)
        XSPRIN='%JNE 1'
        BOUND='East'
        LABEL='JnetEast'
        BNO=CURR
       END SELECT
       WRITE(IPHEL,'()')
       WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
       WRITE(IPHEL,*) 'Labels Array    : PrintVector'
       WRITE(IPHEL,*) 'List Title(s)  1) %XS_PRIN ',XSPRIN
       WRITE(IPHEL,'(16X,3A)')'2) Meaning : (E-.-E-G-.) ',
     >  BOUND,'-Face, G-Group'

       CALL SET_RIEGO(IPHEL)
       WRITE(IPHEL,'(31X,A8,4X,A8)') LABEL,LABEL
       WRITE(IPHEL,'(18X,A)') 'Label E   1-.-E-1-.   1-.-E-2-.'
       DO NB=1,NBU
         WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1   ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
        WRITE(IPHEL,'(ES12.5E2,1X,ES11.4E2)') BNO(:,NB)
        WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
       ENDDO
       ENDDO
      ENDIF
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      END

      SUBROUTINE SETCDF(  IPHEL, BRANCH, ITBRAN,   CDF,    NCDF,   NGP,
     >                      NBU,   BURN, IPRINT                       )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the cross sections for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* CDF        Corner discontinuity factor
* NCDF       number of corner discontinuity factor
* NGP        number of energy groups
* NBU        number of burnup points
* BURN       set of burnup points
* IPRINT  control the printing on screen
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,NGP,NCDF,ITBRAN,NIT,IPRINT
      REAL CDF(NCDF,NGP,NBU),BURN (NBU)
      CHARACTER BRANCH*4
*----
*  LOCAL ARGUMENTS
*----
      INTEGER IT,ITA

      IF(IPRINT>5) WRITE(6,*) 'SETCDF: RECOVER CDF INFO'
      NIT =   NGP*NCDF

      ! LABEL FOR XS TYPE : CDF
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_PRIN %CDF'
      WRITE(IPHEL,*)'2) Meaning : (F-.-E-G-.) F-Face, G-Group'
      IF(NCDF==1) THEN
         WRITE(IPHEL,*)'3)           F=1 denotes average CDF'
      ELSE
         WRITE(IPHEL,*)'3)         F= custom'
      ENDIF

      CALL SET_RIEGO(IPHEL)

      !  LOOP OVER ENERGY GROUPS
      !  CREATION OF LABEL FOR CROSS SECTIONS
      ngrp=1
      nsurf=0
      DO ITA=1,NIT,7  ! ITA
       NITTMP=MIN(NIT-ITA+1,7)
       ngrpb=ngrp
       nsurfb=nsurf
       DO IT=1,NITTMP
          IF((IT==1).AND.(IT.LT.NITTMP)) THEN
            WRITE(IPHEL,'(30X,A6)',advance='no') 'CornDF'
          ELSEIF((IT==1).AND.(IT.EQ.NITTMP)) THEN
            WRITE(IPHEL,'(30X,A6)') 'CornDF'
          ELSE IF(IT.EQ.NITTMP) THEN
            WRITE(IPHEL,'(6X,A6)')'CornDF'
          ELSE
            WRITE(IPHEL,'(6X,A6)',advance='no')'CornDF'
          ENDIF
       ENDDO

       DO IT=1,NITTMP
         nsurf=nsurf+1
         IF(nsurf.GT.NCDF) THEN
           nsurf=1
           ngrp=ngrp+1
         ENDIF
         IF((IT==1).AND.(IT.LT.NITTMP)) THEN
           WRITE(IPHEL,'(6X,A,15X,I1,A,I1,A)',advance='no')
     >        'Label E',nsurf,'-.-E-',ngrp,'-.'
         ELSEIF((IT==1).AND.(IT.EQ.NITTMP)) THEN
           WRITE(IPHEL,'(6X,A,15X,I1,A,I1,A)')
     >        'Label E',nsurf,'-.-E-',ngrp,'-.'
         ELSE IF(IT.EQ.NITTMP) THEN
           WRITE(IPHEL,'(2X,I1,A,I1,A)') nsurf,'-.-E-',ngrp,'-.'
         ELSE
           WRITE(IPHEL,'(2X,I1,A,I1,A)',advance='no')
     >         nsurf,'-.-E-',ngrp,'-.'
         ENDIF
       ENDDO

      ! STORE XS DATA IN HELIOS.DRA FILE
      ! LOOP OVER BURNUP POINTS


       DO NB=1,NBU
         ngrp=ngrpb
         nsurf=nsurfb
         WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1   ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         DO IT=1,NITTMP
           nsurf=nsurf+1
           IF(nsurf.GT.NCDF) THEN
             nsurf=1
             ngrp=ngrp+1
           ENDIF
! in xs_helios_read.f90 l1015    READ(XS_set_unit,hfnF5) rvector(1:RIEGO
! in xs_heliosM.f90 l104         hfnF5='(  X,8F13.5)        '
           IF(IT.EQ.NITTMP) THEN
             WRITE(IPHEL,'(5X,F7.5)') CDF(nsurf,ngrp,NB)
           ELSE
             WRITE(IPHEL,'(5X,F7.5)',advance='no') CDF(nsurf,ngrp,NB)
           ENDIF
         ENDDO
       ENDDO

       WRITE(IPHEL,*)
      ENDDO


  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'

      END

      SUBROUTINE SETGFF(  IPHEL, BRANCH, ITBRAN,   GFF,  NCOLA,  NROWA,
     >                    NPART,    NGP,    NBU,   BURN, NGFF , IPRINT,
     >                     VERS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the cross sections for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* GFF        Group form factor
* NCOLA      number of pin in assembly along x-axis
* NROWA      number of pin in assembly along y-axis
* NPART      symmetry level of assembly
*               0        1        2        3
*               whole    half     quarter  eight
* PARCS Version 32.17 and GenPMAXS 6.1
*               123XXXX  1......  123X...  1......
*               XXXXXXX  23.....  XXXX...  23.....
*               XXXXXXX  XXX....  XXXX...  XXX....
*               XXXXXXX  XXXX...  XXXn...  XXXn...
*               XXXXXXX  XXXXX..  .......  .......
*               XXXXXXX  XXXXXX.  .......  .......
*               XXXXXXn  XXXXXXn  .......  .......
*              Note: Helios format is different from the documentation
*              provided in GenPMAXS.
* Version 32.18 and GenPMAXS 6.2
*               123XXXX  1......  .......  .......
*               XXXXXXX  23.....  .......  .......
*               XXXXXXX  XXX....  .......  .......
*               XXXXXXX  XXXX...  ...123X  ...1...
*               XXXXXXX  XXXXX..  ...XXXX  ...23..
*               XXXXXXX  XXXXXX.  ...XXXX  ...XXX.
*               XXXXXXn  XXXXXXn  ...XXXn  ...XXXn
*              Note: Helios format is the same as in the documentation
*              provided in GenPMAXS.
* NGP        number of energy groups
* NBU        number of burnup points
* BURN       set of burnup points
* IPRINT  control the printing on screen
* VERS       version of PARCS to be used
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,NGP,ITBRAN,NIT,IPRINT,NGFF
      REAL GFF(NCOLA,NROWA,NGP,NBU),BURN (NBU),VERS
      CHARACTER BRANCH*4
*----
*  LOCAL ARGUMENTS
*----
      INTEGER IT,ITA,ipxn,ipyn

      IF(IPRINT>5) WRITE(6,*) 'SETGFF: RECOVER GFF INFO'
      NIT   = NGP*NCOLA*NROWA
      NPIN2 = NCOLA*NROWA
      NCOLA2= 1
      ipxn=1
      ipyn=1
      IF((NPART.GE.1).AND.(NCOLA.NE.NROWA))THEN
        CALL XABORT('@D2PBRA: NPART > 0 and  NCOLA.NE.NROWA')
      ENDIF
      IF(NPART.EQ.1)THEN
        NIT=NGP*NCOLA*(NCOLA+1)/2
        NPIN2 = NCOLA*(NCOLA+1)/2
      ELSEIF(NPART.EQ.2)THEN
        NCOLA2=CEILING(REAL(NCOLA)/2)
        NIT=NGP*NCOLA2*NCOLA2
        NPIN2 = NCOLA2*NCOLA2
      ELSEIF(NPART.EQ.3)THEN
        NCOLA2=CEILING(REAL(NCOLA)/2)
        NIT=NGP*NCOLA2*(NCOLA2+1)/2
        NPIN2 = NCOLA2*(NCOLA2+1)/2
      ENDIF
      IF((VERS.GE.3.2018).AND.(NPART.GE.2))THEN
        ipxn=CEILING(REAL(NCOLA)/2)
        ipyn=CEILING(REAL(NCOLA)/2)
        NCOLA2=NCOLA
      ENDIF
      IF (NGFF.NE.NPIN2) THEN
       WRITE (6,*) '@D2PBRA: INCOHERENT NUMBER OF GFF IN MCO (',
     > NGFF,') AND COMPUTED PART OF ASSEMBLY (PART =',
     > NPART,').'
      CALL XABORT ('')
      ENDIF
      ! LABEL FOR XS TYPE: GFF
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_GFF %GFF 1 '
      WRITE(IPHEL,*)'2) Meaning : (F-.-E-G-.) F-Face, G-Group'
      WRITE(IPHEL,*)'3)         F=1 to NPIN*NPIN average GFF'

      CALL SET_RIEGO(IPHEL)

      !  LOOP OVER ENERGY GROUPS
      !  CREATION OF LABEL FOR CROSS SECTIONS
      ngrp=1
      nsurf=0
      ipx=ipxn-1
      ipy=ipyn
      DO ITA=1,NIT,7  ! ITA
       NITTMP=MIN(NIT-ITA+1,7)
       ngrpb=ngrp
       nsurfb=nsurf
       ipxb=ipx
       ipyb=ipy
       DO IT=1,NITTMP
          IF((IT==1).AND.(IT.LT.NITTMP)) THEN
            WRITE(IPHEL,'(33X,A6)',advance='no') 'GNorRR'
          ELSEIF((IT==1).AND.(IT.EQ.NITTMP)) THEN
            WRITE(IPHEL,'(33X,A6)') 'GNorRR'
          ELSE IF(IT.EQ.NITTMP) THEN
            WRITE(IPHEL,'(8X,A6)')'GNorRR'
          ELSE
            WRITE(IPHEL,'(8X,A6)',advance='no')'GNorRR'
          ENDIF
       ENDDO
       DO IT=1,NITTMP
         nsurf=nsurf+1
         IF(nsurf.GT.NPIN2) THEN
           nsurf=1
           ngrp=ngrp+1
         ENDIF
         IF((IT==1).AND.(IT.LT.NITTMP)) THEN
           WRITE(IPHEL,'(6X,A,12X,I3,A,I1,A)',advance='no')
     >        'Label E',nsurf,'-.-E-',ngrp,'-.'
         ELSEIF((IT==1).AND.(IT.EQ.NITTMP)) THEN
           WRITE(IPHEL,'(6X,A,12X,I3,A,I1,A)')
     >        'Label E',nsurf,'-.-E-',ngrp,'-.'
         ELSE IF(IT.EQ.NITTMP) THEN
           WRITE(IPHEL,'(1X,I3,A,I1,A)') nsurf,'-.-E-',ngrp,'-.'
         ELSE
           WRITE(IPHEL,'(1X,I3,A,I1,A)',advance='no')
     >         nsurf,'-.-E-',ngrp,'-.'
         ENDIF
       ENDDO
      ! STORE XS DATA IN HELIOS.DRA FILE
      ! LOOP OVER BURNUP POINTS


       DO NB=1,NBU
         ngrp=ngrpb
         nsurf=nsurfb
         ipx=ipxb
         ipy=ipyb
         WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1   ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         DO IT=1,NITTMP
           ipx=ipx+1
           IF(((NPART.EQ.0).AND.(ipx.GT.NCOLA)).OR.
     >        ((NPART.EQ.2).AND.(ipx.GT.NCOLA2)).OR.
     >        (((NPART.EQ.1).OR.(NPART.EQ.3)).AND.(ipx.GT.ipy)))THEN
             ipx=ipxn
             ipy=ipy+1
           ENDIF
           nsurf=nsurf+1
           IF(nsurf.GT.NPIN2) THEN
             nsurf=1
             ngrp=ngrp+1
             ipy=ipxn
             ipx=ipyn
           ENDIF
! in xs_helios_read.f90 l1015    READ(XS_set_unit,hfnE4) rvector(1:RIEGO
! in xs_heliosM.f90 l104         hfnF5='(  X,8F13.5)        '
!                   l114         hfnE4=hfnE5
!                   l115         hfnE4(11:11)='4'
           IF(IT.EQ.NITTMP) THEN
             WRITE(IPHEL,'(5X,F7.4)') GFF(ipx,ipy,ngrp,NB)
           ELSE
             WRITE(IPHEL,'(5X,F7.4)',advance='no') GFF(ipx,ipy,ngrp,NB)
           ENDIF
         ENDDO
       ENDDO

       WRITE(IPHEL,*)
      ENDDO


  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'


      END

      SUBROUTINE SETSCT(IPHEL,BRANCH,ITBRAN,SCAT,NGP,NBU,BURN,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the scattering cross sections for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* SCAT       table of elements of scattering matrix
* NGP        number of energy groups
* NBU        number of burnup points
* BURN       set of burnup points
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,NGP,ITBRAN,IPRINT
      REAL SCAT(NGP*NGP,NBU),BURN (NBU)

      CHARACTER BRANCH*2

*----
*  LOCAL ARGUMENTS
*----
      INTEGER IT,G,I
      REAL SCATTMP(8,NBU)
      CHARACTER*45 LABEL
      CHARACTER*45 LABELE
      CHARACTER*210 :: TOTLABELE = ''
      CHARACTER*210 :: TOTLABEL = ''

      IF(IPRINT>5) WRITE(6,*) 'SETSCT: WRITE SCATTERING INFO'

      ! LABEL FOR SCATTERING XS
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,110) '%SCT'
      WRITE(IPHEL,*)'Meaning : (.-.-E-G-O) From O to G-th Group scat'

      CALL SET_RIEGO(IPHEL)
      IT=1
      ITT=1
      I=0
      ! CREATION OF HEADER FOR SCATTERING BLOCK IN HELIOS.DRA FILE

          DO G=1,NGP
           DO J=1, NGP
             IF (IT==1) THEN
              TOTLABELE = ''
              TOTLABEL = ''
              WRITE(LABELE,'(6X,A7,14X)') 'Label E'
              TOTLABELE=TOTLABELE(1:len( trim(TOTLABELE) ))
     1        // LABELE
             ENDIF
             IF (IT==1) THEN
              WRITE(LABEL,'(25X,A)')'ScattMatrix'
              WRITE(LABELE,'(12X,A,I1,A,I1)')
     1        '1-.-E-',G,'-',J
             ELSE
              WRITE(LABEL,'(1X,A)')'ScattMatrix'
              WRITE(LABELE,'(3X,A,I1,A,I1)')
     1        '1-.-E-',G,'-',J
             ENDIF
             SCATTMP(IT,:)=SCAT(ITT,:)
             TOTLABEL=TOTLABEL(1:len( trim(TOTLABEL) ))
     >        //LABEL
             TOTLABELE=TOTLABELE(1:len( trim(TOTLABELE) ))
     >        //LABELE

             IF ((IT==8).OR.(ITT==NGP*NGP)) THEN
              WRITE(IPHEL,'(A)') TOTLABEL
              WRITE(IPHEL,'(A)') TOTLABELE
              DO NB=1, NBU
               WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1         ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
               WRITE(IPHEL,'(8(ES12.5E2))')SCATTMP(1:IT,NB)
              ENDDO
              WRITE (IPHEL,*)
              TOTLABELE = ''
              TOTLABEL = ''
              IT=1
             ELSE
              IT=IT+1
             ENDIF

             ITT=ITT+1
            ENDDO
           ENDDO

  !     DO JT=1, NGP
  !      IF(IT==1 .and. JT==1) THEN
  !         WRITE(IPHEL,'(28X,A)',advance='no') 'ScattMatrix'
  !      ELSE IF((IT==(NGP).and.JT==(NGP)) .OR. JT==7 ) THEN
  !         WRITE(IPHEL,'(3X,A)')'ScattMatrix'
  !      ELSE
  !         WRITE(IPHEL,'(3X,A)',advance='no')'ScattMatrix'
  !      ENDIF
  !     ENDDO
  !   DO IT=1, NGP
  !     DO JT=1, NGP
  !       IF(IT==1 .and. JT==1) THEN
  !         WRITE(IPHEL,'(6X,A,14X,A,I1,A,I1)',advance='no')
  !  1       'Label E','1-.-E-',JT,'-',IT
  !       ELSE IF((IT==(NGP).and.JT==(NGP)) .OR. JT==8 ) THEN
  !         WRITE(IPHEL,'(5X,A,I1,A,I1)')
  !  1       '1-.-E-',JT,'-',IT
  !       ELSE
  !         WRITE(IPHEL,'(5X,A,I1,A,I1)',advance='no')
  !  1       '1-.-E-',JT,'-',IT
  !       ENDIF
  !     ENDDO
  !   ENDDO
  !
  !   DO NB=1, NBU
  !       WRITE(IPHEL,220,advance='no') NB,'t',BRANCH(1:2),
  !  1    ITBRAN,'(s',BRANCH(1:2),ITBRAN,'):',NINT(BURN(NB))
  !       DO IG=1, NGP*NGP
  !         IF(IG == 1) THEN
  !          WRITE(IPHEL,'(3X,ES11.5E2)',advance='no') SCAT(IG,NB)
  !         ELSE IF(IG.EQ.NGP*NGP) THEN
  !           WRITE(IPHEL,'(3X,ES11.5E2)') SCAT(IG,NB)
  !         ELSE IF(IG.EQ.8  ) THEN
  !          WRITE(IPHEL,'(3X,ES11.5E2)') SCAT(IG,NB)
  !         ELSE
  !          WRITE(IPHEL,'(3X,ES11.5E2)',advance='no') SCAT(IG,NB)
  !         ENDIF
  !       ENDDO
  !   ENDDO

  110 FORMAT(28H List Title(s)  1) %XS_SCT  ,A)
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
      END

      SUBROUTINE SETND(IPHEL,BRANCH,ITBRAN,ND,NBU,BURN,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the scattering cross sections for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* NGP        number of energy groups
* NBU        number of burnup points
* ND         number densities for Xenon and samarium : KEFF , KINF, B2
* BURN       set of burnup points
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,ITBRAN,IPRINT
      REAL ND(2,NBU),BURN (NBU)
      CHARACTER BRANCH*2
*----
*  LOCAL ARGUMENTS
*----
      INTEGER NB

      IF(IPRINT>5) WRITE(6,*) 'SETND: WRITE HEADER FOR XENON DENSITY'

      ! CREATION OF HEADER FOR XENON DENSITY
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_XESM %XEND'
      WRITE(IPHEL,*)'Meaning : Xe-135 Number Density [/cm.barn]'

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,'(31X,A)') 'nXe'
      WRITE(IPHEL,'(6X,A,12X,A)') 'Label E','1-1-E-.-.'

      ! LOOP OVER BUNRNUP
      DO NB=1,NBU
              WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1        ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
        WRITE(IPHEL,'(ES12.5E2)') ND(1,NB)
      ENDDO

      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'

      ! CREATION OF HEADER FOR SAMARIUM DENSITY
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_XESM %SMND'
      WRITE(IPHEL,*)'Meaning : SM-149 Number Density [/cm.barn]'

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,'(27X,A)')  'nSm'
      WRITE(IPHEL,'(6X,A,12X,A)') 'Label E','1-1-E-.-.'

      ! LOOP OVER BUNRNUP
      DO NB=1,NBU
              WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1        ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         WRITE(IPHEL,'(ES12.5E2)') ND(2,NB)
      ENDDO
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'

  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      END

      SUBROUTINE SETDIV(IPHEL,BRANCH,ITBRAN,DIV,NBU,BURN,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the scattering cross sections for a banch
*  (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*parameters: input
* IPHEL      file unit of the HELIOS.dra file
* BRANCH     nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN     index of the current branch
* NGP     number of energy groups
* NBU        number of burnup points
* DIV        conttnent of DIV table : KEFF , KINF, B2
* BURN       set of burnup points
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IPHEL,NBU,ITBRAN,IPRINT
      REAL DIV(3,NBU),BURN (NBU)
      CHARACTER BRANCH*2
*----
*  LOCAL ARGUMENTS
*----
      INTEGER NB
      REAL M2

      IF(IPRINT>5) WRITE(6,*) 'SETDIV: WRITE HEADER FOR DIVERS INFO'

      ! CREATION OF HEADER FOR DIVERS INFO (B2, KEFF, KINF)
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_PRIN %KINF'
      WRITE(IPHEL,*)' Meaning : K-eff, K-inf, M^2, B^2 [cm^-2]  '

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,'(27X,A,10X,A,6X,A,6X,A)') 'K-EFF','KINF',
     1  'MigrArea','CritArea'
      WRITE(IPHEL,'(6X,A,12X,A,5X,A,5X,A,5X,A)')
     1       'Label E','.-.-E-.-.','.-.-E-.-.','.-.-E-.-.','.-.-E-.-.'
      ! LOOP OVER BURNUP POINTS
      DO NB=1,NBU
         M2=(DIV(2,NB)-1)/(DIV(3,NB))
              WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1        ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         WRITE(IPHEL,'(5X,F7.5,5X,F7.5,ES12.5E2,ES12.5E2)')
     1    DIV(1,NB),DIV(2,NB),M2,DIV(3,NB)
      ENDDO
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'

  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
       END

       SUBROUTINE SETTH(  IPHEL, BRANCH, ITBRAN,  BURN,     NBU, JOBOPT,
     >                      NGP,  IPDAT, IPRINT                        )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write in the HELIOS.dra file the invaraint TH DATA for a banch
* (including all burnup points).
* This routine write sequentially the HELIOS.dra file, branch after
* branch.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      TYPE(C_PTR) IPDAT,IPTH,KPTH
      INTEGER :: DIM_LAMBDA = 6
      INTEGER NGP,ILONG
      INTEGER ITYLCM,NBU,ITBRAN
      CHARACTER (len=4) BRANCH,DLAY
      CHARACTER JOBOPT(16)
      INTEGER :: BU = 1
      REAL BURN(NBU)
      REAL YLDXe(NBU),YLDPm(NBU),YLDI(NBU)
      REAL OVERV(NGP,NBU),CHI(NGP,NBU),LAMBDA(6,NBU),BETA(6,NBU)
      LOGICAL :: LAMB = .FALSE.
      LOGICAL :: LCHI = .FALSE.
      LOGICAL :: LYLD = .FALSE.
      LOGICAL :: LINV = .FALSE.
      LOGICAL :: LBET = .FALSE.

      IF(IPRINT>5) WRITE(6,*) 'SETTH: WRITE TH DATA'

      ! RECOVER FLAG INFORMATION
      IF(JOBOPT(5)=='T') LCHI = .TRUE.
      IF(JOBOPT(7)=='T') LINV = .TRUE.
      IF(JOBOPT(9)=='T') LYLD = .TRUE.
      IF(JOBOPT(13)=='T') LAMB = .TRUE.
      IF(JOBOPT(12)=='T') LBET = .TRUE.

      IF(NGP>2)THEN
        CALL XABORT('@D2P: NGP > 2 NOT IMPLEMENTED FOR T/H BLOCK')
      ENDIF

      CALL LCMSIX(IPDAT,' ',0)
      IPTH=LCMGID(IPDAT,'TH_DATA')

      DO BU=1,NBU
         KPTH=LCMDIL(IPTH,BU)

         IF(LCHI) THEN
          IF(BU ==1) THEN
           CALL LCMLEN(KPTH,'CHI',ILONG,ITYLCM)
           IF  (ILONG .NE. NGP) THEN
            CALL XABORT (' MORE THAN 2 (NGP) VALUES FOR CHI RECORD')
           ENDIF
          ENDIF
          CALL LCMGET(KPTH,'CHI',CHI(1:NGP,BU))

         ENDIF

         IF(LINV) THEN
          CALL LCMGET(KPTH,'OVERV',OVERV(1:NGP,BU))

         ENDIF

         IF(LYLD) THEN
          CALL LCMGET(KPTH,'YLDPm',YLDPm(BU))
          CALL LCMGET(KPTH,'YLDXe',YLDXe(BU))
          CALL LCMGET(KPTH,'YLDI',YLDI(BU))

         ENDIF

         IF(LAMB)THEN
          IF(BU == 1) THEN
           CALL LCMLEN(KPTH,'LAMBDA',ILONG,ITYLCM)
           IF  (ILONG .NE. DIM_LAMBDA) THEN
            CALL XABORT('MORE THAN 6 (NDLAY) VALUES FOR LAMBDA RECORD')
           ENDIF
          ENDIF
          CALL LCMGET(KPTH,'LAMBDA',LAMBDA(1:DIM_LAMBDA,BU))
         ENDIF

         IF(LBET)THEN
          IF(BU == 1) THEN
           CALL LCMLEN(KPTH,'BETA',ILONG,ITYLCM)
           IF  (ILONG .NE. DIM_LAMBDA) THEN
            CALL XABORT('MORE THAN 6 (NDLAY) VALUES FOR BETA RECORD')
           ENDIF
          ENDIF
          CALL LCMGET(KPTH,'BETA',BETA(1:DIM_LAMBDA,BU))
         ENDIF
      ENDDO

      IF(LCHI) CALL SET_CHI(IPHEL,BRANCH,ITBRAN,BURN,CHI,NGP,NBU)
      IF(LINV) CALL SET_OVERV(IPHEL,BRANCH,ITBRAN,BURN,OVERV,NGP,NBU)
      IF(LYLD) CALL SET_YIELD(IPHEL,BRANCH,ITBRAN,BURN,YLDPm,YLDXe,
     >  YLDI,NBU)
      IF(LAMB) THEN
        DLAY='LAMB'
        CALL SET_DLAY(IPHEL,BRANCH,ITBRAN,BURN,LAMBDA,DIM_LAMBDA,DLAY,
     >  NBU)
      ENDIF
      IF(LBET) THEN
        DLAY='BETA'
        CALL SET_DLAY(IPHEL,BRANCH,ITBRAN,BURN,BETA,DIM_LAMBDA,DLAY,NBU)
      ENDIF
      END

      SUBROUTINE SET_CHI(IPHEL,BRANCH,ITBRAN,BURN,CHI,DIM_CHI,NBU)
      INTEGER DIM_CHI,NBU,ITBRAN,NB
      REAL CHI(DIM_CHI,NBU),BURN(NBU)
      CHARACTER (len=4) BRANCH

      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_PRIN %CHI'
      WRITE(IPHEL,*) 'Meaning :(.-.-E-G-.) G-th Group Fission Spect'

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,'(31X,A3,11X,A3)') 'chi','chi'
      WRITE(IPHEL,'(6X,A,12X,A,5X,A)') 'Label E','1-.-E-1-.','1-.-E-2-.'
      ! LOOP OVER burnup points
      DO NB=1,NBU
         WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1   ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         WRITE(IPHEL,'(ES12.5E2,ES12.5E2)')
     1   CHI(1:DIM_CHI,NB)
      ENDDO
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      END

      SUBROUTINE SET_OVERV(IPHEL,BRANCH,ITBRAN,BURN,OVERV,NG,NBU)
      INTEGER NG,NBU,ITBRAN,NB
      REAL OVERV(NG,NBU),BURN(NBU)
      CHARACTER (len=4) BRANCH

      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      WRITE(IPHEL,*) 'List Title(s)  1) %XS_PRIN %VEL'
      WRITE(IPHEL,*) 'Meaning         :'
      WRITE(IPHEL,*) '(.-.-E-G-.) G-th Group Neutron Velocity [m/s]'

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,'(31X,A3,11X,A3)') 'vel','vel'
      WRITE(IPHEL,'(6X,A,12X,A,5X,A)') 'Label E','.-.-E-1-.','.-.-E-2-.'
      ! LOOP OVER burnup points
      DO NB=1,NBU

         WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1   ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
         WRITE(IPHEL,'(ES12.5E2,ES12.5E2)')
     1   (1/(OVERV(1:NG,NB)))
      ENDDO
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      END

      SUBROUTINE SET_YIELD(IPHEL,BRANCH,ITBRAN,BURN,YLDPm,YLDXe,YLDI,
     1 NBU)
      INTEGER NBU,ITBRAN,NB,I,iXe,iPm,iI
      REAL YLDPm(NBU), YLDXe(NBU),YLDI(NBU),BURN(NBU)
      REAL YLD(NBU)
      CHARACTER (len=4) BRANCH
      CHARACTER (len=5) YIELD
      CHARACTER (len=6) MEANING
      CHARACTER (len=10 ) LABEL

      DO I=1, 3
         SELECT CASE (I)
         CASE(1)
           YIELD='YLDXE'
           MEANING='Xe-135'
           LABEL='YieldXe135'
           DO iXe=1,NBU
            YLD(iXe)=YLDXe(iXe)
           ENDDO
         CASE(2)
           YIELD='YLDID'
           MEANING=' I-135'
           LABEL=' YieldI135'
           DO iI=1,NBU
            YLD(iI)=YLDI(iI)
           ENDDO
         CASE(3)
           YIELD='YLDPM'
           MEANING='Pr-149'
           LABEL='YieldPm149'
           DO iPm=1,NBU
            YLD(iPm)=YLDPm(iPm)
           ENDDO
         END SELECT

         WRITE(IPHEL,'()')
         WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
         WRITE(IPHEL,*) 'Labels Array    : PrintVector'
         WRITE(IPHEL,*) 'List Title(s)  1) %XS_XESM %',YIELD
         WRITE(IPHEL,*) 'Meaning : Effective ,',MEANING,' Yield'

         CALL SET_RIEGO(IPHEL)

         WRITE(IPHEL,'(29X,A10)') LABEL
         WRITE(IPHEL,'(6X,A,12X,A)') 'Label E','1-.-E-1-.'
          ! LOOP OVER burnup points
         DO NB=1,NBU

           WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1     ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
           WRITE(IPHEL,'(ES12.5E2)') YLD(NB)
         ENDDO
         WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      ENDDO
      END

      SUBROUTINE SET_DLAY(IPHEL,BRANCH,ITBRAN,BURN,VECT,DIM_LAMBDA,
     1 DLAY,NBU)
      INTEGER DIM_LAMBDA,NBU,ITBRAN,NB
      REAL VECT(DIM_LAMBDA,NBU),BURN(NBU)
      CHARACTER (len=4) BRANCH,DLAY
      CHARACTER (len=6) LABEL

      IF(DLAY.EQ.'LAMB') THEN
         LABEL="lambda"
      ELSE
         LABEL="beta  "
      ENDIF
      IF(DIM_LAMBDA.GT.8) THEN
         WRITE (6,*) "@D2PBRA: NB OF DELAY NEUTRON GROUPS:",DIM_LAMBDA
         CALL XABORT("MAX EIGHT DELAY NEUTRON GROUPS ARE ALLOWED")
      ENDIF
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : PrintVector'
      IF(LABEL=="lambda")THEN
         WRITE(IPHEL,*) 'List Title(s)  1) %XS_BETA %DCAYB 1'
         WRITE(IPHEL,*) 'Meaning : Decay Cst of the Delayed Neutron /s'
      ELSE
         WRITE(IPHEL,*) 'List Title(s)  1) %XS_BETA %BETA 1  '
         WRITE(IPHEL,*) 'Meaning : Delayed Neutron Fraction'
      ENDIF
      WRITE(IPHEL,*) '          (.-.-E-G-.) From 0 To 6-th Group'

      CALL SET_RIEGO(IPHEL)
      IF(DIM_LAMBDA.EQ.6) THEN
          WRITE(IPHEL,'(31X,A6,6X,A6,6X,A6,6X,A6,6X,A6,6X,A6)')
     >      LABEL,LABEL,LABEL,LABEL,LABEL,LABEL
          WRITE(IPHEL,200)
     >      'Label E','.-.-E-1-.','.-.-E-2-.','.-.-E-3-.','.-.-E-4-.',
     >      '.-.-E-5-.','.-.-E-6-.'
         ! LOOP OVER burnup points
          DO NB=1,NBU
           WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1     ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
           WRITE(IPHEL,'(ES12.5E2,5(ES12.5E2))')
     >      VECT(1:DIM_LAMBDA,NB)
          ENDDO
      ELSE IF(DIM_LAMBDA.EQ.8) THEN
         WRITE(IPHEL,'(26X,A6,6X,A6,6X,A6,6X,A6,6X,A6,6X,A6,6X,A6)')
     >   LABEL,LABEL,LABEL,LABEL,LABEL,LABEL,LABEL
         WRITE(IPHEL,210)
     >     'Label E','.-.-E-1-.','.-.-E-2-.','.-.-E-3-.','.-.-E-4-.',
     >     '.-.-E-5-.','.-.-E-6-.','.-.-E-7-.'
         ! LOOP OVER burnup points
         DO NB=1,NBU
           WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1     ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
           WRITE(IPHEL,'(ES12.5E2,6(ES12.5E2))') VECT(1:7,NB)
         ENDDO

         WRITE(IPHEL,*)

         WRITE(IPHEL,'(26X,A6)') 'lambda'
         WRITE(IPHEL,'(6X,A,12X,A)') 'Label E',LABEL
         DO NB=1, NBU
          WRITE(IPHEL,220,advance='no') NB,BRANCH(1:2),
     1    ITBRAN,' ',BRANCH(1:2),ITBRAN,':',NINT(BURN(NB))
          WRITE(IPHEL,'(ES12.5E2)') VECT(DIM_LAMBDA,NB)
         ENDDO
      ENDIF
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
  200 FORMAT(6X,A,12X,A,3X,A,3X,A,3X,A,3X,A,3X,A)
  210 FORMAT(6X,A,12X,A,3X,A,3X,A,3X,A,3X,A,3X,A,3X,A)
  220 FORMAT(I4,1X,A,I4,A,A,I4,A,I5)
      END
