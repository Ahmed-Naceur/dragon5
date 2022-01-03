*DECK D2PXS
      SUBROUTINE D2PXS (IPDAT,IPMIC,IPSAP,STAVEC,SIGNAT,MIXDIR,
     >                 JOBOPT,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover cross sections from a microlib object and write cross
* sections for one branch at a fixed burnup point in the INFO data
* block.
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of info data block
* IPSAP   address of the saphyb object
* IPMIC   address of the microlib object
* STAVEC  various parameters associated with the IPDAT structure
* SIGNAT  signature of the object containing cross sections
* MIXDIR  directory that contains homogeneous mixture information
* IPRINT  control the printing on screen
*
*Parameters: 
* JOBOPT  
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPMIC,IPSAP
      INTEGER STAVEC(40),IPRINT
      CHARACTER*12 SIGNAT,MIXDIR
*----
*  LOCAL VARIABLES
*----
      ! INDEX OF CURRENT VALUE FOR EACH STATE VARIABLES
      PARAMETER(NSTATE=40)
      INTEGER STAIDX (STAVEC(2)),ISTATE(NSTATE)
      INTEGER DIMSAP(50)
      INTEGER ITBRA,NSF,ITR
      INTEGER ::NREA = 0
      INTEGER :: NISO = 0
      INTEGER ::NMIL = 0
      INTEGER ::NBISO = 0
      INTEGER ::NANI = 0
      INTEGER ::NFISS = 0
      INTEGER :: NADD = 0
      INTEGER :: NBMIX = 0
      INTEGER :: NMAC = 0
      INTEGER :: NADRX = 0
      INTEGER :: NPAR = 0
      INTEGER :: NDEL = 0
      INTEGER :: ISPH = 0
      ! INDICATES THE END OF A BRANCH CALCULATION (REW=1), AND A
      ! DEFAULT MESHING (GRID)
      INTEGER REW,GRID
      ! NUMBER OF STATES VARIABLES
      INTEGER NVAR
      ! NUMBER OF BURNUP POINTS
      INTEGER NBU,NGP
      INTEGER :: NADF = 1
      INTEGER :: NCDF = 1
      INTEGER :: NGFF = 1
      INTEGER :: NPIN = 1
      INTEGER :: NTYPE = 1
      INTEGER FLAG
      INTEGER ICOR
      REAL    STATE(STAVEC(2)),BURN(STAVEC(4)),REFSTA(STAVEC(2)-1)
      ! DATSRC BLOCK OF INFO/GENPMAXS DIRECTORY
      REAL DATSRC(5),FLUX(STAVEC(1))
      ! STATE VARIABLE NAMES
      CHARACTER(len=12) STAVAR(STAVEC(2))
      CHARACTER JOBOPT(16)

      CHARACTER*4 BRANCH
      CHARACTER*3 ADF_T,CDF_T,GFF_T
      LOGICAL LABS(3),SCAT
      LOGICAL :: LADF = .FALSE.
      LOGICAL :: LCDF = .FALSE.
      LOGICAL :: LGFF = .FALSE.
      LOGICAL :: LXES = .FALSE.
      LOGICAL :: LDET = .FALSE.
      LOGICAL :: LTH  = .FALSE.
      LOGICAL :: LCOR = .FALSE.


      ! INITIALIZATION OF PARAMETERS
      NVAR=STAVEC(2)
      NBU=STAVEC(4)
      GRID=STAVEC(5)
      NGP=STAVEC(1)
      NSF=STAVEC(11)
      ICOR=STAVEC(22)

      ! RECOVER INFORMATION FROM INFO date block
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMGET(IPDAT,'FLAG',FLAG)
      CALL LCMGET(IPDAT,'DAT_SRC',DATSRC)


      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)

      IF (ICOR>0) LCOR=.TRUE.
      IF(JOBOPT(1)=='T') THEN
         CALL LCMGTC(IPDAT,'ADF_TYPE',3,1,ADF_T)
         LADF = .TRUE.
         IF((ADF_T.EQ.'SEL').OR.(ADF_T.EQ.'GET')) THEN
          STAVEC(13)=NSF
          STAVEC(14)=1
         ENDIF
         IF((ADF_T.EQ. 'DRA').OR.(ADF_T.EQ. 'GEN'))THEN
          STAVEC(13)=1
          CALL LCMSIX(IPMIC,' ',0)
          CALL LCMSIX(IPMIC,'MACROLIB',1)
          CALL LCMSIX(IPMIC,'ADF',1)
          CALL LCMGET(IPMIC,'NTYPE',STAVEC(14))
         ENDIF
         NADF=STAVEC(13)
         NTYPE=STAVEC(14)
      ENDIF

      IF(JOBOPT(2)=='T') LXES = .TRUE.
      IF(JOBOPT(8)=='T') LDET = .TRUE.
      IF((JOBOPT(5)=='T').OR.(JOBOPT(7)=='T').OR.
     >   (JOBOPT(9)=='T').OR.(JOBOPT(13)=='T')) THEN
        LTH =.TRUE.
      ENDIF

      IF(JOBOPT(10)=='T') THEN
         CALL LCMGTC(IPDAT,'CDF_TYPE',3,1,CDF_T)
         LCDF = .TRUE.
         IF(CDF_T.EQ. 'DRA')THEN
          CALL LCMSIX(IPMIC,' ',0)
          CALL LCMSIX(IPMIC,'MACROLIB',1)
          CALL LCMSIX(IPMIC,'ADF',1)
          CALL LCMGET(IPMIC,'NTYPE',STAVEC(14))
         ENDIF
         NCDF=STAVEC(15)
         NTYPE=STAVEC(14)
      ENDIF
      IF(JOBOPT(11)=='T') THEN
         CALL LCMGTC(IPDAT,'GFF_TYPE',3,1,GFF_T)
         LGFF = .TRUE.
         NGFF=STAVEC(16)
         NPIN=STAVEC(17)
      ENDIF

      IF(DATSRC(3).NE.0.0) THEN
        CALL LCMGET(IPDAT,'LABS',LABS)
        CALL LCMGET(IPDAT,'SCAT',SCAT)
      ENDIF

      CALL LCMGTC(IPDAT,'STATE_VAR',12,NVAR,STAVAR)
      CALL LCMGET(IPDAT,'BURN',BURN)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMGET(IPDAT,'REWIND',REW)
      CALL LCMGTC(IPDAT,'BRANCH',4,1,BRANCH)
      CALL LCMGET(IPDAT,'BRANCH_IT',ITBRA)

      CALL LCMGET(IPDAT,'STATE_INDEX',STAIDX)
      CALL LCMGET(IPDAT,'STATE',STATE)
      CALL LCMSIX(IPMIC,' ',0)
      CALL LCMGET(IPMIC,'STATE-VECTOR',ISTATE)

      NBISO=ISTATE(2)         ! NUMBER OF ISOTOPES
      NDEL=ISTATE(19)         ! NUMBER OF DELAYED NEUTRON GROUPS

      IF(NDEL.NE.STAVEC(7)) THEN
         WRITE(6,*) "@D2PXS: ERROR IN NUMBER OF DELAYED NEUTRON GROUPS"
         WRITE(6,*) "THE NUMBER OF DELAYED NEUTRON GROUPS IN SAP (",
     1   STAVEC(7),") IS DIFFERENT FROM MICROLIB (",NDEL,")"
         CALL XABORT('@D2PXS: DELAYED NEUTRON DATA ERROR')
      ENDIF

      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMSIX(IPMIC,' ',0)
      CALL LCMSIX(IPMIC,'MACROLIB',1)
      CALL LCMGET(IPMIC,'STATE-VECTOR',ISTATE)

      NBMIX=ISTATE(2)          ! NUMBER OF MIXTURESS
      NANI=ISTATE(3)           ! SCATTERING ANISOTROPY
      NADD=ISTATE(5)           ! NUMBER OF ADDITIONAL CROSS SECTIONS
      NFISS=ISTATE(4)          ! NUMBER OF FISSILE ISOTOPES
      ITR=ISTATE(6)            ! TRANSPORT CORRECTION OTPION
      NED=ISTATE(13)           ! NUMBER OF P0 ADDITIONAL XS
      ISPH=ISTATE(14)

      IF(IPRINT > 0)  THEN
         WRITE(6,*)
         WRITE(6,*) "******     BRANCH CHARACTERISTICS     ******"
         WRITE(6,*) "BRANCH TYPE         :",BRANCH
         WRITE(6,*) "BRANCH INDEX        :",ITBRA
         WRITE(6,*) "STATE VARIABLE NAME :",STAVAR
         WRITE(6,*) "BRANCH STATE VALUES :",STATE
      ENDIF

      IF(DATSRC(3)==0.0) THEN
        CALL D2PRFL(  IPDAT, IPMIC , IPRINT,    NBU,   NGP,   NBMIX,
     >                NANI,   NVAR,  STAIDX,   LADF,  NADF,   NTYPE)
      ELSE IF(DATSRC(3) == 1.0) THEN
        ! CASE FOR FUEL CROSS SECTIONS
        CALL LCMSIX(IPSAP,' ',0)
        CALL XDISET(DIMSAP,50,0)
        IF (SIGNAT .EQ. 'L_SAPHYB') THEN
         CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)  ! recover DIMSAP info
         NREA=DIMSAP(4)       ! NUMBER OF REACTIONS
         NISO=DIMSAP(5)       ! NUMBER OF PARTICULARIZED ISOTOPES
         NMAC=DIMSAP(6)       ! NUMBER OF MACROSCOPIC SETS
         NMIL=DIMSAP(7)       ! NUMBER OF MIXTURES
         NPAR=DIMSAP(8)       ! NUMBER OF STATE VARIABLE IN SAPHYB
         NADRX=DIMSAP(18)     ! CONCERN CROSS SECTIONS
                        ! (INCLUDING FLUE AND TIME)
        ELSE
         CALL LCMSIX(IPSAP,' ',0)
         CALL LCMSIX(IPSAP,MIXDIR,1)
         CALL LCMGET(IPSAP,'STATE-VECTOR',DIMSAP)
         NMIL = DIMSAP(1)
        ENDIF
        IF(STAVEC(1).NE.ISTATE(1)) THEN
          CALL XABORT("@D2PBRA: INCOHERENT NUMBER OF ENERGY GROUPS ")
        ENDIF


        IF(NMIL.NE.NBMIX) THEN
         CALL XABORT("@D2PBRA: DIFFERENT NUMBER OF MIX ")
         ENDIF

        ! RECOVER MACROLIB CROSS SECTIONS FROM SAPHYB
        CALL D2PMAC(     IPDAT, IPMIC , IPRINT,    NBU,   NGP,   NBMIX,
     >                    NADD,   NANI,  NVAR,  STAIDX,  LADF,    NADF,
     >                    NTYPE,  LCDF,  NCDF,    LGFF,    NGFF,  NPIN,
     >                     FLUX                                       )

        IF(LTH) THEN
          ICOR=STAVEC(22)
          ! RECOVER THE T/H INVARIANT BLOCK (OPTIONAL IN PMAXS FILES)
          CALL D2PTH(    IPDAT, IPMIC , IPRINT,    NBU,   NGP,   NBMIX,
     >                   NFISS,   NDEL,   NVAR, STAIDX,JOBOPT,    FLAG)
        ENDIF

        IF((LXES).OR.(LDET).OR.(LCOR)) THEN
          ! RECOVER MICROSCOPIC CROSS SECTIONS FROM SAPHYB
          CALL D2PMIC  (  IPDAT, IPMIC , IPRINT,    NGP,  NBMIX, NBISO,
     >                      NED,   NVAR, STAIDX,   LXES,   LDET,  LCOR,
     >                     FLUX                                       )
        ENDIF

        IF((GRID<2).and. (SIGNAT .EQ. 'L_SAPHYB')) THEN

          ! RECOVER THE DIVERS DIRECTORY OF SAPHYB
          CALL D2PDIV(    IPDAT, IPSAP , IPRINT,    NGP,   NBU,   NVAR,
     >                     GRID,   NPAR,   NREA,   NISO,  NMAC,   NMIL,
     >                     NANI,  NADRX, STAIDX,  STATE, STAVAR,   NSF,
     >                     LABS,   SCAT,   LADF                       )
        ENDIF


      ENDIF

      IF(REW.EQ.NBU) THEN
        ! REINITIALIZATION OF INDEX
        IF (FLAG.EQ.-1) THEN
         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
         CALL LCMGET(IPDAT,'REF_STATE',REFSTA)
         STATE(1:NVAR-1)=REFSTA(:)
         FLAG=0
         CALL LCMPUT(IPDAT,'FLAG',1,1,FLAG)
        ENDIF
         STAIDX(NVAR)= 1
         REW = 1
         STATE(NVAR)=BURN(1)


      ELSE
        !  UPDATE THE INDEX FOR THE CALCULATION OF THE NEXT BRANCH
        REW=0
        STAIDX(NVAR)= STAIDX(NVAR)+1
        REW = STAIDX(NVAR)
        STATE(NVAR)=BURN(STAIDX(NVAR))
      ENDIF

      ! STORE NEW VALUES OF BRANCH CALCULATION
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMPUT(IPDAT,'REWIND',1,1,REW)
      CALL LCMPUT(IPDAT,'STATE',NVAR,2,STATE)
      CALL LCMPUT(IPDAT,'STATE_INDEX',NVAR,1,STAIDX)
      END
