*DECK D2PTH
      SUBROUTINE D2PTH(  IPDAT, IPMIC , IPRINT,    NBU,   NGP,    NBMIX,
     >                   NFISS,   NDEL,   NVAR, STAIDX,JOBOPT,     FLAG)

*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover T/H inveariant data block and store in INFO/TH_DATA/
* WARNING: These data are extracted only if the corresponding flag is
* set to T in the GENPMAXS_INP/JOBOPT vector.
* NB 1 : The data for T/H are recovered from the reference state,  the
* branching calculation not includes the TH informations.
* NB 3 : The Helios format cannot recover the CHID (delay neutron
* fission spectrum), it is fixed to default values even if JOBOPT(6)=T.
* NB 4 : The Helios format cannot recover the Decay Heat Data (DBET and
* DLAM in GenPMAXS), it is fixed to default values even if JOBOPT(14)=T.
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of the INFO data block
* IPMIC   address of the MICROLIB object
* IPPRINT control the printing on screen
* NGP     number of energy groups
* NBU     number of burnup point in IPSAP
* NVAR    number of state parameters in INFO data block
* NDEL    number of delaued neutron groups
* NBMIX   number of mixtrures  in IPSAP
* NFISS   number of fissile isotopes
* STAIDX  index of state variables
* FLAG    End of a bran calculation (=-1: branch for yields calculation)
*
*Parameters: 
* IPRINT  
* JOBOPT  
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPMIC
      INTEGER IPRINT,NVAR,NBU, NBMIX,NGP
      INTEGER NFISS,NDEL
      INTEGER STAIDX(NVAR)
      INTEGER FLAG
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMIC,KPMIC,IPTH,KPTH
      PARAMETER(NSTATE=40)
      INTEGER DSTATE(NSTATE)
      INTEGER NDFI,NDFP,MR,MI,MI_REAL,ITYLCM
      INTEGER :: I_PF = 0
      INTEGER :: iso = 1
      REAL YLDI,YLDXe,YLDPm
      REAL CHI(NFISS,NGP)
      REAL OVERV(NGP),BURN(NBU), STATE(NVAR)
      REAL FLX(NGP),NUSIGF_D(NDEL,NGP),NUSIGF(NGP)
      REAL BETA_D(NDEL,NFISS),LAMBDA_D(NDEL,NFISS)
      REAL NUM(NDEL)
      REAL :: DEN = 0.0
      CHARACTER*12 ISOTOPES(4)
      CHARACTER*1 JOBOPT(16)
      CHARACTER*8 NUSID
      CHARACTER*3 YLDOPT
      REAL YLDFIX(3)


      REAL, ALLOCATABLE, DIMENSION(:) :: DEPLETE_ENER,DEPLETE_DECA
      REAL, ALLOCATABLE, DIMENSION(:) :: FISSIONYIELD
      CHARACTER(len=12),ALLOCATABLE, DIMENSION(:) :: ISOTOPERNAME
      CHARACTER(len=12),ALLOCATABLE, DIMENSION(:) :: ISOTOPESDEPL,PF

      IF(IPRINT > 1)  THEN
         WRITE(6,*)
         WRITE(6,*) "**************************************************"
         WRITE(6,*) "*               T/H INVARIANT BLOCK              *"
         WRITE(6,*) "**************************************************"
         WRITE(6,*)
      ENDIF

      CALL LCMSIX(IPMIC,' ',0)
      CALL LCMSIX(IPMIC,'MACROLIB',1)

      IF(JOBOPT(13)=='T') CALL LCMGET(IPMIC,'LAMBDA-D',LAMBDA_D)



      JPMIC=LCMGID(IPMIC,'GROUP')

      IF(NBMIX.NE.1) THEN
        CALL XABORT('@D2PTH: MORE THAN ONE MIXTURE IN SAPHYB')
      ENDIF
      IF(NFISS.NE.1) THEN
        CALL XABORT('@D2PTH: MORE THAN 1 FISSILE ISOTOPE IN MACROLIB')
      ENDIF

      DO IGR=1,NGP
        KPMIC=LCMGIL(JPMIC,IGR)
        IF(JOBOPT(7)=='T')CALL LCMGET(KPMIC,'OVERV',OVERV(IGR))
        IF(JOBOPT(5)=='T')CALL LCMGET(KPMIC,'CHI',CHI(1:NFISS,IGR))
        IF(JOBOPT(12)=='T') THEN
         CALL LCMGET(KPMIC,'NUSIGF',NUSIGF(IGR))
         CALL LCMGET(KPMIC,'FLUX-INTG',FLX(IGR))
         DO ND=1,NDEL
           WRITE(NUSID,' (A6, I2.2)') 'NUSIGF', ND
           CALL LCMGET(KPMIC,NUSID,NUSIGF_D(ND,IGR))
         ENDDO
        ENDIF
      ENDDO

      IF(JOBOPT(12)=='T') THEN
        DO ND=1,NDEL

           DEN=0.
           NUM(ND)=0.0
          DO IGR= 1,NGP
           DEN=DEN+NUSIGF(IGR)*FLX(IGR)
           NUM(ND)=NUM(ND)+NUSIGF_D(ND,IGR)*FLX(IGR)
          ENDDO
           BETA_D(ND,NFISS)=NUM(ND)/DEN
         ENDDO
!        CALL XABORT ('STOP TEST')
      ENDIF
      IF(JOBOPT(9)=='T') THEN

        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
        CALL LCMGTC(IPDAT,'YLD_OPT',3,1,YLDOPT)

        CALL LCMGET(IPDAT,'YLD_FIX',YLDFIX)

        IF ((YLDOPT=='REF').OR.(YLDOPT=='MAN')) THEN
         CALL LCMSIX(IPMIC,' ',0)
         CALL LCMLEN(IPMIC,'ISOTOPESDENS',MI_REAL,ITYLCM)
         CALL LCMLEN(IPMIC,'ISOTOPERNAME',MI,ITYLCM)
         ALLOCATE (ISOTOPERNAME(MI))
         CALL LCMGTC(IPMIC,'ISOTOPERNAME',12,MI,ISOTOPERNAME)
         CALL LCMLEN(IPMIC,'DEPL-CHAIN',ILONG,ITYLCM)
         IF (ILONG.EQ.0) THEN
          YLDI=YLDFIX(1)
          YLDXe=YLDFIX(2)
          YLDPm=YLDFIX(3)
           WRITE(6,*)"@D2PTH : NO RECORD DEPL-CHAIN IN SAP/MCO :"
           WRITE(6,*)"=> DEFAULT VALUES FOR FISSION YLDS CONSIDERED"
         ELSE
          CALL LCMSIX(IPMIC,'DEPL-CHAIN',1)
          CALL LCMGET(IPMIC,'STATE-VECTOR',DSTATE)

          NDEPL = DSTATE(1)
          NDFI = DSTATE(2)
          NDFP = DSTATE(3)
          MR = DSTATE(8)

          ALLOCATE (FISSIONYIELD(NDFI*NDFP), DEPLETE_ENER(NDEPL*MR))
          ALLOCATE (ISOTOPESDEPL(NDEPL), PF(NDEPL),DEPLETE_DECA(NDEPL))
          CALL LCMGET(IPMIC,'DEPLETE-DECA',DEPLETE_DECA)
          CALL LCMGET(IPMIC,'DEPLETE-ENER',DEPLETE_ENER)
          CALL LCMGTC(IPMIC,'ISOTOPESDEPL',12,NDEPL,ISOTOPESDEPL)

          IF ((NDFI.EQ.0 ).OR. (NDFP .EQ. 0)) THEN
           WRITE(6,*) "@D2PTH : NUMBER OF DIRECT FISSILE ISOTOPES",
     1     " OR FISSION FRAGMENT IS ZERO"
           CALL XABORT('=> PLEASE TURN OFF THE LYLD FLAG IN JOB_OPT'
     >    //' OR USE THE "YLD FIX" OPTION' )
          ENDIF
          CALL LCMGET(IPMIC,'FISSIONYIELD',FISSIONYIELD)


          I_PF=0
          iso=1

          CALL LCMSIX(IPDAT,' ',0)
          CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
          CALL LCMGET(IPDAT,'BURN',BURN)
          CALL LCMSIX (IPDAT,'ISOTOPES',1)
          CALL LCMGTC (IPDAT,'XE135',12,1,ISOTOPES(1))
          CALL LCMGTC (IPDAT,'SM149',12,1,ISOTOPES(2))
          CALL LCMGTC (IPDAT,'I135',12,1,ISOTOPES(3))
          CALL LCMGTC (IPDAT,'PM149',12,1,ISOTOPES(4))

          DO iso=1, NDEPL
            IF(INDEX(ISOTOPESDEPL(iso), 'MACR')==0) THEN

             I_PF=I_PF+1
             PF(I_PF)=ISOTOPESDEPL(iso)
             IF(PF(I_PF)==ISOTOPES(3)) YLDI=FISSIONYIELD(I_PF)
             IF(PF(I_PF)==ISOTOPES(1)) YLDXe=FISSIONYIELD(I_PF)
             IF(PF(I_PF)==ISOTOPES(4)) YLDPm=FISSIONYIELD(I_PF)
            ENDIF
          ENDDO

          IF(IPRINT > 1)  THEN
           WRITE(6,*)"********* STATE VECTOR INFORMATION *************"
           WRITE(6,*)
           WRITE(6,*)"Number of isotopes (MI)                  : ",MI
           WRITE(6,*)"Number of groups (NGP)                   : ",NGP
           WRITE(6,*)"Number of fissile isotopes (NFISS)       : ",NFISS
           WRITE(6,*)"Number of delayed neutron groups (NDEL)  : ",NDEL
           WRITE(6,*)"Number of depleted isotopes (NDEPL)      : ",NDEPL
           WRITE(6,*)"Number of direct fissile isotopes (NDFI) : ",NDFI
           WRITE(6,*)"Number of fission fragments (NDFP)       : ",NDFP
           WRITE(6,*)"Maximum number of depleting reactions(MR): ",MR
           WRITE(6,*)
           WRITE(6,*)"**************** ISOTOPE NAME ******************"
           WRITE(6,*)
           WRITE(6,'(10A12)')ISOTOPERNAME(1:MI_REAL)
           WRITE(6,*)
          ENDIF
          DEALLOCATE (ISOTOPERNAME)
          DEALLOCATE (FISSIONYIELD,ISOTOPESDEPL,PF)
          DEALLOCATE (DEPLETE_ENER,DEPLETE_DECA)
         ENDIF
        ELSE IF (YLDOPT=='FIX') THEN
         YLDI=YLDFIX(1)
         YLDXe=YLDFIX(2)
         YLDPm=YLDFIX(3)
        ENDIF
      ENDIF

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMGET(IPDAT,'STATE',STATE)
      CALL LCMSIX(IPDAT,' ',0)

      IF(STAIDX(NVAR)==1) THEN
        IPTH=LCMLID(IPDAT,'TH_DATA',NBU)
      ELSE
        IPTH=LCMGID(IPDAT,'TH_DATA')
      ENDIF

      KPTH=LCMDIL(IPTH,STAIDX(NVAR))

      IF(JOBOPT(13)=='T') THEN
        CALL LCMPUT(KPTH,'LAMBDA',NDEL*NFISS,2,LAMBDA_D)
      ENDIF
      IF(JOBOPT(9)=='T') THEN
       IF((YLDOPT.EQ.'FIX').OR.(YLDOPT.EQ.'REF')) THEN
        CALL LCMPUT(KPTH,'YLDPm',1,2,YLDPm)
        CALL LCMPUT(KPTH,'YLDXe',1,2,YLDXe)
        CALL LCMPUT(KPTH,'YLDI',1,2,YLDI)
       ELSE IF ((YLDOPT.EQ.'MAN').AND.(FLAG.EQ.-1)) THEN
        CALL LCMPUT(KPTH,'YLDPm',1,2,YLDPm)
        CALL LCMPUT(KPTH,'YLDXe',1,2,YLDXe)
        CALL LCMPUT(KPTH,'YLDI',1,2,YLDI)
       ENDIF
      ENDIF

      IF(JOBOPT(7)=='T')CALL LCMPUT(KPTH,'OVERV',NGP,2,OVERV)
      IF(JOBOPT(5)=='T')CALL LCMPUT(KPTH,'CHI',NFISS*NGP,2,CHI)
      IF(JOBOPT(12)=='T')CALL LCMPUT(KPTH,'BETA',NDEL*NFISS,2,BETA_D)
      IF(IPRINT > 1)  THEN
        WRITE(6,*) "**************** T/H INFORMATION *****************"
        IF(JOBOPT(5)=='T') WRITE(6,*) "CHI(NFISS,NGP)     :",CHI
        IF(JOBOPT(7)=='T') WRITE(6,*) "OVERV(NGP)         :",OVERV
        IF(JOBOPT(13)=='T')WRITE(6,*) "LAMBDA(NDEL,NFISS) :",LAMBDA_D
        IF(JOBOPT(12)=='T')WRITE(6,*) "BETA(NDEL,NFISS)   :",BETA_D
        IF(JOBOPT(9)=='T') WRITE(6,*) "PM-149 YIELD       :",YLDPm
        IF(JOBOPT(9)=='T') WRITE(6,*) "XE-135 YIELD       :",YLDXe
        IF(JOBOPT(9)=='T') WRITE(6,*) "I-135 YIELD        :",YLDI
        WRITE(6,*)
      ENDIF

      END
