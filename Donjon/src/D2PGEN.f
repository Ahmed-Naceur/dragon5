*DECK D2PGEN
      SUBROUTINE D2PGEN( IPINP, IPDAT,   STAVEC, JOBTIT, FILNAM,    DER,
     >                    VERS,   COM,   JOBOPT,   IUPS,   FA_K,   SFAC,
     >                    BFAC,   DEB,     XESM,  FC1  ,    FC2,    FC3,
     >                     FC4,   XSC,   IPRINT                        )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create the GENPMAXS input file GENPMAXS.inp at phase 1
* WARNING: 04/2014: the format of this file respects the GENPMAXS format
* (it can't be changed)
* The information is recovered from the input file (.x2m) and stored in
* the INFO DATA block. The user can change any values in the input file
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPINP   file unit of GENPMAXS.inp file
* IPDAT   address of info data block
* VERS    version of PARCS to be used
* SFAC    the scattering cross section factor
* BFAC    the multiplier for betas
* DEB     FLAG to indicate the first call to the D2PGEN subroutine
* FA_K    assembly type
*         =0 reflector
*         =1 assembly
* IUPS    treatment for upscattering
*         =0 keep up scatter XS
*         =1 remove up scatter XS, modify  down scatter with DRAGON
*         spectrum (not available in this version)
*         =2 remove up scatter XS, modify  down scatter with infinite
*         medium spectrum
* STAVEC  various parameters associated with the IPDAT structure
* FILNAM  name of IPINP
* JOBTIT  title of in header of PMAXS file
* COM     comment to be printed in PMAXS file
* DER     partials derivative (T) or row cross section (F) to be stored
*         in PMAXS
* JOBOPT  array of flag to indicate the content option in the HELIOS
*         like file and PMAXS
* XESM    option for comparing k-inf in GenPMAX (1: using Pm/Sm data; 
*         2: using I/Xe data; 3: using I/Xe/Pm/Sm data)
*
*Parameters: 
* FC1     
* FC2     
* FC3     
* FC4     
* XSC     
* IPRINT  
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT
      INTEGER IPINP,STAVEC(40),FA_K,IUPS,DEB,XESM
      REAL SFAC,BFAC,VERS
      CHARACTER FILNAM*12,COM*40
      CHARACTER*16 JOBTIT
      CHARACTER*1 DER
      REAL FC1(5)
      REAL FC2(8)
      REAL FC3(7)
      REAL FC4(3)
      REAL XSC(3)

*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPTH,KPTH
      ! INDEX AND FLAG EXISTENCE OF : TEMPERATURE OF FUEL AND MODERATO
      INTEGER ITCOM, ITMOD,TMOD, TCOM
      ! NUMBER OF STATES VARIABLES, BURNUP EXEPTED
      INTEGER  NVAR
      ! NUMBER OF STATES VARIABLES
      INTEGER  STVARN
      ! IF FLAG=1, END OF A BRANCH CALCULATION
      INTEGER FLAG_PRINT, FLAG
      ! NUMBER OF BRANCH CONTAINED IN THE FINAL PMAXS FILE
      INTEGER NBR
      ! INDEX OF THE CURRENT BRANCH AND NUMBER OF BURNUP POINTS
      INTEGER BR_IT, NBU
      INTEGER CRDINF(STAVEC(6)),STAIDX(STAVEC(2))
      INTEGER IB,PK,GRID,NCRD,NDEL,NLOC,ST
      ! DATA SOURCE INFORMATION (CF GENPMAXS MANUAL)
      REAL DATSRC(5),LOCYLD(5),THCK
      REAL STATE(STAVEC(2)),HISTORY(STAVEC(2)-1), BU(STAVEC(4))
      CHARACTER(len=12) STATE_VAR(STAVEC(2))
      CHARACTER(len=4) STAVAR(STAVEC(2))
      INTEGER PKIDX(STAVEC(2))
      CHARACTER*1 JOBOPT(16)
      CHARACTER*4 BR
      CHARACTER*3 ADF_T
      CHARACTER*12,DIMENSION(6) :: PKNAM
      DATA PKNAM/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/
      LOGICAL LFLAG(6)
      LOGICAL :: LYLD=.FALSE.
      CHARACTER*3 YLDOPT

      ! INITIALIZATION OF ARRAYS
      DATSRC(1)= 2.0
      DATSRC(2)= 1.0
      DATSRC(3)= FA_K
      DATSRC(4)= SFAC
      DATSRC(5)= BFAC
      NGP=STAVEC(1)
      NBU=STAVEC(4)
      STVARN=STAVEC(2)
      NVAR=STVARN-1
      NCRD=STAVEC(6)
      NDEL=STAVEC(7)

      !RECOVER INFORMATION FROM INFO DATA block
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      IF (JOBOPT(9).EQ. 'T') LYLD = .TRUE.
      IF (LYLD) THEN
      CALL LCMGTC(IPDAT,'YLD_OPT',3,1,YLDOPT)
      CALL LCMGET(IPDAT,'YLD_LOC',LOCYLD)
      ENDIF
      CALL LCMGTC(IPDAT,'STATE_VAR',12,STVARN,STATE_VAR)
      CALL LCMGET(IPDAT,'PKIDX',PKIDX)

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

      IF(FA_K==1) THEN
        GRID=STAVEC(5)
      ELSE
        GRID = 2  ! NO XE/SM FOR REFLECTOR CASE
      ENDIF
      ITCOM=0
      ITMOD=0
      TCOM=0
      TMOD=0

      DO IST=1, STVARN
        IF(STATE_VAR(IST)==PKNAM(1)) STAVAR(IST)='CR  '
        IF(STATE_VAR(IST)==PKNAM(2)) STAVAR(IST)='DC  '
        IF(STATE_VAR(IST)==PKNAM(3)) STAVAR(IST)='PC  '
        IF(STATE_VAR(IST)==PKNAM(6)) STAVAR(IST)='BU  '
        IF(STATE_VAR(IST)==PKNAM(4)) THEN
          ITCOM=IST
          TCOM=1
          STAVAR(IST)='TF  '
        ENDIF
        IF(STATE_VAR(IST)==PKNAM(5)) THEN
          ITMOD=IST
          TMOD = 1
          STAVAR(IST)='TC  '
        ENDIF
      ENDDO

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMPTC(IPDAT,'IDEVAR',4,STVARN,STAVAR)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)

      IF(DEB.LE.0) THEN

          CALL LCMSIX(IPDAT,' ',0)
          CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
        CALL LCMPUT(IPDAT,'FLAG',1,1,DEB)
        FLAG=DEB
        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMGET(IPDAT,'BARR_INFO',CRDINF)
        CALL D2PREF( IPDAT, STVARN,   CRDINF,   NCRD,  GRID, PKIDX,
     1               PKNAM, IPRINT                                )
        IF (LYLD.and. (YLDOPT.EQ.'MAN')) THEN
         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
         CALL LCMGET(IPDAT,'STATE',STATE)
         NLOC=0
         ST=0
         DO I=1,5
          IF (LOCYLD(I).NE. -1) THEN
           NLOC=NLOC+1
           IF (LFLAG(I)) THEN
            ST=ST+1
            STATE(I)=LOCYLD(I)
            STATE_VAR(I)=PKNAM(I)
           ELSE IF (I.EQ.1) THEN
            ST=ST+1
            STATE_VAR(I)=PKNAM(I)
           ENDIF
          ENDIF
         ENDDO
         IF ((NLOC.NE.ST).OR.(NLOC.NE.NVAR)) THEN
          WRITE(6,*) '@D2PGEN : INCORRECT NUMBER OF STATE PARAMETERS',
     >    ' SET IN "YLD MAN" CARD : ',NLOC
          CALL XABORT('=> PLEASE FOLLOW THE SAP/MCO OBJECT CONTENT.')
         ENDIF

         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
         CALL LCMPUT(IPDAT,'STATE',STVARN,2,STATE)
        ENDIF

      ELSE
        CALL LCMGET(IPDAT,'FLAG',FLAG)
      ENDIF

      IF(FLAG .LE. 0) THEN
        !FIRST CALL TO D2PGEN SUBROUTINE
        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
        ! CHECK THE JOB_OPT VECTOR

        ! LDED : DIRECT ENERGY DEPOSITION FRACTION NOT IMPLEMENTED
        IF(JOBOPT(3)=='T') JOBOPT(3)='F'
        ! LJ1F : J1 FACTOR FOR MINIMAL CRITICAL POWER RATIO
        IF(JOBOPT(6)=='T') JOBOPT(6)='F'
        ! LCHD : DELAY NEUTRON FISSION SPECTRUM NOT IMPLEMENTED
        ! IF(JOBOPT(8)=='T') JOBOPT(8)='F'
        ! LBET : BETA NOT IMPLEMENTED
        IF((JOBOPT(12)=='T') .and. NDEL > 6) THEN
          JOBOPT(12)='F'
          WRITE(6,*) "@D2PGEN: WARNING "
          WRITE(6,*) "NUMBER OF DELAYED NEUTRON GROUPS > 6 "
          ! HELIOS FORMAT ACCEPTS ONLY NDEL =6
          WRITE(6,*) "lbet (JOBOPT(12)) FLAG FORCED TO FALSE "
        ENDIF
        ! LDEC : DECAY HEAT DATA NOT IMPLEMENTED
        IF(JOBOPT(14)=='T') JOBOPT(14)='F'
        IF((JOBOPT(13)=='T') .and. NDEL > 6) THEN
          JOBOPT(13)='F'
          WRITE(6,*) "@D2PGEN: WARNING "
          WRITE(6,*) "NUMBER OF DELAYED NEUTRON GROUPS > 6 "
          ! HELIOS FORMAT ACCEPTS ONLY NDEL =6
          WRITE(6,*) "lamb (JOBOPT(13)) FLAG FORCED TO FALSE "
        ENDIF

        ! RECOVER INFORMATION FROM GENPMAXS_INP
        CALL LCMPTC(IPDAT,'JOB_TIT',16,1,JOBTIT)
        CALL LCMPTC(IPDAT,'DERIVATIVE',1,1,DER)
        CALL LCMPTC(IPDAT,'JOB_OPT',1,16,JOBOPT(:16))
        CALL LCMPUT(IPDAT,'IUPS',1,1,IUPS)
        CALL LCMPUT(IPDAT,'XESMOPT',1,1,XESM)
        CALL LCMPUT(IPDAT,'DAT_SRC',5,2,DATSRC)
        CALL LCMPTC(IPDAT,'COMMENT',40,1,COM)
        CALL LCMPUT(IPDAT,'VERSION',1,2,VERS)
        CALL LCMPTC(IPDAT,'FILE_NAME',12,1,FILNAM)

        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'HELIOS_HEAD',1)
        CALL LCMPUT(IPDAT,'FILE_CONT_1',2,2,FC1(4:5))
        CALL LCMPUT(IPDAT,'FILE_CONT_2',8,2,FC2)
        CALL LCMPUT(IPDAT,'FILE_CONT_3',7,2,FC3)
        CALL LCMPUT(IPDAT,'FILE_CONT_4',3,2,FC4)
        CALL LCMPUT(IPDAT,'XS_CONT',3,2,XSC)
        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
        ! RECOVER HISTORY STATE AND number of branches
        CALL LCMGET(IPDAT,'HST_STATE',HISTORY)
        CALL LCMGET(IPDAT,'BRANCH_NB',NBR)

        ! WRITING JOBTIT CARD

        IF(IPRINT > 2)  THEN
          WRITE(6,*)
          WRITE(6,*) "******* INFORMATION FOR GENPMAXS INPUT *********"
          WRITE(6,*)
          WRITE(6,*) "JOB_TIT CARD : JOB_TIT,DERIVATIVE,",
     1    " VERSION, COMMENT"
          WRITE(6,*) "VALUES    :",JOBTIT,DER, VERS, COM
          WRITE(6,*)
          WRITE(6,*) "JOB_OPT CARD : ad,xe,de,j1,ch,Xd,iv,dt,yl,cd,gf,",
     1    " be,lb,dc,ups"
          WRITE(6,'(A,14(A,1X))') "VALUES    :",JOBOPT(1:14)
          WRITE(6,*)
          WRITE(6,*) "DAT_SRC CARD : SRC_KIND, NFILE, FA_KIND, SFAC,",
     1    " BFAC"
          WRITE(6,*) "VALUES    :",INT(DATSRC)
          WRITE(6,*)
          WRITE(6,*) "STAVAR CARD :"
          WRITE(6,*) "NUMBER    :",STVARN
          WRITE(6,*) "VALUES    :",STAVAR(1:STVARN)
          WRITE(6,*)
          WRITE(6,*) "HISTORY CARD (IN GENPMAXS FORMALISM):"
          WRITE(6,*) "VALUES OF STATES VARIABLES :",HISTORY(1:NVAR)
          WRITE(6,*)
          WRITE(6,*) "BRANCH CARD :"
          WRITE(6,*) "NUMBER OF BRANCHES : ",NBR
          WRITE(6,*)

        ENDIF
      ELSE IF(FLAG == 1) THEN
        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMGET(IPDAT,'BARR_INFO',CRDINF)
        CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
        CALL LCMGET(IPDAT,'STATE',STATE)
        CALL LCMGET(IPDAT,'STATE_INDEX',STAIDX)
        CALL LCMGET(IPDAT,'PRINT',FLAG_PRINT)
        CALL LCMGTC(IPDAT,'BRANCH',4,1,BR)
        CALL LCMGET(IPDAT,'BRANCH_IT',BR_IT)

        ! REORG  BARR INFORMATION

        DO IB=1 ,NCRD
          IF(STATE(1)==CRDINF(IB)) THEN
           STATE(1)=IB-1
           EXIT
          ENDIF
        ENDDO
        ! TEMPERATURE CONVERSION
        IF (STAVEC(19).EQ.0) THEN
         IF(TCOM==1) STATE(ITCOM)=STATE(ITCOM)+273.15 ! convert C to K
         IF(TMOD==1) STATE(ITMOD)=STATE(ITMOD)+273.15
        ENDIF
        ! CONTINUE WRITING BRANCH CARD
        IF(FLAG_PRINT==1) THEN
          WRITE (IPINP,'(A,A,I4.4,3X,3(F11.5,1X,F11.5,1X))')
     1    'HIST',BR(1:2),BR_IT,(STATE(I), I=1,NVAR)
        ENDIF

        IF(IPRINT > 2)  THEN
          WRITE(6,*)
          WRITE(6,*) "*CONTINUE WRITING BRANCH CARD IN GENPMAXS INPUT*"
          WRITE(6,*)
          WRITE(6,*) "BRANCH TYPE          : ",BR
          WRITE(6,*) "BRANCH INDEX         : ",BR_IT
          WRITE(6,*) "BRANCH STATE VALUES  : ",STATE(1:NVAR)
          WRITE(6,*) "BRANCH STATE INDEX   : ",STAIDX(1:NVAR)
          WRITE(6,*)
        ENDIF
      ELSE IF(FLAG == 2) THEN
        ! RECOVER INFORMATION FROM INFO
        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)

        CALL LCMGET(IPDAT,'BURN',BU)
        IF(JOBOPT(1)=='T') THEN
         CALL LCMGTC(IPDAT,'ADF_TYPE',3,1,ADF_T)
         IF (ADF_T.EQ.'GEN') CALL LCMGET(IPDAT,'THCK',THCK)
        ENDIF
        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
        CALL LCMGET(IPDAT,'BRANCH_NB',NBR)


        CALL LCMSIX(IPDAT,' ',0)
        CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
        CALL LCMGTC(IPDAT,'FILE_NAME',12,1,FILNAM)
        CALL LCMGET(IPDAT,'DAT_SRC',DATSRC)


        ! WRITING BURNUP CARD
        WRITE (IPINP,'(A)') '%BURNUP'
        WRITE (IPINP,'(I1,/,A,3X,I2)') 1, 'set1',NBU
        WRITE (IPINP,'(5(3X,F8.3),/)') (BU(I)/1000, I=1,NBU)
        WRITE (IPINP,'(A,1X,I4,A,I1)') 'HIST01',NBR,'*',1

        IF ((INT(DATSRC(3)).EQ.0).AND.(JOBOPT(1)=='T')
     >  .AND.(ADF_T.EQ.'GEN'))THEN
         WRITE (IPINP,'(A)')'%ADF_1D'
         WRITE (IPINP,'(A)')'ANM   0   1'
         WRITE (IPINP,'(F8.5)')THCK
        ENDIF
        ! WRITING HEL_FMT CARD
        WRITE (IPINP,'(A)') '%HEL_FMT'
        WRITE (IPINP,'(I1,/,I1,1X,I2,1X,I2,1X,I1)') 1,1,24,12,8

        ! WRITING FIL_CNT CARD
        WRITE (IPINP,'(A)') '%FIL_CNT'
        WRITE (IPINP,'(I1,3X,A,3X,I4,3X,I1)') 1,FILNAM,NBR,1
        DO IB=1,NBR
          WRITE (IPINP,'(I4,1X,I1,1X,I4,1X,I1,1X,I2)') IB,1,IB,1,NBU
        ENDDO

        ! WRITING JOB_END CARD
        WRITE (IPINP,'(A)') '%JOB_END'

        IF(IPRINT > 2)  THEN
          WRITE(6,*)
          WRITE(6,*) "***** END OF EDITING THE GENPMAXS INPUT ******"
          WRITE(6,*)
          WRITE(6,*) "BURNUP CARD : "
          WRITE(6,*) "VALUES OF BURNUP POINTS :",BU/1000
          WRITE(6,*)
          WRITE(6,*) "HEL_FMT CARD : NFMT, Index, LABEL, WIDTH, COLUMN"
          WRITE(6,*) "VALUES (FIXED)   :",1,4,24,12,8
          WRITE(6,*)
          WRITE(6,*) "EDIT FIL_CNT CARD "
          WRITE(6,*)
        ENDIF
      ENDIF
      CALL LCMSIX(IPDAT,' ',0)

      END
