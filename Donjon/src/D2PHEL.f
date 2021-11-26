      SUBROUTINE D2PHEL ( IPHEL,  IPDAT,  IPMIC ,  IPINP,     STAVEC,
     >                   JOBOPT, IPRINT                             )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store the header of HELIOS.dra file - (independant data compared with
* branching calculation) at phase 1
* WARNING: 04/2014 : the format of this file respect the HELIOS format
* (it cannot be changed)
* The information is recovered from the input file (.x2m) and stored in
* the INFO DATA block. The user can change any values in the input file
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPHEL   file unit of HELIOS like file
* IPDAT   adress of info data block
* STAVEC  various parameters associated with the IPDAT structure
* FC1     FILE_CONT_1 recovered from D2P: input
* FC2     FILE_CONT_2 recovered from D2P: input
* FC3     FILE_CONT_3 recovered from D2P: input
* FC4     FILE_CONT_4 recovered from D2P: input
* XSC     XS_CONT     recovered from D2P: input
* IPRINT  control the printing on screen
*
*Parameters: 
* IPMIC   
* IPINP   
* JOBOPT  
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPMIC
      INTEGER IPHEL
      INTEGER STAVEC(40)
      ! FILE_CONT DATA BLOC ( CF D2P: DOCUMENTATION)
      REAL FC1(2)
      REAL FC2(8)
      REAL FC3(7)
      REAL FC4(3)
      REAL XSC(3)
      REAL DATSRC(5)
*----
*  LOCAL VARIABLES
*----
      INTEGER NBU,FA_K
      CHARACTER*16 JOBTIT
      CHARACTER*12 FILNAM
      CHARACTER*1 DER
      CHARACTER*40 COM
      CHARACTER*1 JOBOPT(16)
      REAL HISTORY(STAVEC(2)-1)
      CHARACTER*4 STAVAR(STAVEC(2))
      INTEGER IUPS,XESM
      REAL VERS


      NBU=STAVEC(4)
      NPAR=STAVEC(2)
      NVAR=NPAR-1

      ! RECOVER INFORMATION FROM INFO/HELIOS_HEAD DATA BLOCK
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)

      CALL LCMGTC(IPDAT,'IDEVAR',4,NPAR,STAVAR)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMGET(IPDAT,'DAT_SRC',DATSRC)
      CALL LCMGTC(IPDAT,'JOB_TIT',16,1,JOBTIT)
      CALL LCMGTC(IPDAT,'DERIVATIVE',1,1,DER)
      CALL LCMGET(IPDAT,'IUPS',IUPS)
      CALL LCMGET(IPDAT,'XESMOPT',XESM)
      CALL LCMGTC(IPDAT,'COMMENT',40,1,COM)
      CALL LCMGET(IPDAT,'VERSION',VERS)
      CALL LCMGTC(IPDAT,'FILE_NAME',12,1,FILNAM)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'HELIOS_HEAD',1)
      CALL LCMGET(IPDAT,'FILE_CONT_1',FC1)
      CALL LCMGET(IPDAT,'FILE_CONT_2',FC2)
      CALL LCMGET(IPDAT,'FILE_CONT_3',FC3)
      CALL LCMGET(IPDAT,'FILE_CONT_4',FC4)
      CALL LCMGET(IPDAT,'XS_CONT',XSC)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      !RECOVER HISTORY STATE AND number of branches
      CALL LCMGET(IPDAT,'HST_STATE',HISTORY)
      CALL LCMGET(IPDAT,'BRANCH_NB',NBR)

      IF (IUPS.EQ.2) IUPS=0
      FA_K=INT(DATSRC(3))
      IF ((STAVEC(21).EQ.1) .and. (JOBOPT(1).EQ.'T') )THEN
       JOBOPT(1)='F'
      ENDIF
      IF (STAVEC(19).EQ.0) THEN
       DO I=1,NVAR
        IF (STAVAR(I).EQ.'TF  ') THEN

         HISTORY(I)=HISTORY(I)+273.15
        ENDIF
        IF (STAVAR(I).EQ.'TC  ') THEN
         HISTORY(I)=HISTORY(I)+273.15
        ENDIF
       ENDDO
      ENDIF

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMPUT(IPDAT,'HST_STATE',NVAR,2,HISTORY)
        ! WRITING JOBTIT CARD
      WRITE (IPINP,*) '%JOB_TIT'
      WRITE (IPINP,'(A,A,A,1X,A,1X,F3.1,1X,A,A,A)')
     1'"',JOBTIT,'"',DER, VERS, '"',COM,'"'

      ! WRITING JOB_OPT CARD
      WRITE (IPINP,*) '%JOB_OPT'
      WRITE (IPINP,'(14(A,1X),2(I1,1X))',advance="no")
     1 JOBOPT(1:14),IUPS,XESM
      WRITE (IPINP,'(/A)')
     1'!ad,xe,de,j1,ch,Xd,iv,dt,yl,cd,gf,be,lb,dc,ups'

      ! WRITING DAT_SRC CARD
      WRITE (IPINP,*) '%DAT_SRC'
      WRITE(IPINP,'(I2,1X,I2,1X,I2,1X,F3.1,1X,F3.1)')INT(DATSRC(1)),
     1INT(DATSRC(2)),INT(DATSRC(3)),DATSRC(4),DATSRC(5)

      ! WRITING STA_VAR CARD
      WRITE (IPINP,*) '%STA_VAR'
      WRITE (IPINP,'(I2/,3(A,1X,A))') NVAR,(STAVAR(I), I=1,NVAR)

      ! WRITING HISTORY CARD
      ! CONCERN THE CONTROL ROD COMPOSITION
      IF(HISTORY(1)==0) THEN
        HISTORY(1)=1
      ELSE IF(HISTORY(1)==1) THEN
        HISTORY(1)=0
      ELSE IF(HISTORY(1)==2) THEN
        HISTORY(1)=2
      ENDIF

      WRITE (IPINP,*) '%HISTORY'
      WRITE (IPINP,'(I1,1X,I1,/,A,1X,3(F11.5,1X,F11.5,1X))') 1,1,
     1'HIST01',(HISTORY(I), I=1,NVAR)

      ! WRITING BRANCH CARD
      WRITE (IPINP,*) '%BRANCH'
      WRITE (IPINP,'(I4,1X,I1)') NBR, 1


      ! WRITE FILE_CONT DATA in HELIOS.dra file
      IF(IPRINT > 0) WRITE(6,*) "STEP 1 : EDIT THE HEADER "
      CALL SET_INFO(IPHEL)
      IF(IPRINT > 0) WRITE(6,*) "STEP 2 : EDIT THE CONT1 BLOCK "
      IF (FA_K.EQ.0) THEN
       FC1(1)=0.
      ELSE
       IF (FC1(1).EQ.0.) THEN

        CALL LCMSIX(IPMIC,'',0)
        CALL LCMSIX(IPMIC,'MACROLIB',1)
        CALL LCMLEN(IPMIC,'MASL',ILONG,ITYLCM)
        IF (ILONG.GT.1) THEN
         CALL XABORT("@D2PHEL: MORE THAN 1 METAL DENS. IN THE MICROLIB")
        ELSE IF (ILONG.EQ.0) THEN
         WRITE(6,*)"@D2PHEL: RECORD MASL NOT FOUND IN MICROLIB"
         WRITE(6,*)"=> PLEASE USE THE FILE_CONT_1 CARD IN D2P:"
         CALL XABORT(" OR USE THE 'REFLECTOR' KEYWORD")
        ELSE
         CALL LCMGET(IPMIC,'MASL',FC1(1))
        ENDIF
       ELSE  IF (FC1(1).LE.0.) THEN
        CALL XABORT('@D2PHEL: NEGATIVE VALUE FOR HEAVY METAL DENSITY')
       ENDIF
      ENDIF
      CALL LCMPUT(IPDAT,'FILE_CONT_1',2,2,FC1)
      CALL SET_CONT1(IPHEL,STAVEC,FC1,IPRINT)
      ! IF(IPRINT > 0) WRITE(6,*) "STEP 3 : EDIT THE CONT2 BLOCK "
       ! CALL SET_CONT2(IPHEL,FC2,NGP,IPRINT)
      IF(IPRINT > 0) WRITE(6,*) "STEP 4 : EDIT THE CONT3 BLOCK "
      CALL SET_CONT3(IPHEL,FC3,IPRINT)
      IF(IPRINT > 0) WRITE(6,*) "STEP 5 : EDIT THE CONT4 BLOCK "
      CALL SET_CONT4(IPHEL,FC4,IPRINT)

      ! MOVE TO GENPMAXS_INP DIRECTORY
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      IF ((STAVEC(21).EQ.1) .and. (JOBOPT(1).EQ.'F') )THEN
       JOBOPT(1)='T'
      ENDIF
      END

      SUBROUTINE SET_CONT1(IPHEL,STAVEC,FILE_CONT_1,IPRINT)
      INTEGER STAVEC(40)
      REAL FILE_CONT_1(2)

      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s)  1) ==========================='
      WRITE(IPHEL,*) '               2) %FILE_CONT 1'
      WRITE(IPHEL,*) '               3) ==========================='
      WRITE(IPHEL,*) '               4) Meaning : NGROUP, NCOLS, NR'
     1 //'OWS, PART,'
      WRITE(IPHEL,*) '                  HM Density, Bypass Density '
      CALL SET_RIEGO(IPHEL)
      WRITE(IPHEL,120) 'NGROUP','NCOLS','NROWS','PART',
     1 'DenHM','DenByp'
      WRITE(IPHEL,125) 'Label E','.-.-E-.-.','.-.-E-.-.','.-.-E-.-.',
     1 '.-.-E-.-.','1-.-E-.-.','.-.-E-.-.'
      WRITE(IPHEL,130) '   1 HST  1 HST :      0',STAVEC(1),
     1 STAVEC(8),STAVEC(9), STAVEC(10),
     2 FILE_CONT_1(1),FILE_CONT_1(2)
      WRITE(IPHEL,'()')
      IF(IPRINT > 1)  THEN
        WRITE(6,*)
        WRITE(6,*) "CONTENT : NGROUP, NCOLS, NROWS, PART,",
     1  " HM Density, Bypass Density "
        WRITE(6,*) "VALUES  :",STAVEC(1),STAVEC(8:10),FILE_CONT_1
        WRITE(6,*)
      ENDIF
  120 FORMAT(27X,A,9X,A,9X,A,10X,A,9X,A,8X,A)
  125 FORMAT(6X,A,12X,A,3X,A,3X,A,3X,A,3X,A,3X,A)
  130 FORMAT(A,10X,I2,10X,I2,10X,I2,10X,I2,5X,F7.5,5X,F7.5)
      END

      SUBROUTINE  SET_CONT2(IPHEL,FILE_CONT_2,NGROUP,IPRINT)
      INTEGER NGROUP
      CHARACTER*9 LABEL
      REAL FILE_CONT_2(NGROUP)

      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s) 1) ==========================='
      WRITE(IPHEL,*) '              2) %FILE_CONT 2'
      WRITE(IPHEL,*) '              3) ==========================='
      WRITE(IPHEL,*) '              4)Meaning : Lower Energy of Neu'
     1 //'tron Groups'
      CALL SET_RIEGO(IPHEL)

      IF(NGROUP .EQ. 8) THEN
        WRITE(IPHEL,220) 'EMIN','EMIN'
        WRITE(IPHEL,225) 'Label E'
        DO I=1, NGROUP
         WRITE(LABEL,'(A,I1,A)')".-.-E-",I,"-."
         PRINT*,"LABEL",LABEL
         WRITE(IPHEL,'(A9,5X)',advance='no')LABEL
        ENDDO
        WRITE(IPHEL,230) '   1 HST  1 HST :      0',FILE_CONT_2(1),
     1 FILE_CONT_2(2)
      ELSE
        CALL XABORT ("@D2PHEL: NUMBER OF ENERGY GROUPS MUST BE 2")
      ENDIF

      IF(IPRINT > 1)  THEN
        WRITE(6,*)
        WRITE(6,*) "CONTENT : Lower Energy of Neutron Groups"
        WRITE(6,*) "VALUES  :",FILE_CONT_2 (1:NGROUP)
        WRITE(6,*)
      ENDIF

  220 FORMAT(32X,A,10X,A)
  225 FORMAT(6X,A,17X)
  230 FORMAT(A,ES12.5E2,ES12.5E2)
      END

      SUBROUTINE  SET_CONT3(IPHEL,FILE_CONT_3,IPRINT)
      REAL FILE_CONT_3(7)

      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s)  1) ==========================='
      WRITE(IPHEL,*) '               2) %FILE_CONT 3'
      WRITE(IPHEL,*) '               3) ==========================='
      WRITE(IPHEL,*) '               4)Meaning : Regions Volume'

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,320) 'VCool','VWatR','VModr','VCnRd','VFuel',
     1 'VClad','VChan'
      WRITE(IPHEL,310) 'Label E','1-.-E-.-.','.-.-E-.-.','1-.-E-.-.',
     1 '1-.-E-.-.','1-.-E-.-.','1-.-E-.-.','1-.-E-.-.'
      WRITE(IPHEL,390) '   1 HST  1 HST :      0',FILE_CONT_3(1),
     1  FILE_CONT_3(2),FILE_CONT_3(3),FILE_CONT_3(4),FILE_CONT_3(5),
     2  FILE_CONT_3(6),FILE_CONT_3(7)


      IF(IPRINT > 1)  THEN
        WRITE(6,*)
        WRITE(6,*) "CONTENT : VCool, VWatR, VModr, VCnRd, VFuel,",
     1  " VClad, VChan"
        WRITE(6,*) "VALUES  :",FILE_CONT_3
        WRITE(6,*)
      ENDIF

  310 FORMAT(6X,A,12X,A,3X,A,3X,A,3X,A,3X,A,3X,A,3X,A)
  320 FORMAT(27X,A,2X,A,9X,A,9X,A,9X,A,9X,A,9X,A,9X,A)
  390 FORMAT(A,ES12.5E2,ES12.5E2,ES12.5E2,ES12.5E2,
     1 ES12.5E2,ES12.5E2,ES12.5E2,ES12.5E2)
      END

      SUBROUTINE  SET_CONT4(IPHEL,FILE_CONT_4,IPRINT)
      REAL FILE_CONT_4(3)

      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) 'List name       : PMAX_FILE_DATA'
      WRITE(IPHEL,*) 'Labels Array    : KINF'
      WRITE(IPHEL,*) 'List Title(s)  1) ==========================='
      WRITE(IPHEL,*) '               2) %FILE_CONT 4'
      WRITE(IPHEL,*) '               3) ==========================='
      WRITE(IPHEL,*) '               4) Cell Pitch and X,Y Pos of F'
     1 //'irst Cell'

      CALL SET_RIEGO(IPHEL)

      WRITE(IPHEL,320) 'PITCH','XBE','YBE'
      WRITE(IPHEL,410) 'Label E','.-.-E-.-.','.-.-E-.-.','.-.-E-.-.'
      WRITE(IPHEL,390) '   1 HST  1 HST :      0',FILE_CONT_4(1),
     1  FILE_CONT_4(2),FILE_CONT_4(3)

      IF(IPRINT > 1)  THEN
        WRITE(6,*)
        WRITE(6,*) "CONTENT : PITCH ,XBE , YBE"
        WRITE(6,*) "VALUES  :", FILE_CONT_4
        WRITE(6,*)
      ENDIF

  320 FORMAT(24X,A,11X,A,11X,A)
  390 FORMAT(A,ES12.5E2,ES12.5E2,ES12.5E2)
  410 FORMAT(6X,A,12X,A,5X,A,5X,A)
      END

      SUBROUTINE SET_INFO(IPHEL)
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
      WRITE(IPHEL,*) 'Pre-processing for PMAXS Generation'
      DO I=1, 18
        WRITE(IPHEL,*) '*'
      ENDDO
      WRITE(IPHEL,*) '<<< DRAGON >>> Version: 5.0.0 <<< DRAGON >>>'
      WRITE(IPHEL,*) 'DRAGON CALCULATION BY J.TAFOREAU'

      WRITE(IPHEL,*) 'HELIOS Cases Used:'
      WRITE(IPHEL,'()')
      WRITE(IPHEL,*) '        1) IMP-operator name : kkk'
      WRITE(IPHEL,*) '        DRAGON case       : kkk'
      WRITE(IPHEL,*) '        Title(s)      1   : kkk'
      WRITE(IPHEL,'()')
      END

      SUBROUTINE  SET_RIEGO(IPDRA)
      WRITE(IPDRA,'()')
      WRITE(IPDRA,*) '(R) Area/Face names     :   unlabeled'
      WRITE(IPDRA,*) '(I) Isotope Identifiers :   unlabeled'
      WRITE(IPDRA,*) '(E) Path (STATE) idents : *            '
      WRITE(IPDRA,*) '(G) Group name          :   unlabeled'
      WRITE(IPDRA,*) '(O) Originating Group   :   unlabeled'
      WRITE(IPDRA,'()')
      END
