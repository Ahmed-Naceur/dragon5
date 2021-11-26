*DECK D2PMUL
      SUBROUTINE D2PMUL(  IPMUL,  IPDAT, STAVEC,    MIX, IPRINT        )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the global stated variable data contained in the Multicompo
* object (for reflector cross sections)
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of the INFO data block
* IPMUL   address of the MULTICOMPO object
* STAVEC  various parameters associated with the IPDAT structure
* MIX     index of mixture on which XS are to be extracted (only for
*         reflector cases)
* IPRINT  control the printing on screen
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMUL, IPDAT
      INTEGER IPRINT
      INTEGER MIX    ! MIX = 1 (RADIAL); MIX = 2 (LOW) ; MIX = 3 (TOP)
      INTEGER STAVEC(40)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPROOT,IPTH,KPTH
      PARAMETER(NSTATE=40)
      ! DEFAULT CR DC PC TF
      INTEGER :: NPAR = 5
      ! NUMBER OF CROSS SECTIONS TO BE RECOVERED
      INTEGER :: N_XS = 8
      ! NUMBER OF CB VALUES CONTAINED IN MULTICOMPO
      INTEGER CB_NB
      ! NUMBER OF VALUES FOR EACH DEFAULT STATES VARIABLES
      INTEGER NVAL(5)
      ! VALUES FOR EACH DEFAULT STATES VARIABLES
      REAL VALPAR(5,100)
      ! NAME OF PKEY
      CHARACTER (len=4) PKEY(5)
      ! NAME OF PKEY FOR BORON CONCENTRATION (MUST BE C-BORE)
      CHARACTER(LEN=6) CB_name
      ! VALUES FOR BORON CONCENTRATION
      REAL, ALLOCATABLE, DIMENSION(:) :: VAL_CB

      STAVEC(1)=2
      STAVEC(2)=NPAR
      STAVEC(3)=N_XS
      STAVEC(4)=1
      STAVEC(5)=2
      STAVEC(6)=1
      STAVEC(7)=0

      IPROOT=IPMUL
      ! MOVING AND RECOVER INFORMATION FROM MULTICOMPO
      CALL LCMSIX(IPMUL,'default',1)
      CALL LCMSIX(IPMUL,'GLOBAL',1)
      CALL LCMGTC(IPMUL,'PARKEY',6,1,CB_name)
      ! CHECK IF PKEY FOR BORON CONCENTRATION IS C-BORE
      IF(CB_name.NE.'C-BORE') THEN
        CALL XABORT('@D2PMUL: ONLY C-BORE PKEY EXPECTED')
      ENDIF
      ! RECOVER BORON CONCENTRATION VALUES
      CALL LCMLEN(IPMUL,'pval00000001',CB_NB,ITYLCM)
      ALLOCATE (VAL_CB(CB_NB))
      CALL LCMGET(IPMUL,'pval00000001',VAL_CB)

      ! CREATION OF INFO/SAPHYB_INFO/ CONTENT
      CALL LCMPUT(IPDAT,'BARR_INFO',1,1,1)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMPTC(IPDAT,'STATE_VAR',4,5,'BARRDMODCBORTCOMBURN')
      CALL LCMPUT(IPDAT,'MIX',1,1,MIX)

      ! ATTRIBUTION OF DEFAULT VALUES FOR OTHER STATE VARIABLES THAN
      ! C_BORE
      PKEY(1)='BARR'           ! CONTROL ROD
      PKEY(2)='DMOD'           ! MODERATOR DENSITY
      PKEY(3)='CBOR'           ! BORON CONCENTRATION
      PKEY(4)='TCOM'           ! FUEL TEMPERATURE
      PKEY(5)='BURN'           ! BURN UP
      ! ALL STATE VARIABLE (EXCEPT CBOR) ARE FIXED
      NVAL(1)=1
      NVAL(2)=1
      NVAL(3)= CB_NB
      NVAL(4)=1
      NVAL(5)=1
      VALPAR(1,1) = 1       ! NO CONTROL ROD IS INSERTED
      VALPAR(3,1:CB_NB) = VAL_CB
      VALPAR(2,1) = 0.75206 ! DEFAULT MODERATOR DENSITY= 0.75206 G/CM3
      VALPAR(4,1) = 560     ! FUEL TEMPERATURE= 560 Celsius
      VALPAR(5,1) = 0       ! BURN-UP= 0 MWJ/T

      ! CREATION OF INFO/SAPHYB_INFO/SVNAME
      ! LOOP OVER STATE VARIABLE
      DO i=1, NPAR
        CALL  LCMPUT(IPDAT,PKEY(i),NVAL(i),2,VALPAR(i,1:NVAL(i)))
      ENDDO

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      IPTH=LCMLID(IPDAT,'PKEY_INFO',6)
      DO J=1, 6
         KPTH=LCMDIL(IPTH,J)
         IF(J==1) THEN
           CALL LCMPTC(KPTH,"NAME",8,1,"BARR    ")
           CALL LCMPUT(KPTH,"LFLAG",1,5,.TRUE.)
         ELSE IF(J==2)THEN
           CALL LCMPTC(KPTH,"NAME",8,1,"DMOD    ")
           CALL LCMPUT(KPTH,"LFLAG",1,5,.TRUE.)
         ELSE IF(J==3) THEN
           CALL LCMPTC(KPTH,"NAME",8,1,"CBOR    ")
           CALL LCMPUT(KPTH,"LFLAG",1,5,.TRUE.)
         ELSE IF(J==4)THEN
           CALL LCMPTC(KPTH,"NAME",8,1,"TCOM    ")
           CALL LCMPUT(KPTH,"LFLAG",1,5,.TRUE.)
         ELSE IF(J==5)THEN
           CALL LCMPTC(KPTH,"NAME",8,1,"TMOD    ")
           CALL LCMPUT(KPTH,"LFLAG",1,5,.FALSE.)
         ELSE IF(J==6) THEN
           CALL LCMPTC(KPTH,"NAME",8,1,"BURN    ")
           CALL LCMPUT(KPTH,"LFLAG",1,5,.TRUE.)
         ENDIF
      ENDDO
      ! CREATION OF :
      ! INFO/HELIOS_HEAD/ DIRECTORY
      ! INFO/GENPMAXS_INP/ DIRECTORY

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'HELIOS_HEAD',1)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMSIX(IPDAT,' ',0)

      ! EDIT THE LISTING FILE
       IF(IPRINT > 0)  THEN
                   !"**************************************************"
         WRITE(6,*) "******** CONTENT OF MULTICOMPO RECOVERED *********"
         WRITE(6,*)
         WRITE(6,*) "NUMBER OF STATE VARIABLES :", NPAR
         WRITE(6,*) "NAME OF STATE VARIABLES :", PKEY
         WRITE(6,*)
         DO i=1, NPAR
          WRITE(6,*) "NUMBER OF VALUES FOR ",PKEY(i)," PARAMETERS :",
     1    NVAL(i)
          WRITE(6,*) "VALUES FOR ",PKEY(i)," PARAMETERS :",
     1    VALPAR(i,1:NVAL(i))
          WRITE(6,*)
         ENDDO
        WRITE(6,*)
      ENDIF

      ! FREE MEMORY
      DEALLOCATE (VAL_CB)
      END
