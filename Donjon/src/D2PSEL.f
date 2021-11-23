*DECK D2PSEL
      SUBROUTINE D2PSEL( IPDAT,  IPINP, STAVEC,BRANCH,  ITBRAN, STAIDX,
     >                    NVAR, JOBOPT,    DEB,  FC1  ,    FC2,    FC3,
     >                    FC4,   XSC,   IPRINT                        )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Select the next branch calculation . This routine determines also
* when to stop the calculation and updates the INFO data block.
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of info data block
* IPINP   file unit of the GENPMAXS input file
* JOBOPT  array for JOBOPT configuration
* NGP     number of energy groups
* BRANCH  nature of the current branch ( CR, DC, CB, TC, TM etc )
* ITBRAN  index of the current branch
* STAIDX  array of state variables index
* NVAR    number of state variables
* STAVEC  various parameters associated with the IPDAT structure
* DEB     flag for D2PGEN
*
*Parameters: 
* FC1     
* FC2     
* FC3     
* FC4     
* XSC     
* IPRINT  
* X       
* 
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT
      INTEGER IPINP,STAVEC(40),NVAR,ITBRAN,IPRINT,DEB
      INTEGER STAIDX(NVAR)
      CHARACTER*4 BRANCH
      CHARACTER JOBOPT(16)
*----
*  LOCAL VRAIABLES
*----
      TYPE(C_PTR) IPTH,KPTH
      INTEGER CHANGE,ITYLCM,BRAIDX,PK
      INTEGER FA_K
      INTEGER :: IP = 0
      INTEGER NVAL(NVAR),REFIDX(NVAR)
      ! VALUES OF CURRENT STATE VARIABLE ( IE FOR THE CURRENT BRANCH
      ! CALCULATION)
      REAL  STATE(NVAR)
      ! VALUES OF THE CHOOSEN REFERENCE STATE VARIABLES
      REAL  REFSTA(NVAR)
      ! VALUES OF STATES VARIABLES IN SAPHYB
      REAL  VALPAR(NVAR,100)
      REAL SFAC,BFAC,IUPS,VERS,XESM
      CHARACTER*12 BARNAM
      CHARACTER*12 PKEY(NVAR),PKNAM(6)
      CHARACTER FILNAM*12,COM*40
      CHARACTER*16 JOBTIT
      CHARACTER*1 DER
      CHARACTER*12,DIMENSION(6) :: PKREF
      DATA PKREF/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/
      LOGICAL :: BRANCH_STOP = .FALSE.
      LOGICAL :: ONE_VAL = .FALSE.
      LOGICAL LFLAG(6)

      CALL XDISET(VALPAR,NVAR*100,0)
      ! RECOVER INFORMATION FROM INFO data block
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMGTC(IPDAT,'STATE_VAR',12,NVAR,PKEY)
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

      BARNAM=PKNAM(1)
      IP=0
      IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 0'
      ! RECOVER VALUES FOR STATE VARIABLES
      DO i=1,6
         IF (LFLAG(i).OR. i==1 .OR. i==6) THEN
          IP=IP+1
          CALL LCMLEN(IPDAT,PKREF(i),NVAL(IP),ITYLCM)
          CALL LCMGET(IPDAT,PKREF(i),VALPAR(IP,1:NVAL(IP)))
         ENDIF
      ENDDO

      ! RECOVER INFORMATION ABOUT THE CURRENT BRANCH CALCULATION
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMGET(IPDAT,'STATE',STATE)
      CALL LCMGET(IPDAT,'REF_INDEX',REFIDX)
      CALL LCMGET(IPDAT,'REF_STATE',REFSTA)
      CALL LCMGET(IPDAT,'BRANCH_INDEX',BRAIDX)

      DO i=1, NVAR

        IF(BRANCH==PKEY(i)(:4)) THEN
        BRAIDX=i
!        IF (PKEY(i)(:4) == 'C-BO') CALL XABORT( 'STOP BRANCH')
        ENDIF
      ENDDO

      ! initialization of the flag: CHANGE
      CHANGE=1
   30 DO i=1, NVAR
       IF(i<=BRAIDX) THEN

        IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 1'
        ! A NEW BRANCH TYPE MUST BE SET IF THE CURRENT VALUE OF A
        ! GIVEN STATE VARIABLE IS THE LAST OF THE LIST
        IF(STAIDX(i)==NVAL(i)) THEN
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 2'
           ! WE KEEP THE FLAG CHANGE TO 1
           CHANGE=CHANGE*1
        ELSE
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 3'
          ! IF THE BRANCH INDEX CORREPOND TO THE LAST "REAL" STATE
          ! VARIABLE (IE THE STATE VARIABLE BEFORE BURN)
         IF((BRAIDX==NVAR-1)) THEN
            IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 4'
            ! THE CHANGE FLAG MUST BE SET TO FALSE
           CHANGE=0
           IF(NVAL(BRAIDX)==1) THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 5'
             ! EXCEPT IF THERE IS ONLY ONE VALUE FOR THE STATE VARIABLE
             ! IN THIS CASE THE CHANGE FLAG IS RESET TO 1
             CHANGE=1
           ENDIF
         ELSE
           ! IN OTHER CASE WE CONTINUE THE CURRENT BRANCH TYPE
           ! CALCULATION
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 6'
           CHANGE=0
           IF(NVAL(BRAIDX)==1) THEN
             ! EXCEPT IF THERE IS ONLY ONE VALUE FOR THE STATE VARIABLE
             ! IN THIS CASE THE CHANGE FLAG IS RESET TO 1
             CHANGE=1
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 7'
           ENDIF
         ENDIF
        ENDIF
       ENDIF
      ENDDO
      ONE_VAL=.FALSE.

      IF(CHANGE==1) THEN
         IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 8'
         IF(NVAL(BRAIDX+1)==1 .and. (BRAIDX >.1))THEN
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 9'
          IF((BRAIDX+1)<(NVAR)) THEN
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 10'
           BRAIDX=BRAIDX+1
           IF(NVAL(BRAIDX)==1) THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 11'
            IF(BRAIDX==NVAR-1) THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 12'
             BRANCH_STOP=.TRUE.
            ELSE
              IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 13'
              ONE_VAL=.TRUE.
            ENDIF
           ENDIF
          ELSE
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 14'
           BRANCH_STOP=.TRUE.
          ENDIF
         ENDIF

         IF(ONE_VAL) GO TO 30

         IF((BRAIDX+1)<(NVAR)) THEN
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 15'
          ! UPDATE OF THE INDEX OF THE BRANCH TYPE
          BRAIDX=BRAIDX+1
          ! UPDATE OF THE BRANCH TYPE
          BRANCH=PKEY(BRAIDX) (:4)
          ! INITIALIZATION OF THE INDEX OF THE BRANCH TYPE
          ITBRAN=1
          DO i=1,NVAR
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 16'
           IF(i<=BRAIDX) THEN
            IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 17'
            !INITIALIZATION AT THE FIRST VALUE OF STATE PARAMETERS
            STATE(i)=VALPAR(i,1)
            ! INITIALIZATION AT THE FIRST ORDER NUMBERS OF STATE
            ! PARAMETERS
            STAIDX(i)=1
            ! CASE WHERE THE REFERENCE VALUE IS THE FIRST VALUE
            ! (IE WHEN NVAL(BRAIDX) = 2)
            IF(i==BRAIDX) THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 18'
             IF(STAIDX(i)==REFIDX(i)) THEN
              IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 19'
              STAIDX(i)=2
              STATE(i)=VALPAR(i,2)
             ENDIF
            ENDIF
           ELSE
            IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 20'
            ! INITIALIZATION AT REFERENCE VALUES OF STATE PARAMETERS
            STATE(i)=VALPAR(i,REFIDX(i))
            ! INITIALIZATION AT REFERENCE ORDER NUMBERS OF STATE
            ! PARAMETERS
            STAIDX(i)=REFIDX(i)
           ENDIF
          ENDDO
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 21'
          CALL LCMSIX(IPDAT,' ',0)
          CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
          ! THE FLAG STOP IS SET TO FALSE (IE THE BRANCHING CALCULATION
          ! MUST CONTINUE)
          CALL LCMPUT(IPDAT,'STOP',1,1,0)
        ELSE
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 22'
          BRANCH_STOP=.TRUE.
        ENDIF

        IF(BRANCH_STOP) THEN
         IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 23'
         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
         ! THE FLAG STOP IS SET TO TRUE (IE THE BRANCHING CALCULATION
         ! MUST STOP)
         CALL LCMPUT(IPDAT,'STOP',1,1,1)
         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
         ! THE FLAG FOR WRITTING THE GENPMAXS.INP IS SET TO 2
         CALL LCMPUT(IPDAT,'FLAG',1,1,2)
         ! UPDATE OF THE GENPMAXS.INP FILE (MANY ARGUMENTS IN THIS CALL
         ! ARE NOT USED IN D2PGEN)
      CALL      D2PGEN( IPINP, IPDAT,   STAVEC, JOBTIT, FILNAM,    DER,
     >                   VERS,   COM,   JOBOPT,   IUPS,   FA_K,   SFAC,
     >                   BFAC,   DEB,     XESM,  FC1  ,    FC2,    FC3,
     >                    FC4,   XSC,   IPRINT                        )

        ENDIF
      ELSE
        ! update of the index of the branch type
        IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 24'
        ITBRAN=ITBRAN+1
        ! CASE WHERE THE STATE VARIABLE VALUE CORRESPOND TO THE
        ! REFERENCE STATE VALUE
        IF(STATE(BRAIDX)==REFSTA(BRAIDX)) THEN
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 25'
          ! we skip the reference value'
          STAIDX(BRAIDX)=STAIDX(BRAIDX)+1
          IF(NVAL(BRAIDX)>=1) THEN
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 26'
           ! the new value for the state variable is the next in the
           ! list
           STATE(BRAIDX)=VALPAR(BRAIDX,STAIDX(BRAIDX))
          ENDIF
        ELSE

          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 27'
          ! POSITIONNING OF THE LOOP INDEX AT THE CURRENT BRANCH TYPE
          ! CALCULATION
          i=BRAIDX
          ! DECREASE THE INDEX WHILE THE STATE VARIABLE IS BARR
          DO WHILE (i>0)
          IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 28'
          ! IF THE CURRENT VALUE OF STATE VARIABLE IS THE LAST OF THE
          ! LIST
           IF(STAIDX(i)==NVAL(i)) THEN
            IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 29'
            IF(NVAL(i)>2) THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 30'
             ! RESET OF THE ORDER NUMBERS FOR THE STATE VALUE
             STAIDX(i)=1
             ! ATTRIBUTION OF THE FIRST VALUE OF THE LIST TO THE STATE
             STATE(i)=VALPAR(i,STAIDX(i))
            ELSE
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 31'
             j=i-1
             ! INCREASE THE ORDER NUMBERS OF THE VALUE OF THIS STATE
             STAIDX(j)=STAIDX(j)+1
             ! ATTRIBUTION OF THE STATE(J) VALUES
             STATE(j)=VALPAR(j,STAIDX(j))
             ! WHILE J>0 (IE THE STATE VARIABLE EXISTS)
             DO WHILE (STAIDX(j)>NVAL(j).and.j>0)
              IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 32'
              ! IF THE STATE VARAIBLE IS NOT BARR: INITIALIZATION OF THE
              ! ORDER NUMBERS
              IF(j>1)STAIDX(j)=1
              ! IF THE STATE VARAIBLE IS NOT BARR: ATTRIBUTION OF THE
              ! STATE VARIABLE VALUE
              IF(j>1)STATE(j)=VALPAR(j,STAIDX(j))
              ! DECREASE THE J PARAMETERS
              j=j-1
              ! IF THE STATE PRAMETER EXISTS: UPDATE THE ORDER NUMBERS
              IF(j>0)STAIDX(j)=STAIDX(j)+1
              ! IF THE STATE PRAMETER EXISTS: ATTRIBUTION OF THE STATE
              ! VARIABLE VALUE
              IF(j>0)STATE(j)=VALPAR(j,STAIDX(j))
              ! EXIT OF THE IF CONDITION
             ENDDO
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 33'
             EXIT
            ENDIF
           ELSE IF(NVAL(i)==2) THEN
            IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 34'
            IF(PKEY(i).NE.BARNAM)THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 35'
             IF(STAIDX(i-1).NE.NVAL(i-1)) THEN
              IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 36'
              j=i-1
              ! INCREASE THE ORDER NUMBERS OF THE VALUE OF THIS STATE
              STAIDX(j)=STAIDX(j)+1
              ! ATTRIBUTION OF THE STATE(J) VALUES
              STATE(j)=VALPAR(j,STAIDX(j))
              EXIT
             ELSE
              IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 37'
              ! IF THE BRANCH TYPE IS BARR OR THE CURRENT STATE VALUE I$
              STAIDX(i)=STAIDX(i)+1
              IF(i>1)STAIDX(i-1)=1
              STATE(i)=VALPAR(i,STAIDX(i))
              IF(i>1)STATE(i-1)=VALPAR(i-1,STAIDX(i-1))
              EXIT
             ENDIF
            ELSE
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 38'
             IF(STAIDX(i).NE.NVAL(i)) THEN
              IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 39'
              j=i
              ! INCREASE THE ORDER NUMBERS OF THE VALUE OF THIS STATE
              STAIDX(j)=STAIDX(j)+1
              ! ATTRIBUTION OF THE STATE(J) VALUES
              STATE(j)=VALPAR(j,STAIDX(j))
              EXIT
             ENDIF
            ENDIF
           ELSE

            ! IF THE BRANCH TYPE IS BARR OR THE CURRENT STATE VALUE IS
            ! NOT THE LAST OF THE LIST
            STAIDX(i)=STAIDX(i)+1
            IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 40'
            IF((STAIDX(i)==REFIDX(i)).and.(BRANCH.NE.BARNAM)) THEN
             IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 41'
             ! IF IT IS THE REFERENCE VALUE BUT NOT THE BARR REF VALUE
             ! UPDATE THE ORDER NUMBERS OF STATE VARIABLE VALUE
             IF(i==BRAIDX) STAIDX(i)=STAIDX(i)+1
            ENDIF
            ! ATTRIBUTION OF THE STATE VARIABLE VALUE
            STATE(i)=VALPAR(i,STAIDX(i))
            EXIT
           ENDIF
           IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 42'
           i=i-1
          ENDDO
        ENDIF
      ENDIF
      IF (IPRINT>100) WRITE(6,*) '@D2PSEL : STEP 43'
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      IF((BRANCH .NE.BARNAM(:4)).and.NVAL(BRAIDX)==1) THEN
        CALL LCMPUT(IPDAT,'PRINT',1,1,0)
      ELSE
        CALL LCMPUT(IPDAT,'PRINT',1,1,1)
      ENDIF

      CALL LCMPTC(IPDAT,'BRANCH',4,1,BRANCH)
      CALL LCMPUT(IPDAT,'BRANCH_IT',1,1,ITBRAN)
      CALL LCMPUT(IPDAT,'STATE',NVAR,2,STATE)
      CALL LCMPUT(IPDAT,'STATE_INDEX',NVAR,1,STAIDX)
      CALL LCMPUT(IPDAT,'BRANCH_INDEX',1,1,BRAIDX)

      IF(IPRINT > 0)  THEN
         WRITE(6,*)
         WRITE(6,*) "**** SELECTING THE NEXT BRANCH CALCULATION ****"
         WRITE(6,*) "******     NEXT BRANCH CHARACTERISTICS    *****"
         WRITE(6,*) "BRANCH TYPE         :",BRANCH
         WRITE(6,*) "BRANCH INDEX        :",BRAIDX
         WRITE(6,*) "BRANCH ITERATION    :",ITBRAN
         WRITE(6,*) "STATE VARIABLE NAME :",PKEY
         WRITE(6,*) "BRANCH STATE VALUES :",STATE
         WRITE(6,*) "BRANCH STATE INDEX  :",STAIDX
      ENDIF
      CALL LCMSIX(IPDAT,' ',0)

      END
