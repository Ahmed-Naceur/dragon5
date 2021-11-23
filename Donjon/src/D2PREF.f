*DECK D2PREF
      SUBROUTINE D2PREF( IPDAT,   NVAR,   CRDINF,   NCRD,  GRID,  PKIDX,
     >                   PKNAM, IPRINT                                 )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Select the reference state. This routine determine the reference state
* for all cases of meshing
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of info data block
* NVAR    number of state variables
* CRDINF  control rod compostition array
* NCRD    number of crontrol rod comosition
* GRID    type of griddind for branching calculation
*
*Parameters: 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT
      INTEGER NVAR,NBR,NCRD,GRID
      INTEGER CRDINF(NCRD)
      INTEGER PKIDX(NVAR)
      CHARACTER*12 PKNAM(6)
*----
*  LOCAL VARIABLES
*----
      INTEGER ITYLCM,i,IDX
      INTEGER :: IP = 2
      INTEGER STAIDX(NVAR),REFIDX(NVAR)
      INTEGER NVALPA(NVAR)
      REAL STATE(NVAR) ,REFSTA(NVAR-1),HSTSTA(NVAR-1)
      REAL VALPAR(NVAR,100)
      CHARACTER(LEN=12) PKEY(NVAR),BARNAM
      CHARACTER*12,DIMENSION(6) :: PKREF
      DATA PKREF/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/
      ! RECOVER INFORMATION FROM INFO DATA BLOCK
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMGTC(IPDAT,'STATE_VAR',12,NVAR,PKEY)

      !INITIALIZATION OF THE NUMBER OF BRANCHES TO BE CALCULATED

      CALL XDISET(VALPAR,100*NVAR,0)
      NBR=1
      DO i=1, NVAR
         IF (PKIDX(i).EQ.-1) THEN
          IDX=1
         ELSE
          IDX=PKIDX(i)
         ENDIF
         CALL LCMLEN(IPDAT,PKREF(IDX),NVALPA(i),ITYLCM)
         CALL LCMGET(IPDAT,PKREF(IDX),VALPAR(i,1:NVALPA(i)))
      ENDDO

      DO i=1, NVAR
         IF (PKIDX(i).EQ.-1) THEN
          IDX=1
         ELSE
          IDX=PKIDX(i)
         ENDIF
        ! ATTRIBUTION OF VALUES FOR THE BARR PARAMETERS
         IF (PKREF(IDX)==PKREF(1)) THEN
          BARNAM=PKNAM(1)
          REFSTA(1)= CRDINF(1)
          HSTSTA(1)= CRDINF(1)
          REFIDX(1)=1         ! INITIALIZATION OF BARR REFERENCE INDEX
          STATE(1)=CRDINF(1)  ! ATTRIBUTION OF  CONTROL ROD COMPOSITION
          STAIDX(1)=1    ! ATTRIBUTION OF  CONTROL ROD COMPOSITION INDEX
          NBR=NBR*NVALPA(i)   ! CALCULATION OF NUMBER OF BRANCHES
         ! IDEM FOR BURN PARAMETERS
         ELSE IF (PKREF(IDX)==PKREF(6)) THEN
          STATE(NVAR)=VALPAR(i,1)
          STAIDX(NVAR)=1
          REFIDX(NVAR)=1
          !IDEM FOR OTHER PARAMETERS
         ! EXIT

         ELSE

         ! THE REFERENCE STATES IS SET TO THE MIDDLE VALUE IN THE LIST
          REFSTA(IP)=VALPAR(i,NINT(NVALPA(i)/2.0))
          HSTSTA(IP)= VALPAR(i,NINT(NVALPA(i)/2.0))
          REFIDX(IP)=NINT(NVALPA(i)/2.0)
          STATE(IP)=REFSTA(IP)
          STAIDX(IP)=NINT(NVALPA(i)/2.0)
          NBR=NBR*NVALPA(i)
          IP=IP+1
         ENDIF
      ENDDO

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMPUT(IPDAT,'PRINT',1,1,1)

      IF((NBR>9999).OR.(GRID==0)) THEN
         ! IN THE CASE WHERE THE NUMBER OF BRANCHES EXCEED 999, A
         ! DEFAULT BANCHING CALCULATION IS CALLED
         GRID = 0
         CALL D2PDEF( IPDAT,  PKEY, VALPAR,   NVALPA,  STAIDX, REFIDX,
     >                REFSTA,HSTSTA,  STATE,   CRDINF,    NCRD,  NVAR,
     >                PKIDX ,IPRINT                                  )
      ELSE
         ! UPDATE THE INFO DATA BLOCK
         ! WITH THE INITIAL MESHING FROMSAPHYB
         CALL LCMSIX(IPDAT,' ',0)
         CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
         CALL LCMPTC(IPDAT,'BRANCH',12,1,BARNAM)
         CALL LCMPUT(IPDAT,'BRANCH_IT',1,1,1)
         CALL LCMPUT(IPDAT,'REF_STATE',NVAR-1,2,REFSTA)
         CALL LCMPUT(IPDAT,'HST_STATE',NVAR-1,2,REFSTA)
         CALL LCMPUT(IPDAT,'REF_INDEX',NVAR,1,REFIDX)
         CALL LCMPUT(IPDAT,'BRANCH_NB',1,1,NBR)
         CALL LCMPUT(IPDAT,'STATE',NVAR,2,STATE)
         CALL LCMPUT(IPDAT,'STATE_INDEX',NVAR,1,STAIDX)
         CALL LCMPUT(IPDAT,'BRANCH_INDEX',1,1,1)
         CALL LCMPUT(IPDAT,'REWIND',1,1,1)
         CALL LCMPUT(IPDAT,'STOP',1,1,0)

         IF(IPRINT > 1)  THEN
          WRITE(6,*)
          WRITE(6,*) "*** INFORMATION ABOUT BRANCHING CALCULATION  ***"
          WRITE(6,*)
          WRITE(6,*) "DEFAULT MESHING (Y/N) : N"
          IF(GRID==4) WRITE(6,*) "MESHING: NEW GRID WITH ADDITIONAL PTS"
          IF(GRID==3) WRITE(6,*) "MESHING: SAP/MCO WITH ADDITIONAL PTS"
          IF(GRID==2) WRITE(6,*) "MESHING: USER DEFINED "
          IF(GRID==1) WRITE(6,*) "MESHING: SAP/MCO "
          WRITE(6,*) "STATE PARAMETERS : ",PKEY(1:NVAR)
          WRITE(6,*) "REFERENCE STATES VALUES :", REFSTA
          WRITE(6,*) "INITIAL STATES VALUES :", STATE
          WRITE(6,*) "INITIAL STATES INDEX VALUES :", STAIDX
         ENDIF

      ENDIF
      END
