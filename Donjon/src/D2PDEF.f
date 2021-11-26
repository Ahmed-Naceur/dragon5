*DECK D2PDEF
      SUBROUTINE D2PDEF( IPDAT,  PKEY, VALPAR,   NVALPA,  STAIDX,REFIDX,
     >                  REFSTA,HSTSTA,  STATE,   CRDINF,    NCRD,  NVAR,
     >                  PKIDX,  IPRINT                                 )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Select the reference state. This routine determine the reference
* state in both cases: default meshing and initial meshing from Saphyb
* the default meshing is the folllowing :
* For other parameters than BARR and BURN, the subroutine keep three
* values from the list: the first, middle and last of Saphyb. For
* parameters BARR and BURN, all values are kept
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT      address of info data block
* NVAR       number of state variables
* NCRD       number of control rod composotion
* CRDINF     control rod compositions array
* VALPAR     array of values taken for each state variables
* STATE      state values for current branch calculation
* STAIDX     index of state values for current branch calculation
* REFSTA     values for each state variables of reference branch
* HSTSTA     values for each state variables of history branch
*
*Parameters: 
* IPDAT    
* PKEY     
* VALPAR   
* NVALPA   
* STAIDX   
* REFIDX   
* REFSTA   
* HSTSTA   
* STATE    
* CRDINF   
* NCRD     
* NVAR     
* PKIDX    
* IPRINT   
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT

      INTEGER NVAR,NCRD
      INTEGER NVALPA(NVAR),CRDINF(NCRD)
      INTEGER STAIDX(NVAR),REFIDX(NVAR)
      REAL STATE(NVAR),VALPAR(NVAR,100)
      REAL REFSTA(NVAR-1), HSTSTA(NVAR-1)
      CHARACTER*12 PKEY(NVAR)
      INTEGER PKIDX(NVAR)

*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPTH,KPTH
      INTEGER ITYLCM,i,PK,IDX(NVAR),j
      INTEGER :: NBR = 1
      ! number of values for each default state variable ( 1 if the
      ! initial number of values is less than 3, 3 otherwise)
      ! 1 : DMOD ; 2 : CBOR ; 3 : TCOM ; 4 : TMOD
      INTEGER :: DMS(5) = 0  ! NB OF VALUE FOR PARAMETER
      REAL :: REF(5) = -999.9 ! REFERENC VALUE
      REAL :: STA(5) = -999.9 ! INITIAL VALUE
      REAL :: HST(5) = -999.9 ! HISTORY VALUE
      CHARACTER*12 PKNAM(6)
      LOGICAL LFLAG(6)
      CHARACTER*12,DIMENSION(6) :: PKREF
      DATA PKREF/ "BARR","DMOD","CBOR","TCOM","TMOD","BURN"/
      REAL DEF(5,3)
      CALL XDISET(DEF,12,0)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)


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


      DO i=1, NVAR
       IF (PKIDX(i).EQ.-1) THEN
        IDX(i)=1
       ELSE
        IDX(i)=PKIDX(i)
       ENDIF
       CALL LCMLEN(IPDAT,PKREF(IDX(i)),NVALPA(i),ITYLCM)
       CALL LCMGET(IPDAT,PKREF(IDX(i)),VALPAR(i,1:NVALPA(i)))
       IF (PKREF(IDX(i)).EQ.PKREF(1)) THEN
        NBR=NBR*NVALPA(i)
        PKEY(1)=PKNAM(1)
        REFSTA(1)=CRDINF(1)
        HSTSTA(1)= CRDINF(1)
        STATE(1)= CRDINF(1)
        STAIDX(1)=1
        REFIDX(1)=1
       ENDIF
      ENDDO

      DO i=1, NVAR
       DO j=2,5
        IF (PKREF(IDX(i)).EQ.PKREF(j)) THEN
          IF(NVALPA(i)>2) THEN
           DMS(j)=3
           DEF(j,2)=VALPAR(i,NINT(NVALPA(i)/2.0))
           DEF(j,3)=VALPAR(i,NVALPA(i))
           STAIDX(j)=2   ! DMOD INDEX OF INITIAL DEFAULT VALUE
           REFIDX(j)=2   ! DMOD INDEX OF REFERENCE DEFAULT VALUE
           NBR=NBR*3
          ELSE
           DMS(j)=1
           STAIDX(j)=1   ! DMOD INDEX OF INITIAL DEFAULT VALUE
           REFIDX(j)=1   ! DMOD INDEX OF REFERENCE DEFAULT VALUE
          ENDIF
          DEF(j,1)=VALPAR(i,1)
          STA(j)=VALPAR(i,NINT(NVALPA(i)/2.0))
          HST(j)=VALPAR(i,NINT(NVALPA(i)/2.0))
          REF(j)=HST(j)
        ENDIF
       ENDDO

      ENDDO


      DO k=2,5
       IF (k==6) THEN
        PKEY(k)=PKNAM(6)
        STATE(k)= VALPAR(NVAR,1)
        STAIDX(k)=1
        REFIDX(k)=1
        CALL LCMDEL(IPDAT,PKREF(k))
        CALL LCMPUT(IPDAT,PKREF(k),NVALPA(NVAR),2,
     1  VALPAR(NVAR,1:NVALPA(NVAR)) )
       ELSE IF (LFLAG(k)) THEN
        l=k-1
        DO WHILE ((.NOT.(LFLAG(l)).and. (l.GT.1)))
        l=l-1
        ENDDO
        PKEY(l+1)=PKNAM(k)
        REFSTA(l+1)=REF(k)
        HSTSTA(l+1)=HST(k)
        STATE(l+1)=STA(k)
        CALL LCMPUT(IPDAT,PKREF(k),DMS(k),2,DEF(k,1:DMS(k)))
       ENDIF
      ENDDO

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      CALL LCMPTC(IPDAT,'STATE_VAR',12,NVAR,PKEY)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      CALL LCMPTC(IPDAT,'BRANCH',12,1,PKNAM(1))
      CALL LCMPUT(IPDAT,'BRANCH_IT',1,1,1)
      CALL LCMPUT(IPDAT,'REF_STATE',NVAR-1,2,REFSTA)
      CALL LCMPUT(IPDAT,'HST_STATE',NVAR-1,2,REFSTA)
      CALL LCMPUT(IPDAT,'BRANCH_NB',1,1,NBR)
      CALL LCMPUT(IPDAT,'STATE_INDEX',NVAR,1,STAIDX)
      CALL LCMPUT(IPDAT,'REF_INDEX',NVAR,1,REFIDX)
      CALL LCMPUT(IPDAT,'BRANCH_INDEX',1,1,1)
      CALL LCMPUT(IPDAT,'REWIND',1,1,1)
      CALL LCMPUT(IPDAT,'STATE',NVAR,2,STATE)
      CALL LCMPUT(IPDAT,'STOP',1,1,0)

      IF(IPRINT > 1)  THEN
          WRITE(6,*)
          WRITE(6,*) "*** INFORMATION ABOUT BRANCHING CALCULATION  ***"
          WRITE(6,*)
          WRITE(6,*) "DEFAULT MESHING (Y/N) : Y"
          WRITE(6,*) "==> NEW VALUES FOR PARAMTERS",
     1    " OTHER THAN BARR AND BURN :"
          WRITE(6,*) "  DMOD : ", DEF(2,1:DMS(2))
          WRITE(6,*) "  CBOR : ", DEF(3,1:DMS(3))
          WRITE(6,*) "  TCOM : ", DEF(4,1:DMS(4))
          IF(LFLAG(5)) THEN
          WRITE(6,*) "  TMOD : ", DEF(5,1:DMS(5))
          ENDIF
          WRITE(6,*)
          WRITE(6,*) "NUMBER OF BRANCHES : ", NBR
          WRITE(6,*)
          WRITE(6,*) "STATE PARAMETERS : ",PKEY(1:NVAR)
          WRITE(6,*) "REFERENCE STATES VALUES :", REFSTA
          WRITE(6,*)
          WRITE(6,*) "INITIAL STATES VALUES :", STATE
          WRITE(6,*) "INITIAL STATES INDEX VALUES :", STAIDX
          WRITE(6,*)
      ENDIF

      END
