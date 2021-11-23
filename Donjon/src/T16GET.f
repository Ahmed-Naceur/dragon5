*DECK T16GET
      SUBROUTINE T16GET(MAXMIX,MNLOCP,MNCPLP,MNPERT,NALOCP,IDLCPL,
     >                  NCMIXS,MNBURN,NAMMIX,MIXRCI,PARRCI,MIXPER,
     >                  PARPER,MIXREG)
*
*----
*
*Purpose:
*  Read from input T16CPO processing options.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* MAXMIX  maximum number of mixtures.
* MNLOCP  maximum number of local parameters.
* MNCPLP  maximum number of coupled parameters.
* MNPERT  maximum number of perturbations per local parameter.
* NALOCP  local parameter names allowed.
* IDLCPL  local ID for perturbation parameters.
*
*Parameters: input/output
* NCMIXS  number of current mixtures.
* MNBURN  current and final number of burnup steps
* NAMMIX  names of mixtures.
* MIXRCI  reference information for mixtures where:
*         =0 no information for mixture;
*         >0 information not updated; 
*         <0 information to be updated. 
* PARRCI  reference local parameters for mixtures.
* MIXPER  perturbation information for mixtures.
*         =0 no information for mixture;
*         >0 information not updated; 
*         <0 information to be updated. 
* PARPER  perturbation parameters for mixtures.
*
*Parameters: output
* MIXREG  mixture update identifier where:
*          =0 do not update;
*          =-1 update using CELLAV information;
*          > 0 update using specified region number.
*
*Comments:
* Input data is of the form:
*    [[ MIXNAM [ { CELLAV | REGION noreg } ]
*    [ RC [ nburn ] frstrec ]
*    [[ NAMPER valref npert
*       (valper(i),frstrec(i),i=1,npert)]]
*    ]]
*    [ MTMD [ valreft valrefd ] npert
*       (valpert(i), valperd(i), frstrec(i),i=1,npert)]]
*    ]
*
*----
*
      IMPLICIT         NONE
      INTEGER          MAXMIX,MNLOCP,MNCPLP,MNPERT,NCMIXS,MNBURN
      CHARACTER        NALOCP(MNLOCP+MNCPLP)*4
      INTEGER          IDLCPL(2,MNLOCP+MNCPLP),NAMMIX(2,MAXMIX),
     >                 MIXRCI(2+MNLOCP+MNCPLP,MAXMIX),
     >                 MIXPER(MNPERT,MNLOCP+MNCPLP,MAXMIX),
     >                 MIXREG(MAXMIX)
      REAL             PARRCI(MNLOCP,MAXMIX),
     >                 PARPER(MNPERT,2,MNLOCP+MNCPLP,MAXMIX)
*----
*  READ VARIABLES
*----
      CHARACTER        TEXT12*12
      INTEGER          ITYPE,NITMA
      REAL             FLOTT
      DOUBLE PRECISION DFLOTT
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,NTC
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NTC=3,NAMSBR='T16GET')
      INTEGER          KCHAR(NTC),INEXTM,ILOCP,ILOCL,NLPAR,
     >                 ILPAR,NBRCI,IPAR,IMIX,IRLOC
*----
*  READ INPUT DATA.
*----
      INEXTM=0
 100  CONTINUE
      CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
 101  CONTINUE
      IF(ITYPE .NE. 3) CALL XABORT(NAMSBR//
     >  ': KEYWORD EXPECTED')
      IF(TEXT12 .EQ. ';') THEN
*----
*  END OF INPUT REACHED
*  EXIT READ
*----
        GO TO 105
      ELSE IF(TEXT12 .EQ. 'CELLAV') THEN
*----
*  CELLAV KEYWORD FOUND
*----
        IF(INEXTM .EQ. 0) CALL XABORT(NAMSBR//
     >  ': MIXTURE NAME MUST BE DEFINED BEFORE CELLAV')
        MIXREG(INEXTM)=-1
      ELSE IF(TEXT12 .EQ. 'REGION') THEN
*----
*  REGION KEYWORD FOUND
*----
        IF(INEXTM .EQ. 0) CALL XABORT(NAMSBR//
     >  ': MIXTURE NAME MUST BE DEFINED BEFORE REGION')
        CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(ITYPE .NE. 1) CALL XABORT(NAMSBR//
     >  ': REGION NUMBER MUST FOLLOW REGION KEYWORD')
        IF(NITMA .LT. 1) CALL XABORT(NAMSBR//
     >  ': REGION NUMBER MUST BE > 0')
        MIXREG(INEXTM)=NITMA
      ELSE IF(TEXT12 .EQ. 'RC') THEN
*----
*  REFERENCE CASE INFORMATION
*----
        IF(INEXTM .EQ. 0) CALL XABORT(NAMSBR//
     >  ': MIXTURE NAME MUST BE DEFINED RC')
        CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(ITYPE .NE. 1) CALL XABORT(NAMSBR//
     >  ': DATA TYPE FOLLOWING RC MUST BE INTEGER')
        IF(NITMA .LT. 1) CALL XABORT(NAMSBR//
     >  ': FIRST INTEGER VALUE FOLLOWING RC MUST BE > 0')
        MIXRCI(1,INEXTM)=NITMA
        IF(MIXRCI(2,INEXTM) .EQ. 0) THEN
          MIXRCI(2,INEXTM)=-1
        ENDIF
        CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
        IF(ITYPE .NE. 1) GO TO 101
        IF(NITMA .LT. 1) CALL XABORT(NAMSBR//
     >  ': SECOND INTEGER VALUE FOLLOWING RC MUST BE > 0')
        MNBURN=MAX(MNBURN,MIXRCI(1,INEXTM))
        MIXRCI(2,INEXTM)=-MIXRCI(1,INEXTM)
        MIXRCI(1,INEXTM)=NITMA
      ELSE
*----
*  EITHER PERTURBATION OR NEW MIXTURE
*  1) IF PERTURBATION
*     TREAT INPUT AND RETURN TO READ NEXT KEYWORD
*     OTHERWISE TEXT12 IS NEW MIXTURE NAME
*----
        IRLOC=2
        DO ILOCP=1,MNLOCP+MNCPLP
          NLPAR=1
          IF(ILOCP .GT. MNLOCP) NLPAR=2
          IF(TEXT12 .EQ. NALOCP(ILOCP)) THEN
            IF(INEXTM .EQ. 0) CALL XABORT(NAMSBR//
     >      ': MIXTURE NAME REQUIRED FOR PERTURBATIONS')
*----
*  SAVE REFERENCE PARAMETER AND TEST FOR COHERENCE
*----
            CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
            IF(ITYPE .EQ. 2) THEN
              DO ILPAR=1,NLPAR
                IF(ITYPE .NE. 2) CALL XABORT(NAMSBR//
     >          ': REFERENCES EXPECTED FOR PERTURBATIONS')
                ILOCL=IDLCPL(ILPAR,ILOCP)
                IF(MIXRCI(IRLOC+ILOCL,INEXTM) .EQ. 0) THEN
                  PARRCI(ILOCL,INEXTM)=FLOTT
                ELSE IF(PARRCI(ILOCL,INEXTM) .NE. FLOTT) THEN
                  CALL XABORT(NAMSBR//
     >            ': REFERENCE PARAMETER NOT COHERENT FOR '//
     >            NALOCP(ILOCP)//
     >            ' PERTURBATION INITIALIZATION')
                ENDIF
                CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
              ENDDO
            ELSE IF( MIXRCI(IRLOC+ILOCP,INEXTM) .EQ. 0) THEN
              CALL XABORT(NAMSBR//
     >        ': REFERENCE CASE NOT INITIALIZED FOR '//
     >        NALOCP(ILOCP)//' PERTURBATION')
            ENDIF
*----
*  READ NUMBER OF PERTURBATIONS
*----
            IF(ITYPE .NE. 1) CALL XABORT(NAMSBR//
     >      ': INVALID RECORD FOLLOWING PERTURBATION')
            IF(NITMA .LT. 0) CALL XABORT(NAMSBR//
     >      ': NUMBER OF PERTURBATION MUST BE >= 0')
            NBRCI=NITMA
            MIXRCI(IRLOC+ILOCP,INEXTM)=-NITMA
*----
*  READ PERTURBATIONS PARAMETERS
*----
            DO IPAR=1,NBRCI
              DO ILPAR=1,NLPAR
                ILOCL=IDLCPL(ILPAR,ILOCP)
                CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
                IF(ITYPE .NE. 2) CALL XABORT(NAMSBR//
     >          ': INVALID RECORD FOR REFERENCE PARAMETER')
                PARPER(IPAR,ILPAR,ILOCP,INEXTM)=FLOTT
              ENDDO
              CALL REDGET(ITYPE,NITMA,FLOTT,TEXT12,DFLOTT)
              IF(ITYPE .NE. 1) CALL XABORT(NAMSBR//
     >        ': INVALID RECORD FOLLOWING PERTURBATION')
              IF(NITMA .LT. 0) CALL XABORT(NAMSBR//
     >        ': NUMBER OF PERTURBATION MUST BE >= 0')
              MIXPER(IPAR,ILOCP,INEXTM)=NITMA
            ENDDO
            GO TO 100
          ENDIF
        ENDDO
*----
*  3) TEXT12 IS A NEW MIXTURE NAME
*     TREAT INPUT AND RETURN TO READ NEXT KEYWORD
*----
        READ(TEXT12,'(A4,A2)') KCHAR(1),KCHAR(2)
        DO IMIX=1,NCMIXS
          IF(KCHAR(1) .EQ. NAMMIX(1,IMIX) .AND.
     >       KCHAR(2) .EQ. NAMMIX(2,IMIX) ) THEN
            INEXTM=IMIX
            GO TO 145
          ENDIF
        ENDDO
        NCMIXS=NCMIXS+1
        NAMMIX(1,NCMIXS)=KCHAR(1)
        NAMMIX(2,NCMIXS)=KCHAR(2)
        IF(NCMIXS .GT. MAXMIX) CALL XABORT(NAMSBR//
     >  ': TOO MANY MIXTURES READ')
        INEXTM=NCMIXS
 145    CONTINUE
*----
*  ASSUME CELLAV BY DEFAULT
*----
        MIXREG(INEXTM)=-1
      ENDIF
      GO TO 100
 105  CONTINUE
*----
*  ALL THE REQUIRED INFORMATION READ
*  RETURN
*----
      RETURN
      END
