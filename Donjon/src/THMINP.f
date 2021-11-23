*DECK THMINP
      SUBROUTINE THMINP(HNAME,NCH,VECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read channel-dependent data.
*
*Copyright:
* Copyright (C) 2018 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* HNAME   character*8 name of the data
* NCH     number of channels
*
*Parameters: output
* VECT    data vector
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*(*) HNAME
      INTEGER NCH
      REAL VECT(NCH)
*----
*  LOCAL VARIABLES
*----
      INTEGER ITYP,NITMA,ICH
      REAL FLOT
      CHARACTER TEXT*8,HSMG*131
      DOUBLE PRECISION DFLOT
*
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.EQ.1) THEN
        FLOT=REAL(NITMA)
        CALL XDRSET(VECT,NCH,FLOT)
      ELSE IF(ITYP.EQ.2) THEN
        CALL XDRSET(VECT,NCH,FLOT)
      ELSE IF((ITYP.EQ.3).AND.(TEXT.EQ.'CHAN')) THEN
        DO ICH=1,NCH
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
          IF(ITYP.EQ.1) THEN
            VECT(ICH)=REAL(NITMA)
          ELSE IF(ITYP.EQ.2) THEN
            VECT(ICH)=FLOT
          ELSE
             WRITE(HSMG,'(14H@THMINP: NAME=,A,21H. INTEGER OR REAL VAL,
     >       12HUE EXPECTED.)') HNAME
             CALL XABORT(HSMG)
          ENDIF
        ENDDO
      ELSE
        WRITE(HSMG,'(14H@THMINP: NAME=,A,26H. SINGLE INTEGER OR REAL V,
     >  30HALUE OR CHAN KEYWORD EXPECTED.)') HNAME
        CALL XABORT(HSMG)
      ENDIF
      RETURN
      END
