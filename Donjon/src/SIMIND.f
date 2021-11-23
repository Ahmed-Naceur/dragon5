*DECK SIMIND
      INTEGER FUNCTION SIMIND(IPMAP,IMPX,HCYCLE,INDCY,BURNCY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Return the list index of an existing fuel cycle.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAP   fuel map object
* IMPX    print parameter
* HCYCLE  cycle list directory in IPMAP
* INDCY   integer index in directory HCYCLE. INDCY=-1 if undefined at
*         input.
* BURNCY  average burnup in directory HCYCLE. BURNCY=-999.0 if undefined
*         at input.
*
*Parameters: output
* SIMIND  list index in HCYCLE
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER IMPX,INDCY
      CHARACTER HCYCLE*12,HSMG*131
      REAL BURNCY
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAP,KPMAP
*
      SIMIND=0
      CALL LCMLEN(IPMAP,HCYCLE,ILONG,ITYLCM)
      IF(ILONG.EQ.0) THEN
        CALL LCMLIB(IPMAP)
        WRITE(HSMG,'(24H@SIMIND: NO CYCLE NAMED ,A12,1H.)') HCYCLE
        CALL XABORT(HSMG)
      ENDIF
      IF(INDCY.NE.-1) THEN
        IF(INDCY.GT.ILONG) CALL XABORT('@SIMIND: INDCY.GT.ILONG')
        IF(BURNCY.NE.-999.0) CALL XABORT('@SIMIND: BURNCY.NE.-999.0')
        SIMIND=INDCY
        IF(IMPX.GT.0) THEN
          JPMAP=LCMGID(IPMAP,HCYCLE)
          KPMAP=LCMGIL(JPMAP,INDCY)
          CALL LCMGET(KPMAP,'BURNAVG',BURNAVG)
          WRITE(6,100) INDCY,HCYCLE,BURNAVG
        ENDIF
        RETURN
      ELSE IF((INDCY.EQ.-1).AND.(BURNCY.EQ.-999.0)) THEN
        SIMIND=ILONG
        IF(IMPX.GT.0) THEN
          JPMAP=LCMGID(IPMAP,HCYCLE)
          KPMAP=LCMGIL(JPMAP,ILONG)
          CALL LCMGET(KPMAP,'BURNAVG',BURNAVG)
          WRITE(6,100) ILONG,HCYCLE,BURNAVG
        ENDIF
        RETURN
      ENDIF
*
      DELTA=1.0E10
      BURNK=0.0
      JPMAP=LCMGID(IPMAP,HCYCLE)
      DO I=1,ILONG
        KPMAP=LCMGIL(JPMAP,I)
        CALL LCMLEN(KPMAP,'BURNAVG',ILONG,ITYLCM)
        IF(ILONG.EQ.0) CYCLE
        CALL LCMGET(KPMAP,'BURNAVG',BURNAVG)
        IF(ABS(BURNAVG-BURNCY).LT.DELTA) THEN
          SIMIND=I
          DELTA=ABS(BURNAVG-BURNCY)
          BURNK=BURNAVG
        ENDIF
      ENDDO
      IF(DELTA.GT.2.0) THEN
        WRITE(HSMG,'(47H@SIMIND: UNABLE TO FIND AN EXISTING AVERAGE BUR,
     >  12HNUP EQUAL TO,1P,E12.4,10H IN CYCLE ,A12,1H.)') BURNCY,HCYCLE
        CALL XABORT(HSMG)
      ENDIF
      IF(IMPX.GT.0) WRITE(6,100) SIMIND,HCYCLE,BURNK
      RETURN
*
  100 FORMAT(/40H SIMIND: RECOVER LIST INDEX IN LIST ITEM,I3,7H OF CYC,
     > 3HLE ,A12,10H AT BURNUP,1P,E12.4,8H MW-D/T./)
      END
