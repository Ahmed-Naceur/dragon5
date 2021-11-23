*DECK MCCMOD
      SUBROUTINE MCCMOD(IMPX,IPMAP,NPARAM,NCH,NB,REC1,VAL,MODTYPE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modification of the data stored in the PARAM folder of a fuel map
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M. Cordiez
*
*Parameters: input
* IMPX    printing index (=0 for no print).
* IPMAP   pointer of the fuel map
* NPARAM  number of parameters in the PARAM folder
* NCH     number of fuel channels
* NB      number of fuel bundles per channel
* REC1    record to be updated
* VAL    uniform value (real) that is to be set to REC1
* MODTYPE type of modification (0: new uniform value, 2: value added)
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NPARAM,NCH,NB,MODTYPE
      REAL VAL
      CHARACTER REC1*40
      TYPE(C_PTR) IPMAP
*----
*  LOCAL VARIABLES
*----
      INTEGER VALSIZE,PTYPE
      CHARACTER PNAME*12
      REAL, ALLOCATABLE, DIMENSION(:) :: VALMOD
      TYPE(C_PTR) JPMAP,KPMAP
*----
*  RECOVERY OF L_MAP PARAMETERS
*----
      JPMAP=LCMGID(IPMAP,'PARAM')
      DO IPAR=1,NPARAM,1
        KPMAP=LCMGIL(JPMAP,IPAR)
        CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
        CALL LCMGET(KPMAP,'P-TYPE',PTYPE)
        IF(PNAME.EQ.REC1) THEN
          IF(IMPX.GE.3) CALL LCMLIB(KPMAP)
          EXIT
        ENDIF
      ENDDO

* Checking of the type (local or global) of REC1
      IF((PTYPE.NE.1).AND.(PTYPE.NE.2)) THEN
        CALL XABORT('@MCCMOD: '//PNAME//'IS NEITHER LOCAL NOR GLOBAL'
     >         //'AND THAT IS IMPOSSIBLE.') 
      ENDIF
      IF(PTYPE.EQ.1) THEN
        IF(IMPX.GE.1) WRITE(6,210) 'PARAMETER ',PNAME,' IS GLOBAL'
        VALSIZE=1
        ALLOCATE(VALMOD(VALSIZE))
      ELSE
        IF(IMPX.GE.1) WRITE(6,210) 'PARAMETER ',PNAME,' IS LOCAL'
        VALSIZE=NCH*NB
        ALLOCATE(VALMOD(VALSIZE))
      ENDIF

* Modification of REC1
      IF(MODTYPE.EQ.0) THEN
        VALMOD=VAL
      ELSE IF(MODTYPE.EQ.1) THEN
        CALL LCMGET(KPMAP,'P-VALUE',VALMOD)
        VALMOD=VALMOD+VAL
      ENDIF

      CALL LCMPUT(KPMAP,'P-VALUE',VALSIZE,2,VALMOD)

      IF((MODTYPE.EQ.0).AND.(IMPX.GT.0)) THEN
        WRITE(6,220) 'EVERY VALUE OF THE RECORD ',REC1,' HAS BEEN '
     >             //'SET TO ',VAL,'.'
      ELSE IF((MODTYPE.EQ.1).AND.(IMPX.GT.0)) THEN
        WRITE(6,220) 'EVERY VALUE OF THE RECORD ',REC1,' HAS BEEN '
     >             //'INCREASED OF ',VAL,'.'
      ENDIF

* Array deallocation
      DEALLOCATE(VALMOD)

      RETURN

  210 FORMAT(1X,A,A6,A/)
  220 FORMAT(1X,A,A6,A,F7.2,A/)
      END
