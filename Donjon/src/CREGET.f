*DECK CREGET
      SUBROUTINE CREGET(IPMAP,NCH,NB,IBTYP,IMPX,BRN0,BRN1,FMIX,ZONEDP,
     1 IVARTY,VARVAL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* recover the necessary information from the fuel-map object.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki, A. Hebert
*
*Parameters: input
* IPMAP   pointer to the fuel-map information.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* IBTYP   type of interpolation:
*         =0 not provided;
*         =1 time-average;
*         =2 instantaneous;
*         =3 derivative with respect to a single exit burnup.
* IMPX    printing index (=0 for no print).
* IVARTY  index of the exit burnup used to compute derivatives;
*         used if IBTYP=3.
*
*Parameters: output
* FMIX    fuel mixture indices per fuel bundle.
* BRN0    contains either low burnup integration limits or 
*         instantaneous burnups per fuel bundle.
* BRN1    upper burnup integration limits per fuel bundle.
* VARVAL  single exit burnup; used if IBTYP=3.
* ZONEDP  switch related to Chambon formula.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,IBTYP,IMPX,FMIX(NCH,NB),ZONEDP(NCH,NB),IVARTY
      REAL BRN0(NCH,NB),BRN1(NCH,NB),VARVAL
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6,NSTATE=40)
      INTEGER ISTATE(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IZONE
      REAL, ALLOCATABLE, DIMENSION(:) :: VARC
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IZONE(NCH))
*
      CALL XDISET(FMIX,NCH*NB,0)
      CALL LCMGET(IPMAP,'FLMIX',FMIX)
      CALL XDRSET(BRN0,NCH*NB,0.0)
      CALL XDRSET(BRN1,NCH*NB,0.0)
      IF(IBTYP.EQ.0) THEN
         CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
         IBTYP=ISTATE(5)
      ENDIF
*----
*  TIME-AVERAGE
*----
      IF(IBTYP.EQ.1)THEN
*       LOW BURNUP LIMITS
        CALL LCMLEN(IPMAP,'BURN-BEG',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@CREGET: MISSI'
     1   //'NG BURN-BEG VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-BEG',BRN0)
*       UPPER BURNUP LIMITS
        CALL LCMLEN(IPMAP,'BURN-END',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@CREGET: MISSI'
     1   //'NG BURN-END VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-END',BRN1)
        IF(IMPX.GT.0)WRITE(IOUT,1000)
*----
*  INSTANTANEOUS
*----
      ELSEIF(IBTYP.EQ.2)THEN
        CALL LCMLEN(IPMAP,'BURN-INST',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@CREGET: MISSI'
     1   //'NG BURN-INST VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-INST',BRN0)
        IF(IMPX.GT.0)WRITE(IOUT,1001)
*----
*  SINGLE EXIT BURNUP
*----
      ELSEIF(IBTYP.EQ.3)THEN
        IF(IVARTY.EQ.0)CALL XABORT('@CREGET: IVARTY NOT SET.')
*       LOW BURNUP LIMITS
        CALL LCMLEN(IPMAP,'BURN-BEG',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@CREGET: MISSI'
     1   //'NG BRN0 VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-BEG',BRN0)
*       UPPER BURNUP LIMITS
        CALL LCMLEN(IPMAP,'BURN-END',LENGT,ITYP)
        IF(LENGT.EQ.0)CALL XABORT('@CREGET: MISSI'
     1   //'NG BRN1 VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-END',BRN1)
        IF(IMPX.GT.0)WRITE(IOUT,1000)
        CALL LCMGET(IPMAP,'B-ZONE',IZONE)
        DO 35 ICH=1,NCH
        DO 30 IB=1,NB
        IF(IZONE(ICH).EQ.IVARTY)THEN
          ZONEDP(ICH,IB)=1
        ELSE
          ZONEDP(ICH,IB)=0
        ENDIF
   30   CONTINUE
   35   CONTINUE
        CALL LCMLEN(IPMAP,'BURN-AVG',ILONG,ITYP)
        IF (ILONG.EQ.0)CALL XABORT('@CREGET: NO SAVED VA'
     1  //'LUES FOR THIS TYPE OF VARIABLE IN FUEL MAP')
        ALLOCATE(VARC(ILONG))
        CALL LCMGET(IPMAP,'BURN-AVG',VARC)
        VARVAL=VARC(IVARTY)
        DEALLOCATE(VARC)
      ELSE
        CALL XABORT('@CREGET: INVALID OPTION IBTYP.')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IZONE)
      RETURN
*
 1000 FORMAT(/1X,'** PERFORMING THE TIME-AVERAGE',
     1 1X,'INTEGRATION OVER THE FUEL LATTICE **'/)
 1001 FORMAT(/1X,'** PERFORMING THE INSTANTANEOU',
     1'S INTERPOLATION OVER THE FUEL LATTICE **'/)
      END
