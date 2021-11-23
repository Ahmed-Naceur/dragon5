*DECK NCRMAP
      SUBROUTINE NCRMAP(IPMAP,NPARM,HPARM,NCH,NB,IBTYP,IMPX,BURN0,
     1 BURN1,WPAR,LPARM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* rRcover global parameter values from the fuel-map object.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki, R. Chambon
*
*Parameters: input
* IPMAP   pointer to the fuel-map information.
* NPARM   number of expected global parameters to be recovered from
*         the fuel-map (burnup not included).
* HPARM   names of these global parameters.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* IBTYP   type of interpolation:
*         =0 not provided; =1 time-average; =2 instantaneous;
*         =3 derivative with respect to a single exit burnup.
* IMPX    printing index (=0 for no print).
*
*Parameters: output
* BURN0   contains either low burnup integration limits or
*         instantaneous burnups per fuel bundle.
* BURN1   upper burnup integration limits per fuel bundle.
* WPAR    values of the other global parameters in each bundle.
* LPARM   existence flag for each expected global parameters.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NPARM,NCH,NB,IBTYP,IMPX
      REAL BURN0(NCH,NB),BURN1(NCH,NB),WPAR(NCH,NB,NPARM)
      LOGICAL LPARM(NPARM+1)
      CHARACTER HPARM(NPARM+1)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::NSTATE=40
      INTEGER  ISTATE(NSTATE)
      INTEGER IB, ICH, ILONG, ITYLCM, ITYPEP, JPARM
      REAL VARTMP
      CHARACTER HSMG*131
      TYPE(C_PTR) JPMAP,KPMAP
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: BURNB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(BURNB(NCH,NB))
*----
*  TIME-AVERAGE BURNUP
*----
      BURN0(:NCH,:NB)=0.0
      BURN1(:NCH,:NB)=0.0
      WPAR(:NCH,:NB,:NPARM)=0.0
      LPARM(:NPARM+1)=.FALSE.
      IF(IBTYP.EQ.0) THEN
         CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
         IBTYP=ISTATE(5)
      ENDIF
      IF((IBTYP.EQ.1).OR.(IBTYP.EQ.3))THEN
*       LOW BURNUP LIMITS
        CALL LCMLEN(IPMAP,'BURN-BEG',ILONG,ITYLCM)
        IF(ILONG.EQ.0)CALL XABORT('@NCRMAP: MISSING'
     1   //' BURN0 VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-BEG',BURN0)
*       UPPER BURNUP LIMITS
        CALL LCMLEN(IPMAP,'BURN-END',ILONG,ITYLCM)
        IF(ILONG.EQ.0)CALL XABORT('@NCRMAP: MISSING'
     1   //' BURN1 VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-END',BURN1)
        IF(IMPX.GT.0)WRITE(IOUT,1000)
        LPARM(NPARM+1)=.TRUE.
*----
*  INSTANTANEOUS BURNUP
*----
      ELSEIF(IBTYP.EQ.2)THEN
        CALL LCMLEN(IPMAP,'BURN-INST',ILONG,ITYLCM)
        IF(ILONG.EQ.0)CALL XABORT('@NCRMAP: MISSING'
     1   //' BURN-INST VALUES IN FUEL MAP.')
        CALL LCMGET(IPMAP,'BURN-INST',BURNB)
        DO ICH=1,NCH
          DO IB=1,NB
            BURN0(ICH,IB)=BURNB(ICH,IB)
            BURN1(ICH,IB)=BURNB(ICH,IB)
          ENDDO
        ENDDO
        IF(IMPX.GT.0)WRITE(IOUT,1001)
        LPARM(NPARM+1)=.TRUE.
      ELSEIF(IBTYP.NE.0)THEN
        CALL XABORT('@NCRMAP: INVALID BURNUP INTERPOLATION OPTION '
     1  //'IBTYP IN FUEL MAP.')
      ENDIF
*----
*  RECOVER OTHER PARAMETERS
*----
      IF(NPARM.GT.0) THEN
        JPMAP=LCMGID(IPMAP,'PARAM')
        DO 30 JPARM=1,NPARM
          KPMAP=LCMGIL(JPMAP,JPARM)
          CALL LCMGTC(KPMAP,'PARKEY',12,1,HPARM(JPARM))
          CALL LCMGET(KPMAP,'P-TYPE',ITYPEP)
          LPARM(JPARM)=.TRUE.
*       Global parameter
          IF(ITYPEP.EQ.1) THEN
            CALL LCMLEN(KPMAP,'P-VALUE',ILONG,ITYLCM)
            IF(ILONG.NE.1) THEN
              WRITE(HSMG,'(37H@NCRMAP: P-VALUE LENGTH OF PARAMETER ,A,
     1        12H IS EQUAL TO,I6,13H (MUST BE 1).)') HPARM(JPARM),ILONG
              CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(KPMAP,'P-VALUE',VARTMP)
            WPAR(:NCH,:NB,JPARM)=VARTMP
*       Local parameter
          ELSEIF (ITYPEP.EQ.2) THEN
            CALL LCMGET(KPMAP,'P-VALUE',WPAR(1,1,JPARM))
          ENDIF
   30   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(BURNB)
      RETURN
*
 1000 FORMAT(/1X,'** PERFORMING THE TIME-AVERAGE',
     1 1X,'INTEGRATION OVER THE FUEL LATTICE **'/)
 1001 FORMAT(/1X,'** PERFORMING THE INSTANTANEOU',
     1'S INTERPOLATION OVER THE FUEL LATTICE **'/)
      END
