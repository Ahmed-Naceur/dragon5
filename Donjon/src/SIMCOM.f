*DECK SIMCOM
      SUBROUTINE SIMCOM(IPMAP,IMPX,IMODE,NCH,NB,HC1,HC2,INDCY1,INDCY2,
     1 BURNCY1,BURNCY2,ERROR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compare two fields of values, corresponding to two cycles.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAP   fuel map
* IMPX    print parameter
* IMODE   type of field (1: burnup distribution; 2: power distribution)
* NCH     number of assemblies or number of quart-of-assemblies.
* NB      number of axial burnup subdivisions in an assembly.
* HC1     first cycle list directory in IPMAP
* HC2     first cycle list directory in IPMAP
* INDCY1  integer index in directory HCY1. INDCY1=-1 if undefined at
*         input.
* INDCY2  integer index in directory HCY2. INDCY2=-1 if undefined at
*         input.
* BURNCY1 average burnup in directory HCY1. BURNCY1=-999.0 if undefined
*         at input.
* BURNCY2 average burnup in directory HCY2. BURNCY2=-999.0 if undefined
*         at input.
*
*Parameters: output
* ERROR   discrepancy between the two distributions
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER IMPX,IMODE,NCH,NB,INDCY1,INDCY2
      REAL BURNCY1,BURNCY2,ERROR
      CHARACTER HC1*12,HC2*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAP,KPMAP
      INTEGER SIMIND
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DIST1,DIST2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DIST1(NCH,NB),DIST2(NCH,NB))
*----
*  RECOVER INFORMATION FROM THE FIRST CYCLE
*----
      JNDCY=SIMIND(IPMAP,IMPX,HC1,INDCY1,BURNCY1)
      JPMAP=LCMGID(IPMAP,HC1)
      KPMAP=LCMGIL(JPMAP,JNDCY)
      IF(IMODE.EQ.1) THEN
        CALL LCMGET(KPMAP,'BURN-INST',DIST1)
      ELSE IF(IMODE.EQ.2) THEN
        CALL LCMGET(KPMAP,'POWER-BUND',DIST1)
      ENDIF
*----
*  RECOVER INFORMATION FROM THE SECOND CYCLE
*----
      JNDCY=SIMIND(IPMAP,IMPX,HC2,INDCY2,BURNCY2)
      JPMAP=LCMGID(IPMAP,HC2)
      KPMAP=LCMGIL(JPMAP,JNDCY)
      IF(IMODE.EQ.1) THEN
        CALL LCMGET(KPMAP,'BURN-INST',DIST2)
      ELSE IF(IMODE.EQ.2) THEN
        CALL LCMGET(KPMAP,'POWER-BUND',DIST2)
      ENDIF
*----
*  COMPUTE DISCREPANCY
*----
      ERROR=0.0
      ICHMAX=0
      IBMAX=0
      IF(IMODE.EQ.1) THEN
        DO ICH=1,NCH
          DO IB=1,NB
            FLOT=ABS(DIST1(ICH,IB)-DIST2(ICH,IB))
            IF(FLOT.GE.ERROR) THEN
              ERROR=FLOT
              ICHMAX=ICH
              IBMAX=IB
            ENDIF
          ENDDO
        ENDDO
        IF(IMPX.GT.1) WRITE(6,100) ERROR,ICHMAX,IBMAX
      ELSE IF(IMODE.EQ.2) THEN
        DO ICH=1,NCH
          DO IB=1,NB
            FLOT=ABS(DIST1(ICH,IB)-DIST2(ICH,IB))/ABS(DIST2(ICH,IB))
            IF(FLOT.GE.ERROR) THEN
              ERROR=FLOT
              ICHMAX=ICH
              IBMAX=IB
            ENDIF
          ENDDO
        ENDDO
        IF(IMPX.GT.1) WRITE(6,110) ERROR,ICHMAX,IBMAX
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIST2,DIST1)
      RETURN
*
  100 FORMAT(/49H SIM: MAXIMUM DISCREPANCY ON BURNUP DISTRIBUTION=,1P,
     > E11.4,18H MW-D/T IN CHANNEL,I4,22H AND AXIAL SUBDIVISION,I4,1H./)
  110 FORMAT(/51H SIM: MAXIMUM RELATIVE ERROR ON POWER DISTRIBUTION=,1P,
     > E11.4,18H MW-D/T IN CHANNEL,I4,22H AND AXIAL SUBDIVISION,I4,1H./)
      END
