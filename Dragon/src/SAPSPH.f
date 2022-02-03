*DECK SAPSPH
      SUBROUTINE SAPSPH(IPEDIT,NG,NMIL,ILOC,NLOC,RVALOC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover a set of sph equivalence factors and store them as local
* variables.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPEDIT  pointer to the edition object (L_EDIT signature).
* NG      number of condensed energy groups.
* NMIL    number of mixtures in the Saphyb.
* ILOC    position of local parameter in RVALOC.
* NLOC    first dimension of matrix RVALOC.
*
*Parameters: output
* RVALOC  local variable values in mixtures located in RVALOC(ILOC,:).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT
      INTEGER NG,NMIL,ILOC,NLOC
      REAL RVALOC(NLOC,NMIL)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPEDIT,KPEDIT
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT12*12
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NMIL))
*
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,TEXT12)
      CALL LCMSIX(IPEDIT,TEXT12,1)
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      CALL LCMGET(IPEDIT,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NG) CALL XABORT('SAPSPH: BAD VALUE OF NG.')
      IF(ISTATE(2).NE.NMIL) CALL XABORT('SAPSPH: BAD VALUE OF NMIL.')
*----
*  RECOVER SPH EQUIVALENCE FACTORS.
*----
      JPEDIT=LCMGID(IPEDIT,'GROUP')
      DO 30 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
      CALL LCMLEN(KPEDIT,'NSPH',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
         CALL LCMGET(KPEDIT,'NSPH',WORK)
         DO 10 IMIL=1,NMIL
         RVALOC(ILOC+IGR-1,IMIL)=WORK(IMIL)
   10    CONTINUE
      ELSE
         DO 20 IMIL=1,NMIL
         RVALOC(ILOC+IGR-1,IMIL)=1.0
   20    CONTINUE
      ENDIF
   30 CONTINUE
      CALL LCMSIX(IPEDIT,' ',2)
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
