*DECK FLDREL
      SUBROUTINE FLDREL(RELAX,IPLIST,NGRP,NUN,ARRAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Relaxation procedure for flux distribution information.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): A. Hebert
*
*Parameters: input
* RELAX   relaxation factor
* IPLIST  pointer to object information.
* NGRP    number of energy groups
* NUN     number of unknowns per energy group
* ARRAY   real record to relax 
*
*Parameters: output
* ARRAY   real record after relaxation
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER NGRP,NUN
      REAL RELAX,ARRAY(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      INTEGER IGR,IUN,ILONG,ITYLCM
      REAL, ALLOCATABLE, DIMENSION(:) :: ARRAY0
*
      IF(RELAX.EQ.1.0) RETURN
      ALLOCATE(ARRAY0(NUN))
      DO IGR=1,NGRP
        CALL LCMLEL(IPLIST,IGR,ILONG,ITYLCM)
        IF(ILONG.NE.NUN) CALL XABORT('FLDREL: UNABLE TO RELAX.')
        CALL LCMGDL(IPLIST,IGR,ARRAY0)
        DO IUN=1,NUN
          ARRAY(IUN,IGR)=RELAX*ARRAY(IUN,IGR)+(1.0-RELAX)*ARRAY0(IUN)
        ENDDO
      ENDDO
      DEALLOCATE(ARRAY0)
      RETURN
      END
