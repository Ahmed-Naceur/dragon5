*DECK ERRABS
      SUBROUTINE ERRABS(IPMAC,NREG2,NREG,NGRP,XABS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover absorption cross sections from the macrolib.
*
*Copyright:
* Copyright (C) 2016 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC   pointer to the macrolib.
* NREG2   number of regions in the absorption array.
* NREG    number of regions in the macrolib.
* NGRP    number of energy groups in the macrolib.
* XABS    absorption cross sections.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC
      INTEGER NREG2,NREG,NGRP
      REAL XABS(NREG,NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC,KPMAC
      INTEGER, DIMENSION(:), ALLOCATABLE :: NJJ,IJJ,IPOS
      REAL, DIMENSION(:), ALLOCATABLE :: TOTAL,XSIGS,XSCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NREG),IJJ(NREG),IPOS(NREG))
      ALLOCATE(TOTAL(NREG),XSIGS(NREG),XSCAT(NREG*NGRP))
*
      CALL XDRSET(XABS,NREG*NGRP,0.0)
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO IGR=1,NGRP
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'NTOT0',TOTAL)
        CALL LCMLEN(KPMAC,'SIGS00',ILONG,ITYLCM)
        IF(ILONG.NE.0) THEN
          CALL LCMGET(KPMAC,'SIGS00',XSIGS)
          DO I=1,NREG2
            XABS(I,IGR)=XABS(I,IGR)+TOTAL(I)-XSIGS(I)
          ENDDO
        ELSE
          CALL LCMGET(KPMAC,'NJJS00',NJJ)
          CALL LCMGET(KPMAC,'IJJS00',IJJ)
          CALL LCMGET(KPMAC,'IPOS00',IPOS)
          CALL LCMGET(KPMAC,'SCAT00',XSCAT)
          DO I=1,NREG2
            XABS(I,IGR)=XABS(I,IGR)+TOTAL(I)
            IPO=IPOS(I)
            J2=IJJ(I)
            J1=IJJ(I)-NJJ(I)+1
            DO JGR=J2,J1,-1
              XABS(I,JGR)=XABS(I,JGR)-XSCAT(IPO)
              IPO=IPO+1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSCAT,XSIGS,TOTAL)
      DEALLOCATE(IPOS,IJJ,NJJ)
      RETURN
      END
