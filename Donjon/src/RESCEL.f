*DECK RESCEL
      SUBROUTINE RESCEL(IPMAP,NCH,NK,ALCH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute fuel bundle burnups from the age pattern ALCH between 
* begin-of-cyle burnups BINI and end-of-cycle burnups BFIN
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* routine partly recovered from OPTEX-4 (coef3e)
*
*Parameters: input
* IPMAP  address of the MAP linked list or xsm file
* NCH    number of channels
* NK     number of bundles per channel
* ALCH   integer representing channel age. 
*
*Parameters: output
* IPMAP  address of the MAP linked list or xsm file
*
*Reference:
* J. Tajmouati, "Optimisation de la gestion du combustible enrichi d'un
* reacteur CANDU avec prise en compte des parametres locaux", These
* Ph. D., Ecole Polytechnique de Montreal (1993). Voir Eq. (4.7).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER     NCH,NK,ALCH(NCH)
      REAL, ALLOCATABLE, DIMENSION(:) :: F
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WINT,BINI,BFIN
*----
*  LOCAL VARIABLES
*----
      INTEGER     I,J,ILONG,ITYP
*----
*  SCRATCH STORAGE ALLOCATION
*   BINI   initial burnup map
*   BFIN   final burnup map
*   WINT   instantaneous burnup
*   F      age values in real
*----
      ALLOCATE(WINT(NCH,NK),BINI(NCH,NK),BFIN(NCH,NK),F(NCH))
*
* RECOVER FUEL BURNUPS
      CALL LCMLEN(IPMAP,'BURN-BEG',ILONG,ITYP)
      IF(ILONG.EQ.0) THEN
         CALL XABORT('SHIFTB: INITIAL BURNUP REQUIRED')
      ENDIF
      CALL LCMGET(IPMAP,'BURN-BEG',BINI)
      CALL LCMLEN(IPMAP,'BURN-END',ILONG,ITYP)
      IF(ILONG.EQ.0) THEN
         CALL XABORT('SHIFTB: FINAL BURNUP REQUIRED')
      ENDIF
      CALL LCMGET(IPMAP,'BURN-END',BFIN)
*
      DO 10 I=1,NCH
        F(I) = (FLOAT(ALCH(I)) - 0.5) / FLOAT(NCH)
        IF( ALCH(I).EQ.0 ) F(I) = 0.0
        DO 11 J=1,NK
          WINT(I,J) = BINI(I,J) + F(I) * (BFIN(I,J) - BINI(I,J))
 11     CONTINUE
 10   CONTINUE
      CALL LCMPUT(IPMAP,'BURN-INST',NCH*NK,2,WINT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(F,BFIN,BINI,WINT)
      RETURN
      END
