*DECK RANDD
      SUBROUTINE RANDD(ISEED,NRAND,DRAND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* This subroutine returns a pseudo-random number for each invocation.
* It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
* standard number generator whose Pascal code appears in reference.
* This is the double precision version of the single precision
* routine RANDF.
* 
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input
* ISEED   the seed for the generation of random numbers. If ISEED=0
*         use ISEED=3141592654
* NRAND   number of random number requested.
*                 
*Parameters: ouput
* DRAND   random numbers between 0 and 1.
*
*Reference:                     
* Park, Steven K. and Miller, Keith W., "Random Number Generators:
* Good Ones are Hard to Find", Communications of the ACM, October 1988.
*                          
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          ISEED,NRAND
      DOUBLE PRECISION DRAND(NRAND)
*----
*  LOCAL VARIABLES
*----
      INTEGER          MPLIER,MODLUS,MOBYMP,MOMDMP
      PARAMETER       (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +                 MOMDMP=2836)
      INTEGER          IRAND,HVLUE,LVLUE,TESTV,NEXTN
*
      IF(ISEED .EQ. 0) THEN
        NEXTN = 314159265
      ELSE
        NEXTN = ISEED
      ENDIF
*
      DO IRAND=1,NRAND
        HVLUE = NEXTN / MOBYMP
        LVLUE = MOD(NEXTN, MOBYMP)
        TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
        IF (TESTV .GT. 0) THEN
          NEXTN = TESTV
        ELSE
          NEXTN = TESTV + MODLUS
        ENDIF
        DRAND(IRAND) = DBLE(NEXTN)/DBLE(MODLUS)
      ENDDO
      ISEED= NEXTN
      RETURN
      END
