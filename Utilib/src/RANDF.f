*DECK RANDF
      SUBROUTINE RANDF(ISEED,IFIRST,RAND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* This subroutine returns a pseudo-random number for each invocation.
* It is a FORTRAN 77 adaptation of the "Integer Version 2" minimal
* standard number generator whose Pascal code appears in reference.
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
* ISEED   the seed for the generation of random numbers.
* IFIRST  set to 1 to indicate that the seed is being generated.
*                 
*Parameters: ouput
* RAND    random number between 0 and 1.
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
      INTEGER ISEED,IFIRST
      REAL    RAND
      DOUBLE PRECISION D1,D2
*----
*  LOCAL VARIABLES
*----
      INTEGER  MPLIER,MODLUS,MOBYMP,MOMDMP
      PARAMETER (MPLIER=16807,MODLUS=2147483647,MOBYMP=127773,
     +           MOMDMP=2836)
      INTEGER HVLUE, LVLUE, TESTV, NEXTN
*
      IF (IFIRST .EQ. 0) IFIRST = 1
      NEXTN = ISEED
*
      HVLUE = NEXTN / MOBYMP
      LVLUE = MOD(NEXTN, MOBYMP)
      TESTV = MPLIER*LVLUE - MOMDMP*HVLUE
      IF (TESTV .GT. 0) THEN
        NEXTN = TESTV
      ELSE
        NEXTN = TESTV + MODLUS
      ENDIF
      D1=DBLE(NEXTN)
      D2=DBLE(MODLUS)
      RAND=REAL(D1/D2)
      ISEED= NEXTN
      RETURN
      END
