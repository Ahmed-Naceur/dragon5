*DECK SYB7VO
      SUBROUTINE SYB7VO(NR,HSIDES,RAD,VOLINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the volumes of an hexagonal cell.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NR      one thus the number of tubes.
* HSIDES  hexagon side.
* RAD     radius of the tubes.
*
*Parameters: output
* VOLINT  volumes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NR
      REAL HSIDES,RAD(NR-1),VOLINT(NR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (SQRT32=1.7320508075689/2.0)
*
      TETA0  = 6.0
      DO 10 IR=1,NR
      VOLINT(IR)=0.0
   10 CONTINUE
      DDXY = HSIDES * SQRT32
*
      TETAC = ACOS(DDXY / HSIDES)
      TETA2   = 0.
      DO 20 IR = 1, NR-1
      IF (RAD(IR) .GT. DDXY) THEN
         IRNEXT = IR
         GOTO 30
      ENDIF
   20 CONTINUE
      IRNEXT = NR
*----
*  IRNEXT : NEXT RADIUS INTERCEPTING A SIDE.
*  DDY2   : LAST PROCESSED COORDINATE.
*  ISNEXT : NEXT SECTOR.
*----
   30 DDY2   = 0.
      ISNEXT = 1
*----
*  NEXT SECTOR
*----
   40 IS = ISNEXT
      IF (IS .NE. 0) THEN
         IR0   = IRNEXT
         DDY1  = DDY2
         TETA1 = TETA2
         TETA2 = IS * TETA0
*
*        THE ANGLE IS LIMITED BY THE DIAGONAL.
         IF (TETA2 .GE. (TETAC - 1.E-6)) THEN
            ISNEXT = 0
            TETA2  = TETAC
            RAD2 = HSIDES
         ELSE
            ISNEXT = IS + 1
            RAD2 = DDXY / COS(TETA2)
         ENDIF
*
*        THE NEXT RADIUS IS INTERCEPTING THE SECTOR.
         IF (IR0 .LT. NR) THEN
            RADIR = RAD(IR0)
            IF (RADIR .LT. (RAD2 * (1. - 1.E-6))) THEN
               RAD2 = RADIR
               TETA2 = ACOS(DDXY / RAD2)
               IRNEXT = IR0 + 1
               ISNEXT = IS
            ELSE IF (RADIR .LE. (RAD2 * (1. + 1.E-6))) THEN
*              THE NEXT RADIUS IS EQUAL TO THE SECTOR.
               IRNEXT = IR0 + 1
            ENDIF
         ENDIF
*
*        DDY2 IS THE NEXT COORDINATE AND DT IS HALF THE ANGLE
*        INCREMENT FOR THE SECTOR.
         DDY2 = RAD2 * SIN(TETA2)
         DT = (TETA2 - TETA1) * 0.5
*
*        COMPLETE TUBES.
         R1 = 0.
         DO 50 IR = 1, IR0 - 1
         R0 = R1
         R1 = RAD(IR)
         DR = (R1 - R0) * (R1 + R0)
         VOLINT(IR) = VOLINT(IR) + DT * DR
   50    CONTINUE
*
*        LAST SIDE-INTERCEPTED TUBE.
         VOLMAX = DDXY * 0.5 * (DDY2 - DDY1)
         VOLMAX = VOLMAX - DT * R1 * R1
         VOLINT(IR0) = VOLINT(IR0) + VOLMAX
         GOTO 40
      ENDIF
*
      DO 60 I = 1, NR
      VOLINT(I) = 12.0 * VOLINT(I)
   60 CONTINUE
      RETURN
      END
