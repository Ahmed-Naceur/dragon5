*DECK SYB4VO
      SUBROUTINE SYB4VO(NSECT,NR,DX,DY,RAD,VOLINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the volumes of a Cartesian sectorized cell.
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
* NSECT   number of sectors (multiple of 4).
* NR      one plus the number of tubes.
* DX      X-oriented side.
* DY      Y-oriented side.
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
      INTEGER   NSECT,NR
      REAL      DX,DY,RAD(NR-1),VOLINT(NSECT,NR)
*
      IF(MOD(4+NSECT,4).NE.0) CALL XABORT('SYB4VO: INVALID NSECT.')
      DDX = DX / 2.
      DDY = DY / 2.
      RADDD = SQRT(DDX*DDX + DDY*DDY)
      TETA0 = 3.14159265358979 * 2. / REAL(NSECT)
*----
*  ITER=1 IS FOR (0, PI/4) AND ITER=2 IS FOR (PI/2, PI/4)
*----
      DO 15 IR=1,NR
      DO 10 IS=1,NSECT
      VOLINT(IS,IR)=0.0
   10 CONTINUE
   15 CONTINUE
      DO 60 ITER=1,2
      IF(ITER.EQ.1) THEN
         ISA = 1
         ISB = 0
         DDXY=DDX
      ELSE
         ISA = - 1
         ISB = NSECT/4 + 1
         DDXY=DDY
      ENDIF
*
      TETAC = ACOS(DDXY / RADDD)
      TETA2   = 0.
      DO 20 IR = 1, NR-1
      IF (RAD(IR) .GT. DDXY) THEN
         IRNEXT = IR
         GOTO 30
      ENDIF
   20 CONTINUE
      IRNEXT = NR
*----
*  IRNEXT : NEXT RADIUS INTERCEPTING A SIDE
*  DDY2   : LAST PROCESSED COORDINATE
*  ISNEXT : NEXT SECTOR
*----
   30 DDY2   = 0.
      ISNEXT = 1
*----
*  NEXT SECTOR
*----
   40 IS = ISNEXT
      IF (IS .NE. 0) THEN
         ISV   = ISA * IS + ISB
         IR0   = IRNEXT
         DDY1  = DDY2
         TETA1 = TETA2
         TETA2 = IS * TETA0
*
*        THE ANGLE IS LIMITED BY THE DIAGONAL.
         IF (TETA2 .GE. (TETAC - 1.E-6)) THEN
            ISNEXT = 0
            TETA2  = TETAC
            RAD2 = RADDD
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
         VOLINT(ISV, IR) = VOLINT(ISV, IR) + DT * DR
   50    CONTINUE
*
*        LAST SIDE-INTERCEPTED TUBE.
         VOLMAX = DDXY * 0.5 * (DDY2 - DDY1)
         VOLMAX = VOLMAX - DT * R1 * R1
         VOLINT(ISV, IR0) = VOLINT(ISV, IR0) + VOLMAX
         GOTO 40
      ENDIF
   60 CONTINUE
*
      DO 90 IR=1,NR
      DO 70 IS=NSECT/4+1,NSECT/2
      VOLINT(IS,IR)=VOLINT(NSECT/2-IS+1,IR)
   70 CONTINUE
      DO 80 IS=NSECT/2+1,NSECT
      VOLINT(IS,IR)=VOLINT(NSECT-IS+1,IR)
   80 CONTINUE
   90 CONTINUE
      RETURN
      END
