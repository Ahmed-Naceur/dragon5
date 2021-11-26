*DECK NXTTRS
      FUNCTION NXTTRS(ITRCUR,ISYM)
*
*----------
*
*Purpose:
* Find new DRAGON TURN factor associated after a Cartesian symmetry
* is applied on an old DRAGON TURN factor.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* ITRCUR  initial turn factor.
* ISYM    symmetry to consider where:
*         =-1 indicates $Z$ reflection symmetry;
*         = 1 indicates $X$ reflection symmetry;
*         = 2 indicates $X=Y$ diagonal symmetry;
*         = 3 indicates $Y$ reflection symmetry;
*         = 4 indicates $X=-Y$ diagonal symmetry.
*
*Parameters: output
* NXTTRS  turn factor after symmetry is applied.
*
*----
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          ITRCUR,ISYM
      INTEGER          NXTTRS
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTRS')
      INTEGER          MAXTUR
      PARAMETER       (MAXTUR=12)
*----
*  Local variables
*----
      INTEGER          ISZ,ICTP,IRTR
*----
*  Test input data
*  ITRCUR can be 1-8 or 13-20
*  ISYM can be 1-4
*----
      IF((ITRCUR .LE.  0) .OR.
     >   (ITRCUR .GE.  9 .AND. ITRCUR .LE. 12) .OR.
     >   (ITRCUR .GE. 21)) CALL XABORT(NAMSBR//
     >  ': Invalid TURN')
      IF(ISYM .LT. -1 .OR.
     >   ISYM .EQ.  0 .OR.
     >   ISYM .GT. 4 ) CALL XABORT(NAMSBR//
     >  ': Invalid symmetry')
*----
*  Find current symmetry factor in Z (ISZ) and current turn
*  number in plane X-Y (ICTP)
*----
      ISZ=((ITRCUR-1)/MAXTUR)
      ICTP=MOD(ITRCUR-1,MAXTUR)+1
      IF(ISYM .EQ. -1) THEN
*----
*  Z symmetry
*----
        ISZ=(1-ISZ)
      ELSE
*----
*  X, X-Y AND Y symmetry
*----
        IRTR=((ICTP-1)/4)*4
        ICTP=MOD(4-ICTP+IRTR+ISYM,4)+5-IRTR
      ENDIF
      NXTTRS=ICTP+MAXTUR*ISZ
      RETURN
      END
