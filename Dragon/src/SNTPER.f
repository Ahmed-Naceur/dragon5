*DECK SNTPER
      SUBROUTINE SNTPER (JP,IDIR,LX,LT4,IP,IENV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the permutation vectors in 2-D hexagonal geometry (heavily
* adapted from BIVPER.f).
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud, A. A. Calloo
*
*Parameters: input
* JP      index of first hexagon.
* IDIR    choice of direction (=1 W axis ; =2 X axis ; =3 Y axis).
* LX      number of hexagons, including virtual hexagons.
* LT4     number of non virtual hexagons.
* IENV    index of non virtual hexagon corresponding to each hexagon.
*
*Parameters: output
* IP      permutation vector.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER JP,IDIR,LX,LT4,IP(LX),IENV(LX)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LPAS
*
      LPAS = .TRUE.
      DO 10 I=1,LT4
      IP(I)=0
 10   CONTINUE
      NC = INT((SQRT(REAL((4*LX-1)/3))+1.)/2.)
      IFACE1 = 0
      IFACE2 = 0
      IFACE3 = 0
      IF(IDIR.EQ.1) THEN
         IFACE1 = 6
         IFACE2 = 2
         IFACE3 = 1
      ELSE IF(IDIR.EQ.2) THEN
         IFACE1 = 1
         IFACE2 = 3
         IFACE3 = 2
      ELSE IF(IDIR.EQ.3) THEN
         IFACE1 = 2
         IFACE2 = 4
         IFACE3 = 3
      ELSE IF(IDIR.EQ.4) THEN
         IFACE1 = 3
         IFACE2 = 5
         IFACE3 = 4
      ELSE IF(IDIR.EQ.5) THEN
         IFACE1 = 4
         IFACE2 = 6
         IFACE3 = 5
      ELSE IF(IDIR.EQ.6) THEN
         IFACE1 = 5
         IFACE2 = 1
         IFACE3 = 6
      ELSE
         CALL XABORT('SNTPER: INVALID DATA')
      ENDIF
*
      JI = JP
      JS = JP
      KEL = 0
      M = JI + 1
      IF(JI.EQ.LX) M = JI -(((NC -1)*6) -1)
*
      IF(IENV(JI).GT.0) THEN
         KEL = KEL + 1
         IP(IENV(JI)) = 1
      ENDIF
*
      IC = JP + NC - 1
*
 20   IF(KEL.EQ.LT4) RETURN
*
      IF(JI.EQ.IC) THEN
         IFACE2 = IDIR + 3
         CALL XABORT('SNTPER: JI == IC! UNSURE OF THE IF CONDITION. '
     1 //'CHECK CODE!')
      ENDIF
*
 30   IF(M.LE.LX) THEN
         IF(IENV(M).GT.0) THEN
            KEL = KEL + 1
            IP(IENV(M)) = KEL
         ENDIF
         JI = M
         M = NEIGHB(JI,IFACE1,9,LX,POIDS)
         GOTO 30
      ELSE
 40      JI = NEIGHB(JS,IFACE2,9,LX,POIDS)
         IF(JI.GT.LX.AND.LPAS) THEN
            IFACE2 = IFACE3
            LPAS = .FALSE.
            GO TO 40
         ENDIF
         M  = JI
         JS = JI
         GOTO 20
      ENDIF
      END
