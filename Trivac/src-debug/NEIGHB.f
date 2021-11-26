*DECK NEIGHB
      FUNCTION NEIGHB (J,K,IHEX,NH,POIDS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the index of a neighbour hexagon taking into account the
* symmetries.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud
*
*Parameters: input
* J       index of the considered hexagon.
* K       index of the side.
* IHEX    type of symmetry:
*         =1: S30;   =2: SA60;   =3: SB60;   =4: S90;   =5: R120;
*         =6: R180;  =7: SA180;  =8: SB180;  =9: complete;
*         =10: S30 with HBC SYME; =11: sa60 with HBC SYME.
* NH      total number of hexagons.
* POIDS   weight of the hexagon.
*
*Parameters: output
* NEIGHB  index of the neighbour hexagon. Note that:
*         ABS(NEIGHB).GT.NH:  external boundary;
*         NEIGHB=J:           reflection on side k;
*         NEIGHB.LT.0:        axial symmetry or rotation.
*
*-----------------------------------------------------------------------
*
*                             side 2
*                             xxxxxx
*                  side 3   x        x   side 1
*                         x            x
*                       x                x
*                         x            x
*                  side 4   x        x   side 6
*                             xxxxxx
*                             side 5
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER J,K,IHEX,NH
      REAL POIDS
*
      IF ((IHEX.EQ.1).OR.(IHEX.EQ.10)) THEN
         VI = 2.* SQRT(REAL(NH)) - 1.
         VP = SQRT(REAL(4*NH+1)) - 1.
         IF (AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF (AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(1).')
         ENDIF
      ELSE IF ((IHEX.EQ.2).OR.(IHEX.EQ.11)) THEN
         VA = (SQRT(REAL(8*NH+1)) - 1.)/2.
         IF (AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(2).')
         ENDIF
      ELSE IF (IHEX.EQ.3) THEN
         VI = SQRT(REAL(2*NH-1))
         VP = SQRT(REAL(2*NH))
         IF (AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF (AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(3).')
         ENDIF
      ELSE IF (IHEX.EQ.4) THEN
         VI = SQRT(REAL((4*NH-1)/3))
         VP = SQRT(REAL(4*NH/3))
         IF (AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF (AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(4).')
         ENDIF
      ELSE IF (IHEX.EQ.5) THEN
         VA = (SQRT(REAL(4*(NH-1)+1)) + 1.)/2.
         IF (AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(5).')
         ENDIF
      ELSE IF (IHEX.EQ.6) THEN
         VA = (SQRT(REAL(8*(NH-1)/3+1)) + 1)/2
         IF (AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(6).')
         ENDIF
      ELSE IF (IHEX.EQ.7) THEN
         VA = (SQRT(REAL(24*NH+1)) + 1.)/6.
         IF (AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(7).')
         ENDIF
      ELSE IF (IHEX.EQ.8) THEN
         VI = (1.+SQRT(REAL(3*(2*NH-1)+1)))/3.
         VP = (1.+SQRT(REAL(6*NH+1)))/3.
         IF (AINT(VI).EQ.VI) THEN
            NC = INT(VI)
         ELSE IF (AINT(VP).EQ.VP) THEN
            NC = INT(VP)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(8).')
         ENDIF
      ELSE IF (IHEX.EQ.9) THEN
         VA = (SQRT(REAL((4*NH-1)/3)) + 1.)/2.
         IF (AINT(VA).EQ.VA) THEN
            NC = INT(VA)
         ELSE
            CALL XABORT('NEIGHB: INVALID NUMBER OF HEXAGONS(9).')
         ENDIF
      ELSE
         CALL XABORT('NEIGHB: INVALID TYPE OF SYMMETRY.')
      ENDIF
*
      IF (IHEX.EQ.1) THEN
         CALL NEIGH1 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.2) THEN
         CALL NEIGH2 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.3) THEN
         CALL NEIGH3 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.4) THEN
         CALL NEIGH4 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.5) THEN
         CALL NEIGH5 (NC,J,K,N,POIDS)
         IF (-N.GT.NH) N=-N
      ELSE IF (IHEX.EQ.6) THEN
         CALL NEIGH6 (NC,J,K,N,POIDS)
         IF (-N.GT.NH) N=-N
      ELSE IF (IHEX.EQ.7) THEN
         CALL NEIGH7 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.8) THEN
         CALL NEIGH8 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.9) THEN
         CALL NEIGH9 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.10) THEN
         CALL NEIG10 (NC,J,K,N,POIDS)
      ELSE IF (IHEX.EQ.11) THEN
         CALL NEIG11 (NC,J,K,N,POIDS)
      ENDIF
      NEIGHB=N
      RETURN
      END
