*DECK LELCRN
      FUNCTION LELCRN( CENTEC, RAYONC, X, Y)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Decide if the crown intersect a rectangular mesh.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* CENTEC  coordinates of center.                   
* RAYONC  inner and outer radius**2 of the crown.  
* X       X of the square.                         
* Y       Y of the square.                         
*
*Parameters: output
* LELCRN  checking flag: =.true. if interaction exists and
*         =.false. otherwise.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
      LOGICAL  LELCRN
*
      DOUBLE PRECISION CENTEC(2), RAYONC(2), X(2), Y(2), R
      INTEGER NBEXT, NBINT, IX, IY
*
      NBEXT=0
      NBINT=0
      DO 20 IX=1, 2
      DO 10 IY=1, 2
         R= (X(IX)-CENTEC(1))*(X(IX)-CENTEC(1))
     >    + (Y(IY)-CENTEC(2))*(Y(IY)-CENTEC(2))
         IF( R.LE.RAYONC(1) ) NBINT= NBINT+1
         IF( R.GE.RAYONC(2) ) NBEXT= NBEXT+1
   10 CONTINUE
   20 CONTINUE
      IF( NBINT.EQ.4 )THEN
*
*        RECTANGLE IS CONTAINED INSIDE THE INTERNAL RADIUS
         LELCRN=.FALSE.
      ELSEIF( NBEXT.EQ.4 )THEN
         IF( Y(1).LT.CENTEC(2).AND.CENTEC(2).LT.Y(2) )THEN
            IF( CENTEC(1).LT.X(1) )THEN
               LELCRN= (X(1)-CENTEC(1))*(X(1)-CENTEC(1)).LT.RAYONC(2)
            ELSEIF( X(2).LT.CENTEC(1) )THEN
               LELCRN= (X(2)-CENTEC(1))*(X(2)-CENTEC(1)).LT.RAYONC(2)
            ELSE
               LELCRN=.TRUE.
            ENDIF
         ELSEIF( X(1).LT.CENTEC(1).AND.CENTEC(1).LT.X(2) )THEN
            IF( CENTEC(2).LT.Y(1) )THEN
               LELCRN= (Y(1)-CENTEC(2))*(Y(1)-CENTEC(2)).LT.RAYONC(2)
            ELSEIF( Y(2).LT.CENTEC(2) )THEN
               LELCRN= (Y(2)-CENTEC(2))*(Y(2)-CENTEC(2)).LT.RAYONC(2)
            ELSE
               LELCRN=.TRUE.
            ENDIF
         ELSE
            LELCRN=.FALSE.
         ENDIF
      ELSE
         LELCRN=.TRUE.
      ENDIF
*
      RETURN
      END
