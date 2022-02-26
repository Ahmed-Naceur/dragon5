      FUNCTION           AXGTRS(ITRCUR,ISYM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transform turn according to symmetry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* ITRCUR  turn number.
* ISYM    symmetry option.
*
*Parameters: output
* AXGTRS  turn after symmetry is applied. 
*
*-----------------------------------------------------------------------
*
      IMPLICIT           NONE
*----
*  Local parameters
*----
      INTEGER            MAXTUR,MAXS
      CHARACTER          NAMSBR*6
      PARAMETER         (MAXTUR=12,MAXS=3,NAMSBR='AXGTRS')
*----
*  Routine input and output variables 
*----     
      INTEGER            ITRCUR,ISYM
      INTEGER            AXGTRS
*----
*  Local variables
*----
      INTEGER            ITURN(2*MAXTUR,MAXS)
      SAVE               ITURN
*----
*  Definition of turns
*----
      DATA        ITURN /
*----
*                         SYMMETRY IN *X*
*----
     >                     5 ,  8 ,  7 ,  6 ,  1 ,  4 ,  3 ,  2 ,
     >                     0 ,  0 ,  0 ,  0 ,
     >                    17 , 20 , 19 , 18 , 13 , 16 , 15 , 14 ,
     >                     0 ,  0 ,  0 ,  0 ,
*----
*                         SYMMETRY IN *Y*
*----
     >                     7 ,  6 ,  5 ,  8 ,  3 ,  2 ,  1 ,  4 ,
     >                     0 ,  0 ,  0 ,  0 ,
     >                    19 , 18 , 17 , 20 , 15 , 14 , 13 , 16 ,
     >                     0 ,  0 ,  0 ,  0 ,
*----
*                         TSYMMETRY IN *X-Y*
*----
     >                     6 ,  5 ,  8 ,  7 ,  2 ,  1 ,  4 ,  3 ,
     >                     0 ,  0 ,  0 ,  0 ,
     >                    18 , 17 , 20 , 19 , 14 , 13 , 16 , 15 ,
     >                     0 ,  0 ,  0 ,  0 /   
      IF(ITRCUR .LE. 0 .OR. ITRCUR .GT. 2*MAXTUR) CALL XABORT(NAMSBR//
     >  ': INVALID TURN NUMBER')
      IF(ISYM .LE. 0 .OR. ISYM .GT. MAXS+1 ) CALL XABORT(NAMSBR//
     >  ': INVALID SYMMETRY')
      IF(ISYM .LE. MAXS) THEN 
        AXGTRS=ITURN(ITRCUR,ISYM) 
      ELSE
        AXGTRS=MOD(MAXTUR+ITRCUR,2*MAXTUR) 
      ENDIF
      RETURN
      END
