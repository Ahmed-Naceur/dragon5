      FUNCTION           AXGTRN(ITRCUR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Associate to TURN number a DRAGON name.
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
*
*Parameters: output
* AXGTRN  DRAGON turn name. 
*
*-----------------------------------------------------------------------
*
      IMPLICIT           NONE
*----
*  Local parameters
*----
      INTEGER            MAXTUR
      CHARACTER          NAMSBR*6
      PARAMETER         (MAXTUR=12,NAMSBR='AXGTRN')
*----
*  Routine input and output variables 
*----     
      INTEGER            ITRCUR
      CHARACTER          AXGTRN*(*)
*----
*  local variables
*----
      CHARACTER*2        CTURN(2*MAXTUR)
      SAVE               CTURN
*----
*  DEFINITION OF TURNS
*----
      DATA CTURN        /' A',' B',' C',' D',' E',' F',' G',' H',
     >                   ' I',' J',' K',' L',
     >                   '-A','-B','-C','-D','-E','-F','-G','-H',
     >                   '-I','-J','-K','-L'/
      IF(ITRCUR .LE. 0 .OR. ITRCUR .GT. 2*MAXTUR) CALL XABORT(NAMSBR//
     >  ': INVALID TURN NUMBER')
      AXGTRN=CTURN(ITRCUR)
      RETURN
      END
