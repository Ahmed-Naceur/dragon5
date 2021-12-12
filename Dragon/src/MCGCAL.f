*DECK MCGCAL
      SUBROUTINE MCGCAL(N,NOMCEL,NREG,MCUW,MCUI,LMCU,LMXMCU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of connection matrices.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov
*
*Parameters: input
* N       number of segments on this track.
* NOMCEL  integer tracking elements.
* NREG    number of volumes.
* LMCU    dimension (used) of MCUW.
* LMXMCU  real dimension of MCUW MCUI.
*
*Parameters: input/output
* MCUW    cell connection matrix.
* MCUI    cell connection matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,NOMCEL(N),NREG,MCUW(LMXMCU),MCUI(LMXMCU),LMCU,LMXMCU
*
      CHARACTER HSMG*131
*
      DO 10 I=1,N
      ICEL=NOMCEL(I)
      IF (I.EQ.N) THEN
         ICEL1=-1
      ELSE
         ICEL1=NOMCEL(I+1)
      ENDIF
      IF(ICEL.EQ.ICEL1) THEN
         IF (ICEL1.GT.NREG) THEN
            ICEL1=-1
         ELSE
            GOTO 6
         ENDIF
      ENDIF
*     IS THERE AREADY AN ELEMENT IN MATRIX FOR CELL ICEL ?
      IF (MCUW(ICEL).NE.0) GOTO 5
*      NO :
      MCUW(ICEL)=ICEL1
      GOTO 6
*      YES :
    5 II=ICEL
      IF(MCUW(II).EQ.ICEL1) GOTO 6
      ICEL=MCUI(II)
      IF(ICEL.NE.0) GOTO 5
*     ADD NEW ELEMENT 
      LMCU=LMCU+1
      IF(LMCU.GT.LMXMCU) THEN
         WRITE(HSMG,'(46HMCGCAL: MEMORY OVERFLOW. INCREASE MCU. LMXMCU=
     1              ,I10,1H.)') LMXMCU
         CALL XABORT(HSMG)
      ENDIF
      MCUW(LMCU)=ICEL1
      MCUI(II)=LMCU
    6 CONTINUE
   10 CONTINUE
*
      RETURN
      END
