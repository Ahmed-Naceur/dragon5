*DECK S2MGET
      SUBROUTINE S2MGET(IFIN,HNAME,IDX,LRMS,NGRP,XS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To read a record in the Matlab-formatted SERPENT output file.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IFIN    SERPENT unit number.
* HNAME   record name.
* IDX     instance index on the record (generally equal to the burnup
*         step).
* LRMS    standard deviation flag (.TRUE. if the standard deviation is
*         present).
* NGRP    number of energy groups.
*
*Parameters: output
* XS      record recovered from SERPENT file.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IFIN,IDX,NGRP
      LOGICAL LRMS
      CHARACTER*(*) HNAME
      REAL XS(NGRP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HLINE*50000,PREFIX*26,HSMG*131
*
      REWIND(IFIN)
      IF(LEN(HNAME).GT.26) CALL XABORT('LCMGET: PREFIX OVERFLOW.')
      PREFIX=HNAME
      N=0
      DO
        READ(IFIN,'(A)',END=10) HLINE
        IND1=INDEX(HLINE,PREFIX)
        IF(IND1.GT.0) N=N+1
        IF(N.EQ.IDX) GO TO 20
      ENDDO
   10 WRITE(HSMG,'(22HS2MGET: UNABLE TO FIND,I5,14H INSTANCES OF ,A,
     > 9H RECORD (,I5,8H FOUND).)') IDX,HNAME,N
      CALL XABORT(HSMG)
   20 IND1=47
      DO IGR=1,NGRP
        READ(HLINE(IND1:IND1+11),'(E12.0)') XS(IGR)
        IF(LRMS) THEN
          IND1=IND1+21
        ELSE
          IND1=IND1+13
        ENDIF
      ENDDO
      RETURN
      END
