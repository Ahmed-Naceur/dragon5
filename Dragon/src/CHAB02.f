*DECK CHAB02
      SUBROUTINE CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,GAR1,DELTA,FMULT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modify an array of cross section values.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NGRP    number of energy groups.
* IMOD    type of modification: =1: complete replacement; =2: replace
*         specific values by VALUE; =3: increase by VALUE; =4: multiply
*         by VALUE.
* VALUE   value used in modification operation.
* IGM     first energy group to modify.
* IGP     last energy group to modify.
* VAL     array of values used if IMOD=1.
*
*Parameters: input/output
* GAR1    cross section array to modify on input and
*         mofified cross section array at output.
*
*Parameters: output
* DELTA   difference in cross section value.
* FMULT   modification factors.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,IMOD,IGM,IGP
      REAL VALUE,VAL(NGRP),GAR1(NGRP),DELTA(NGRP),FMULT(NGRP)
*
      CALL XDRSET(DELTA,NGRP,0.0)
      CALL XDRSET(FMULT,NGRP,1.0)
      IF(IMOD.EQ.1) THEN
         DO 10 IG=IGM,IGP
         IF(GAR1(IG).EQ.0.0) THEN
            FMULT(IG)=1.0
         ELSE
            FMULT(IG)=VAL(IG)/GAR1(IG)
         ENDIF
         DELTA(IG)=VAL(IG)-GAR1(IG)
         GAR1(IG)=VAL(IG)
   10    CONTINUE
      ELSE IF(IMOD.EQ.2) THEN
         DO 20 IG=IGM,IGP
         IF(GAR1(IG).EQ.0.0) THEN
            FMULT(IG)=1.0
         ELSE
            FMULT(IG)=VALUE/GAR1(IG)
         ENDIF
         DELTA(IG)=VALUE-GAR1(IG)
         GAR1(IG)=VALUE
   20    CONTINUE
      ELSE IF(IMOD.EQ.3) THEN
         DO 30 IG=IGM,IGP
         IF(GAR1(IG).EQ.0.0) THEN
            FMULT(IG)=1.0
         ELSE
            FMULT(IG)=1.0+VALUE/GAR1(IG)
         ENDIF
         DELTA(IG)=VALUE
         GAR1(IG)=GAR1(IG)+VALUE
   30    CONTINUE
      ELSE IF(IMOD.EQ.4) THEN
         DO 40 IG=IGM,IGP
         FMULT(IG)=VALUE
         DELTA(IG)=GAR1(IG)*(VALUE-1.0)
         GAR1(IG)=GAR1(IG)*VALUE
   40    CONTINUE
      ENDIF
      RETURN
      END
