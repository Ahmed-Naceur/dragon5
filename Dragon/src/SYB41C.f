*DECK SYB41C
      SUBROUTINE SYB41C(PP,TAU0,TAUI,TAUJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $C_{ij}$ function in 1D spherical geometry.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* TAU0   side to side optical path.
* TAUI   optical path in volume i (or volume j).
* TAUJ   optical path in volume j (or volume i).
*
*Parameters: output
* PP   value of the expression.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL PP,TAU0,TAUI,TAUJ
*----
*  SYMMETRIC FORMULA IN (I,J). WE SET POPI <= POPJ
*----
      IF(TAUJ.LT.TAUI) THEN
         POPI = TAUJ
         POPJ = TAUI
      ELSE
         POPI = TAUI
         POPJ = TAUJ
      ENDIF
      IF(POPI.GT.0.004) THEN
         PP=EXP(-TAU0)-EXP(-(TAU0+POPI))-EXP(-(TAU0+POPJ))+
     1      EXP(-(TAU0+POPI+POPJ))
      ELSE IF((POPI.GT.0.002).AND.(POPJ.GT.0.004)) THEN
         PP=EXP(-TAU0)-EXP(-(TAU0+POPI))-EXP(-(TAU0+POPJ))+
     1      EXP(-(TAU0+POPI+POPJ))
         PQ=POPI*(EXP(-TAU0)-EXP(-(TAU0+POPJ)))
         FACTI=500.0*(POPI-0.002)
         PP=PP*FACTI+PQ*(1.0-FACTI)
      ELSE IF((POPI.GT.0.002).AND.(POPJ.GT.0.002)) THEN
         PP=EXP(-TAU0)-EXP(-(TAU0+POPI))-EXP(-(TAU0+POPJ))+
     1      EXP(-(TAU0+POPI+POPJ))
         PQ=POPI*(EXP(-TAU0)-EXP(-(TAU0+POPJ)))
         PR=POPJ*(EXP(-TAU0)-EXP(-(TAU0+POPI)))
         PS=POPI*POPJ*EXP(-TAU0)
         FACTI=500.0*(POPI-0.002)
         FACTJ=500.0*(POPJ-0.002)
         PP=PP*FACTI*FACTJ+PQ*(1.0-FACTI)*FACTJ+PR*FACTI*(1.0-FACTJ)
     1     +PS*(1.0-FACTI)*(1.0-FACTJ)
      ELSE IF(POPJ.GT.0.004) THEN
         PP=POPI*(EXP(-TAU0)-EXP(-(TAU0+POPJ)))
      ELSE IF(POPJ.GT.0.002) THEN
         PP=POPI*(EXP(-TAU0)-EXP(-(TAU0+POPJ)))
         PS=POPI*POPJ*EXP(-TAU0)
         FACTJ=500.0*(POPJ-0.002)
         PP=PP*FACTJ+PS*(1.0-FACTJ)
      ELSE
         PP=POPI*POPJ*EXP(-TAU0)
      ENDIF
      RETURN
      END
