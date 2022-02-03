*DECK SYB31C
      SUBROUTINE SYB31C (PPLUS,TAUP,XOPJ,XOPI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $C_{ij}$ function in 1D cylindrical and 2D geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* TAUP   side to side optical path.
* XOPJ   optical path in volume j (or volume i).
* XOPI   optical path in volume i (or volume j).
*
*Parameters: output
* PPLUS   value of the probability.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL PPLUS,TAUP,XOPJ,XOPI
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600)
      REAL T(2,2),B(2,2)
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
*----
*  ASYMPTOTIC VALUE
*----
      IF(TAUP.GE.XLIM3) THEN
        PPLUS = 0.
        RETURN
      ENDIF
*----
*  SYMMETRIC FORMULA IN (I,J). WE SET POPI <= POPJ
*----
      IF(XOPJ.LT.XOPI) THEN
        POPI = XOPJ
        POPJ = XOPI
      ELSE
        POPI = XOPI
        POPJ = XOPJ
      ENDIF
*----
*  GENERAL CASE
*
*  POPI < POPJ < 1.
*  VOID (I,J) UP TO THE END OF PROGRAM
*----
      T(2,2) = TAUP + POPJ + POPI
      T(1,2) = TAUP + POPJ
      T(2,1) = TAUP + POPI
      T(1,1) = TAUP
*
      DO 15 I=1,2
      DO 10 J=1,2
      B(I,J) = TABKI(3,T(I,J))
   10 CONTINUE
   15 CONTINUE
*----
*  GENERAL DIFFERENCE FORMULA. THIS FORMULAS SHOULD NOT BE APPLIED TO
*  VOIDED (I,J) VOLUMES
*----
      IF(POPI.GE.0.004) THEN
*       LARGE POPI AND POPJ (MOST GENERAL CASE)
        PNLUS = B(2,2) + B(1,1) - (B(1,2) + B(2,1))
        PPLUS = PNLUS
      ELSE IF(POPJ.GT.0.002) THEN
*       SMALL POPI, LARGE POPJ. USE DERIVATIVE DIFFERENCES.
        IF(TAUP.GE.XLIM3) THEN
          APLUS=0.
        ELSE IF(TAUP+POPI.GE.XLIM3) THEN
          APLUS=TABKI(3,TAUP)
        ELSE IF(POPI.LE.0.002) THEN
          APLUS=(TABKI(2,TAUP)+TABKI(2,TAUP+POPI))*POPI*0.5
        ELSE IF(POPI.LT.0.004) THEN
          PQLUS=(TABKI(2,TAUP)+TABKI(2,TAUP+POPI))*POPI*0.5
          PRLUS=TABKI(3,TAUP)-TABKI(3,TAUP+POPI)
          FACT=500.0*(POPI-0.002)
          APLUS=PRLUS*FACT+PQLUS*(1.0-FACT)
        ELSE
          APLUS=TABKI(3,TAUP)-TABKI(3,TAUP+POPI)
        ENDIF
        IF(TAUP+POPJ.GE.XLIM3) THEN
          BPLUS=0.
        ELSE IF(TAUP+POPI+POPJ.GE.XLIM3) THEN
          BPLUS=TABKI(3,TAUP+POPJ)
        ELSE IF(POPI.LE.0.002) THEN
          BPLUS=(TABKI(2,TAUP+POPJ)+TABKI(2,TAUP+POPI+POPJ))*POPI*0.5
        ELSE IF(POPI.LT.0.004) THEN
          PQLUS=(TABKI(2,TAUP+POPJ)+TABKI(2,TAUP+POPI+POPJ))*POPI*0.5
          PRLUS=TABKI(3,TAUP+POPJ)-TABKI(3,TAUP+POPI+POPJ)
          FACT=500.0*(POPI-0.002)
          BPLUS=PRLUS*FACT+PQLUS*(1.0-FACT)
        ELSE
          BPLUS=TABKI(3,TAUP+POPJ)-TABKI(3,TAUP+POPI+POPJ)
        ENDIF
        PPLUS = APLUS - BPLUS
      ELSE IF(T(2,2).GE.XLIM3) THEN
*       SIMILAR TO A SECOND DERIVATIVE. ASYMPTOTIC TAUP.
        PNLUS = B(2,2) + B(1,1) - (B(1,2) + B(2,1))
        PPLUS = PNLUS
      ELSE
*       SIMILAR TO A SECOND DERIVATIVE.
        IF(T(2,2).EQ.0.) THEN
          PPLUS = 0.
        ELSE IF(TAUP.EQ.0.) THEN
          PPLUS = 1.57079632679489 + TABKI(1,TAUP+POPI+POPJ)
          PPLUS = PPLUS * 0.5 * POPI*POPJ
        ELSE IF(TAUP.LE.0.002) THEN
          PPLUS = TABKI(1,TAUP) + TABKI(1,TAUP+POPI+POPJ)
          PPLUS = PPLUS * 0.5 * POPI*POPJ
        ELSE IF(TAUP.LT.0.004) THEN
          PQLUS = TABKI(1,TAUP) + TABKI(1,TAUP+POPI+POPJ)
          PRLUS = TABKI(1,TAUP)
          TAUP1 = TAUP + POPI + POPJ
          PRLUS = PRLUS + TABKI(1,TAUP1)
          FACT=500.0*(TAUP-0.002)
          PPLUS = PRLUS * FACT + PQLUS * (1.0-FACT)
          PPLUS = PPLUS * 0.5 * POPI*POPJ
        ELSE
*         VOIDED VOLUMES (I,J) SEPARATED WITH NON-VOIDED VOLUMES.
          PRLUS = TABKI(1,TAUP)
          TAUP1 = TAUP + POPI + POPJ
          PRLUS = PRLUS + TABKI(1,TAUP1)
          PPLUS = PRLUS * 0.5 * POPI*POPJ
        ENDIF
      ENDIF
      RETURN
      END
