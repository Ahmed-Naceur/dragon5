*DECK SYB32C
      SUBROUTINE SYB32C (PPLUS,TAUP,POPI,M)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluation of the $E_i$ function in 2D geometry.
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
* TAUP   initial optical path.
* POPI   delta optical path.
* M      order of the Bickley function (equal to M+1).
*
*Parameters: output
* PPLUS   value of the difference.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL PPLUS,TAUP,POPI
      INTEGER M
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600)
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
*
      TAUQ=TAUP+POPI
      IF(TAUP.GE.XLIM3) THEN
        PPLUS=0.
      ELSE IF(TAUQ.GE.XLIM3) THEN
        PPLUS=TABKI(M+1,TAUP)
      ELSE IF(POPI.LE.0.002) THEN
        PPLUS=(TABKI(M,TAUP)+TABKI(M,TAUQ))*POPI*0.5
      ELSE IF(POPI.LT.0.004) THEN
        PQLUS=(TABKI(M,TAUP)+TABKI(M,TAUQ))*POPI*0.5
        PRLUS=TABKI(M+1,TAUP)-TABKI(M+1,TAUQ)
        FACT=500.0*(POPI-0.002)
        PPLUS=PRLUS*FACT+PQLUS*(1.0-FACT)
      ELSE
        PPLUS=TABKI(M+1,TAUP)-TABKI(M+1,TAUQ)
      ENDIF
      RETURN
      END
