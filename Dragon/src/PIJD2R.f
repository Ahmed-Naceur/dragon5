*DECK PIJD2R
      SUBROUTINE PIJD2R(NREG,NSOUT,PROB,FACTOR,LPIJK,NELPIJ,N2PROB,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Charge PIJ matrices in the DRAGON symmetrized format.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NREG    number of zones for geometry.
* NSOUT   number of surfaces for geometry. 
* PROB    collision probabilities.
* FACTOR  one over total xs.
* LPIJK   PIJK flag.
* NELPIJ  number of terms in PIJ.
* N2PROB  number of terms in PROB.
*
*Parameters: output
* PIJ     symmetric probability matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* VARIABLES
*----
      INTEGER          NREG,NSOUT,NELPIJ,N2PROB,IUN,JUN,KPRB,IVV
      DOUBLE PRECISION PROB(N2PROB)
      REAL             FACTOR(NREG),PIJ(NELPIJ),COEF
      LOGICAL          LPIJK
*----
*  STORE IN SYMMETRIC FORMAT
*----
      IVV=0
      COEF=1.0
      IF(LPIJK) COEF=1.5
      KPRB=(NSOUT+1)*(NSOUT+2)/2+NSOUT+1
      DO 20 IUN=1,NREG
         DO 10 JUN=1,IUN
            KPRB=KPRB+1
            IVV=IVV+1
            PIJ(IVV)=COEF*REAL(PROB(KPRB))*FACTOR(IUN)*FACTOR(JUN)
   10    CONTINUE
         KPRB=KPRB+NSOUT+1
   20 CONTINUE
*
      RETURN
      END
