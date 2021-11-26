*DECK PIJD2S
      SUBROUTINE PIJD2S(NREG,NSOUT,PROB,PROBKS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Charge PROBKS matrices in the DRAGON square format.
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
*
*Parameters: output
* PROBKS  square probability matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
C----
C VARIABLES
C----
      INTEGER          NREG,NSOUT,KPRB,IIU,IIL,II,JJ
      DOUBLE PRECISION PROB(((NSOUT+NREG+2)*(NSOUT+NREG+1))/2)
      REAL             PROBKS(NREG*NREG)
C----
C  STORE IN SQUARE FORMAT
C----
      KPRB=(NSOUT+1)*(NSOUT+2)/2+NSOUT+1
      DO 20 JJ=1,NREG
         IIU=JJ
         IIL=(JJ-1)*NREG+1
         DO 10 II=1,JJ
            KPRB=KPRB+1
            PROBKS(IIL)=REAL(PROB(KPRB))
            PROBKS(IIU)=PROBKS(IIL)
            IIU=JJ+II*NREG
            IIL=IIL+1
   10    CONTINUE
         KPRB=KPRB+NSOUT+1
   20 CONTINUE
C
      RETURN
      END
