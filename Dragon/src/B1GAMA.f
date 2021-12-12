*DECK B1GAMA
      DOUBLE PRECISION FUNCTION B1GAMA(IAPROX,B2,SIG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the gamma function for a P1 or B1 calculation.
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
* IAPROX  type of beta function calculation:
*         =0: LKRD or RHS; =1: P0 or P1; =2: B0 or B1.
* B2      buckling.
* SIG     total macroscopic cross section.
*
*Parameters: output
* B1GAMA  value of the gamma function.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IAPROX
      DOUBLE PRECISION B2,SIG
      PARAMETER(B2LIM=0.5D0)
*
      SIG2=SIG*SIG
      B1GAMA=0.0D0
      IF(IAPROX.LE.1) THEN
*        P0 OR P1 APPROXIMATION
         B1GAMA=1.0D0
      ELSE IF(IAPROX.EQ.2) THEN
*        B0 OR B1 APPROXIMATION
         IF(B2.EQ.0.0) THEN
            B1GAMA=1.0D0
         ELSE IF(B2.GT.0.05D0*SIG2) THEN
            TMP=SQRT(B2)/SIG
            BBS=SIG/SQRT(B2)
            ATG=ATAN(TMP)
            B1GAMA=ATG/(3.0D0*BBS*(1.0D0-BBS*ATG))
         ELSE IF((B2.LE.0.05D0*SIG2).AND.(B2.GE.-0.05D0*SIG2)) THEN
            TMP=B2/SIG2
            B1GAMA=1.0D0+TMP*(4.0D0/15.0D0-12.0D0*TMP/175.0D0
     1      +92.0D0*TMP*TMP/2625.0D0)
         ELSE IF((B2.LT.-0.05D0*SIG2).AND.(B2.GE.-B2LIM*SIG2)) THEN
            TMP=SQRT(-B2)
            SB=SIG/TMP
            BLN=0.5D0*LOG((SIG+TMP)/(SIG-TMP))
            B1GAMA=BLN/(3.0D0*SB*(SB*BLN-1.0D0))
         ELSE IF(B2.LT.-B2LIM*SIG2) THEN
*           Pn-type fundamental mode extension for extreme subcritical
*           cases
            TMP2=SQRT(B2LIM*SIG2)
            ALPHA1=0.5D0*LOG((SIG+TMP2)/(SIG-TMP2))/TMP2
            ALPHA2=3.0D0*SIG/(3.0D0*SIG2-B2LIM*SIG2)
            ALPHA3=3.0D0*SIG/(3.0D0*SIG2+B2)
            B1BETA=(1.0D0-(ALPHA1-ALPHA2+ALPHA3)*SIG)/B2
            B1GAMA=(ALPHA1-ALPHA2+ALPHA3)/(3.0D0*SIG*B1BETA)
         ENDIF
      ELSE
         CALL XABORT('B1GAMA: INVALID VALUE OF IAPROX.')
      ENDIF
      RETURN
      END
