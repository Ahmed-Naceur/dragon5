*DECK LIBA28
      SUBROUTINE LIBA28(X,XP,NXP,LL,FP,L,IPROX,I0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the interpolation weights FP for the Lagrange interpolation.
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
* X       abscissa.
* XP      tabulated abscissa.
* NXP     number of tabulated points.
* LL      requested interpolation order.
*
*Parameters: output
* FP      weights for Lagrange interpolation.
* L       interpolation limit (number of non-zero weights).
* IPROX   index of closest tabulated point.
* I0      number of leading zero weights.
*
*Comments:
*  Evaluation method.
*  F(X) = sum for I = 1 to L of F(I+I0)*FP(I)
*  for LL.le.NXP uses an LL-point "centered" Lagrange interpolation,
*  otherwise it uses a linear interpolation formula.
*  Attention: it is assumed that XP(I+1) > XP(I) for all I
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NXP,LL,L,IPROX,I0
      REAL X
      DOUBLE PRECISION XP(NXP),FP(LL)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION PP,XI,XX
*
      L=LL
      IF(L.LE.NXP) GO TO 20
*----
*  CONSTANT FUNCTION DEFINED BY A SINGLE POINT
*----
      IF(NXP.LE.1)THEN
       I0=0
       L=1
       IPROX=1
       FP(1)=1.0D0
       GO TO 100
      ENDIF
*----
*  NXP < L: SWITCH TO LINEAR INTERPOLATION
*----
      L=2
   20 L2=(L+1)/2
*----
*  LOCATE FIRST POINT TO THE RIGHT OF XX (PT IXP)
*----
      DO 30 IXP=1,NXP
       IF(XP(IXP).GT.X)GO TO 40
   30 CONTINUE
*----
*  X IS TO THE RIGHT OF EVERY POINT XP
*----
      IXP=1
      IPROX=NXP
      GO TO 60
*----
*  XP(IMIN) IS THE FIRST POINT FOR THE INTERPOLATION
*----
   40 IMIN=IXP-L2
      IPROX=IXP
      IF(IXP.GT.1)THEN
       IF((X-XP(IXP-1)).LT.(XP(IXP)-X))IPROX=IXP-1
      ENDIF
      IF(L.EQ.1)THEN
       FP(1)=1.0D0
       I0=IPROX-1
       GO TO 100
      ENDIF
      IF(IMIN.GE.1)GO TO 50
      IMIN=1
      IMAX=L
      GO TO 70
*----
*  XP(IMAX) IS THE LAST POINT FOR THE INTERPOLATION
*----
   50 IMAX=IMIN+L-1
      IF(IMAX.LE.NXP)GO TO 70
   60 IMAX=NXP
      IMIN=NXP-L+1
*----
*  CENTERED POLYNOMIAL INTERPOLATION OF DEGRE L
*----
   70 I0=IMIN-1
      XX=X
      DO 90 I=IMIN,IMAX
      PP=1.0D0
      XI=XP(I)
      DO 80 J=IMIN,IMAX
      IF(I.NE.J)PP=PP*((XX-XP(J))/(XI-XP(J)))
   80 CONTINUE
      FP(I+1-IMIN)=PP
   90 CONTINUE
  100 RETURN
      END
