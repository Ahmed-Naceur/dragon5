*DECK ALDFIT
      SUBROUTINE ALDFIT(N,MA,X,Y,W,PARAM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* performs linear least squares fitting to a polynomial of a specified
* order in one independent variable using the Forsythe method.
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* N       number of data points i.e., number of X,Y values.
* MA      integer specifying the order of the polynomial.
* X       array of values of indep. variable.
* Y       array of values of dependent variable.
* W       array of weights.
*
*Parameters: output
* PARAM   real array of coefficients of the fitted polynomial.
*         PARAM(I)=coeff. of X**I.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,MA
      DOUBLE PRECISION X(N),Y(N),W(N),PARAM(0:MA)
*----
*  ALLOCATABLE ARRAYS
*----
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GAMMA,POLY,PP
*
      IF(MA.GE.N) CALL XABORT('ALDFIT: UNDER-DETERMINED SYSTEM.')
      AA=0.0D0
      BB=0.0D0
      CC=0.0D0
      DO 10 I=1,N
      AA=AA+W(I)*X(I)
      BB=BB+W(I)*Y(I)
      CC=CC+W(I)
   10 CONTINUE
      PARAM(0)=BB/CC
      IF(MA.EQ.0) RETURN
      ALLOCATE(GAMMA(MA,3),POLY(N,0:2),PP(0:MA,0:2))
      CALL XDDSET(POLY(1,0),N,1.0D0)
      GAMMA(1,1)=AA/CC
      GAMMA(1,2)=0.0D0
      AA=0.0D0
      BB=0.0D0
      DO 20 I=1,N
      POLY(I,1)=X(I)-GAMMA(1,1)
      AA=AA+W(I)*POLY(I,1)*Y(I)
      BB=BB+W(I)*POLY(I,1)**2
   20 CONTINUE
      GAMMA(1,3)=AA/BB
      DO 50 J=2,MA
      AA=0.0D0
      BB=0.0D0
      CC=0.0D0
      DD=0.0D0
      DO 30 I=1,N
      AA=AA+W(I)*X(I)*POLY(I,MOD(J-1,3))**2
      BB=BB+W(I)*POLY(I,MOD(J-1,3))**2
      CC=CC+W(I)*X(I)*POLY(I,MOD(J-1,3))*POLY(I,MOD(J-2,3))
      DD=DD+W(I)*POLY(I,MOD(J-2,3))**2
   30 CONTINUE
      GAMMA(J,1)=AA/BB
      GAMMA(J,2)=CC/DD
      AA=0.0D0
      BB=0.0D0
      DO 40 I=1,N
      POLY(I,MOD(J,3))=(X(I)-GAMMA(J,1))*POLY(I,MOD(J-1,3))-GAMMA(J,2)*
     1 POLY(I,MOD(J-2,3))
      AA=AA+W(I)*POLY(I,MOD(J,3))*Y(I)
      BB=BB+W(I)*POLY(I,MOD(J,3))**2
   40 CONTINUE
      GAMMA(J,3)=AA/BB
   50 CONTINUE
*
      DO 60 I=1,MA
      PP(I,0)=0.0D0
      PARAM(I)=0.0D0
   60 CONTINUE
      PP(0,0)=1.0D0
      DO 90 J=1,MA
      DO 70 I=0,MA
      PP(I,MOD(J,3))=0.0D0
   70 CONTINUE
      DO 80 I=0,J
      IF(I.LT.J) PP(I+1,MOD(J,3))=PP(I,MOD(J-1,3))
      PP(I,MOD(J,3))=PP(I,MOD(J,3))-PP(I,MOD(J-1,3))*GAMMA(J,1)
      IF(J.GT.1) PP(I,MOD(J,3))=PP(I,MOD(J,3))-PP(I,MOD(J-2,3))*
     1 GAMMA(J,2)
      PARAM(I)=PARAM(I)+PP(I,MOD(J,3))*GAMMA(J,3)
   80 CONTINUE
   90 CONTINUE
      DEALLOCATE(GAMMA,POLY,PP)
      RETURN
      END
