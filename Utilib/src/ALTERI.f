*DECK ALTERI
      SUBROUTINE ALTERI(LCUBIC,N,X,VAL0,VAL1,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* determination of the TERP integration components using the order 4
* Ceschino method with cubic Hermite polynomials.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* LCUBIC  =.TRUE.: cubic Ceschino interpolation; =.FALSE: linear
*         Lagrange interpolation.
* N       number of points.
* X       abscissas
* VAL0    left integration limit.
* VAL1    right integration limit.
*
*Parameters: output
* TERP    integration components.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL LCUBIC
      INTEGER N
      REAL X(N),VAL0,VAL1,TERP(N)
*----
*  LOCAL VARIABLES
*----
      REAL UU(2)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WK
*
      IF(N.LE.1) CALL XABORT('ALTERI: INVALID NUMBER OF POINTS.')
      IF(VAL1.LE.VAL0) CALL XABORT('ALTERI: INVALID LIMITS.')
      IF((VAL0.LT.X(1)).OR.(VAL1.GT.X(N))) CALL XABORT('ALTERI: UNABLE'
     1 //' TO INTEGRATE.')
      IF(N.EQ.2) GO TO 110
      DO 10 I=1,N
      TERP(I)=0.0
   10 CONTINUE
*----
*  LINEAR LAGRANGE POLYNOMIALS.
*----
      IF(.NOT.LCUBIC) THEN
         DO 15 I0=1,N-1
         IF((VAL0.LT.X(I0+1)).AND.(VAL1.GT.X(I0))) THEN
            A=MAX(VAL0,X(I0))
            B=MIN(VAL1,X(I0+1))
            DX=X(I0+1)-X(I0)
            TERP(I0)=TERP(I0)+(X(I0+1)-0.5*(A+B))*(B-A)/DX
            TERP(I0+1)=TERP(I0+1)+(0.5*(A+B)-X(I0))*(B-A)/DX
         ENDIF
   15    CONTINUE
         RETURN
      ENDIF
*----
*  CESCHINO CUBIC POLYNOMIALS.
*----
      ALLOCATE(WK(3,N))
      DO 16 I=1,N
      WK(3,I)=0.0
   16 CONTINUE
      DO 30 I0=1,N-1
      IF((VAL0.LT.X(I0+1)).AND.(VAL1.GT.X(I0))) THEN
         A=MAX(VAL0,X(I0))
         B=MIN(VAL1,X(I0+1))
         CC=0.5*(B-A)
         DX=X(I0+1)-X(I0)
         U1=(A-0.5*(X(I0)+X(I0+1)))/DX
         U2=(B-0.5*(X(I0)+X(I0+1)))/DX
         UU(1)=0.5*(-(U2-U1)*0.577350269189626+U1+U2)
         UU(2)=0.5*((U2-U1)*0.577350269189626+U1+U2)
         DO 20 JS=1,2
         H1=(3.0*(0.5-UU(JS))**2-2.0*(0.5-UU(JS))**3)*CC
         H2=((0.5-UU(JS))**2-(0.5-UU(JS))**3)*CC
         H3=(3.0*(0.5+UU(JS))**2-2.0*(0.5+UU(JS))**3)*CC
         H4=(-(0.5+UU(JS))**2+(0.5+UU(JS))**3)*CC
         TERP(I0)=TERP(I0)+H1
         TERP(I0+1)=TERP(I0+1)+H3
         WK(3,I0)=WK(3,I0)+H2*DX
         WK(3,I0+1)=WK(3,I0+1)+H4*DX
   20    CONTINUE
      ENDIF
   30 CONTINUE
*----
*  COMPUTE THE COEFFICIENT MATRIX.
*----
      HP=1.0/(X(2)-X(1))
      WK(1,1)=HP
      WK(2,1)=HP
      DO 40 I=2,N-1
      HM=HP
      HP=1.0/(X(I+1)-X(I))
      WK(1,I)=2.0*(HM+HP)
      WK(2,I)=HP
   40 CONTINUE
      WK(1,N)=HP
      WK(2,N)=HP
*----
*  FORWARD ELIMINATION.
*----
      PMX=WK(1,1)
      WK(3,1)=WK(3,1)/PMX
      DO 50 I=2,N
      GAR=WK(2,I-1)
      WK(2,I-1)=WK(2,I-1)/PMX
      PMX=WK(1,I)-GAR*WK(2,I-1)
      WK(3,I)=(WK(3,I)-GAR*WK(3,I-1))/PMX
   50 CONTINUE
*----
*  BACK SUBSTITUTION.
*----
      DO 60 I=N-1,1,-1
      WK(3,I)=WK(3,I)-WK(2,I)*WK(3,I+1)
   60 CONTINUE
*----
*  COMPUTE THE INTERPOLATION FACTORS.
*----
      TEST=1.0
      DO 100 J=1,N
      IMIN=MAX(2,J-1)
      IMAX=MIN(N-1,J+1)
      DO 70 I=1,N
      WK(1,I)=0.0
   70 CONTINUE
      WK(1,J)=1.0
      HP=1.0/(X(2)-X(1))
      YLAST=WK(1,IMIN-1)
      WK(1,IMIN-1)=2.0*HP*HP*(WK(1,IMIN)-WK(1,IMIN-1))
      DO 80 I=IMIN,IMAX
      HM=HP
      HP=1.0/(X(I+1)-X(I))
      PMX=3.0*(HM*HM*(WK(1,I)-YLAST)+HP*HP*(WK(1,I+1)-WK(1,I)))
      YLAST=WK(1,I)
      WK(1,I)=PMX
   80 CONTINUE
      WK(1,IMAX+1)=2.0*HP*HP*(WK(1,IMAX+1)-YLAST)
      DO 90 I=IMIN-1,IMAX+1
      TERP(J)=TERP(J)+WK(1,I)*WK(3,I)
   90 CONTINUE
      IF(ABS(TERP(J)).LE.1.0E-7) TERP(J)=0.0
      TEST=TEST-TERP(J)/(VAL1-VAL0)
  100 CONTINUE
      IF(ABS(TEST).GT.1.0E-5) CALL XABORT('ALTERI: WRONG TERP FACTORS.')
      DEALLOCATE(WK)
      RETURN
*
  110 TERP(1)=(X(2)-0.5*(VAL0+VAL1))*(VAL1-VAL0)/(X(2)-X(1))
      TERP(2)=(0.5*(VAL0+VAL1)-X(1))*(VAL1-VAL0)/(X(2)-X(1))
      RETURN
      END
