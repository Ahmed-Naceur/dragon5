*DECK ALDERV
      SUBROUTINE ALDERV(N,X,Y)
*
*-----------------------------------------------------------------------
*
*Purpose:
* numerical derivation of an array of values using the order 4 Ceschino
* method (compatible with cubic splines).
*
*Copyright:
* Copyright (C) 1981 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* N       number of points.
* X       abscissas.
* Y       ordinates.
*
*Parameters: output
* Y       derivatives.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N
      REAL X(N),Y(N)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, DIMENSION(:,:), ALLOCATABLE :: WK
*
      IF(N.LE.1) CALL XABORT('ALDERV: INVALID NUMBER OF POINTS.')
      IF(N.EQ.2) GO TO 40
      ALLOCATE(WK(2,N))
      HP=1.0/(X(2)-X(1))
      WK(1,1)=HP
      WK(2,1)=HP
      YLAST=Y(1)
      Y(1)=2.0*HP*HP*(Y(2)-Y(1))
      DO 10 I=2,N-1
      HM=HP
      HP=1.0/(X(I+1)-X(I))
      WK(1,I)=2.0*(HM+HP)
      WK(2,I)=HP
      PMX=3.0*(HM*HM*(Y(I)-YLAST)+HP*HP*(Y(I+1)-Y(I)))
      YLAST=Y(I)
      Y(I)=PMX
   10 CONTINUE
      HM=HP
      WK(1,N)=HM
      WK(2,N)=HM
      Y(N)=2.0*HM*HM*(Y(N)-YLAST)
*
* FORWARD ELIMINATION.
      PMX=WK(1,1)
      Y(1)=Y(1)/PMX
      DO 20 I=2,N
      GAR=WK(2,I-1)
      WK(2,I-1)=WK(2,I-1)/PMX
      PMX=WK(1,I)-GAR*WK(2,I-1)
      Y(I)=(Y(I)-GAR*Y(I-1))/PMX
   20 CONTINUE
*
* BACK SUBSTITUTION.
      DO 30 I=N-1,1,-1
      Y(I)=Y(I)-WK(2,I)*Y(I+1)
   30 CONTINUE
      DEALLOCATE(WK)
      RETURN
*
   40 Y(1)=(Y(2)-Y(1))/(X(2)-X(1))
      Y(2)=Y(1)
      RETURN
      END
