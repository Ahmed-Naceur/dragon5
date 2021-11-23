*DECK MSRLUS
      SUBROUTINE MSRLUS1(LFORW,N,LC,IM,MCU,JU,ILUDF,CF,XIN,XOUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve LU y = x where LU is stored in the "L\U-I" form in MSR format.
* Can be use "in-place" i.e. XOUT=XIN.
* Special case for which the non-diagonal elements of U 
* are the same as the original matrix.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* LFORW   flag set to .false. to transpose the coefficient matrix.
* N       order of the system.
* LC      dimension of MCU.
* IM      elements of row i in MCU(IM(i):IM(i+1)-1)
* MCU
* JU      MCU(JU(i):IM(i+1)) corresponds to U.
*         MCU(IM(i)+1:JU(i)-1) correspond to L.
* ILUDF   diagonal elements of U (inversed diagonal).
* CF      non-diagonal elements of the original matrix. 
* XIN     input vector x. 
*    
*Parameters: output
* XOUT    output vector y.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N,LC,IM(N+1),MCU(LC),JU(N)
      REAL ILUDF(N),CF(LC)
      DOUBLE PRECISION XIN(N),XOUT(N)
      LOGICAL LFORW
*---
* LOCAL VARIABLES
*---
      INTEGER I,J,IJ
      IF(LFORW) THEN
*---
* FORWARD SOLVE
*---
        DO I=1,N
           XOUT(I)=XIN(I)
           DO IJ=IM(I)+1,JU(I)-1
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) XOUT(I)=XOUT(I)
     1                                          -ILUDF(J)*CF(IJ)*XOUT(J)
           ENDDO
        ENDDO
*---
* BACKWARD SOLVE
*---
        DO I=N,1,-1
           XOUT(I)=ILUDF(I)*XOUT(I)
           DO IJ=JU(I),IM(I+1)
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) XOUT(I)=XOUT(I)
     1                                          -ILUDF(I)*CF(IJ)*XOUT(J)
           ENDDO
        ENDDO
      ELSE
*---
* FORWARD SOLVE (TRANSPOSED LINEAR SYSTEM)
*---
        DO I=1,N
           XOUT(I)=XIN(I)
        ENDDO
        DO I=1,N
           DO IJ=JU(I),IM(I+1)
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) XOUT(J)=XOUT(J)
     1                                          -ILUDF(I)*CF(IJ)*XOUT(I)
           ENDDO
        ENDDO
*---
* BACKWARD SOLVE (TRANSPOSED LINEAR SYSTEM)
*---
        DO I=N,1,-1
           XOUT(I)=ILUDF(I)*XOUT(I)
           DO IJ=IM(I)+1,JU(I)-1
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) XOUT(J)=XOUT(J)
     1                                          -ILUDF(J)*CF(IJ)*XOUT(I)
           ENDDO
        ENDDO
      ENDIF
      RETURN
*
      END
