*DECK MSRLUS
      SUBROUTINE MSRLUS2(LFORW,N,LC,LC0,IM,MCU,IM0,MCU0,JU,ILUDF,ILUCF,
     1                   CF,XIN,XOUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve LU y = x where LU is stored in the "L\U-I" form in MSR format.
* Can be use "in-place" i.e. XOUT=XIN.
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
* LC0     dimension of MCU0.
* IM      elements of row i in MCU(IM(i):IM(i+1)-1) for CF.
* MCU
* IM0     elements of row i in MCU0(IM0(i):IM0(i+1)-1) for ILUCF.
* MCU0
* JU      MCU(JU(i):IM(i+1)) corresponds to U.
*         MCU(IM(i)+1:JU(i)-1) correspond to L.
* ILUDF   diagonal elements of U (inversed diagonal).
* ILUCF   non-diagonal elements of U which differs from CF.
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
      INTEGER N,LC,LC0,IM(N+1),MCU(LC),IM0(N+1),MCU0(LC0),JU(N)
      REAL ILUDF(N),ILUCF(LC0),CF(LC)
      DOUBLE PRECISION XIN(N),XOUT(N)
      LOGICAL LFORW
*---
* LOCAL VARIABLES
*---
      INTEGER I,J,IJ,IK
      REAL ICFIJ
      IF(LFORW) THEN
*---
* FORWARD SOLVE
*---
        DO I=1,N
           XOUT(I)=XIN(I)
           DO IJ=IM(I)+1,JU(I)-1
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) THEN
                 DO IK=IM0(I)+1,IM0(I+1)
                    IF (MCU0(IK).EQ.J) THEN
                       ICFIJ=ILUCF(IK)
                       GOTO 10
                    ENDIF
                 ENDDO
                 ICFIJ=ILUDF(J)*CF(IJ)
 10              CONTINUE
                 XOUT(I)=XOUT(I)-ICFIJ*XOUT(J)
              ENDIF
           ENDDO
        ENDDO
*---
* BACKWARD SOLVE
*---
        DO I=N,1,-1
           DO IJ=JU(I),IM(I+1)
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) THEN
                 DO IK=IM0(I)+1,IM0(I+1)
                    IF (MCU0(IK).EQ.J) THEN
                       ICFIJ=ILUCF(IK)
                       GOTO 20
                    ENDIF
                 ENDDO
                 ICFIJ=CF(IJ)
 20              CONTINUE
                 XOUT(I)=XOUT(I)-ICFIJ*XOUT(J)
              ENDIF
           ENDDO
           XOUT(I)=ILUDF(I)*XOUT(I)
        ENDDO
      ELSE
*---
* FORWARD SOLVE
*---
        DO I=1,N
           XOUT(I)=XIN(I)
        ENDDO
        DO I=1,N
           DO IJ=JU(I),IM(I+1)
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) THEN
                 DO IK=IM0(I)+1,IM0(I+1)
                    IF (MCU0(IK).EQ.J) THEN
                       ICFIJ=ILUCF(IK)
                       GOTO 30
                    ENDIF
                 ENDDO
                 ICFIJ=ILUDF(I)*CF(IJ)
 30              CONTINUE
                 XOUT(J)=XOUT(J)-ICFIJ*XOUT(I)
              ENDIF
           ENDDO
        ENDDO
*---
* BACKWARD SOLVE
*---
        DO I=N,1,-1
           XOUT(I)=ILUDF(I)*XOUT(I)
           DO IJ=IM(I)+1,JU(I)-1
              J=MCU(IJ)
              IF ((J.GT.0).AND.(J.LE.N)) THEN
                 DO IK=IM0(I)+1,IM0(I+1)
                    IF (MCU0(IK).EQ.J) THEN
                       ICFIJ=ILUCF(IK)
                       GOTO 40
                    ENDIF
                 ENDDO
                 ICFIJ=CF(IJ)
 40              CONTINUE
                 XOUT(J)=XOUT(J)-ILUDF(J)*ICFIJ*XOUT(I)
              ENDIF
           ENDDO
        ENDDO
      ENDIF
*
      RETURN
      END
