*DECK MCGPRA
      SUBROUTINE MCGPRA(LFORW,OPT,PACA,FLOUT,NLONG,LC,IM,MCU,JU,DIAGM,
     1                  CM,ILUDF,ILUCF,DIAGF,XIN,XOUT,LC0,IM0,MCU0,CF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiply a matrix stored in MSR format with a vector
* and apply preconditioner (diagonal / ILUO stored in MSR format) 
* to this vector.
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
* OPT     1 product matrix-vector only;
*         2 preconditioner only;
*         3 both.
* PACA    type of preconditioner.
* FLOUT   flag for output:
*         .false. matrix-vector product to -preconditioned vector
*         to XOUT;
*         .true.  matrix-vector product to XOUT and
*         preconditioned vector to XIN.
* NLONG   size of the corrective system.
* LC      dimension of profiled matrices MCU and CM.
* IM      MSR indexes vector.
* MCU     MSR indexes vector.
* JU      used in ACA  acceleration for ilu0.
* DIAGM   diagonal of matrix to multiply by.
* CM      non-diagonal elements of matrix to multiply by.
* ILUDF   diagonal of ilu0 matrix.
* ILUCF   non-diagonal elements of ilu0 matrix.
* DIAGF   vector of diagonal preconditioning.
* LC0     used in ILU0-ACA acceleration.
* IM0     used in ILU0-ACA acceleration.
* MCU0    used in ILU0-ACA acceleration.
* CF      non-diagonal elements of the matrix
*         corresponding to the ILU0 decomposition ILUCF.
*
*Parameters: input/output
* XIN     input vector to be precondition.
*
*Parameters: output
* XOUT    output vector preconditioned.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* SUBROUTINE ARGUMENTS
*----
      INTEGER OPT,PACA,NLONG,LC,IM(NLONG),MCU(LC),JU(NLONG),LC0,IM0(*),
     1 MCU0(*)
      REAL DIAGM(NLONG),CM(LC),ILUDF(NLONG),ILUCF(LC),DIAGF(NLONG),
     1 CF(LC)
      DOUBLE PRECISION XIN(NLONG),XOUT(NLONG)
      LOGICAL LFORW,FLOUT
*---
* LOCAL VARIABLES
*---
      INTEGER I,J,IJ
      DOUBLE PRECISION FF
*
      IF(MOD(OPT,2).EQ.1) THEN
*---
* MATRIX VECTOR PRODUCT
*---
        IF(LFORW) THEN
*         direct calculation
          DO I=1,NLONG
             FF=XIN(I)*DIAGM(I)
             DO IJ=IM(I)+1,IM(I+1)
                J=MCU(IJ)
                IF((J.GT.0).AND.(J.LE.NLONG)) FF=FF+CM(IJ)*XIN(J)
             ENDDO
             XOUT(I)=FF
          ENDDO
        ELSE
*         adjoint calculation
          DO I=1,NLONG
             XOUT(I)=XIN(I)*DIAGM(I)
          ENDDO
          DO I=1,NLONG
             DO IJ=IM(I)+1,IM(I+1)
                J=MCU(IJ)
                IF((J.GT.0).AND.(J.LE.NLONG)) XOUT(J)=XOUT(J)
     1                                               +CM(IJ)*XIN(I)
             ENDDO
          ENDDO
        ENDIF
      ENDIF
*
      IF((OPT/2).EQ.1) THEN
*---
* APPLY PRECONDITIONER
*---
      IF(PACA.EQ.4) THEN
*     apply ILU0 preconditioner (optimized storage = no extra-storage)
         IF(FLOUT) THEN
            CALL MSRLUS1(LFORW,NLONG,LC,IM,MCU,JU,ILUDF,CF,XOUT,XIN)
         ELSE
            CALL MSRLUS1(LFORW,NLONG,LC,IM,MCU,JU,ILUDF,CF,XOUT,XOUT)
         ENDIF
      ELSE IF(PACA.EQ.3) THEN
*     apply ILU0 preconditioner (optimized storage)
         IF(FLOUT) THEN
            CALL MSRLUS2(LFORW,NLONG,LC,LC0,IM,MCU,IM0,MCU0,JU,ILUDF,
     1           ILUCF,CF,XOUT,XIN)
         ELSE
            CALL MSRLUS2(LFORW,NLONG,LC,LC0,IM,MCU,IM0,MCU0,JU,ILUDF,
     1           ILUCF,CF,XOUT,XOUT)
         ENDIF
      ELSEIF(PACA.EQ.2) THEN
*     apply ILU0 preconditioner (complete storage)
         IF(FLOUT) THEN
            CALL MSRLUS(LFORW,NLONG,LC,IM,MCU,JU,ILUDF,ILUCF,XOUT,XIN)
         ELSE
            CALL MSRLUS(LFORW,NLONG,LC,IM,MCU,JU,ILUDF,ILUCF,XOUT,XOUT)
         ENDIF
      ELSEIF(PACA.EQ.1) THEN
*     apply Diagonal preconditioner         
         IF(FLOUT) THEN
            DO I=1,NLONG
               XIN(I)=XOUT(I)/DIAGF(I)
            ENDDO
         ELSE
            DO I=1,NLONG
               XOUT(I)=XOUT(I)/DIAGF(I)
            ENDDO
         ENDIF
      ELSEIF(PACA.EQ.0) THEN
*     no preconditioner
         IF(FLOUT) THEN
            DO I=1,NLONG
               XIN(I)=XOUT(I)
            ENDDO
         ENDIF
      ENDIF
      ELSE
*---
* DO NOT APPLY PRECONDITIONER
*---
         IF(FLOUT) THEN
            DO I=1,NLONG
               XIN(I)=XOUT(I)
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
      END
