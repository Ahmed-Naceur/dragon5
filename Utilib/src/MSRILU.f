*DECK MSRILU
      SUBROUTINE MSRILU(N,LC,IM,MCU,JU,DIAGF,CF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* In-place ILU0 decomposition of a matrix store in MSR format.
* The resulting matrix is in "L\U-I" format with inversed diagonal elements
* It is assumed that each line of the matrix is sorted by increasing
* column index.
*
*Reference
* Y. Saad, "Iterative Methods for Sparse Linear Systems",
* PWS Publishing Company, Boston, 1996
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
* N       size of the matrix.
* LC      dimension of IM.
* IM      elements of row i in MCU(IM(i):IM(i+1)-1)
* MCU
* JU      MCU(JU(i):IM(i+1)) corresponds to U.
*         MCU(IM(i)+1:JU(i)-1) correspond to L.
*
* Parameters: input/output
* DIAGF   diagonal elements of the matrix to process / U (inversed diagonal).
* CF      non-diagonal elements of the matrix to process / L and U.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER N,LC,IM(N+1),MCU(LC),JU(N)
      DOUBLE PRECISION DIAGF(N),CF(LC)
*---
* LOCAL VARIABLES
*---
      DOUBLE PRECISION EPS
      PARAMETER(EPS=1D-7)
      INTEGER I,K,J,H,IK,IJ,KJ,IH
      DOUBLE PRECISION E
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IW
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IW(N))
*
      CALL XDISET(IW,N,0)
*---
*  LOOP OVER LINE INDEX I
*---
      DO I=1,N
*     store component H in relation with I
         DO IH=IM(I)+1,IM(I+1)
            H=MCU(IH)
            IF (H.GT.0) IW(H)=IH
         ENDDO
*       ---
*        LOOP OVER LINE K < I IN RELATION WITH I
*       ---
         DO IK=IM(I)+1,JU(I)-1
            K=MCU(IK)
            IF (K.GT.0) THEN
*              entry (I>K,K)
               E=DIAGF(K)*CF(IK)
               CF(IK)=E
*             ---
*              LOOP OVER COLUMN J > K IN RELATION WITH K
*             ---
               DO KJ=JU(K),IM(K+1)
                  J=MCU(KJ)
*                 entry (I>K,J>K)
                  IF (I.EQ.J) THEN ! diagonal term
                     DIAGF(J)=DIAGF(J)-E*CF(KJ)
                  ELSEIF(IW(J).GT.0) THEN ! if J is in relation with I
                     IJ=IW(J)
                     CF(IJ)=CF(IJ)-E*CF(KJ)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         IF (ABS(DIAGF(I)).GT.EPS) THEN
            DIAGF(I)=1.0/DIAGF(I)
         ELSE
            WRITE(6,*) 'I=',I,' PIVOT=',DIAGF(I)
            CALL PRINDM('DIAGF    ',DIAGF,N)
            CALL PRINDM('CF       ',CF,LC)
            CALL XABORT('MSRILU: ZERO PIVOT')
         ENDIF
*        reset IW
         DO IH=IM(I)+1,IM(I+1)
            H=MCU(IH)
            IF (H.GT.0) THEN
               IW(H)=0
            ENDIF
         ENDDO
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IW)
      RETURN
      END
