*DECK XELTSW
      SUBROUTINE XELTSW(ABSC,NANGLE,PTSANG,WGTANG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To compute the integration weights for cyclic tracking.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* ABSC   multidimensional width of the cell.
* NANGLE number of angles.
* PTSANG integration points.
*
*Parameters: output
* WGTANG integration weights.
*
*Reference:
* R. Roy, G. Marleau, A. Hebert and D. Rozon,
* A Cyclic Tracking Procedure for Collision Probability Calculations
* in 2-D Lattices, Advances in Mathematics, Computations and 
* Reactor Physics, Pittsburgh, PA, April 28 - May 2 (1991).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      DOUBLE PRECISION ABSC(2)
      INTEGER          NANGLE
      DOUBLE PRECISION PTSANG(NANGLE)
      DOUBLE PRECISION WGTANG(NANGLE)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,MXANGL
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,MXANGL=30,NAMSBR='XELTSW')
      DOUBLE PRECISION XDRCST,PI
*----
*  LOCAL VARIABLES
*----
      INTEGER          I,K,N,P
      DOUBLE PRECISION X(0:MXANGL-1),W(0:MXANGL-1),ACCUW,A,B
      INTEGER          INEGW
      DOUBLE PRECISION AC0,AC1,AC2,AC3
*----
*  Test for validity of NANGLE
*----
      PI=XDRCST('Pi',' ')
      IF(NANGLE.GT. MXANGL) CALL XABORT(NAMSBR//
     >': Number of specular azimuthal points too large')
      A= ABSC(1)
      B= ABSC(2)
      P= NANGLE-1
      ACCUW= 1.D0
      DO 10 I= 0, P
         X(I)= PTSANG(NANGLE-I)*PTSANG(NANGLE-I)
         W(I)= ACCUW
         ACCUW= ACCUW*DBLE(2*I+1)/DBLE(2*(I+1))
   10 CONTINUE
      N= NANGLE-1
      DO 30 K= 0, N-1
         DO 20 I= N, K+1, -1
            W(I)=  W(I) - W(I-1) * X(K)
   20    CONTINUE
   30 CONTINUE
      DO 60 K= N-1, 0, -1
         DO 40 I= K+1, N
            W(I)= W(I)  / ( X(I) - X(I-K-1) )
   40    CONTINUE
         DO 50 I= K, N-1
            W(I)= W(I)  - W(I+1)
   50    CONTINUE
   60 CONTINUE
      INEGW=0
      DO 70 I= 0, P
         WGTANG(NANGLE-I)= W(I)
         IF(WGTANG(NANGLE-I) .LT. 0.0D0) INEGW=INEGW+1
   70 CONTINUE
*----
*  If some weights are negative write warning
*  and use Sanchez weighting
*  R. Sanchez, L. Mao, S. Santandrea
*  Nucl. Sci. Eng. 140, 23-50 (2002).
*----
      IF(INEGW .GT. 0) THEN
        WRITE(IOUT,7000) NAMSBR
        I=1
        AC0=PI
        AC1=ACOS(PTSANG(I))
        AC2=ACOS(PTSANG(I+1))
        WGTANG(I)=ABS(AC2-AC1)/AC0
        ACCUW=WGTANG(1)
        AC3=0
        DO 80 I=2,NANGLE-1
          AC3=ACOS(PTSANG(I+1))
          WGTANG(I)=ABS(AC3-AC1)/AC0
          ACCUW=ACCUW+WGTANG(I)
          AC1=AC2
          AC2=AC3
  80    CONTINUE
        I=NANGLE
        WGTANG(I)=ABS(AC3-AC1)/AC0
        ACCUW=ACCUW+WGTANG(I)
        DO 90 I=1,NANGLE
        WGTANG(I)=WGTANG(I)/ACCUW
  90    CONTINUE
      ENDIF
*----
*  Processing finished, return
*----
      RETURN
*----
*  FORMATS
*----
 7000 FORMAT(' ******  WARNING in : ',A6,' *****'/
     >10X,'Some of the integration weights are negative'/
     >10X,'this may result in invalid integration of CP'/
     >10X,'Use Sanchez instead of Roy weighting')
      END
