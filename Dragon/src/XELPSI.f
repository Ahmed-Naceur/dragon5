*DECK XELPSI
      FUNCTION XELPSI(ITYP,RANN2,XYPOS,XYPOS2,SPXY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute intersection surface between part of annular region to the 
* left of X-plane and either the part of the annular region above 
* Y-plane or the part below the Y-plane.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* ITYP    type of calculation:
*         =1 above y-plane;
*         =2 below y-plane.
* RANN2   pin radius squared.
* XYPOS   cartesian plane location:
*         (1,1)  left   x-plane;
*         (1,2)  right  x-plane;
*         (2,1)  bottom y-plane;
*         (2,2)  top    y-plane.
* XYPOS2  cartesian mesh squared with same notation as for XYPOS.
* SPXY    annular surface outside of planes with same notation 
*         as for XYPOS.
*
*Parameters: output
* XELPSI  intersection surface.
*
*-----------------------------------------------------------------------
*
      IMPLICIT     NONE
      INTEGER      ITYP
      DOUBLE PRECISION XELPSI,RANN2,XYPOS(2,2),XYPOS2(2,2),SPXY(2,2)
*----
*  LOCAL PARAMETERS
*----
      DOUBLE PRECISION   PI,DZERO
      PARAMETER         (PI=3.14159265358979323846D0,DZERO=0.0D0)
      DOUBLE PRECISION   SQANN,YFC,XFC
*----
*  TEST IF POINT OF INTEREST IS LOCATED INSIDE
*  ANNULAR REGION
*----
      XELPSI=0.0D0
      IF(XYPOS2(2,ITYP)+XYPOS2(1,1).LT.RANN2) THEN
*----
*  FOR POINT INSIDE ANNULAR REGION
*  1) (SUM OF ANNULAR INTERSECTION SURFACES)/2
*     -INTERSECTION SURFACE
*     +(INTERNAL REGION CARTESIAN SURFACE)/4
*     =(ANNULAR SURFACE)/4
*----
        SQANN=0.25D0*PI*RANN2
        YFC=-XYPOS(2,ITYP)
        IF(ITYP.EQ.1) THEN
          XFC=-XYPOS(1,1)
        ELSE
          XFC=XYPOS(1,1)
        ENDIF
        XELPSI=0.5D0*(SPXY(1,1)+SPXY(2,ITYP))+XFC*YFC-SQANN
      ELSE
        IF(ITYP.EQ.1) THEN
          IF(XYPOS(2,ITYP).LT.DZERO.AND.XYPOS(1,1).LT.DZERO) THEN
            XELPSI=DZERO
          ELSE IF(XYPOS(2,ITYP).GT.DZERO.AND.XYPOS(1,1).GT.DZERO) THEN
            XELPSI=SPXY(2,ITYP)+SPXY(1,1)-PI*RANN2
          ELSE IF(XYPOS(2,ITYP).GT.DZERO.AND.XYPOS(1,1).LT.DZERO) THEN
            XELPSI=SPXY(1,1)
          ELSE IF(XYPOS(2,ITYP).LT.DZERO.AND.XYPOS(1,1).GT.DZERO) THEN
            XELPSI=SPXY(2,ITYP)
          ENDIF
        ELSE
          IF(XYPOS(2,ITYP).LT.DZERO.AND.XYPOS(1,1).LT.DZERO) THEN
            XELPSI=SPXY(1,1)
          ELSE IF(XYPOS(2,ITYP).GT.DZERO.AND.XYPOS(1,1).GT.DZERO) THEN
            XELPSI=SPXY(2,ITYP)
          ELSE IF(XYPOS(2,ITYP).GT.DZERO.AND.XYPOS(1,1).LT.DZERO) THEN
            XELPSI=DZERO
          ELSE IF(XYPOS(2,ITYP).LT.DZERO.AND.XYPOS(1,1).GT.DZERO) THEN
            XELPSI=SPXY(2,ITYP)+SPXY(1,1)-PI*RANN2
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
