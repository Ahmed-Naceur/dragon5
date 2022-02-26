*DECK NXTQPS
      SUBROUTINE NXTQPS(NDIM  ,DANGLT,DNPDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate direction defining the planes normal to a solid angle.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
*  G. Marleau, R. Roy, M. Hampartzounian
*
*Parameters: input
* NDIM    number of dimensions for geometry.
* DANGLT  direction of track.
*
*Parameters: output
* DNPDIR  directions defining plane normal to track.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  Extracted from the subroutine XELEQN of EXCELL.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NDIM
      DOUBLE PRECISION DANGLT(3)
      DOUBLE PRECISION DNPDIR(3,2,3)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTQPS')
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Local variables
*----
      INTEGER          IPL
      DOUBLE PRECISION X,Y,Z,SUPX,SUPY,SUPZ,OOSUPX,OOSUPY,OOSUPZ,
     >                 XOSUPX,YOSUPY,ZOSUPZ
*----
*  Define reference position
*----
      X  = DANGLT(1)
      Y  = DANGLT(2)
      IF(NDIM .EQ. 2) THEN
        DNPDIR(1,1,1)=-Y
        DNPDIR(2,1,1)=X
      ELSE IF(NDIM .EQ. 3) THEN
        SUPX  = SQRT( DONE - X * X )
        SUPY  = SQRT( DONE - Y * Y )
        OOSUPX= DONE / SUPX
        OOSUPY= DONE / SUPY
        XOSUPX=  X  / SUPX
        YOSUPY=  Y  / SUPY
        Z  = DANGLT(3)
        SUPZ  = SQRT( DONE - Z * Z )
        OOSUPZ= DONE / SUPZ
        ZOSUPZ=  Z  / SUPZ
        DO IPL=1,2*NDIM-3
          IF(IPL .EQ. 1) THEN
            DNPDIR( 1, 1 ,IPL)= -Y * OOSUPZ
            DNPDIR( 2, 1 ,IPL)=  X * OOSUPZ
            DNPDIR( 3, 1 ,IPL)=  DZERO
            DNPDIR( 1, 2 ,IPL)=  X * ZOSUPZ
            DNPDIR( 2, 2 ,IPL)=  Y * ZOSUPZ
            DNPDIR( 3, 2 ,IPL)=      - SUPZ
          ELSE IF(IPL .EQ. 2) THEN
            DNPDIR( 1, 1 ,IPL)= -Z * OOSUPY
            DNPDIR( 2, 1 ,IPL)=  DZERO
            DNPDIR( 3, 1 ,IPL)=  X * OOSUPY
            DNPDIR( 1, 2 ,IPL)=  X * YOSUPY
            DNPDIR( 2, 2 ,IPL)=      - SUPY
            DNPDIR( 3, 2 ,IPL)=  Z * YOSUPY
          ELSE IF(IPL .EQ. 3) THEN
            DNPDIR( 1, 1 ,IPL)=  DZERO
            DNPDIR( 2, 1 ,IPL)= -Z * OOSUPX
            DNPDIR( 3, 1 ,IPL)=  Y * OOSUPX
            DNPDIR( 1, 2 ,IPL)=      - SUPX
            DNPDIR( 2, 2 ,IPL)=  Y * XOSUPX
            DNPDIR( 3, 2 ,IPL)=  Z * XOSUPX
          ENDIF
        ENDDO
      ENDIF
*----
*  Processing finished: return
*----
      RETURN
*----
*  Output formats
*----
      END
