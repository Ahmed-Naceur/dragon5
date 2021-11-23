*DECK NXTIRR
      FUNCTION NXTIRR(XYCAR ,XYPIN ,VOLINT)
*
*----------
*
*Purpose:
* Compute the volume of intersection between
* a rectangular region and a Cartesian pin.
* centered at the origin.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* XYCAR   spatial description of the Cartesian region with:
*         XYCAR(1) for left face; XYCAR(2) for right face;
*         XYCAR(3) for bottom face, XYCAR(4) for top face
*         positions.
* XYPIN   spatial description of the Cartesian pin region with
*         XYPIN(1) for left face; XYPIN(2) for right face;
*         XYPIN(3) for bottom face; XYPIN(4) for top face
*         positions.
*
*Parameters: output
* NXTIRR  type of intersection between Cartesian region and
*         annular pin or annular region and Cartesian pin, where:
*         =0  means that there is no intersection
*         between the two regions;
*         = 1 means that the Cartesian region
*         is all located inside the Cartesian pin;
*         = 2 means that the Cartesian pin
*         is all located inside the Cartesian region;
*         =-1 means that the intersection between
*         the Cartesian region and the Cartesian pin is partial.
* VOLINT  2-D volume of intersection (area) between Cartesian region and
*         Cartesian pin.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NXTIRR
      DOUBLE PRECISION XYCAR(4),XYPIN(4)
      DOUBLE PRECISION VOLINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTIRR')
      INTEGER          IPRINT
      PARAMETER       (IPRINT=100)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO
      PARAMETER       (DZERO=0.0D0)
*----
*  Functions
*----
      INTEGER          NXTPRR,ITYPRR
*----
*  Local variables
*----
      INTEGER          IFACE
      DOUBLE PRECISION VOLCAR,VOLPIN,XYINT(4)
      DOUBLE PRECISION DT1,DT2,DT3
*----
*  Initialize NXTIRR and VOLINT
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) (XYCAR(IFACE),IFACE=1,4)
        WRITE(IOUT,6011) (XYPIN(IFACE),IFACE=1,4)
      ENDIF
      NXTIRR=0
      VOLINT=DZERO
      VOLCAR=(XYCAR(2)-XYCAR(1))*(XYCAR(4)-XYCAR(3))
      VOLPIN=(XYPIN(2)-XYPIN(1))*(XYPIN(4)-XYPIN(3))
*----
*  Find rectangle of intersection between the two rectangles.
*----
      ITYPRR=NXTPRR(XYCAR ,XYPIN ,XYINT )
*----
*  For cases with intersection, compute volume of intersection
*  and type of intersection
*----
      IF(ITYPRR .NE. 0) THEN
        VOLINT=(XYINT(2)-XYINT(1))*(XYINT(4)-XYINT(3))
        DT1=ABS(VOLINT-VOLPIN)
        DT2=ABS(VOLINT-VOLCAR)
        DT3=ABS(VOLINT)
        IF(DT1 .LT. DCUTOF) THEN
          VOLINT=VOLPIN
          NXTIRR=2
        ELSE IF(DT2 .LT. DCUTOF) THEN
          VOLINT=VOLCAR
          NXTIRR=1
        ELSE IF(DT3 .LT. DCUTOF) THEN
          VOLINT=DZERO
          NXTIRR=0
        ELSE
          NXTIRR=-1
        ENDIF
      ENDIF
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6012) NAMSBR,NXTIRR,VOLINT
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('XYCAR ={',3(F20.10,','),F20.10,'};')
 6011 FORMAT('XYPIN ={',3(F20.10,','),F20.10,'};')
 6012 FORMAT(A6,'={',I5,',',F20.10,'};')
      END
