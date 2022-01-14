*DECK NXTIAA
      FUNCTION NXTIAA(POSANN ,POSPIN,VOLINT)
*
*----------
*
*Purpose:
* Compute the volume of intersection between
* a 2--D annular region and an annular pin.
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
* POSANN  spatial description of the annular region with
*         POSANN(0) the radius, POSANN(1) the $X$ position
*         of center and POSANN(2) the $Y$ position
*         of center.
* POSPIN  spatial description of the annular pin region with
*         POSPIN(0) the radius, POSPIN(1) the $X$ position
*         of center and POSPIN(2) the $Y$ position
*         of center.
*
*Parameters: output
* NXTIAA  type of intersection between annular region and
*         annular pin, where:
*         = 0 means that there is no intersection
*         between the two regions;
*         = 1 means that the annular region
*         is all located inside the annular pin;
*         = 2 means that the annular pin
*         is all located inside the annular region;
*         =-1 means that the intersection between
*         the annular region and the annular pin is partial.
* VOLINT  2-D volume of intersection (area) between annular region and
*         annular pin.
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
      INTEGER          NXTIAA
      DOUBLE PRECISION POSANN(0:2),POSPIN(0:2)
      DOUBLE PRECISION VOLINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTIAA')
      INTEGER          IPRINT
      PARAMETER       (IPRINT=100)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO
      PARAMETER       (DZERO=0.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          IFACE
      DOUBLE PRECISION PX,PY,PX2,PY2,RA2,RP2,DAP2,DAP,XINT,YINT,
     >                 ANGLEP,ANGLEA,SP,SA,VOLANN,VOLPIN,ACARG
      DOUBLE PRECISION DT1,DT2,DT3
*----
*  Initialize NXTIAA and VOLINT
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) (POSANN(IFACE),IFACE=0,2)
        WRITE(IOUT,6011) (POSPIN(IFACE),IFACE=0,2)
      ENDIF
      PI=XDRCST('Pi',' ')
      NXTIAA=0
      VOLINT=DZERO
*----
*  Find distance from center of annulus to pin center
*  $d=\sqrt(y^{2}+y^{2})$
*  annular region/pin intersection if $d<r_{p}+r_{a}$
*  where $r_{a}$ is annular region radius
*  and $r_{p}$ is pin radius
*  Annular region in pin if $r_{a}<r_{p}$ $d<r_{p}-r_{a}$
*  Pin in annular region if $r_{p}<r_{a}$ $d<r_{a}-r_{p}$
*----
      PX=POSANN(1)-POSPIN(1)
      PY=POSANN(2)-POSPIN(2)
      PX2=PX*PX
      PY2=PY*PY
      RA2=POSANN(0)*POSANN(0)
      RP2=POSPIN(0)*POSPIN(0)
      VOLANN=PI*RA2
      VOLPIN=PI*RP2
      DAP2=PX2+PY2
      DAP=SQRT(DAP2)
      IF(DAP .LE. ABS(POSANN(0)-POSPIN(0)) ) THEN
        IF(POSANN(0) .LE. POSPIN(0)) THEN
          VOLINT=VOLANN
          NXTIAA=1
        ELSE
          VOLINT=VOLPIN
          NXTIAA=2
        ENDIF
      ELSE IF(DAP .LT. POSANN(0)+POSPIN(0)) THEN
*----
*  Find intersection points $(x_{i},y_{i})$ assuming
*  the annular region is rotated in such a way that its center
*  is located at $(d,0)$.
*  Equation for annular region is :
*    $(x-d)^{2}+y^{2}=r^{2}_{a}$
*  Equation for pin region is :
*    $x^{2}+y^{2}=r^{2}_{p}$
*  Intersection point is :
*    $x_{i}=(d^{2}+r^{2}_{p}-r^{2}_{s})/(2*d)$
*    $y_{i}=\sqrt{r^{2}_{p}-x^{2}_{i}}$
*----
        XINT=(DAP2+RP2-RA2)/(2*DAP)
        YINT=SQRT(RP2-XINT*XINT)
*----
*  Find angular opening for pin and annulus
*  Pin surface to the right of $x=x_{i}$ plane
*  $S_{p}=r_{p}^{2}\arccos(x_{i}/r_{p})-x_{i}y_{i}$
*  Annular surface to the left of $x=x_{i}$ plane
*  $S_{a}=r_{a}^{2}\arccos((d-x_{i})/r_{a})-(d-x_{i})y_{i}$
*----
        ACARG=XINT/POSPIN(0)
        IF(ACARG .GE. 1.0D0) THEN
          ANGLEP=ACOS(1.0D0)
        ELSE IF(ACARG .LE. -1.0D0) THEN
          ANGLEP=ACOS(-1.0D0)
        ELSE
          ANGLEP=ACOS(ACARG)
        ENDIF
        ACARG=(DAP-XINT)/POSANN(0)
        IF(ACARG .GE. 1.0D0) THEN
          ANGLEA=ACOS(1.0D0)
        ELSE IF(ACARG .LE. -1.0D0) THEN
          ANGLEA=ACOS(-1.0D0)
        ELSE
          ANGLEA=ACOS(ACARG)
        ENDIF
        SP=RP2*ANGLEP-XINT*YINT
        SA=RA2*ANGLEA-(DAP-XINT)*YINT
        VOLINT=SP+SA
        DT1=ABS(VOLINT-VOLPIN)
        DT2=ABS(VOLINT-VOLANN)
        DT3=ABS(VOLINT)
        IF(DT1 .LT. DCUTOF) THEN
          VOLINT=VOLPIN
          NXTIAA=2
        ELSE IF(DT2 .LT. DCUTOF) THEN
          VOLINT=VOLANN
          NXTIAA=1
        ELSE IF(DT3 .GT. DCUTOF) THEN
          NXTIAA=-1
        ELSE
          VOLINT=DZERO
        ENDIF
      ENDIF
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6012) NAMSBR,NXTIAA,VOLINT
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('POSANN={',2(F20.10,','),F20.10,'};')
 6011 FORMAT('POSPIN={',2(F20.10,','),F20.10,'};')
 6012 FORMAT(A6,'={',I5,',',F20.10,'};')
      END
