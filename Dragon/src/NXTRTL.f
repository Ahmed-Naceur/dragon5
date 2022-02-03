*DECK NXTRTL
      SUBROUTINE NXTRTL(IPRINT,NDIM  ,ITRN  ,TRKORI,ANGLES,
     >                  TRKORR,ANGROT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Rotate tracking line according to reference turn.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* NDIM    dimensions of problem.
* ITRN    geometry original turn number.
* TRKORI  original track origin.
* ANGLES  original track direction.
*
*Parameters: output
* TRKORR  rotated geometry track origin.
* ANGROT  rotated geometry track direction.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,ITRN
      DOUBLE PRECISION TRKORI(NDIM),ANGLES(NDIM),
     >                 TRKORR(NDIM),ANGROT(NDIM)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTRTL')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IKT,IDIR
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6011) 'Initial starting point  ',
     >                   (TRKORI(IDIR),IDIR=1,NDIM)
        WRITE(IOUT,6011) 'Initial direction       ',
     >                   (ANGLES(IDIR),IDIR=1,NDIM)
      ENDIF
*----
*  Z axis reflection for 3-D problems
*----
      IKT=ITRN
      IF(NDIM .EQ. 3) THEN
        IF(ITRN .GT. 12 ) THEN
          IKT=IKT-12
          TRKORR(NDIM)=-TRKORI(NDIM)
          ANGROT(NDIM)=-ANGLES(NDIM)
        ELSE
          TRKORR(NDIM)=TRKORI(NDIM)
          ANGROT(NDIM)=ANGLES(NDIM)
        ENDIF
      ENDIF
      IF(IKT .EQ. 1) THEN
*----
* no turn in $X-Y$ plane
*----
        DO IDIR=1,2
          TRKORR(IDIR)=TRKORI(IDIR)
          ANGROT(IDIR)=ANGLES(IDIR)
        ENDDO
      ELSE IF(IKT .EQ. 2) THEN
*----
*  ROTATION OF -PI/2 OF GEOMETRY IMPLIES A ROTATION
*  OF PI/2 OF LINE.
*----
        TRKORR(1)=-TRKORI(2)
        TRKORR(2)= TRKORI(1)
        ANGROT(1)=-ANGLES(2)
        ANGROT(2)= ANGLES(1)
      ELSE IF(IKT .EQ. 3) THEN
*----
*  ROTATION OF PI OF GEOMETRY IMPLIES A ROTATION
*  OF -PI OF LINE.
*----
        TRKORR(1)=-TRKORI(1)
        TRKORR(2)=-TRKORI(2)
        ANGROT(1)=-ANGLES(1)
        ANGROT(2)=-ANGLES(2)
      ELSE IF(IKT .EQ. 4) THEN
*----
*  ROTATION OF -3*PI/2 OF GEOMETRY IMPLIES A ROTATION
*  OF 3PI/2 OF LINE.
*----
        TRKORR(1)= TRKORI(2)
        TRKORR(2)=-TRKORI(1)
        ANGROT(1)= ANGLES(2)
        ANGROT(2)=-ANGLES(1)
      ELSE IF(IKT .EQ. 5) THEN
*----
*  REFLECTION WITH RESPECT TO AXIS  // TO Y
*----
        TRKORR(1)=-TRKORI(1)
        TRKORR(2)= TRKORI(2)
        ANGROT(1)=-ANGLES(1)
        ANGROT(2)= ANGLES(2)
      ELSE IF(IKT .EQ. 6) THEN
*----
*  ROTATION OF PI/2 FOLLOWED BY
*  REFLECTION WITH RESPECT TO AXIS  // TO Y
*  IMPLIES REFLECTION WITH RESPECT TO AXIS  // TO Y
*  FOLLOWED BY A ROTATION OF -PI/2 OF LINE.
*----
        TRKORR(1)= TRKORI(2)
        TRKORR(2)= TRKORI(1)
        ANGROT(1)= ANGLES(2)
        ANGROT(2)= ANGLES(1)
      ELSE IF(IKT .EQ. 7) THEN
*----
*  REFLECTION WITH RESPECT TO AXIS // TO X
*----
        TRKORR(1)= TRKORI(1)
        TRKORR(2)=-TRKORI(2)
        ANGROT(1)= ANGLES(1)
        ANGROT(2)=-ANGLES(2)
      ELSE IF(IKT .EQ. 8) THEN
*----
*  ROTATION OF PI/2 FOLLOWED BY
*  REFLECTION WITH RESPECT TO AXIS // TO X
*  IMPLIES REFLECTION WITH RESPECT TO AXIS  // TO X
*  FOLLOWED BY A ROTATION OF -PI/2 OF LINE.
*----
        TRKORR(1)=-TRKORI(2)
        TRKORR(2)=-TRKORI(1)
        ANGROT(1)=-ANGLES(2)
        ANGROT(2)=-ANGLES(1)
      ENDIF
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6011) 'Final starting point    ',
     >                   (TRKORR(IDIR),IDIR=1,NDIM)
        WRITE(IOUT,6011) 'Final direction         ',
     >                   (ANGROT(IDIR),IDIR=1,NDIM)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6011 FORMAT(2X,A24,1P,3E20.12)
      END
