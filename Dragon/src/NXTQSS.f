*DECK NXTQSS
      SUBROUTINE NXTQSS(IPRINT,NDIM  ,ITYPBC,MAXMSH,NUCELL,DENUSR,
     >                  DGMESH,NPLANE,NPOINT,DENLIN,SPACLN,
     >                  WEIGHT,RADIUS,CENTER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define standard spatial quadrature.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
*  G. Marleau, R.Roy
*
*Parameters: input
* IPRINT  print level.
* NDIM    number of dimensions for geometry.
* ITYPBC  type of boundary conditions (
*         +/- 2 for hexagones;
*         +/- 1 for annular;
*         0     for Cartesian).
* MAXMSH  maximum number of elements in mesh vector for
*         each directions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* DENUSR  user defined track density.
* DGMESH  meshing vector for global geometry.
*
*Parameters: output
* NPLANE  number of normal planes considered.
* NPOINT  number of integration points along each axis
*         in a plane mormal to track direction.
* DENLIN  effective track density in plane.
* SPACLN  linear track spacing in the plane.
* WEIGHT  weight associated with each line in the plane.
* RADIUS  radius of circle (2-D) or sphere (3-D) surrounding
*         the geometry.
* CENTER  center of circle (2-D) or sphere (3-D) surrounding
*         the geometry.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  Extracted from the subroutine XELTI2 and XELTI3.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,ITYPBC,MAXMSH,NUCELL(3)
      DOUBLE PRECISION DENUSR,DGMESH(-1:MAXMSH,3)
      INTEGER          NPLANE,NPOINT
      DOUBLE PRECISION DENLIN,SPACLN,WEIGHT
      DOUBLE PRECISION RADIUS,CENTER(NDIM)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTQSS')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IDIR,LSTCEL,NPO2,IX
      DOUBLE PRECISION DM,XD
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Find the radius of the sphere (3-D) or circle surrounding
*  the cell.
*  Also find the true center of the cell
*----
      IF(NDIM .EQ. 3) THEN
        NPLANE=3
        DENLIN=SQRT(DENUSR)
      ELSE IF(NDIM .EQ. 2) THEN
        NPLANE=1
        DENLIN=DENUSR
      ENDIF
      IF(ABS(ITYPBC) .EQ. 2) THEN
*----
*  HEXAGONAL
*----
        RADIUS=DZERO
        DO IDIR=1,2
          LSTCEL=NUCELL(IDIR)
          CENTER(IDIR)=DGMESH(1,IDIR)
          DM=DZERO
          DO IX=1,LSTCEL
            XD=DGMESH(IX,IDIR)-DGMESH(1,IDIR)+DGMESH(0,IDIR)
            DM=MAX(DM,XD*XD)
          ENDDO
          RADIUS=RADIUS+DM
        ENDDO
        DM=DZERO
        DO IDIR=3,NDIM
          LSTCEL=NUCELL(IDIR)
          CENTER(IDIR)=(DGMESH(LSTCEL,IDIR)+DGMESH(0,IDIR))/DTWO
          DM=((DGMESH(LSTCEL,IDIR)-DGMESH(0,IDIR))/DTWO)**2
        ENDDO
        RADIUS=SQRT(RADIUS+DM)
      ELSE IF(ABS(ITYPBC) .EQ. 1) THEN
*----
*  RADIAL
*----
        RADIUS=DZERO
        DO IDIR=1,NDIM
          LSTCEL=NUCELL(IDIR)
          CENTER(IDIR)=(DGMESH(LSTCEL,IDIR)+DGMESH(0,IDIR))/DTWO
          RADIUS=RADIUS+(DGMESH(LSTCEL,IDIR)-DGMESH(0,IDIR))**2
        ENDDO
        RADIUS=SQRT(RADIUS)/DTWO
      ELSE
*----
*  CARTESIAN
*----
        RADIUS=DZERO
        DO IDIR=1,NDIM
          LSTCEL=NUCELL(IDIR)
          CENTER(IDIR)=(DGMESH(LSTCEL,IDIR)+DGMESH(0,IDIR))/DTWO
          RADIUS=RADIUS+(DGMESH(LSTCEL,IDIR)-DGMESH(0,IDIR))**2
        ENDDO
        RADIUS=SQRT(RADIUS)/DTWO
      ENDIF
      NPOINT=INT(DTWO*RADIUS*DENLIN)+1
      NPO2=NPOINT/2
      NPOINT=2*NPO2+1
      SPACLN=DTWO*RADIUS/DBLE(NPOINT)
      DENLIN=DONE/SPACLN
      IF(NDIM .EQ. 3) DENLIN=DENLIN*DENLIN
      WEIGHT=DONE/(DBLE(NPLANE)*DENLIN)
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6002) DENUSR,DENLIN,WEIGHT,NPOINT,SPACLN
        WRITE(IOUT,6010) RADIUS,(CENTER(IDIR),IDIR=1,NDIM)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Processing finished: return
*----
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(' Spatial tracking information :'/
     >       ' Minimum tracking density requested =',1P,E15.7/
     >       ' Effective tracking density selected=',E15.7/
     >       ' Effective tracking weight          =',E15.7/
     >       ' Number of points per direction     =',I15/
     >       ' Linear track spacing               =',E15.7)
 6010 FORMAT(' RADIUS = ',1P,E20.12/
     >       ' CENTER = ',3E20.12)
      END
