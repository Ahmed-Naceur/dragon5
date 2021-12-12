*DECK NXTXYZ
      SUBROUTINE NXTXYZ(IPTRK ,IPRINT,NDIM  ,ITYPBC,MAXMSH,NUCELL,
     >                  ABSC,DGMESH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find global cell limits.
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
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* IPRINT  print level.
* NDIM    number of dimensions for geometry.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesian boundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* MAXMSH  maximum number of elements in mesh vector for
*         each directions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
*
*Parameters: output
* ABSC    cell width and upper limit.
* DGMESH  meshing vector for global geometry.
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
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPRINT,NDIM,ITYPBC,MAXMSH,NUCELL(3)
      DOUBLE PRECISION ABSC(3,2),DGMESH(-1:MAXMSH,4)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTXYZ')
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Local variables
*----
      INTEGER          IDIR,ICELL
      CHARACTER        NAMREC*12
      DOUBLE PRECISION SIDEH,CENTH,DHMAX
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 20) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL XDDSET(DGMESH,(MAXMSH+2)*4,DZERO)
      IF(ITYPBC .EQ. 0) THEN
        DO IDIR=1,NDIM
          NAMREC='G00000001SM'//CDIR(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          ABSC(IDIR,2)=DGMESH(NUCELL(IDIR),IDIR)
          ABSC(IDIR,1)=ABSC(IDIR,2)-DGMESH(0,IDIR)
          IF(IPRINT .GE. 20) THEN
            WRITE(IOUT,6010) CDIR(IDIR),ABSC(IDIR,1)
          ENDIF
        ENDDO
        DO IDIR=NDIM+1,3
          ABSC(IDIR,1)=DONE
          ABSC(IDIR,2)=DONE
        ENDDO
      ELSE IF(ITYPBC .EQ. 1) THEN
*----
*  Find Cartesian box surrounding circle in plane
*----
        IDIR=4
        NAMREC='G00000001SM'//CDIR(IDIR)
        CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
        ABSC(1,2)=DGMESH(NUCELL(IDIR),IDIR)
        ABSC(1,1)=ABSC(1,2)
        ABSC(2,1)=ABSC(1,1)
        ABSC(2,2)=ABSC(1,2)
        IF(IPRINT .GE. 20) THEN
          WRITE(IOUT,6010) CDIR(IDIR),ABSC(IDIR,1)
        ENDIF
        IDIR=3
        IF(NDIM .EQ. 3) THEN
          NAMREC='G00000001SM'//CDIR(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          ABSC(IDIR,2)=DGMESH(NUCELL(IDIR),IDIR)
          ABSC(IDIR,1)=ABSC(IDIR,2)-DGMESH(0,IDIR)
          IF(IPRINT .GE. 20) THEN
            WRITE(IOUT,6010) CDIR(IDIR),ABSC(IDIR,1)
          ENDIF
        ELSE
          ABSC(IDIR,1)=DONE
          ABSC(IDIR,2)=DONE
        ENDIF
      ELSE
*----
*  Find Cartesian box surrounding hexagons in plane
*----
        DO IDIR=1,2
          NAMREC='G00000001SM'//CDIR(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          SIDEH=DGMESH(0,IDIR)
          CENTH=DGMESH(1,IDIR)
          DHMAX=DZERO
          DO ICELL=2,NUCELL(IDIR)
            DHMAX=MAX(DHMAX,ABS(DGMESH(NUCELL(IDIR),IDIR)-CENTH))
          ENDDO
          ABSC(IDIR,2)=DHMAX+SIDEH
          ABSC(IDIR,1)=2.0*ABSC(IDIR,2)
        ENDDO
        IDIR=3
        IF(NDIM .EQ. 3) THEN
          NAMREC='G00000001SM'//CDIR(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          ABSC(IDIR,2)=DGMESH(NUCELL(IDIR),IDIR)
          ABSC(IDIR,1)=ABSC(IDIR,2)-DGMESH(0,IDIR)
          IF(IPRINT .GE. 20) THEN
            WRITE(IOUT,6010) CDIR(IDIR),ABSC(IDIR,1)
          ENDIF
        ELSE
          ABSC(IDIR,1)=DONE
          ABSC(IDIR,2)=DONE
        ENDIF
      ENDIF
      IF(IPRINT .GE. 20) THEN
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
 6010 FORMAT(' Geometry width in ',A1,' = ',F20.15)
      END
