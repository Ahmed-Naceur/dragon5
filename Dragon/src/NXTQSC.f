*DECK NXTQSC
      SUBROUTINE NXTQSC(IPRINT,NDIM  ,NBANGL,MAXMSH,NUCELL,
     >                  DGMESH,DANGLT,DDENWT,DNSANG,NBSANG,DEPART)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define spatial quadrature for cyclic tracking.
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
* NBANGL  number of angles.
* MAXMSH  maximum number of elements in mesh vector for
*         each directions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* DGMESH  meshing vector for global geometry.
* DANGLT  director cosines of angles.
* DDENWT  angular density for each angle.
*
*Parameters: input/output
* DNSANG  spatial density required.
* NBSANG  number of segments for each angles.
*
*Parameters: output
* DEPART  track starting point.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  Extracted from the subroutine XELTS2.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,NBANGL,MAXMSH,NUCELL(3)
      DOUBLE PRECISION DGMESH(-1:MAXMSH,3),DANGLT(NDIM,NBANGL),
     >                 DDENWT(NBANGL),DNSANG(NBANGL)
      INTEGER          NBSANG(5,NBANGL)
      DOUBLE PRECISION DEPART(NDIM,2,NBANGL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTQSC')
      DOUBLE PRECISION DZERO,DONE,DTWO,DHALF
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0,
     >                 DHALF=DONE/DTWO)
*----
*  Local variables
*----
      INTEGER          IANG,IGEN,IX,IY
      DOUBLE PRECISION PROJ(4),PMIN,PMAX,DP
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6002)
      ENDIF
*----
*  Find the radius of the sphere (3-D) or circle surrounding
*  the cell.
*  Also find the true center of the cell
*----
      DO IANG=1,NBANGL
        IGEN=0
        DP=DONE/DNSANG(IANG)
        DO IX=0,NUCELL(1),NUCELL(1)
          DO IY=0,NUCELL(2),NUCELL(2)
            IGEN=IGEN+1
            PROJ(IGEN)=DGMESH(IX,1)*DANGLT(2,IANG)
     >                -DGMESH(IY,2)*DANGLT(1,IANG)
          ENDDO
        ENDDO
        PMIN=PROJ(1)
        PMAX=PROJ(1)
        DO IGEN=2,4
          PMIN=MIN(PMIN,PROJ(IGEN))
          PMAX=MAX(PMAX,PROJ(IGEN))
        ENDDO
        NBSANG(5,IANG)=NINT((PMAX-PMIN)*DNSANG(IANG))+1
        PMIN=PMIN+DHALF*DP
        DEPART(1,1,IANG)=PMIN*DANGLT(2,IANG)
        DEPART(2,1,IANG)=-PMIN*DANGLT(1,IANG)
        DEPART(1,2,IANG)=DP*DANGLT(2,IANG)
        DEPART(2,2,IANG)=-DP*DANGLT(1,IANG)
        DNSANG(IANG)=DP/(DTWO*DDENWT(IANG))
        IF(IPRINT .GE. 10) THEN
          WRITE(IOUT,6003) IANG,DNSANG(IANG),NBSANG(5,IANG)
          WRITE(IOUT,6004) (DANGLT(IGEN,IANG),IGEN=1,NDIM)
          WRITE(IOUT,6005) (DEPART(IGEN,1,IANG),IGEN=1,NDIM)
          WRITE(IOUT,6006) (DEPART(IGEN,2,IANG),IGEN=1,NDIM)
        ENDIF
      ENDDO
      IF(IPRINT .GE. 10) THEN
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
 6002 FORMAT(' Spatial tracking information :')
 6003 FORMAT(' Tracking density and number of points=',
     >        I10,1X,F20.15,1X,I10)
 6004 FORMAT(' Track direction         =',3(F20.15,2X))
 6005 FORMAT(' Track starting point    =',3(F20.15,2X))
 6006 FORMAT(' Track displacement      =',3(F20.15,2X))
      END
