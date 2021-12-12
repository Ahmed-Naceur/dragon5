*DECK TLMGEO
      SUBROUTINE TLMGEO(IPTRK,IPMAT,IPRINT,ITGEO,MAXMDH,NDIM,NUCELL,
     >                  DGMESH,XYZL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the Matlab instruction for drawing the global geometry.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* C. Plamondon, G. Marleau
*
*Parameters: input
* IPTRK   tracking data structure.
* IPMAT   pointer to Matlab-m file.
* IPRINT  print level.
* ITGEO   type of geometry (0 for annular; 1 for Cartesian;
*         2 for hexagonal).
* MAXMDH  maximum dimensions of DGMESH.
* NDIM    number of dimensions for problem.
* NUCELL  cell dimensions for each direction.
* DGMESH  mesh of global beometry.
*
*Parameters: output
* XYZL    mesh limits.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPMAT
      INTEGER          IPRINT,ITGEO,MAXMDH,NDIM
      INTEGER          NUCELL(NDIM)
      DOUBLE PRECISION DGMESH(-1:MAXMDH,4),XYZL(2,3)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMGEO')
      DOUBLE PRECISION DZERO
      PARAMETER       (DZERO=0.0D0)
*----
*  Other local variables
*----
      INTEGER          IDIR,NX,NY,NZ,IPT,IH
      CHARACTER        NAMREC*12
      DOUBLE PRECISION POSHRX(7),POSHRY(7),POSHDX(7),POSHDY(7),
     >                 ROTMAT(2,2)
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      DO IDIR=1,NDIM
        NAMREC='G00000001SM'//CDIR(IDIR)
        IF(ITGEO .EQ. 2) THEN
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
        ELSE
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          XYZL(1,IDIR)=DGMESH(0,IDIR)
          XYZL(2,IDIR)=DGMESH(NUCELL(IDIR),IDIR)
        ENDIF
      ENDDO
      WRITE(IPMAT,7000)
      NX=NUCELL(1)
      NY=NUCELL(2)
*      write(6,*) 'ITGEO,NDIM=',ITGEO,NDIM,NX,NY
*      write(6,*) 'DGMESHx=',(DGMESH(IH,1),IH=-1,NX)
*      write(6,*) 'DGMESHy=',(DGMESH(IH,2),IH=-1,NX)
*----
*  For hexagon, find reference corner positions
*----
      IF(ABS(ITGEO) .EQ. 2) THEN
*        IF(ITGEO .EQ. 2) THEN
*----
*  One side parallel to x-axis
*----
*          POSHRX(1)=DGMESH(0,1)
*          POSHRY(1)=DZERO
*        ELSE
*----
*  One side parallel to y-axis
*----
          POSHRX(1)=DZERO
          POSHRY(1)=DGMESH(0,2)
*        ENDIF
        ROTMAT(1,1)=0.5D0
        ROTMAT(2,1)=SQRT(3.0D0)/2.0D0
        ROTMAT(1,2)=-ROTMAT(2,1)
        ROTMAT(2,2)=ROTMAT(1,1)
        DO IPT=2,7
         POSHRX(IPT)=ROTMAT(1,1)*POSHRX(IPT-1)
     >              +ROTMAT(1,2)*POSHRY(IPT-1)
         POSHRY(IPT)=ROTMAT(2,1)*POSHRX(IPT-1)
     >              +ROTMAT(2,2)*POSHRY(IPT-1)
        ENDDO
        WRITE(IPMAT,7040) CDIR(1),NX
        WRITE(IPMAT,7040) CDIR(2),NX
        IF(NDIM .EQ. 3) THEN
          WRITE(IPMAT,7040) CDIR(3),NZ
        ENDIF
        DO IH=1,NX
          DO IPT=1,7
            POSHDX(IPT)=POSHRX(IPT)+DGMESH(IH,1)
            POSHDY(IPT)=POSHRY(IPT)+DGMESH(IH,2)
          ENDDO
          WRITE(IPMAT,7041) CDIR(1),IH,(POSHDX(IPT),IPT=1,7)
          WRITE(IPMAT,7041) CDIR(2),IH,(POSHDY(IPT),IPT=1,7)
        ENDDO
      ENDIF
*----
* Print IPMAT header
*----
      IF(NDIM .EQ. 2) THEN
        IF(ITGEO .EQ. 2) THEN
          DO IH=1,NX
            WRITE(IPMAT,7042) IH,IH
          ENDDO
        ELSE
          WRITE(IPMAT,7020)
     >    DGMESH(0,1),DGMESH(NX,1),
     >    DGMESH(0,2),DGMESH(NY,2)
          WRITE(IPMAT,7021)
        ENDIF
      ELSE IF(NDIM .EQ. 3) THEN
        NZ=NUCELL(3)
        WRITE(IPMAT,7001)
        IF(ITGEO .EQ. 2) THEN
          WRITE(IPMAT,7043)
     >    DGMESH(1,3),DGMESH(NZ,3)
          DO IH=1,NX
            WRITE(IPMAT,7042) IH,IH
          ENDDO
        ELSE
          WRITE(IPMAT,7030)
     >    DGMESH(1,1),DGMESH(NX,1),
     >    DGMESH(1,2),DGMESH(NY,2),
     >    DGMESH(1,3),DGMESH(NZ,3)
          WRITE(IPMAT,7031)
        ENDIF
      ENDIF
*----
*  Processing finished, return
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
*----
*  Matlab .m file format
*----
 7000 FORMAT('figure;'/'hold on;'/'axis equal;'/
     >12Hxlabel('x');/12Hylabel('y');)
 7001 FORMAT(12Hzlabel('z');)
 7020 FORMAT('xmin=',F18.10,';'/'xmax=',F18.10,';'/
     >       'ymin=',F18.10,';'/'ymax=',F18.10,';'/)
 7021 FORMAT('line([xmin,xmin],[ymin,ymax]);'/
     >       'line([xmin,xmax],[ymin,ymin]);'/
     >       'line([xmax,xmax],[ymin,ymax]);'/
     >       'line([xmin,xmax],[ymax,ymax]);')
 7030 FORMAT('xmin=',F18.10,';'/'xmax=',F18.10,';'/
     >       'ymin=',F18.10,';'/'ymax=',F18.10,';'/
     >       'zmin=',F18.10,';'/'zmax=',F18.10,';'/)
 7031 FORMAT('line([xmin,xmin],[ymin,ymax],[zmin,zmin]);'/
     >       'line([xmin,xmax],[ymin,ymin],[zmin,zmin]);'/
     >       'line([xmax,xmax],[ymin,ymax],[zmin,zmin]);'/
     >       'line([xmin,xmax],[ymax,ymax],[zmin,zmin]);'/
     >       'line([xmin,xmin],[ymin,ymin],[zmin,zmax]);'/
     >       'line([xmax,xmax],[ymin,ymin],[zmin,zmax]);'/
     >       'line([xmax,xmax],[ymax,ymax],[zmin,zmax]);'/
     >       'line([xmin,xmin],[ymax,ymax],[zmin,zmax]);'/
     >       'line([xmin,xmin],[ymin,ymax],[zmax,zmax]);'/
     >       'line([xmin,xmax],[ymin,ymin],[zmax,zmax]);'/
     >       'line([xmax,xmax],[ymin,ymax],[zmax,zmax]);'/
     >       'line([xmin,xmax],[ymax,ymax],[zmax,zmax]);')
 7040 FORMAT('HexM',A1,'=zeros(',I4,',7);')
 7041 FORMAT('HexM',A1,'(',I4,',:)=[',7F18.10,'];')
 7042 FORMAT('line(HexMX(',I4,',:),HexMY(',I4,',:));')
 7043 FORMAT('zmin=',F18.10,';'/'zmax=',F18.10,';'/)
      END
