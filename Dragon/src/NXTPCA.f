*DECK NXTPCA
      SUBROUTINE NXTPCA(IPRINT,NDIM  ,IDIRC ,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DMESH ,NPIN  ,ITPIN ,DPIN  ,
     >                  NBSUR ,NBREG ,INDXSR,SURVOL)
*
*----------
*
*Purpose:
* Remove from the volumes or surfaces
* associated with a Cartesian 2-D or 3-D geometry
* the volumes or surfaces of the overlapping pins.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* NDIM    dimension of problem.
* IDIRC   the direction of the first axis of a Cartesian geometry
*         assuming the axis are in a cyclic rotation.
* MXMESH  maximum number of spatial subdivision in
*         $X$, $Y$ and $Z$.
* MAXSUR  maximum number of surfaces in the geometry.
* MAXREG  maximum number of regions in the geometry.
* MESH    effective number of spatial subdivision in
*         each direction ($X$, $Y$ and $Z$).
* DMESH   spatial description of the Cartesian geometry.
* NPIN    number of pins to superimpose on geometry.
* ITPIN   type of pin.
* DPIN    pin location and dimensions.
* NBSUR   number of surfaces in the geometry.
* NBREG   final number of non void regions in the geometry.
*
*Parameters: input/output
* INDXSR  local indexing of regions.
* SURVOL  volume of regions.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*Comments:
*  1- Meaning of IDIRC:
*         IDIRC      axes in 1-D   axes in 2-D   axes in 3-D
*         1          x             (x,y)         (x,y,z)
*         2          y             (y,z)         (y,z,x)
*         3          z             (z,x)         (z,x,y)
*  2- Contents of the DMESH array:
*     mesh in $X$  is x(i)=DMESH(i,1) for i=0,MESH(1);
*     mesh in $Y$  is y(j)=DMESH(j,2) for j=0,MESH(2);
*     mesh in $Z$  is z(k)=DMESH(k,3) for k=0,MESH(3);
*  3- Contents of the DPIN array for pin IPIN:
*     if(ITPIN(3,IPIN) = 1) then
*     -> annular pin
*       if(IDIRC = 1) then
*       ->annular regions in the $X-Y$ plane
*         centre (x,y,z)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         outer pin radius      r=DPIN(4,IPIN)
*         pin height            dz(iz)=DPIN(3,IPIN)
*       else if(IDIRC = 2) then
*       ->annular regions in the $Y-Z$ plane
*         centre (y,z,x)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         outer pin radius      r=DPIN(4,IPIN)
*         pin height            dx(ix)=DPIN(1,IPIN)
*       else if(IDIRC = 3) then
*       ->annular regions in the $Z-X$ plane
*         centre (z,x,y)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         outer pin radius      r=DPIN(4,IPIN)
*         pin height            dy(iy)=DPIN(2,IPIN)
*       endif
*     else if(ITPIN(3,IPIN) = 2) then
*     -> Cartesian pin
*       if(IDIRC = 1) then
*       ->Cartesian region in the $X-Y$ plane
*         centre (x,y,z)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         pin width in $X$   dx=DPIN(1,IPIN)
*         pin width in $Y$   dy=DPIN(2,IPIN)
*         pin height         dz(iz)=DPIN(3,IPIN)
*       else if(IDIRC = 2) then
*       ->Cartesian region in the $Y-Z$ plane
*         centre (y,z,x)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         pin width in $Y$   dy=DPIN(2,IPIN)
*         pin width in $Z$   dz=DPIN(3,IPIN)
*         pin height         dx(ix)=DPIN(1,IPIN)
*       else if(IDIRC = 3) then
*       ->Cartesian region in the $Z-X$ plane
*         centre (z,x,y)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         pin width in $Z$   dz=DPIN(3,IPIN)
*         pin width in $X$   dx=DPIN(1,IPIN)
*         pin height         dy(iy)=DPIN(2,IPIN)
*       endif
*     endif
*  4- Contents of the INDXSR array:
*     For i>0
*       INDXSR(1,i)= ix is the $X$ location of region i
*       INDXSR(2,i)= iy is the $Y$ location of region i
*       INDXSR(3,i)= iz is the $Z$ location of region i
*       INDXSR(4,i)= ir =0 is the $R$ location of region i.
*       INDXSR(5,i)= not used.
*     For i<0
*       INDXSR(1,i)= ix is the $X$ location of surface i
*       INDXSR(2,i)= iy is the $Y$ location of surface i
*       INDXSR(3,i)= iz is the $Z$ location of surface i
*       INDXSR(4,i)= ir =0 is the $R$ location of surface i.
*       INDXSR(5,i)= not used.
*       with INDXSR(n,i)=-1 for surface associated with
*                           location 0 in direction n.
*       with INDXSR(n,i)=-2 for surface associated with
*                           location MESH(n) in direction n.
*       Note that for radial regions INDXSR(n,i)=-1 does not
*       exists.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,IDIRC,MXMESH,MAXSUR,MAXREG
      INTEGER          MESH(4)
      DOUBLE PRECISION DMESH(-1:MXMESH,4)
      INTEGER          NPIN,ITPIN(3,NPIN)
      DOUBLE PRECISION DPIN(-1:4,NPIN)
      INTEGER          NBSUR,NBREG,INDXSR(5,-MAXSUR:MAXREG)
      DOUBLE PRECISION SURVOL(-MAXSUR:MAXREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTPCA')
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
      INTEGER          NXTIRA,NXTIRR
      INTEGER          ITYIRP
      DOUBLE PRECISION VOLINT
*----
*  Local variables
*----
      INTEGER          ID,IDG,IDIR(MAXDIM),NM(MAXDIM),IDM(MAXDIM),
     >                 IDGP1,IDGP2,IDGPP,NSOZB,NSOZT
      INTEGER          IPIN,IXYZ,IZ,IY,IX,ISUR,IVOL,IR,ILOC
      DOUBLE PRECISION ZB,ZT
      DOUBLE PRECISION POSPIN(0:2),XYPIN(4),XYCAR(4)
      DOUBLE PRECISION VOLPIN,PPRMIN,PPRMAX,PPPMIN,PPPMAX,DPP
*----
*  Data
*----
      CHARACTER        CDIR(MAXDIM)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Prepare loops over spatial directions as a function
*  of IDIRC and NDIM for Cartesian mesh.
*----
      IF(NDIM .LT. 2) CALL XABORT(NAMSBR//
     >': Problem must be 2-D or 3-D')
      ITYIRP=0
      VOLPIN=DZERO
      ZB=DZERO
      ZT=DZERO
      PI=XDRCST('Pi',' ')
      DO 100 ID=1,NDIM
        IDG=MOD(IDIRC+ID-2,3)+1
        IDIR(ID)=IDG
        NM(IDG)=MESH(IDG)
        IDM(IDG)=1
 100  CONTINUE
      DO 101 ID=NDIM+1,3
        IDG=MOD(IDIRC+ID-2,3)+1
        IDIR(ID)=IDG
        NM(IDG)=1
        IDM(IDG)=0
 101  CONTINUE
      IDGP1=IDIR(1)
      IDGP2=IDIR(2)
      IDGPP=IDIR(3)
      NSOZB=2*NM(IDGPP)*(NM(IDGP1)+NM(IDGP2))
      NSOZT=NSOZB+NM(IDGP1)*NM(IDGP2)
*----
*  Print mesh if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
*----
*  Cartesian mesh
*----
        DO 600 ID=1,NDIM
          IDG=IDIR(ID)
          WRITE(IOUT,6002) 'MESH'//CDIR(IDG)
          WRITE(IOUT,6006) (DMESH(IXYZ,IDG),IXYZ=0,NM(IDG))
          WRITE(IOUT,6003)
 600    CONTINUE
*----
* Pin description
*----
        WRITE(IOUT,6007) 'PINDIR',IDGPP
        DO 610 IPIN=1,NPIN
          WRITE(IOUT,6017) 'PinTyp',IPIN,ITPIN(3,IPIN)
          IF(ITPIN(3,IPIN) .GE. 1) THEN
            WRITE(IOUT,6019) 'PinC'//CDIR(IDGP1)//CDIR(IDGP2),IPIN,
     >                     DPIN(0,IPIN)*COS(DPIN(-1,IPIN)),
     >                     DPIN(0,IPIN)*SIN(DPIN(-1,IPIN))
            WRITE(IOUT,6018) 'PinRad',IPIN,DPIN(4,IPIN)
          ELSE IF(ITPIN(3,IPIN) .LE. -1) THEN
            WRITE(IOUT,6019) 'PinC'//CDIR(IDGP1)//CDIR(IDGP2),IPIN,
     >                     DPIN(0,IPIN)*COS(DPIN(-1,IPIN)),
     >                     DPIN(0,IPIN)*SIN(DPIN(-1,IPIN))
            WRITE(IOUT,6019) 'PinW'//CDIR(IDGP1)//CDIR(IDGP2),IPIN,
     >                    DPIN(IDGP1,IPIN),DPIN(IDGP2,IPIN)
          ENDIF
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6019) 'PinPo'//CDIR(IDGPP),IPIN,
     >                    DMESH(-1,IDGPP)-DPIN(IDGPP,IPIN)/DTWO,
     >                    DMESH(-1,IDGPP)+DPIN(IDGPP,IPIN)/DTWO
          ENDIF
 610    CONTINUE
      ENDIF
*----
*  Loop over pins
*----
      DO 110 IPIN=1,NPIN
*----
*  for 3-D problem,
*  Find pin bottom (ZB) and top (ZT) z location.
*----
        IF(NDIM .EQ. 3) THEN
          ZB=DMESH(-1,IDGPP)-DPIN(IDGPP,IPIN)/DTWO
          ZT=ZB+DPIN(IDGPP,IPIN)
        ENDIF
        IF(ITPIN(3,IPIN) .GE. 1) THEN
*----
*  Annular pin properties
*----
          POSPIN(0)=DPIN(4,IPIN)
          POSPIN(1)=DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
          POSPIN(2)=DPIN(0,IPIN)*SIN(DPIN(-1,IPIN))
          VOLPIN=PI*POSPIN(0)*POSPIN(0)
        ELSE IF (ITPIN(3,IPIN) .LE. -1) THEN
*----
*  Rectangular pin properties
*----
          XYPIN(1)=DMESH(-1,IDGP1)-DPIN(IDGP1,IPIN)/DTWO
          XYPIN(2)=DMESH(-1,IDGP1)+DPIN(IDGP1,IPIN)/DTWO
          XYPIN(3)=DMESH(-1,IDGP2)-DPIN(IDGP2,IPIN)/DTWO
          XYPIN(4)=DMESH(-1,IDGP2)+DPIN(IDGP2,IPIN)/DTWO
          VOLPIN=DPIN(IDGP2,IPIN)*DPIN(IDGP1,IPIN)
        ENDIF
*----
*  1- Loop over second normal direction
*----
        DO 120 IY=1,NM(IDGP2)
          XYCAR(3)=DMESH(IY-1,IDGP2)
          XYCAR(4)=DMESH(IY,IDGP2)
*----
*  2- Loop over first normal direction
*----
          DO 121 IX=1,NM(IDGP1)
            XYCAR(1)=DMESH(IX-1,IDGP1)
            XYCAR(2)=DMESH(IX,IDGP1)
            IF(ITPIN(3,IPIN) .GE. 1) THEN
*----
*  Find rectangle/annular pin intersection
*----
              ITYIRP=NXTIRA(XYCAR ,POSPIN,VOLINT)
            ELSE IF (ITPIN(3,IPIN) .LE. -1) THEN
*----
*  Find rectangle/rectangular pin intersection
*----
              ITYIRP=NXTIRR(XYCAR ,XYPIN ,VOLINT)
            ENDIF
            IF(ITYIRP .NE. 0) THEN
              VOLPIN=VOLPIN-VOLINT
*----
*  There is an intersection possible between the pin and
*  the rectangle.
*  1- Look for botton and top surface in 3-D
*----
              IF(NDIM .EQ. 3) THEN
                ILOC=(IY-1)*NM(IDGP1)+IX
                IF(ZB .LE. DMESH(0,IDGPP) .AND.
     >             ZT .GE. DMESH(0,IDGPP)) THEN
*----
*  Remove area contribution from bottom surface
*----
                  ISUR=NSOZB+ILOC
                  SURVOL(-ISUR)=SURVOL(-ISUR)-VOLINT
                ENDIF
                IF(ZB .LE. DMESH(NM(IDGPP),IDGPP) .AND.
     >             ZT .GE. DMESH(NM(IDGPP),IDGPP)) THEN
*----
*  Remove area contribution from top surface
*----
                  ISUR=NSOZT+ILOC
                  SURVOL(-ISUR)=SURVOL(-ISUR)-VOLINT
                ENDIF
                PPPMIN=ZB
                PPPMAX=PPPMIN+DPIN(IDGPP,IPIN)
                DO 130 IZ=1,NM(IDGPP)
                  PPRMIN=MAX(DMESH(IZ-1,IDGPP),PPPMIN)
                  PPRMAX=MIN(DMESH(IZ,IDGPP),PPPMAX)
                  IF(PPRMIN .LT. PPRMAX) THEN
                    DPP=VOLINT*(PPRMAX-PPRMIN)
                    IVOL=((IZ-1)*NM(IDGP2)+(IY-1))*NM(IDGP1)+IX
                    SURVOL(IVOL)=SURVOL(IVOL)-DPP
                  ENDIF
 130            CONTINUE
              ELSE
                IVOL=(IY-1)*NM(IDGP1)+IX
                SURVOL(IVOL)=SURVOL(IVOL)-VOLINT
              ENDIF
            ENDIF
*----
*  If pin all extracted, go to next pin
*----
            IF(VOLPIN .LE. DZERO) GO TO 115
 121      CONTINUE
 120    CONTINUE
 115    CONTINUE
 110  CONTINUE
*----
*  Test for negative surface area and volumes
*----
      DO 140 ISUR=1,NBSUR
        IF(SURVOL(-ISUR) .LT. -DCUTOF) THEN
          WRITE(IOUT,9000) NAMSBR,-ISUR
          WRITE(IOUT,9002) (INDXSR(IR,-ISUR),IR=1,5),SURVOL(-ISUR)
          CALL XABORT(NAMSBR//
     >    ': Region with negative surface area detected')
        ELSE IF(SURVOL(-ISUR) .LT. DCUTOF) THEN
          SURVOL(-ISUR)=DZERO
        ENDIF
 140  CONTINUE
      DO 141 IVOL=1,NBREG
        IF(SURVOL(IVOL) .LT. -DCUTOF) THEN
          WRITE(IOUT,9001) NAMSBR,IVOL
          WRITE(IOUT,9002) (INDXSR(IR,IVOL),IR=1,5),SURVOL(IVOL)
          CALL XABORT(NAMSBR//
     >    ': Region with negative volume detected')
        ELSE IF(SURVOL(IVOL) .LT. DCUTOF) THEN
          SURVOL(IVOL)=DZERO
        ENDIF
 141  CONTINUE
*----
*  Print volumes if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVol'
        WRITE(IOUT,6005) (IVOL,(INDXSR(IR,IVOL),IR=1,5),SURVOL(IVOL),
     >                     IVOL=-NBSUR,NBREG)
        WRITE(IOUT,6003)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(A12,'={')
 6003 FORMAT('};')
 6005 FORMAT((6(I10,','),F20.10,:,','))
 6006 FORMAT(6(F20.10,:,','))
 6007 FORMAT(A6,'=',I10,';')
 6017 FORMAT(A6,I4.4,'=',I10,';')
 6018 FORMAT(A6,I4.4,'=',F20.10,';')
 6019 FORMAT(A6,I4.4,'={',F20.10,',',F20.10,'};')
*----
*  Error and Warning formats
*----
 9000 FORMAT('**** ERROR in -- ',A6,'-- found'/
     >       '     Area of region ',I5,' is negative')
 9001 FORMAT('**** ERROR in -- ',A6,'-- found'/
     >       '     Volume of region ',I5,' is negative')
 9002 FORMAT(5I10,F20.10)
      END
