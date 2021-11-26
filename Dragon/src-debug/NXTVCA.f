*DECK NXTVCA
      SUBROUTINE NXTVCA(IPRINT,NDIM  ,IDIRC ,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DMESH ,NBSUR ,NBREG ,INDXSR,SURVOL)
*
*----------
*
*Purpose:
* Compute the volume and area associated with each region
* or surface for a Cartesian 1-D, 2-D or 3-D geometry.
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
* DMESH   spatial description of the parallepiped.
*
*Parameters: output
* NBSUR   number of surfaces in the geometry.
* NBREG   number of regions in the geometry.
* INDXSR  local indexing of surfaces/regions.
* SURVOL  area/volume of regions.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*Comments:
*  1- Contents of IDIRC:
*         IDIRC      axes in 1-D   axes in 2-D   axes in 3-D
*         1          x             (x,y)         (x,y,z)
*         2          y             (y,z)         (y,z,x)
*         3          z             (z,x)         (z,x,y)
*  2- Contents of the DMESH array:
*     mesh in $X$  is x(i)=DMESH(i,1) for i=0,MESH(1);
*     mesh in $Y$  is y(j)=DMESH(j,2) for j=0,MESH(2);
*     mesh in $Z$  is z(k)=DMESH(k,3) for k=0,MESH(3);
*     if(IDIRC = 1) then
*     ->annular regions in the $X-Y$ plane
*       centre of cylinder in (x,y)=(DMESH(-1,1),DMESH(-1,2))
*       radius of shells      r(l)=DMESH(l,4), l=1,MESH(4)
*     else if(IDIRC = 2) then
*     ->annular regions in the $Y-Z$ plane
*       centre of cylinder in (y,z)=(DMESH(-1,2),DMESH(-1,3))
*       radius of shells      r(l)=DMESH(l,4), l=1,MESH(4)
*     else if(IDIRC = 3) then
*     ->annular regions in the $Z-X$ plane
*       centre of cylinder in (z,x)=(DMESH(-1,3),DMESH(-1,1))
*       radius of shells      r(l)=DMESH(l,4), l=1,MESH(4)
*     endif
*  3- Contents of the INDXSR array:
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
*       INDXSR(4,i)= ir is the $R$ location of surface i.
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
      INTEGER          NBSUR,NBREG,INDXSR(5,-MAXSUR:MAXREG)
      DOUBLE PRECISION SURVOL(-MAXSUR:MAXREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTVCA')
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Local variables
*----
      INTEGER          IDIR(MAXDIM),NM(MAXDIM),IDM(MAXDIM)
      INTEGER          ID,IDG,IDGP1,IDGP2,IDGPP,IX,IY,IZ,IXYZ,IR,
     >                 ISUR,ISUR2,IVOL,IDIRCX
      DOUBLE PRECISION DX,DY,DZ
*----
*  Data
*----
      CHARACTER        CDIR(MAXDIM)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Prepare loops over spatial directions as a function
*  of IDIRC and NDIM.
*----
      CALL XDDSET(SURVOL,(MAXSUR+MAXREG+1),DZERO)
      IDG=0
      IDGP1=0
      IDGP2=0
      IR=0
      NBREG=1
      DO 100 ID=1,NDIM
        IDG=MOD(IDIRC+ID-2,3)+1
        IDIR(ID)=IDG
        NM(IDG)=MESH(IDG)
        NBREG=NBREG*NM(IDG)
        IDM(IDG)=1
 100  CONTINUE
      DO 101 ID=NDIM+1,3
        IDG=MOD(IDIRC+ID-2,3)+1
        IDIR(ID)=IDG
        NM(IDG)=1
        IDM(IDG)=0
 101  CONTINUE
      IF(MAXREG .LT. NBREG) CALL XABORT(NAMSBR//
     >': Insufficient space to store region volumes')
*----
*  number of surfaces
*----
      NBSUR=0
      DO 102 ID=1,NDIM
        IDG=MOD(IDIRC+ID-2,3)+1
        NBSUR=NBSUR+2*NBREG/NM(IDG)
 102  CONTINUE
      IF(MAXSUR .LT. NBSUR) CALL XABORT(NAMSBR//
     >': Insufficient space to store surface areas')
*----
*  Print mesh if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        DO 600 ID=1,NDIM
          IDG=IDIR(ID)
          WRITE(IOUT,6002) 'MESH'//CDIR(IDG)
          WRITE(IOUT,6006) (DMESH(IXYZ,IDG),IXYZ=0,NM(IDG))
          WRITE(IOUT,6003)
 600    CONTINUE
      ENDIF
*----
*  Compute surface area
*  1- Loop over directions (110)
*----
      ISUR=0
      IDIRCX=1
      DO 110 ID=1,NDIM
        IDG=MOD(IDIRCX+ID-2,3)+1
        IF(IDM(IDG) .EQ. 1) THEN
          IDGP1=MOD(IDIRCX+ID-1,3)+1
          IDGP2=MOD(IDIRCX+ID,3)+1
*----
*  2- Loop over second normal direction (111)
*----
          ISUR2=ISUR+NM(IDGP2)*NM(IDGP1)
          DO 111 IY=1,NM(IDGP2)
            DY=DONE
            IF(IDM(IDGP2) .EQ. 1)
     >      DY=DMESH(IY,IDGP2)-DMESH(IY-1,IDGP2)
*----
*  3- Loop over first normal direction (112)
*----
            DO 112 IX=1,NM(IDGP1)
              DX=DONE
              IF(IDM(IDGP1) .EQ. 1)
     >        DX=DMESH(IX,IDGP1)-DMESH(IX-1,IDGP1)
              ISUR=ISUR+1
              SURVOL(-ISUR)=DX*DY
              INDXSR(IDGP1,-ISUR)=IX
              INDXSR(IDGP2,-ISUR)=IY
              INDXSR(IDG,-ISUR)=-1
              INDXSR(4,-ISUR)=IR
              ISUR2=ISUR2+1
              SURVOL(-ISUR2)=DX*DY
              INDXSR(IDGP1,-ISUR2)=IX
              INDXSR(IDGP2,-ISUR2)=IY
              INDXSR(IDG,-ISUR2)=-2
              INDXSR(4,-ISUR2)=IR
 112        CONTINUE
 111      CONTINUE
          ISUR=ISUR2
        ENDIF
 110  CONTINUE
*----
*  Computes regional volumes
*  1- Loop on $Z$ (120)
*----
      IR=0
      IVOL=0
      SURVOL(IVOL)=DZERO
      INDXSR(IDGP1,IVOL)=0
      INDXSR(IDGP2,IVOL)=0
      INDXSR(IDG,IVOL)=0
      INDXSR(4,IVOL)=0
      IDGPP=IDIR(3)
      IDGP2=IDIR(2)
      IDGP1=IDIR(1)
      IDGPP=3
      IDGP2=2
      IDGP1=1
      DO 120 IZ=1,NM(IDGPP)
        DZ=DONE
        IF(IDM(IDGPP) .EQ. 1)
     >  DZ=DMESH(IZ,IDGPP)-DMESH(IZ-1,IDGPP)
*----
*  2- Loop on $Y$ (121)
*----
        DO 121 IY=1,NM(IDGP2)
          DY=DONE
          IF(IDM(IDGP2) .EQ. 1)
     >    DY=DMESH(IY,IDGP2)-DMESH(IY-1,IDGP2)
*----
*  3- Loop on $X$ (122)
*----
          DO 122 IX=1,NM(IDGP1)
            DX=DONE
            IF(IDM(IDGP1) .EQ. 1)
     >      DX=DMESH(IX,IDGP1)-DMESH(IX-1,IDGP1)
            IVOL=IVOL+1
            SURVOL(IVOL)=DX*DY*DZ
            INDXSR(IDGP1,IVOL)=IX
            INDXSR(IDGP2,IVOL)=IY
            INDXSR(IDGPP,IVOL)=IZ
            INDXSR(4,IVOL)=IR
 122      CONTINUE
 121    CONTINUE
 120  CONTINUE
*----
*  Print volumes if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVol'
        WRITE(IOUT,6005) (IVOL,(INDXSR(IR,IVOL),IR=1,4),SURVOL(IVOL),
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
 6005 FORMAT((5(I10,','),F20.10,:,','))
 6006 FORMAT(4(F20.10,:,','))
      END
