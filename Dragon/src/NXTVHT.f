*DECK NXTVHT
      SUBROUTINE NXTVHT(IPRINT,NDIM  ,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DMESH ,NBSUR ,NBREG ,INDXSR,SURVOL)
*
*----------
*
*Purpose:
* Compute the volume and area associated with each region
* or surface for an hexagon with triangular mesh 
* in 2-D or 3-D geometry.
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
*Comments:
*  1- Contents of the DMESH array:
*     hexagonal mesh is DMESH(i,1) for i=0,MESH(1);
*     mesh in $Z$  is z(k)=DMESH(k,3) for k=0,MESH(3);
*     annular regions in the $X-Y$ plane
*     centre of cylinder in (x,y)=(DMESH(-1,1),DMESH(-1,2))
*     radius of shells      r(l)=DMESH(l,4), l=1,MESH(4)
*  2- Contents of the INDXSR array:
*     For i>0
*       INDXSR(1,i)= iu is the $U$ location of region i
*       INDXSR(2,i)= iv is the $V$ location of region i
*       INDXSR(3,i)= iz is the $Z$ location of region i
*       INDXSR(4,i)= ir is the $R$ location of region i
*       INDXSR(5,i)= iw is the $W$ location of region i
*     For i<0
*       INDXSR(1,i)= iu is the $U$ location of surface i
*       INDXSR(2,i)= iv is the $V$ location of surface i
*       INDXSR(3,i)= iz is the $Z$ location of surface i
*       INDXSR(4,i)= ir is the $R$ location of surface i
*       INDXSR(5,i)= iw is the $W$ location of surface i
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
      INTEGER          IPRINT,NDIM,MXMESH,MAXSUR,MAXREG
      INTEGER          MESH(4)
      DOUBLE PRECISION DMESH(-1:MXMESH,4)
      INTEGER          NBSUR,NBREG,INDXSR(5,-MAXSUR:MAXREG)
      DOUBLE PRECISION SURVOL(-MAXSUR:MAXREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTVHT')
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
      DOUBLE PRECISION DZERO,DONE,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,
     >                 DHALF=0.5D0,DSQ3O2=0.86602540378444D0)
*----
*  Local variables
*----
      INTEGER          NX,NZ,NRTP,NRP,NSTP,NSP,IX,IZ,ISECT,IR,IS,
     >                 IOFS,IVSI,JVSI,IVSIR,NVSIR,ISUVW,IZ1
      DOUBLE PRECISION SIDEU,SIDEUP,SIDER,SIDERP,SIDEL,SIDELP,
     >                 AREAT,AREAL,AREAR,DZ
*----
*  Prepare loop over Z-direction.
*----
      NX=MESH(1)
      NZ=MESH(3)
      NRTP=NX**2
      NRP=6*NRTP
      NSTP=2*NX-1
      NSP=6*NSTP
      IF(NDIM .EQ. 3) THEN
        NBREG=NRP*NZ
        NBSUR=NSP*NZ+2*NRP
      ELSE
        NBREG=NRP
        NBSUR=NSP
      ENDIF
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) NX,NZ
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DMESH(IX,1),IX=-1,2*NX)
        IF(NDIM .EQ. 3) THEN
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DMESH(IZ,3),IZ=-1,NZ)
        ENDIF
        WRITE(IOUT,6002) 'SurVol init'
        WRITE(IOUT,6005) (IVSI,(INDXSR(IR,IVSI),IR=1,5),SURVOL(IVSI),
     >                     IVSI=-NBSUR,NBREG)
        WRITE(IOUT,6003)
      ENDIF
      CALL XDDSET(SURVOL,MAXSUR+1+MAXREG,DZERO)
*----
*  All the triangles are equilateral and have the same surface area
*  except possibly those in the last crown because this last crown
*  may be thinner than the other
*  1) Find first the area of generic triangle
*----
      SIDEU=DMESH(NX+1,1)-DMESH(NX,1)
      SIDEUP=SIDEU/DSQ3O2
      AREAT=DHALF*SIDEU*SIDEUP
      CALL XDDSET(SURVOL(1),NRP,AREAT)
*----
*  Process last crown
*  for right triangle remove volume past last cell boundary
*  for left triangle volume extends only to cell boundary
*  Add surfaces
*----
      SIDEL=DMESH(1,1)-DMESH(0,1)
      SIDELP=SIDEL/DSQ3O2
      SIDER=SIDEU-SIDEL
      SIDERP=SIDER/DSQ3O2
      AREAL=DHALF*SIDEL*SIDELP
      AREAR=AREAT-DHALF*SIDER*SIDERP
*----
*  Finish filling first plan with AREA
*----
      IS=0
      IOFS=NRTP-NSTP
      DO ISECT=1,6
*----
*  right triangles
*----
        DO IX=1,NX-1
          IR=IX+IOFS
          IS=IS-1
          SURVOL(IR)=AREAR
          SURVOL(IS)=SIDERP
        ENDDO
*----
*  Left triangles
*----
        DO IX=NX,NSTP
          IR=IX+IOFS
          IS=IS-1
          SURVOL(IR)=AREAL
          SURVOL(IS)=SIDELP
        ENDDO
        IOFS=IOFS+NRTP
      ENDDO
      IF(NDIM .EQ. 3) THEN
*----
* Hexagonal Surfaces
*----
        DO IZ=NZ,1,-1
          DZ=DMESH(IZ,3)-DMESH(IZ-1,3)
          IX=-NSP*(IZ-1)
          DO IS=-1,-NSP,-1
            SURVOL(IX+IS)=DZ*SURVOL(IS)
          ENDDO
        ENDDO
*----
*  Fill bottom face and top face
*----
        IS=-NSP*NZ
        DO IR=1,NRP
          IS=IS-1
          SURVOL(IS)=SURVOL(IR)
          SURVOL(IS-NRP)=SURVOL(IR)
        ENDDO
*----
*  Fill all planes with volumes
*----
        DO IZ=NZ,1,-1
          DZ=DMESH(IZ,3)-DMESH(IZ-1,3)
          IX=NRP*(IZ-1)
*----
*  Volumes
*----
          DO IR=1,NRP
            SURVOL(IX+IR)=DZ*SURVOL(IR)
          ENDDO
        ENDDO
      ENDIF
*----
*  Fill INDXSR to identify regions and volumes
*  Process first sector (-pi/6 to pi/6) for volume
*  First plane or 2-D
*----
      IF(NDIM .EQ. 3) THEN
        IZ1=1
      ELSE
        IZ1=0
      ENDIF
      IVSI=0
      DO IX=1,NX
*----
*  First line of triangle in the crown
*----
        DO IR=1,IX-1
          IVSI=IVSI+1
          INDXSR(1,IVSI)=NX+IX
          INDXSR(5,IVSI)=NX+IX-IR
          INDXSR(3,IVSI)=IZ1
          INDXSR(2,IVSI)=NX+IR
        ENDDO
*----
*  Second line of triangle in the crown
*----
        DO IR=1,IX
          IVSI=IVSI+1
          INDXSR(1,IVSI)=NX+IX
          INDXSR(5,IVSI)=NX+IX-IR+1
          INDXSR(3,IVSI)=IZ1
          INDXSR(2,IVSI)=NX+IR
        ENDDO
      ENDDO
*----
*  Complete for sectors 2-6 on first plane
*----
      NVSIR=IVSI
      IVSIR=NVSIR
      ISUVW=2*NX+1
      DO ISECT=2,6
        IF(ISECT .EQ. 6) THEN
          DO JVSI=1,NVSIR
            IVSIR=IVSIR+1
            INDXSR(1,IVSIR)=INDXSR(2,JVSI)
            INDXSR(2,IVSIR)=ISUVW-INDXSR(5,JVSI)
            INDXSR(3,IVSIR)=IZ1
            INDXSR(5,IVSIR)=INDXSR(1,JVSI)
          ENDDO
        ELSE IF (ISECT.EQ.5) THEN
          DO JVSI=1,NVSIR
            IVSIR=IVSIR+1
            INDXSR(1,IVSIR)=ISUVW-INDXSR(5,JVSI)
            INDXSR(2,IVSIR)=ISUVW-INDXSR(1,JVSI)
            INDXSR(3,IVSIR)=IZ1
            INDXSR(5,IVSIR)=INDXSR(2,JVSI)
          ENDDO
        ELSE IF (ISECT.EQ.4) THEN
          DO JVSI=1,NVSIR
            IVSIR=IVSIR+1
            INDXSR(1,IVSIR)=ISUVW-INDXSR(1,JVSI)
            INDXSR(2,IVSIR)=ISUVW-INDXSR(2,JVSI)
            INDXSR(3,IVSIR)=IZ1
            INDXSR(5,IVSIR)=ISUVW-INDXSR(5,JVSI)
          ENDDO
        ELSE IF (ISECT.EQ.3) THEN
          DO JVSI=1,NVSIR
            IVSIR=IVSIR+1
            INDXSR(1,IVSIR)=ISUVW-INDXSR(2,JVSI)
            INDXSR(2,IVSIR)=INDXSR(5,JVSI)
            INDXSR(3,IVSIR)=IZ1
            INDXSR(5,IVSIR)=ISUVW-INDXSR(1,JVSI)
          ENDDO
        ELSE IF (ISECT.EQ.2) THEN
          DO JVSI=1,NVSIR
            IVSIR=IVSIR+1
            INDXSR(1,IVSIR)=INDXSR(5,JVSI)
            INDXSR(2,IVSIR)=INDXSR(1,JVSI)
            INDXSR(3,IVSIR)=IZ1
            INDXSR(5,IVSIR)=ISUVW-INDXSR(2,JVSI)
          ENDDO
        ENDIF
      ENDDO
*----
*  Process other planes in 3-D
*----
      IVSI=NRP
      DO IZ=2,NZ
        DO IR=1,NRP
          IVSI=IVSI+1
          INDXSR(1,IVSI)=INDXSR(1,IR)
          INDXSR(2,IVSI)=INDXSR(2,IR)
          INDXSR(3,IVSI)=IZ
          INDXSR(5,IVSI)=INDXSR(5,IR)
        ENDDO
      ENDDO
*----
*  Process surfaces
*  First plane or 2-D
*----
      IVSI=0
      IOFS=3-6*NX
      IX=NX
*----
*  Sector 1 and 4
*  First line of triangle on the last crown
*----
      DO IR=1,IX-1
        IVSI=IVSI-1
        INDXSR(1,IVSI+IOFS)=-1
        INDXSR(5,IVSI+IOFS)=IR+1
        INDXSR(3,IVSI+IOFS)=IZ1
        INDXSR(2,IVSI+IOFS)=NX+1-IR
        INDXSR(1,IVSI)=-2
        INDXSR(5,IVSI)=NX+IX-IR
        INDXSR(3,IVSI)=IZ1
        INDXSR(2,IVSI)=NX+IR
      ENDDO
*----
*  Second line of triangle on the last crown
*----
      DO IR=1,IX
        IVSI=IVSI-1
        INDXSR(1,IVSI+IOFS)=-1
        INDXSR(5,IVSI+IOFS)=IR
        INDXSR(3,IVSI+IOFS)=IZ1
        INDXSR(2,IVSI+IOFS)=NX+1-IR
        INDXSR(1,IVSI)=-2
        INDXSR(5,IVSI)=NX+IX-IR+1
        INDXSR(3,IVSI)=IZ1
        INDXSR(2,IVSI)=NX+IR
      ENDDO
*----
*  Sector 2 and 5
*  First line of triangle on the last crown
*----
      DO IR=1,IX-1
        IVSI=IVSI-1
        INDXSR(1,IVSI+IOFS)=IR+1
        INDXSR(5,IVSI+IOFS)=NX+IR
        INDXSR(3,IVSI+IOFS)=IZ1
        INDXSR(2,IVSI+IOFS)=-1
        INDXSR(1,IVSI)=2*NX-IR
        INDXSR(5,IVSI)=NX+1-IR
        INDXSR(3,IVSI)=IZ1
        INDXSR(2,IVSI)=-2
      ENDDO
*----
*  Second line of triangle on the last crown
*----
      DO IR=1,IX
        IVSI=IVSI-1
        INDXSR(1,IVSI+IOFS)=IR
        INDXSR(5,IVSI+IOFS)=NX+IR
        INDXSR(3,IVSI+IOFS)=IZ1
        INDXSR(2,IVSI+IOFS)=-1
        INDXSR(1,IVSI)=2*NX+1-IR
        INDXSR(5,IVSI)=NX+1-IR
        INDXSR(3,IVSI)=IZ1
        INDXSR(2,IVSI)=-2
      ENDDO
*----
*  Sector 3 and 6
*  First line of triangle on the last crown
*----
      DO IR=1,IX-1
        IVSI=IVSI-1
        INDXSR(1,IVSI)=NX+1-IR
        INDXSR(5,IVSI)=-1
        INDXSR(3,IVSI)=IZ1
        INDXSR(2,IVSI)=2*NX-IR
        INDXSR(1,IVSI+IOFS)=NX+IR
        INDXSR(5,IVSI+IOFS)=-2
        INDXSR(3,IVSI+IOFS)=IZ1
        INDXSR(2,IVSI+IOFS)=IR+1
      ENDDO
*----
*  Second line of triangle on the last crown
*----
      DO IR=1,IX
        IVSI=IVSI-1
        INDXSR(1,IVSI)=NX+1-IR
        INDXSR(5,IVSI)=-1
        INDXSR(3,IVSI)=IZ1
        INDXSR(2,IVSI)=2*NX+1-IR
        INDXSR(1,IVSI+IOFS)=NX+IR
        INDXSR(5,IVSI+IOFS)=-2
        INDXSR(3,IVSI+IOFS)=IZ1
        INDXSR(2,IVSI+IOFS)=IR
      ENDDO
      IF(NDIM .EQ. 3) THEN
*----
*  Process other planes in 3-D
*----
        IVSI=IVSI-6*NX+3
        DO IZ=2,NZ
          DO IR=-1,-NSP,-1
            IVSI=IVSI-1
            INDXSR(1,IVSI)=INDXSR(1,IR)
            INDXSR(2,IVSI)=INDXSR(2,IR)
            INDXSR(3,IVSI)=IZ
            INDXSR(5,IVSI)=INDXSR(5,IR)
          ENDDO
        ENDDO
*----
*  Process bottom and top faces if required
*----
        DO IR=1,NRP
          IVSI=IVSI-1
          IVSIR=IVSI-NRP
          INDXSR(1,IVSI)=INDXSR(1,IR)
          INDXSR(2,IVSI)=INDXSR(2,IR)
          INDXSR(3,IVSI)=-1
          INDXSR(5,IVSI)=INDXSR(5,IR)
          INDXSR(1,IVSIR)=INDXSR(1,IR)
          INDXSR(2,IVSIR)=INDXSR(2,IR)
          INDXSR(3,IVSIR)=-2
          INDXSR(5,IVSIR)=INDXSR(5,IR)
        ENDDO
      ENDIF
*----
*  Print surfaces and volumes if required and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVol'
        WRITE(IOUT,6005) (IVSI,(INDXSR(IR,IVSI),IR=1,5),SURVOL(IVSI),
     >                     IVSI=-NBSUR,NBREG)
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
 6010 FORMAT(1X,'MESH DIMENSIONS =',2I10)
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F20.10)
      END
