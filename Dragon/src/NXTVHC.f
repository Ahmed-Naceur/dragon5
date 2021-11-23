*DECK NXTVHC
      SUBROUTINE NXTVHC(IPRINT,NDIM  ,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DMESH ,NBSUR ,NBREG ,INDXSR,SURVOL,
     >                  POSTRI)
*
*----------
*
*Purpose:
* Compute the volume and area associated with each region
* or surface for a annular/hexagon with triangular mesh 
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
* POSTRI  triangle position:
*         POSTRI(1,*,*,*) is X position;
*         POSTRI(2,*,*,*) is Y position;
*         POSTRI(*,1,*,*) is location of first corner;
*         POSTRI(*,2,*,*) is location of second corner;
*         POSTRI(*,3,*,*) is location of third corner;
*         POSTRI(*,*,i,j) is location of triangle i in cector j.
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
      DOUBLE PRECISION POSTRI(2,3,MXMESH*MXMESH,6)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTVHC')
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
      DOUBLE PRECISION DZERO
      PARAMETER       (DZERO=0.0D0)
*----
*  Functions
*----
      INTEGER          NXTITA,ITYITA
      DOUBLE PRECISION VOLINT
*----
*  Local variables
*----
      INTEGER          NX,NZ,NR,NRTP,NRP,NSTP,NSP,NRTPP,NRP1,
     >                 IT,IZ,IR,ISECT,IOFS,IVSI,NBVHEX,NBSHEX
      DOUBLE PRECISION POSANN(0:2),VOLZ
      INTEGER          ILVI,ILVT,IDIR,IOFZ,IS1,ISURTF,ISURTI,
     >                 IVOLF,IVOLI,NBREGT,
     >                 NSKPF,NSKPI,NBSURT
*----
*  Prepare loop over Z-direction.
*----
      NX=MESH(1)
      NZ=MESH(3)
      NR=MESH(4)
      NRTP=NX**2
      NRP=6*NRTP
      NSTP=2*NX-1
      NSP=6*NSTP
      NRP1=NR+1
      NRTPP=NRP*NRP1
*----
*  Compute number of surfaces
*  1- Surface parallel to cylinder axis
*----
      NBVHEX=NRP
      NBREG=NRTPP
      NBSUR=NSP
      NBSHEX=NSP
      IF(NDIM .EQ. 3) THEN
        NBREG=NBREG*NZ
        NBSUR=NBSHEX*NZ+2*NRTPP
        NBVHEX=NBVHEX*NZ
        NBSHEX=NBSHEX*NZ+2*NRP
      ENDIF
      NBREGT=NBVHEX
      NBSURT=NBSHEX
*----
*  2- Surface normal to cylinder axis (if any)
*----
      NSKPI=0
      NSKPF=0
      IF(NDIM .NE. 3) THEN
        NSKPI=NRP
        NSKPF=(NR+1)*NSKPI
      ENDIF
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) 'H',NX
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DMESH(IT,1),IT=-1,2*NX)
        IF(NDIM .EQ. 3) THEN
          WRITE(IOUT,6010) 'Z',NZ
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DMESH(IZ,3),IZ=-1,NZ)
        ENDIF
          WRITE(IOUT,6010) 'R',NR
        WRITE(IOUT,6011) 'MESHR ='
        WRITE(IOUT,6012) (DMESH(IR,4),IR=-1,NR)
      ENDIF
      CALL XDDSET(SURVOL,MAXSUR+1+MAXREG,DZERO)      
*----
*  Call NXTVHT to obtain the triangles volumes and external surfaces
*----
      CALL NXTVHT(IPRINT,NDIM  ,MXMESH,NBSHEX,NBVHEX,
     >            MESH  ,DMESH ,NBSURT,NBREGT,
     >            INDXSR(1,-NBSHEX),SURVOL(-NBSHEX))
*----
*  For 3-D case, displace hexagonal surfaces towards
*  the end of SURVOL leaving space for the possible
*  annular sub-surfaces (skip NR spaces for radial volumes).
*----
      ISURTF=NBSUR
      ISURTI=NBSURT
      IF(NDIM .EQ. 3) THEN
        DO IZ=1,2
          DO IT=1,NRP
            SURVOL(-ISURTF)=SURVOL(-ISURTI)
            SURVOL(-ISURTI)=DZERO
            DO IDIR=1,5
              INDXSR(IDIR,-ISURTF)=INDXSR(IDIR,-ISURTI)
              INDXSR(IDIR,-ISURTI)=0
            ENDDO
            ISURTF=ISURTF-NRP1
            ISURTI=ISURTI-1
          ENDDO
        ENDDO
      ENDIF
*----
*  Displace triangular volumes towards the end of vector
*  SURVOL leaving space for the possible annular sub-regions.
*----
      IVOLF=NBREG
      DO IVOLI=NBREGT,1,-1
        SURVOL(IVOLF)=SURVOL(IVOLI)
        SURVOL(IVOLI)=DZERO
        DO IDIR=1,5
          INDXSR(IDIR,IVOLF)=INDXSR(IDIR,IVOLI)
          INDXSR(IDIR,IVOLI)=0
        ENDDO
        IVOLF=IVOLF-NRP1
      ENDDO
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVolTri'
        WRITE(IOUT,6005) (ILVT,(INDXSR(IDIR,ILVT),IDIR=1,5),
     >                    SURVOL(ILVT),ILVT=-NBSUR,NBREG)
        WRITE(IOUT,6003)
      ENDIF
*----
*  Loop over radial regions for triangular region
*  intersection.
*----
      POSANN(1)=0.0
      POSANN(2)=0.0
      DO IR=NR,1,-1
        POSANN(0)=DMESH(IR,4)
*----
*  Loop over a single sector of triangle since centered annular
*  region
*----
        ISECT=1
        DO IT=1,NRTP
          ITYITA=NXTITA(POSTRI(1,1,IT,ISECT),POSANN,VOLINT)
          IF(ITYITA .NE. 0) THEN
*            WRITE(IOUT,*) 'IR,ISECT,IT,ITYITA,VOLINT',
*     >      IR,ISECT,IT,ITYITA,VOLINT
            IF(NDIM .EQ. 3) THEN
*----
*  Remove contribution from top and bottom surfaces
*----
              IOFZ=NSP*NZ
              DO IZ=1,2
                IOFS=(IT-1)*NRP1+IR
                DO IS1=1,6
                  ISURTF=IOFZ+IOFS
                  ISURTI=ISURTF+1
                  SURVOL(-ISURTI)=SURVOL(-ISURTI)-VOLINT
                  SURVOL(-ISURTF)=VOLINT
                  INDXSR(1,-ISURTF)=INDXSR(1,-ISURTI)
                  INDXSR(2,-ISURTF)=INDXSR(2,-ISURTI)
                  INDXSR(3,-ISURTF)=INDXSR(3,-ISURTI)
                  INDXSR(4,-ISURTF)=IR
                  INDXSR(5,-ISURTF)=INDXSR(5,-ISURTI)
                  IOFS=IOFS+NRP1*NRTP
                ENDDO
                IOFZ=IOFZ+NRTPP
              ENDDO
*----
*  Remove contributions to volumes
*----
              IOFZ=0
              DO IZ=1,NZ
                IOFS=(IT-1)*NRP1+IR
                VOLZ=VOLINT*(DMESH(IZ,3)-DMESH(IZ-1,3))
                DO IS1=1,6
                  ILVT=IOFZ+IOFS
                  ILVI=ILVT+1
                  SURVOL(ILVI)=SURVOL(ILVI)-VOLZ
                  SURVOL(ILVT)=VOLZ
                  INDXSR(1,ILVT)=INDXSR(1,ILVI)
                  INDXSR(2,ILVT)=INDXSR(2,ILVI)
                  INDXSR(3,ILVT)=INDXSR(3,ILVI)
                  INDXSR(4,ILVT)=IR
                  INDXSR(5,ILVT)=INDXSR(5,ILVI)
                  IOFS=IOFS+NRP1*NRTP
                ENDDO
                IOFZ=IOFZ+NRTPP
              ENDDO
            ELSE
              IOFS=(IT-1)*NRP1+IR
              VOLZ=VOLINT
              DO IS1=1,6
                ILVT=IOFS
                ILVI=ILVT+1
                SURVOL(ILVI)=SURVOL(ILVI)-VOLZ
                SURVOL(ILVT)=VOLZ
                INDXSR(1,ILVT)=INDXSR(1,ILVI)
                INDXSR(2,ILVT)=INDXSR(2,ILVI)
                INDXSR(3,ILVT)=INDXSR(3,ILVI)
                INDXSR(4,ILVT)=IR
                INDXSR(5,ILVT)=INDXSR(5,ILVI)
                IOFS=IOFS+NRP1*NRTP
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
*----
*  Print surfaces and volumes if required and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVolTriAnn'
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
 6010 FORMAT(1X,'MESH DIMENSIONS IN ',A1,' =',I10)
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F20.10)
      END
