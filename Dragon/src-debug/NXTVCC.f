*DECK NXTVCC
      SUBROUTINE NXTVCC(IPRINT,NDIM  ,IDIRCX,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DMESH ,NBSUR ,NBREG ,INDXSR,SURVOL)
*
*----------
*
*Purpose:
* Compute the volume of each region for a mixed annular/Cartesian
* 2-D or 3-D geometry using the NXT tracking procedure.
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
* IDIRCX  the direction of the first axis of an annular geometry
*         assuming the axis are in a cyclic rotation.
*         A negative value means that the external boundary is
*         annular while a positive boundary implies that the
*         external boundaries are Cartesian.
* MXMESH  maximum number of spatial subdivision in
*         $R$ and $X$, $Y$ or $Z$.
* MAXSUR  maximum number of surfaces in the geometry.
* MAXREG  maximum number of regions in the geometry.
* MESH    effective number of spatial subdivision in $R$
*         and $X$, $Y$ or $Z$.
* DMESH   spatial description of the cylinder.
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
*  1- Contents of IDIRCX:
*     IDIRCX      Annulus in 2-D plane Cylinder directions in 3-D
*     +/- 1       (x,y)                z
*     +/- 2       (y,z)                x
*     +/- 3       (z,x)                y
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
*       INDXSR(4,i)= ir is the $R$ location of region i.
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
      INTEGER          IPRINT,NDIM,IDIRCX,MXMESH,MAXSUR,MAXREG
      INTEGER          MESH(4)
      DOUBLE PRECISION DMESH(-1:MXMESH,4)
      INTEGER          NBSUR,NBREG,INDXSR(5,-MAXSUR:MAXREG)
      DOUBLE PRECISION SURVOL(-MAXSUR:MAXREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTVCC')
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
      INTEGER          NXTIRA,ITYIRA
      DOUBLE PRECISION VOLINT
*----
*  Local variables
*----
      INTEGER          NBVCAR,ID,IDG,IDIR(MAXDIM),
     >                 NM(MAXDIM),IDM(MAXDIM)
      INTEGER          IDIRC,IDGP1,IDGP2,IDGPP,NBSCAR,NSKPI,NSKPF
      INTEGER          NBSURT,NBREGT,NRP1
      INTEGER          ISURTF,ISURBF,ISURTI,ISURBI,IS,IN,IVOLF,IVOLI
      INTEGER          IR,IX,IY,IZ,IXYZ,ISURTN,ISURBN,ILVI,ILVT,IOF
      INTEGER          IDX,IDY,IDZ,IDR,NOFSUR
      DOUBLE PRECISION XYCAR(4),POSANN(0:2),DZ
      INTEGER          IPRNT2
*----
*  Data
*----
      CHARACTER        CDIR(MAXDIM)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Prepare radial and axial loops over spatial directions
*  as a function of IDIRC and NDIM.
*----
      DZ=DONE
      IDIRC=ABS(IDIRCX)
      IF(IDIRC .EQ. 1) THEN
        IDX=1
        IDY=2
        IDZ=3
      ELSE IF(IDIRC .EQ. 2) THEN
        IDX=2
        IDY=3
        IDZ=1
      ELSE
        IDX=3
        IDY=1
        IDZ=2
      ENDIF
      IDR=4
      IPRNT2=IPRINT/2
      PI=XDRCST('Pi',' ')
      IF(NDIM .EQ. 1) CALL XABORT(NAMSBR//
     >': Only 2-D and 3-D problems permitted')
      CALL XDDSET(SURVOL,(MAXSUR+MAXREG+1),DZERO)
*----
*  Prepare loops over spatial directions as a function
*  of IDIRC and NDIM.
*  Compute number of Cartesian surfaces.
*----
      NBVCAR=1
      DO 100 ID=1,NDIM
        IDG=MOD(IDIRC+ID-2,3)+1
        IDIR(ID)=IDG
        NM(IDG)=MESH(IDG)
        NBVCAR=NBVCAR*NM(IDG)
        IDM(IDG)=1
 100  CONTINUE
      DO 101 ID=NDIM+1,3
        IDG=MOD(IDIRC+ID-2,3)+1
        IDIR(ID)=IDG
        NM(IDG)=1
        IDM(IDG)=0
 101  CONTINUE
      IDG=4
      IDIR(4)=IDG
      NM(IDG)=MESH(IDG)
      IDM(IDG)=1
      NBREG=NBVCAR*NM(IDG)
      NRP1=NM(4)
      IF(IDIRCX .GT. 0) THEN
        NRP1=NRP1+1
        NBREG=NBREG+NBVCAR
      ENDIF
      IDGP1=IDIR(1)
      IDGP2=IDIR(2)
      IDGPP=IDIR(3)
      IF(MAXREG .LT. NBREG) CALL XABORT(NAMSBR//
     >': Insufficient space to store region volumes')
*----
*  Compute number of Cartesian surfaces
*  1- Surface parallel to cylinder axis
*----
      IF(IDIRCX .GT. 0) THEN
        NBSCAR=0
        NBSUR=0
        DO 102 ID=1,2
          IDG=IDIR(ID)
          NBSCAR=NBSCAR+2*NBVCAR/NM(IDG)
          NBSUR=NBSUR+2*NBVCAR/NM(IDG)
 102    CONTINUE
      ELSE
        NBSUR=0
        NBSCAR=0
      ENDIF
*----
*  2- Surface normal to cylinder axis (if any)
*----
      NSKPI=0
      NSKPF=0
      DO 103 ID=3,NDIM
        IDG=IDIR(ID)
        NSKPI=NBVCAR/NM(IDG)
        NSKPF=NM(4)*NSKPI
        IF(IDIRCX .GT. 0) THEN
          NSKPF=NSKPF+NSKPI
        ENDIF
        NBSCAR=NBSCAR+2*NSKPI
        NBSUR=NBSUR+2*NSKPF
 103  CONTINUE
      IF(MAXSUR .LT. NBSUR) CALL XABORT(NAMSBR//
     >': Insufficient space to store surface areas')
*----
*  Print mesh if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6002) 'CENTER'//CDIR(IDGP1)//CDIR(IDGP2)
        WRITE(IOUT,6006) (DMESH(-1,IDIR(ID)),ID=1,2)
        WRITE(IOUT,6003)
        WRITE(IOUT,6002) 'RADIAL'
        WRITE(IOUT,6006) (DMESH(IR,4),IR=1,MESH(4))
        WRITE(IOUT,6003)
        DO 600 ID=1,NDIM
          IDG=IDIR(ID)
          WRITE(IOUT,6002) 'MESH'//CDIR(IDG)
          WRITE(IOUT,6006) (DMESH(IXYZ,IDG),IXYZ=0,NM(IDG))
          WRITE(IOUT,6003)
 600    CONTINUE
      ENDIF
      NOFSUR=NBSUR
      IF(IDIRCX .GT. 0) THEN
*----
*  Compute volumes and surfaces associated with
*  Cartesian regions.
*----
        CALL NXTVCA(IPRNT2,NDIM  ,IDIRC ,MXMESH,
     >              NBSCAR,NBVCAR,MESH  ,DMESH ,
     >              NBSURT,NBREGT,
     >              INDXSR(1,-NBSCAR),SURVOL(-NBSCAR))
*----
*  For 3-D case, displace Cartesian surfaces towards
*  the end of SURVOL leaving space for the possible
*  annular sub-surfaces.
*----
        IF(NDIM .EQ. 3) THEN
          ISURTF=NBSUR
          ISURTI=NBSURT
          IF(IDIRC .GT. 1) THEN
*----
*  Displace Z faces
*----
            DO IZ=1,2
              DO IY=1,NM(2)
                DO IX=1,NM(1)
                  SURVOL(-ISURTF)=SURVOL(-ISURTI)
                  SURVOL(-ISURTI)=DZERO
                  DO IN=1,5
                    INDXSR(IN,-ISURTF)=INDXSR(IN,-ISURTI)
                    INDXSR(IN,-ISURTI)=0
                  ENDDO
                  ISURTF=ISURTF-1
                  ISURTI=ISURTI-1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF(IDIRC .EQ. 2) THEN
*----
*  Displace Y faces
*----
            DO IY=1,2
              DO IX=1,NM(1)
                DO IZ=1,NM(3)
                  SURVOL(-ISURTF)=SURVOL(-ISURTI)
                  SURVOL(-ISURTI)=DZERO
                  DO IN=1,5
                    INDXSR(IN,-ISURTF)=INDXSR(IN,-ISURTI)
                    INDXSR(IN,-ISURTI)=0
                  ENDDO
                  ISURTF=ISURTF-1
                  ISURTI=ISURTI-1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
*----
*  Displace (X (2), Y (3) , or Z (1)) and leave space for annular
*----
          NOFSUR=ISURTF
          DO IS=NSKPI,1,-1
            SURVOL(-ISURTF)=SURVOL(-ISURTI)
            SURVOL(-ISURTI)=DZERO
            DO IN=1,5
              INDXSR(IN,-ISURTF)=INDXSR(IN,-ISURTI)
              INDXSR(IN,-ISURTI)=0
            ENDDO
            ISURTF=ISURTF-NRP1
            ISURTI=ISURTI-1
          ENDDO
          DO IS=NSKPI,1,-1
            SURVOL(-ISURTF)=SURVOL(-ISURTI)
            SURVOL(-ISURTI)=DZERO
            DO IN=1,5
              INDXSR(IN,-ISURTF)=INDXSR(IN,-ISURTI)
              INDXSR(IN,-ISURTI)=0
            ENDDO
            ISURTF=ISURTF-NRP1
            ISURTI=ISURTI-1
          ENDDO
        ENDIF
*----
*  Displace Cartesian volumes towards the end of vector
*  SURVOL leaving space for the possible annular sub-regions.
*----
        IVOLF=NBREG
        DO IVOLI=NBREGT,1,-1
          SURVOL(IVOLF)=SURVOL(IVOLI)
          SURVOL(IVOLI)=DZERO
          DO IN=1,5
            INDXSR(IN,IVOLF)=INDXSR(IN,IVOLI)
            INDXSR(IN,IVOLI)=0
          ENDDO
          IVOLF=IVOLF-NRP1
        ENDDO
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6002) 'SurVol'
          WRITE(IOUT,6005) (ILVT,(INDXSR(IR,ILVT),IR=1,5),SURVOL(ILVT),
     >                       ILVT=-NBSUR,NBREG)
          WRITE(IOUT,6003)
        ENDIF
      ENDIF
*----
*  Loop over radial regions for Cartesian/annular region
*  intersection.
*----
      ISURTF=NOFSUR-NSKPF+1
      ISURBF=NOFSUR-2*NSKPF+1
      DO 120 IR=NM(4),1,-1
        POSANN(0)=DMESH(IR,4)
        POSANN(1)=DMESH(-1,IDGP1)
        POSANN(2)=DMESH(-1,IDGP2)
*----
*  Loop over second normal direction
*----
        DO 121 IY=1,NM(IDGP2)
          XYCAR(3)=DMESH(IY-1,IDGP2)
          XYCAR(4)=DMESH(IY,IDGP2)
*----
*  Loop over first normal direction
*----
          DO 122 IX=1,NM(IDGP1)
            XYCAR(1)=DMESH(IX-1,IDGP1)
            XYCAR(2)=DMESH(IX,IDGP1)
*----
*  Rectangle/annular region intersection
*----
            ITYIRA=NXTIRA(XYCAR,POSANN,VOLINT)
            IF(ITYIRA .NE. 0) THEN
              IF(NDIM .EQ. 3) THEN
*----
*  For 3-D problem when
*  rectangle and annular regions intersect:
*  Correct top and bottom surfaces
*----
                ISURTI=ISURTF+NRP1*(IX-1+NM(IDGP1)*(IY-1))
                ISURBI=ISURBF+NRP1*(IX-1+NM(IDGP1)*(IY-1))
                ISURTN=ISURTI+IR-1
                ISURTI=ISURTI+IR
                ISURBN=ISURBI+IR-1
                ISURBI=ISURBI+IR
                IF(IR .NE. NM(4) .OR. IDIRCX .GT. 0) THEN
                  SURVOL(-ISURTI)=SURVOL(-ISURTI)-VOLINT
                  SURVOL(-ISURBI)=SURVOL(-ISURBI)-VOLINT
                ENDIF
                SURVOL(-ISURTN)=VOLINT
                SURVOL(-ISURBN)=VOLINT
                INDXSR(IDX,-ISURTN)=IX
                INDXSR(IDX,-ISURBN)=IX
                INDXSR(IDY,-ISURTN)=IY
                INDXSR(IDY,-ISURBN)=IY
                INDXSR(IDZ,-ISURTN)=-2
                INDXSR(IDZ,-ISURBN)=-1
                INDXSR(IDR,-ISURTN)=IR
                INDXSR(IDR,-ISURBN)=IR
              ENDIF
*----
*  2- Volumes
*----
              DO 124 IZ=1,NM(IDGPP)
                IF(IDIRC .EQ. 1) THEN
                  IOF=NRP1*(IX-1+NM(IDGP1)*((IY-1)+(IZ-1)*NM(IDGP2)))
                ELSE IF(IDIRC .EQ. 2) THEN
                  IOF=NRP1*(IZ-1+NM(IDGPP)*((IX-1)+(IY-1)*NM(IDGP1)))
                ELSE
                  IOF=NRP1*(IY-1+NM(IDGP2)*((IZ-1)+(IX-1)*NM(IDGPP)))
                ENDIF
                ILVT=IR+IOF
                ILVI=ILVT+1
                DZ=DONE
                IF(IDM(IDGPP) .EQ. 1)
     >          DZ=DMESH(IZ,IDGPP)-DMESH(IZ-1,IDGPP)
                IF(IR .NE. NM(4) .OR. IDIRCX .GT. 0) THEN
                  SURVOL(ILVI)=SURVOL(ILVI)-VOLINT*DZ
                ENDIF
                SURVOL(ILVT)=VOLINT*DZ
                INDXSR(IDX,ILVT)=IX
                INDXSR(IDY,ILVT)=IY
                INDXSR(IDZ,ILVT)=IZ
                INDXSR(IDR,ILVT)=IR
 124          CONTINUE
            ENDIF
 122      CONTINUE
 121    CONTINUE
 120  CONTINUE
      IF(IDIRCX .LT. 0) THEN
*----
*  Add radial surfaces
*----
        ISURTN=NBSUR
        DO 130 IZ=1,NM(IDGPP)
          ISURTN=ISURTN+1
          IF(IDM(IDGPP) .EQ. 1) DZ=DMESH(IZ,IDGPP)-DMESH(IZ-1,IDGPP)
          SURVOL(-ISURTN)=DZ*DTWO*PI*DMESH(NM(4),4)
          INDXSR(IDX,-ISURTN)=0
          INDXSR(IDY,-ISURTN)=0
          INDXSR(IDZ,-ISURTN)=IZ
          INDXSR(IDR,-ISURTN)=-2
 130    CONTINUE
        NBSUR=ISURTN
      ENDIF
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVol'
        WRITE(IOUT,6005) (ILVT,(INDXSR(IR,ILVT),IR=1,5),SURVOL(ILVT),
     >                     ILVT=-NBSUR,NBREG)
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
 6005 FORMAT((6(I10,','),D20.10,:,','))
 6006 FORMAT(4(F20.10,:,','))
      END
