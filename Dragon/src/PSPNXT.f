*DECK PSPNXT
      SUBROUTINE PSPNXT(IPRINT,ISPSP ,ICOLR ,IPTRK ,ITYPBC,MAXMSH,
     >                  NDIM  ,NFSUR ,NFREG ,NUCELL,NBUCEL,
     >                  MXGREG,MAXPIN,COLREG,IUNFLD,MATALB,DGMESH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the graphics for a 2-D NXT geometry.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPRINT  print level.
* ISPSP   pointer to the POSTSCRIPT file.
* ICOLR   color set used where:
*         =-4 HSB filling with no contour;
*         =-3 CYMK filling with no contour;
*         =-2 RGB filling with no contour;
*         =-1 BW filling with no contour;
*         = 0 no filling with contour;
*         = 1 BW filling with contour;
*         = 2 RGB filling with contour;
*         = 3 CMYK filling with contour;
*         = 4 HSB filling with  contour.
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* ITYPBC  type of cell boundary.
* MAXMSH  maximum number of elements in MESH array.
* NDIM    dimension of the problem.
* NFSUR   number of surfaces.
* NFREG   number of regions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* NBUCEL  number of cells in unfolded geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* COLREG  region color.
* IUNFLD  description of unfolded geometry.
* MATALB  global mixture/albedo identification vector.
* DGMESH  meshing vector for global geometry.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPRINT,ISPSP,ICOLR,ITYPBC,MAXMSH,
     >                 NDIM,NFSUR,NFREG,NUCELL(3),NBUCEL,
     >                 MXGREG,MAXPIN
      REAL             COLREG(4,NFREG)
      INTEGER          IUNFLD(2,NBUCEL),
     >                 MATALB(-NFSUR:NFREG)
      DOUBLE PRECISION DGMESH(0:MAXMSH,4)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPNXT')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
      DOUBLE PRECISION DIMX,DIMY
      PARAMETER       (DIMX=3.5D0,DIMY=3.5D0)
*----
*  Local variables
*----
      INTEGER          KPSP(7)
      CHARACTER        NAMREC*12
      INTEGER          IDIR,ILPD,IX,IY,ICELL,ICEL,ITRN,ILONG,ITYLCM
      DOUBLE PRECISION RCIRC,ABSC(2),OFFC(2),FACT,CELLPO(2,2)
      REAL             XYPOS(2,4)
      DOUBLE PRECISION SIDEH,CENTH,DHMAX
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDREG,ITPIN,NBPTS,REGI,EVENT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DCMESH,DRAPIN,
     > POSTRI,COOR
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
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Initialize ICOL for color treatment
*  and ICONT for contour
*  KPSP(1)=ICONT
*  KPSP(2)=ICOL
*  KPSP(3)=IWLFAC (0->1.0,1->2.5)
*  KPSP(4)=KFS
*  KPSP(5)=KFR
*  KPSP(6)=KSS
*  KPSP(7)=KSR
*-----
      ILPD=MATALB(0)
      KPSP(1)=1
      KPSP(2)=ABS(ICOLR)
      KPSP(3)=0
      KPSP(4)=0
      KPSP(5)=0
      KPSP(6)=0
      KPSP(7)=0
      IF(ICOLR .EQ. 0) THEN
        KPSP(3)=1
      ELSE IF(ICOLR .LT. 0) THEN
        KPSP(1)=0
      ELSE
        KPSP(4)=1
        KPSP(7)=1
      ENDIF
      RCIRC=1.0D0
*----
*  Read global mesh for geometry
*  and determine graphics size
*----
      IF(ITYPBC .EQ. 0) THEN
*----
*  Cartesian
*----
        DO IDIR=1,NDIM
          NAMREC='G00000001SM'//CDIR(IDIR)
          ILPD=NUCELL(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          ABSC(IDIR)=0.5D0*(DGMESH(ILPD,IDIR)-DGMESH(0,IDIR))
          OFFC(IDIR)=0.5D0*(DGMESH(ILPD,IDIR)+DGMESH(0,IDIR))
        ENDDO
        RCIRC=DZERO
        DO IDIR=1,NDIM
          RCIRC=MAX(RCIRC,ABSC(IDIR))
        ENDDO
      ELSE IF(ITYPBC .EQ. 1) THEN
*----
*  Annular
*----
        DO IDIR=1,NDIM
          NAMREC='G00000001SM'//CDIR(IDIR)
          ILPD=NUCELL(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          ABSC(IDIR)=0.5D0*(DGMESH(ILPD,IDIR)-DGMESH(0,IDIR))
          OFFC(IDIR)=0.5D0*(DGMESH(ILPD,IDIR)+DGMESH(0,IDIR))
        ENDDO
        IDIR=4
        NAMREC='G00000001SM'//CDIR(IDIR)
        CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
        RCIRC=DGMESH(1,IDIR)
      ELSE IF(ITYPBC .EQ. 2) THEN
*----
*  Hexagonal
*----
        DO IDIR=1,2
          NAMREC='G00000001SM'//CDIR(IDIR)
          CALL LCMGET(IPTRK,NAMREC,DGMESH(0,IDIR))
          SIDEH=DGMESH(0,IDIR)
          CENTH=DGMESH(1,IDIR)
          OFFC(IDIR)=CENTH
          DHMAX=DZERO
          DO ICELL=2,NUCELL(IDIR)
            DHMAX=MAX(DHMAX,ABS(DGMESH(ICELL,IDIR)-CENTH))
          ENDDO
          ABSC(IDIR)=DHMAX+SIDEH
        ENDDO
        RCIRC=DZERO
        DO IDIR=1,NDIM
          RCIRC=MAX(RCIRC,ABSC(IDIR))
        ENDDO
      ELSE
        CALL XABORT(NAMSBR//
     >  ': Invalid geometry boundary types for PSP')
      ENDIF
*----
*  Locate pen at center of page
*----
      XYPOS(1,1)=DIMX
      XYPOS(2,1)=DIMY
      CALL PSMOVE(ISPSP,XYPOS,-3)
      FACT=DIMX/RCIRC
      ALLOCATE(IDREG(MXGREG),ITPIN(3*(MAXPIN)))
      ALLOCATE(DCMESH(4*(MAXMSH+2)),DRAPIN(6*(MAXPIN)))
*----
*  Scan over all Cartesian cells
*  1) Mesh in $Y$ direction
*----
      IF(ITYPBC .EQ. 0) THEN
        ICELL=0
        DO IY=1,NUCELL(2)
          CELLPO(2,1)=(DGMESH(IY-1,2)-OFFC(2))
          CELLPO(2,2)=(DGMESH(IY,2)-OFFC(2))
*----
*  2) Mesh in $X$ direction
*----
          DO IX=1,NUCELL(1)
            CELLPO(1,1)=(DGMESH(IX-1,1)-OFFC(1))
            CELLPO(1,2)=(DGMESH(IX,1)-OFFC(1))
            ICELL=ICELL+1
            ICEL=IUNFLD(1,ICELL)
            ITRN=IUNFLD(2,ICELL)
            IF(ITRN .EQ. 1) THEN
              CALL PSPTCR(IPTRK ,ISPSP ,IPRINT,ICEL  ,NDIM  ,NFREG ,
     >                    MAXMSH,MXGREG,MAXPIN,KPSP  ,COLREG,FACT  ,
     >                    CELLPO,IDREG ,ITPIN ,DCMESH,DRAPIN)
            ENDIF
          ENDDO
        ENDDO
      ELSE IF(ITYPBC .EQ. 2) THEN
        ALLOCATE(NBPTS(MXGREG),POSTRI(2*4*MXGREG))
        DO ICELL=1,NUCELL(1)
*----
*  Scan over all hexagonal cells
*----
          CELLPO(2,1)=(DGMESH(ICELL,2)-OFFC(2))
          CELLPO(2,2)=(DGMESH(ICELL,2)-OFFC(2))
          CELLPO(1,1)=(DGMESH(ICELL,1)-OFFC(1))
          CELLPO(1,2)=(DGMESH(ICELL,1)-OFFC(1))
          ICEL=IUNFLD(1,ICELL)
          ITRN=IUNFLD(2,ICELL)
          IF(ITRN .EQ. 1) THEN
            CALL PSPTHR(IPTRK ,ISPSP ,IPRINT,ICEL  ,NDIM  ,NFREG ,
     >                  MAXMSH,MXGREG,MAXPIN,KPSP  ,COLREG,FACT  ,
     >                  CELLPO,IDREG ,ITPIN ,
     >                  DCMESH,DRAPIN,NBPTS ,POSTRI)
          ENDIF
        ENDDO
        DEALLOCATE(POSTRI,NBPTS)
      ENDIF 
      DEALLOCATE(DRAPIN,ITPIN,DCMESH,IDREG)
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  plot MC: neutron paths if present
*----
      CALL LCMLEN(IPTRK,'MCpoints',ILONG,ITYLCM)
      IF (ITYLCM.EQ.0) THEN
         CALL LCMSIX(IPTRK,'MCpoints',1)
         CALL LCMLEN(IPTRK,'REGI',ILONG,ITYLCM)
         ALLOCATE(REGI(ILONG),EVENT(ILONG))
         ALLOCATE(COOR(3*ILONG))
         CALL LCMGET(IPTRK,'COORD',COOR)
         CALL LCMGET(IPTRK,'REGI',REGI)
         CALL LCMGET(IPTRK,'EVENT',EVENT)
         CALL PSPMCP(ISPSP,OFFC,FACT,ILONG,COOR,REGI,EVENT)
         DEALLOCATE(EVENT,REGI) 
         DEALLOCATE(COOR)
         CALL LCMSIX(IPTRK,' ',2)
      ENDIF
*----
*  Save track normalisation vector
*----
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
