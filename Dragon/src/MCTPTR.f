*DECK MCTPTR
      SUBROUTINE MCTPTR(IPTRK,IPRINT,NDIM,MAXMSH,ITYPBC,NUCELL,MXGSUR,
     1                  MXGREG,MAXPIN,IUNFLD,DGMESH,XYZL,NBIND,POS,
     2                  LENGTH,VDIR,ODIR,IDS,IDSO,IREG,INDX,IDIRC,
     3                  MESHC,NSURC,NREGC,NTPIN,CELLPO,PINCEN,INDEX,
     4                  IDREG,DCMESH,ITPIN,DRAPIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find region/surface index from an initial position and a path to
* travel.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NDIM    problem dimensions.
* MAXMSH  maximum number of elements in MESH array.
* ITYPBC  type of boundary.
* NUCELL  number of cell after unfolding in 
*         $X$, $Y$ and $Z$ directions.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* IUNFLD  description of unfolded geometry.
* DGMESH  meshing vector for global geometry.
* XYZL    undefined.
* NBIND   first dimension of INDX.
* VDIR    travel direction (unit vector).
*
*Parameters: output
* IDS     outer surface orientation index. ABS(IDSO).
* IDSO    outer surface oreientation index with +-1:X+-; +-2:Y+-;
*         +-3:Z+- faces.
* IREG    region/surface index.
* ODIR    search (octant) direction. 
*
*Parameters: input/output
* POS     initial/final position.
* LENGTH  length to travel/remaining length on the path.
* INDX    location index of the initial/final position in the 
*         geometry structure.
*
*Parameters: scratch
* IDIRC   undefined.
* MESHC   undefined.
* NSURC   undefined.
* NREGC   undefined.
* NTPIN   undefined.
* CELLPO  undefined.
* PINCEN  undefined.
* INDEX   undefined.
* IDREG   undefined.
* DCMESH  undefined.
* ITPIN   undefined.
* DRAPIN  undefined.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER IPRINT,NDIM,MAXMSH,ITYPBC,NUCELL(3),MXGSUR,MXGREG,MAXPIN,
     1 IUNFLD(2,NUCELL(1),NUCELL(2),NUCELL(3)),NBIND,ODIR(3),IREG,
     2 INDX(NBIND,0:2),IDIRC(2),MESHC(4,2),NSURC(2),NREGC(2),NTPIN,
     3 INDEX(5,-MXGSUR:MXGREG,2),IDREG(MXGREG,2),ITPIN(3,MAXPIN)
      DOUBLE PRECISION DGMESH(-1:MAXMSH,4),XYZL(2,NDIM),POS(3),LENGTH,
     1 VDIR(3),CELLPO(3,2),PINCEN(3),DCMESH(-1:MAXMSH,4,2),
     2 DRAPIN(-1:4,MAXPIN)
*----
*  LOCAL VARIABLES
*----
      INTEGER IDIR,IBCO(3),IDS,IDSO,JJ
      DOUBLE PRECISION VCOR(3),LEN,LENM
      DOUBLE PRECISION EPS
      PARAMETER(EPS=1.D-8)
      INTEGER INDOS(2,3)
      DATA INDOS / 2,3,
     1             3,1,
     2             1,2 /
*---- 
*  VERIFY IF A BOUNDARY IS REACHED
*----
      IDSO=0
      LENM=0.0D0
      IF(ITYPBC.EQ.0) THEN
*       CARTESIAN BOUNDARY
*       find corresponding geometry corner corresponding to travel direction 
*       check if a boundary is reached and if so, which one it is
        IDS=0
        LENM=LENGTH
        DO IDIR=1,NDIM
          IF(VDIR(IDIR).GT.EPS) THEN
            ODIR(IDIR)=1
            IBCO(IDIR)=2
            VCOR(IDIR)=XYZL(2,IDIR)-POS(IDIR)
            LEN=VCOR(IDIR)/VDIR(IDIR)
            IF(LEN.LT.LENM) THEN
               LENM=LEN
               IDS=IDIR
               IF(IDIR.EQ.1) THEN
                 IREG=-2 ;
               ELSEIF(IDIR.EQ.2) THEN
                 IREG=-4 ;
               ELSEIF(IDIR.EQ.3) THEN
                 IREG=-6 ;
               ENDIF
            ENDIF
          ELSEIF(VDIR(IDIR).LT.-EPS) THEN
            ODIR(IDIR)=-1
            IBCO(IDIR)=1
            VCOR(IDIR)=XYZL(1,IDIR)-POS(IDIR)
            LEN=VCOR(IDIR)/VDIR(IDIR)
            IF(LEN.LT.LENM) THEN
               LENM=LEN
               IDS=IDIR
               IF(IDIR.EQ.1) THEN
                 IREG=-1 ;
               ELSEIF(IDIR.EQ.2) THEN
                 IREG=-3 ;
               ELSEIF(IDIR.EQ.3) THEN
                 IREG=-5 ;
               ENDIF
            ENDIF
          ELSE
            ODIR(IDIR)=0
          ENDIF
        ENDDO
      ELSE
*       CYLINDRICAL BOUNDARY
        CALL XABORT('MCTPTR: CYLINDRICAL/HEXAGONAL BOUNDARY NOT IMPLEM'
     1  //'ENTED YET.')
      ENDIF
      IF(IDS.EQ.0) THEN
*----
*  NO: LOCATE POINT WITHIN THE GEOMETRY
*----
        DO IDIR=1,NDIM
          POS(IDIR)=POS(IDIR)+LENGTH*VDIR(IDIR)
        ENDDO
        DO IDIR=NDIM+1,3
          POS(IDIR)=1.D0
        ENDDO
!!!!!! for the time being it is the same routine as for a starting point
!!!!!! an optimized version should use the info of the previous point 
!!!!!! to start its search
        CALL MCTCTR(IPTRK,IPRINT,NDIM,MAXMSH,NUCELL,MXGSUR,MXGREG,
     1     MAXPIN,IUNFLD,DGMESH,NBIND,ODIR,POS,IREG,INDX,IDIRC,MESHC,
     2     NSURC,NREGC,NTPIN,CELLPO,PINCEN,INDEX,IDREG,DCMESH,ITPIN,
     3     DRAPIN)
      ELSE
*----
*  YES: LOCATE POINT ON THE BOUNDARY
*----
        POS(IDS)=XYZL(IBCO(IDS),IDS)
        DO JJ=1,2
          IDIR=INDOS(JJ,IDS)
          POS(IDIR)=POS(IDIR)+LENM*VDIR(IDIR)
        ENDDO
        LENGTH=LENGTH-LENM
        IDSO=IDS
        IF(IBCO(IDS).EQ.1) IDSO=-IDSO
      ENDIF
*
      RETURN
      END
