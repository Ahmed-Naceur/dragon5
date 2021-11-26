*DECK NXTTGS
      SUBROUTINE NXTTGS(IPTRK ,IPRINT,NDIM  ,NBTDIR,NPOINT,NTRK  ,
     >                  ITRAK ,MAXMSH,NFSUR ,NFREG ,NUCELL,NBUCEL,
     >                  MXGSUR,MXGREG,MAXPIN,LINMAX,ITYPBC,IUNFLD,
     >                  MATALB,SURVOL,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                  DORITR,NSLINE,NCOR  ,WEIGHT,NUMERO,LENGTH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate a specific standard tracking line (isotropic tracking)
* for a geometry. This routine is used for line by line integration of
* the collision probability matrix without tracking file.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* IPRINT  print level.
* NDIM    problem dimensions.
* NBTDIR  number of tracks directions considered.
* NPOINT  number of integration points along each axis
*         in a plane mormal to track direction.
* NTRK    maximum number of track that can be generated.
* ITRAK   track number considered. For 3-D problems
*         ITRAK=(ITDIR-1)*(3*NPOINT**2)
*              +(IPLANE-1)*NPOINT**2
*              +(IPA2-1)*NPOINT
*              +IPA3
*         while for 2-D problems
*         ITRAK=(ITDIR-1)*(NPOINT)
*              +IPA3
*         will be used.
* MAXMSH  maximum number of elements in MESH array.
* NFSUR   number of surfaces.
* NFREG   number of regions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* NBUCEL  number of cells in unfolded geometry.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* MAXPIN  maximum number of pins in a cell.
* LINMAX  maximum number of segments in a track.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesian boundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* IUNFLD  description of unfolded geometry.
* MATALB  global mixture/albedo identification vector.
* SURVOL  global surface volume vector.
* DGMESH  meshing vector for global geometry.
* DANGLT  angles.
* DVNOR   ratio of analytic to tracked volume.
* DWGTRK  weight of tracking lines.
* DORITR  origin of tracking lines.
*
*Parameters: output
* NSLINE  number of segments for this track.
* NCOR    number of start/end surfaces.
* WEIGHT  weight associated with each line.
* NUMERO  region/surface identification number
*         for segment.
* LENGTH  segment length.
*
*References:
*  This routine represent a simplified version of NXTTLS where
*  the track origin, direction and weights are already known.
*  Moreover, instead of computing the track normalization
*  fators, one assume that they are available and the
*  track are renormalized directly in this module.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
*      INTEGER          IPTRK
      INTEGER          IPRINT,NDIM,NBTDIR,NPOINT,NTRK,ITRAK,MAXMSH,
     >                 NFSUR,NFREG,NUCELL(3),NBUCEL,MXGSUR,MXGREG,
     >                 MAXPIN,LINMAX,ITYPBC
      INTEGER          IUNFLD(2,NBUCEL),MATALB(-NFSUR:NFREG)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG),DGMESH(-1:MAXMSH,4),
     >                 DANGLT(NDIM,NBTDIR),DVNOR(NFREG)
      DOUBLE PRECISION DWGTRK(NBTDIR),
     >                 DORITR(NDIM*(NDIM+1),2*NDIM-3,NBTDIR)
      INTEGER          NSLINE,NCOR
      DOUBLE PRECISION WEIGHT
      INTEGER          NUMERO(LINMAX)
      DOUBLE PRECISION LENGTH(LINMAX)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTGS')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-8,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
      INTEGER          NXTLCA,NXTLHA,IRLA
*----
*  Local variables
*----
      INTEGER          ITDIR,NPLANE,IPLANE,IPTA2,IPTA3,IDIR,ITST,
     >                 NBCOR(2),NBSINT,ITLOC
      INTEGER          ITRN,ICEL,ICI,JLINE,IBLIN,IELIN,ISURF,
     >                 IX,IY,IZ,IOX,IOY,IOZ,IOC,LMAXT
      DOUBLE PRECISION ANGLES(3,3),
     >                 TRKORI(3),TRKOR2(3),TRKLIM(2),
     >                 CELLPO(3,2),DSTART,DERR
      INTEGER          ISEG,JSEG,IREG,ILREG,IPRINL,NBREG,IOFF
      INTEGER          NSDEB,NSFIN
      REAL             FACSC
      INTEGER          NICSS,ICISS,ICSR,ICSRR
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DLENGT,DCINT
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      NSLINE=0
      NCOR=0
      ITST=NTRK
      ITST=MATALB(0)
      FACSC=REAL(SURVOL(0))
      PI=XDRCST('Pi',' ')
      LMAXT=4*(NBUCEL+4)
      ITST=1
      NPLANE=2*NDIM-3
      IOC=0
*----
*  Scratch storage allocation
*   DLENGT  spatial location of each line segment.
*   ICINT   identification of spatial position for each
*           line segment in cell description of geometry.
*   DCINT   position of each intersection point for each
*           line segment in cell description of geometry.
*----
      ALLOCATE(ICINT(0:5,LMAXT))
      ALLOCATE(DLENGT(LINMAX),DCINT(LMAXT))
*----
*  Identify track direction plane number and point.
*  For 3-D:
*    ITRAK=(ITDIR-1)*(3*NPOINT**2)
*         +(IPLANE-1)*NPOINT**2
*         +(IPA2-1)*NPOINT
*         +IPA3
*  For 2-D:
*    ITRAK=(ITDIR-1)*(NPOINT)
*         +IPA2
*         will be used.
*----
      ITLOC=ITRAK
      IF(NDIM .EQ. 3) THEN
        IPTA3=MOD(ITLOC-1,NPOINT)+1
        ITLOC=ITLOC/NPOINT
        IPTA2=MOD(ITLOC,NPOINT)+1
        ITLOC=ITLOC/NPOINT
        IPLANE=MOD(ITLOC,NPLANE)+1
        ITLOC=ITLOC/NPLANE
        ITDIR=MOD(ITLOC,NBTDIR)+1
      ELSE
        IPTA3=1
        IPTA2=MOD(ITLOC-1,NPOINT)+1
        IPLANE=1
        ITLOC=ITLOC/NPOINT
        ITDIR=MOD(ITLOC,NBTDIR)+1
      ENDIF
      WEIGHT=DWGTRK(ITDIR)
*----
*  Find planes mormal to selected direction
*----
      CALL XDDSET(ANGLES,3*3,DZERO)
      DO IDIR=1,NDIM
        ANGLES(IDIR,1)=DANGLT(IDIR,ITDIR)
        TRKOR2(IDIR)=DORITR(IDIR,IPLANE,ITDIR)
      ENDDO
      IOFF=NDIM
      DO JLINE=1,NDIM
        DO IDIR=1,NDIM
          ANGLES(IDIR,JLINE)=DORITR(IDIR+IOFF,IPLANE,ITDIR)
        ENDDO
        IOFF=IOFF+NDIM
      ENDDO
*----
*  Position TRKORI with respect to IPTA2
*----
      DO IDIR=1,NDIM
        TRKORI(IDIR)=TRKOR2(IDIR)+DBLE(IPTA2)*ANGLES(IDIR,2)
      ENDDO
*----
*  Position TRKORI with respect to IPTA3
*----
      DO IDIR=1,NDIM
        TRKORI(IDIR)=TRKORI(IDIR)+DBLE(IPTA3)*ANGLES(IDIR,3)
      ENDDO
      IPRINL=IPRINT
      IRLA=-1
      IF(ITYPBC .EQ. 0) THEN
        IRLA=NXTLCA(IPRINL,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >               NUCELL,TRKORI,ANGLES,DGMESH,
     >               NBCOR ,NBSINT,ICINT ,DCINT)
      ELSE IF(ITYPBC .EQ. 1) THEN
*----
*  Annular assembly
*----
        CALL XABORT(NAMSBR//': Circular BC not implemented')
      ELSE IF(ITYPBC .EQ. 2) THEN
*----
*  Hexagonal assembly
*----
        IRLA=NXTLHA(IPRINL,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >              NUCELL,TRKORI,ANGLES,DGMESH,
     >              NBCOR ,NBSINT,ICINT ,DCINT)
      ENDIF
*----
*  When no external face crossed go to next line
*----
      IF(IRLA .EQ. -1) CALL XABORT(NAMSBR//
     >': This type of cell cannot be tracked by NXT:')
      IF(IRLA .EQ. 0) GO TO 115
*----
*  Test for multiple line segments in hexagonal assemblies
*
*----
      NICSS=1                  
      IF(ITYPBC .EQ. 2) THEN
        DO ICI=2,NBSINT
          IF(ICINT(0,ICI) .EQ. 0) THEN
            NICSS=NICSS+1                  
          ENDIF
        ENDDO
        ICSRR=NBCOR(1)+1
      ELSE
        ICSRR=NBCOR(1)+1
      ENDIF
*----
*  For each region crossed loop track geometry
*  present in this region
*----
      DO ICISS=1,NICSS
      ICSR=ICSRR
      IF(IPRINL .GE. 500) THEN
        WRITE(IOUT,6045)
        WRITE(IOUT,6040) ITRAK,ITDIR,IPLANE,
     >  IPTA2,IPTA3,WEIGHT
        WRITE(IOUT,6041) ITRAK,TRKORI
        WRITE(IOUT,6042) ITRAK,ANGLES
        WRITE(IOUT,6043) ITRAK,DCINT(1)
      ELSE IF(IPRINL .GE. 50) THEN
        WRITE(IOUT,6040) ITRAK,ITDIR,IPLANE,
     >  IPTA2,IPTA3,WEIGHT
      ENDIF
      IBLIN=1
      IELIN=0
      NUMERO(IBLIN)=0
      DLENGT(IBLIN)=DCINT(ICSR-1)
      DSTART=DCINT(ICSR-1)
      IBLIN=IBLIN+1
      DO ICI=ICSR,NBSINT-NBCOR(2)+1
        ICSRR=ICSRR+1
        IF(ITYPBC .EQ. 0) THEN
          IX=ICINT(1,ICI)
          IY=ICINT(2,ICI)
          IOX=IX
          CELLPO(1,2)=DGMESH(IX,1)
          CELLPO(1,1)=DGMESH(IX-1,1)
          CELLPO(2,2)=DGMESH(IY,2)
          CELLPO(2,1)=DGMESH(IY-1,2)
          IOY=(IY-1)*NUCELL(1)
          IOZ=0
          IF(NDIM .EQ. 3) THEN
            IZ=ICINT(3,ICI)
            IOZ=(IZ-1)*NUCELL(1)*NUCELL(2)
            CELLPO(3,2)=DGMESH(IZ,3)
            CELLPO(3,1)=DGMESH(IZ-1,3)
          ENDIF
          IOC=IOX+IOY+IOZ
        ELSE IF(ITYPBC .EQ. 2) THEN
          IOC=ICINT(0,ICI+1)
*----
*  For multiple track segment IOC=0 indicates that the current
*  segment is completed and that a new track segment should be
*  started at ICI+2
*----
          IF(IOC .EQ. 0) THEN
            ICSRR=ICSRR+2
            GO TO 125
          ENDIF        
          IOC=ICINT(0,ICI)
          IX=ICINT(1,ICI)
          IY=ICINT(1,ICI)
          CELLPO(1,2)=DGMESH(IX,1)
          CELLPO(1,1)=DGMESH(IX,1)
          CELLPO(2,2)=DGMESH(IY,2)
          CELLPO(2,1)=DGMESH(IY,2)
          IOZ=0
          IF(NDIM .EQ. 3) THEN
            IZ=ICINT(3,ICI)
            IOZ=(IZ-1)*NUCELL(1)
            CELLPO(3,2)=DGMESH(IZ,3)
            CELLPO(3,1)=DGMESH(IZ-1,3)
          ENDIF 
          IOC=IOZ+IX                  
        ENDIF
        ICEL=IUNFLD(1,IOC)
        ITRN=IUNFLD(2,IOC)
        TRKLIM(1)=DSTART
        DSTART=DSTART+DCINT(ICI)
        TRKLIM(2)=DSTART
        IF(ICI .EQ. NBCOR(1)+1) THEN
*----
*  initial surfaces (at TRKLIM(1)) considered
*----
          ISURF=-1
        ELSE IF(ICI .EQ. NBSINT-NBCOR(2)) THEN
*----
*  final surfaces (at TRKLIM(2)) considered
*----
          ISURF=1
        ELSE
*----
*  no surface considered
*----
          ISURF=0
        ENDIF
*----
*  Track turned Cell
*----
        CALL NXTTCR(IPTRK ,IPRINL,ICEL  ,ITRN  ,ISURF ,
     >              NDIM  ,MAXMSH,LINMAX,MXGSUR,MXGREG,
     >              MAXPIN,CELLPO,TRKLIM,TRKORI,ANGLES,
     >              IBLIN ,IELIN ,NUMERO,DLENGT)
        DERR=MAX(ABS(TRKLIM(1)),ABS(TRKLIM(2)))
        DERR=(DLENGT(IELIN)-DSTART)/DERR
        IF(DERR .GT. DCUTOF) THEN
          WRITE(IOUT,9100) NAMSBR,IELIN,
     >    DLENGT(IELIN),DSTART,DERR
          CALL XABORT(NAMSBR//
     >': End of track does not coincide with end of cell')
        ENDIF
        IBLIN=IELIN+1
      ENDDO
      NSLINE=IELIN
*----
*  Compress tracking vector for region with DLENGT=0.0 and
*  for successive segments in the same region.
*----
      JSEG=0
      ILREG=-1
      NBREG=0
      DO ISEG=1,NSLINE
        IREG=NUMERO(ISEG)
        IF(DLENGT(ISEG) .GT. DZERO) THEN
          IF(IREG .GT. 0) THEN
            IF(IREG .EQ. ILREG) THEN
              DLENGT(JSEG)=DLENGT(JSEG)+DLENGT(ISEG)
            ELSE
              JSEG=JSEG+1
              DLENGT(JSEG)=DLENGT(ISEG)
              NUMERO(JSEG)=IREG
              ILREG=IREG
            ENDIF
          ELSE
            JSEG=JSEG+1
            DLENGT(JSEG)=DLENGT(ISEG)
            NUMERO(JSEG)=IREG
            ILREG=-1
          ENDIF
        ENDIF
      ENDDO
      NSLINE=JSEG
*----
*  Add contribution of track to volume integration for this angle
*  in this quadrant and compress tracking line by removing
*  segments with NUMERO=0.
*----
      JSEG=0
      NSDEB=0
      NSFIN=0
      NBREG=0
      DO ISEG=1,NSLINE
        IREG=NUMERO(ISEG)
        IF(DLENGT(ISEG) .GT. DZERO) THEN
          IF(IREG .GT. 0) THEN
            NBREG=NBREG+1
            JSEG=JSEG+1
            LENGTH(JSEG)=DLENGT(ISEG)*DVNOR(IREG)
            NUMERO(JSEG)=IREG
          ELSE IF(IREG .LT. 0) THEN
            IF(NBREG .EQ. 0) THEN
              NSDEB=NSDEB+1
            ELSE
              NSFIN=NSFIN+1
            ENDIF
            JSEG=JSEG+1
            LENGTH(JSEG)=DLENGT(ISEG)
            NUMERO(JSEG)=IREG
          ENDIF
        ENDIF
      ENDDO
      NSLINE=JSEG
      NCOR=MAX(NSDEB,NSFIN)
      FACSC=1.0/REAL(NCOR)
      IF(NCOR .GT. 1) THEN
        IF(NSDEB .EQ. 2*NSFIN) THEN
*----
*  Duplicate final surfaces
*  And change their weight
*----
          DO ISEG=NSLINE-NSFIN+1,NSLINE
            LENGTH(ISEG)=LENGTH(ISEG)*FACSC
            LENGTH(ISEG+NSFIN)=LENGTH(ISEG)
            NUMERO(ISEG+NSFIN)=NUMERO(ISEG)
          ENDDO
          NSLINE=NSLINE+NSFIN
          DO ISEG=1,NSDEB
            LENGTH(ISEG)=LENGTH(ISEG)*FACSC
          ENDDO
        ELSE IF(2*NSDEB .EQ. NSFIN) THEN
*----
*  Displace tracks by NSDEB places
*  This automatically double the NSDEB first faces.
*----
          DO ISEG=NSLINE,1,-1
            LENGTH(ISEG+NSDEB)=LENGTH(ISEG)
            NUMERO(ISEG+NSDEB)=NUMERO(ISEG)
          ENDDO
          NSLINE=NSLINE+NSDEB
*----
*  Duplicate surface weights
*----
          DO ISEG=NSLINE-NSFIN+1,NSLINE
            LENGTH(ISEG)=LENGTH(ISEG)*FACSC
          ENDDO
          DO ISEG=1,NSFIN
            LENGTH(ISEG)=LENGTH(ISEG)*FACSC
          ENDDO
        ELSE IF(NSDEB .EQ. NSFIN) THEN
*----
*  Duplicate surface weights
*----
          DO ISEG=NSLINE-NSFIN+1,NSLINE
            LENGTH(ISEG)=LENGTH(ISEG)*FACSC
          ENDDO
          DO ISEG=1,NSDEB
            LENGTH(ISEG)=LENGTH(ISEG)*FACSC
          ENDDO
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Number of begin and end surfaces not compatible')
        ENDIF
      ENDIF
 125  CONTINUE
      ENDDO               
*----
*  Exit because line is outside circle or sphere surrounding geometry
*----
 115  CONTINUE
*----
*  Processing finished:
*  print track normalization vector.
*  and routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DCINT,DLENGT)
      DEALLOCATE(ICINT)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6040 FORMAT('Track',I10.10,'={',4(I10,','),
     >F15.8,'};')
 6041 FORMAT('Torig',I10.10,'={',2(F15.8,','),F15.8,'};')
 6042 FORMAT('Tdire',I10.10,'={',8(F15.8,','),F15.8,'};')
 6043 FORMAT('Tstrt',I10.10,'={',F15.8,'};')
 6045 FORMAT(1X)
 9100 FORMAT(1X,' ***** Error in ',A6,'***** for line ',I8/
     >       7X,'Positions (current and reference ) =',1P,2D21.14/
     >       7X,'Relative error = ',D21.14)
      END
