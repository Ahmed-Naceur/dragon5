*DECK NXTTLS
      SUBROUTINE NXTTLS(IPTRK ,IFTEMP,IPRINT,IGTRK ,NDIM  ,MAXMSH,
     >                  NFSUR ,NFREG ,NUCELL,NBUCEL,NBANGL,NQUAD ,
     >                  NPLANE,NPOINT,LINMAX,MXGSUR,MXGREG,RENO  ,
     >                  MAXPIN,NBTDIR,NBDR  ,ITYPBC,IFMT  ,
     >                  RCUTOF,SPACLN,WEIGHT,RADIUS,CENTER,
     >                  IUNFLD,SURVOL,DGMESH,DANGLT,DDENWT,
     >                  MAXSGL,NTLINE,DVNOR ,DSNOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the standard tracking lines (isotropic tracking)
* for a geometry.
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
* IFTEMP  pointer to a temporary TRACKING data structure in
*         update or creation mode.
* IPRINT  print level.
* IGTRK   flag to generate the tracking file. In the case where
*         IGTRK=1, the tracking is performed and
*         used to evaluate the track normalisation factor and the
*         tracking file is generated. When IGTRK=0, the tracking is
*         still performed and used to evaluate the
*         track normalisation factor but the tracking file is not
*         generated.
* NDIM    problem dimensions.
* MAXMSH  maximum number of elements in MESH array.
* NFSUR   number of surfaces.
* NFREG   number of regions.
* NUCELL  number of cell after unfolding in
*         $X$, $Y$ and $Z$ directions.
* NBUCEL  number of cells in unfolded geometry.
* NBANGL  number of angles.
* NQUAD   number of quadrant (in 3-D) and quarter (in 2-D).
* NPLANE  number of normal planes considered.
* NPOINT  number of integration points along each axis
*         in a plane mormal to track direction.
* LINMAX  maximum number of segments in a track.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* RENO    track normalisation option. A value RENO=-1 implies
*         a direction dependent normalization of the tracks
*         for the volume while a value RENO=0, implies
*         a global normalisation.
* MAXPIN  maximum number of pins in a cell.
* NBTDIR  number of tracks directions considered.
* NBDR    number of directions for track normalization.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesianb oundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* IFMT    tracking file format:
*         =0 for short file;
*         =1 long file required by TLM:.
* RCUTOF  corner cutoff.
* SPACLN  linear track spacing in the plane.
* WEIGHT  weight associated with each line in the plane.
* RADIUS  radius of circle (2-D) or sphere (3-D) surrounding
*         the geometry.
* CENTER  center of circle (2-D) or sphere (3-D) surrounding
*         the geometry.
* IUNFLD  description of unfolded geometry.
* SURVOL  global surface volume vector.
* DGMESH  meshing vector for global geometry.
* DANGLT  angles.
* DDENWT  angular density for each angle.
*
*Parameters: output
* MAXSGL  maximum number of segments in a line.
* NTLINE  total number of lines generated.
* DVNOR   ratio of analytic to tracked volume.
* DSNOR   ratio of analytic to tracked surface area.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  Based on the XELTI2 and XELTI3 routines.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IFTEMP
      INTEGER          IPRINT,IGTRK,NDIM,MAXMSH,NFSUR,NFREG,
     >                 NUCELL(3),NQUAD,NBANGL,NBUCEL,
     >                 NPLANE,NPOINT,LINMAX,MXGSUR,MXGREG,RENO,
     >                 MAXPIN,NTLINE,NBTDIR,NBDR,ITYPBC,IFMT,MAXSGL
      DOUBLE PRECISION RCUTOF,SPACLN,WEIGHT
      DOUBLE PRECISION RADIUS,CENTER(NDIM)
      INTEGER          IUNFLD(2,NBUCEL)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG),DGMESH(-1:MAXMSH,5),
     >                 DANGLT(NDIM,NQUAD,NBANGL),DDENWT(NQUAD,NBANGL),
     >                 DVNOR(NFREG,NBDR),
     >                 DSNOR(NFSUR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTLS')
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
      INTEGER          NPO2,NCUTOF,
     >                 ITDIR,IANGL,IQUAD,IPLANE,IDIR,
     >                 NPTA2,NPTA3,IPTA2,IPTA3,ITST,NBCOR(2),NBSINT,
     >                 ITRN,ICEL,ICI,JLINE,IBLIN,IELIN,NBSEG,ISURF,
     >                 IX,IY,IZ,IOX,IOY,IOZ,IOC,LMAXT
      DOUBLE PRECISION RAD2G,RAD2T,ANGLES(3,3),DNPDIR(3,2,3),
     >                 TRKORI(3),TRKOR2(3),TCUTOF(3,2,2),TRKLIM(2),
     >                 CELLPO(3,2),DSTART,DCERR,FACVOL,FACSUR,
     >                 DSVERR,DMVERR,DAVERR,DSSERR,DMSERR,DASERR,DERR
      DOUBLE PRECISION DWGT,DAWGT,VCONT,VCONTA,WGTFAC,TORIG(3)
      INTEGER          ISEG,JSEG,IREG,ILREG,ISUR,NBVERR,NBSERR,
     >                 NBV0,NBV1,NBS0,IPRINL
      INTEGER          ISD,NSDEB,ISF,NSFIN,NBREG,ISBL,NTSEG,IST
      INTEGER          IOFF,NGLINE
      DOUBLE PRECISION TOTVE,TOTVA,TOTVD,TOTVDR,DELV,
     >                 TOTSE,TOTSA,TOTSD,TOTSDR,DELS
      INTEGER          NICSS,ICISS,ICSR,ICSRR
*----
*  Allocatable arrays (local)
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMERO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LENGTH,DLENGT,
     >                                               DCINT,DWGTRK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DDIRET
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DORITR
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
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      LMAXT=4*(NBUCEL+4)
      WGTFAC=DONE
      PI=XDRCST('Pi',' ')
      MAXSGL=0
      ITST=RENO
      ITST=1
      NPTA2=NPOINT
      NPTA3=1
      RAD2G=RADIUS*RADIUS
      IF(NDIM .EQ. 3) NPTA3=NPOINT
      CALL XDDSET(DVNOR,NFREG*NBDR,DZERO)
      CALL XDDSET(DSNOR,NFSUR,DZERO)
      NPO2=NPOINT/2
      IF(RCUTOF .EQ. DZERO) THEN
        NCUTOF= 1
      ELSE
        NCUTOF= 4
      ENDIF
      IOC=0
*----
* Allocate: temporary storage (local)
*   NUMERO  region/surface identification number for segment.
*   LENGTH  segment length.
*   DLENGT  spatial location of each line segment.
*   ICINT   identification of spatial position for each
*           line segment in cell description of geometry.
*   DCINT   position of each intersection point for each
*           line segment in cell description of geometry.
*   DDIRET  direction of tracking lines.
*   DWGTRK  weight of tracking lines.
*   DORITR  origin of tracking lines.
*----
      ALLOCATE(NUMERO(LINMAX),ICINT(0:5,LMAXT))
      ALLOCATE(LENGTH(LINMAX),DLENGT(LINMAX),DCINT(LMAXT),
     >         DDIRET(NDIM,NQUAD*NBANGL),DWGTRK(NQUAD*NBANGL),
     >         DORITR(NDIM*(NDIM+1),NPLANE,NQUAD*NBANGL))
*----
*  Compute number of track directions
*----
      NBTDIR=0
      DO IANGL=1,NBANGL
        DO IQUAD=1,NQUAD
          IF(DDENWT(IQUAD,IANGL) .GT. DZERO) THEN
            NBTDIR=NBTDIR+1
          ENDIF
        ENDDO
      ENDDO
      FACVOL=DTWO
      IF(NDIM .EQ. 2) THEN
        FACSUR=PI
      ELSE
        FACSUR=DTWO*DTWO
      ENDIF
      IF(IPRINT .GT. 1000) THEN
        IF(ITYPBC .EQ. 2) THEN
            WRITE(IOUT,6070) 'DGMESH X',NUCELL(1)
            WRITE(IOUT,6071) (DGMESH(IX,1),IX=-1,NUCELL(1))
            WRITE(IOUT,6070) 'DGMESH Y',NUCELL(1)
            WRITE(IOUT,6071) (DGMESH(IX,2),IX=-1,NUCELL(1))
            IF(NDIM .EQ.3) THEN
              WRITE(IOUT,6070) 'DGMESH Z',NUCELL(3)
              WRITE(IOUT,6071) (DGMESH(IX,3),IX=-1,NUCELL(3))
            ENDIF
        ELSE
          DO IDIR=1,NDIM
            WRITE(IOUT,6070) 'DGMESH '//CDIR(IDIR),NUCELL(IDIR)
            WRITE(IOUT,6071) (DGMESH(IX,IDIR),IX=-1,NUCELL(IDIR))
          ENDDO
        ENDIF
      ENDIF
      CALL XDDSET(DDIRET,NDIM*NQUAD*NBANGL,DZERO)
      CALL XDDSET(DWGTRK,NQUAD*NBANGL,DZERO)
      CALL XDDSET(DORITR,NDIM*(NDIM+1)*NPLANE*NQUAD*NBANGL,DZERO)
*----
*  Loop over angles in a quarter (2-D) or a quadrant (3-D)
*----
      NGLINE=0
      NTLINE=0
      ITDIR=0
      DO IANGL=1,NBANGL
*----
*  Loop over 2 quarters (2-D) or 4 quadrants (3-D)
*----
        DO IQUAD=1,NQUAD
*----
*  Do not track angle with 0 density
*  because of the problem symmetry
*----
          IF(DDENWT(IQUAD,IANGL) .EQ. DZERO) GO TO 105
          DWGT=WEIGHT/DDENWT(IQUAD,IANGL)
          DAWGT=WEIGHT
*----
*  Track this angle
*----
          ITDIR=ITDIR+1
          DWGTRK(ITDIR)=DWGT
*----
*  Find planes mormal to selected direction
*----
          CALL XDDSET(TRKORI,3,DZERO)
          CALL XDDSET(TRKOR2,3,DZERO)
          CALL XDDSET(ANGLES,3*3,DZERO)
          CALL XDDSET(DNPDIR,3*2*3,DZERO)
          DO IDIR=1,NDIM
            ANGLES(IDIR,1)=DANGLT(IDIR,IQUAD,IANGL)
            DDIRET(IDIR,ITDIR)=ANGLES(IDIR,1)
          ENDDO
          CALL NXTQPS(NDIM,ANGLES,DNPDIR)
*----
*  Loop over planes normal to direction
*  1 in 2-D and 3 in 3-D
*----
          DO IPLANE=1,NPLANE
            DO IDIR=1,NDIM
              DO JLINE=2,NDIM
                ANGLES(IDIR,JLINE)=SPACLN*DNPDIR(IDIR,JLINE-1,IPLANE)
              ENDDO
              IF(NCUTOF .NE. 1)THEN
                TCUTOF(IDIR,1,1)=RCUTOF*(ANGLES(IDIR,2)+ANGLES(IDIR,3))
                TCUTOF(IDIR,1,2)=RCUTOF*(ANGLES(IDIR,2)-ANGLES(IDIR,3))
                TCUTOF(IDIR,2,1)=-TCUTOF(IDIR,1,2)
                TCUTOF(IDIR,2,2)=-TCUTOF(IDIR,1,1)
              ENDIF
              TRKOR2(IDIR)=CENTER(IDIR)
     >                    -DBLE(NPO2+1)*(ANGLES(IDIR,2)+ANGLES(IDIR,3))
            ENDDO
*----
*  Fill array for localisation of integration lines
*----
            DO IDIR=1,NDIM
              DORITR(IDIR,IPLANE,ITDIR)=TRKOR2(IDIR)
            ENDDO
            IOFF=NDIM
            DO JLINE=1,NDIM
              DO IDIR=1,NDIM
                DORITR(IDIR+IOFF,IPLANE,ITDIR)=ANGLES(IDIR,JLINE)
              ENDDO
              IOFF=IOFF+NDIM
            ENDDO
*----
*  Loop over lines on first normal axis (direction 2)
*----
            DO IPTA2=1,NPTA2
*----
*  Displace starting point by an additional value of SPALLN
*  along first normal axis (direction 2)
*----
              DO IDIR=1,NDIM
                TRKOR2(IDIR)=TRKOR2(IDIR)+ANGLES(IDIR,2)
                TRKORI(IDIR)=TRKOR2(IDIR)
              ENDDO
*----
*  For 3-D models:
*  Loop over lines on second normal axis (direction 3)
*  For 2-D models:
*  NPTA3=1
*----
              DO IPTA3=1,NPTA3
                RAD2T=0.0
*----
*  For 3-D models:
*  Displace starting point by an additional value of SPALLN
*  along second normal axis (direction 3)
*  For 2-D models:
*  No displacement since ANGLES(IDIR,3)=0
*----
                NGLINE=NGLINE+1
                DO IDIR=1,NDIM
                  TRKORI(IDIR)=TRKORI(IDIR)+ANGLES(IDIR,3)
                  RAD2T=RAD2T+(TRKORI(IDIR)-CENTER(IDIR))**2
                ENDDO
*----
*  Eliminate tracks outside circle or sphere surrounding geometry
*----
                IF(RAD2T .GT. RAD2G) GO TO 115
*----
*  Find cells crossed by track
*----
                IPRINL=IPRINT
*                IF(NTLINE+1 .EQ. 676
*     >        .OR. NTLINE+1 .EQ. 1814
*     >        .OR. NTLINE+1 .EQ. 1831
*     >                           ) IPRINL=IPRINT+4000
                IRLA=-1
                IF(ITYPBC .EQ. 0) THEN
*----
*  Cartesian assembly
*----
                  IRLA=NXTLCA(IPRINL,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >                        NUCELL,TRKORI,ANGLES,DGMESH,
     >                        NBCOR ,NBSINT,ICINT ,DCINT)
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
     >                        NUCELL,TRKORI,ANGLES,DGMESH,
     >                        NBCOR ,NBSINT,ICINT ,DCINT)
                ENDIF
                IF(IRLA .EQ. -1) CALL XABORT(NAMSBR//
     >': This type of cell cannot be tracked by NXT:')
*----
*  When no external face crossed go to next line
*----
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
*                write(6,*) 'NICSS,ICSRR=',NICSS,ICSRR
                DO ICISS=1,NICSS
                ICSR=ICSRR
                NTLINE=NTLINE+1
                IF(IPRINL .GE. 500) THEN
                  WRITE(IOUT,6039)
                  WRITE(IOUT,6040) NTLINE,ITDIR,IPLANE,
     >            IPTA2,IPTA3,WEIGHT,SPACLN
                  WRITE(IOUT,6041) NTLINE,TRKORI
                  WRITE(IOUT,6042) NTLINE,ANGLES
                  WRITE(IOUT,6043) NTLINE,DCINT(ICSR-1)
                ELSE IF(IPRINL .GE. 50) THEN
                  WRITE(IOUT,6040) NTLINE,ITDIR,IPLANE,
     >            IPTA2,IPTA3,WEIGHT,SPACLN
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
                    IF(NDIM .EQ. 3) THEN
                      IZ=ICINT(3,ICI)
                      IOZ=(IZ-1)*NUCELL(1)
                      CELLPO(3,2)=DGMESH(IZ,3)
                      CELLPO(3,1)=DGMESH(IZ-1,3)
                      IOC=IOZ+IX                  
                    ENDIF
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
     >                        NDIM  ,MAXMSH,LINMAX,MXGSUR,MXGREG,
     >                        MAXPIN,CELLPO,TRKLIM,TRKORI,ANGLES,
     >                        IBLIN ,IELIN ,NUMERO,DLENGT)
                  DERR=MAX(ABS(TRKLIM(1)),ABS(TRKLIM(2)))
                  DERR=(DLENGT(IELIN)-DSTART)/DERR
                  IF(DERR .GT. DCUTOF) THEN
                    WRITE(IOUT,9100) NAMSBR,NTLINE,
     >              DLENGT(IELIN),DSTART,DERR
                    CALL XABORT(NAMSBR//
     >': End of track does not coincide with end of cell')
                  ENDIF
                  IBLIN=IELIN+1
                ENDDO
 125            CONTINUE               
                NBSEG=IELIN
*----
*  Compress tracking vector for region with DLENGT=0.0 and
*  for successive segments in the same region.
*----
                JSEG=0
                ILREG=-1
                NBREG=0
                DO ISEG=1,NBSEG
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
                NBSEG=JSEG
*----
*  Add contribution of track to volume integration for this angle
*  in this quadrant and compress tracking line by removing
*  segments with NUMERO=0.
*----
                JSEG=0
                NSDEB=0
                NSFIN=0
                NBREG=0
                DO ISEG=1,NBSEG
                  IREG=NUMERO(ISEG)
                  IF(DLENGT(ISEG) .GT. DZERO) THEN
                    IF(IREG .GT. 0) THEN
                      NBREG=NBREG+1
                      VCONT=DLENGT(ISEG)*DWGT*FACVOL
                      DVNOR(IREG,1)=DVNOR(IREG,1)+VCONT
                      IF(NBDR .GT. 1) THEN
                        VCONTA=DLENGT(ISEG)*DAWGT
                        DVNOR(IREG,ITDIR+1)=DVNOR(IREG,ITDIR+1)+VCONTA
                      ENDIF
                      JSEG=JSEG+1
                      LENGTH(JSEG)=DLENGT(ISEG)
                      NUMERO(JSEG)=IREG
                    ELSE IF(IREG .LT. 0) THEN
                      IF(NBREG .EQ. 0) THEN
                        NSDEB=NSDEB+1
                      ELSE
                        NSFIN=NSFIN+1
                      ENDIF
                      DSNOR(-IREG)=DSNOR(-IREG)
     >                +DLENGT(ISEG)*DWGT*FACSUR
                      JSEG=JSEG+1
                      LENGTH(JSEG)=DLENGT(ISEG)
                      NUMERO(JSEG)=IREG
                    ENDIF
                  ENDIF
                ENDDO
                NBSEG=JSEG
                IF(NSDEB .GT. 1 .OR. NSFIN .GT. 1) THEN
                  IF(IPRINL .GE. 500) THEN
                    WRITE(IOUT,6025) NTLINE,ITDIR,NBSEG,
     >                               NSDEB,NSFIN,NBREG,
     >                               DLENGT(1),DLENGT(NBSEG)
                    WRITE(IOUT,6023) 'LineReg',NTLINE
                    WRITE(IOUT,6021) (NUMERO(ISEG),ISEG=1,NBSEG)
                    WRITE(IOUT,6024)
                    WRITE(IOUT,6023) 'LinePos',NTLINE
                    WRITE(IOUT,6022) (LENGTH(ISEG),ISEG=1,NBSEG)
                    WRITE(IOUT,6024)
                  ENDIF
                ENDIF
                MAXSGL=MAX(MAXSGL,NBSEG)
                IF(NSDEB*NSFIN .EQ. 0) THEN
*----
*  Missing outer or inner surface
*  Skip track and return warning
*----
                  WRITE(IOUT,9006) NAMSBR
                  WRITE(IOUT,9026) NTLINE,ITDIR,NBSEG,NSDEB,NSFIN,
     >                             DLENGT(1),DLENGT(NBSEG)
                ELSE
*----
*  Store line on temporary tracking file if required
*----
                  WGTFAC=DONE/DBLE(NSDEB*NSFIN)
                  IF(IPRINL .GE. 500) THEN
                    WRITE(IOUT,6020) NTLINE,ITDIR,NBSEG,DWGT*WGTFAC,
     >                               DLENGT(1),DLENGT(NBSEG)
*                    WRITE(6,*) IPLANE,IPTA2,IPTA3,
*     >                    (TRKORI(IST),IST=1,NDIM),DCINT(ICSR-1)
                    WRITE(IOUT,6023) 'LineReg',NTLINE
                    WRITE(IOUT,6021) (NUMERO(ISEG),ISEG=1,NBSEG)
                    WRITE(IOUT,6024)
                    WRITE(IOUT,6023) 'LinePos',NTLINE
                    WRITE(IOUT,6022) (LENGTH(ISEG),ISEG=1,NBSEG)
                    WRITE(IOUT,6024)
                  ENDIF
                  NTSEG=NBSEG-NSDEB-NSFIN+2
                  ISBL=0
                  DO ISD=1,NSDEB
                    DO ISF=0,NSFIN-1
                      ISBL=ISBL+1
                      IF(NSDEB*NSFIN .GT. 1 .AND. IPRINL .GE. 500)
     >                WRITE(IOUT,6026) NTLINE-1+ISBL,NTLINE,ISBL
                      IF(IGTRK .EQ. 1) THEN
                        IF(IFMT .EQ. 1) THEN
                          DO IST=1,NDIM
                            TORIG(IST)=TRKORI(IST)+DCINT(ICSR-1)*
     >                      DANGLT(IST,IQUAD,IANGL)
                          ENDDO
                          WRITE(IFTEMP) 1,NTSEG,DWGT*WGTFAC,ITDIR,
     >                    NUMERO(ISD),
     >                    (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                    NUMERO(NBSEG-ISF),
     >                    DONE,
     >                    (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                    DONE,
     >                    NTLINE-1+ISBL,IPLANE,IPTA2,IPTA3,
     >                    (TORIG(IST),IST=1,NDIM)
                        ELSE
                          WRITE(IFTEMP) 1,NTSEG,DWGT*WGTFAC,ITDIR,
     >                    NUMERO(ISD),
     >                    (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                    NUMERO(NBSEG-ISF),
     >                    DONE,
     >                    (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                    DONE
                        ENDIF
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF
                IF(ISBL .EQ. 0) THEN
                  WRITE(IOUT,6027) NTLINE,ITDIR,NBSEG,
     >                             NSDEB,NSFIN,NBREG,
     >                             NUMERO(1),DLENGT(1),
     >                             NUMERO(NBSEG),DLENGT(NBSEG)
                  ISBL=1
                  IF(IGTRK .EQ. 1) THEN
                    WRITE(IFTEMP) 1,NTSEG,DWGT*WGTFAC,ITDIR,
     >              NUMERO(1),
     >              (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >              NUMERO(NBSEG),
     >              DONE,
     >              (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >              DONE
                  ENDIF
                ENDIF
                NTLINE=NTLINE-1+ISBL
                ENDDO
*----
*  Exit because line is outside circle or sphere surrounding geometry
*----
 115            CONTINUE
*----
*  END loop over points on second normal axis
*----
              ENDDO
*----
*  END loop over points on first normal axis
*----
            ENDDO
*----
*  END loop over planes
*----
          ENDDO
*----
*  Exit because angle with 0 density not tracked
*----
 105      CONTINUE
*----
*  END loop over quarter or quadrant
*----
        ENDDO
*----
*  END loop over angles
*----
      ENDDO
*----
*  Save general tracking information
*----
      CALL LCMPUT(IPTRK,'TrackingOrig',
     >            NDIM*(NDIM+1)*NPLANE*NBTDIR,4,DORITR)
      CALL LCMPUT(IPTRK,'TrackingWgtD',NBTDIR,4,DWGTRK)
      CALL LCMPUT(IPTRK,'TrackingDirc',NDIM*NBTDIR,4,DDIRET)
*----
*  Compute DVNOR and DSNOR by comparing ratio of analytical
*  to numerically integrated volume or surfaces.
*----
      NBVERR=0
      DSVERR=DZERO
      DMVERR=DZERO
      DAVERR=DZERO
      NBV0=0
      NBV1=0
      TOTVE=DZERO
      TOTVA=DZERO
      TOTVD=DZERO
      TOTVDR=DZERO
      DO IREG=1,NFREG
        IF(IPRINT .GE. 20 ) THEN
          WRITE(IOUT,6030) IREG,
     >    SURVOL(IREG),DVNOR(IREG,1)
        ENDIF
        TOTVE=TOTVE+SURVOL(IREG)
        TOTVA=TOTVA+DVNOR(IREG,1)
        DELV=SURVOL(IREG)-DVNOR(IREG,1)
        TOTVDR=TOTVDR+DELV*DELV
        DO IDIR=1,NBDR
          IF(DVNOR(IREG,IDIR) .EQ. DZERO) THEN
            IF(IPRINT .GE. 10) THEN
              WRITE(IOUT,9000) NAMSBR,IREG,ITDIR
            ENDIF
            DVNOR(IREG,IDIR)=DONE
            IF(IDIR .EQ. 1) THEN
              NBV0=NBV0+1
            ELSE
              NBV1=NBV1+1
            ENDIF
          ELSE
            DVNOR(IREG,IDIR)=SURVOL(IREG)
     >                       /DVNOR(IREG,IDIR)
            IF(IDIR .EQ. 1) THEN
              NBVERR=NBVERR+1
            ENDIF
          ENDIF
        ENDDO
        DCERR=100.0D0*(DONE-DVNOR(IREG,1))
        DMVERR=MAX(DMVERR,ABS(DCERR))
        DSVERR=DSVERR+DCERR*DCERR
        DAVERR=DAVERR+DCERR
      ENDDO
      TOTVD=100.0D0*(TOTVE-TOTVA)/TOTVE
      TOTVDR=100.0D0*SQRT(TOTVDR/DBLE(NBVERR))/TOTVE
      IF(NBV0 .GT. 0) THEN
        WRITE(IOUT,9002) NAMSBR
      ENDIF
      IF(NBV1 .GT. 0) THEN
        WRITE(IOUT,9005) NAMSBR
      ENDIF
      DSVERR=SQRT(DSVERR/DBLE(NBVERR))
      DAVERR=DAVERR/DBLE(NBVERR)
      NBSERR=0
      DSSERR=DZERO
      DMSERR=DZERO
      DASERR=DZERO
      NBS0=0
      TOTSE=DZERO
      TOTSA=DZERO
      TOTSD=DZERO
      TOTSDR=DZERO
      DO ISUR=1,NFSUR
        IF(IPRINT .GE. 20 ) THEN
          WRITE(IOUT,6031) ISUR,
     >    SURVOL(-ISUR),DSNOR(ISUR)
        ENDIF
        TOTSE=TOTSE+SURVOL(-ISUR)
        TOTSA=TOTSA+DSNOR(ISUR)
        DELS=SURVOL(-ISUR)-DSNOR(ISUR)
        TOTSDR=TOTSDR+DELS*DELS
        IF(DSNOR(ISUR) .EQ. DZERO) THEN
          IF(IPRINT .GE. 10) THEN
            WRITE(IOUT,9001) NAMSBR,-ISUR
          ENDIF
          NBS0=NBS0+1
          DSNOR(ISUR)=DONE
        ELSE
          DSNOR(ISUR)=SURVOL(-ISUR)
     >                /DSNOR(ISUR)
          NBSERR=NBSERR+1
        ENDIF
        DCERR=100.0D0*(DONE-DSNOR(ISUR))
        DMSERR=MAX(DMSERR,ABS(DCERR))
        DSSERR=DSSERR+DCERR*DCERR
        DASERR=DASERR+DCERR
      ENDDO
      TOTSD=100.0D0*(TOTSE-TOTSA)/TOTSE
      TOTSDR=100.0D0*SQRT(TOTSDR/DBLE(NBSERR))/TOTSE
      IF(NBS0 .GT. 0) THEN
        WRITE(IOUT,9003) NAMSBR,NBS0
      ENDIF
      DSSERR=SQRT(DSSERR/DBLE(NBSERR))
      DASERR=DASERR/DBLE(NBSERR)
*----
*  Processing finished:
*  print track normalization vector.
*  and routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6007) TOTVE,TOTVA,TOTVD,TOTVDR
        WRITE(IOUT,6005) DSVERR,DMVERR,DAVERR
        IF(IPRINT .GE. 10) THEN
          DO IREG=1,NFREG
            WRITE(IOUT,6010) IREG,SURVOL(IREG)
            WRITE(IOUT,6012) DVNOR(IREG,1),
     >      100.0D0*(DONE-DVNOR(IREG,1))
          ENDDO
        ENDIF
        WRITE(IOUT,6008) TOTSE,TOTSA,TOTSD,TOTSDR
        WRITE(IOUT,6006) DSSERR,DMSERR,DASERR
        IF(IPRINT .GE. 10) THEN
          DO ISUR=1,NFSUR
            WRITE(IOUT,6011) -ISUR,SURVOL(-ISUR)
            WRITE(IOUT,6012) DSNOR(ISUR),
     >      100.0D0*(DONE-DSNOR(ISUR))
          ENDDO
          WRITE(IOUT,6001) NAMSBR
        ENDIF
      ENDIF
*----
*  Save track normalisation vector
*----
      IF(MAXSGL .EQ. 0) THEN
        WRITE(IOUT,9004) NAMSBR
        MAXSGL=LINMAX
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DORITR,DWGTRK,DDIRET,DCINT,DLENGT,LENGTH)
      DEALLOCATE(ICINT,NUMERO)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6005 FORMAT(' Global RMS, maximum and average errors (%) ',
     >       'on region volumes :',3(2X,F10.5))
 6006 FORMAT(' Global RMS, maximum and average errors (%) ',
     >       'on surface areas  :',3(2X,F10.5))
 6007 FORMAT(' Total exact volume       = ',F19.10/
     >       ' Total approximate volume = ',F19.10/
     >       ' Error on total volume    = ',4X,F10.5/
     >       ' RMS Error on total volume= ',4X,F10.5)
 6008 FORMAT(' Total exact surface area       = ',F19.10/
     >       ' Total approximate surface area = ',F19.10/
     >       ' Error on total surface area    = ',4X,F10.5/
     >       ' RMS Error on total surface area= ',4X,F10.5)
 6010 FORMAT(' Normalisation factors and relative errors (%) ',
     >       'for region ',I8,' with volume ',F19.10)
 6011 FORMAT(' Normalisation factors and relative error (%) ',
     >       'for surface ',I8,' with area   ',F19.10)
 6012 FORMAT((2X,F15.10,2X,F10.5))
 6020 FORMAT('Line',I10.10,'={',2(I10,','),
     >       F18.10,',',F18.10,',',F18.10,'};')
 6021 FORMAT(6(I10,:,',',9X))
 6022 FORMAT(6(F18.10,:,','))
 6023 FORMAT(A7,I10.10,'={')
 6024 FORMAT(18X,'};')
 6025 FORMAT('Line',I10.10,'={',5(I10,','),
     >       F18.10,',',F18.10,'};')
 6026 FORMAT('Line',I10.10,'={',I10,',',I10,'};')
 6027 FORMAT('Problem with Line',I10.10,'={',5(I10,','),
     >       I10,',',F18.10,',',I10,',',F18.10,'};')
 6030 FORMAT(' Normalization volumes  =',I10,1P,2D20.10)
 6031 FORMAT(' Normalization surfaces =',I10,1P,2D20.10)
 6039 FORMAT(1X)
 6040 FORMAT('Track',I10.10,'={',4(I10,','),
     >F15.8,',',F15.8,'};')
 6041 FORMAT('Torig',I10.10,'={',2(F15.8,','),F15.8,'};')
 6042 FORMAT('Tdire',I10.10,'={',8(F15.8,','),F15.8,'};')
 6043 FORMAT('Tstrt',I10.10,'={',F15.8,'};')
 6070 FORMAT(1X,A8,5X,I10)
 6071 FORMAT(5F18.10)
 9000 FORMAT(1X,'***** Warning in ',A6,'*****'/
     >       7X,'For region ',I8,
     >       1X,'no crossing by angle ',I8)
 9001 FORMAT(1X,'***** Warning in ',A6,'*****'/
     >       7X,'For surface ',I8,
     >       1X,'no crossing by any angle ')
 9002 FORMAT(1X,'***** Warning in ',A6,'*****'/
     >       7X,' regions not tracked for any direction  ')
 9003 FORMAT(1X,' ***** Warning in ',A6,'*****'/
     >       7X,I8,' surfaces not tracked for direction ',I8)
 9004 FORMAT(1X,' ***** Warning in ',A6,'*****'/
     >       7X,' no line segments detected in tracking ')
 9005 FORMAT(1X,'***** Warning in ',A6,'*****'/
     >       7X,' regions not tracked for some directions')
 9006 FORMAT(1X,'***** Warning in ',A6,'*****'/
     >       7X,' Final or initial surface could not be identified')
 9026 FORMAT('Line',I10.10,'={',4(I10,','),
     >       F18.10,',',F18.10,'};')
 9100 FORMAT(1X,' ***** Error in ',A6,'***** for line ',I8/
     >       7X,'Positions (current and reference ) =',1P,2D21.14/
     >       7X,'Relative error = ',D21.14)
      END
