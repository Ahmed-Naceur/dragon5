*DECK NXTTLC
      SUBROUTINE NXTTLC(IPTRK ,IFTEMP,IPRINT,IGTRK ,NDIM  ,MAXMSH,
     >                  NFSUR ,NFREG ,NUCELL,NBUCEL,NBANGL,
     >                  LINMAX,MXGSUR,MXGREG,RENO  ,NBDR  ,ITYPBC,
     >                  IFMT  ,MAXPIN,AZMQUA,IPER  ,IUNFLD,SURVOL,
     >                  DGMESH,DANGLT,DDENWT,DNSANG,NBSANG,DEPART,
     >                  MAXSUB,MAXSGL,NTLINE,DVNOR ,DSNOR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To generate the cyclic tracking lines (specular tracking)
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
* IFTEMP  pointer to a temporary TRACKING file in
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
* LINMAX  maximum number of segments in a track.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
* RENO    track normalisation option. A value RENO=-1 implies
*         a direction dependent normalization of the tracks
*         for the volume while a value RENO=0, implies
*         a global normalisation.
* NBDR    number of directions for track normalization.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesianb oundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* IFMT    tracking file format:
*         =0 for short file;
*         =1 long file required by TLM:.
* MAXPIN  maximum number of pins in a cell.
* AZMQUA  tracking type.
* IPER    cell periodicity factor in each direction.
* IUNFLD  description of unfolded geometry.
* SURVOL  global surface volume vector.
* DGMESH  meshing vector for global geometry.
* DANGLT  angles.
* DDENWT  angular density for each angle.
* DNSANG  spatial density required.
* NBSANG  number of segments for each angles.
* DEPART  track starting point.
*
*Parameters: output
* MAXSUB  maximum number of subtracks in a line.
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
     >                 NUCELL(3),NBUCEL,NBANGL,LINMAX,MXGSUR,MXGREG,
     >                 RENO,NBDR,ITYPBC,IFMT,MAXPIN,AZMQUA,IPER(3),
     >                 MAXSUB,MAXSGL,NTLINE
      INTEGER          IUNFLD(2,NBUCEL)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG),DGMESH(-1:MAXMSH,4),
     >                 DANGLT(NDIM,NBANGL),DDENWT(NBANGL),
     >                 DNSANG(NBANGL)
      INTEGER          NBSANG(5,NBANGL)
      DOUBLE PRECISION DEPART(NDIM,2,NBANGL)
      DOUBLE PRECISION DVNOR(NFREG,NBDR),
     >                 DSNOR(NFSUR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTLC')
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
      INTEGER          LMAXT,ITST,IANGL,NGLINE,IDIR,ITXY(3),
     >                 NSCAN,NUMANG,NPOINT,NBPTS,ISCAN,
     >                 ISTART,ISUM,NBCOR(2),IFREF
      INTEGER          NBSINT,ISDIR,IBLIN,IELIN,ICI,IX,IY,IZ,
     >                 IOX,IOY,IOZ,IOC,ICEL,ITRN,ISURF
      DOUBLE PRECISION FACVOL,FACSUR,DWGT,WGTFAC,DAWGT,PMAX,PMIN,
     >                 DSTART,RADT2,ANGT2
      DOUBLE PRECISION ANGLED(3),ANGLES(3),ANGLEN(3),ANGLER(3),
     >                 TRKORI(3),TRKOR2(3),TRKORD(3),TRKORN(3),
     >                 TRKLIM(2),CELLPO(3,2)
      DOUBLE PRECISION DERR,DCERR,DSVERR,DMVERR,DAVERR,DSSERR,
     >                 DMSERR,DASERR,VCONTA,VCONT
      INTEGER          ISEG,JSEG,IREG,ILREG,ISUR,NBVERR,NBSERR,
     >                 NBV0,NBV1,NBS0,IPRINL,IPRINC
      INTEGER          ISD,NSDEB,ISF,NSFIN,NBREG,ISBL,NTSEG
      INTEGER          ISINT,NBSEG,IPT,JSDIR,ISDIRX,IPTN
      DOUBLE PRECISION DHALF
      DOUBLE PRECISION AAA
      INTEGER          IND,II
      LOGICAL          LNEW
*----
*  Allocatable arrays (local)
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NUMERO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ICINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LENGTH,DLENGT,
     >                                               DCINT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KANGL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: TORIG
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
* Allocate: temporary storage (local)
* NUMERO  region/surface identification number
*         for segment.
* LENGTH  segment length.
* DLENGT  spatial location of each line segment.
* ICINT   identification of spatial position for each
*         line segment in cell description of geometry.
* DCINT   position of each intersection point for each
*         line segment in cell description of geometry.
* TORIG   track origin
* KANGL   angle index (quadrant)
*----
      ALLOCATE(NUMERO(LINMAX),LENGTH(LINMAX),DLENGT(LINMAX),
     > ICINT(0:5,2*(NBUCEL+4)),DCINT(2*(NBUCEL+4)))
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 2) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      PI=XDRCST('Pi',' ')
      DHALF=DONE/DTWO
      MAXSUB=0
      MAXSGL=0
      NTLINE=0
      LMAXT=2*(NBUCEL+4)
      ITST=RENO
      AAA=AZMQUA
      ITST=1
      CALL XDDSET(DVNOR,NFREG*NBDR,DZERO)
      CALL XDDSET(DSNOR,NFSUR,DZERO)
      FACVOL=DTWO
      IF(NDIM .EQ. 2) THEN
        FACSUR=DTWO*PI
      ELSE
        FACSUR=DTWO*DTWO
      ENDIF
      IF(IPRINT .GT. 1000) THEN
        DO IDIR=1,NDIM
          WRITE(IOUT,6070) 'DGMESH '//CDIR(IDIR),NUCELL(IDIR)
          WRITE(IOUT,6071) (DGMESH(IX,IDIR),IX=1,NUCELL(IDIR))
        ENDDO
      ENDIF
*----
*  Loop over directions
*----
      IPRINL=IPRINT
      IPRINC=IPRINT
      NGLINE=0
      DO IANGL=1,NBANGL
        DO IDIR=1,NDIM
          ITXY(IDIR)=NBSANG(IDIR,IANGL)
        ENDDO
        ITXY(3)=ITXY(1)*ITXY(2)
        NSCAN=NBSANG(3,IANGL)
        NUMANG=NBSANG(4,IANGL)
        NPOINT=NBSANG(5,IANGL)
        ALLOCATE(TORIG(NDIM,NUMANG),KANGL(NUMANG))
*----
*  Loop over quadrant if both quadrant not scanned simultaneously
*  If NSCAN = 2  0 TO PI/2 AND PI/2 TO PI are scanned simultaneously
*  If NSCAN = 1  0 TO PI/2 AND PI/2 TO PI are scanned independently
*----
        NBPTS=NSCAN*NPOINT
        DO ISCAN=2,NSCAN,-1
*----
*  Initialize weight
*  Store track direction in ANGLES
*  and track starting point in TRKOR2
*----
          DWGT=DNSANG(IANGL)
          DAWGT=DWGT*DDENWT(IANGL)
          DO IDIR=1,NDIM
            ANGLED(IDIR)=DANGLT(IDIR,IANGL)
            TRKOR2(IDIR)=DEPART(IDIR,1,IANGL)
          ENDDO
          IF(ISCAN .EQ. 1) THEN
            IF(ITXY(1) .EQ. 0) THEN
              ANGLED(2)=-ANGLED(2)
            ELSE IF(ITXY(2) .EQ. 0) THEN
              ANGLED(1)=-ANGLED(1)
            ENDIF
          ENDIF
          DO IDIR=1,NDIM
            ANGLES(IDIR)=ANGLED(IDIR)
            ANGLEN(IDIR)=ANGLED(IDIR)
            ANGLER(IDIR)=ANGLED(IDIR)
            TRKORD(IDIR)=TRKOR2(IDIR)
            TRKORI(IDIR)=TRKORD(IDIR)
            TRKORN(IDIR)=TRKORD(IDIR)
          ENDDO
*----
*  Loop over points
*----
          IF(IPRINT .GE. 10) THEN
            WRITE(IOUT,6050) IANGL,ISCAN,ITXY,DWGT,
     >      (DANGLT(IDIR,IANGL),IDIR=1,NDIM),
     >      (DEPART(IDIR,1,IANGL),IDIR=1,NDIM)
          ENDIF
          ISTART=0
          ISUM=ISTART
          IFREF=1
          DO IPT=1,NBPTS
*----
*  Find if this track crosses geometry
*----
            IF(ISTART .EQ. 0) THEN
              IF(IFREF .EQ. 1) THEN
                DO IDIR=1,NDIM
                  ANGLES(IDIR)=ANGLER(IDIR)
                  ANGLED(IDIR)=ANGLER(IDIR)
                ENDDO
              ENDIF
            ENDIF
            ISINT=1
            DCINT(ISINT)=DZERO
            IRLA=-1
            IF(ITYPBC .EQ. 0) THEN
*----
*  Cartesian assembly
*----
              IRLA=NXTLCA(IPRINC,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >                     NUCELL,TRKORI,ANGLES,DGMESH,
     >                     NBCOR ,NBSINT,ICINT ,DCINT)
            ELSE IF(ITYPBC .EQ. 1) THEN
*----
*  Annular assembly
*----
              CALL XABORT(NAMSBR//': Circular BC not implemented')
            ELSE IF(ITYPBC .EQ. 2) THEN
*----
*  Hexagonal assembly
*----
              IRLA=NXTLHA(IPRINC,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >                     NUCELL,TRKORI,ANGLES,DGMESH,
     >                     NBCOR ,NBSINT,ICINT ,DCINT)
            ENDIF
            IF(IRLA .EQ. -1) CALL XABORT(NAMSBR//
     >': This type of cell cannot be tracked by NXT:')
*----
*  When no external face crossed go to next line
*----
            IF(IRLA .EQ. 0) GO TO 105
            ISDIRX=ABS(ICINT(0,1))
            ISDIR=ABS(ICINT(0,NBSINT))
            JSDIR=MOD(ISDIR,2)+1
            IF(ISTART .EQ. 0) THEN
*----
*  Define position of cyclic starting point
*----
              NGLINE=NGLINE+1
              DO IDIR=1,NDIM
                TRKORD(IDIR)=TRKORD(IDIR)+DCINT(ISINT)*ANGLES(IDIR)
                TRKORI(IDIR)=TRKORD(IDIR)
              ENDDO
              DO IDIR=1,NBCOR(1)
                DCINT(IDIR)=DZERO
              ENDDO
              IF(ISCAN .EQ. 1 .AND. ITXY(3) .NE. 0) THEN
                ANGLED(ISDIRX)=-ANGLED(ISDIRX)
                DO IDIR=1,NDIM
                  ANGLES(IDIR)=ANGLED(IDIR)
                ENDDO
                IF(ITYPBC .EQ. 0) THEN
*----
*  Cartesian assembly
*----
                  IRLA=NXTLCA(IPRINL,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >                         NUCELL,TRKORI,ANGLES,DGMESH,
     >                         NBCOR ,NBSINT,ICINT ,DCINT)
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
     >                         NUCELL,TRKORI,ANGLES,DGMESH,
     >                         NBCOR ,NBSINT,ICINT ,DCINT)
                ELSE
                  CALL XABORT(NAMSBR//
     >': This type of cell cannot be tracked by NXT:')
                ENDIF
                IF(IRLA .EQ. 0) GO TO 105
                DO IDIR=1,NDIM
                  TRKORD(IDIR)=TRKORD(IDIR)+DCINT(ISINT)*ANGLES(IDIR)
                  TRKORI(IDIR)=TRKORD(IDIR)
                ENDDO
                DO IDIR=1,NBCOR(1)
                  DCINT(IDIR)=DZERO
                ENDDO
                ISDIR=ABS(ICINT(0,NBSINT))
                JSDIR=MOD(ISDIR,2)+1
              ENDIF
              DO IDIR=1,NDIM
                ANGLES(IDIR)=ANGLED(IDIR)
                ANGLEN(IDIR)=ANGLED(IDIR)
              ENDDO
            ENDIF
*----
*  For track crossing geometry, end find surface direction
*  and end of track location
*----
            DO IDIR=1,NDIM
              AAA=DZERO
              DO ISINT=NBCOR(1),NBSINT-NBCOR(2)+1
                AAA=AAA+DCINT(ISINT)*ANGLES(IDIR)
              ENDDO
              TRKORN(IDIR)=TRKORI(IDIR)+AAA
              TORIG(IDIR,ISTART+1)=TRKORI(IDIR)
            ENDDO
            KANGL(ISTART+1)=0
            DO II=1,4*NBANGL
              IF((ANGLEN(1).EQ.DANGLT(1,II)).AND.
     >           (ANGLEN(2).EQ.DANGLT(2,II))) THEN
                  KANGL(ISTART+1)=II
                  GO TO 10
              ENDIF
            ENDDO
            CALL XABORT(NAMSBR//
     >': Unable to find an angular index for a subtrack')
 10         CONTINUE
            IPRINL=IPRINT
*            IF(NTLINE+1 .EQ. 1411
*     >    .OR. NTLINE+1 .EQ. 1414
*     >    .OR. NTLINE+1 .EQ. 1415
*     >                            ) IPRINL=IPRINT+4000
            IF(ISTART.EQ.0) THEN
              IF(IPRINL .GE. 500) THEN
                WRITE(IOUT,6051) NTLINE+1
              ENDIF
            ENDIF
            IF(IPRINL .GE. 500) THEN
              WRITE(IOUT,6052) TRKORI,ANGLEN
            ENDIF
*----
*  Find location with respect to end surfaces positions
*----
            PMAX=ABS(TRKORN(ISDIR)-DGMESH(NUCELL(ISDIR),ISDIR))
            PMIN=ABS(TRKORN(ISDIR)-DGMESH(0,ISDIR))
            IF(IPER(ISDIR) .EQ. 1) THEN
*----
*  Translate next initial starting point if possible
*----
              IF(PMAX .LT. DCUTOF) THEN
                TRKORN(ISDIR)=DGMESH(0,ISDIR)
              ELSE IF(PMIN .LT. DCUTOF) THEN
                TRKORN(ISDIR)=DGMESH(NUCELL(ISDIR),ISDIR)
              ELSE
                WRITE(IOUT,9100) NAMSBR,ISDIR,NUCELL(ISDIR),
     >          TRKORI(ISDIR),DGMESH(NUCELL(ISDIR),ISDIR),PMAX
                CALL XABORT(NAMSBR//
     >          ': Translation before cell limit attempted')
              ENDIF
            ELSE
*----
*  Change track direction upon reflexion if possible
*----
              IF(PMAX .GT. DCUTOF .AND. PMIN .GT. DCUTOF) THEN
                CALL XABORT(NAMSBR//
     >          ': Reflexion before cell limit attempted')
              ENDIF
              ANGLEN(ISDIR)=-ANGLEN(ISDIR)
            ENDIF
*----
*  Track seems ok, locate track segments for geometry
*----
            IF(ISTART .EQ. 0) THEN
*----
*  For first series of line segments initialize track information
*----
              IBLIN=1
              IELIN=0
            ELSE
*----
*  Next series of line segments, add at the end of
*  current line
*----
              IBLIN=IELIN+1
              DLENGT(IBLIN)=DZERO
            ENDIF
            ISTART=ISTART+1
            ISUM=ISUM+1
            NUMERO(IBLIN)=0
            DLENGT(IBLIN)=DCINT(1)
            DSTART=DCINT(1)
            IBLIN=IBLIN+1
*----
*  Find line segments for this track series
*----
            DO ICI=NBCOR(1)+1,NBSINT-NBCOR(2)+1
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
     >                    NDIM  ,MAXMSH,LINMAX,MXGSUR,MXGREG,
     >                    MAXPIN,CELLPO,TRKLIM,TRKORI,ANGLES,
     >                    IBLIN ,IELIN ,NUMERO,DLENGT)
              IPRINL=IPRINT
              AAA=MAX(ABS(TRKLIM(1)),ABS(TRKLIM(2)))
              DERR=ABS(DLENGT(IELIN)-DSTART)
              IF(AAA .GT. DONE) DERR=DERR/AAA
              IF(DERR .GT. DCUTOF) THEN
                IF(AAA .GT. DONE) THEN
                  WRITE(IOUT,9101) NAMSBR,NGLINE,
     >            DLENGT(IELIN),DSTART,DERR
                ELSE
                  WRITE(IOUT,9102) NAMSBR,NGLINE,
     >            DLENGT(IELIN),DSTART,DERR
                ENDIF
                CALL XABORT(NAMSBR//
     >': End of track does not coincide with end of cell')
              ENDIF
              IBLIN=IELIN+1
*----
*  End of tracking for this line segment
*----
            ENDDO
            IPRINL=IPRINT
*----
*  This series of line segments completed
*  Test if the cycle is back to starting point
*----
            RADT2=DZERO
            ANGT2=DZERO
            DO IDIR=1,NDIM
              RADT2=RADT2+(TRKORD(IDIR)-TRKORN(IDIR))**2
              ANGT2=ANGT2+(ANGLED(IDIR)-ANGLEN(IDIR))**2
            ENDDO
            RADT2=SQRT(RADT2)
            ANGT2=SQRT(ANGT2)
            IF(RADT2 .GT. DCUTOF .OR. ANGT2 .GT. DCUTOF) THEN
*----
*  Cycle is incomplete
*  Check if period is higher than predicted
*  reset starting point and direction
*----
              IF(ISTART .GE. NUMANG) CALL XABORT(NAMSBR//
     >        ': Cyclic period is too long')
                IF(IPRINL .GE. 500) THEN
                  IF(IPT .EQ. NBPTS) WRITE(IOUT,6064)
                ENDIF
              DO IDIR=1,NDIM
                TRKORI(IDIR)=TRKORN(IDIR)
                ANGLES(IDIR)=ANGLEN(IDIR)
              ENDDO
            ELSE
*----
*  Cycle is complete
*  Check if period is ok
*----
              NTLINE=NTLINE+1
*----
*  Process tracking line and save
*----
              NBSEG=IELIN
*----
*  Compress tracking vector for region with DLENGT=0.0,
*  for multiple surface of intersection and
*  for successive segments in the same region.
*----
              ISD=-1
              JSEG=0
              ILREG=-1
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
                    ISD=IREG
                  ELSE
                    ISF=NUMERO(ISEG+1)
                    IF(ISD .LT. 0) THEN
                      IF(ISF .GT. 0) THEN
                        JSEG=JSEG+1
                        DLENGT(JSEG)=DLENGT(ISEG)*DHALF
                        NUMERO(JSEG)=IREG
                        ILREG=-1
                      ENDIF
                    ELSE
                      JSEG=JSEG+1
                      DLENGT(JSEG)=DLENGT(ISEG)*DHALF
                      NUMERO(JSEG)=IREG
                      ILREG=-1
                    ENDIF
                    ISD=IREG
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
              NBREG=0
              IND=0
              LNEW=.TRUE.
              DO ISEG=1,NBSEG
                IREG=NUMERO(ISEG)
                IF(DLENGT(ISEG) .GT. DZERO) THEN
                  IF(IREG .GT. 0) THEN
                    IF(LNEW) THEN
                      IND=IND+1
                      IF(IND.GT.ISTART) CALL XABORT(NAMSBR//
     >                ': ISTART overflow')
                      LNEW=.FALSE.
                    ENDIF
                    NBREG=NBREG+1
                    VCONT=DLENGT(ISEG)*DWGT*FACVOL
                    DVNOR(IREG,1)=DVNOR(IREG,1)+VCONT
                    IF(NBDR .GT. 1) THEN
                      II=KANGL(IND)
                      IF(II .GT. 2*NBANGL) II=II-2*NBANGL
                      IF(DANGLT(1,II).EQ.DZERO) THEN
                        VCONTA=DHALF*DLENGT(ISEG)*DAWGT
                      ELSE
                        VCONTA=DLENGT(ISEG)*DAWGT
                      ENDIF
                      DVNOR(IREG,II+1)=DVNOR(IREG,II+1)+VCONTA
                    ENDIF
                    JSEG=JSEG+1
                    LENGTH(JSEG)=DLENGT(ISEG)
                    NUMERO(JSEG)=IREG
                  ELSE IF(IREG .LT. 0) THEN
                    LNEW=.TRUE.
                    IF(NBREG .EQ. 0) THEN
                      NSDEB=NSDEB+1
                    ENDIF
                    DSNOR(-IREG)=DSNOR(-IREG)+DLENGT(ISEG)*DWGT*FACSUR
                    JSEG=JSEG+1
                    LENGTH(JSEG)=DLENGT(ISEG)
                    NUMERO(JSEG)=IREG
                  ENDIF
                ENDIF
              ENDDO
              NBSEG=JSEG
              NSFIN=0
              DO ISEG=NBSEG,1,-1
                IREG=NUMERO(ISEG)
                IF(DLENGT(ISEG) .GT. DZERO) THEN
                  IF(IREG .GT. 0) THEN
                    GO TO 115
                  ELSE IF(IREG .LT. 0) THEN
                    NSFIN=NSFIN+1
                  ENDIF
                ENDIF
              ENDDO
 115          CONTINUE
*----
*  Matlab commands to create a simili-TLM plot
*----
              IF(IPRINL .GE. 2000) THEN
                WRITE(IOUT,6064)
                WRITE(IOUT,6063) NTLINE,'LineReg'
                WRITE(IOUT,6061) (NUMERO(ISEG),ISEG=1,NBSEG)
                WRITE(IOUT,6064)
                WRITE(IOUT,6063) NTLINE,'LinePos'
                WRITE(IOUT,6062) (LENGTH(ISEG),ISEG=1,NBSEG)
                WRITE(IOUT,6064)
                WRITE(IOUT,6065)
              ENDIF
*----
*  Store line on temporary tracking file if required
*----
              WGTFAC=DONE/DBLE(NSDEB*NSFIN)
              MAXSUB=MAX(MAXSUB,ISTART)
              MAXSGL=MAX(MAXSGL,NBSEG)
              NTSEG=NBSEG-NSDEB-NSFIN+2
              ISBL=0
              DO ISD=1,NSDEB
                DO ISF=0,NSFIN-1
                  ISBL=ISBL+1
                  IF(NSDEB*NSFIN .GT. 1 .AND. IPRINL .GE. 1000)
     >            WRITE(IOUT,6026) NTLINE-1+ISBL,NTLINE,ISBL
                  IF(IGTRK .EQ. 1) THEN
                    IF(IFMT == 1) THEN
                      WRITE(IFTEMP) ISTART,NTSEG,DWGT*WGTFAC,
     >                (KANGL(II),II=1,ISTART),
     >                NUMERO(ISD),
     >                (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                NUMERO(NBSEG-ISF),
     >                LENGTH(ISD),
     >                (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                LENGTH(NBSEG-ISF),
     >                NTLINE,1,1,1,
     >                ((TORIG(IDIR,II),IDIR=1,NDIM),II=1,ISTART)
                    ELSE
                      WRITE(IFTEMP) ISTART,NTSEG,DWGT*WGTFAC,
     >                (KANGL(II),II=1,ISTART),
     >                NUMERO(ISD),
     >                (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                NUMERO(NBSEG-ISF),
     >                LENGTH(ISD),
     >                (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >                LENGTH(NBSEG-ISF)
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
              IF(ISBL .EQ. 0) THEN
                WRITE(IOUT,6027) NTLINE,IANGL,NBSEG,
     >                           NSDEB,NSFIN,NBREG,
     >                           NUMERO(1),DLENGT(1),
     >                           NUMERO(NBSEG),DLENGT(NBSEG)
                ISBL=1
                IF(IFMT == 1) THEN
                  WRITE(IFTEMP) ISTART,NTSEG,DWGT*WGTFAC,
     >            (KANGL(II),II=1,ISTART),
     >            NUMERO(1),
     >            (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >            NUMERO(NBSEG),
     >            LENGTH(1),
     >            (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >            LENGTH(NBSEG),
     >            NTLINE,1,1,1,
     >            ((TORIG(IDIR,II),IDIR=1,NDIM),II=1,ISTART)
                ELSE
                  WRITE(IFTEMP) ISTART,NTSEG,DWGT*WGTFAC,
     >            (KANGL(II),II=1,ISTART),
     >            NUMERO(1),
     >            (NUMERO(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >            NUMERO(NBSEG),
     >            LENGTH(1),
     >            (LENGTH(ISEG),ISEG=NSDEB+1,NBSEG-NSFIN),
     >            LENGTH(NBSEG)
                ENDIF
              ENDIF
              NTLINE=NTLINE-1+ISBL
              IF(ISTART .EQ. NUMANG) THEN
*----
*  Displace line to next starting point
*----
                DO IDIR=1,NDIM
                  ANGLES(IDIR)=ANGLED(IDIR)
                  ANGLEN(IDIR)=ANGLED(IDIR)
                  TRKOR2(IDIR)=TRKOR2(IDIR)+DEPART(IDIR,2,IANGL)
                  TRKORD(IDIR)=TRKOR2(IDIR)
                  TRKORI(IDIR)=TRKORD(IDIR)
                  TRKORN(IDIR)=TRKORD(IDIR)
                ENDDO
              ELSE
                IF(NSCAN .EQ. 2) THEN
                  IF(IPER(ISDIR) .EQ. 1) THEN
                    IF(IFREF .EQ. -1) THEN
                      IFREF=1
                      DO IDIR=1,NDIM
                        ANGLES(IDIR)=ANGLED(IDIR)
                        ANGLEN(IDIR)=ANGLED(IDIR)
                        TRKOR2(IDIR)=TRKOR2(IDIR)+DEPART(IDIR,2,IANGL)
                        TRKORD(IDIR)=TRKOR2(IDIR)
                        TRKORI(IDIR)=TRKORD(IDIR)
                        TRKORN(IDIR)=TRKORD(IDIR)
                      ENDDO
                    ELSE
                      ANGLED(ISDIR)=-ANGLED(ISDIR)
                      IFREF=-1
                    ENDIF
                  ELSE IF(IPER(JSDIR) .EQ. 1) THEN
                    IF(IFREF .EQ. 1) THEN
                      DO IPTN=1,NBPTS
                        DO IDIR=1,NDIM
                          TRKORI(IDIR)=TRKORN(IDIR)
                          ANGLES(IDIR)=ANGLEN(IDIR)
                        ENDDO
*----
*  Find next surface intersection track crosses geometry
*----
                        IF(ITYPBC .EQ. 0) THEN
*----
*  Cartesian assembly
*----
                          IRLA=NXTLCA(IPRINL,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >                               NUCELL,TRKORI,ANGLES,DGMESH,
     >                               NBCOR ,NBSINT,ICINT ,DCINT)
                        ELSE IF(ITYPBC .EQ. 1) THEN
                          CALL XABORT(NAMSBR//
     >                    ': Circular BC not implemented')
                        ELSE IF(ITYPBC .EQ. 2) THEN
*----
*  Hexagonal assembly
*----
                          IRLA=NXTLHA(IPRINL,ITST  ,NDIM  ,MAXMSH,LMAXT,
     >                               NUCELL,TRKORI,ANGLES,DGMESH,
     >                               NBCOR ,NBSINT,ICINT ,DCINT)
                        ELSE
                          CALL XABORT(NAMSBR//
     >': This type of cell cannot be tracked by NXT:')
                        ENDIF
*----
*  Define position of cyclic starting point
*----
                        DO IDIR=1,NDIM
                          AAA=DZERO
                          DO ISINT=NBCOR(1),NBSINT-NBCOR(2)+1
                            AAA=AAA+DCINT(ISINT)*ANGLES(IDIR)
                          ENDDO
                          TRKORN(IDIR)=TRKORI(IDIR)+AAA
                        ENDDO
                        ISDIRX=ABS(ICINT(0,NBSINT))
                        IF(ISDIRX .EQ. JSDIR) THEN
                          DO IDIR=1,NDIM
                            ANGLED(IDIR)=ANGLES(IDIR)
                          ENDDO
                          ANGLED(JSDIR)=-ANGLED(JSDIR)
                          GO TO 125
                        ENDIF
                        IF(IPER(ISDIRX) .EQ. 1) THEN
                          ANGLEN(ISDIRX)=-ANGLEN(ISDIRX)
                        ENDIF
                      ENDDO
 125                  CONTINUE
                      IFREF=-1
                    ELSE
                      DO IDIR=1,NDIM
                        ANGLES(IDIR)=ANGLED(IDIR)
                        ANGLEN(IDIR)=ANGLED(IDIR)
                        TRKOR2(IDIR)=TRKOR2(IDIR)+DEPART(IDIR,2,IANGL)
                        TRKORD(IDIR)=TRKOR2(IDIR)
                        TRKORI(IDIR)=TRKORD(IDIR)
                        TRKORN(IDIR)=TRKORD(IDIR)
                      ENDDO
                      IFREF=1
                    ENDIF
                  ENDIF
                ENDIF
                DO IDIR=1,NDIM
                  ANGLES(IDIR)=ANGLED(IDIR)
                  ANGLEN(IDIR)=ANGLED(IDIR)
                  TRKORD(IDIR)=TRKORN(IDIR)
                  TRKORI(IDIR)=TRKORD(IDIR)
                  TRKORN(IDIR)=TRKORD(IDIR)
                ENDDO
              ENDIF
              ISTART=0
            ENDIF
*----
*  End of tracking for this direction
*----
 105        CONTINUE
          ENDDO
*----
*  End of tracking for this set of quadrant
*----
        ENDDO
        DEALLOCATE(KANGL,TORIG)
*----
*  End of tracking for this set of angles
*----
      ENDDO
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
      DO IREG=1,NFREG
        IF(IPRINT .GE. 20 ) THEN
          WRITE(IOUT,6030) IREG,SURVOL(IREG),DVNOR(IREG,1)
        ENDIF
        DO IDIR=1,NBDR
          IF(DVNOR(IREG,IDIR) .EQ. DZERO) THEN
            IF(IPRINT .GE. 10) THEN
              IANGL=NBDR-1
              WRITE(IOUT,9000) NAMSBR,IREG,IANGL
            ENDIF
            DVNOR(IREG,IDIR)=DONE
            IF(IDIR .EQ. 1) THEN
              NBV0=NBV0+1
            ELSE
              NBV1=NBV1+1
            ENDIF
          ELSE
            DVNOR(IREG,IDIR)=SURVOL(IREG)/DVNOR(IREG,IDIR)
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
      DO ISUR=1,NFSUR
        IF(IPRINT .GE. 20 ) THEN
          WRITE(IOUT,6031) ISUR,
     >    SURVOL(-ISUR),DSNOR(ISUR)
        ENDIF
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
      IF(ABS(IPRINT) .GE. 1) THEN
        WRITE(IOUT,6005) DSVERR,DMVERR,DAVERR
        IF(IPRINT .GE. 10) THEN
          DO IREG=1,NFREG
            WRITE(IOUT,6010) IREG,SURVOL(IREG)
            WRITE(IOUT,6012) DVNOR(IREG,1),
     >      100.0D0*(DONE-DVNOR(IREG,1))
          ENDDO
        ENDIF
        WRITE(IOUT,6006) DSSERR,DMSERR,DASERR
        IF(IPRINT .GE. 10) THEN
          DO ISUR=1,NFSUR
            WRITE(IOUT,6011) -ISUR,SURVOL(-ISUR)
            WRITE(IOUT,6012) DSNOR(ISUR),
     >      100.0D0*(DONE-DSNOR(ISUR))
          ENDDO
        ENDIF
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Save track normalisation vector
*----
      IF(MAXSGL .EQ. 0) THEN
        WRITE(IOUT,9004) NAMSBR
        MAXSGL=LINMAX
      ENDIF
*----
*  Deallocate temporary records
*----
      DEALLOCATE(DCINT,ICINT,DLENGT,LENGTH,NUMERO)
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
 6010 FORMAT(' Normalisation factors and relative errors (%) ',
     >       'for region ',I8,' with volume ',F19.10)
 6011 FORMAT(' Normalisation factors and relative error (%) ',
     >       'for surface ',I8,' with area   ',F19.10)
 6012 FORMAT((2X,F15.10,2X,F10.5))
 6026 FORMAT('Line',I10.10,'={',I10,',',I10,'};')
 6027 FORMAT('Problem with Line',I10.10,'={',5(I10,','),
     >       I10,',',F18.10,',',I10,',',F18.10,'};')
 6030 FORMAT(' Normalization volumes  =',I10,1P,2D20.10)
 6031 FORMAT(' Normalization surfaces =',I10,1P,2D20.10)
 6050 FORMAT(' Direction = ',I5,2X,'Scan =',4I5,2x,'Weight =',F20.15/
     >' Directions and starting point =',6F20.15)
 6051 FORMAT('% Torig',I10.10,/'Torig=[')
 6052 FORMAT(9F15.8)
 6061 FORMAT(I10)
 6062 FORMAT(F18.10)
 6063 FORMAT('% ',I10.10,/A7,'=[')
 6064 FORMAT(18X,'];')
 6065 FORMAT(
     >'xcol=jet(    2);'/
     >'lls=length(Torig);'/'lli=length(LinePos);'/
     >'is=-1;'/'nums=0;'/
     >'for i=1:lli'/'  if LineReg(i) < 0'/'    if is == -1'/
     >'      nums=nums+1;'/'      is=0;'/
     >'      xi=Torig(nums,1);'/'      yi=Torig(nums,2);'/
     >'      tx=Torig(nums,4);'/'      ty=Torig(nums,5);'/
     >'    elseif is == 0'/'      is=-1;'/'    end'/'  else'/
     >'    xf=xi+tx*LinePos(i);'/'    yf=yi+ty*LinePos(i);'/
     >'    xxx=line([xi,xf],[yi,yf]);'/
     >'    xi=xf;'/'    yi=yf;'/'    reg=LineReg(i);'/
     >'    set(xxx,''Color'',[xcol(reg,:)]);'/
     >'  end'/'end'/'clear Torig LinePosLineReg'/'pause ;')
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
 9100 FORMAT(1X,' ***** Error in ',A6,'***** for line ',2I10/
     >       3F20.15)
 9101 FORMAT(1X,' ***** Error in ',A6,'***** for line ',I8/
     >       7X,'Positions (current and reference ) =',1P,2D21.14/
     >       7X,'Relative error = ',D21.14)
 9102 FORMAT(1X,' ***** Error in ',A6,'***** for line ',I8/
     >       7X,'Positions (current and reference ) =',1P,2D20.12/
     >       7X,'Absolute error = ',D21.14)
      END
