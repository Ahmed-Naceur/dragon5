*DECK XCWSCL
      SUBROUTINE XCWSCL(  NDIM, NSURX,  NVOL,  NBAN,   NRT, MSROD,
     >                   MAROD, NANGL,  DENS,IFTEMP,  IPRT, NCODE,
     >                  SWZERO,NRINFO,   RAN,  COTE, NRODS,  RODS,
     >                   NRODR,  RODR, MXSUB, MXSEG,  NXRI,   IMS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform specular tracking for 2-D square cluster.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G.Marleau
*
*Parameters: input
* NDIM    dimension of problem.
* NSURX   number of initial surface.
* NVOL    total number of regions.
* NBAN    number of concentric regions.
* NRT     number of rod types.
* MSROD   maximum number of subrod per rods.
* MAROD   maximum number of rod in any cluster.
* NANGL   number of integration angles.
* DENS    minimum parallel line trak density.
* IFTEMP  temporary tracking file unit.
* SWZERO  logical for specular tracking.
* IPRT    print level.
* NCODE   boundary type.
* NRINFO  type of concentric region:
*         NRINFO(1,IAN) = new region number;
*         NRINFO(2,IAN) = associated cluster;
*                       = 0 no cluster.
* RAN     radius/lattice side of region.
* COTE    Y dimension for rectangle.
* NRODS   integer description of rod  type:
*         NRODS(1,IRT) = number of rod;
*         NRODS(2,IRT) = number of subrods in rod;
*         NRODS(3,IRT) = associated region.
* RODS    description of rod of a given type:
*         RODS(1,IRT) = rod center radius;
*         RODS(2,IRT) = angle position of one rod.
* NRODR   subrod region.
* RODR    subrod radius.
* MXSUB   current maximum number of subtracks.
* MXSEG   current maximum track length.
* NXRI    annular region content multi-rod.
* IMS     surface merge.
*
*----------------------------------------------------------------------
*
      PARAMETER       (IUNOUT=6,PI=3.1415926535897932,EPS=1.E-5)
      CHARACTER        TEDATA*13
      INTEGER          NDIM,NSURX,NVOL,NBAN,NRT,MSROD,MAROD,NANGL,
     >                 IFTEMP,IPRT,NCODE(6),NRINFO(2,NBAN),
     >                 NRODS(3,NRT),NRODR(NRT),MXSUB,MXSEG,
     >                 INDS(2),NXRI(NRT,NBAN),IMS(6),IPER(2)
      LOGICAL          LINTER,LNEWP,SWZERO
      REAL             DENS,RAN(NBAN),COTE,RODS(2,NRT),RODR(MSROD,NRT)
      DOUBLE PRECISION DFACX,DFACY,SIDE(2),RCIRC,DENSP,DENLIN,
     >                 PROJ,PMAX,PMIN,DEPART,TRKPOS(2,2),ROTPOS(2,2),
     >                 TRKBEG(2,2),DIRBEG(2),RONEPS,ANGD,ANGC,RADC,
     >                 RADC2,WEIGHT,XPO
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NRSEG,NNSEG,KANGL
      REAL, ALLOCATABLE, DIMENSION(:) :: ATOP
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: RODP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SEGLEN,WGTANG,
     > DNSANG,PTSANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DANGLE
*----
*  SCRATCH STORAGE ALLOCATION
*   NRSEG : region crossed by track
*   NNSEG : region crossed by track (left)
*   SEGLEN: length of track
*   RODP  : rod position in cartesian geometry
*   ATOP  : number of rod between origin and rod1
*   DANGLE: integration angles
*   WGTANG: integration weight
*   DNSANG: integration densities
*   PTSANG: principal integration angles
*----
      ALLOCATE(NRSEG(MXSEG),NNSEG(MXSEG),KANGL(MXSUB))
      ALLOCATE(SEGLEN(MXSEG),RODP(2,MAROD,NRT,2),ATOP(NRT))
      ALLOCATE(DANGLE(NDIM,2,4*NANGL),WGTANG(4*NANGL),DNSANG(NANGL),
     > PTSANG(NANGL))
*----
*  DETERMINE INTEGRATION LIMITS FOR CLUSTER REGIONS
*----
      IF(IPRT.GE.1) THEN
        WRITE(IUNOUT,'(//1X,A20)') 'SPECULAR TRACKING   '
      ENDIF
      MLSEG=MXSEG/(2*NANGL)
      SIDE(1)=DBLE(RAN(NBAN))
      SIDE(2)=DBLE(COTE)
      IF( ABS((SIDE(1)-SIDE(2))/ABS(1)).GT.10.*EPS )THEN
         CALL XABORT('XCWSCL: AVAILABLE ONLY FOR SQUARE GEOMETRIES')
      ENDIF
      RCIRC=SQRT(SIDE(1)**2+SIDE(2)**2)
      SIDE(1)= SIDE(1)/RCIRC
      SIDE(2)= SIDE(2)/RCIRC
      NTAN=NBAN-1
      IF(IPRT.GT.0) THEN
        WRITE(IUNOUT,6000) NVOL,NSURX,NBAN,NRT
        WRITE(IUNOUT,6001)
        WRITE(IUNOUT,6002) (II,NRODS(1,II),NRODS(2,II),
     >                      NRODS(3,II),II=1,NRT)
        WRITE(IUNOUT,6003) NANGL,DENS
      ENDIF
*----
*  SET FLAG FOR SURFACE CROSSING
*    IPER(1) = X-PERIOD
*    IPER(2) = Y-PERIOD
*  VALUES ARE
*    IPER(I) = 1 FOR PERIODIC BC
*    IPER(I) = 2 FOR OTHER BC
*----
      IPER(1)=2
      IPER(2)=2
      IF( (NCODE(1) .EQ. 4) .AND. (NCODE(2) .EQ. 4) ) THEN
        IPER(1)=1
      ENDIF
      IF( (NCODE(3) .EQ. 4) .AND. (NCODE(4) .EQ. 4) ) THEN
        IPER(2)=1
      ENDIF
      IPERG=MIN(IPER(1),IPER(2))
      IF( SWZERO )THEN
        IFIN= NANGL-1
        IDEB= 0
        ISTRID=1
      ELSE
        IFIN= 2*NANGL
        IDEB= 0
        ISTRID=2
      ENDIF
      IANG=0
      DO 100 ITX= IDEB, IFIN, ISTRID
        INDS(1)= ITX
        ITY=IFIN-ITX
        INDS(2)=ITY
        IANG= IANG+1
        CALL XELTSA( NDIM, SIDE, INDS, DENSP, DANGLE(1,1,IANG))
*----
*  CHANGE DENLIN FOR HORIZONTAL & VERTICAL ANGLES
*----
        IF( (ITX .EQ. 0) .OR. (ITY .EQ. 0) )THEN
          DNSANG(IANG)=DBLE(DENS)
        ELSE
          DENLIN= DENSP / RCIRC
          NTRAC=MAX(1,INT(DBLE(DENS)/DENLIN+0.5D0))
          DNSANG(IANG)= DBLE(NTRAC) * DENLIN
        ENDIF
*----
*  COMPUTE NTRAK AND CHANGE DENS ACCORDING TO INPUT
*----
        PTSANG(IANG)= DANGLE(1,1,IANG)
 100  CONTINUE
      CALL XELTSW(SIDE,NANGL,PTSANG,WGTANG)
      IF( IPRT.GT.2 )THEN
        DO 110 IANG= 1, NANGL
          WRITE(IUNOUT,6004) IANG, DANGLE(1,1,IANG),WGTANG(IANG),
     >                     DNSANG(IANG),WGTANG(IANG)/DNSANG(IANG)
  110   CONTINUE
       ENDIF
*----
*  LOCALIZE CENTER OF REFERENCE ROD WITH RESPECT TO X-Y AXIS
*----
      DO 120 IRT=1,NRT
        IF(NRODS(3,IRT).GT.0) THEN
          NBROD=NRODS(2,IRT)
          DANGR=2.*PI/FLOAT(NRODS(1,IRT))
          IF(RODR(NBROD,IRT).GT.RODS(1,IRT)) THEN
            ATOP(IRT)=0.0
          ELSE
            ATOP(IRT)=(RODS(2,IRT)
     >             +ASIN(RODR(NBROD,IRT)/RODS(1,IRT)))/DANGR
          ENDIF
        ENDIF
 120  CONTINUE
      SIDE(1)=DBLE(RAN(NBAN))
      SIDE(2)=DBLE(COTE)
*----
*  COPY ANGLES AND DENSITIES ON TEMPORARY TRACKING FILE
*----
      DO 125 IANG=1,NANGL
        DANGLE(1,1,2*NANGL-IANG+1)=-DANGLE(1,1,IANG)
        DANGLE(2,1,2*NANGL-IANG+1)=DANGLE(2,1,IANG)
        WGTANG(2*NANGL-IANG+1)=WGTANG(IANG)
 125  CONTINUE
      DO 126 IANG=1,2*NANGL
        DANGLE(1,1,4*NANGL-IANG+1)=DANGLE(1,1,IANG)
        DANGLE(2,1,4*NANGL-IANG+1)=-DANGLE(2,1,IANG)
        WGTANG(4*NANGL-IANG+1)=WGTANG(IANG)
 126  CONTINUE
      WRITE(IFTEMP) ((DANGLE(IDIM,1,IANG),IDIM=1,NDIM),IANG=1,4*NANGL)
      WRITE(IFTEMP) (2.0D0/WGTANG(IANG),IANG=1,4*NANGL)
*----
*  PRINT TRACKING INFORMATION IF REQUIRED
*----
      NSOLMX=0
      IF((IPRT.GT.1).AND.(IPRT.LT.100))THEN
        WRITE(IUNOUT,'(/8H0ECHO = ,I3,27H SOLID ANGLES TO BE TRACKED)')
     >                 NANGL
        NSOLMX= MIN(9, NANGL/10)
        IREF1=0
        WRITE(IUNOUT,'(1X,10(I1,9X))') (IREF1, IZZ=0,NSOLMX)
        WRITE(IUNOUT,'(1X,10(I1,9X))') (MOD(IZZ,10), IZZ=0,NSOLMX)
        WRITE(IUNOUT,'(2H 0)')
        TEDATA='(1H+,TXXX,I1)'
      ENDIF
      NOTRAK= 0
*----
*  ANGULAR TRACK SWEEP
*----
      IXYN=0
      IXYR=0
      N0LSEG=0
      DO 130 IANG=1,NANGL
        DENLIN = DNSANG(IANG)
        DENSP   = 1.D0 / DENLIN
*----
*  PRINT TRACKING INFORMATION IF REQUIRED
*----
        IF((IPRT.GT.1).AND.(IPRT.LT.100))THEN
          IF( MOD(IANG,100) .EQ. 0 )THEN
            IREF1=IREF1+1
            NDEBS= NSOLMX+1
            NSOLMX=MIN(NDEBS+9, NANGL/10)
            WRITE(IUNOUT,'(1X,10(I1,9X))')(IREF1,IZZ=NDEBS,NSOLMX)
            WRITE(IUNOUT,'(1X,10(I1,9X))')
     >           (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
            WRITE(IUNOUT,'(2H 0)')
          ELSE
            WRITE(TEDATA(7:9),'(I3.3)') MOD(IANG,100) + 2
            WRITE(IUNOUT,TEDATA) MOD(IANG,10)
          ENDIF
        ENDIF
*----
*  LOCALIZE ROD POSITIONS WITH RESPECT TO 2 DIFFERENT ANGLES
*  POSSIBLE (+-COS(THETA),SIN(THETA))
*----
        ANGD=ATAN2(DANGLE(2,1,IANG),DANGLE(1,1,IANG))
        DO 300 IA=1,2
          DO 310 IRT=1,NRT
            IF(NRODS(3,IRT).GT.0) THEN
              DANGR=2.*PI/FLOAT(NRODS(1,IRT))
              ANGC=(ANGD/DANGR)-ATOP(IRT)
              IF(ANGC.GT.0.0) THEN
                IRDEP=INT(ANGC+0.9999)
              ELSE
                IRDEP=INT(ANGC)
              ENDIF
              ANGC=RODS(2,IRT)-ANGD+IRDEP*DANGR
              DO 320 IRD=1,NRODS(1,IRT)
                RODP(1,IRD,IRT,IA)=RODS(1,IRT)*REAL(COS(ANGC))
                RODP(2,IRD,IRT,IA)=RODS(1,IRT)*REAL(SIN(ANGC))
                ANGC=ANGC+DANGR
 320          CONTINUE
            ENDIF
 310      CONTINUE
          ANGD=PI-ANGD
 300    CONTINUE
*----
*  PROJECT THE 4 CORNERS OF SQUARE LOCATED AT
*  -SIDE(1)/2 < X < SIDE(1)/2 AND -SIDE(2)/2 < Y < SIDE(2)/2
*  ON LINE NORMAL TO TRACK DIRECTION
*----
        PMIN = +1.0D+50
        PMAX = -1.0D+50
        DFACX=1.0D0
        DO 150 IX=1,2
          DFACY=1.0D0
          DO 160 IY=1,2
            PROJ  = (SIDE(1)*DFACX*DANGLE(1,2,IANG)
     >            + SIDE(2)*DFACY*DANGLE(2,2,IANG))/2.0
            IF( PROJ.LT.PMIN ) PMIN = PROJ
            IF( PROJ.GT.PMAX ) PMAX = PROJ
            DFACY=-1.0D0*DFACY
 160      CONTINUE
          DFACX=-1.0D0*DFACX
 150    CONTINUE
*----
*  FIND NUMBER OF PARALLEL TRACK: NEAREST INTEGER +1 FOR SECURITY
*----
        NPOINT =NINT((PMAX-PMIN)*DENLIN)+1
        DEPART =0.5D0*(PMAX+PMIN-DBLE(NPOINT)*DENSP)
        DO 170 J = 1, 2
          TRKPOS(J,1)= DEPART*DANGLE(J,2,IANG)
          DANGLE(J,2,IANG)= DANGLE(J,2,IANG)*DENSP
 170    CONTINUE
        LNEWP=.TRUE.
*----
*  TRACK OVER 2*NPOINT PARALLEL TRACK FOR DIRECTION
*  TRACK AND REFLECTION
*----
        IXYF=0
        DO 180 IPOINT = 1,2*NPOINT
          NRIN=0
          IF(LNEWP)THEN
            NSUB=0
            IA=1
            NOTRAK=NOTRAK+1
            NSEG=0
            N0FSEG=1
            N0LSEG=MLSEG
            IF(IXYF.EQ.0) THEN
              DO 181 J=1,2
                TRKPOS(J,1)= TRKPOS(J,1) +DANGLE(J,2,IANG)
 181          CONTINUE
              IXYN=0
            ENDIF
          ELSE
            N0FSEG=N0LSEG+1
            N0LSEG=N0LSEG+MLSEG
          ENDIF
          DO 182 ISEG=N0FSEG,N0LSEG
            NRSEG(ISEG)=0
            NNSEG(ISEG)=0
            SEGLEN(ISEG)=0.0D0
 182      CONTINUE
          NLSEG=N0LSEG
          NFSEG=N0FSEG
*----
*  FIND EXTERNAL SURFACES CROSSED BY THIS TRACK
*----
          CALL XCWREC(DANGLE(1,1,IANG),SIDE,TRKPOS,LINTER,ROTPOS,
     >                INDS,IMS)
*----
* REJECT TRACK IF LINTER IS FALSE
*----
          IF(.NOT.LINTER) GO TO 183
*----
*  KEEP THE TRACK IF LINTER IS TRUE
*  A) SAVE INITIAL AND FINAL SURFACE INFORMATION
*----
          NRSEG(NFSEG)=-INDS(1)
          SEGLEN(NFSEG)=0.5D0
          NRSEG(NLSEG)=-INDS(2)
          SEGLEN(NLSEG)=0.5D0
          NFSEG=NFSEG+1
          NLSEG=NLSEG-1
*----
*  SAVE INFORMATION FOR INITIAL AND FINAL ANNULAR TRACKING
*----
          NRSEG(NFSEG)=NRINFO(1,NBAN)
          NNSEG(NFSEG)=NRIN
          SEGLEN(NFSEG)=ROTPOS(1,1)
          NRSEG(NLSEG)=NRIN
          NNSEG(NLSEG)=NRINFO(1,NBAN)
          SEGLEN(NLSEG)=ROTPOS(1,2)
          NLSEG=NLSEG-1
          NFSEG=NFSEG+1
          NRIN=NRINFO(1,NBAN)
*----
*  TRACK INSIDE ANNULAR REGIONS
*----
          RADC=ABS(ROTPOS(2,1))
          RADC2=RADC**2
          DO 210 IAN=NTAN,1,-1
            IF(RADC.GE.RAN(IAN)) GO TO 211
*----
*  LINE INTERSECT ANNULUS IAN
*----
            XPO=SQRT(RAN(IAN)**2-RADC2)
            NRSEG(NLSEG)=NRIN
            NNSEG(NFSEG+1)=NRIN
            SEGLEN(NLSEG)=XPO
            NLSEG=NLSEG-1
            NRIN=NRINFO(1,IAN)
            NFSEG=NFSEG+1
            NRSEG(NFSEG)=NRIN
            NNSEG(NLSEG+1)=NRIN
            SEGLEN(NFSEG)=-XPO
            IF(NRINFO(2,IAN).NE.0) THEN
*----
*  TRACK INSIDE RODS
*----
              DO 146 KRT=1,NRT
                JRT=NXRI(KRT,IAN)
                LRT=MOD(JRT,1000000)
                IF((JRT.GT.3000000).OR.
     >            ((JRT.GT.0).AND.(JRT.LT.1000000)) ) THEN
                  CALL XCWROD(NRIN,NRODS(1,LRT),NRODR(LRT),
     >                        RODR(1,LRT),RODP(1,1,LRT,IA),
     >                        ROTPOS(2,1),NFSEG,NLSEG,SEGLEN,NRSEG,
     >                        NNSEG)
                ELSE IF(JRT.EQ.0) THEN
                  GO TO 147
                ENDIF
 146          CONTINUE
 147          CONTINUE
              DO 143 KRT=1,NRT
                JRT=NXRI(KRT,IAN)
                IF(JRT.LT.0) THEN
                  IRT=-JRT
                  NXTR=NRODR(IRT)
                  DO 144 IRD=NRODS(2,IRT),1,-1
                    IF(RADC.GT.RODR(IRD,IRT)) GO TO 211
*----
*  LINE INTERSECT CENTERED ROD IRD
*----
                    XPO=SQRT(RODR(IRD,IRT)*RODR(IRD,IRT)-RADC2)
                    NRSEG(NLSEG)=NRIN
                    NNSEG(NFSEG+1)=NRIN
                    SEGLEN(NLSEG)=XPO
                    NLSEG=NLSEG-1
                    NRIN=NXTR
                    NXTR=NXTR-1
                    NFSEG=NFSEG+1
                    NRSEG(NFSEG)=NRIN
                    NNSEG(NLSEG+1)=NRIN
                    SEGLEN(NFSEG)=-XPO
 144              CONTINUE
                  GO TO 211
                ENDIF
 143          CONTINUE
            ENDIF
 210      CONTINUE
 211      CONTINUE
          IF( LNEWP )THEN
            IF(IXYF .EQ. 0) THEN
              IXYF=MOD(INDS(1)+1,2)+1
            ENDIF
            DO 250 J= 1, 2
              TRKBEG(J,IXYF)= TRKPOS(J,1)
              DIRBEG(J)= DANGLE(J,1,IANG)
 250        CONTINUE
          ELSE IF(IXYN .EQ. 0) THEN
            IXY=MOD(INDS(1)+1,2)+1
            IF(IXY.NE.IXYF) THEN
              IXYN=IXY
              DO 251 J= 1, 2
                TRKBEG(J,IXYN)= TRKPOS(J,1)
 251          CONTINUE
            ENDIF
          ENDIF
          IF(IPRT.GE.100) THEN
            WRITE(IUNOUT,6100) IANG,DANGLE(1,1,IANG),DANGLE(2,1,IANG),
     >                         IPOINT,INDS(1),(TRKPOS(II,1),II=1,2),
     >                         IPOINT,INDS(2),(TRKPOS(II,2),II=1,2)
          ENDIF
          NSUB=NSUB+1
          IF(NSUB.GT.MXSUB) CALL XABORT('XCWSCL: MXSUB OVERFLOW.')
          KANGL(NSUB)=IANG
*----
*  COMPRESS AND SORT TRACK VECTOR
*----
          ISRT=N0FSEG+1
          NSRT=MLSEG-2
          CALL XCWSRT(IPRT,NSRT,SEGLEN(ISRT),NRSEG(ISRT),
     >                NNSEG(ISRT),NTSEG)
          NOSEG=NSEG+NTSEG+2
          IF(IPRT.GE.200) THEN
            WRITE(IUNOUT,6101) ROTPOS(2,1),
     >        (SEGLEN(IIJJ),NRSEG(IIJJ),IIJJ=NSEG+2,NOSEG),
     >         SEGLEN(N0LSEG),NRSEG(N0LSEG)
          ENDIF
*----
*  CONVERT SEGMENT DIVISION TO SEGMENT LENGTH
*----
          DO 240 ISEG=NSEG+2,NOSEG-1
            SEGLEN(ISEG)=SEGLEN(ISEG+1)-SEGLEN(ISEG)
 240      CONTINUE
          SEGLEN(NOSEG)=SEGLEN(N0LSEG)
          NRSEG(NOSEG)=NRSEG(N0LSEG)
          IF(IPRT.GE.200) THEN
            WRITE(IUNOUT,6102) NOSEG-NSEG,
     >        (SEGLEN(IIJJ),NRSEG(IIJJ),IIJJ=NSEG+1,NOSEG)
          ENDIF
          NSEG=NOSEG
          N0LSEG=NSEG
*----
*  FOR TRANSLATION -> CHANGE TRACK STARTUP POINT
*  FOR REFLECTION  -> CHANGE TRACK DIRECTION
*----
          JINT=MOD(INDS(2)+1,2)+1
          KINT=MOD(INDS(2),2)+1
          IF(IPER(JINT) .EQ. 1)  THEN
            TRKPOS(JINT,2)=-TRKPOS(JINT,2)
          ELSE
            DANGLE(JINT,1,IANG)=-DANGLE(JINT,1,IANG)
            IA=MOD(IA,2)+1
          ENDIF
          RONEPS= 0.0D0
          DO 260 J= 1, 2
            TRKPOS(J,1)= TRKPOS(J,2)
            RONEPS= RONEPS + (TRKPOS(J,1)-TRKBEG(J,IXYF))**2
     >                     + (DANGLE(J,1,IANG)-DIRBEG(J))**2
 260      CONTINUE
          LNEWP= RONEPS.LT.EPS
          IF(LNEWP)THEN
*----
*  NOW, WRITE THE TRACK
*----
            WEIGHT= 0.25*WGTANG(IANG)/DNSANG(IANG)
            WRITE(IFTEMP) NSUB,NSEG, WEIGHT,
     >                   (KANGL(I),I=1,NSUB),
     >                   (NRSEG(I),I=1,NSEG),
     >                   (SEGLEN(I),I=1,NSEG)
            IF(IPRT.GE.300) THEN
               WRITE(IUNOUT,6103) NOTRAK,IANG,NSEG,
     >            (SEGLEN(I),NRSEG(I),I=1,NSEG)
            ENDIF
            IF(IPERG .EQ. 1) THEN
              IF(IXYN .EQ. IXYF) THEN
                TRKPOS(1,1)=TRKBEG(1,IXYR)
                TRKPOS(2,1)=TRKBEG(2,IXYR)
                DANGLE(IXYF,1,IANG)=ABS(DANGLE(IXYF,1,IANG))
                IXYF=0
              ELSE IF(IXYN .EQ. 0) THEN
                IF(IPER(IXYF).EQ.1) THEN
                  DANGLE(IXYF,1,IANG)=-DANGLE(IXYF,1,IANG)
                  IXYN=IXYF
                  IXYR=IXYF
                ELSE
                  IXYF=0
                ENDIF
              ELSE IF(IXYN .NE. IXYF) THEN
                IF(IPER(IXYF).EQ.1) THEN
                  IXYR=IXYN
                  TRKBEG(1,IXYR)=TRKBEG(1,IXYF)
                  TRKBEG(2,IXYR)=TRKBEG(2,IXYF)
                  DANGLE(IXYF,1,IANG)=-DANGLE(IXYF,1,IANG)
                  IXYN=IXYF
                ELSE IF(IPER(IXYN).EQ.1) THEN
                  DANGLE(IXYN,1,IANG)=-DANGLE(IXYN,1,IANG)
                  TRKPOS(1,1)=TRKBEG(1,IXYN)
                  TRKPOS(2,1)=TRKBEG(2,IXYN)
                  IXYR=IXYF
                  IXYF=IXYN
                ENDIF
              ENDIF
            ELSE
              IXYF=0
            ENDIF
          ENDIF
  183     CONTINUE
  180   CONTINUE
  130 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PTSANG,DNSANG,WGTANG,DANGLE)
      DEALLOCATE(ATOP,RODP,SEGLEN)
      DEALLOCATE(KANGL,NNSEG,NRSEG)
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(1X,'          TOTAL NUMBER OF REGIONS =',I10/
     >       1X,'       NUMBER OF INITIAL SURFACES =',I10/
     >       1X,'        NUMBER OF ANNULAR REGIONS =',I10/
     >       1X,'             NUMBER OF RODS TYPES =',I10)
 6001 FORMAT(1X,'  ROD TYPE',10X,'  NB. RODS',10X,
     >       'NB. SUBROD',10X,'IN ANNULUS')
 6002 FORMAT((1X,I10,10X,I10,10X,I10,10X,I10))
 6003 FORMAT(1X,'INTEGRATION PARAMETERS',/
     >       1X,'        NUMBER OF ANGLES =',I10,/
     >       1X,'  MINIMUM TRACK DENSITY  =',1P,E15.7)
 6004 FORMAT( 1X,I4,': COS=',F10.6,' WGT=',F10.6,' DNS=',F10.6,
     >       ' WGT/DEN=',F10.6)
 6100 FORMAT(//' *** TRACKING INFORMATION ***'/
     >       ' ANGLE(',I5,')           :',1P,2E15.7/
     >       ' START SURFACE  (',I5,') :',5X,I10,5X,1P,2E15.7/
     >       ' FINISH SURFACE (',I5,') :',5X,I10,5X,1P,2E15.7)
 6101 FORMAT(' INTERSECTION OF REGION AT NORMAL DISTANCE =',1P,E15.7/
     >       3(5X,E15.7,1X,I5))
 6102 FORMAT(' NUMBER OF SEGMENTS  ',I10/1P,3(5X,E15.7,1X,I5))
 6103 FORMAT(/' INFORMATION TO TRACKING FILE: ',
     >       ' TRACK NUMBER =',I5,2X,'IANG =',I5,2X,'NSEG =',I7/1P,
     >       3(5X,E15.7,1X,I5))
      END
