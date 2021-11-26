*DECK NXTLCY
      FUNCTION NXTLCY(IPRINT,ITST  ,NDIM  ,MXMESH,LINMAX,
     >                MESH  ,ORITRK,DIRTRK,DCMESH,IDIRC ,
     >                NBCOR ,NBSINT,ISINT ,TRKLSI)
*
*----------
*
*Purpose:
* To track an annular 2-D or 3-D geometry
* using the NXT tracking procedure.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* ITST    type of tracking, where:
*         =-1   only the exact geometry
*               is considered taking into account the
*               submesh in each direction;
*         = 0   only the global geometry
*               is considered without taking into account the
*               submesh in each direction;
*         = 1   both the global
*               geometry (as a first step) and the exact geometry
*               are considered taking into account the
*               submesh in each direction.
* NDIM    dimension of problem.
* MXMESH  maximum number of spatial subdivision in
*         $R$ and $X$, $Y$ or $Z$.
* LINMAX  maximum number of segments in a track.
* MESH    effective number of spatial subdivision in $X$
*         $Y$, $Z$ and $R$.
* ORITRK  a point on the track (origin).
* DIRTRK  the track direction (director cosines).
* DCMESH  spatial description of the cylinder.
* IDIRC   the direction of the cylinder axis.
*
*Parameters: output
* NXTLCY  number of surfaces intersected.
* NBCOR   number of corner found for each external faces.
* NBSINT  number of surface crossed by track.
* ISINT   the surfaces crossed by the track.
* TRKLSI  the surface intersection distance.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,ITST,NDIM,MXMESH,LINMAX
      INTEGER          MESH(5)
      DOUBLE PRECISION ORITRK(NDIM),
     >                 DIRTRK(NDIM),
     >                 DCMESH(-1:MXMESH,5)
      INTEGER          IDIRC
      INTEGER          NBCOR(2),NBSINT
      INTEGER          ISINT(0:5,LINMAX)
      DOUBLE PRECISION TRKLSI(LINMAX)
      INTEGER          NXTLCY
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLCY')
      INTEGER          MXDIM
      PARAMETER       (MXDIM=3)
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
      DOUBLE PRECISION DCUTOL
      PARAMETER       (DCUTOL=1.0D-11)
*----
*  Local variables
*----
      INTEGER          IDIR,INDIR,INEXT,IFRST,ILAST,IFACE,
     >                 IDGR,IDG1,IDG2,IDGP,NR,NP
      INTEGER          JFACE,JSUR,KNEXT,KSFRST,KSLAST,
     >                 KVFRST,KVLAST,KFACE,KSUR,KSUB
      INTEGER          IC1,IDIRS,IDIRF,IDIRB,IRFIN,IDIRFS
      DOUBLE PRECISION DP1,DP2,TRKK
      DOUBLE PRECISION TRKDIS,TRKOLD,TRKTMP
      DOUBLE PRECISION TORI(MXDIM),ROTANG(MXDIM)
      DOUBLE PRECISION XYCEN(2),RADIUS,CYLLIM(2),PROJ,PROJN,
     >                 YLOC,XORI,XLOC(2),CYLINT,TCINTD,XXINTD
      DOUBLE PRECISION DCUTOF
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Verify ITST option and reset to default value if invalid
*----
      IF(ITST .LT. -1 .OR. ITST .GT. 1) THEN
*----
*  Reset ITST=1 (complete analysis) if the value of ITST is invalid
*----
        ITST=1
      ENDIF
*----
*  Initialise output vectors
*----
      KSFRST=1
      KSLAST=1
      KVFRST=1
      KVLAST=1
      KNEXT=1
      KSUB=1
      TRKOLD=DZERO
      NBCOR(1)=0
      NBCOR(2)=0
      NBSINT=0
      CALL XDISET(ISINT,6*LINMAX,0)
      CALL XDDSET(TRKLSI,LINMAX,DZERO)
*----
*  IDG1 is first direction of plane perpendicular
*       to main direction ($Y, $Z$ or $X$).
*  IDG2 is second direction of plane perpendicular
*       to main direction ($Z$, $X$ or $Y$).
*  IDGP is main cylinder direction ($X$, $Y$ or $Z$)
*       for 2D case IDGP=3
*----
      IDGR=4
      IDGP=IDIRC
      IDG1=MOD(IDGP,3)+1
      IDG2=MOD(IDGP+1,3)+1
      NR=MESH(4)
      NP=MESH(IDGP)
*----
*  Print header if required
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6011) 'radial'//CDIR(IDG1)//CDIR(IDG2)//
     >  '={          '
        WRITE(IOUT,6012)  (DCMESH(IC1,IDGR),IC1=1,NR)
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'center'//CDIR(IDG1)//CDIR(IDG2)//
     >  '={          '
        WRITE(IOUT,6012)  DCMESH(-1,IDG1),DCMESH(-1,IDG2)
        WRITE(IOUT,6013)
        IF(NDIM .EQ. 3) THEN
          WRITE(IOUT,6011) 'axial'//CDIR(IDGP)//'={            '
          WRITE(IOUT,6012)  (DCMESH(IC1,IDGP),IC1=0,NP)
          WRITE(IOUT,6013)
        ENDIF
        WRITE(IOUT,6011) 'trackorigin={       '
        WRITE(IOUT,6012) ORITRK
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'trackdirection={    '
        WRITE(IOUT,6012) DIRTRK
        WRITE(IOUT,6013)
      ENDIF
*----
*  Select planes order in direction IDIR
*----
      IFRST=NR+1
      INEXT=-1
      ILAST=2
*----
*  Select planes order in direction KDIR (if required)
*----
      DCUTOF=DCUTOL*DCUTOL
      IF(NDIM .EQ. 3) THEN
        DCUTOF=DCUTOF*DCUTOL
        IF(DIRTRK(IDGP) .LT. 0) THEN
          KNEXT=-1
          KSFRST=NP+1
          KSLAST=1
          KVFRST=NP+1
          KVLAST=2
          KSUB=-1
        ELSE
          KNEXT=1
          KSFRST=1
          KSLAST=NP+1
          KVFRST=1
          KVLAST=NP
          KSUB=0
        ENDIF
      ENDIF
*----
*  Projection on 2-D plane
*  and translation to origin of circle
*----
      XYCEN(1)=DCMESH(-1,IDG1)
      XYCEN(2)=DCMESH(-1,IDG2)
      TORI(1)=ORITRK(IDG1)
      TORI(2)=ORITRK(IDG2)
      IF(NDIM .EQ. 3) THEN
        TORI(3)=ORITRK(IDGP)
        ROTANG(3)=DIRTRK(IDGP)
      ELSE
        ROTANG(3)=DZERO
        TORI(3)=DZERO
      ENDIF
      PROJ=DONE/SQRT(DONE-ROTANG(3)**2)
      ROTANG(1)=DIRTRK(IDG1)*PROJ
      ROTANG(2)=DIRTRK(IDG2)*PROJ
      PROJN=DONE/SQRT(ROTANG(1)*ROTANG(1)+ROTANG(2)*ROTANG(2))
      ROTANG(1)=ROTANG(1)*PROJN
      ROTANG(2)=ROTANG(2)*PROJN
*----
*  Translate (x,y) coordinates to circle center and
*  rotate y coordinate in such a way that tracking
*  line parallel to the x-axis
*----
      XORI= (TORI(1)-XYCEN(1))*ROTANG(1)
     >     +(TORI(2)-XYCEN(2))*ROTANG(2)
      YLOC=-(TORI(1)-XYCEN(1))*ROTANG(2)
     >     +(TORI(2)-XYCEN(2))*ROTANG(1)
*      IF(IPRINT .GT. 1000) THEN
*        write(6,*) 'XYCEN =',XYCEN
*        write(6,*) 'TORI  =',TORI
*        write(6,*) 'ROTANG=',ROTANG
*        write(6,*) 'RotOri=',XORI,YLOC
*      ENDIF
*----
*  For option ITST=0, 1, consider global geometry
*  and test if line crosses this geometry.
*----
      IF(ITST .EQ. 0 .OR. ITST .EQ. 1) THEN
        NBSINT=0
        RADIUS=DCMESH(NR,IDGR)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'YLOC, RADIUS=',YLOC,RADIUS
        IF(NDIM .EQ. 3) THEN
          CYLLIM(1)=DCMESH(0,IDGP)
          CYLLIM(2)=DCMESH(NP,IDGP)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*)'CYLLIM ',CYLLIM
        ENDIF
*----
*  Test if coordinate is inside circle.
*----
        IF(ABS(YLOC) .LT. RADIUS) THEN
*----
*  Line crosses circle (infinite cylinder)
*----
          XLOC(2)=SQRT(RADIUS**2-YLOC**2)
          XLOC(1)=-XLOC(2)
          IF(NDIM .EQ. 3) THEN
*----
*  Find intersection points for finite cylinders in 3D
*----
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'XLOC ',XLOC(1),XLOC(2)
              DO 100 IFACE=1,2
              TCINTD=YLOC*ROTANG(1)+
     >               XLOC(IFACE)*ROTANG(2)+XYCEN(2)-TORI(2)
              XXINTD=XLOC(IFACE)*ROTANG(1)-
     >               YLOC*ROTANG(2)+XYCEN(1)-TORI(1)
              TRKTMP=(TCINTD*PROJ)/ROTANG(2)
              TRKDIS=SQRT(XXINTD*XXINTD+TCINTD*TCINTD)*PROJ
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'Changedir=',TRKDIS,TCINTD,ROTANG(2)
              IF(ROTANG(2)*TCINTD .LT. DZERO) TRKDIS=-TRKDIS
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'TRKDIS=',TRKDIS,TRKTMP,ROTANG(2)
              CYLINT=TRKDIS*ROTANG(3)+TORI(3)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'CYLINT= ',CYLINT
              IF(CYLINT .GE. CYLLIM(1) .AND. CYLINT .LE. CYLLIM(2)) THEN
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'iface r ',IFACE
                IF(NBSINT .EQ. 0) THEN
                  NBSINT=1
                  TRKOLD=TRKDIS
                ELSE IF(TRKOLD .NE. TRKDIS) THEN
                  NBSINT=2
                  GO TO 105
                ENDIF
              ENDIF
 100        CONTINUE
          ELSE
            NBSINT=2
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'Intersected in 2D'
            GO TO 105
          ENDIF
        ENDIF
*----
*  For 3-D geometry look for intersection with bottom
*  and top of cylinder
*----
        IF(NDIM .EQ. 3) THEN
          DO 102 IFACE=1,2
            TRKDIS=(CYLLIM(IFACE)-TORI(3))/ROTANG(3)
            XLOC(1)=TORI(1)+TRKDIS*ROTANG(1)/PROJ-XYCEN(1)
            XLOC(2)=TORI(2)+TRKDIS*ROTANG(2)/PROJ-XYCEN(2)
            TCINTD=SQRT(XLOC(1)**2+XLOC(2)**2)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'dir 3 int ', XLOC(1),XLOC(2),TRKDIS
            IF(TCINTD .LT. RADIUS) THEN
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'iface z ',IFACE
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) TORI(1)+TRKDIS*ROTANG(1)/PROJ,
*     >                 TORI(2)+TRKDIS*ROTANG(2)/PROJ,
*     >                 TORI(3)+TRKDIS*ROTANG(3)
              IF(NBSINT .EQ. 0) THEN
                NBSINT=1
                TRKOLD=TRKDIS
              ELSE IF(TRKOLD .NE. TRKDIS) THEN
                NBSINT=2
                GO TO 105
              ENDIF
            ENDIF
 102      CONTINUE
        ENDIF
 105    CONTINUE
      ELSE
*----
*  If ITST=-1, assume that 2 surface are intersected
*  and continue
*----
        NBSINT=2
      ENDIF
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'Nb sur=',NBSINT
      NXTLCY=NBSINT
*----
*  If line does not intersect any surfaces, exit
*----
      IF(NXTLCY .EQ. 0) THEN
        IF(IPRINT .GT. 1000) THEN
          WRITE(IOUT,6002) NAMSBR
        ENDIF
        RETURN
      ENDIF
*----
*  If line does not intersect 2 surfaces, abort
*----
      IF(NBSINT .NE. 2) CALL XABORT(NAMSBR//
     >': Invalid line since it did not intersect 0 or 2 surfaces')
*----
*  If ITST=0, return number of surfaces intersected by line
*----
      IF(ITST .EQ. 0) RETURN
*----
*  If ITST=1 or -1
*  track geometry with submesh
*----
      NBSINT=0
*----
*  Loop over radial regions from maximal radius
*  from outer to inner
*----
      DO 110 IFACE=IFRST,ILAST,INEXT
        RADIUS=DCMESH(IFACE-1,IDGR)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'YLOC, RADIUS=',YLOC,RADIUS
        IF(ABS(YLOC) .LE. RADIUS) THEN
          XLOC(2)=SQRT(RADIUS**2-YLOC**2)
          XLOC(1)=-XLOC(2)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'XLOC',XLOC
*----
*  Scan over two possible cylinder intersections
*----
          INDIR=1
          DO 111 JFACE=1,2
            TCINTD=YLOC*ROTANG(1)+XLOC(JFACE)*ROTANG(2)+XYCEN(2)-TORI(2)
            XXINTD=XLOC(JFACE)*ROTANG(1)-YLOC*ROTANG(2)+XYCEN(1)-TORI(1)
            TRKDIS=SQRT(XXINTD*XXINTD+TCINTD*TCINTD)*PROJ
*      IF(IPRINT .GT. 1000)
*     >  TRKTMP=(TCINTD*PROJ)/ROTANG(2)
            IF(ROTANG(2)*TCINTD .LT. -DCUTOF) TRKDIS=-TRKDIS
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'Distances',XXINTD,TCINTD,ROTANG(2),TRKDIS,TRKTMP
            IF(NDIM .EQ. 3) THEN
*----
*  Scan over planes in axial cylinder direction
*----
              TRKK=DBLE(KNEXT)*(TRKDIS*ROTANG(3)+TORI(3))
              DO 112 KFACE=KVFRST,KVLAST,KNEXT
                DP1=DBLE(KNEXT)*DCMESH(KFACE-1,IDGP)
                DP2=DBLE(KNEXT)*DCMESH(KFACE-1+KNEXT,IDGP)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'DP1,DP2 ',DP1,DP2,TRKK
                IF(TRKK .GE. DP1 .AND.
     >             TRKK .LE. DP2 ) THEN
                  NBSINT=NBSINT+1
                  IF(IPRINT .GT. 1000) THEN
                    WRITE(IOUT,6011) 'Intersection at     '
                      WRITE(IOUT,6014) IFACE,JFACE,KFACE
                      WRITE(IOUT,6015) TRKDIS
                  ENDIF
                  IF(NBSINT .EQ. 1) THEN
                    ISINT(0,NBSINT)=(IDG1+3)*INDIR
                    ISINT(IDGR,NBSINT)=IFACE-1
                    ISINT(IDGP,NBSINT)=KFACE+KSUB
                    TRKLSI(NBSINT)=TRKDIS
                  ELSE
*----
*  Scan previous distances and reorder if necessary
*----
                    DO 120 JSUR=NBSINT-1,1,-1
                      IF(TRKDIS .GE. TRKLSI(JSUR)) THEN
*----
*  For more than one intersection at a given distance:
*  a) if this is an external surface, store the distance
*     for corner duplication
*  b) if this is an internal surface select adequate next region
*     to cross and reset NBSINT=NBSINT-1
*----
                        IF(TRKDIS .EQ. TRKLSI(JSUR) .AND.
     >                     IFACE .NE. IFRST             ) THEN
                          ISINT(IDGR,JSUR)=IFACE-1
                          ISINT(IDGP,JSUR)=KFACE+KSUB
                          NBSINT=NBSINT-1
                          GO TO 125
                        ENDIF
                        DO 121 KSUR=NBSINT,JSUR+2,-1
                          ISINT(0,KSUR)=ISINT(0,KSUR-1)
                          ISINT(IDGR,KSUR)=ISINT(IDGR,KSUR-1)
                          ISINT(IDGP,KSUR)=ISINT(IDGP,KSUR-1)
                          TRKLSI(KSUR)=TRKLSI(KSUR-1)
 121                    CONTINUE
                        ISINT(0,JSUR+1)=(IDG1+3)*INDIR
                        ISINT(IDGR,JSUR+1)=IFACE-1
                        ISINT(IDGP,JSUR+1)=KFACE+KSUB
                        TRKLSI(JSUR+1)=TRKDIS
                        GO TO 125
                      ENDIF
 120                CONTINUE
                    DO 122 KSUR=NBSINT,2,-1
                      ISINT(0,KSUR)=ISINT(0,KSUR-1)
                      ISINT(IDGR,KSUR)=ISINT(IDGR,KSUR-1)
                      ISINT(IDGP,KSUR)=ISINT(IDGP,KSUR-1)
                      TRKLSI(KSUR)=TRKLSI(KSUR-1)
 122                CONTINUE
                    ISINT(0,1)=(IDG1+3)*INDIR
                    ISINT(IDGR,1)=IFACE-1
                    ISINT(IDGP,1)=KFACE+KSUB
                    TRKLSI(1)=TRKDIS
 125                CONTINUE
                  ENDIF
                ENDIF
 112          CONTINUE
            ELSE
*----
*  For 2-D, intersection with a circle
*----
              NBSINT=NBSINT+1
              IF(IPRINT .GT. 1000) THEN
                WRITE(IOUT,6011) 'Intersection at     '
                WRITE(IOUT,6014) IFACE,JFACE
                WRITE(IOUT,6015) TRKDIS
              ENDIF
              IF(NBSINT .EQ. 1) THEN
                ISINT(0,NBSINT)=(IDG1+3)*INDIR
                ISINT(IDGR,NBSINT)=IFACE-1
                TRKLSI(NBSINT)=TRKDIS
              ELSE
*----
*  Scan previous distances and reorder if necessary
*----
                DO 130 JSUR=NBSINT-1,1,-1
                  IF(TRKDIS .GE. TRKLSI(JSUR)) THEN
*----
*  For more than one intersection at a given distance:
*  a) if this is an external surface, store the distance
*     for corner duplication
*  b) if this is an internal surface select adequate next region
*     to cross and reset NBSINT=NBSINT-1
*----
                    IF(TRKDIS .EQ. TRKLSI(JSUR) .AND.
     >                 IFACE .NE. IFRST              ) THEN
                      ISINT(IDGR,JSUR)=IFACE-1
                      NBSINT=NBSINT-1
                      GO TO 135
                    ENDIF
                    DO 131 KSUR=NBSINT,JSUR+2,-1
                      ISINT(0,KSUR)=ISINT(0,KSUR-1)
                      ISINT(IDGR,KSUR)=ISINT(IDGR,KSUR-1)
                      TRKLSI(KSUR)=TRKLSI(KSUR-1)
 131                CONTINUE
                    ISINT(0,JSUR+1)=(IDG1+3)*INDIR
                    ISINT(IDGR,JSUR+1)=IFACE-1
                    TRKLSI(JSUR+1)=TRKDIS
                    GO TO 135
                  ENDIF
 130            CONTINUE
                DO 132 KSUR=NBSINT,2,-1
                  ISINT(0,KSUR)=ISINT(0,KSUR-1)
                  ISINT(IDGR,KSUR)=ISINT(IDGR,KSUR-1)
                  TRKLSI(KSUR)=TRKLSI(KSUR-1)
 132            CONTINUE
                ISINT(0,1)=(IDG1+3)*INDIR
                ISINT(IDGR,1)=IFACE-1
                TRKLSI(1)=TRKDIS
 135            CONTINUE
              ENDIF
            ENDIF
            INDIR=-1
 111      CONTINUE
        ENDIF
 110  CONTINUE
      IF(NDIM .EQ. 3) THEN
*----
*  Loop over axial regions from Top to bottom or
*  vice versa
*----
        DO 140 KFACE=KSFRST,KSLAST,KNEXT
          TRKDIS=(DCMESH(KFACE-1,IDGP)-TORI(3))/ROTANG(3)
          DP1=TORI(1)+TRKDIS*ROTANG(1)/PROJ-XYCEN(1)
          DP2=TORI(2)+TRKDIS*ROTANG(2)/PROJ-XYCEN(2)
          TCINTD=SQRT(DP1**2+DP2**2)
*----
*  Loop over radial regions from inner to outer
*----
          DO 141 IFACE=ILAST,IFRST,-INEXT
            RADIUS=DCMESH(IFACE-1,IDGR)
*      IF(IPRINT .GT. 1000)
*     >  write(6,*) 'zz= ',DP1,DP2,TCINTD,RADIUS
            IF(TCINTD .LE. RADIUS) THEN
              IF(IPRINT .GT. 1000) THEN
                WRITE(IOUT,6011) 'Intersection at     '
                WRITE(IOUT,6014) KFACE,IFACE
                WRITE(IOUT,6015) TRKDIS
              ENDIF
              NBSINT=NBSINT+1
              IF(NBSINT .EQ. 1) THEN
                ISINT(0,NBSINT)=IDGP
                ISINT(IDGR,NBSINT)=IFACE-1
                ISINT(IDGP,NBSINT)=KFACE
                TRKLSI(NBSINT)=TRKDIS
              ELSE
*----
*  Scan previous distances and reorder if necessary
*----
                DO 150 JSUR=NBSINT-1,1,-1
                  IF(TRKDIS .GE. TRKLSI(JSUR)) THEN
*----
*  For more than one intersection at a given distance:
*  a) if this is an external surface, store the distance
*     for corner duplication
*  b) if this is an internal surface select adequate next region
*     to cross and reset NBSINT=NBSINT-1
*----
                    IF(TRKDIS .EQ. TRKLSI(JSUR) .AND.
     >                 KFACE .NE. KSFRST         .AND.
     >                 KFACE .NE. KSLAST              ) THEN
                      ISINT(IDGR,JSUR)=IFACE-1
                      ISINT(IDGP,JSUR)=KFACE
                      NBSINT=NBSINT-1
                      GO TO 155
                    ENDIF
                    DO 151 KSUR=NBSINT,JSUR+2,-1
                      ISINT(0,KSUR)=ISINT(0,KSUR-1)
                      ISINT(IDGR,KSUR)=ISINT(IDGR,KSUR-1)
                      ISINT(IDGP,KSUR)=ISINT(IDGP,KSUR-1)
                      TRKLSI(KSUR)=TRKLSI(KSUR-1)
 151                CONTINUE
                    ISINT(0,JSUR+1)=IDGP
                    ISINT(IDGR,JSUR+1)=IFACE-1
                    ISINT(IDGP,JSUR+1)=KFACE
                    TRKLSI(JSUR+1)=TRKDIS
                    GO TO 155
                  ENDIF
 150            CONTINUE
                DO 152 KSUR=NBSINT,2,-1
                  ISINT(0,KSUR)=ISINT(0,KSUR-1)
                  ISINT(IDGR,KSUR)=ISINT(IDGR,KSUR-1)
                  ISINT(IDGP,KSUR)=ISINT(IDGP,KSUR-1)
                  TRKLSI(KSUR)=TRKLSI(KSUR-1)
 152            CONTINUE
                ISINT(0,1)=IDGP
                ISINT(IDGR,1)=IFACE-1
                ISINT(IDGP,1)=KFACE
                TRKLSI(1)=TRKDIS
 155            CONTINUE
              ENDIF
              GO TO 145
            ENDIF
 141      CONTINUE
 145      CONTINUE
 140    CONTINUE
      ENDIF
      IF(NBSINT .GE. 2) THEN
*----
*  Find number of corners for initial surface intersections
*----
        TRKDIS=TRKLSI(1)
        NBCOR(1)=1
        DO 160 JSUR=2,NBSINT
          IF(TRKLSI(JSUR) .GT. TRKDIS) GO TO 165
          NBCOR(1)=NBCOR(1)+1
 160    CONTINUE
 165    CONTINUE
*----
*  Find number of corners for final surface intersections
*----
        TRKDIS=TRKLSI(NBSINT)
        NBCOR(2)=1
        DO 170 JSUR=NBSINT-1,1,-1
          IF(TRKLSI(JSUR) .LT. TRKDIS) GO TO 175
          NBCOR(2)=NBCOR(2)+1
 170    CONTINUE
 175    CONTINUE
      ENDIF
*----
*  Print intersection points
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6011) 'Initial planes      '
        DO 200 JSUR=1,NBCOR(1)
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
 200    CONTINUE
        IF(NBCOR(1)+1 .LE. NBSINT-NBCOR(2)) THEN
          WRITE(IOUT,6011) 'Intermediate planes '
          DO 201 JSUR=NBCOR(1)+1,NBSINT-NBCOR(2)
            WRITE(IOUT,6010) TRKLSI(JSUR),
     >      (ISINT(IDIR,JSUR),IDIR=0,5)
 201      CONTINUE
        ENDIF
        WRITE(IOUT,6011) 'Final planes        '
        DO 202 JSUR=NBSINT-NBCOR(2)+1,NBSINT
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
 202    CONTINUE
      ENDIF
*----
*  Identify final faces
*----
      DO JSUR=NBSINT,NBSINT-NBCOR(2)+1,-1
        DO IDIRS=0,5
          ISINT(IDIRS,JSUR+1)=ISINT(IDIRS,JSUR)
        ENDDO
        TRKLSI(JSUR+1)=TRKLSI(JSUR)
        IDIRS=ABS(ISINT(0,JSUR+1))
        IF(IDIRS .GE. 4) THEN
          IDIRF=ISINT(4,JSUR+1)
          IF(IDIRF .NE. MESH(4)) THEN
            WRITE(IOUT,9001) NAMSBR,'final   ',
     >                    JSUR,4,IDIRF,MESH(4)
            CALL XABORT(NAMSBR//
     >      ': Invalid '//CDIR(4)//' directed surface')
          ENDIF
          ISINT(4,JSUR+1)=-2
        ELSE
          IDIRF=ISINT(IDIRS,JSUR+1)
          IF(IDIRF .EQ. 1) THEN
            ISINT(IDIRS,JSUR+1)=-1
          ELSE IF(IDIRF .EQ. MESH(IDIRS)+1) THEN
            ISINT(IDIRS,JSUR+1)=-2
          ELSE
            WRITE(IOUT,9001) NAMSBR,'final   ',
     >                    JSUR,IDIRS,IDIRF,MESH(IDIRS)
            CALL XABORT(NAMSBR//': Invalid '//CDIR(IDIRS)//
     >      ' directed surface')
          ENDIF
        ENDIF
      ENDDO
*----
*  Identify regions
*----
      DO JSUR=NBSINT-NBCOR(2)+1,NBCOR(1)+1,-1
        IDIRFS=ISINT(0,JSUR)
        IF(IDIRFS .GE. 4) THEN
          IDIRFS=1
        ELSE
          IDIRFS=0
        ENDIF
        IDIRB=ABS(ISINT(0,JSUR-1))
        IF(IDIRB .GE. 4) IDIRB=4
        IDIRF=ABS(ISINT(0,JSUR))
        IF(IDIRF .GE. 4) IDIRF=4
        IF(IDIRB .EQ. IDIRF) THEN
          IF(IDIRB .EQ. 4) THEN
            IRFIN=MAX(ISINT(4,JSUR-1),ISINT(4,JSUR))
          ELSE
            IRFIN=MIN(ISINT(IDIRF,JSUR-1),ISINT(IDIRF,JSUR))
          ENDIF
        ELSE
          IF(IDIRF .EQ. 4) THEN
            IF(IDIRFS .EQ. 1) THEN
              IRFIN=ISINT(4,JSUR-1)
            ELSE
              IRFIN=ISINT(4,JSUR)
            ENDIF
          ELSE
            IRFIN=ISINT(IDIRF,JSUR-1)
          ENDIF
        ENDIF
        IF(IPRINT .GT. 1000) THEN
          WRITE(IOUT,*) 'JSUR =',JSUR,ISINT(0,JSUR-1),ISINT(0,JSUR),
     >                  IDIRB,IDIRF,IRFIN
        ENDIF
        TRKLSI(JSUR)=TRKLSI(JSUR)-TRKLSI(JSUR-1)
        DO IDIRS=1,5
          IF(IDIRS .EQ. IDIRF) THEN
            ISINT(IDIRS,JSUR)=IRFIN
          ELSE
            ISINT(IDIRS,JSUR)=ISINT(IDIRS,JSUR)
          ENDIF
        ENDDO
      ENDDO
*----
*  Identify initial face
*----
      DO JSUR=1,NBCOR(1)
        IDIRS=ABS(ISINT(0,JSUR))
        IF(IDIRS .GE. 4) THEN
          IDIRB=ISINT(4,JSUR)
          IF(IDIRB .NE. MESH(4)) THEN
            WRITE(IOUT,9001) NAMSBR,'initial ',
     >                    JSUR,4,IDIRB,MESH(4)
            CALL XABORT(NAMSBR//
     >      ': Invalid '//CDIR(4)//' directed surface')
          ENDIF
          ISINT(4,JSUR)=-2
        ELSE
          IDIRB=ISINT(IDIRS,JSUR)
          IF(IDIRB .EQ. 1) THEN
            ISINT(IDIRS,JSUR)=-1
          ELSE IF(IDIRB .EQ. MESH(IDIRS)+1) THEN
            ISINT(IDIRS,JSUR)=-2
          ELSE
            WRITE(IOUT,9001) NAMSBR,'initial ',
     >                    JSUR,IDIRS,IDIRB,MESH(IDIRS)
            CALL XABORT(NAMSBR//': Invalid '//CDIR(IDIRS)//
     >      ' directed surface')
          ENDIF
        ENDIF
      ENDDO
*----
*  Print final track information
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6011) 'Initial face        '
        DO JSUR=1,NBCOR(1)
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
        ENDDO
        WRITE(IOUT,6011) 'Regions             '
        DO JSUR=NBCOR(1)+1,NBSINT-NBCOR(2)+1
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
        ENDDO
        WRITE(IOUT,6011) 'Final face          '
        DO JSUR=NBSINT-NBCOR(2)+2,NBSINT+1
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT('   No region crossed'/
     >       '   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,F25.16,6I10)
 6011 FORMAT(A20)
 6012 FORMAT(6(1X,F25.16,:,','))
 6013 FORMAT('};')
 6014 FORMAT(6I10)
 6015 FORMAT(6(1X,F25.16))
*----
*  Warning formats
*----
 9001 FORMAT(1X,'Error in --',A6,' --'/
     >       1X,'invalid ',A8,' surface identification ',4I10)
      END
