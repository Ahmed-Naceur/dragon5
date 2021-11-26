*DECK NXTLCA
      FUNCTION NXTLCA(IPRINT,ITST  ,NDIM  ,MXMESH,LINMAX,
     >                MESH  ,ORITRK,DIRTRK,DCMESH,
     >                NBCOR ,NBSINT,ISINT ,TRKLSI)
*
*----------
*
*Purpose:
* To track a Cartesian 2-D or 3-D geometry
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
*         $X$, $Y$ or $Z$.
* LINMAX  maximum number of segments in a track.
* MESH    effective number of spatial subdivision in
*         each direction ($X$, $Y$ and $Z$).
* ORITRK  a point on the track (origin).
* DIRTRK  the track direction (director cosines).
* DCMESH  spatial description of the parallepiped.
*
*Parameters: output
* NXTLCA  number of surfaces intersected.
* NBCOR   number of corner found for each external faces.
* NBSINT  number of surface crossed by track.
* ISINT   direction of plane intersected and
*         the surfaces crossed by the track.
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
      INTEGER          MESH(NDIM)
      DOUBLE PRECISION ORITRK(NDIM),
     >                 DIRTRK(NDIM),
     >                 DCMESH(-1:MXMESH,5)
      INTEGER          NBCOR(2),NBSINT
      INTEGER          ISINT(0:5,LINMAX)
      DOUBLE PRECISION TRKLSI(LINMAX)
      INTEGER          NXTLCA
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLCA')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-8,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IDIR,INEXT,IFRST,ILAST,IFACE,IOFF,NBTINT
      INTEGER          JDIR,JNEXT,JFRST,JLAST,JFACE,JOFF,JSUR
      INTEGER          KDIR,KNEXT,KFRST,KLAST,KFACE,KOFF,KSUR
      INTEGER          IC1,IDIRS,IDIRF,IDIRB,IRFIN(5)
      INTEGER          ISDIR,JSDIR,KSDIR
      DOUBLE PRECISION DP1,DP2,TRKK,TRKJ
      DOUBLE PRECISION TRKDIS,TRKOLD,DELTKD
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
      KFRST=1
      KLAST=1
      KNEXT=1
      TRKK=DZERO
      TRKOLD=DZERO
      KOFF=1
*----
*  Initialise output vectors
*----
      NBCOR(1)=0
      NBCOR(2)=0
      NBSINT=0
      CALL XDISET(ISINT ,6*LINMAX,0)
      CALL XDDSET(TRKLSI,LINMAX,DZERO)
*----
*  Print header if required
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6011) 'meshx={             '
        WRITE(IOUT,6012) (DCMESH(IC1,1),IC1=0,MESH(1))
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'meshy={             '
        WRITE(IOUT,6012) (DCMESH(IC1,2),IC1=0,MESH(2))
        WRITE(IOUT,6013)
        IF(NDIM .EQ. 3) THEN
          WRITE(IOUT,6011) 'meshz={             '
          WRITE(IOUT,6012) (DCMESH(IC1,3),IC1=0,MESH(3))
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
*  For option ITST=0, 1, consider global geometry
*  and test if line crosses this geometry.
*----
      IF(ITST .EQ. 0 .OR. ITST .EQ. 1) THEN
*----
*  Scan over Cartesian directions
*----
        DO 100 IDIR=1,NDIM
*----
*  IDIR is main planes direction ($X$, $Y$ or $Z$)
*  JDIR is first direction of plane perpendicular
*       to main direction ($Y, $Z$ or $X$).
*  KDIR is second direction of plane perpendicular
*       to main direction ($Z$, $X$ or $Y$).
*----
          JDIR=MOD(IDIR,NDIM)+1
          KDIR=MOD(IDIR+1,NDIM)+1
*----
*  Select planes order in direction IDIR
*----
          IF(DIRTRK(IDIR) .LT. 0) THEN
            INEXT=-1
            IFRST=MESH(IDIR)+1
            ILAST=1
          ELSE
            INEXT=1
            IFRST=1
            ILAST=MESH(IDIR)+1
          ENDIF
*----
*  Select planes order in direction JDIR
*----
          IF(DIRTRK(JDIR) .LT. 0) THEN
            JNEXT=-1
            JFRST=MESH(JDIR)+1
            JLAST=1
          ELSE
            JNEXT=1
            JFRST=1
            JLAST=MESH(JDIR)+1
          ENDIF
*----
*  Select planes order in direction KDIR (if required)
*----
          IF(NDIM .EQ. 3) THEN
            IF(DIRTRK(KDIR) .LT. 0) THEN
              KNEXT=-1
              KFRST=MESH(KDIR)+1
              KLAST=1
            ELSE
              KNEXT=1
              KFRST=1
              KLAST=MESH(KDIR)+1
            ENDIF
          ENDIF
*----
*  Scan over planes in direction IDIR
*----
          DO 101 IFACE=IFRST,ILAST,ILAST-IFRST
*----
*  Compute track length required to reach a face
*  and intersection point with infinite plane
*----
            TRKDIS=(DCMESH(IFACE-1,IDIR)-ORITRK(IDIR))/DIRTRK(IDIR)
            TRKJ=DBLE(JNEXT)*(TRKDIS*DIRTRK(JDIR)+ORITRK(JDIR))
            IF(NDIM .EQ. 3) THEN
              TRKK=DBLE(KNEXT)*(TRKDIS*DIRTRK(KDIR)+ORITRK(KDIR))
            ENDIF
*----
*  Test if point is in finite Cartesian plane (j,k)
*----
            DP1=DBLE(JNEXT)*DCMESH(JFRST-1,JDIR)
            DP2=DBLE(JNEXT)*DCMESH(JLAST-1,JDIR)
            IF(TRKJ .GE. DP1-DCUTOF
     >   .AND. TRKJ .LE. DP2+DCUTOF ) THEN
              IF(NDIM .EQ. 3) THEN
*----
*  For 3-D, consider intersection with plane
*----
                DP1=DBLE(KNEXT)*DCMESH(KFRST-1,KDIR)
                DP2=DBLE(KNEXT)*DCMESH(KLAST-1,KDIR)
                IF(TRKK .GE. DP1-DCUTOF
     >       .AND. TRKK .LE. DP2+DCUTOF ) THEN
*----
*  If no point, save distance set NBSINT to 1 and continue
*  if one point already known, check if this point is
*  at a different distance exit, otherwise, this point
*  is already identified (intersection of 2 or more planes)
*  just continue
*----
                  IF(NBSINT .EQ. 0) THEN
                    NBSINT=1
                    TRKOLD=TRKDIS
                  ELSE IF(TRKOLD .NE. TRKDIS) THEN
                    NBSINT=2
                    GO TO 105
                 ENDIF
                ENDIF
              ELSE
*----
*  For 2-D, intersection with a line is sufficient
*----
                IF(NBSINT .EQ. 0) THEN
                  NBSINT=1
                  TRKOLD=TRKDIS
                ELSE IF(TRKOLD .NE. TRKDIS) THEN
                  NBSINT=2
                  GO TO 105
                ENDIF
              ENDIF
            ENDIF
 101      CONTINUE
 100    CONTINUE
 105    CONTINUE
      ELSE
*----
*  If ITST=-1, assume that 2 surface are intersected
*  and continue
*----
        NBSINT=2
      ENDIF
      NXTLCA=NBSINT
*----
*  If line does not intersect any surfaces, exit
*----
      IF(NXTLCA .EQ. 0) THEN
        IF(IPRINT .GT. 1000) THEN
          WRITE(IOUT,6002) NAMSBR
        ENDIF
        RETURN
      ENDIF
*----
*  If line does not intersect 2 surfaces, abort
*----
      IF(NBSINT .NE. 2) THEN
        IF(NBSINT .EQ. 1) CALL XABORT(NAMSBR//
     >  ': Invalid line since it intersects only 1 surface')
        IF(IPRINT .GT. 1000)
     >   WRITE(IOUT,9000) NBSINT
      ENDIF
*----
*  If ITST=0, return number of surfaces intersected by line
*----
      IF(ITST .EQ. 0) RETURN
*----
*  If ITST=1 and NXTLCA=2 or ITST=-1
*  track geometry with submesh
*----
      NBSINT=0
*----
*  Scan over Cartesian directions
*----
      DO 110 IDIR=1,NDIM
*----
*  IDIR is main planes direction ($X$, $Y$ or $Z$)
*  JDIR is first direction of plane perpendicular
*       to main direction ($Y, $Z$ or $X$).
*  KDIR is second direction of plane perpendicular
*       to main direction ($Z$, $X$ or $Y$).
*----
        JDIR=MOD(IDIR,NDIM)+1
        KDIR=MOD(IDIR+1,NDIM)+1
*----
*  Select planes order in direction IDIR
*----
        IF(DIRTRK(IDIR) .LT. 0) THEN
          INEXT=-1
          IFRST=MESH(IDIR)+1
          ILAST=1
          IOFF=-1
        ELSE
          INEXT=1
          IFRST=1
          ILAST=MESH(IDIR)+1
          IOFF=0
        ENDIF
*----
*  Select planes order in direction JDIR
*----
        IF(DIRTRK(JDIR) .LT. 0) THEN
          JNEXT=-1
          JFRST=MESH(JDIR)+1
          JLAST=2
          JOFF=-1
        ELSE
          JNEXT=1
          JFRST=1
          JLAST=MESH(JDIR)
          JOFF=0
        ENDIF
*----
*  Select planes order in direction KDIR (if required)
*----
        IF(NDIM .EQ. 3) THEN
          IF(DIRTRK(KDIR) .LT. 0) THEN
            KNEXT=-1
            KFRST=MESH(KDIR)+1
            KLAST=2
            KOFF=-1
          ELSE
            KNEXT=1
            KFRST=1
            KLAST=MESH(KDIR)
            KOFF=0
          ENDIF
        ENDIF
*----
*  Scan over planes in direction IDIR
*----
        DO 111 IFACE=IFRST,ILAST,INEXT
          ISDIR=IDIR
          IF(IFACE .EQ. IFRST .OR. IFACE .EQ. ILAST)
     >    ISDIR=-ISDIR
*----
*  Compute track length required to reach a face
*  and intersection point with infinite plane
*----
          TRKDIS=(DCMESH(IFACE-1,IDIR)-ORITRK(IDIR))/DIRTRK(IDIR)
          TRKJ=DBLE(JNEXT)*(TRKDIS*DIRTRK(JDIR)+ORITRK(JDIR))
          IF(NDIM .EQ. 3) THEN
            TRKK=DBLE(KNEXT)*(TRKDIS*DIRTRK(KDIR)+ORITRK(KDIR))
          ENDIF
*----
*  Test if point is in finite Cartesian plane (j,k) and add to number
*  of surfaces crossed if this is the case
*  (j=1 or NPLJ+1 and k=1 or NPLK+1)
*  For values of j=1,NPLJ and k=1,NPLK a volume was crossed
*----
          DO 112 JFACE=JFRST,JLAST,JNEXT
            JSDIR=JDIR
            IF(JFACE .EQ. JFRST .OR. JFACE .EQ. JLAST)
     >      JSDIR=-JSDIR
            DP1=DBLE(JNEXT)*DCMESH(JFACE-1,JDIR)
            DP2=DBLE(JNEXT)*DCMESH(JFACE+JNEXT-1,JDIR)
            IF(TRKJ .GE. DP1-DCUTOF
     >   .AND. TRKJ .LE. DP2+DCUTOF ) THEN
              IF(NDIM .EQ. 3) THEN
*----
*  For 3-D, consider intersection with plane
*----
                DO 113 KFACE=KFRST,KLAST,KNEXT
                  KSDIR=KDIR
                  IF(KFACE .EQ. KFRST .OR. KFACE .EQ. KLAST)
     >            KSDIR=-KSDIR
                  DP1=DBLE(KNEXT)*DCMESH(KFACE-1,KDIR)
                  DP2=DBLE(KNEXT)*DCMESH(KFACE+KNEXT-1,KDIR)
                  IF(TRKK .GE. DP1-DCUTOF
     >         .AND. TRKK .LE. DP2+DCUTOF ) THEN
                    NBSINT=NBSINT+1
                    IF(IPRINT .GT. 1000) THEN
                      WRITE(IOUT,6011) 'Intersection at     '
                      WRITE(IOUT,6014) IDIR,IFACE,
     >                JDIR,JFACE+JOFF,KDIR,KFACE+KOFF
                      WRITE(IOUT,6015) TRKDIS,DIRTRK
                    ENDIF
                    IF(NBSINT .EQ. 1) THEN
                      ISINT(0,NBSINT)=ISDIR
                      ISINT(IDIR,NBSINT)=-IFACE
                      ISINT(JDIR,NBSINT)=JFACE+JOFF
                      ISINT(KDIR,NBSINT)=KFACE+KOFF
                      TRKLSI(NBSINT)=TRKDIS
                    ELSE
*----
*  Scan previous distances and reorder if necessary
*----
                      NBTINT=NBSINT
                      DO 120 JSUR=NBTINT-1,1,-1
                        DELTKD=TRKDIS-TRKLSI(JSUR)
                        IF(DELTKD .GE. -DCUTOF) THEN
*----
*  For more than one intersection at a given distance:
*  a) if this is an external surface, store the distance
*     for corner duplication
*  b) if this is an internal surface select adequate next region
*     to cross and reset NBSINT=NBSINT-1
*----
                          IF(ABS(DELTKD) .LE. DCUTOF) THEN
                            IF(ISDIR .GT. 0 ) THEN
                              IF(ISINT(0,JSUR) .GT. 0) THEN
                                ISINT(0,JSUR)=ISINT(0,JSUR)+10*ISDIR
                                ISINT(IDIR,JSUR)=-IFACE
                              ENDIF
                              NBSINT=NBSINT-1
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'SurInt A =',JSUR,IFACE,JFACE,KFACE,
*     >ISINT(0,JSUR),ISINT(IDIR,JSUR),ISINT(JDIR,JSUR),
*     >ISINT(KDIR,JSUR),
*     >TRKDIS,TRKLSI(JSUR)
                              GO TO 125
                            ELSE
                              IF(ISINT(0,JSUR) .GT. 0) THEN
                                ISINT(0,JSUR)=ISDIR
                                ISINT(IDIR,JSUR)=-IFACE
                                ISINT(JDIR,JSUR)=JFACE+JOFF
                                ISINT(KDIR,JSUR)=KFACE+KOFF
                                NBSINT=NBSINT-1
                                GO TO 125
                              ENDIF
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'SurInt B =',JSUR,IFACE,JFACE,KFACE,
*     >ISINT(0,JSUR),ISINT(IDIR,JSUR),ISINT(JDIR,JSUR),
*     >ISINT(KDIR,JSUR),
*     >TRKDIS,TRKLSI(JSUR)
                            ENDIF
                          ENDIF
                          DO 121 KSUR=NBTINT,JSUR+2,-1
                            ISINT(0,KSUR)=ISINT(0,KSUR-1)
                            ISINT(IDIR,KSUR)=ISINT(IDIR,KSUR-1)
                            ISINT(JDIR,KSUR)=ISINT(JDIR,KSUR-1)
                            ISINT(KDIR,KSUR)=ISINT(KDIR,KSUR-1)
                            TRKLSI(KSUR)=TRKLSI(KSUR-1)
 121                      CONTINUE
                          ISINT(0,JSUR+1)=ISDIR
                          ISINT(IDIR,JSUR+1)=-IFACE
                          ISINT(JDIR,JSUR+1)=JFACE+JOFF
                          ISINT(KDIR,JSUR+1)=KFACE+KOFF
                          TRKLSI(JSUR+1)=TRKDIS
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'SurInt C =',JSUR,IFACE,JFACE,KFACE,
*     >ISINT(0,JSUR),ISINT(IDIR,JSUR),ISINT(JDIR,JSUR),
*     >ISINT(KDIR,JSUR),
*     >TRKDIS,TRKLSI(JSUR)
                          GO TO 125
                        ENDIF
 120                  CONTINUE
                      DO 122 KSUR=NBTINT,2,-1
                        ISINT(0,KSUR)=ISINT(0,KSUR-1)
                        ISINT(IDIR,KSUR)=ISINT(IDIR,KSUR-1)
                        ISINT(JDIR,KSUR)=ISINT(JDIR,KSUR-1)
                        ISINT(KDIR,KSUR)=ISINT(KDIR,KSUR-1)
                        TRKLSI(KSUR)=TRKLSI(KSUR-1)
 122                  CONTINUE
                      ISINT(0,1)=ISDIR
                      ISINT(IDIR,1)=-IFACE
                      ISINT(JDIR,1)=JFACE+JOFF
                      ISINT(KDIR,1)=KFACE+KOFF
                      TRKLSI(1)=TRKDIS
 125                  CONTINUE
                    ENDIF
                  ENDIF
 113            CONTINUE
              ELSE
*----
*  For 2-D, intersection with a line is sufficient
*----
                NBSINT=NBSINT+1
                IF(IPRINT .GT. 1000) THEN
                  WRITE(IOUT,6011) 'Intersection at     '
                  WRITE(IOUT,6014) IDIR,IFACE,JDIR,JFACE
                  WRITE(IOUT,6015) TRKDIS,DIRTRK
                ENDIF
                IF(NBSINT .EQ. 1) THEN
                  ISINT(0,NBSINT)=ISDIR
                  ISINT(IDIR,NBSINT)=-IFACE
                  ISINT(JDIR,NBSINT)=JFACE+JOFF
                  TRKLSI(NBSINT)=TRKDIS
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'Int sur A =',NBSINT,
*     >ISINT(IDIR,NBSINT),ISINT(JDIR,NBSINT),TRKLSI(NBSINT)
                ELSE
*----
*  Scan previous distances and reorder if necessary
*----
                  NBTINT=NBSINT
                  DO 130 JSUR=NBTINT-1,1,-1
                    DELTKD=TRKDIS-TRKLSI(JSUR)
                    IF(DELTKD .GE. -DCUTOF) THEN
*----
*  For more than one intersection at a given distance:
*  a) if this is an external surface, store the distance
*     for corner duplication
*  b) if this is an internal surface select adequate next region
*     to cross and reset NBSINT=NBSINT-1
*----
                      IF(ABS(DELTKD) .LE. DCUTOF) THEN
                        IF(ISDIR .GT. 0 ) THEN
                          IF(ISINT(0,JSUR) .GT. 0) THEN
                            ISINT(0,JSUR)=ISINT(0,JSUR)+10*ISDIR
                            ISINT(IDIR,JSUR)=-IFACE
                          ENDIF
                          NBSINT=NBSINT-1
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'Int sur A =',JSUR,
*     >ISINT(IDIR,JSUR),ISINT(JDIR,JSUR),TRKLSI(JSUR)
                          GO TO 135
                        ELSE
                          IF(ISINT(0,JSUR) .GT. 0) THEN
                            ISINT(0,JSUR)=ISDIR
                            ISINT(IDIR,JSUR)=-IFACE
                            ISINT(JDIR,JSUR)=JFACE+JOFF
                            NBSINT=NBSINT-1
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'Int sur B =',JSUR,
*     >ISINT(IDIR,JSUR),ISINT(JDIR,JSUR),TRKLSI(JSUR)
                            GO TO 135
                          ENDIF
                        ENDIF
                      ENDIF
                      DO 131 KSUR=NBTINT,JSUR+2,-1
                        ISINT(0,KSUR)=ISINT(0,KSUR-1)
                        ISINT(IDIR,KSUR)=ISINT(IDIR,KSUR-1)
                        ISINT(JDIR,KSUR)=ISINT(JDIR,KSUR-1)
                        TRKLSI(KSUR)=TRKLSI(KSUR-1)
 131                  CONTINUE
                      ISINT(0,JSUR+1)=ISDIR
                      ISINT(IDIR,JSUR+1)=-IFACE
                      ISINT(JDIR,JSUR+1)=JFACE+JOFF
                      TRKLSI(JSUR+1)=TRKDIS
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'Int sur C =',JSUR+1,
*     >ISINT(IDIR,JSUR+1),ISINT(JDIR,JSUR+1),TRKLSI(JSUR+1)
                      GO TO 135
                    ENDIF
 130              CONTINUE
                  DO 132 KSUR=NBTINT,2,-1
                    ISINT(0,KSUR)=ISINT(0,KSUR-1)
                    ISINT(IDIR,KSUR)=ISINT(IDIR,KSUR-1)
                    ISINT(JDIR,KSUR)=ISINT(JDIR,KSUR-1)
                    TRKLSI(KSUR)=TRKLSI(KSUR-1)
 132              CONTINUE
                  ISINT(0,1)=ISDIR
                  ISINT(IDIR,1)=-IFACE
                  ISINT(JDIR,1)=JFACE+JOFF
                  TRKLSI(1)=TRKDIS
*      IF(IPRINT .GT. 1000)
*     >write(6,*) 'Int sur D =',1,
*     >ISINT(IDIR,1),ISINT(JDIR,1),TRKLSI(1)
 135              CONTINUE
                ENDIF
              ENDIF
            ENDIF
 112      CONTINUE
 111    CONTINUE
 110  CONTINUE
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6011) 'Intersections'
        DO 210 JSUR=1,NBSINT
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
 210    CONTINUE
      ENDIF
      IF(NBSINT .GE. 2) THEN
*----
*  Find number of corners for initial surface intersections
*----
        TRKDIS=TRKLSI(1)
        NBCOR(1)=1
        DO 140 JSUR=2,NBSINT
          DELTKD=TRKLSI(JSUR)-TRKDIS
          IF(DELTKD .GT. DCUTOF) GO TO 145
          NBCOR(1)=NBCOR(1)+1
 140    CONTINUE
 145    CONTINUE
*----
*  Find number of corners for final surface intersections
*----
        TRKDIS=TRKLSI(NBSINT)
        NBCOR(2)=1
        DO 150 JSUR=NBSINT-1,1,-1
          DELTKD=TRKDIS-TRKLSI(JSUR)
          IF(DELTKD .GT. DCUTOF) GO TO 155
          NBCOR(2)=NBCOR(2)+1
 150    CONTINUE
 155    CONTINUE
      ENDIF
*----
*  Print intersection points
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6011) 'Initial planes      '
        DO JSUR=1,NBCOR(1)
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
        ENDDO
        IF(NBCOR(1)+1 .LE. NBSINT-NBCOR(2)) THEN
          WRITE(IOUT,6011) 'Intermediate planes '
          DO JSUR=NBCOR(1)+1,NBSINT-NBCOR(2)
            WRITE(IOUT,6010) TRKLSI(JSUR),
     >      (ISINT(IDIR,JSUR),IDIR=0,5)
          ENDDO
        ENDIF
        WRITE(IOUT,6011) 'Final planes        '
        DO JSUR=NBSINT-NBCOR(2)+1,NBSINT
          WRITE(IOUT,6010) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5)
        ENDDO
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
        IF(IDIRS .GT. 10) CALL XABORT(NAMSBR//
     >': Final outer surface identified as an internal corner')
        IDIRF=ABS(ISINT(IDIRS,JSUR+1))
        IF(IDIRF .EQ. 1) THEN
          ISINT(IDIRS,JSUR+1)=-1
        ELSE IF(IDIRF .EQ. MESH(IDIRS)+1) THEN
          ISINT(IDIRS,JSUR+1)=-2
        ELSE
*----
* Print information for debug
*----
          WRITE(IOUT,6011) 'meshx={             '
          WRITE(IOUT,6012) (DCMESH(IC1,1),IC1=0,MESH(1))
          WRITE(IOUT,6013)
          WRITE(IOUT,6011) 'meshy={             '
          WRITE(IOUT,6012) (DCMESH(IC1,2),IC1=0,MESH(2))
          WRITE(IOUT,6013)
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6011) 'meshz={             '
            WRITE(IOUT,6012) (DCMESH(IC1,3),IC1=0,MESH(3))
            WRITE(IOUT,6013)
          ENDIF
          WRITE(IOUT,6011) 'trackorigin={       '
          WRITE(IOUT,6012) ORITRK
          WRITE(IOUT,6013)
          WRITE(IOUT,6011) 'trackdirection={    '
          WRITE(IOUT,6012) DIRTRK
          WRITE(IOUT,6013)
          WRITE(IOUT,6011) 'Initial planes      '
          DO KSUR=1,NBCOR(1)
            WRITE(IOUT,6010) TRKLSI(KSUR),
     >      (ISINT(IDIR,KSUR),IDIR=0,5)
          ENDDO
          IF(NBCOR(1)+1 .LE. NBSINT-NBCOR(2)) THEN
            WRITE(IOUT,6011) 'Intermediate planes '
            DO KSUR=NBCOR(1)+1,NBSINT-NBCOR(2)
              WRITE(IOUT,6010) TRKLSI(KSUR),
     >        (ISINT(IDIR,KSUR),IDIR=0,5)
            ENDDO
          ENDIF
          WRITE(IOUT,6011) 'Final planes        '
          DO KSUR=NBSINT-NBCOR(2)+1,NBSINT
            WRITE(IOUT,6010) TRKLSI(KSUR),
     >      (ISINT(IDIR,KSUR),IDIR=0,5)
          ENDDO
          WRITE(IOUT,9001) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5),IDIRF
          CALL XABORT(NAMSBR//': Invalid final '//CDIR(IDIRS)//
     >' directed surface')
        ENDIF
      ENDDO
*----
*  Identify regions
*----
      DO JSUR=NBSINT-NBCOR(2)+1,NBCOR(1)+1,-1
        DO IDIRS=1,5
          IRFIN(IDIRS)=ABS(ISINT(IDIRS,JSUR-1))
          IF(ISINT(IDIRS,JSUR)   .LT. 0 ) THEN
            IF(ISINT(IDIRS,JSUR-1) .LT. 0 ) THEN
              IRFIN(IDIRS)=-MAX(ISINT(IDIRS,JSUR-1),ISINT(IDIRS,JSUR))
            ENDIF
          ELSE
            IF(ISINT(IDIRS,JSUR-1) .LT. 0 ) THEN
              IRFIN(IDIRS)=ISINT(IDIRS,JSUR)
            ENDIF
          ENDIF
        ENDDO
        TRKLSI(JSUR)=TRKLSI(JSUR)-TRKLSI(JSUR-1)
        DO IDIRS=1,5
          IF(ISINT(IDIRS,JSUR)   .LT. 0 ) THEN
            ISINT(IDIRS,JSUR)=IRFIN(IDIRS)
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
        IF(IDIRS .GT. 10) CALL XABORT(NAMSBR//
     >': Initial outer surface identified as an internal corner')
        IDIRB=ABS(ISINT(IDIRS,JSUR))
        IF(IDIRB .EQ. 1) THEN
          ISINT(IDIRS,JSUR)=-1
        ELSE IF(IDIRB .EQ. MESH(IDIRS)+1) THEN
          ISINT(IDIRS,JSUR)=-2
        ELSE
*----
* Print information for debug
*----
          WRITE(IOUT,6011) 'meshx={             '
          WRITE(IOUT,6012) (DCMESH(IC1,1),IC1=0,MESH(1))
          WRITE(IOUT,6013)
          WRITE(IOUT,6011) 'meshy={             '
          WRITE(IOUT,6012) (DCMESH(IC1,2),IC1=0,MESH(2))
          WRITE(IOUT,6013)
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6011) 'meshz={             '
            WRITE(IOUT,6012) (DCMESH(IC1,3),IC1=0,MESH(3))
            WRITE(IOUT,6013)
          ENDIF
          WRITE(IOUT,6011) 'trackorigin={       '
          WRITE(IOUT,6012) ORITRK
          WRITE(IOUT,6013)
          WRITE(IOUT,6011) 'trackdirection={    '
          WRITE(IOUT,6012) DIRTRK
          WRITE(IOUT,6013)
          WRITE(IOUT,6011) 'Initial planes      '
          DO KSUR=1,NBCOR(1)
            WRITE(IOUT,6010) TRKLSI(KSUR),
     >      (ISINT(IDIR,KSUR),IDIR=0,5)
          ENDDO
          IF(NBCOR(1)+1 .LE. NBSINT-NBCOR(2)) THEN
            WRITE(IOUT,6011) 'Intermediate planes '
            DO KSUR=NBCOR(1)+1,NBSINT-NBCOR(2)
              WRITE(IOUT,6010) TRKLSI(KSUR),
     >        (ISINT(IDIR,KSUR),IDIR=0,5)
            ENDDO
          ENDIF
          WRITE(IOUT,6011) 'Final planes        '
          DO KSUR=NBSINT-NBCOR(2)+1,NBSINT
            WRITE(IOUT,6010) TRKLSI(KSUR),
     >      (ISINT(IDIR,KSUR),IDIR=0,5)
          ENDDO
         WRITE(IOUT,9001) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=0,5),IDIRB
          CALL XABORT(NAMSBR//': Invalid initial '//CDIR(IDIRS)//
     >' directed surface')
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
 6015 FORMAT(6(F25.16,2X))
 9000 FORMAT(' Warning : ',I10,' surfaces crossed')
 9001 FORMAT(1X,F25.16,7I10)
      END
