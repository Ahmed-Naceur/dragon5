*DECK NXTPRA
      FUNCTION NXTPRA(NFACES,POSCAR,POSANN,POSPIN,VOLINT)
*
*----------
*
*Purpose:
* Compute the volume of intersection between
* a 2-D Cartesian region defined by N planes
* an 2-D annular region and
* an annular pin centered at the origin.
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
* NFACES  number of planes for Cartesian geometre (3 for triangles,
*         4 for rectangle s and 6 for hexagones).
* POSCAR  Cartesian region corner definition:
*         POSCAR(1,*) is X position;
*         POSCAR(2,*) is Y position;
*         POSCAR(*,IPLANE) is location of first corner of
*         plane IPLANE;
*         POSCAR(*,IPLANE+1) is location of second corner of
*         plane IPLANE. 
*         For last plane, the position of the
*         second corner is POSCAR(*,1)
* POSANN  spatial description of the annular region with
*         POSANN(0) the radius, POSANN(1) the $X$ position
*         of center and POSANN(2) the $Y$ position
*         of center.
* POSPIN  spatial description of the annular pin region with
*         POSPIN(0) the radius, POSPIN(1) the $X$ position
*         of center and POSPIN(2) the $Y$ position
*         of center.
*
*Parameters: output
* NXTPRA  type of intersection between the three regions, where:
*         = 0 means that the volume of intersection
*         between the three regions vanishes;
*         =-1 means that the volume of intersection
*         between the three regions was computed.
* VOLINT  2-D volume of intersection (area) between the three regions.
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
      INTEGER          NXTPRA
      INTEGER          NFACES
      DOUBLE PRECISION POSCAR(2,NFACES),POSANN(0:2),POSPIN(0:2),VOLINT
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI,PIO2,TPIO2,FPIO2
*----
*  Local parameters
*----
      INTEGER          IOUT,IPRINT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,IPRINT=1,NAMSBR='NXTPRA')
      DOUBLE PRECISION DZERO,DONE,DTWO,DHALF,CUTOFF
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0,DHALF=0.5D0,
     >                 CUTOFF=1.0D-10)
*----
*  Local variables
*----
      DOUBLE PRECISION PX,PY,RADAN2,RACEN2,RACEN,RADPI2,
     >                 COSA,SINA,XIAPR,YIAPR,TRANX,TRANY,
     >                 XYIAPD(2,2)
      INTEGER          IPOINT,IFACE,IPLIN,JPLIN,NPLIN,ICYL,
     >                 INT,INTT
      DOUBLE PRECISION XLOC,YLOC,XMIN,XMAX,COSR,SINR,ACARG
      DOUBLE PRECISION POSCYL(2),RADCYL,R2CYL,YBOT,SOL,
     >                 XINT(2),XLOCR,YLOCR
      DOUBLE PRECISION XYBEG(2),RADBEG,THBEG,XYEND(2),RADEND,THEND,
     >                 XYADD(2,2),HDT,FACT,DVOL
      INTEGER          NPOINT,NBPTS,NADSEG,TYADD(2,2),ISEG
      DOUBLE PRECISION PNTINT(6),POINTS(2,24)
      INTEGER          TYPINT(6),TYPES(2,24)
      DOUBLE PRECISION X1,Y1,X2,Y2,VLEN,DP1,DP2,DDP
      INTEGER          NPIN,IP1,IP2
      DOUBLE PRECISION ANGR

*----
*  Initialize NXTPRA to no intersection
*  Initialize PI and multiples
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
        write(IOUT,6100)
        WRITE(IOUT,6101) (POSCAR(1,IP1),POSCAR(2,IP1),IP1=1,NFACES)
        WRITE(IOUT,6102) POSCAR(1,1),POSCAR(2,1)
        WRITE(IOUT,6103) POSANN(0),
     >                   POSANN(1),POSANN(2)
        WRITE(IOUT,6104) POSPIN(0),
     >                   POSPIN(1),POSPIN(2)
      ENDIF
      COSR=DZERO
      SINR=DZERO
      XMAX=DZERO
      XMIN=DZERO
      YLOC=DZERO
      NXTPRA=0
      PI=XDRCST('Pi',' ')
      PIO2=PI/DTWO
      TPIO2=3.0D0*PIO2
      FPIO2=5.0D0*PIO2
*----
*  Locate annular region/annular pin intersection points
*  and find transformation matrix to locate intersection points with
*  respect to center of annular region/annular pin intersection.
*  In this system of reference the intersection points are located at
*  $(0,y_{i})$ and $(0,y_{i})$ or $\theta=\pm \pi/2$ respectively.
*  This will become usefull when all the intersection points between
*  the rectangular region, the annular region and the annular pins 
*  must be classified according to a counter-clockwise order.
*  See validation in mathematica file: PRA.nb
*----
      RADAN2=POSANN(0)*POSANN(0)
      PX=POSANN(1)-POSPIN(1)
      PY=POSANN(2)-POSPIN(2)
      RACEN2=PX*PX+PY*PY
      RACEN=SQRT(RACEN2)
      COSA=PX/RACEN
      SINA=PY/RACEN
      RADPI2=POSPIN(0)*POSPIN(0)
      XIAPR=(RADPI2+RACEN2-RADAN2)/(DTWO*RACEN)
      YIAPR=SQRT(RADPI2-XIAPR*XIAPR)
      TRANX=-XIAPR
      TRANY=DZERO
*----
*  Rotate points back to local frame of reference
*  XYIAPD(*,1) is at $\theta = \pi/2$
*  XYIAPD(*,2) is at $\theta = -\pi/2$
*----
      XYIAPD(1,1)=COSA*XIAPR-SINA*YIAPR+POSPIN(1)
      XYIAPD(2,1)=SINA*XIAPR+COSA*YIAPR+POSPIN(2)
      XYIAPD(1,2)=COSA*XIAPR+SINA*YIAPR+POSPIN(1)
      XYIAPD(2,2)=SINA*XIAPR-COSA*YIAPR+POSPIN(2)
*----
*  For each of the faces of the Cartesian region, 
*  classify the intersection between the annular region,
*  annular pins and corners in
*  increasing order in a counter-clockwise fashion.
*  A maximum of 6 intersection points can be found.
*  Note that the first and last corners are identified by $\pm 3$
*  respectively, the first and last intersection point with the
*  annular region are indicated by $\pm 1$ respectively
*  and the first and last intersection point with the
*  annular pin are indicated by $\pm 2$ respectively.
*  Procedure to classify the intersection points:
*  1 - Rotate faces in such a way that they are parallel to the
*      $X$ axis and below the Cartesian region. 
*  2 - Fill in corner locations in increasing order
*  3 - Locate intersection with annular region (after rotation) and
*      insert at adequate location in intersection point vector.
*  4 - Locate intersection with annular pin (after rotation) and
*      insert at adequate location in intersection point vector.
*----
      IPOINT=0
      DO IFACE=1,NFACES
        IP1=IFACE
        IP2=MOD(IFACE,NFACES)+1
        ANGR=ATAN2(POSCAR(2,IP2)-POSCAR(2,IP1),
     >             POSCAR(1,IP2)-POSCAR(1,IP1))
        COSR=COS(-ANGR)
        SINR=SIN(-ANGR)
*----
*  Left triangles
*----
*        write(6,'(A20,4F20.15)')
*     >  'Angles de rotation  ',ANGR,180.0*ANGR/PI,COSR,SINR
        XMIN=COSR*POSCAR(1,IP1)-SINR*POSCAR(2,IP1)
        XMAX=COSR*POSCAR(1,IP2)-SINR*POSCAR(2,IP2)
        YLOC=SINR*POSCAR(1,IP1)+COSR*POSCAR(2,IP1)
*        write(6,'(3F20.15)') XMIN,XMAX,YLOC
*----
*  Save corner location
*----
        IPLIN=1
        PNTINT(IPLIN)=XMIN
        TYPINT(IPLIN)=4
        IPLIN=2
        PNTINT(IPLIN)=XMAX
        TYPINT(IPLIN)=4
        NPLIN=IPLIN
*----
*  Loop over cylinder
*  1-  annular region
*  2-  annular pin
*----
        DO ICYL=1,2
*----
*  Extract cylinder information
*----
          IF(ICYL .EQ. 1) THEN
*----
*  Cylinder is annular region
*  Rotate as required.
*----
            POSCYL(1)=COSR*POSANN(1)-SINR*POSANN(2)
            POSCYL(2)=SINR*POSANN(1)+COSR*POSANN(2)
            RADCYL=POSANN(0)
          ELSE
*----
*  Cylinder is annular pin
*----
            POSCYL(1)=COSR*POSPIN(1)-SINR*POSPIN(2)
            POSCYL(2)=SINR*POSPIN(1)+COSR*POSPIN(2)
            RADCYL=POSPIN(0)
          ENDIF
*----
*  Find intersection points between Cartesian face and
*  cylindrical region
*----
          R2CYL=RADCYL*RADCYL
          YBOT=YLOC-POSCYL(2)
          SOL=R2CYL-YBOT*YBOT
          IF(SOL .GE. DZERO) THEN
            XINT(1)=POSCYL(1)-SQRT(SOL)
            XINT(2)=POSCYL(1)+SQRT(SOL)
*----
*  Classify intersection points per order of increasing
*  x location for annular region and annular pin
*----
            DO INT=1,2
              DO IPLIN=1,NPLIN
                IF(XINT(INT) .LE. PNTINT(IPLIN)) THEN
                  DO JPLIN=NPLIN,IPLIN,-1
                    PNTINT(JPLIN+1)=PNTINT(JPLIN)
                    TYPINT(JPLIN+1)=TYPINT(JPLIN)
                  ENDDO
                  PNTINT(IPLIN)=XINT(INT)
                  TYPINT(IPLIN)=ICYL
                  GO TO 100
                ENDIF
              ENDDO
              IPLIN=NPLIN+1
              PNTINT(IPLIN)=XINT(INT)
              TYPINT(IPLIN)=ICYL
 100          CONTINUE
              NPLIN=NPLIN+1
            ENDDO
          ENDIF
        ENDDO
*----
*  All intersection points located and ordered for this face
*  of the Cartesian region.
*  Scan and locate those defining the intersection of three
*  regions
*  sum of TYPINT = 7 namely:
*   +1 -> one annular region
*   +2 -> one pin crossing
*   +4 -> inside rectangle
*  TYPES(1,*) is type of line segment before point
*  TYPES(2,*) is type of line segment after point
*  Here
*  TYPES=1 means annular region,
*  TYPES=2 means annular pin and
*  TYPES=4 means rectangle side
*----
*        write(6,*) 'AVANT  NPLIN',NPLIN
*        write(6,'(F20.15,I10)')
*     >  (PNTINT(IPLIN),TYPINT(IPLIN),IPLIN=1,NPLIN)
        INTT=0
        DO IPLIN=1,NPLIN
*----
*  Rotate back line to original location.
*----
          XLOC=PNTINT(IPLIN)
          XLOCR=XLOC*COSR+YLOC*SINR
          YLOCR=-XLOC*SINR+YLOC*COSR
          IF(INTT .EQ. 7) THEN
*----
*  Already in 3 region intersection
*  find the point at which one leaves this region
*----
            IPOINT=IPOINT+1
            POINTS(1,IPOINT)=XLOCR
            POINTS(2,IPOINT)=YLOCR
            TYPES(1,IPOINT)=4
            TYPES(2,IPOINT)=TYPINT(IPLIN)
          ENDIF
          INTT=INTT+TYPINT(IPLIN)
          IF(INTT .EQ. 7) THEN
            IPOINT=IPOINT+1
            POINTS(1,IPOINT)=XLOCR
            POINTS(2,IPOINT)=YLOCR
            TYPES(1,IPOINT)=TYPINT(IPLIN)
            TYPES(2,IPOINT)=4
*----
*  Test if new point is at the same location as previous point
*  for rectangle corners and get rid of duplicates
*----
            IF(IPOINT .GE. 2) THEN
*            write(6,'(A8,I10)') 'IPOINT  ',IPOINT
*            write(6,'(A8,2I10,2F20.15)') 'CURRENT ',
*     >      TYPES(1,IPOINT),TYPES(2,IPOINT),
*     >      POINTS(1,IPOINT),POINTS(2,IPOINT)
*            write(6,'(A8,2I10,2F20.15)') 'PREVIOUS',
*     >      TYPES(1,IPOINT-1),TYPES(2,IPOINT-1),
*     >      POINTS(1,IPOINT-1),POINTS(2,IPOINT-1)
              IF(TYPES(1,IPOINT) .EQ. 4) THEN
                IF(TYPES(1,IPOINT-1) .EQ. 4 .AND.
     >             TYPES(2,IPOINT-1) .EQ. 4) THEN
                  DP1=POINTS(1,IPOINT-1)-POINTS(1,IPOINT)
                  DP2=POINTS(2,IPOINT-1)-POINTS(2,IPOINT)
                  DDP=SQRT(DP1*DP1+DP2*DP2)
*                   write(6,*) 'DP1,DP2,DDP',DP1,DP2,DDP
                  IF(DDP .LT. CUTOFF) THEN
*                  IF(POINTS(1,IPOINT-1) .EQ.  POINTS(1,IPOINT) .AND.
*     >               POINTS(2,IPOINT-1) .EQ.  POINTS(2,IPOINT) ) THEN
                    IPOINT=IPOINT-1
                  ELSE
                    CALL XABORT(NAMSBR//
     >              ': Problem with corner position')
                  ENDIF
                ELSE
                  CALL XABORT(NAMSBR//
     >            ': Problem with corner order')
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
*        write(6,*) 'APRES  NPLIN',NPLIN
*        write(6,'(F20.15,I10)')
*     >  (PNTINT(IPLIN),TYPINT(IPLIN),IPLIN=1,NPLIN)
      ENDDO
      NPOINT=IPOINT
*----
*  Complete segment for geometry
*----
*      write(6,*) 'NPOINT',NPOINT
      IF(NPOINT .LE. 1) THEN
*----
*  Path is empty or contain a single point for intersections with
*  sides of the Cartesian region:
*  A- Add annular/pin intersection points (XYIAPD) if both inside
*     Cartesian region and create path with two arc segments
*     1- From pin to annulus with arc in annular region (1)
*     2- From annulus to pin with arc in pin region (2)
*     3- Closed loop (0)
*  B- otherwise, there is no intersection
*----
        NPIN=0
        DO IPOINT=1,2
          DO IFACE=1,NFACES
            IP1=IFACE
            IP2=MOD(IP1,NFACES)+1
            X1=POSCAR(1,IP2)-POSCAR(1,IP1)
            Y1=POSCAR(2,IP2)-POSCAR(2,IP1)
            VLEN=SQRT(X1*X1+Y1*Y1)
            X2=-Y1/VLEN
            Y2=X1/VLEN
            VLEN=(XYIAPD(1,IPOINT)-POSCAR(1,IP1))*X2+
     >           (XYIAPD(2,IPOINT)-POSCAR(2,IP1))*Y2
            IF(VLEN. LT. DZERO) GO TO 101
          ENDDO
          NPIN=NPIN+1
 101      CONTINUE
        ENDDO  
        IF(NPIN .EQ. 2) THEN
*----
*  TYPES(1,*) is type of line segment before point
*  TYPES(2,*) is type of line segment after point
*  where TYPES=1 means annular region and
*  TYPES=2 means annular pin
*----
          NPOINT=NPOINT+3
          POINTS(1,1)=XYIAPD(1,1)
          POINTS(2,1)=XYIAPD(2,1)
          TYPES(1,1)  =2
          TYPES(2,1)  =1
          POINTS(1,2)=XYIAPD(1,2)
          POINTS(2,2)=XYIAPD(2,2)
          TYPES(1,2)  =1
          TYPES(2,2)  =2
          POINTS(1,3)=XYIAPD(1,1)
          POINTS(2,3)=XYIAPD(2,1)
          TYPES(1,2)  =2
          TYPES(2,3)  =1
        ELSE
          NPOINT=0
        ENDIF
      ELSE
*----
*  Test for cyclic track if first point is a corner
*----
        IPOINT=1
*           write(6,'(A8,I10)') 'IPOINT  ',IPOINT
*            write(6,'(A8,2I10,2F20.15)') 'CURRENT ',
*     >      TYPES(1,IPOINT),TYPES(2,IPOINT),
*     >      POINTS(1,IPOINT),POINTS(2,IPOINT)
*            write(6,'(A8,2I10,2F20.15)') 'LAST    ',
*     >      TYPES(1,NPOINT),TYPES(2,NPOINT),
*     >      POINTS(1,NPOINT),POINTS(2,NPOINT)
        IF(TYPES(1,IPOINT) .EQ. 4 .AND.
     >     TYPES(2,IPOINT) .EQ. 4) THEN
          IF(TYPES(1,NPOINT) .EQ. 4 .AND.
     >       TYPES(2,NPOINT) .EQ. 4) THEN
            DP1=POINTS(1,NPOINT)-POINTS(1,IPOINT)
            DP2=POINTS(2,NPOINT)-POINTS(2,IPOINT)
            DDP=SQRT(DP1*DP1+DP2*DP2)
*            write(6,*) 'DP1,DP2,DDP',DP1,DP2,DDP
            IF(DDP .GE. CUTOFF) THEN
*            IF(POINTS(1,NPOINT) .NE.  POINTS(1,IPOINT) .OR.
*     >         POINTS(2,NPOINT) .NE.  POINTS(2,IPOINT) ) THEN
              CALL XABORT(NAMSBR//
     >        ': Problem with end corner position')
            ENDIF
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Problem with end corner order')
          ENDIF
        ELSE
*----
*  Duplicate first point for cyclic track
*----
          NPOINT=NPOINT+1
          POINTS(1,NPOINT)=POINTS(1,IPOINT)
          POINTS(2,NPOINT)=POINTS(2,IPOINT)
          TYPES(1,NPOINT)=TYPES(1,IPOINT)
          TYPES(2,NPOINT)=TYPES(2,IPOINT)
        ENDIF
        IF(IPRINT .GE. 200) THEN
          WRITE(IOUT,6015)
          DO IPOINT=1,NPOINT
            IF(IPOINT .EQ. NPOINT) THEN
              WRITE(IOUT,6011) POINTS(1,IPOINT),POINTS(2,IPOINT),
     >                         TYPES(1,IPOINT),TYPES(2,IPOINT)
            ELSE
              WRITE(IOUT,6012) POINTS(1,IPOINT),POINTS(2,IPOINT),
     >                         TYPES(1,IPOINT),TYPES(2,IPOINT)
            ENDIF
          ENDDO
        ENDIF
*----
*  Add missing arc segment if required
*----
        NBPTS=NPOINT
        DO IPOINT=NPOINT,2,-1
          NADSEG=0
          IF(TYPES(1,IPOINT) .NE. 4) THEN
*----
*  This point finishes an arc segment
*  previous point must begin an arc
*----
            IF(TYPES(2,IPOINT-1) .EQ. 4) CALL XABORT(NAMSBR//
     >      ': Starting point for arc not found')
*----
*  Find position of intersection points with respect to
*  annular/pin center location and angular location
*  Rotate to center annular region on $X_{+}$ axis (COSA,SINA)
*  and translate by (-XIAPR,0) to center
*  annular region/annular pin at $x=0$
*----
            X1=POINTS(1,IPOINT-1)-POSPIN(1)
            Y1=POINTS(2,IPOINT-1)-POSPIN(2)
*            XYBEG(1)=COSA*POINTS(1,IPOINT-1)+SINA*POINTS(2,IPOINT-1)
*     >              -XIAPR
*            XYBEG(2)=-SINA*POINTS(1,IPOINT-1)+COSA*POINTS(2,IPOINT-1)
            XYBEG(1)=COSA*X1+SINA*Y1-XIAPR
            XYBEG(2)=-SINA*X1+COSA*Y1
            RADBEG=SQRT(XYBEG(1)*XYBEG(1)+XYBEG(2)*XYBEG(2))
            X2=POINTS(1,IPOINT)-POSPIN(1)
            Y2=POINTS(2,IPOINT)-POSPIN(2)
*            XYEND(1)=COSA*POINTS(1,IPOINT)+SINA*POINTS(2,IPOINT)
*     >              -XIAPR
*            XYEND(2)=-SINA*POINTS(1,IPOINT)+COSA*POINTS(2,IPOINT)
            XYEND(1)=COSA*X2+SINA*Y2-XIAPR
            XYEND(2)=-SINA*X2+COSA*Y2
            RADEND=SQRT(XYEND(1)*XYEND(1)+XYEND(2)*XYEND(2))
*----
*  Find angular location of points
*----
            ACARG=XYBEG(1)/RADBEG
            IF(ACARG .GE. 1.0D0) THEN
              THBEG=ACOS(1.0D0)
            ELSE IF(ACARG .LE. -1.0D0) THEN
              THBEG=ACOS(-1.0D0)
            ELSE
              THBEG=ACOS(ACARG)
            ENDIF
            IF(XYBEG(2) .LT. DZERO) THBEG=-THBEG
            ACARG=XYEND(1)/RADEND
            IF(ACARG .GE. 1.0D0) THEN
              THEND=ACOS(1.0D0)
            ELSE IF(ACARG .LE. -1.0D0) THEN
              THEND=ACOS(-1.0D0)
            ELSE
              THEND=ACOS(ACARG)
            ENDIF
            IF(XYEND(2) .LT. DZERO) THEND=-THEND
            IF(THEND .LT. THBEG) THEND=DTWO*PI+THEND
            IF(THBEG .LT. -PIO2) THEN
*----
*  For $\theta_{i}\le -\pi/2$ the segment must be of
*  type 1 (annular region)
*----
              IF(TYPES(2,IPOINT-1) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Initial line segment must be an annular region')
*----
*  For $\theta_{f}\le -\pi/2$ the segment must be of
*  type 1 (annular region) and there is no segment
*  to add
*  For $-\pi/2 < \theta_{f}\le \pi/2$ the segment must be of
*  type 2 (annular region) and there is 1 segment
*  to add
*  For $\pi/2 < \theta_{f}$ the segment must be of
*  type 1 (annular region) and there are 2 segments
*  to add
*----
              IF(THEND .LT. -PIO2) THEN
                IF(TYPES(1,IPOINT) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be an annular region')
                NADSEG=0
              ELSE IF(THEND .LT. PIO2) THEN
                IF(TYPES(1,IPOINT) .NE. 2) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be a pin region')
                NADSEG=1
                XYADD(1,1)=XYIAPD(1,2)
                XYADD(2,1)=XYIAPD(2,2)
                TYADD(1,1)=1
                TYADD(2,1)=2
              ELSE
                IF(TYPES(1,IPOINT) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be an annular region')
                NADSEG=2
                XYADD(1,1)=XYIAPD(1,2)
                XYADD(2,1)=XYIAPD(2,2)
                TYADD(1,1)=1
                TYADD(2,1)=2
                XYADD(1,2)=XYIAPD(1,1)
                XYADD(2,2)=XYIAPD(2,1)
                TYADD(1,2)=2
                TYADD(2,2)=1
              ENDIF
            ELSE IF(THBEG .LT. PIO2) THEN
*----
*  For $-\pi/2 < \theta_{i}\le \pi/2$ the segment must be of
*  type 2 (pin region)
*----
              IF(TYPES(2,IPOINT-1) .NE. 2) CALL XABORT(NAMSBR//
     >': Error -> Initial line segment must be a pin region')
*----
*  For $-\pi/2 < \theta_{f}\le \pi/2$ the segment must be of
*  type 2 (pin region) and there is no segment
*  to add
*  For $\pi/2 < \theta_{f}\le 3\pi/2$ the segment must be of
*  type 1 (annular region) and there is 1 segment
*  to add
*  For $3\pi/2< \theta_{f}$ the segment must be of
*  type 2 (pin region) and there are 2 segments
*  to add
*----
              IF(THEND .LT. PIO2) THEN
                IF(TYPES(1,IPOINT) .NE. 2) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be a pin region')
                NADSEG=0
              ELSE IF(THEND .LT. TPIO2) THEN
                IF(TYPES(1,IPOINT) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be an annular region')
                NADSEG=1
                XYADD(1,1)=XYIAPD(1,1)
                XYADD(2,1)=XYIAPD(2,1)
                TYADD(1,1)=2
                TYADD(2,1)=1
              ELSE
                IF(TYPES(1,IPOINT) .NE. 2) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be a pin region')
                NADSEG=2
                XYADD(1,1)=XYIAPD(1,1)
                XYADD(2,1)=XYIAPD(2,1)
                TYADD(1,1)=2
                TYADD(2,1)=1
                XYADD(1,2)=XYIAPD(1,2)
                XYADD(2,2)=XYIAPD(2,2)
                TYADD(1,2)=1
                TYADD(2,2)=2
              ENDIF
            ELSE
*----
*  For $\pi/2 < \theta_{i}$ the segment must be of
*  type 1 (annular region)
*----
              IF(TYPES(2,IPOINT-1) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Initial line segment must be an annular region')
*----
*  For $\pi/2 < \theta_{f}\le 3\pi/2$ the segment must be of
*  type 1 (annular region) and there is no segment
*  to add
*  For $3\pi/2< \theta_{f}\le 5*\pi/2$ the segment must be of
*  type 2 (pin region) and there is 1 segment
*  to add
*  For $5*\pi/2 < \theta_{f}$ the segment must be of
*  type 1 (annular region) and there are 2 segments
*  to add
*----
              IF(THEND .LT. TPIO2) THEN
                IF(TYPES(1,IPOINT) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be an annular region')
                NADSEG=0
              ELSE IF(THEND .LT. FPIO2) THEN
                IF(TYPES(1,IPOINT) .NE. 2) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be a pin region')
                NADSEG=1
                XYADD(1,1)=XYIAPD(1,2)
                XYADD(2,1)=XYIAPD(2,2)
                TYADD(1,1)=1
                TYADD(2,1)=2
              ELSE
                IF(TYPES(1,IPOINT) .NE. 1) CALL XABORT(NAMSBR//
     >': Error -> Final line segment must be an annular region')
                NADSEG=2
                XYADD(1,1)=XYIAPD(1,2)
                XYADD(2,1)=XYIAPD(2,2)
                TYADD(1,1)=1
                TYADD(2,1)=2
                XYADD(1,2)=XYIAPD(1,1)
                XYADD(2,2)=XYIAPD(2,1)
                TYADD(1,2)=2
                TYADD(2,2)=1
              ENDIF
            ENDIF
          ENDIF
*----
*  Move end segments to create place for new segments
*----
          IF(NADSEG .GT. 0) THEN
            DO ISEG=NBPTS,IPOINT,-1
              POINTS(1,ISEG+NADSEG)=POINTS(1,ISEG)
              POINTS(2,ISEG+NADSEG)=POINTS(2,ISEG)
              TYPES(1,ISEG+NADSEG)=TYPES(1,ISEG)
              TYPES(2,ISEG+NADSEG)=TYPES(2,ISEG)
            ENDDO
*----
*  Insert new segments
*----
            DO ISEG=NADSEG,1,-1
              POINTS(1,IPOINT+ISEG-1)=XYADD(1,ISEG)
              POINTS(2,IPOINT+ISEG-1)=XYADD(2,ISEG)
              TYPES(1,IPOINT+ISEG-1)=TYADD(1,ISEG)
              TYPES(2,IPOINT+ISEG-1)=TYADD(2,ISEG)
            ENDDO
            NBPTS=NBPTS+NADSEG
          ENDIF
        ENDDO
        NPOINT=NBPTS
      ENDIF
      IF(NPOINT .EQ. 0) THEN
        NXTPRA=0
        VOLINT=DZERO
      ELSE
        VOLINT=DZERO
        IF(IPRINT .GE. 200) THEN
*----
*  Print cell description if required
*----
          WRITE(IOUT,6010)
          DO IPOINT=1,NPOINT
            IF(IPOINT .EQ. NPOINT) THEN
              WRITE(IOUT,6011) POINTS(1,IPOINT),POINTS(2,IPOINT),
     >                         TYPES(1,IPOINT),TYPES(2,IPOINT)
            ELSE
              WRITE(IOUT,6012) POINTS(1,IPOINT),POINTS(2,IPOINT),
     >                         TYPES(1,IPOINT),TYPES(2,IPOINT)
            ENDIF
          ENDDO
        ENDIF
        DO IPOINT=1,NPOINT-1
*----
*  Add contribution under line segments
*----
          DVOL=(POINTS(1,IPOINT)-POINTS(1,IPOINT+1))
     >        *(POINTS(2,IPOINT)+POINTS(2,IPOINT+1))/DTWO
          VOLINT=VOLINT+DVOL
          IF(TYPES(2,IPOINT) .EQ. 1) THEN
*----
*  Add annular region contribution (annular region is not centered)
*  1- Find angular width for two points
*  2- Compute volume above line joining the two points
*----
            XYBEG(1)=POINTS(1,IPOINT)-POSANN(1)
            XYBEG(2)=POINTS(2,IPOINT)-POSANN(2)
            XYEND(1)=POINTS(1,IPOINT+1)-POSANN(1)
            XYEND(2)=POINTS(2,IPOINT+1)-POSANN(2)
*----
*  Find angular location of points
*----
            ACARG=XYBEG(1)/POSANN(0)
            IF(ACARG .GE. 1.0D0) THEN
              THBEG=ACOS(1.0D0)
            ELSE IF(ACARG .LE. -1.0D0) THEN
              THBEG=ACOS(-1.0D0)
            ELSE
              THBEG=ACOS(ACARG)
            ENDIF
            IF(XYBEG(2) .LT. DZERO) THBEG=-THBEG
            ACARG=XYEND(1)/POSANN(0)
            IF(ACARG .GE. 1.0D0) THEN
              THEND=ACOS(1.0D0)
            ELSE IF(ACARG .LE. -1.0D0) THEN
              THEND=ACOS(-1.0D0)
            ELSE
              THEND=ACOS(ACARG)
            ENDIF
            IF(XYEND(2) .LT. DZERO) THEND=-THEND
            IF(THEND .LT. THBEG) THEND=DTWO*PI+THEND
            HDT=(THEND-THBEG)/DTWO
            FACT=COS(HDT)*SIN(HDT)
            DVOL=RADAN2*(HDT-FACT)
            VOLINT=VOLINT+DVOL
          ELSE IF (TYPES(2,IPOINT) .EQ. 2) THEN
*----
*  Add pin region contribution (pin is centered)
*  1- Find angular width for the two points
*  2- Compute volume above line joining the two points
*----
            XYBEG(1)=POINTS(1,IPOINT)-POSPIN(1)
            XYBEG(2)=POINTS(2,IPOINT)-POSPIN(2)
            XYEND(1)=POINTS(1,IPOINT+1)-POSPIN(1)
            XYEND(2)=POINTS(2,IPOINT+1)-POSPIN(2)
*----
*  Find angular location of points
*----
            ACARG=XYBEG(1)/POSPIN(0)
            IF(ACARG .GE. 1.0D0) THEN
              THBEG=ACOS(1.0D0)
            ELSE IF(ACARG .LE. -1.0D0) THEN
              THBEG=ACOS(-1.0D0)
            ELSE
              THBEG=ACOS(ACARG)
            ENDIF
            IF(XYBEG(2) .LT. DZERO) THBEG=-THBEG
            ACARG=XYEND(1)/POSPIN(0)
            IF(ACARG .GE. 1.0D0) THEN
              THEND=ACOS(1.0D0)
            ELSE IF(ACARG .LE. -1.0D0) THEN
              THEND=ACOS(-1.0D0)
            ELSE
              THEND=ACOS(ACARG)
            ENDIF
            IF(XYEND(2) .LT. DZERO) THEND=-THEND
            IF(THEND .LT. THBEG) THEND=DTWO*PI+THEND
            HDT=(THEND-THBEG)/DTWO
            FACT=COS(HDT)*SIN(HDT)
            DVOL=RADPI2*(HDT-FACT)
            VOLINT=VOLINT+DVOL
          ENDIF
        ENDDO
      ENDIF
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6020) VOLINT
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('FinalSegments={')
 6011 FORMAT('{',F20.10,',',F20.10,',',I10,',',I10,'}};')
 6012 FORMAT('{',F20.10,',',F20.10,',',I10,',',I10,'},')
 6015 FORMAT('OriginalSegments={')
 6020 FORMAT('Volint=',F20.10,';')
 6100 FORMAT('CartesianRegion={')
 6101 FORMAT(('{',F15.10,',',F15.10,'}',','/))
 6102 FORMAT('{',F15.10,',',F15.10,'}','};')
 6103 FORMAT('RADAN = ',F15.10,';'/
     >       'POSANN={',F15.10,',',F15.10,'};')
 6104 FORMAT('RADIUS= ',F15.10,';'/
     >       'xypin={',F15.10,',',F15.10,'};')
      END
