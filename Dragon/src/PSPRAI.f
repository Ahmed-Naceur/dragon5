*DECK PSPRAI
      SUBROUTINE PSPRAI(MXSEG ,NPTS  ,XYPOS ,CENTER,RCIRC ,
     >                  NSEG  ,IORDER,RADANG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find general closed Cartesian region intersections
* with annular region and order points for plotting.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* MXSEG   maximum number of segments.
* NPTS    number of corners.
* XYPOS   X and Y position of corners.
* CENTER  X and Y position of annulus center.
* RCIRC   annulus radius.
*
*Parameters: output
* NSEG    number of region intersection.
*         number of segments is NSEG-1
* IORDER  type of region:
*         = -2 arc segment begins;
*         = -1 arc segment ends;
*         =  0 close path;
*         >  0 corner.
* RADANG  segments intersection points
*         with respect to annular region center:
*         RADANG(1) = radial position;
*         RADANG(2) = angular position.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      REAL             PI
      PARAMETER       (IOUT=6,PI=3.1415926535897932,NAMSBR='PSPRAI')
*----
*  ROUTINE PARAMETERS
*----
      INTEGER          MXSEG,NPTS,NSEG
      INTEGER          IORDER(MXSEG)
      REAL             XYPOS(2,NPTS),CENTER(2),RCIRC,
     >                 RADANG(2,MXSEG)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IPT,ICUR,INXT
      REAL             XYCUR(2),XYNXT(2),RADCUR,RADNXT,
     >                 XANN,XANNI,XANNF,XANNT,YANNI,YANNF,YANNT,
     >                 DELX,DELY,DELL,RCIRCM
*----
*  SCAN OVER CORNERS AND SIDES
*----
*      write(6,*) ' Circle ',MXSEG ,NPTS  ,CENTER(1),CENTER(2),RCIRC
      NSEG=1
      IORDER(NSEG)=0
      RCIRCM=0.0
      DO IPT=1,NPTS
        ICUR=IPT
        INXT=MOD(IPT,NPTS)+1
        XYCUR(1)=XYPOS(1,ICUR)-CENTER(1)
        XYCUR(2)=XYPOS(2,ICUR)-CENTER(2)
        XYNXT(1)=XYPOS(1,INXT)-CENTER(1)
        XYNXT(2)=XYPOS(2,INXT)-CENTER(2)
        RADCUR=SQRT(XYCUR(1)*XYCUR(1)+XYCUR(2)*XYCUR(2))
        RADNXT=SQRT(XYNXT(1)*XYNXT(1)+XYNXT(2)*XYNXT(2))
*        write(6,*) ' Point ',IPT ,ICUR  ,XYCUR(1),XYCUR(2),
*     >  XYNXT(1),XYNXT(2),RADCUR,RADNXT,RADCUR-RCIRC
        RCIRCM=RCIRC
        IF(RADCUR .EQ. RCIRC) RCIRCM=RCIRCM+0.00001
*----
*        WRITE(6,6002) IPT,XYCUR,XYNXT
* 6002 FORMAT('Line ',I5,5X,'Starts =',2F20.10,5X,'Ends =',2F20.10)
*----
*  CHECK IF CURRENT CORNER IS LOCATED INSIDE
*  ANNULAR REGIONS
*----
        IF(RADCUR .LE. RCIRCM) THEN
          IF(IORDER(NSEG) .NE. ICUR) THEN
*----
*  IT IS LOCATED INSIDE
*  SET IORDER TO IPT TO SPECIFY THIS POINT TO CORRESPOND TO
*  CORNER IPT
*----
            NSEG=NSEG+1
            IORDER(NSEG)=ICUR
            RADANG(1,NSEG)=RADCUR
            IF(RADCUR .EQ. 0.0) THEN
              RADANG(2,NSEG)=0.0
            ELSE
              RADANG(2,NSEG)=ATAN2(XYCUR(2),XYCUR(1))
            ENDIF
          ENDIF
        ENDIF
*----
*  Find line direction
*----
        DELY=XYNXT(2)-XYCUR(2)
        DELX=XYNXT(1)-XYCUR(1)
        DELL=SQRT(DELY*DELY+DELX*DELX)
        DELY=DELY/DELL
        DELX=DELX/DELL
        XANNI=XYCUR(1)*DELX+XYCUR(2)*DELY
        YANNI=-XYCUR(1)*DELY+XYCUR(2)*DELX
        XANNF=XYNXT(1)*DELX+XYNXT(2)*DELY
        YANNF=-XYNXT(1)*DELY+XYNXT(2)*DELX
*----
*        WRITE(6,6003) DELX,DELY,XANNI,YANNI,XANNF,YANNF
* 6003 FORMAT('Rotation ',2F20.10,5X,
*     >       'Starts =',2F20.10,5X,'Ends =',2F20.10)
*----
        IF(YANNI .GE. -RCIRCM .AND. YANNI .LE. RCIRCM) THEN
          XANN=-SQRT(RCIRCM*RCIRCM-YANNI*YANNI)
          IF(XANN .GE. XANNI .AND.
     >       XANN .LE. XANNF) THEN
            NSEG=NSEG+1
            IORDER(NSEG)=-1
            RADANG(1,NSEG)=RCIRCM
            XANNT=XANN*DELX-YANNI*DELY
            YANNT=XANN*DELY+YANNI*DELX
            RADANG(2,NSEG)=ATAN2(YANNT,XANNT)
          ENDIF
          XANN=-XANN
          IF(XANN .GE. XANNI .AND.
     >       XANN .LE. XANNF) THEN
            NSEG=NSEG+1
            IORDER(NSEG)=-2
            RADANG(1,NSEG)=RCIRCM
            XANNT=XANN*DELX-YANNI*DELY
            YANNT=XANN*DELY+YANNI*DELX
            RADANG(2,NSEG)=ATAN2(YANNT,XANNT)
          ENDIF
        ENDIF
*----
*  CHECK IF NEXT CORNER OF THE RECTANGLE IS LOCATED INSIDE
*  ANNULAR REGIONS
*----
        IF(RADNXT .LE. RCIRCM) THEN
*----
*  IT IS LOCATED INSIDE
*  SET IORDER TO IPT TO SPECIFY THIS POINT TO CORRESPOND TO
*  CORNER IPT
*----
          NSEG=NSEG+1
          IORDER(NSEG) =INXT
          RADANG(1,NSEG)=RADNXT
          IF(RADNXT .EQ. 0.0) THEN
            RADANG(2,NSEG)=0.0
          ELSE
            RADANG(2,NSEG)=ATAN2(XYNXT(2),XYNXT(1))
          ENDIF
        ENDIF
      ENDDO
*----
*  STORE LAST SEGMENT ALSO AT FIRST POSITION
*  FOR CYCLIC TRACKING
*----
      IF(NSEG .EQ. 1) THEN
        NSEG=2
        IORDER(1)=-2
        RADANG(1,1)=RCIRCM
        RADANG(2,1)=0.0
        IORDER(2)=-1
        RADANG(1,2)=RCIRCM
        RADANG(2,2)=2.0*PI
      ELSE
        IORDER(1) =IORDER(NSEG)
        RADANG(1,1)=RADANG(1,NSEG)
        RADANG(2,1)=RADANG(2,NSEG)
      ENDIF
*      write(6,'(I5,2F20.10)') (IORDER(IPT),RADANG(1,IPT),RADANG(2,IPT),
*     >            IPT=1,NSEG)
      RETURN
      END
