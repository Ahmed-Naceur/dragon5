*DECK XELCRN
      SUBROUTINE XELCRN(IPRINT,RANN2,NRSPX,NRSPY,SPAT,AREAI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find 2-D surface of intersection between annular region and 
* Cartesian plane.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPRINT  print level (active if >=10).
* RANN2   annular region  radius**2.
* NRSPX   number of mesh in x- direction.
* NRSPY   number of mesh in x- direction.
* SPAT    spatial mesh x-direction:
*         SPAT(1,1)       = lower X - position;
*         SPAT(NRSPX+1,1) = upper X - position;
*         SPAT(1,2)       = lower Y - position;
*         SPAT(NRSPY+1,2) = upper Y - position.
*
*Parameters: output
* AREAI   area of intersection.
*
*--------------------------    XELCRN    -------------------------------
*
      IMPLICIT          NONE
      INTEGER           IPRINT,NRSPX,NRSPY
      DOUBLE PRECISION  RANN2,SPAT(NRSPX+1,NRSPY+1),
     >                  AREAI(NRSPX,NRSPY)
*----
*  INTERNAL PARAMETERS
*----
      INTEGER      IOUT
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NAMSBR='XELCRN')
      DOUBLE PRECISION   PI,DZERO
      PARAMETER   (PI=3.14159265358979323846D0,DZERO=0.0D0)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IRP(2,2),IX,NMX,IY,NMY
      DOUBLE PRECISION XELPSC,XELPSI,XYPOS(2,2),XYPOS2(2,2),
     >                 SPXY(2,2),SIXY(2,2),RANN,SURANN 
*----
*  COMPUTE GENERAL ANNULAR REGION INFORMATIONS
*    RANN   = ANNULAR REGION RADIUS
*    SURANN = ANNULAR SURFACE
*  COMPUTE CARTESIAN PARAMETERS
*    NMX =NRSPX+1
*    NMY =NRSPY+1
*  INITIALIZE AREAI TO 0.0D0
*----
      RANN=SQRT(RANN2)
      SURANN=PI*RANN2
      NMX=NRSPX+1
      NMY=NRSPY+1
      CALL XDDSET(AREAI,NRSPX*NRSPY,DZERO)
      CALL XDISET(IRP,4,0)
*----
*  PRINT INITIAL MESH IF REQUIRED
*-----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000)
        WRITE(IOUT,6002) 'ANNULAR RADIUS  '
        WRITE(IOUT,6003) RANN
        WRITE(IOUT,6002) 'ANNULAR SURFACE '
        WRITE(IOUT,6003) SURANN
        WRITE(IOUT,6002) 'X-DIRECTED MESH '
        WRITE(IOUT,6003) (SPAT(IX,1),IX=1,NRSPX+1)
        WRITE(IOUT,6002) 'Y-DIRECTED MESH '
        WRITE(IOUT,6003) (SPAT(IY,2),IY=1,NRSPY+1)
        WRITE(IOUT,6002) 'X-Y SURFACES    '
        WRITE(IOUT,6003) (( (SPAT(IX+1,1)-SPAT(IX,1))
     >                     *(SPAT(IY+1,2)-SPAT(IY,2)),
     >                     IX=1,NRSPX),IY=1,NRSPY)
      ENDIF
*----
*  CYCLE OVER CARTESIAN Y-DIRECTIONS STARTING FROM THE END
*  AND LOCATE Y-MESH POSITION WITH RESPECT TO ANNULUS CENTER
*----
      SPXY(2,2)=DZERO
      DO 110 IY=NMY,1,-1
        XYPOS(2,1)=SPAT(IY,2)
        XYPOS2(2,1)=XYPOS(2,1)*XYPOS(2,1)
*----
*  FIND IF ANNULUS ABOVE, BELOW OR INTERSECT CURRENT Y-PLANE
*  AND COMPUTE
*    SPXY = ANNULAR SURFACE BELOW CURRENT PLANE
*  IF ANNULUS BELOW CURRENT PLANE (XYPOS(2,1)>= RANN)
*    IRPY(2,1)=-1
*    SPXY(2,1)=SURANN
*  IF ANNULUS ABOVE CURRENT PLANE (XYPOS(2,1)<= -RANN)
*    IRPY(2,1)= 1
*    SPXY(2,1)=0.0
*  IF ANNULUS INTERSECT CURRENT ( -RANN < XYPOS(2,1) < RANN)
*    IRPY(2,1)= 0
*    SPXY=XELPSC(RANN,XYPOS(2,1))
*----
        IF(XYPOS(2,1) .GE. RANN) THEN
          IRP(2,1)=-1
          SPXY(2,1)=SURANN
        ELSE IF(XYPOS(2,1) .LE. -RANN) THEN
          IRP(2,1)=1
          SPXY(2,1)=DZERO
        ELSE
          IRP(2,1)=0
          SPXY(2,1)=XELPSC(RANN,XYPOS(2,1))
        ENDIF
*----
*  FOR LAST PLANE IN Y DIRECTION OR
*  Y-PLANE ABOVE ANNULAR VOLUME
*  GO TO LABEL 111
*----
        IF(IY .EQ. NMY .OR. IRP(2,1) .EQ. -1) GO TO 111
*----
*  CYCLE OVER CARTESIAN X-DIRECTIONS STARTING FROM THE END
*  AND LOCATE X-MESH POSITION WITH RESPECT TO ANN CENTER
*----
        SPXY(1,2)=DZERO
        SIXY(2,1)=DZERO
        SIXY(2,2)=DZERO
        DO 120 IX=NMX,1,-1
          XYPOS(1,1)=SPAT(IX,1)
          XYPOS2(1,1)=XYPOS(1,1)*XYPOS(1,1)
*----
*  FIND IF ANNULUS LEFT, RIGHT OR INTERSECT CURRENT X-PLANE
*  AND COMPUTE
*    SPXY      THE  ANNULAR SURFACE LEFT OF CURRENT PLANE
*    SIXY(1,1) THE INTERSECTION BETWEEN THE PART OF THE ANNULUS
*              THE LEFT OF X-PLANE
*              AND THE PART OF THE ANNULUS AT
*              THE BOTTOM OF CURRENT Y-PLANE
*    SIXY(1,2) THE INTERSECTION BETWEEN THE PART OF THE ANNULUS
*              THE LEFT OF X-PLANE
*              AND THE PART OF THE ANNULUS AT
*              THE TOP OF PREVIOUS Y-PLANE
*  IF ANNULUS TO THE LEFT OF CURRENT PLANE (XYPOS(1,1)>= RANN)
*    IRPY(1,1)=-1
*    SPXY(1,1)=SURANN
*    SIXY(1,1)=SPXY(2,1)
*    SIXY(1,2)=SPXY(2,2)
*  IF ANNULUS TO THE RIGHT OF CURRENT (XYPOS(1,1)<= -RANN)
*    IRPY(1,1)= 1
*    SPXY(1,1)=0.0
*    SIXY(1,1)=0.0
*    SIXY(1,2)=0.0
*  IF ANNULUS INTERSECT CURRENT PLANE ( -RANN < XYPOS(1,1) < RANN)
*    IRPY(1,1)= 0
*    SPXY=XELPSC(RANN,XYPOS(1,1))
*    SIXY(1,1)=GEOPSI(1,RANN2,XYPOS,XYPOS2,SPXY)
*    SIXY(1,2)=GEOPSI(2,RANN2,XYPOS,XYPOS2,SPXY)
*----
          SPXY(1,1)=DZERO
          SIXY(1,1)=DZERO
          SIXY(1,2)=DZERO
          IF(XYPOS(1,1) .GE. RANN) THEN
            IRP(1,1)=-1
            SPXY(1,1)=SURANN
            SIXY(1,1)=SPXY(2,1)
            SIXY(1,2)=SPXY(2,2)
          ELSE IF(XYPOS(1,1) .LE. -RANN) THEN
            IRP(1,1)=1
          ELSE
            IRP(1,1)=0
            SPXY(1,1)=XELPSC(RANN,XYPOS(1,1))
            IF(IRP(2,1) .EQ. 0)
     >        SIXY(1,1)=XELPSI(1,RANN2,XYPOS,XYPOS2,SPXY)
            IF(IRP(2,2) .EQ. 0)
     >        SIXY(1,2)=XELPSI(2,RANN2,XYPOS,XYPOS2,SPXY)
          ENDIF
*----
*  FOR LAST PLANE IN X DIRECTION OR
*  X-PLANE TO THE RIGHT OF ANNULAR VOLUME
*  GO TO LABEL 121
*----
          IF(IX .EQ. NMX .OR. IRP(1,1) .EQ. -1) GO TO 121
*----
*  GET SURFACE INTERSECTION BETWEEN ANNULUS AND CARTESIAN REGION
*  LOCATED BETWEEN X-PLANES (IX-> IX+1) AND Y-PLANES (IX -> IY+1)
*  AND STORE IN AREAI(IX,IY)
*----
            AREAI(IX,IY)=SURANN
     >       -SPXY(1,1)-SPXY(1,2)-SPXY(2,1)-SPXY(2,2)
     >       +SIXY(1,1)+SIXY(2,1)+SIXY(1,2)+SIXY(2,2)
*----
*  WHEN ANNULUS ALL LOCATED TO THE RIGHT OF CURRENT X-PLANE
*  EXIT FROM IX LOOP BY GOING TO LABLE 122
*---
          IF(IRP(1,1) .EQ. 1) GO TO 122
 121      CONTINUE
*----
*  RESET IN LOCATION 2 VALUES COMPUTED WITH LOCATION 1
*  WITH ADEQUATE CHANGE OF SIGN FOR SURFACE DIRECTION
*  NAMELY SURFACES LOCATED ON THE LEFT OF X-PLANE BECOME SURFACES
*  LOCATED ON THE RIGHT OF X-PLANE
*----
          SIXY(2,1)=SPXY(2,1)-SIXY(1,1)
          SIXY(2,2)=SPXY(2,2)-SIXY(1,2)
          XYPOS(1,2)=XYPOS(1,1)
          XYPOS2(1,2)=XYPOS2(1,1)
          IRP(1,2)=-IRP(1,1)
          SPXY(1,2)=SURANN-SPXY(1,1)
 120    CONTINUE
 122    CONTINUE
*----
*  WHEN ANNULUS ALL LOCATED ABOVE CURRENT Y-PLANE
*  EXIT FROM IY LOOP BY GOING TO LABLE 112
*---
        IF(IRP(2,1) .EQ. 1) GO TO 112
 111    CONTINUE
*----
*  RESET IN LOCATION 2 VALUES COMPUTED WITH LOCATION 1
*  WITH ADEQUATE CHANGE OF SIGN FOR SURFACE DIRECTION
*  NAMELY SURFACES LOCATED ON THE BELOW Y-PLANE BECOME SURFACES
*  LOCATED ABOVE Y-PLANE
*----
        XYPOS(2,2)=XYPOS(2,1)
        XYPOS2(2,2)=XYPOS2(2,1)
        IRP(2,2)=-IRP(2,1)
        SPXY(2,2)=SURANN-SPXY(2,1)
 110  CONTINUE
 112  CONTINUE
*----
*  PRINT SURFACE INTERSECTIONS IF REQUIRED
*-----
      IF(IPRINT.GE.10) THEN
        WRITE(IOUT,6002) 'CART-ANN AREA   '
        WRITE(IOUT,6003) ((AREAI(IX,IY),IX=1,NRSPX),IY=1,NRSPY)
        WRITE(IOUT,6001)
      ENDIF
*----
*  RETURN
*----
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(/5X,'------ OUTPUT FROM XELCRN ------ ')
 6001 FORMAT(5X,' -------------------------------- '/)
 6002 FORMAT(5X,A16)
 6003 FORMAT(1P,5E16.6)
      END
