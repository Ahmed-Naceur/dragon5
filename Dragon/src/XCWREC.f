*DECK XCWREC
      SUBROUTINE XCWREC(ANGD,SIDE,TRKPOS,LINTER,ROTPOS,INDS,IMS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Track outer rectangle for 2-D cluster.
*
*Copyright:
* Copyright (C) 1992 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G.Marleau
*
*Parameters: input
* ANGD    track director cosines (cos(a),sin(a)).
* SIDE    side of rectangle.
* IMS     surface merge.
*
*Parameters: input/output
* TRKPOS  one track point at input  (*,1).
*         Track origin at output    (*,1).
*         Track origin at input     (*,2).
*
*Parameters: output
* LINTER  intersection logical.
* ROTPOS  position wrt rotated axis.
* INDS    surface of intersection.
*
*----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION ANGD(2),SIDE(2),TRKPOS(2,2),ROTPOS(2,2)
      INTEGER    IMS(6),INDS(2)
      LOGICAL    LINTER
*----
*  EQUATIONS FOR SIDES
*     SIDE 1: XR= SIDE(1)/2   (-SIDE(2)/2<=YR<=SIDE(2)/2)
*     SIDE 2: YR= SIDE(2)/2   (-SIDE(1)/2<=XR<=SIDE(1)/2)
*     SIDE 3: XR=-SIDE(1)/2   (-SIDE(2)/2<=YR<=SIDE(2)/2)
*     SIDE 4: YR=-SIDE(2)/2   (-SIDE(1)/2<=XR<=SIDE(1)/2)
*  TRACK EQUATION
*             YR=  TAN(ANGD)*(XR-TRKPOS(1,1))+TRKPOS(2,1)
*          OR XR=COTAN(ANGD)*(YR-TRKPOS(2,1))+TRKPOS(1,1)
*----
      YTOP=0.5D0*SIDE(2)
      XTOP=0.5D0*SIDE(1)
      LINTER=.FALSE.
      IF(ANGD(1).EQ.0.0D0) THEN
*----
*  TRACK PARALLEL TO Y
*  TRACK INTERCEPT SURFACE 4 AND 2
*----
        IF(ABS(TRKPOS(1,1)).LT.XTOP) THEN
          TRKPOS(1,2)=TRKPOS(1,1)
           IF(ANGD(2).LT.0.0) THEN
            INDS(2)=IMS(4)
            INDS(1)=IMS(2)
            TRKPOS(2,2)=-YTOP
            TRKPOS(2,1)=YTOP
          ELSE
            INDS(2)=IMS(2)
            INDS(1)=IMS(4)
            TRKPOS(2,2)=YTOP
            TRKPOS(2,1)=-YTOP
          ENDIF
          LINTER=.TRUE.
        ENDIF
      ELSE IF(ANGD(2).EQ.0.0D0) THEN
*----
*  TRACK PARALLEL TO X
*  TRACK INTERCEPT SURFACE 3 AND 1
*----
        IF(ABS(TRKPOS(2,1)).LT.YTOP) THEN
          TRKPOS(2,2)=TRKPOS(2,1)
          IF(ANGD(1).LT.0.0D0) THEN
            INDS(2)=IMS(3)
            INDS(1)=IMS(1)
            TRKPOS(1,2)=-XTOP
            TRKPOS(1,1)=XTOP
          ELSE
            INDS(2)=IMS(1)
            INDS(1)=IMS(3)
            TRKPOS(1,2)=XTOP
            TRKPOS(1,1)=-XTOP
          ENDIF
          LINTER=.TRUE.
        ENDIF
      ELSE
        NSEG=1
        COSAI=1.0/ANGD(1)
        SINAI=1.0/ANGD(2)
*----
*    SLOPEY=TAN(ANGD)
*    SLOPEX=COTAN(ANGD)
*    RINTY=TRKPOS(2,1)-SLOPEY*TRKPOS(1,1)
*    RINTX=TRKPOS(1,1)-SOLPEX*TRKPOS(2,1)
*----
        SLOPEY=ANGD(2)*COSAI
        SLOPEX=ANGD(1)*SINAI
        RINTY=TRKPOS(2,1)-SLOPEY*TRKPOS(1,1)
        RINTX=TRKPOS(1,1)-SLOPEX*TRKPOS(2,1)
*----
*  SURFACE 3: YR=RINTY-SLOPEY*XTOP
*            (-YTOP <=YR<= YTOP)
*----
        TRKPOS(2,NSEG)=RINTY-SLOPEY*XTOP
        IF( ABS(TRKPOS(2,NSEG)).LE.YTOP ) THEN
*----
*  TRACK INTERSEPT SURFACE 3
*----
          INDS(NSEG)=IMS(3)
          TRKPOS(1,NSEG)=-XTOP
          NSEG=NSEG+1
        ENDIF
*----
*  SURFACE 1: YR=RINTY+SLOPEY*XTOP
*            (-YTOP <=YR<= YTOP)
*----
        TRKPOS(2,NSEG)=RINTY+SLOPEY*XTOP
        IF( ABS(TRKPOS(2,NSEG)).LE.YTOP ) THEN
*----
*  TRACK INTERSEPT SURFACE 1
*----
          INDS(NSEG)=IMS(1)
          TRKPOS(1,NSEG)=XTOP
          IF(NSEG.EQ.2) GO TO 100
          NSEG=NSEG+1
        ENDIF
*----
*  SURFACE 4: XR=RINTX-SLOPEX*YTOP
*            (-XTOP <=XR<= XTOP)
*----
        TRKPOS(1,NSEG)=RINTX-SLOPEX*YTOP
        IF( ABS(TRKPOS(1,NSEG)).LE.XTOP ) THEN
*----
*  TRACK INTERSEPT SURFACE 4
*----
          INDS(NSEG)=IMS(4)
          TRKPOS(2,NSEG)=-YTOP
          IF(NSEG.EQ.2) GO TO 100
          NSEG=NSEG+1
        ENDIF
*----
*  SURFACE 2: XR=RINTX+SLOPEX*YTOP
*            (-XTOP <=XR<= XTOP)
*----
        TRKPOS(1,NSEG)=RINTX+SLOPEX*YTOP
        IF( ABS(TRKPOS(1,NSEG)).LE.XTOP ) THEN
*----
*  TRACK INTERSEPT SURFACE 2
*----
          INDS(NSEG)=IMS(2)
          TRKPOS(2,NSEG)=YTOP
          IF(NSEG.EQ.2) GO TO 100
          NSEG=NSEG+1
        ENDIF
 100    CONTINUE
        IF(NSEG.EQ.2) THEN
          LINTER=.TRUE.
*----
*  REORDER INTERSECTION POINTS FOR DIRECTION
*----
          IF(ANGD(1).LT.0.0D0) THEN
            IF(TRKPOS(1,1).GT.TRKPOS(1,2)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ELSE
            IF(TRKPOS(1,1).GT.TRKPOS(1,2)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ENDIF
          IF(ANGD(2).LT.0.0D0) THEN
            IF(TRKPOS(2,2).GT.TRKPOS(2,1)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ELSE
            IF(TRKPOS(2,1).GT.TRKPOS(2,2)) THEN
              TRKTMP=TRKPOS(1,2)
              TRKPOS(1,2)=TRKPOS(1,1)
              TRKPOS(1,1)=TRKTMP
              TRKTMP=TRKPOS(2,2)
              TRKPOS(2,2)=TRKPOS(2,1)
              TRKPOS(2,1)=TRKTMP
              INDT=INDS(2)
              INDS(2)=INDS(1)
              INDS(1)=INDT
            ENDIF
          ENDIF
        ENDIF
      ENDIF
*----
*  ROTATE RECTANGLE BY ANGD
*----
      IF(LINTER) THEN
        DO 110 II=1,2
          ROTPOS(1,II)=ANGD(1)*TRKPOS(1,II)+ANGD(2)*TRKPOS(2,II)
          ROTPOS(2,II)=-ANGD(2)*TRKPOS(1,II)+ANGD(1)*TRKPOS(2,II)
 110    CONTINUE
      ENDIF
      RETURN
      END
