*DECK XCWHEX
      SUBROUTINE XCWHEX(ANGD,RADC,SIDE,LINTER,XPOS,INDS,IMS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Track outer hexagone for 2-D cluster geometry.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G.Marleau
*
*Parameters: input
* ANGD    track angle.
* RADC    Y-position of track > 0.
* SIDE    side of hexagone.
* IMS     surface merge.
*
*Parameters: output
* LINTER  intersection logical.
* XPOS    points of intersection.
* INDS    surface of intersection.
*
*Comments:
*  Equations for sides of hexagone YR(XR)
*     SIDE 1: YR=-SQ3*XR+SQ3*SIDE   (          0 <=YR<= SQ3*SIDE/2)
*                                   (     SIDE/2 <=XR<= SIDE      )
*     SIDE 2: YR= SQ3*SIDE/2        (    -SIDE/2 <=XR<=  SIDE/2   )
*     SIDE 3: YR= SQ3*XR+SQ3*SIDE   (          0 <=YR<= SQ3*SIDE/2)
*                                   (      -SIDE <=XR<= -SIDE/2   )
*     SIDE 4: YR=-SQ3*XR-SQ3*SIDE   (-SQ3*SIDE/2 <=YR<= 0         )
*                                   (      -SIDE <=XR<= -SIDE/2   )
*     SIDE 5: YR=-SQ3*SIDE/2        (    -SIDE/2 <=XR<= SIDE/2    )
*     SIDE 6: YR= SQ3*XR-SQ3*SIDE   (-SQ3*SIDE/2 <=YR<= 0         )
*                                  (     SIDE/2 <=XR<= SIDE      )
*  Equations for sides of hexagone XR(YR)
*     SIDE 1: XR=-OSQ3*YR+SIDE      (          0 <=YR<= SQ3*SIDE/2)
*                                   (     SIDE/2 <=XR<= SIDE      )
*     SIDE 3: XR= OSQ3*YR-SIDE      (          0 <=YR<= SQ3*SIDE/2)
*                                   (      -SIDE <=XR<= -SIDE/2   )
*     SIDE 4: XR=-OSQ3*YR-SIDE      (-SQ3*SIDE/2 <=YR<= 0         )
*                                   (      -SIDE <=XR<= -SIDE/2   )
*     SIDE 6: XR= OSQ3*YR+SIDE      (-SQ3*SIDE/2 <=YR<= 0         )
*                                   (     SIDE/2 <=XR<= SIDE      )
*  TRACK EQUATION:
*             YR= SQ3*(SLOPEY*XR+RINTY)
*          OR XR= OSQ3*SLOPEX*YR-RINTX
*
*----------------------------------------------------------------------
*
      PARAMETER (SQ3=1.73205080756887729,OSQ3=0.577350269189625795)
      INTEGER    IMS(6),INDS(2)
      LOGICAL    LINTER
      REAL       ANGD,RADC,SIDE,XPOS(2)
*----
      YRINT=SQ3*SIDE
      YLIM=0.5*YRINT
      XLIM=0.5*SIDE
      SINA=SIN(ANGD)
      COSA=COS(ANGD)
      LINTER=.FALSE.
      IF(COSA.EQ.0.0) THEN
*----
*  TRACK PARALLEL TO Y
*----
        IF( RADC.LT. XLIM ) THEN
*----
*  TRACK INTERCEPT SURFACE 5 AND 2
*----
          IF(SINA.LT.0.0) THEN
            INDS(2)=IMS(5)
            INDS(1)=IMS(2)
          ELSE
            INDS(2)=IMS(2)
            INDS(1)=IMS(5)
          ENDIF
          XPOS(2)=YLIM
          XPOS(1)=-XPOS(2)
          LINTER=.TRUE.
        ELSE IF(RADC.LE.SIDE) THEN
*----
*  TRACK INTERCEPT SURFACE 3 AND 4 OR 6 AND 1
*----
          IF(SINA.LT.0.0) THEN
            INDS(2)=IMS(3)
            INDS(1)=IMS(4)
          ELSE
            INDS(2)=IMS(6)
            INDS(1)=IMS(1)
          ENDIF
          XPOS(2)=YRINT-SQ3*RADC
          XPOS(1)=-XPOS(2)
          LINTER=.TRUE.
        ENDIF
      ELSE IF(SINA.EQ.0.0) THEN
*----
*  TRACK PARALLEL TO X
*----
        IF(RADC.LE.YLIM ) THEN
*----
*  TRACK INTERCEPT SURFACE 6 AND 4
*----
          INDS(2)=IMS(4)
          INDS(1)=IMS(6)
          XPOS(2)= OSQ3*RADC+SIDE
          XPOS(1)=-XPOS(2)
          LINTER=.TRUE.
        ENDIF
      ELSE
        NSEG=0
        COSAI=1.0/COSA
        SINAI=1.0/SINA
        SLOPEY=OSQ3*SINA*COSAI
        SLOPEX=SQ3*COSA*SINAI
        RINTY=OSQ3*RADC*COSAI
        RINTX=RADC*SINAI
        XREF=RADC*COSAI*SINA
        OPSY=1.0/(1+SLOPEY)
        OMSY=1.0/(1-SLOPEY)
        XLSX=SLOPEX*XLIM
        SPRY=SIDE+RINTY
        SMRY=SIDE-RINTY
*----
*  SURFACE 1: XR=(SIDE-RINTY)/(1+SLOPEY)
*            (SIDE/2 <=XR<= SIDE)
*----
        XR=SMRY*OPSY
        IF( (XLIM.LE.XR) .AND. (XR.LE.SIDE) ) THEN
*----
*  TRACK INTERSEPT SURFACE 1
*----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(1)
          XPOS(NSEG)=XR
        ENDIF
*----
*  SURFACE 2: XR= SLOPEX*SIDE/2-RINTX
*            (-SIDE/2 <=XR<= SIDE/2)
*----
        XR=XLSX-RINTX
        IF( ABS(XR).LE.XLIM ) THEN
*----
*  TRACK INTERSEPT SURFACE 2
*----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(2)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
*----
*  SURFACE 3: XR=-(SIDE-RINTY)/(1-SLOPEY)
*            (-SIDE <=XR<= -SIDE/2)
*----
        XR=-SMRY*OMSY
        IF( (-SIDE.LE.XR) .AND. (XR.LE.-XLIM) )THEN
*----
*  TRACK INTERSEPT SURFACE 3
*----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(3)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
*----
*  SURFACE 4: XR=-(SIDE+RINTY)/(1+SLOPEY)
*            (-SIDE <=XR<= -SIDE/2)
*----
        XR=-SPRY*OPSY
        IF( (-SIDE.LE.XR) .AND. (XR.LE.-XLIM) ) THEN
*----
*  TRACK INTERSEPT SURFACE 4
*----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(4)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
*----
*  SURFACE 5: XR=-SLOPEX*SIDE/2-RINTX
*            (-SIDE/2 <=XR<= SIDE/2)
*----
        XR=-XLSX-RINTX
        IF( ABS(XR).LE.XLIM ) THEN
*----
*  TRACK INTERSEPT SURFACE 5
*----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(5)
          XPOS(NSEG)=XR
          IF(NSEG.EQ.2) GO TO 100
        ENDIF
*----
*  SURFACE 6: XR=(RINTY+SIDE)/(1-SLOPEY)
*            (SIDE/2 <=XR<= SIDE)
*----
        XR=SPRY*OMSY
        IF( (XLIM.LE.XR) .AND. (XR.LE.SIDE) ) THEN
*----
*  TRACK INTERSEPT SURFACE 6
*----
          NSEG=NSEG+1
          INDS(NSEG)=IMS(6)
          XPOS(NSEG)=XR
        ENDIF
 100    CONTINUE
        IF(NSEG.EQ.2) THEN
          LINTER=.TRUE.
*----
*  ROTATE HEXAGONE BY -ANGD
*----
          XPOS(1)=XREF+XPOS(1)*COSAI
          XPOS(2)=XREF+XPOS(2)*COSAI
          IF( XPOS(1).GT.XPOS(2) ) THEN
            INDT=INDS(2)
            INDS(2)=INDS(1)
            INDS(1)=INDT
            XPOST=XPOS(2)
            XPOS(2)=XPOS(1)
            XPOS(1)=XPOST
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
