*DECK NXTSQD
      SUBROUTINE NXTSQD(IFTRK ,IPRINT,NDIM  ,NQUAD ,NBANGL,
     >                  DANGLT,DDENWT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To transform double precision to simple precision
* quadrature and save on IFTRK.
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
* IFTRK   pointer to the TRACKING file in creation mode.
* IPRINT  print level.
* NDIM    number of dimensions for geometry.
* NQUAD   number of quadrant (in 3-D) and quarter (in 2-D).
* NBANGL  number of angles.
* DANGLT  angles (double precision).
* DDENWT  angular density for each angle (double precision).
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
      INTEGER          IFTRK,IPRINT
      INTEGER          NDIM,NQUAD,NBANGL
      DOUBLE PRECISION DANGLT(NDIM,NQUAD,NBANGL),DDENWT(NQUAD,NBANGL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTSQD')
*----
*  Local variables
*----
      INTEGER          II,IJ,IK,JJ
*----
*  Allocatable arrays
*----
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ANGLT,DENWT
*----
*  Scratch storage allocation
*   ANGLT   angles.
*   DENWT   angular density for each angle.
*----
      ALLOCATE(ANGLT(NDIM*NQUAD*NBANGL),DENWT(NQUAD*NBANGL))
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      JJ=0
      DO IK=1,NBANGL
        DO IJ=1,NQUAD
          DO II=1,NDIM
            JJ=JJ+1
            ANGLT(JJ)=DANGLT(II,IJ,IK)
          ENDDO
        ENDDO
      ENDDO
      JJ=0
      DO IK=1,NBANGL
        DO IJ=1,NQUAD
          JJ=JJ+1
          DENWT(JJ)=DDENWT(IJ,IK)
        ENDDO
      ENDDO
      WRITE(IFTRK) (ANGLT(JJ),JJ=1,NQUAD*NBANGL*NDIM)
      WRITE(IFTRK) (DENWT(JJ),JJ=1,NQUAD*NBANGL)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(DENWT,ANGLT)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
