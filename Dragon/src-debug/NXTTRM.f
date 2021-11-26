*DECK NXTTRM
      SUBROUTINE NXTTRM(ICTRN ,INTRN ,DRW   ,DNW   )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To determine the final mesh of a cell after turn.
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
* ICTRN   turn of geometry.
* INTRN   turn of cell.
* DRW     mesh of geometry before turn.
*
*Parameters: output
* DNW     mesh of cell after turns.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  This routine is based on the LELCSY routine written by
*  R. Roy and G. Marleau for the EXCELT: module.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          ICTRN,INTRN
      DOUBLE PRECISION DRW(3),DNW(3)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTRM')
*----
*  Local variables
*----
      INTEGER          IDIR,IKT,IRXY
      DOUBLE PRECISION DTW(3)
*----
*  1) turn geometry
*----
      DO IDIR=1,3
        DTW(IDIR)=DRW(IDIR)
      ENDDO
      IKT=MOD(ICTRN-1,12)+1
      IRXY=MOD(IKT,2)
      IF(IRXY .EQ. 0) THEN
*----
*  These rotations inply interchange of $X$ and $Y$
*----
        DTW(2)=DRW(1)
        DTW(1)=DRW(2)
      ENDIF
      DO IDIR=1,3
        DNW(IDIR)=DTW(IDIR)
      ENDDO
*----
*  2) turn cell
*----
      IKT=MOD(INTRN-1,12)+1
      IRXY=MOD(IKT,2)
      IF(IRXY .EQ. 0) THEN
*----
*  These rotations inply interchange of $X$ and $Y$
*----
        DTW(2)=DNW(1)
        DTW(1)=DNW(2)
      ENDIF
      DO IDIR=1,3
        DNW(IDIR)=DTW(IDIR)
      ENDDO
*----
*  Processing finished:
*  and return
*----
      RETURN
      END
