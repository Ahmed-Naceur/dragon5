*DECK NXTHCL
      SUBROUTINE NXTHCL(IPRINT,IR    ,IS    ,ISS   ,
     >                  SIDEH ,XLOC  ,YLOC  )
*
*----------
*
*Purpose:
* Locate spatial position of hexagon in assembly of cells.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* IR      crown number.
* IS      cell sector.
* ISS     cell in sector.
* SIDEH   hexagon width.
*
*Parameters: output
* XLOC    X location of cell center in assembly.
* YLOC    Y location of cell center in assembly.
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
      INTEGER          IPRINT,IR,IS,ISS
      DOUBLE PRECISION SIDEH,XLOC,YLOC
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTHCL')
      DOUBLE PRECISION DZERO,DONE,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,
     >                 DHALF=0.5D0,DSQ3O2=0.86602540378444D0)
*----
*  Local variables
*----
      DOUBLE PRECISION SQ32H,H3O2,XLOCR,YLOCR
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6010) IR,IS,ISS,SIDEH
        ENDIF
      ENDIF
      SQ32H=SIDEH*DSQ3O2
      H3O2=3.0D0*DHALF*SIDEH
      IF(IR.EQ.1) THEN
        XLOC=DZERO
        YLOC=DZERO
      ELSE
        XLOCR=SQ32H*(2.0D0*DBLE(IR-1)-DBLE(ISS-1))
        YLOCR=H3O2*DBLE(ISS-1)
        IF(IS .EQ. 1) THEN
*----
*  No rotation
*----
          XLOC=XLOCR
          YLOC=YLOCR
        ELSE IF(IS .EQ. 2) THEN
*----
*  Rotate by Pi/3
*----
          XLOC=DHALF*XLOCR-DSQ3O2*YLOCR
          YLOC=DSQ3O2*XLOCR+DHALF*YLOCR
        ELSE IF(IS .EQ. 3) THEN
*----
*  Rotate by 2*Pi/3
*----
          XLOC=-DHALF*XLOCR-DSQ3O2*YLOCR
          YLOC=DSQ3O2*XLOCR-DHALF*YLOCR
        ELSE IF(IS .EQ. 4) THEN
*----
*  Rotate by Pi
*----
          XLOC=-XLOCR
          YLOC=-YLOCR
        ELSE IF(IS .EQ. 5) THEN
*----
*  Rotate by 4*Pi/3
*----
          XLOC=-DHALF*XLOCR+DSQ3O2*YLOCR
          YLOC=-DSQ3O2*XLOCR-DHALF*YLOCR
        ELSE IF(IS .EQ. 6) THEN
*----
*  Rotate by 5*Pi/3
*----
          XLOC=DHALF*XLOCR+DSQ3O2*YLOCR
          YLOC=-DSQ3O2*XLOCR+DHALF*YLOCR
        ENDIF
      ENDIF
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 10) THEN
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6011) XLOC,YLOC
        ENDIF
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Location of cell in ',
     >' Crown = ',I8,5X,'Sector =',I8,5X,' Cell =',I8,5X,
     >' SIDE  = ',F20.10)
 6011 FORMAT(' X =',F20.10,10X,'Y =',F20.10)
      END
