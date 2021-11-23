*DECK MCTCCC
      SUBROUTINE MCTCCC(NDIM,ITRN,CELLPO,ODIR,POS,ODIRC,POSC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Change global coordinates to turned cell coordinates 
* (adapted from NXTRTL.f).
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* NDIM    dimensions of problem.
* ITRN    geometry original turn number.
* CELLPO  cell global coordinates.
* ODIR    search (octant) direction in global geometry.
* POS     global coordinates.
*
*Parameters: output
* POSC    final coordinates.
* ODIRC   search (octant) direction in cell.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NDIM,ITRN,ODIR(3),ODIRC(3)
      DOUBLE PRECISION CELLPO(3,2),POS(3),POSC(3)
*----
*  LOCAL VARIABLES
*----
      INTEGER IDIR,IKT
      DOUBLE PRECISION POSO(3)
*----
*  CHANGE COORDINATES TO LOCAL CELL COORDINATES (ORIGIN=CELL CENTER)
*----
      DO IDIR=1,NDIM
        POSO(IDIR)=POS(IDIR)
     1            -0.5D0*(CELLPO(IDIR,1)+CELLPO(IDIR,2))
      ENDDO
*----
*  CHANGE COORDINATES ACCORDING TO TURN
*----
*     Z AXIS REFLECTION FOR 3-D GEOMETRY
      IKT=ITRN
      IF(NDIM .EQ. 3) THEN
        IF(ITRN .GT. 12 ) THEN
          IKT=IKT-12
          POSC(NDIM)=-POSO(NDIM)
          ODIRC(NDIM)=-ODIR(NDIM)
        ELSE
          POSC(NDIM)=POSO(NDIM)
          ODIRC(NDIM)=ODIR(NDIM)
        ENDIF
      ENDIF
      IF(IKT .EQ. 1) THEN
*     NO TURN IN X-Y PLANE
        DO IDIR=1,2
          POSC(IDIR)=POSO(IDIR)
          ODIRC(IDIR)=ODIR(IDIR)
        ENDDO
      ELSE IF(IKT .EQ. 2) THEN 
*     ROTATION OF -PI/2 OF GEOMETRY IMPLIES A ROTATION
*     OF PI/2 OF LINE. 
        POSC(1)=-POSO(2)
        POSC(2)= POSO(1)
        ODIRC(1)=-ODIR(2)
        ODIRC(2)= ODIR(1)
      ELSE IF(IKT .EQ. 3) THEN
*     ROTATION OF PI OF GEOMETRY IMPLIES A ROTATION
*     OF -PI OF LINE. 
        POSC(1)=-POSO(1)
        POSC(2)=-POSO(2)
        ODIRC(1)=-ODIR(1)
        ODIRC(2)=-ODIR(2)
      ELSE IF(IKT .EQ. 4) THEN
*     ROTATION OF -3*PI/2 OF GEOMETRY IMPLIES A ROTATION
*     OF 3PI/2 OF LINE. 
        POSC(1)= POSO(2)
        POSC(2)=-POSO(1)
        ODIRC(1)= ODIR(2)
        ODIRC(2)=-ODIR(1)
      ELSE IF(IKT .EQ. 5) THEN
*     REFLECTION WITH RESPECT TO AXIS  // TO Y
        POSC(1)=-POSO(1)
        POSC(2)= POSO(2)
        ODIRC(1)=-ODIR(1)
        ODIRC(2)= ODIR(2)
      ELSE IF(IKT .EQ. 6) THEN
*     ROTATION OF PI/2 FOLLOWED BY 
*     REFLECTION WITH RESPECT TO AXIS  // TO Y
*     IMPLIES REFLECTION WITH RESPECT TO AXIS  // TO Y
*     FOLLOWED BY A ROTATION OF -PI/2 OF LINE. 
        POSC(1)= POSO(2)
        POSC(2)= POSO(1)
        ODIRC(1)= ODIR(2)
        ODIRC(2)= ODIR(1)
      ELSE IF(IKT .EQ. 7) THEN
*     REFLECTION WITH RESPECT TO AXIS // TO X
        POSC(1)= POSO(1)
        POSC(2)=-POSO(2)
        ODIRC(1)= ODIR(1)
        ODIRC(2)=-ODIR(2)
      ELSE IF(IKT .EQ. 8) THEN
*     ROTATION OF PI/2 FOLLOWED BY
*     REFLECTION WITH RESPECT TO AXIS // TO X
*     IMPLIES REFLECTION WITH RESPECT TO AXIS  // TO X
*     FOLLOWED BY A ROTATION OF -PI/2 OF LINE.
        POSC(1)=-POSO(2)
        POSC(2)=-POSO(1)
        ODIRC(1)=-ODIR(2)
        ODIRC(2)=-ODIR(1)
      ENDIF
*
      RETURN
      END
