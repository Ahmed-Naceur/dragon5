*DECK NXTFID
      SUBROUTINE NXTFID(IX,IY,IZ,NSURC,NREGC,MESHCZ,IMX,IMY,IMR,INDEX,
     1                  IDSUR,IDREG,IDZ,ITYP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Locate all the regions/surfaces within a cell/pin corresponding
* to a certain x, y and r position along the projection axis.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IX      first direction perpendicular to the projection axis.
* IY      second direction perpendicular to the projection axis.
* IZ      projection axis.
* NSURC   number of surfaces for the cells/pins.
* NREGC   number of regions for the cells/pins.
* MESHCZ  number of meshes along the projection axis.
* IMX     x index to locate.
* IMY     y index to locate.
* IMR     r index to locate.
* INDEX   cell/pin index vector.
* IDSUR   surface index array.
* IDREG   region index array.
*
*Parameters: output
* IDZ     regions/surfaces encountered along the projection axis.
* ITYP    flag for "what was encountered?":
*          = 0 non existing IMX,IMY,IMR combination;
*          = 1 non-vanishing top/bottom surfaces and regions;
*          = 2 vanishing region;
*          =-1 non-vanishing lateral surface;
*          =-2 vanishing lateral surface.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IX,IY,IZ,NSURC,NREGC,MESHCZ,IMX,IMY,IMR,
     1 INDEX(5,-NSURC:NREGC),IDSUR(NSURC),IDREG(NREGC),IDZ(0:MESHCZ+1),
     2 ITYP
*----
*  LOCAL VARIABLES
*----
      INTEGER II,IND
*
      ITYP=0
      DO II=-NSURC,-1
*     surfaces
         IF ((INDEX(IX,II).EQ.IMX).AND.
     1       (INDEX(IY,II).EQ.IMY).AND.
     2       (INDEX( 4,II).EQ.IMR)) THEN
            IND=INDEX(IZ,II)
            IF (IND.EQ.-1) THEN
*           a top surface
               IDZ(0)=-ABS(IDSUR(-II))
               IF (IDZ(0).NE.0) ITYP=1
            ELSEIF(IND.EQ.-2) THEN
*           a bottom surface
               IDZ(MESHCZ+1)=-ABS(IDSUR(-II))
               IF (IDZ(MESHCZ+1).NE.0) ITYP=1
            ELSE
*           a lateral surface
               IDZ(IND)=-ABS(IDSUR(-II))
               ITYP=-2
               IF (IDZ(IND).NE.0) ITYP=-1
            ENDIF
         ENDIF
      ENDDO
      DO II=1,NREGC
*     regions
         IF ((INDEX(IX,II).EQ.IMX).AND.
     1       (INDEX(IY,II).EQ.IMY).AND.
     2       (INDEX( 4,II).EQ.IMR)) THEN
            IND=INDEX(IZ,II)
            IDZ(IND)=ABS(IDREG(II))
            ITYP=2
            IF (IDZ(IND).NE.0) ITYP=1
         ENDIF
      ENDDO
*
      RETURN
      END
