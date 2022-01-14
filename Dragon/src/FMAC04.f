      SUBROUTINE FMAC04(NGPRT,NPART,NGP,IEMAX,IEMIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the index in a FMAC-M array.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NGPRT   number of energy groups per particle.
* NPART   number of particles in FMAC-M file.
* NGP     number of energy groups for all particles.
* IEMAX   index of upper energy information in FMAC-M array.
*
*Parameters: output
* IEMIN   index of lower energy information in FMAC-M array.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGPRT(NPART),NPART,NGP,IEMAX,IEMIN
*
      IEMIN=IEMAX+1
      II=NGP+NPART+1
      NN=0
      DO I=1,NPART
        II=II-1
        NN=NN+NGPRT(I)
        IF(IEMAX.LT.NN) THEN
          EXIT
        ELSE IF(IEMAX.EQ.NN) THEN
          IEMIN=II
          EXIT
        ENDIF
      ENDDO
      RETURN
      END
