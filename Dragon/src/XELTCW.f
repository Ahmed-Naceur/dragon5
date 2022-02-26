*DECK XELTCW
      SUBROUTINE XELTCW(NANGLE,PTSANG,WGTANG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the integration weights for cyclic tracking.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NANGLE number of angles.
* PTSANG integration points.
*
*Parameters: output
* WGTANG integration weights.
*
*Reference:
* M.J.Halsall, Cactus, a characteritics solution to the neutron
* transport equations in complicated geometries, UK Atomic Energy
* Authority, 1980.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      INTEGER          NANGLE
      DOUBLE PRECISION PTSANG(NANGLE),WGTANG(NANGLE)
*----
*  Local variables
*----
      INTEGER          LL
      DOUBLE PRECISION XDRCST,PI
      DOUBLE PRECISION ACCUW
*----
*  Computes weights associated with azimuthal angles (PTSANG)
*----
      PI=XDRCST('Pi',' ')
      IF(NANGLE>1) THEN
         WGTANG(1)=1.0D0-(ACOS(PTSANG(2))+ACOS(PTSANG(1)))/PI
         ACCUW=WGTANG(1)
         DO LL=2,NANGLE-1
*           WGTANG(LL)=(PTSANG(LL+1)+PTSANG(LL))/PI
*                     -(PTSANG(LL)+PTSANG(LL-1))/PI
            WGTANG(LL)=(ACOS(PTSANG(LL-1))-ACOS(PTSANG(LL+1)))/PI
            ACCUW=ACCUW+WGTANG(LL)
         ENDDO
         WGTANG(NANGLE)=(ACOS(PTSANG(NANGLE))+ACOS(PTSANG(NANGLE-1)))/PI
         ACCUW=ACCUW+WGTANG(NANGLE)
      ELSE
         WGTANG(1)=1.0D0
      ENDIF
      RETURN
      END
