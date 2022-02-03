*DECK PSPFCD
      SUBROUTINE PSPFCD(IPFL,NGROUP,NUNKNO,IGR,ICOND,
     >                  FLUXC,FLUXR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To condense flux for PSP.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPFL    pointer to the FLUX data structure.
* NGROUP  number of groups.
* NUNKNO  number of flux unknowns.
* IGR     condensed group number.
* ICOND   group limit for condensed group.
*
*Parameters: output
* FLUXC   condensed flux.
* FLUXR   multigroup flux read.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='PSPFCD')
*----
*  ROUTINE PARAMETERS
*----
      TYPE(C_PTR)      IPFL
      INTEGER          NGROUP,NUNKNO,IGR
      INTEGER          ICOND(NGROUP)
      REAL             FLUXC(NUNKNO),FLUXR(NUNKNO)
*----
*  LOCAL PARAMETERS
*----
      TYPE(C_PTR)      JPFL
      INTEGER          IGC,IGD,IGF,IUN
*
      CALL XDRSET(FLUXC,NUNKNO,0.0)
      IGF=ICOND(IGR)
      IF(IGR .EQ. 1) THEN
        IGD=1
      ELSE
        IGD=ICOND(IGR-1)+1
      ENDIF
      JPFL=LCMGID(IPFL,'FLUX')
      DO IGC=IGD,IGF
        CALL LCMGDL(JPFL,IGC,FLUXR)
        DO IUN=1,NUNKNO
           FLUXC(IUN)=FLUXC(IUN)+FLUXR(IUN)
        ENDDO
      ENDDO
*----
*  RETURN
*----
      RETURN
      END
