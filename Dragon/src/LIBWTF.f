*DECK LIBWTF
      SUBROUTINE LIBWTF(NGROUP,NTMP,TERP,SCAT,SIGS,TMPXS,TMPSC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform temperature interpolation for WIMS-E P1 scattering matrices.
*
*Copyright:
* Copyright (C) 2016 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NGROUP  number of energy groups.
* NTMP    number of temperatures.
* TERP    temperature coefficients.
* TMPXS   temperature dependent vectorial scattering cross sections.
* TMPSC   temperature dependent scattering matrix.
*
*Parameters: output
* SCAT    complete scattering matrix from ig to jg.
* SIGS    scattering cross sections.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* INTERFACE VARIABLES
*----
      INTEGER          NGROUP,NTMP
      DOUBLE PRECISION TERP(NTMP)
      REAL             SCAT(NGROUP,NGROUP),SIGS(NGROUP),
     >                 TMPXS(NGROUP,NTMP),TMPSC(NGROUP,NGROUP,NTMP)
*----
* LOCAL VARIABLES
*----
      INTEGER          IGF,ITM,IGD
      REAL             RTERP
*----
*  INTERPOLATE SCATTERING CROSS SECTIONS IN TEMPERATURE
*----
      CALL XDRSET(SIGS,NGROUP,0.0)
      CALL XDRSET(SCAT,NGROUP*NGROUP,0.0)
      DO 130 ITM=1,NTMP
        RTERP=REAL(TERP(ITM))
        IF(RTERP.NE.0.0D0) THEN
          DO 131 IGD=1,NGROUP
            SIGS(IGD)=SIGS(IGD)+RTERP*TMPXS(IGD,ITM)
            DO 132 IGF=1,NGROUP
              SCAT(IGF,IGD)=SCAT(IGF,IGD)+RTERP*TMPSC(IGF,IGD,ITM)
 132        CONTINUE
 131      CONTINUE
        ENDIF
 130  CONTINUE
      RETURN
      END
