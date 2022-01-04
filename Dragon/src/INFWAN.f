*DECK INFWAN
      SUBROUTINE INFWAN(TEMPK,PURWGT,PRES,DENSITY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute water density as a function of temperature and pressure.
*
*Copyright:
* Copyright (C) 2016 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): R. Roy and C. Kieffer
*
*Parameters: input
* TEMPK   temperature (kelvin).
* PURWGT  D2O purity (in wgt%).
* PRES    pressure (Pa).
*
*Parameters: output
* DENSITY density (G/CM**3).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
      REAL TEMPK,PRES,PURWGT,DENSITY,R2,R3,R4,R5
      REAL TEMPC,h, zk, zmu, cp
      REAL DEND2O, DENH2O, WGTD2O, WGTH2O
*----------------------------------------------------------
*  Light water density
*----------------------------------------------------------
      CALL THMPT(PRES,TEMPK,DENH2O,h, zk, zmu, cp)
      IF(DENH2O .EQ.0 ) THEN
        DENH2O=0.00000001
      ENDIF
*----------------------------------------------------------
*  Heavy water density
*----------------------------------------------------------
      TEMPC=TEMPK-273.15
      IF(TEMPC.GT.358.5 .OR. TEMPC.LT.90.5 .OR. PRES .GT. 22.0E6) THEN
           DEND2O = 1.11 * DENH2O
      ELSE
           CALL THMHPT(PRES,TEMPK,DEND2O,R2,R3,R4,R5)
      ENDIF
*----------------------------------------------------------
*  Global density for the mixture
*----------------------------------------------------------
      WGTD2O = 0.01 * PURWGT
      WGTH2O = 1.00 - WGTD2O
      IF(PURWGT .EQ. 1.0) THEN
        DENSITY=DEND2O
      ELSE IF(PURWGT .EQ. 0.0) THEN
        DENSITY=DENH2O
      ELSE
        DENSITY=DENH2O*DEND2O/(WGTH2O*DEND2O+WGTD2O*DENH2O)
      ENDIF
      DENSITY=DENSITY/1000.0
      RETURN
      END
