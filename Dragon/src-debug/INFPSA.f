*DECK INFPSA
      SUBROUTINE INFPSA(IPRINT,TEMPK,PURWGT,PRES)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute pressure at saturation for a mixture of light and heavy water.
*
*Copyright:
* Copyright (C) 2016 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): C. Kieffer and G. Marleau
*
*Parameters: input
* IPRINT  print parameter.
* TEMPK   temperature (kelvin).
* PURWGT  D2O purity (in wgt%).
*
*Parameters: output
* PRES    pressure (Pa).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER  IPRINT
      REAL     TEMPK,PURWGT,PRES
*----
*  LOCAL VARIABLES
*----
      INTEGER  IOUT
      CHARACTER NAMSBR*6
      PARAMETER (IOUT=6,NAMSBR='INFPSA')
      REAL      PD2O,PH2O,WGTD2O,WGTH2O,WGID2O,WGIH2O
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  For heavy water, T limited to 90.5 C < T < 358.5 C
*----
      WGTD2O = 0.01 * PURWGT
      WGTH2O = 1.00 - WGTD2O
      IF(PURWGT .LT. 1.0) THEN
        CALL THMSAP(PH2O,TEMPK)
        WGIH2O=WGTH2O/PH2O
      ELSE
        PH2O=0.0
        WGIH2O=0.0
      ENDIF
      IF(PURWGT .GT. 0.0) THEN
        CALL THMHSP(PD2O,TEMPK)
        IF(PD2O .GT. 0.0) THEN
          WGID2O=WGTD2O/PD2O
        ELSE
          WGID2O=WGTD2O/PH2O
        ENDIF
      ELSE
        WGID2O=0.0
      ENDIF
      PRES = 1.0 /( WGIH2O +WGID2O ) + 0.01
*----
*  Processing finished, return
*----
      IF(IPRINT .GE. -10) THEN
        WRITE(IOUT,*) 'Saturation pressure (Pa)', PRES
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
