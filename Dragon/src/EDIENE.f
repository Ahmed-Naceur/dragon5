*DECK EDIENE
      SUBROUTINE EDIENE(NGROUP,NGCR  ,NGCOND,NTENER,
     >                  IGCR  ,EGCR  ,IGCOND,ENERGY,ENERV )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate energy limits for condensation.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input/output
* NGROUP  number of energy groups.
* NGCR    number of condensed groups read on input.
* NGCOND  number of condensed groups read on EDI.
* NTENER  number of energy found on library.
* IGCR    new group limits.
* EGCR    new energy limits.
* IGCOND  old group limits.
* ENERGY  energy/lethargy/average energy.
* ENERV   average group energy.
*
*-----------------------------------------------------------------------
*
      IMPLICIT    NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     NGROUP,NGCR,NGCOND,NTENER
      INTEGER     IGCR(NGROUP+1),IGCOND(NGROUP+1)
      REAL        EGCR(NGROUP+1),ENERGY(2*NGROUP+1),ENERV(NGROUP)
*----
*  LOCAL VARIABLES
*----
      INTEGER     IOUT
      CHARACTER   NAMSBR*6
      PARAMETER  (IOUT=6,NAMSBR='EDIENE')
      INTEGER     IGC,KDGRP,IGRP,JGRP,IGLIM
*----
*  FIND IF NEW ENERGY OR GROUP SPECIFICATIONS FROM INPUT
*----
      IF(NGCR .GT. 0) THEN
        IGC=0
        IF(EGCR(1) .NE. 0.0) THEN
          IF(NTENER .EQ. 0)  CALL XABORT(NAMSBR//
     >    ': CONDENSATION NOT PERMITTED - NO GROUP STRUCTURE')
          KDGRP=1
          DO 100 IGRP=1,NGROUP+1
            IF(EGCR(IGRP) .LT. ENERGY(NGROUP+1)) THEN
              KDGRP=NGROUP
              IGC=IGC+1
              IGCOND(IGC)=KDGRP
            ELSE IF(EGCR(IGRP) .LT. ENERGY(KDGRP)) THEN
              DO 110 JGRP=KDGRP,NGROUP
                IF(EGCR(IGRP) .GE. ENERGY(JGRP+1)) THEN
                  KDGRP=JGRP
                  IGC=IGC+1
                  IGCOND(IGC)=KDGRP
                  GO TO 115
                ENDIF
 110          CONTINUE
 115          CONTINUE
            ENDIF
            IF(KDGRP .EQ. NGROUP) GO TO 105
 100      CONTINUE
 105      CONTINUE
        ELSE
          DO 120 IGRP=1,NGROUP+1
            IGCOND(IGRP)=IGCR(IGRP)
            IF(IGCR(IGRP) .EQ. NGROUP) THEN
              IGC=IGRP
              GO TO 125
            ENDIF
 120      CONTINUE
 125      CONTINUE
        ENDIF
        NGCOND=IGC
      ENDIF
      IF(NTENER .GT. 0) THEN
*----
*  FIND ENERGY LIMITS, LETHARGY AND AVERAGE ENERGY
*----
        DO 130 IGRP=1,NGROUP
          ENERV(IGRP)=SQRT(ENERGY(IGRP)*ENERGY(IGRP+1))
 130    CONTINUE
        DO 140 IGC=1,NGCOND
          IGLIM=IGCOND(IGC)+1
          ENERGY(IGC+1)=ENERGY(IGLIM)
 140    CONTINUE
        IGLIM=NGCOND+1
        IF(ENERGY(IGLIM) .EQ. 0.0) ENERGY(IGLIM)=1.0E-5
        DO 150 IGC=1,NGCOND
          IGLIM=IGLIM+1
          ENERGY(IGLIM)=LOG(ENERGY(IGC)/ENERGY(IGC+1))
 150    CONTINUE
      ENDIF
      RETURN
      END
