*DECK EDIDEL
      SUBROUTINE EDIDEL(IPRINT,NGCOND,NMERGE,IGR,SCATTS,INGSCT,IFGSCT,
     >                  IPOSCT,XSCAT,DELSCT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Scattering delta sigma calculation.
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
*Parameters: input
* IPRINT  print level (=1 to print delta sigmas).
* NGCOND  number of condensed groups.
* NMERGE  number of merged regions.
* IGR     current group treated.
* SCATTS  new scattering matrix.
*
*Parameters: scratch
* INGSCT  number of scattering groups (old).
* IFGSCT  first old scattering group.
* IPOSCT  position of region in scattering matrix.
* XSCAT   new scattering matrix.
* DELSCT  scattering delta sigma.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     IPRINT,NGCOND,NMERGE,IGR,INGSCT(NMERGE),
     >            IFGSCT(NMERGE),IPOSCT(NMERGE)
      REAL        SCATTS(NMERGE,NGCOND,NGCOND),XSCAT(NMERGE*NGCOND),
     >            DELSCT(NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6)
*----
*  COMPUTE DELTA SIGMA BY MATERIAL FOR ALL GROUPS
*----
      DO 100 IREG=1,NMERGE
        JLASTG=IFGSCT(IREG)
        JFRSTG=JLASTG-INGSCT(IREG)+1
        JPOSCT=IPOSCT(IREG)+INGSCT(IREG)
        IF(IPRINT.GE.1) THEN
          WRITE(IUNOUT,6000) IREG,NGCOND
        ENDIF
        DO 110 JGR=1,NGCOND
          IF( (JGR.LT.JFRSTG) .OR. (JGR.GT.JLASTG) ) THEN
            DELSCT(JGR)=SCATTS(IREG,IGR,JGR)
          ELSE
            JPOSCT=JPOSCT-1
            DELSCT(JGR)=SCATTS(IREG,IGR,JGR)-XSCAT(JPOSCT)
          ENDIF
 110    CONTINUE
 100  CONTINUE
      WRITE(IUNOUT,6001) (DELSCT(JGR),JGR=1,NGCOND)
*----
*  FORMAT
*----
 6000 FORMAT(/' REGION=',I5,10X,'NUMBER OF GROUPS ',I10)
 6001 FORMAT(1P,7(3X,E15.7))
      RETURN
      END
