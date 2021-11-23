*DECK RODTYP
      SUBROUTINE RODTYP(IPMAP,NGRP,MAXMIX,RNAME,INS,HMIX,RMIX,NUMMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store rod parameters in the fuelmap
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Tixier
*
*Parameters: input
* IPMAP   pointer to the fuel map
* NGRP    number of rod groups
* MAXMIX  maximum number of rod mix
* RNAME   name of rod group
* INS     rod insertion for each rod group
* HMIX    height of rod mix (if more than 2 rod mix are defined)
* RMIX    number associated to rod mix
* NUMMIX  number of rod mix for each rod group
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NGRP,MAXMIX,INS(NGRP),NUMMIX(NGRP)
      REAL HMIX(NGRP*MAXMIX),RMIX(NGRP*MAXMIX)
      CHARACTER(LEN=3) RNAME(NGRP)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,J,INSM
      REAL FLOT
      TYPE(C_PTR) MPMAP
      CHARACTER TEXT*3
      DOUBLE PRECISION DFLOT
*
      MPMAP=LCMGID(IPMAP,'ROD-INFO')
      CALL LCMGET(MPMAP,'INS-MAX',INSM)
      I=1
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@RODTYP: CHARACTER DATA EXPECTED.')
      DO WHILE (I.LE.NGRP)
        J=1
        RNAME(I)=TEXT
        CALL REDGET(ITYP,INS(I),FLOT,TEXT,DFLOT)
        IF(ITYP.NE.1)CALL XABORT('@RODTYP: INTEGER DA'
     1   //'TA FOR INS EXPECTED.')
        IF(INS(I).GT.INSM)CALL XABORT('@RODTYP: ROD INSERTION IS '
     1   //'LARGER THAN MAXIMUM INSERTION.')
        IF(INS(I).LT.0)CALL XABORT('@RODTYP: ROD INSERTION MUST BE '
     1   //'POSITIVE.')
        CALL REDGET(ITYP,NITMA,RMIX(I+(J-1)*NGRP),TEXT,DFLOT)
        IF(ITYP.NE.2)CALL XABORT('@RODTYP: REAL DA'
     1   //'TA FOR MIX EXPECTED.')
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        DO WHILE (ITYP.NE.3)
          HMIX(I+(J-1)*NGRP)=FLOT
          IF((J.GE.2).AND.(HMIX(I+(J-2)*NGRP).GT.HMIX(I+(J-1)*NGRP)))
     >    CALL XABORT('@RODTYP: THE LENGTH OF THE TOP MATERIAL MUST' 
     1   //'BE HIGHER THAN THE BOTTOM MATERIAL.')
          J=J+1
          CALL REDGET(ITYP,NITMA,RMIX(I+(J-1)*NGRP),TEXT,DFLOT)
          IF(ITYP.NE.2)CALL XABORT('@RODTYP: REAL DATA FOR MIX '
     1   //'EXPECTED.')
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        END DO
          NUMMIX(I)=J
          I=I+1
      END DO
      CALL LCMPTC(MPMAP,'ROD-NAME',3,NGRP,RNAME)
      CALL LCMPUT(MPMAP,'ROD-INSERT',NGRP,1,INS)
      CALL LCMPUT(MPMAP,'ROD-RIN',MAXMIX*NGRP,2,RMIX)
      CALL LCMPUT(MPMAP,'ROD-NBZONE',NGRP,1,NUMMIX)
      CALL LCMPUT(MPMAP,'ROD-HEIGHT',MAXMIX*NGRP,2,HMIX)
      RETURN
      END
