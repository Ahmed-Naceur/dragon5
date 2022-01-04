*DECK RODMOV
      SUBROUTINE RODMOV(IPMAP,NGRP,RNAME,INS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modify rod insertion (second call)
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
* RNAME   name of rod group
* INS     rod insertion for each rod group
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NGRP,INS(NGRP)
      CHARACTER(LEN=3) RNAME(NGRP)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,J,NRMV,INS2(NGRP)
      REAL FLOT
      CHARACTER(LEN=3) RNAME2(NGRP)
      TYPE(C_PTR) MPMAP
      CHARACTER TEXT*3
      DOUBLE PRECISION DFLOT
*
      MPMAP=LCMGID(IPMAP,'ROD-INFO')
      CALL LCMGTC(MPMAP,'ROD-NAME',3,NGRP,RNAME)
      CALL LCMGET(MPMAP,'ROD-INSERT',INS)
      CALL REDGET(ITYP,NRMV,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1)CALL XABORT('@RODMOV: INTEGER'
     1 //' DATA FOR GROUP NUMBER EXPECTED.')
      J=1
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@RODMOV: CHARACTER DATA EXPECTED.')
      DO WHILE(J.LE.NGRP)
        RNAME2(J)=TEXT
        I=1
        DO WHILE (I.LE.NRMV)
          RNAME2(I)=TEXT
          IF(RNAME2(I).EQ.RNAME(J)) THEN
            CALL REDGET(ITYP,INS2(J),FLOT,TEXT,DFLOT)
            IF(ITYP.NE.1)CALL XABORT('@RODMOV: INTEGER DA'
     1       //'TA FOR INS EXPECTED.')
          CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@RODMOV: CHARACTER DATA EXPECTED.')
            GO TO 10
          ELSE
            I=I+1
          ENDIF
        END DO
        INS2(J)=INS(J)
   10   J=J+1
      END DO
      CALL LCMPUT(MPMAP,'ROD-INSERT',NGRP,1,INS2)
      RETURN
      END
