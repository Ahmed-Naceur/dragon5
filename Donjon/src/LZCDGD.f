*DECK LZCDGD
      SUBROUTINE LZCDGD(IPDEV,NLZC,LGRP,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create the liquid-zone-controllers group directories on the device
* data structure.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPDEV  pointer to device information.
* NLZC   total number of liquid zone controllers.
* LGRP   total number of lzc-groups.
* IMPX   printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDEV
      INTEGER NLZC,LGRP,IMPX
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      CHARACTER TEXT*12
      INTEGER LZCID(NLZC)
      DOUBLE PRECISION DFLOT
      TYPE(C_PTR) JPDEV,KPDEV
*----
*  CREATE GROUPS
*----
      JPDEV=LCMLID(IPDEV,'LZC_GROUP',LGRP)
      IGRP=0
      IF(IMPX.GT.0)WRITE(IOUT,1001)
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@LZCDGD: KEYWORD GROUP-ID EXPECTED.')
      IF(TEXT.NE.'GROUP-ID')CALL XABORT('@LZCDGD: KEYWORD GROUP-'
     1 //'ID EXPECTED.')
   10 IGRP=IGRP+1
      CALL REDGET(ITYP,JGRP,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.1)CALL XABORT('@LZCDGD: INTEGER GROUP-ID NUMBER'
     1 //' EXPECTED.')
      IF(JGRP.NE.IGRP)THEN
        WRITE(IOUT,*)'@LZCDGD: READ GROUP-ID NUMBER #',JGRP
        WRITE(IOUT,*)'@LZCDGD: EXPECTED GROUP-ID NUMBER #',IGRP
        CALL XABORT('@LZCDGD: WRONG GROUP-ID NUMBER.')
      ENDIF
      IF(JGRP.GT.LGRP)THEN
        WRITE(IOUT,*)'@LZCDGD: GIVEN TOTAL NUMBER OF GROUPS ',LGRP
        WRITE(IOUT,*)'@LZCDGD: READ GROUP-ID NUMBER #',JGRP
        CALL XABORT('@LZCDGD: WRONG GROUP-ID NUMBER.')
      ENDIF
      CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
      IF(ITYP.NE.3)CALL XABORT('@LZCDGD: KEYWORD EXPECTED.')
*----
*  OPTION ALL
*----
      IF(TEXT.EQ.'ALL')THEN
        KPDEV=LCMDIL(JPDEV,IGRP)
        DO 30 ID=1,NLZC
        LZCID(ID)=ID
   30   CONTINUE
        CALL LCMPUT(KPDEV,'GROUP-ID',1,1,IGRP)
        CALL LCMPUT(KPDEV,'NUM-LZC',1,1,NLZC)
        CALL LCMPUT(KPDEV,'LZC-ID',NLZC,1,LZCID)
*
        CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.NE.3)CALL XABORT('@LZCDGD: WRONG INPUT DATA.')
        IF(TEXT.EQ.';')THEN
          IF(IGRP.EQ.LGRP)THEN
            NDG=NLZC
            GOTO 100
          ENDIF
          WRITE(IOUT,*)'@LZCDGD: GIVEN TOTAL NUMBER OF GROUPS ',LGRP
          WRITE(IOUT,*)'@LZCDGD: CREATED ONLY NUMBER OF GROUPS ',IGRP
          CALL XABORT('@LZCDGD: WRONG NUMBER OF GROUPS.')
        ELSEIF(TEXT.EQ.'GROUP-ID')THEN
          IF(IMPX.GT.0)WRITE(IOUT,1000)IGRP,NLZC
          GOTO 10
        ELSE
          CALL XABORT('@LZCDGD: WRONG KEYWORD '//TEXT)
        ENDIF
*----
*  OPTION LZC-ID
*----
      ELSEIF(TEXT.EQ.'LZC-ID')THEN
        NDG=0
        CALL XDISET(LZCID,NLZC,0)
        KPDEV=LCMDIL(JPDEV,IGRP)
*
   50   CALL REDGET(ITYP,NITMA,FLOT,TEXT,DFLOT)
        IF(ITYP.EQ.3)THEN
          IF(TEXT.EQ.';')THEN
            IF(IGRP.EQ.LGRP)GOTO 100
            WRITE(IOUT,*)'@LZCDGD: GIVEN TOTAL NUMBER OF GROUPS ',LGRP
            WRITE(IOUT,*)'@LZCDGD: CREATED ONLY NUMBER OF GROUPS ',IGRP
            CALL XABORT('@LZCDGD: WRONG NUMBER OF GROUPS.')
          ELSEIF(TEXT.EQ.'GROUP-ID')THEN
            IF(IMPX.GT.0)WRITE(IOUT,1000)IGRP,NDG
            GOTO 10
          ELSE
            CALL XABORT('@LZCDGD: WRONG KEYWORD '//TEXT)
          ENDIF
*----
*  LZC-ID NUMBERS
*----
        ELSEIF(ITYP.EQ.1)THEN
          ID=NITMA
          IF((ID.GT.NLZC).OR.(ID.LE.0))THEN
            WRITE(IOUT,*)'@LZCDGD: FOR THE GROUP #',IGRP
            WRITE(IOUT,*)'@LZCDGD: READ LZC-ID #',ID
            CALL XABORT('@LZCDGD: WRONG LZC-ID NUMBER.')
          ENDIF
          DO I=1,NLZC
            IF(ID.EQ.LZCID(I))THEN
              WRITE(IOUT,*)'@LZCDGD: FOR THE GROUP #',IGRP
              WRITE(IOUT,*)'@LZCDGD: REPEATED LZC-ID #',ID
              CALL XABORT('@LZCDGD: WRONG LZC-ID NUMBER.')
            ENDIF
          ENDDO
*
          NDG=NDG+1
          IF(NDG.GT.NLZC)THEN
            WRITE(IOUT,*)'@LZCDGD: FOR THE GROUP #',IGRP
            WRITE(IOUT,*)'@LZCDGD: WRONG TOTAL NUMBER OF LZC ',NDG
            CALL XABORT('@LZCDGD: INVALID INPUT OF LZC-DEVICES.')
          ENDIF
          LZCID(NDG)=ID
          CALL LCMPUT(KPDEV,'GROUP-ID',1,1,IGRP)
          CALL LCMPUT(KPDEV,'NUM-LZC',1,1,NDG)
          CALL LCMPUT(KPDEV,'LZC-ID',NDG,1,LZCID)
        ELSE
          CALL XABORT('@LZCDGD: WRONG INPUT DATA.')
        ENDIF
        GOTO 50
      ELSE
        CALL XABORT('@LZCDGD: WRONG KEYWORD '//TEXT)
      ENDIF
  100 IF(IMPX.GT.0)WRITE(IOUT,1000)IGRP,NDG
      IF(IMPX.GT.0)WRITE(IOUT,1002)LGRP
      RETURN
*
 1000 FORMAT(/1X,'CREATED A GROUP #',I2.2,
     1        4X,'INCLUDES TOTAL NUMBER OF LZC:',1X,I2)
 1001 FORMAT(/1X,'**  CREATING GROUPS FOR LZC-DEVICES  **')
 1002 FORMAT(/1X,39('-')/1X,'TOTAL NUMBER OF GROUPS CREATED:',I2)
      END
