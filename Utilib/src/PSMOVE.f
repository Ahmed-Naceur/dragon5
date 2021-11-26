*DECK PSMOVE
      SUBROUTINE PSMOVE(ISPSP,XYPOS,ITMOVE)
C
C---------------------------  PSMOVE  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSMOVE
C      USE      : MOVE PLOT REFERENCE POINT
C                 REPLACES PSPLOT ROUTINE PLOT
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      YXPOS    : FINAL (X,Y) POSITION TO REACH          R(2)
C      ITMOVE   : TYPE OF DISPLACEMENT
C
C---------------------------   PSMOVE  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP,ITMOVE
      REAL             XYPOS(2)
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      REAL             CONVER
      PARAMETER       (NAMSBR='PSMOVE',CONVER=72.0)
      CHARACTER        CMDSTR*132
C----
C  IF ITMOVE=999 TERMINATE PLOT SESSION
C----
      IF(ITMOVE .EQ. 999) THEN
        CMDSTR='stroke showpage'
        CALL PSCPUT(ISPSP,CMDSTR)
        RETURN
      ENDIF
C----
C  MOVE WITH PEN UP (3) OR DOWN (OTHER VALUES) AS REQUESTED
C----
      CMDSTR=' '
      IF(ABS(ITMOVE).EQ.3) THEN
        WRITE(CMDSTR,'(2(F8.2,1X),A2)')
     >    XYPOS(1)*CONVER,XYPOS(2)*CONVER,'SM'
      ELSE
        WRITE(CMDSTR,'(2(F8.2,1X),A3)')
     >    XYPOS(1)*CONVER,XYPOS(2)*CONVER,'LSM'
      ENDIF
      CALL PSCPUT(ISPSP,CMDSTR)
C----
C  RESET ORIGIN IF ITMOVE < 0
C----
      CMDSTR=' '
      IF(ITMOVE.LT.0) THEN
        WRITE(CMDSTR,'(2(F8.2,1X),A9)')
     >    XYPOS(1)*CONVER,XYPOS(2)*CONVER,'translate'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      RETURN
      END
