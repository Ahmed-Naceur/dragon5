*DECK PSPAGE
      SUBROUTINE PSPAGE(ISPSP,NPAGE,XYPOS)
C
C---------------------------  PSPAGE  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSPAGE
C      USE      : SET POSTSCRIPT NEW PAGE
C                 REPLACES PART OF PSPLOT ROUTINE CHOPIT
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      NPAGE    : PAGE NUMBER                            I
C      XYPOS    : ORIGIN OF PAGE                         R(2)
C---------------------------   PSPAGE  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP,NPAGE
      REAL             XYPOS(2)
C----
C  LOCAL VARIABLES
C----
      INTEGER          MCOLR
      CHARACTER        NAMSBR*6
      PARAMETER       (MCOLR=4,NAMSBR='PSPAGE')
      INTEGER          ICOLR
      CHARACTER        CMDSTR*132
      REAL             COLR(MCOLR)
C----
C  INITIALIZE COLOR, XYPOS TO 0.0
C----
      DO 100 ICOLR=1,MCOLR
        COLR(ICOLR)=0.0
 100  CONTINUE
      CMDSTR=' '
      WRITE(CMDSTR,'(A8,2(I4,1X))') '%%Page: ',NPAGE,NPAGE
      CALL PSCPUT(ISPSP,CMDSTR)
      IF(NPAGE .GT. 1) THEN
        CMDSTR='newpath 0 0 moveto'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      CMDSTR='/Curfnt /Times-Italic findfont def'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' 12 Setf'
      CALL PSCPUT(ISPSP,CMDSTR)
      IF(NPAGE .EQ. 1) THEN
        CMDSTR=' 1.000 1.000  scale'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      CALL PSLINW(ISPSP,COLR)
      CMDSTR=' '
      WRITE(CMDSTR,'(F7.3,1X,A5)') 0.0,'Sgray'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' '
      WRITE(CMDSTR,'(3(F7.3,1X),A5)') 0.0,0.0,0.0,'Scrgb'
      CALL PSCPUT(ISPSP,CMDSTR)
      CALL PSMOVE(ISPSP,XYPOS,-3)
      RETURN
      END
