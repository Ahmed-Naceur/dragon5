*DECK PSLINW
      SUBROUTINE PSLINW(ISPSP,WLINE)
C
C---------------------------  PSLINW  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSLINW
C      USE      : SET LINE WIDTH
C                 REPLACES PSPLOT ROUTINE SETLW
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      WLINE    : WIDTH OF LINE                          R
C
C---------------------------   PSLINW  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      REAL             WLINE
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      PARAMETER       (NAMSBR='PSLINW')
      CHARACTER        CMDSTR*132
C----
C  SET LINE WIDTH
C  MINIMUM IS 0.00001
C----
      CMDSTR=' '
      WRITE(CMDSTR,'(F7.3,1X,A5)') 0.0,'Sgray'
      CALL PSCPUT(ISPSP,CMDSTR)
      IF(ABS(WLINE).LT.1.E-5) THEN
        CMDSTR='Slw0'
      ELSE
        CMDSTR=' '
        WRITE(CMDSTR,'(F7.3,1X,A3)') WLINE,'Slw'
      ENDIF
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
