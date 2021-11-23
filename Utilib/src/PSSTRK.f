*DECK PSSTRK
      SUBROUTINE PSSTRK(ISPSP,WLINE,KSS,KSR)
C
C---------------------------  PSSTRK  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSSTRK
C      USE      : SET LINE WIDTH AND STROKE PATH
C                 REPLACES PSPLOT ROUTINE SETLW
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      WLINE    : WIDTH OF LINE                          R
C      KSS      : FLAG TO SAVE DRAWING BEFORE FILLING    I
C                 = 0 : NO SAVE
C                 = 1 : SAVE  
C      KSR      : FLAG TO RESTORE DRAWING BEFORE FILLING I
C                 = 0 : NO RESTORE
C                 = 1 : RESTORE 
C
C---------------------------   PSSTRK  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      REAL             WLINE
      INTEGER          KSS,KSR
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      PARAMETER       (NAMSBR='PSSTRK')
      CHARACTER        CMDSTR*132
C----
C  SET LINE WIDTH
C  MINIMUM IS 0.00001
C----
      IF(KSR .EQ. 1) THEN
        CMDSTR='grestore'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      IF(KSS .EQ. 1) THEN
        CMDSTR='gsave'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      CMDSTR=' '
      WRITE(CMDSTR,'(F7.3,1X,A5)') 0.0,'Sgray'
      CALL PSCPUT(ISPSP,CMDSTR)
      IF(ABS(WLINE).LT.1.E-5) THEN
        CMDSTR='SSlw0'
      ELSE
        CMDSTR=' '
        WRITE(CMDSTR,'(F7.3,1X,A4)') WLINE,'SSlw'
      ENDIF
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
