*DECK PSCUTP
      SUBROUTINE PSCUTP(ISPSP)
C
C---------------------------  PSCUTP  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSCUTP
C      USE      : CUT POSTSCRIPT PAGE
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C
C---------------------------   PSCUTP  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      PARAMETER       (NAMSBR='PSCUTP')
      CHARACTER        CMDSTR*132
      CMDSTR='stroke showpage'
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
