*DECK PSCPUT
      SUBROUTINE PSCPUT(ISPSP,CMDSTR)
C
C---------------------------  PSCPUT  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSCPUT
C      USE      : TRANSFER COMMAND LINE TO FILE
C                 REPLACES PSPLOT ROUTINE FILLER
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      CMDSTR   : COMMAND LINE                           C*132
C    LOCAL
C      IBSL     : ASCII REPRESENTATION OF BACKSLASH      I
C
C---------------------------   PSCPUT  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      CHARACTER        CMDSTR*132
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      INTEGER          IBSL
      PARAMETER       (IBSL=92,NAMSBR='PSCPUT')
      INTEGER          LCMD,IBSLH,ISPACE,IPAREN,IC
      CHARACTER        CBSL*1
      CBSL=CHAR(IBSL)
      LCMD=0
      IBSLH=0
      ISPACE=0
      IPAREN=0
C----
C  COMPRESS COMMAND LINE TO REMOVE USELESS BLANKS
C----
      DO 100 IC=1,132
        IF(CMDSTR(IC:IC) .EQ. ' ' ) THEN
C----
C  REMOVE BLANK IF NOT INSERTED BETWEEN () OR
C  2 OR MORE IN SUCCESSION
C----
          IF(IPAREN .EQ. 0) THEN
            ISPACE=ISPACE+1
          ENDIF
          IF(ISPACE .LE. 1 ) THEN
            LCMD=LCMD+1
            CMDSTR(LCMD:LCMD)=CMDSTR(IC:IC)
          ENDIF
        ELSE
          ISPACE=0
          LCMD=LCMD+1
          CMDSTR(LCMD:LCMD)=CMDSTR(IC:IC)
C----
C  TEST FOR SET OF PARENTHESIS
C  "Backslash"( AND "Backslash") ARE CONSIDERED AS COMMENTED PARENTHESIS
C  AND NOT TREATED
C----
          IF(IBSLH .EQ. 0) THEN
            IF(CMDSTR(IC:IC) .EQ. '(') THEN
              IPAREN=IPAREN+1
            ELSE IF(CMDSTR(IC:IC) .EQ. ')') THEN
              IPAREN=IPAREN-1
            ENDIF
          ENDIF
          IBSLH=0
          IF(CMDSTR(IC:IC) .EQ. CBSL) THEN
            IBSLH=1
          ENDIF
        ENDIF
 100  CONTINUE
C----
C  TEST IF LAST CHARACTER IS A BLANK
C----
      IF(CMDSTR(LCMD:LCMD).EQ. ' ') THEN
        LCMD=LCMD-1
      ENDIF
C----
C  CLEAR REST OF COMMAND STRING AFTER COMPRESSION
C  OF BLANK CHARACTERS
C----
      IF(LCMD .LT. 132) THEN
        CMDSTR(LCMD+1:132)=' '
      ENDIF
C----
C  TRANSFER COMPRESSED COMMAND LINE TO FILE
C----
      IF(LCMD .GT. 0) THEN
        WRITE(ISPSP,'(132A1)')(CMDSTR(IC:IC),IC=1,LCMD)
      ENDIF
      RETURN
      END
