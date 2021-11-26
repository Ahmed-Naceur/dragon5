*DECK PSFILL
      SUBROUTINE PSFILL(ISPSP,IFILL,GRYCOL,KFS,KFR)
C
C---------------------------  PSFILL  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSFILL
C      USE      : SET GRAY LEVEL OR COLOR AND FILL PATERN
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      IFILL    : FILL TYPE                              I
C                 = 0 SET TO COLOR(4)
C                 = 1 SET TO GRAY LEVEL
C                 = 2 SET RGB COLLOR PATTERN
C                 = 3 SET GRYCOL COLLOR PATTERN
C                 = 4 SET HSB COLLOR PATTERN
C      GRYCOL   : GRAY LEVEL OF COLOR INTENSITY          R(4)
C      KFS      : FLAG TO SAVE DRAWING BEFORE FILLING    I
C                 = 0 : NO SAVE
C                 = 1 : SAVE  
C      KFR      : FLAG TO RESTORE DRAWING BEFORE FILLING I
C                 = 0 : NO RESTORE
C                 = 1 : RESTORE 
C
C---------------------------   PSFILL  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP,IFILL
      REAL             GRYCOL(4)
      INTEGER          KFS,KFR
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      PARAMETER       (NAMSBR='PSFILL')
      REAL             COLOR(4)
      CHARACTER        CMDSTR*132
C----
C  TAKE COLOR LEVEL BETWEEN 0.0 AND 1.0
C----
      IF(KFR .EQ. 1) THEN
        CMDSTR='grestore'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      IF(KFS .EQ. 1) THEN
        CMDSTR='gsave'
        CALL PSCPUT(ISPSP,CMDSTR)
      ENDIF
      COLOR(1)=MIN(1.0,ABS(GRYCOL(1)))
      COLOR(2)=MIN(1.0,ABS(GRYCOL(2)))
      COLOR(3)=MIN(1.0,ABS(GRYCOL(3)))
      COLOR(4)=MIN(1.0,ABS(GRYCOL(4)))
      COLOR(1)=MAX(0.0,COLOR(1))
      COLOR(2)=MAX(0.0,COLOR(2))
      COLOR(3)=MAX(0.0,COLOR(3))
      COLOR(4)=MAX(0.0,COLOR(4))
      CMDSTR=' '
      IF(IFILL .EQ.4) THEN
        WRITE(CMDSTR,'(3(F7.3,1X),A6)')
     >    COLOR(1),COLOR(2),COLOR(3),'FSchsb'
      ELSE IF(IFILL.EQ.3) THEN
        WRITE(CMDSTR,'(4(F7.3,1X),A6)')
     >    COLOR(1),COLOR(2),COLOR(3),COLOR(4),'FScmyk'
      ELSE IF(IFILL.EQ.2) THEN
        WRITE(CMDSTR,'(3(F7.3,1X),A6)')
     >    COLOR(1),COLOR(2),COLOR(3),'FScrgb'
      ELSE IF(IFILL.EQ.1) THEN
        WRITE(CMDSTR,'(1(F7.3,1X),A6)')
     >    COLOR(1),'FSgray'
      ELSE
        WRITE(CMDSTR,'(1(F7.3,1X),A6)')
     >    0.0,'FSgray'
      ENDIF
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
