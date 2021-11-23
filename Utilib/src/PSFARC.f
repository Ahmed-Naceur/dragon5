*DECK PSFARC
      SUBROUTINE PSFARC(ISPSP,XYCENT,RADIUS,ANGR)
C
C---------------------------  PSFARC  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSFARC
C      USE      : FILL ARC
C                 REPLACES PSPLOT ROUTINE ARC
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      YXCENT   : POSITION (X,Y) OF LINE INTERSECTION    R(2)
C      RADIUS   : ARC RADIUS                             R
C      ANGR     : ARC ANGLE RANGE                        R(2)
C
C---------------------------   PSFARC  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      REAL             XYCENT(2),RADIUS,ANGR(2)
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      REAL             CONVER
      PARAMETER       (NAMSBR='PSFARC',CONVER=72.0)
      CHARACTER        CMDSTR*132
      CMDSTR='Np'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' '
      WRITE(CMDSTR,'(5(F8.2,1X),A3)')
     >  XYCENT(1)*CONVER,XYCENT(2)*CONVER,RADIUS*CONVER,
     >  ANGR(1),ANGR(2),'arc'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR='cf'
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
