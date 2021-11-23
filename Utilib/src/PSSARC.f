*DECK PSSARC
      SUBROUTINE PSSARC(ISPSP,XYCENT,RADIUS,ANGR)
C
C---------------------------  PSSARC  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSSARC
C      USE      : DRAW CIRCLE
C                 ADAPTED FROM PSPLOT ROUTINE ARC
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      YXCENT   : POSITION (X,Y) OF LINE INTERSECTION    R(2)
C      RADIUS   : ARC RADIUS                             R
C      ANGR     : ARC ANGLE RANGE                        R(2)
C
C---------------------------   PSSARC  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      REAL             XYCENT(2),RADIUS,ANGR(2)
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      REAL             CONVER
      PARAMETER       (NAMSBR='PSSARC',CONVER=72.0)
      CHARACTER        CMDSTR*132
      CMDSTR=' '
      WRITE(CMDSTR,'(5(F8.2,1X),A5)')
     >  XYCENT(1)*CONVER,XYCENT(2)*CONVER,RADIUS*CONVER,
     >  ANGR(1),ANGR(2),'arcit'
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
