*DECK PSFCIR
      SUBROUTINE PSFCIR(ISPSP,XYCENT,RADIUS)
C
C---------------------------  PSFCIR  ---------------------------------
C
C  1- PROGRAMME STATISTICS:
C      NAME     : PSFCIR
C      USE      : FILL CIRCLE
C                 REPLACES PSPLOT ROUTINE CIRCLE
C
C  2- ROUTINE PARAMETERS:
C    INPUT/OUTPUT
C      ISPSP    : PSP FILE UNIT                          I
C      YXCENT   : POSITION (X,Y) OF LINE INTERSECTION    R(2)
C      RADIUS   : ARC RADIUS                             R
C
C---------------------------   PSFCIR  --------------------------------
C
      IMPLICIT         NONE
      INTEGER          ISPSP
      REAL             XYCENT(2),RADIUS
C----
C  LOCAL VARIABLES
C----
      CHARACTER        NAMSBR*6
      REAL             CONVER
      PARAMETER       (NAMSBR='PSFCIR',CONVER=72.0)
      CHARACTER        CMDSTR*132
      CMDSTR='Np'
      CALL PSCPUT(ISPSP,CMDSTR)
      CMDSTR=' '
      WRITE(CMDSTR,'(3(F8.2,1X),A8)')
     >  XYCENT(1)*CONVER,XYCENT(2)*CONVER,RADIUS*CONVER,'C fill  '
      CALL PSCPUT(ISPSP,CMDSTR)
      RETURN
      END
