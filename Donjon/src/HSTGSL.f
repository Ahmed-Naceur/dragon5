*DECK HSTGSL
      SUBROUTINE HSTGSL(IPHST, MAXL,   IOK,    TIMPOW, PARAML)
*
*----------
*
*Purpose:
* To read from or save to history file the local parameters 
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IPHST   address of the \dds{history} data structure.
* MAXL    maximum number of local parameters.                   
*
*Parameters: input/output
* IOK     processing option where:
*         --> on input, a negative value indicates that the 
*         information is to be extracted from the \dds{history} data
*         structure and a positive value indicates that the information
*         is to be stored on the \dds{history} data structure 
*         (-1 and 1 for before refueling and -2, 2 for after refueling);             
*         --> on output, a value of 0 indicates that the required
*         processing took place successfully while a negative 
*         value indicates a failure of the processing. 
* TIMPOW  burnup time and power density.
* PARAML  local parameters.
*
*----------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST
      INTEGER          MAXL,IOK
      REAL             PARAML(0:MAXL)
      REAL             TIMPOW(2)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='HSTGSL')
*----
*  LOCAL VARIABLES
*----                  
      INTEGER          ILCMLN,ILCMTY 
*----
*  Local parameters after refuel
*---- 
      IF(IOK .EQ. -2) THEN
*----
*  Get local parameters after refuel
*---- 
        CALL LCMLEN(IPHST,'PARAMLOCALAR',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0 .OR. ILCMLN .GT. MAXL) THEN
          IOK=-1
        ELSE 
          CALL LCMGET(IPHST,'PARAMLOCALAR',PARAML(1))
          IOK=0
        ENDIF
        CALL LCMLEN(IPHST,'PARAMBURNTAR',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0 .OR. ILCMLN .GT. 2) THEN
          IOK=-1
        ELSE 
          CALL LCMGET(IPHST,'PARAMBURNTAR',TIMPOW)
          IOK=0
        ENDIF
      ELSE IF(IOK .EQ. -1) THEN
*----
*  Get local parameters before refuel
*---- 
        PARAML(0)=0
        CALL LCMLEN(IPHST,'PARAMLOCALBR',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0 .OR. ILCMLN .GT. MAXL) THEN
          IOK=-1
        ELSE 
          CALL LCMGET(IPHST,'PARAMLOCALBR',PARAML(1))
          IOK=0
        ENDIF
        CALL LCMLEN(IPHST,'PARAMBURNTBR',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0 .OR. ILCMLN .GT. 2) THEN
          IOK=-1
        ELSE 
          CALL LCMGET(IPHST,'PARAMBURNTBR',TIMPOW)
          IOK=0
        ENDIF
      ELSE IF(IOK .EQ. 1) THEN
*----
*  Save local parameters before refuel
*---- 
        CALL LCMPUT(IPHST,'PARAMLOCALBR',MAXL,2,PARAML(1))
        CALL LCMPUT(IPHST,'PARAMBURNTBR',2,2,TIMPOW)
        IOK=0
      ELSE IF(IOK .EQ. 2) THEN
*----
*  Save local parameters after refuel
*---- 
        CALL LCMPUT(IPHST,'PARAMLOCALAR',MAXL,2,PARAML(1))
        CALL LCMPUT(IPHST,'PARAMBURNTAR',2,2,TIMPOW)
        IOK=0 
      ELSE
        IOK=-2
      ENDIF
      RETURN
      END 
