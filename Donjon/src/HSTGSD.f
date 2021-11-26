*DECK HSTGSD
      SUBROUTINE HSTGSD(IPHST,  MAXI,   IOK,    DENI,   FDEN  )
*
*----------
*
*Purpose:
* To read from or write to to history file
* isotopic and fuel densities. 
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IPHST   address of the \dds{history} data structure.
* MAXI    maximum number of isotopes.            
*
*Parameters: input/output
* IOK     processing option where:
*         --> on input, a negative value indicates 
*         that the information is to be extracted 
*         from the \dds{history} data structure and a 
*         positive value indicates that the information is to be
*         stored on the \dds{history} data structure;              
*         --> on output, a value of 0 indicates that 
*         the required processing took place
*         successfully while a negative value indicates 
*         a failure of the processing.
* DENI    isotopic concentration.
* FDEN    average fuel density and weight.
* IOK     status of read. 
*         On input ->  IOK< 0 means get densities
*                             densities
*         On input ->  IOK> 0 means save densities
*         On output -> IOK= 0 success
*                      IOK=-1 error: density missing 
*                      IOK=-2 error: involid processing option
* DENI    initial and final isotopic concentration.
* FDEN    initial fuel density and heavy element mass.
*
*----------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST
      INTEGER          MAXI,IOK
      REAL             DENI(0:MAXI)
      REAL             FDEN(2)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='HSTGSD')
*----
*  LOCAL VARIABLES
*----                  
      INTEGER          ILCMLN,ILCMTY 
*----
*  Local parameters after refuel
*---- 
      IF(IOK .LT. 0) THEN
        IOK=0
*----
*  Get isotopes concentration
*----
        CALL LCMLEN(IPHST,'ISOTOPESDENS',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0 .OR. ILCMLN .GT. MAXI) THEN
          IOK=-1
        ELSE
          CALL LCMGET(IPHST,'ISOTOPESDENS',DENI(1))
        ENDIF
*----
*  Get fuel density
*----
        CALL LCMLEN(IPHST,'FUELDEN-INIT',ILCMLN,ILCMTY)
        IF(ILCMLN .LE. 0 .OR. ILCMLN .GT. 2) THEN
          IOK=-1
        ELSE
          CALL LCMGET(IPHST,'FUELDEN-INIT',FDEN)
        ENDIF
      ELSE IF(IOK .GT. 0) THEN
        IOK=0
*----
*  Put isotopes concentration
*----
        CALL LCMPUT(IPHST,'ISOTOPESDENS',MAXI,2,DENI(1))
*----
*  Put fuel density
*----
        CALL LCMPUT(IPHST,'FUELDEN-INIT',2,2,FDEN)
      ELSE
        IOK=-2
      ENDIF
      RETURN
      END 
