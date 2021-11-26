*DECK HSTGMA
      SUBROUTINE HSTGMA(IPMAP,  NCHA,   NBUN,  DELTAT,   POWER,
     >                  BURNP,  IREFUS, REFUT, NBFUEL)
*
*----------
*
*Purpose:
* To read from the MAP data structure the power and 
* burnup distribution for each cell as well as the refueling 
* option for each channel.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau, E. Varin
*
*Parameters: input
* IPMAP   address of the \dds{map} data structure.
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.            
* IPMAP   pointer to the MAP data structure
* NCHA    number of fuel channels.
* NBUN    number of bundles per channels.
*
*Parameters: input/output
* DELTAT  last character string read.
* POWER   power for each fuel bundle in each channel.
* BURNP   burnup for each fuel bundle in each channel.
* IREFUS  refueling strategy for each channel.
* REFUT   refueling time for each channel.
* NBFUEL  number of fueled channels.            
* DELTAT  next time steps for burnup.
* POWER   values of local powers.
* IREFUS  fuels shift for each channel. 
*         A channel is refueled using a NBS bundle 
*         shift procedure if IREFUS(I)=NBS. 
*         In the case where NBS $>$ 0,
*         bundles 1 to NBUN-NBS are displaced to position NBS+1 to
*         NBUN while locations 1 to NBS are filled with new fuel. 
*         In the case where NBS $<$ 0,
*         bundles -NBS+1 to NBUN are displaced to position 1 to
*         NBUN+NBS while locations NBUN+NBS+1 to NBUN are filled 
*         with new fuel.
* REFUT   channel refueling time.
* NBFUEL  number of fueled channels 
*
*----------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPMAP
      INTEGER          NCHA,NBUN
      REAL             DELTAT
      REAL             POWER(NCHA,NBUN),BURNP(NCHA,NBUN)
      INTEGER          IREFUS(NCHA)
      REAL             REFUT(NCHA)
      INTEGER          NBFUEL
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT,NTC,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NTC=3,ILCMUP=1,ILCMDN=2,
     >                 NAMSBR='HSTGMA')
*----
*  LOCAL VARIABLES
*----                  
      INTEGER          ILCMLN,ILCMTY
      INTEGER          IC 
*----
*  Read DEPL-TIME
*----
      CALL LCMLEN(IPMAP,'DEPL-TIME   ',ILCMLN,ILCMTY)
      IF(ILCMLN .GT. 0) THEN 
        IF(ILCMLN .GT. 1) CALL XABORT(NAMSBR//
     >  ': Space to store next time step is too small')
        CALL LCMGET(IPMAP,'DEPL-TIME   ',DELTAT)
      ENDIF
*----
*  Read bundle powers
*----
      CALL LCMLEN(IPMAP,'BUND-PW',ILCMLN,ILCMTY)
      IF(ILCMLN .GT. 0) THEN 
        IF(ILCMLN .GT. NCHA*NBUN) CALL XABORT(NAMSBR//
     >  ': Space to store power is too small')
        CALL LCMGET(IPMAP,'BUND-PW',POWER)
      ENDIF
*----
*  Read BURNUP IF DELTAT=0.0
*----
      CALL XDRSET(BURNP,NCHA*NBUN,0.0)
      IF(DELTAT.EQ.0.0) THEN
        CALL LCMLEN(IPMAP,'BURN-INST',ILCMLN,ILCMTY)
        IF(ILCMLN .GT. 0) THEN 
          IF(ILCMLN .GT. NCHA*NBUN) CALL XABORT(NAMSBR//
     >    ': Space to store burnup is too small')
          CALL LCMGET(IPMAP,'BURN-INST',BURNP)
        ENDIF
      ENDIF
*----
*  Read refueling scheme
*----
      CALL LCMLEN(IPMAP,'REF-SCHEME',ILCMLN,ILCMTY)
      IF(ILCMLN .GT. 0) THEN 
        IF(ILCMLN .GT. NCHA) CALL XABORT(NAMSBR//
     >  ': Space to store REF-SCHEME is too small')
        CALL LCMGET(IPMAP,'REF-SCHEME',IREFUS)
      ENDIF
      CALL LCMLEN(IPMAP,'REF-CHAN',ILCMLN,ILCMTY)
      IF(ILCMLN .GT. 0) THEN 
        IF(ILCMLN .GT. NCHA) CALL XABORT(NAMSBR//
     >  ': Space to store REF-CHAN is too small')
        CALL LCMGET(IPMAP,'REF-CHAN',REFUT)
      ENDIF
*----
*  Compute number of channels refueled
*----
      DO 100 IC=1,NCHA
        IF(REFUT(IC) .GT. 0.0) NBFUEL=NBFUEL+1
 100  CONTINUE        
      RETURN
      END 
