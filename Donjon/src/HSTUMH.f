*DECK HSTUMH                                 
      SUBROUTINE HSTUMH(IPMAP, IPHST,  IPRINT, NCHA,   NBUN,   IDCELL,
     >                  BURNUP )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To update the MAP data structure using the information
* provided on the HISTORY data structure.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin
*
*Parameters: input
* IPMAP   address of the \dds{map} data structure.
* IPHST   address of the \dds{history} data structure.
* IPRINT  print level.
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.
* IDCELL  cell identifier for each fuel bundle in each channel.
*
*Parameters: work        
* BURNUP   burnup for each fuel bundle in each channel.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPHST,IPMAP
      INTEGER          IPRINT
      INTEGER          NCHA,NBUN
      INTEGER          IDCELL(NBUN,NCHA)
      REAL             BURNUP(NCHA,NBUN)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT
      INTEGER          ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,ILCMUP=1,ILCMDN=2,NAMSBR='HSTUMH')
*----
*  LOCAL VARIABLES
*---- 
      CHARACTER        NAMP*12 
      INTEGER          ILCMLN,ILCMTY 
      INTEGER          IBT,ICT,ICCT
      REAL             BITH(3)
*----
*  Read isotope densities on CELL TYPE
*  if available
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL LCMGET(IPMAP,'BURN-DEB',BURNUP)
      DO 10 ICT=1,NCHA
       DO 20 IBT=1,NBUN
        ICCT=IDCELL(IBT,ICT)
        WRITE(NAMP,'(A4,I8.8)') 'CELL',ICCT
        CALL LCMLEN(IPHST,NAMP,ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 0) THEN
          CALL XABORT(' HSTUMH: BAD CELL TYPE')
        ELSE
          CALL LCMSIX(IPHST,NAMP,ILCMUP)
          CALL LCMGET(IPHST,'DEPL-PARAM  ',BITH)
          CALL LCMSIX(IPHST,NAMP,ILCMDN)
        ENDIF
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6001) 'CELL TYPE',ICCT
          WRITE(IOUT,'(A6,1X,F8.3,2X,F8.3)') 'BURNUP',
     >       BITH(2),BURNUP(ICT,IBT)
        ENDIF 
        BURNUP(ICT,IBT) = BITH(2)
 20   CONTINUE
 10   CONTINUE
*----
* Store burnup record in MAP data structure
*----
      CALL LCMPUT(IPMAP,'BURN-DEB',NBUN*NCHA,2,BURNUP)
*----
*  Return
*----
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(' ****** OUTPUT FROM ',A6)
 6001 FORMAT(' Contents of ',A9,1X,I8) 
      END 
