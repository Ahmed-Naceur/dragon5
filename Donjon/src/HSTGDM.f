*DECK HSTGDM
      SUBROUTINE HSTGDM(IPRINT, NGLO,   NLOC,   NCHA,   NBUN  ,
     >                  BUNLEN, ITYRED, CARRED)
*
*----------
*
*Purpose:
* To read the editing level and general dimensioning parameters
* for the \dds{history} data structure.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal.
*
*Author(s): 
* G. Marleau
*
*Parameters: input/output
* IPRINT  print level.
* NGLO    number of global parameters.               
* NLOC    number of local parameters.                
* NCHA    number of fuel channels.                   
* NBUN    number of bundles per channel.            
* BUNLEN  length (cm) of a bundle.            
* ITYRED  type of the last variable read.                
* CARRED  last character string read.
*                                                   
*----------                                         
*                                                   
      USE GANLIB
      IMPLICIT         NONE                         
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          IPRINT,NGLO,NLOC,NCHA,NBUN
      REAL             BUNLEN
      INTEGER          ITYRED
      CHARACTER*12     CARRED
*----
*  LOCAL PARAMETERS
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='HSTGDM')
*----
*  INPUT VARIABLES 
*  Input data is of the form
*  [ EDIT iprint ]
*  [ DIMENSIONS 
*    [ GLOBAL  nglo  ]
*    [ LOCAL   nloc   ]
*    [ BUNDL   nbun  bunl ]
*    [ CHANNEL ncha ]
*----
      INTEGER          ITYPLU,INTLIR
      CHARACTER        CARLIR*12
      REAL             REALIR
      DOUBLE PRECISION DBLLIR
*----
*  Initialize output variables  variables
*----
      ITYPLU= 0
      CARLIR='            '
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
 101  CONTINUE
      IF(ITYPLU .EQ. 10) GO TO 105
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': Read error -- Character variable expected')
      IF(CARLIR .EQ. ';') THEN
        GO TO 105
      ELSE IF(CARLIR(1:4) .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) THEN
          IPRINT=1
          GO TO 101
        ENDIF
        IPRINT=INTLIR
        GO TO 100
      ELSE IF(CARLIR(1:4) .EQ. 'DIME') THEN
 110    CONTINUE
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': Read error -- Dimension type expected')
        IF(CARLIR(1:4) .EQ. 'GLOB') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >    ': Number of global parameters expected')
          NGLO=INTLIR
          GO TO 110
        ELSE IF(CARLIR(1:4) .EQ. 'LOCA') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >    ': Number of local parameters expected')
          NLOC=INTLIR
          GO TO 110
        ELSE IF(CARLIR(1:4) .EQ. 'BUND') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >    ': Number of bundles expected')
          NBUN=INTLIR
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >    ': Bundles length (cm) expected')
          IF(REALIR .GT. 0.0) BUNLEN=REALIR
          GO TO 110
        ELSE IF(CARLIR(1:4) .EQ. 'CHAN') THEN
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >    ': Number of channels expected')
          NCHA=INTLIR
          GO TO 110
        ELSE
          GO TO 105
        ENDIF
      ENDIF
 105  CONTINUE
      IF(NGLO .LT. 0) THEN
        NGLO=0
        WRITE(IOUT,8000) NAMSBR,'nglo'
      ENDIF
      IF(NLOC .LT. 0) THEN
        NLOC=0
        WRITE(IOUT,8000) NAMSBR,'nloc'
      ENDIF
      IF(NBUN .LT. 0) THEN
        NBUN=0
        WRITE(IOUT,8000) NAMSBR,'nbun'
      ENDIF
      IF(NCHA .LT. 0) THEN
        NCHA=0
        WRITE(IOUT,8000) NAMSBR,'ncha'
      ENDIF
      ITYRED=ITYPLU
      CARRED=CARLIR
*----
*  Format
*----
 8000 FORMAT(' ****** WARNING in ',A6,' ****** '/
     >       ' Problem  : ',A4,1X,' < 0'/
     >       ' Solution : assume this parameter is not read'/
     >       ' ******************************')
      RETURN
      END 
