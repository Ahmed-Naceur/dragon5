*DECK CFCGET
      SUBROUTINE CFCGET(TINFO,DBNAME,IPRINT,NBPARA,DBPARA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read CFC options.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* TINFO   title of database.
* DBNAME  database name.
* IPRINT  print level.
* NBPARA  number of parameters.
* DBPARA  database parameters: 
*         DBPARA( 1) = Nominal power level (PWR);
*         DBPARA( 2) = Nominal T cool (TCR);
*         DBPARA( 3) = Nominal T mode (TMR);
*         DBPARA( 4) = Nominal T fuel (TFR);
*         DBPARA( 5) = Nominal Density cool;
*         DBPARA( 6) = Nominal Density mode;
*         DBPARA( 7) = Nominal purity mode (XIR);
*         DBPARA( 8) = Perturbed T fuel 1 (TFU);
*         DBPARA( 9) = Perturbed T cool 1 (TCU);
*         DBPARA(10) = Perturbed P1 (PWUL);
*         DBPARA(11) = Perturbed P2 (PWDL);
*         DBPARA(12) = Perturbed P3 (PWU);
*         DBPARA(13) = Perturbed P4 (PWD);
*         DBPARA(14) = Perturbed P mode (XI).
*
*-----------------------------------------------------------------------
*
      USE         GANLIB
      IMPLICIT    NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER   TINFO*72,DBNAME*9
      INTEGER     IPRINT,NBPARA
      REAL        DBPARA(NBPARA)
*----
*  LOCAL PARAMETERS
*----
      INTEGER     IOUT
      CHARACTER   NAMSBR*6
      PARAMETER  (IOUT=6,NAMSBR='CFCGET')
*----
*  INPUT PARAMETERS
*----
      INTEGER     ITYPLU,INTLIR
      CHARACTER   CARLIR*12
      REAL        REALIR
      DOUBLE PRECISION DBLLIR
*----
*  INITIALIZE PARAMETERS
*----
      IPRINT=1
      CALL XDRSET(DBPARA,NBPARA,0.0)
      DBPARA(1)=615.0
      DBPARA(2)=560.66
      DBPARA(3)=345.66
      DBPARA(4)=941.29
      DBPARA(5)=1.08288
      DBPARA(6)=0.81212
      DBPARA(7)=0.99911
      DBPARA(8)=1241.29
      DBPARA(9)=660.66
      DBPARA(10)=878.57143
      DBPARA(11)=307.5
      DBPARA(12)=307.5
      DBPARA(13)=100.0
      DBPARA(14)=0.985
      DBPARA(15)=375.66
      DBPARA(16)=541.29
      DBPARA(17)=300.66
      DBPARA(18)=295.66
*----
*  READ TITLE
*----
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': CHARACTER DATA EXPECTED.')
      IF(CARLIR .EQ. 'EDIT') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 1) THEN
          IPRINT=INTLIR
        ELSE
          CALL XABORT(NAMSBR//': EDIT LEVEL EXPECTED.')
        ENDIF
        CALL REDGET (ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': CHARACTER DATA EXPECTED.')
      ENDIF
      IF(CARLIR .EQ. 'INFOR') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,TINFO,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': CHARACTER DATA EXPECTED.')
      ELSE
        CALL XABORT(NAMSBR//': INFOR KEY WORD EXPECTED.')
      ENDIF
*-----
*  READ THE RECORD NAME FOR THE DATABASE
*-----
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': CHARACTER DATA EXPECTED.')
      IF(CARLIR .EQ. 'DNAME') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': CHARACTER DATA EXPECTED.')
        DBNAME=CARLIR(1:9)
      ENDIF
*----
*  LOOP OVER PARAMETERS TO READ
*----
 100  CONTINUE
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >': CHARACTER DATA EXPECTED.')
      IF(CARLIR .EQ. ';') THEN
        GO TO 105
      ELSE IF(CARLIR .EQ. 'PWR') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal power EXPECTED.')
        DBPARA(1)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Power UP EXPECTED.')
        DBPARA(10)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Power DOWN 1 EXPECTED.')
        DBPARA(11)=REALIR
        DBPARA(12)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Power DOWN 2 EXPECTED.')
        DBPARA(13)=REALIR
      ELSE IF(CARLIR .EQ. 'TCOOL') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal coolant temperature expected.')
        DBPARA(2)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed up coolant temperature expected.')
        DBPARA(9)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed down coolant temperature expected.')
        DBPARA(17)=REALIR
      ELSE IF(CARLIR .EQ. 'TMODE') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal moderator temperature expected.')
        DBPARA(3)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed up moderator temperature expected.')
        DBPARA(15)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed down moderator temperature expected.')
        DBPARA(18)=REALIR
      ELSE IF(CARLIR .EQ. 'TFUEL') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal fuel temperature expected.')
        DBPARA(4)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed up fuel temperature expected.')
        DBPARA(8)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed down fuel temperature expected.')
        DBPARA(16)=REALIR
      ELSE IF(CARLIR .EQ. 'RHOM') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal moderator density expected.')
        DBPARA(5)=REALIR
      ELSE IF(CARLIR .EQ. 'RHOC') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal coolant density expected.')
        DBPARA(6)=REALIR
      ELSE IF(CARLIR .EQ. 'XIR') THEN
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Nominal Xenon.')
        DBPARA(7)=REALIR
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 2) CALL XABORT(NAMSBR//
     >  ': Perturbed Xenon.')
        DBPARA(14)=REALIR
      ENDIF
      GO TO 100
 105  CONTINUE
      RETURN
      END
