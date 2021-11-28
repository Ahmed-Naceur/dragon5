*DECK KTRDRV
      INTEGER FUNCTION KTRDRV(HMODUL,NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Code dependent operator driver for TRIVAC.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* HMODUL  name of the operator.
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Parameters: output
* KTRDRV  completion flag (=0: operator HMODUL exists; =1: does not
*         exists).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER HMODUL*(*),HENTRY(NENTRY)*12
      INTEGER IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      REAL TBEG,TEND
      DOUBLE PRECISION DMEMB,DMEMD
      LOGICAL :: TRIMOD
*
      KTRDRV=0
      TRIMOD=.TRUE.
      CALL KDRCPU(TBEG)
      CALL KDRMEM(DMEMB)
      IF(HMODUL.EQ.'BIVACA:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL BIVACA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'BIVACT:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL BIVACT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'FLUD:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL FLD(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'GEO:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL GEOD(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'MAC:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL MACD(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'TRIVAT:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL TRIVAT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'TRIVAA:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL TRIVAA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'OUT:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL OUT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'ERROR:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL ERROR(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'DELTA:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL DELTA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'GPTFLU:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL GPTFLU(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'INIKIN:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'D. SEKKI'
         CALL KININI(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'KINSOL:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'D. SEKKI'
         CALL KINSOL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'VAL:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'R. CHAMBON'
         CALL VAL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'NSST:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL NSST(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE IF(HMODUL.EQ.'NSSF:') THEN
         WRITE(6,'(//7H EXEC: ,A,4H BY ,A)') HMODUL,'A. HEBERT'
         CALL NSSF(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ELSE
         TRIMOD=.FALSE.
         KTRDRV=GANDRV(HMODUL,NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
      ENDIF
      IF(TRIMOD)THEN
         CALL KDRCPU(TEND)
         CALL KDRMEM(DMEMD)
         WRITE(6,5000) HMODUL,(TEND-TBEG),REAL(DMEMD-DMEMB)
      ENDIF
      RETURN
*
 5000 FORMAT('-->>MODULE ',A12,': TIME SPENT=',F13.3,' MEMORY USAGE=',
     1 1P,E10.3)
      END
