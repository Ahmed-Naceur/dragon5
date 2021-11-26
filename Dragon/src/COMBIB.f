*DECK COMBIB
      SUBROUTINE COMBIB(IPLB1,IPLB2,TYPE,IMILI,HBIB,HISO,MAXISO,VALPAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover a global parameter or a local variable from a microlib object.
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
*Parameters: input
* IPLB1   pointer to the first microlib object.
* IPLB2   pointer to the second (optional) microlib object.
* TYPE    ='temp' or 'conc'.
* IMILI   get the value in mixture imili.
* HBIB    character*12 name of the microlib.
* HISO    character*8 name of the isotope.
* MAXISO  allocated storage for isotopes.
*
*Parameters: output
* VALPAR  global parameter or local variable.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLB1,IPLB2
      INTEGER IMILI,MAXISO
      REAL VALPAR
      CHARACTER TYPE*4,HBIB*(*),HISO*(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) IPLIB
      CHARACTER NAMLCM*12,NAMMY*12,TEXT8*8
      INTEGER ISTATE(NSTATE)
      LOGICAL EMPTY,LCM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DEN,TN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISONAM(3,MAXISO),MIX(MAXISO))
      ALLOCATE(DEN(MAXISO),TN(MAXISO))
*
      IPLIB=C_NULL_PTR
      CALL LCMINF(IPLB1,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
      IF(NAMLCM.EQ.HBIB) THEN
         IPLIB=IPLB1
      ELSE IF(C_ASSOCIATED(IPLB2)) THEN
         CALL LCMINF(IPLB2,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
         IF(NAMLCM.EQ.HBIB) IPLIB=IPLB2
      ENDIF
      IF(.NOT.C_ASSOCIATED(IPLIB)) THEN
         NAMLCM=HBIB
         CALL XABORT('COMBIB: UNABLE TO FIND A MICROLIB NAMED '//
     1   NAMLCM//'.')
      ENDIF
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NBISOT=ISTATE(2)
      IF(NBISOT.GT.MAXISO) CALL XABORT('COMBIB: MAXISO OVERFLOW.')
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DEN)
      CALL LCMGET(IPLIB,'ISOTOPESTEMP',TN)
      IF(TYPE.EQ.'TEMP') THEN
         VALPAR=99999.0
         DO 10 I=1,NBISOT
         IF(MIX(I).EQ.IMILI) VALPAR=MIN(VALPAR,TN(I))
   10    CONTINUE
         IF(VALPAR.EQ.99999.0) CALL XABORT('COMBIB: UNABLE TO FIND A'//
     1   ' TEMP-TYPE PARAMETER OR LOCAL VARIABLE.')
      ELSE IF(TYPE.EQ.'CONC') THEN
         DO 20 I=1,NBISOT
         IF(MIX(I).EQ.IMILI) THEN
            WRITE(TEXT8,'(2A4)') (ISONAM(I0,I),I0=1,2)
            IF(TEXT8.EQ.HISO) THEN
               VALPAR=DEN(I)
               GO TO 30
            ENDIF
         ENDIF
   20    CONTINUE
         VALPAR=0.0
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   30 DEALLOCATE(TN,DEN)
      DEALLOCATE(MIX,ISONAM)
      RETURN
      END
