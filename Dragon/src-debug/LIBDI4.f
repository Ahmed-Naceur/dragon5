*DECK LIBDI4
      SUBROUTINE LIBDI4 (MAXDIL,NAMFIL,HSHI,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in Apolib-1 format.
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
* MAXDIL  maximum number of dilutions.
* NAMFIL  name of the Apolib file.
* HSHI    library name of the self-shielding data.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXDIL,NDIL
      CHARACTER NAMFIL*(*),HSHI*12
      REAL DILUT(MAXDIL)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXIT=1000)
      CHARACTER FORM*4,HSMG*131
      INTEGER IT(MAXIT),NTETA(3)
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(MAXDIL))
*
      NIN=KDROPN(NAMFIL,2,2,0)
      IF(NIN.LE.0) THEN
         WRITE (HSMG,'(35HLIBDI4: UNABLE TO OPEN LIBRARY FILE,1X,A16,
     1   6H. NIN=,I4,1H.)') NAMFIL,NIN
         CALL XABORT(HSMG)
      ENDIF
      I=INDEX(HSHI,' ')
      IF(HSHI.EQ.' ') THEN
         NISBEF=0
      ELSE IF(I.EQ.0) THEN
         READ(HSHI,'(I8)') NISBEF
      ELSE
         WRITE(FORM,'(2H(I,I1,1H))') I-1
         READ(HSHI,FORM) NISBEF
      ENDIF
   10 READ(NIN) INDLOR,NR,NIT,(IT(I),I=1,NIT),(DUMMY,I=1,18)
      IF(NIT.GT.MAXIT) CALL XABORT('LIBDI4: INVALID MAXIT.')
      IF(INDLOR.EQ.9999) THEN
         CALL XABORT('LIBDI4: UNABLE TO FIND ISOTOPE '//HSHI//'.')
      ELSE IF(INDLOR.EQ.NISBEF) THEN
         NTYPE=0
         JTYSEC=0
         DO 20 IK=1,IT(4)
         IF(IT(IK+4).NE.JTYSEC) THEN
            NTYPE=NTYPE+1
            NTETA(NTYPE)=1
            JTYSEC=IT(IK+4)
         ELSE
            NTETA(NTYPE)=NTETA(NTYPE)+1
         ENDIF
   20    CONTINUE
         DO 55 I=1,NTYPE
         READ (NIN) TEMP,NSEI,(WORK(K),K=1,NSEI)
         IF(NSEI.GT.MAXDIL) CALL XABORT('LIBDI4: INVALID MAXDIL.')
         IF(I.EQ.1) THEN
            NDIL=NSEI
            DO 30 K=NSEI,1,-1
            IF(WORK(K).GE.1.0E10) THEN
               NDIL=NDIL-1
            ELSE
               DILUT(K)=WORK(K)
            ENDIF
   30       CONTINUE
            DILUT(NDIL+1)=1.0E10
         ELSE
            DO 40 K=NSEI,1,-1
            IF((WORK(K).LT.1.0E10).AND.(WORK(K).NE.DILUT(K))) THEN
               WRITE(HSMG,'(26HLIBDI4: INVALID DILUTION (,1P,E12.4,
     1         9H) ON TYPE,I2,11H REACTIONS.,E12.4,10H EXPECTED.)')
     2         WORK(K),I,DILUT(K)
               CALL XABORT(HSMG)
            ENDIF
   40       CONTINUE
         ENDIF
         DO 50 ITET=2,NTETA(I)
         READ(NIN)
   50    CONTINUE
   55    CONTINUE
      ELSE
         DO 60 K=1,NR
         READ(NIN)
   60    CONTINUE
         GO TO 10
      ENDIF
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE (HSMG,'(36HLIBDI4: UNABLE TO CLOSE LIBRARY FILE,1X,A16,
     1   1H.)') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
