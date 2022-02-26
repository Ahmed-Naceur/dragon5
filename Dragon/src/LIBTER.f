*DECK LIBTER
      SUBROUTINE LIBTER (NGRO,NSUBM,TMIX,SMIX,TN,SN,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the Lagrange interpolation factors (TERP) for temperature and
* dilution interpolation of cross sections. TRANSX CTR algorithm.
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
* NGRO    number of energy groups.
* NSUBM   number of submaterials (number of temperature/dilution
*         collocation points).
* TMIX    temperature of each submaterial in the library.
* SMIX    dilution of each submaterial in the library.
*         The submaterials are ordered by decreasing dilution and
*         then by increasing temperature.
* TN      temperature of the isotope.
* SN      dilution cross section in each energy group of the isotope.
*         A value of 1.0E10 is used for infinite dilution.
*
*Parameters: output
* TERP    Lagrange interpolation factors.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRO,NSUBM
      REAL TMIX(NSUBM),SMIX(NSUBM),TN,SN(NGRO),TERP(NGRO,NSUBM)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6)
      DOUBLE PRECISION CHECK
*----
*  TEMPERATURE AND BACKGROUND LOOP
*----
      IF (NSUBM.EQ.1) THEN
         DO 10 I=1,NGRO
         TERP(I,1)=1.0
   10    CONTINUE
         RETURN
      ENDIF
      BREAK=0.0
      SMIN=1.0E10
      DO 20 ISUBM=1,NSUBM
      ST=SMIX(ISUBM)
      IF ((ST.LT.1.0E8).AND.(ST.GT.BREAK)) BREAK=ST
      IF (ST.LT.SMIN) SMIN=ST
   20 CONTINUE
      DO 70 ISUBM=1,NSUBM
      DO 30 JJ=1,NGRO
      TERP(JJ,ISUBM)=0.0
   30 CONTINUE
      TMAT=TMIX(ISUBM)
      SMAT=SMIX(ISUBM)
*----
*  COMPUTE TERP FACTORS
*----
      TERPT=1.0
      DO 40 ISM=1,NSUBM
      TT=TMIX(ISM)
      ST=SMIX(ISM)
      IF (ST.LE.0.99E10) GO TO 40
      IF (ABS(TMAT-TT).LT.1.0E-5*TMAT+1.0E-5) GO TO 40
      TERPT=TERPT*(TN-TT)/(TMAT-TT)
      IF (ABS(TERPT).LT.1.0E-5) GO TO 70
   40 CONTINUE
      DO 60 JJ=1,NGRO
      TERPS=TERPT
      IF ((SN(JJ).LT.SMIN).AND.(SMAT.GT.1.01*SMIN)) GO TO 60
      IF ((SN(JJ).LT.BREAK).AND.(SMAT.GT.1.01*BREAK)) GO TO 60
      IF ((SN(JJ).GE.BREAK).AND.(SMAT.LT.0.99*BREAK)) GO TO 60
      TLAST=-1.0
      DO 50 ISM=1,NSUBM
      TT=TMIX(ISM)
      ST=SMIX(ISM)
      IF (TLAST.LT.0.) TLAST=TT
      IF (TT.NE.TLAST) GO TO 50
      IF (ABS(SMAT-ST).LT.1.0E-5*SMAT) GO TO 50
      IF ((SN(JJ).GE.SMIN).AND.(SN(JJ).LT.BREAK)) THEN
         IF (ST.GT.1.01*BREAK) GO TO 50
         TERPS=TERPS*LOG(SN(JJ)/ST)/LOG(SMAT/ST)
      ELSE IF ((SN(JJ).GE.SMIN).AND.(SN(JJ).GE.BREAK)) THEN
         IF (ST.LT.0.99*BREAK) GO TO 50
         TERPS=TERPS*((ST/SN(JJ))-1.)/((ST/SMAT)-1.)
      ELSE
         IF (ST.GT.1.01*SMIN) GO TO 50
         TERPS=TERPS*(SN(JJ)**2-ST**2)/(SMAT**2-ST**2)
      ENDIF
      IF (ABS(TERPS).LE.1.0E-5) GO TO 60
   50 CONTINUE
      TERP(JJ,ISUBM)=TERPS
   60 CONTINUE
   70 CONTINUE
*----
*  CHECK FOR CONSISTENCY OF THE TERP FACTORS.
*----
      DO 90 JJ=1,NGRO
      CHECK=0.0D0
      DO 80 ISUBM=1,NSUBM
      CHECK=CHECK+TERP(JJ,ISUBM)
   80 CONTINUE
      IF (ABS(CHECK-1.0D0).GT.5.0D-3) THEN
         WRITE (IOUT,100) JJ,CHECK,(TERP(JJ,ISUBM),ISUBM=1,NSUBM)
         CALL XABORT('LIBTER: INTERPOLATION FAILURE.')
      ENDIF
   90 CONTINUE
      RETURN
  100 FORMAT (/51H LIBTER: INCONSISTENT LAGRANGE INTERPOLATION FACTOR,
     1 9H IN GROUP,I4,8H. CHECK=,1P,E13.3/9H FACTORS=,1P,9E13.3/
     2 (9X,1P,9E13.3))
      END
