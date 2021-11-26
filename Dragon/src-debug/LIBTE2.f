*DECK LIBTE2
      SUBROUTINE LIBTE2 (NGRO,NSUBM,TMIX,SMIX,TN,SN,TERP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the Lagrange interpolation factors (TERP) for temperature and
* dilution interpolation of cross sections. TRANSX 2.0 algorithm.
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
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6)
      DOUBLE PRECISION CHECK
*
      IF(NSUBM.EQ.1) THEN
         DO 10 I=1,NGRO
         TERP(I,1)=1.0
   10    CONTINUE
         RETURN
      ENDIF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NSUBM))
*
      DO 115 ISUBM=1,NSUBM
      TMAT=TMIX(ISUBM)
      SMAT=SMIX(ISUBM)
*----
*  FIND TEMPERATURES VALUES
*----
      NTEMP=0
      DO 20 JSUBM=1,NSUBM
      IF(SMIX(JSUBM).LT.0.9E10) GO TO 20
      NTEMP=NTEMP+1
      WORK(NTEMP)=TMIX(JSUBM)
   20 CONTINUE
      TERPT=0.0
*----
*  LIMIT TEMPERATURE INTERPOLATION TO LINEAR IS TN A GRID TEMPERATURE?
*----
      DO 30 ITMP=1,NTEMP
      TT=WORK(ITMP)
      IF(ABS(TN-TT).LT.1.E-3*TN+1.E-3) THEN
         IF(ABS(TN-TMAT).LT.1.E-3*TMAT+1.E-3) TERPT=1.0
         GO TO 70
      ENDIF
   30 CONTINUE
*----
*  IF TEMP OUT OF RANGE USE ENDPOINTS
*----
      IF((TN.GT.WORK(NTEMP).AND.ABS(TMAT-WORK(NTEMP)).LT.1.E-3*TMAT
     > +1.E-3).OR.(TN.LT.WORK(1).AND.ABS(TMAT-WORK(1)).LT.
     > 1.E-3*TMAT+1.E-3)) THEN
         TERPT=1.0
         GO TO 70
      ENDIF
*----
*  FIND BRACKETING TEMPS
*----
      IF(NTEMP.EQ.1) THEN
         TERPT=1.0
         GO TO 70
      ENDIF
      DO 40 ITMP=1,NTEMP-1
      IF(WORK(ITMP).LT.TN.AND.WORK(ITMP+1).GT.TN) THEN
         ILOW=ITMP
         IHIGH=ITMP+1
         IF(ABS(TMAT-WORK(ITMP)).LT.1.E-3*TMAT+1.E-3.OR.ABS(TMAT-
     >    WORK(ITMP+1)).LT.1.E-3*TMAT+1.E-3) GO TO 50
      ENDIF
   40 CONTINUE
      GO TO 70
   50 TERPT=1.0
      DO 60 ITMP=ILOW,IHIGH
      TT=WORK(ITMP)
      IF(ABS(TMAT-TT).LT.1.E-3*TMAT+1.E-3) GO TO 60
      TERPT=TERPT*(TN-TT)/(TMAT-TT)
      IF(ABS(TERPT).LT.1.E-3) GO TO 70
   60 CONTINUE
*
*----
*  FIND SIGMA-ZERO VALUES
*----
   70 NTEMP=0
      NSIGZ=0
      DO 80 JSUBM=1,NSUBM
      IF(SMIX(JSUBM).GE.0.9E10) NTEMP=NTEMP+1
      IF(NTEMP.GT.1) GO TO 80
      NSIGZ=NSIGZ+1
      WORK(NSIGZ)=SMIX(JSUBM)
   80 CONTINUE
*----
*  FIND TERP FACTOR FOR SIGMA-ZERO
*----
      DO 110 JJ=1,NGRO
      TERPS=0.0
      IF((SN(JJ).GE.WORK(1)).OR.(NSIGZ.EQ.1)) THEN
         IF(SMAT.EQ.WORK(1)) TERPS=1.0
      ELSE IF(SN(JJ).LE.WORK(NSIGZ)) THEN
         IF(SMAT.EQ.WORK(NSIGZ)) TERPS=1.0
      ELSE IF((SN(JJ).GT.WORK(2)).OR.(NSIGZ.EQ.2)) THEN
         IF(SMAT.EQ.WORK(1)) TERPS=1.0-WORK(2)/SN(JJ)
         IF(SMAT.EQ.WORK(2)) TERPS=WORK(2)/SN(JJ)
      ELSE
         DO 90 I=2,NSIGZ-1
         IF(SN(JJ).LT.WORK(I+1)) GO TO 90
         IF(SMAT.EQ.WORK(I+1)) TERPS=(ALOG10(WORK(I))-ALOG10(SN(JJ)))
     X    /(ALOG10(WORK(I))-ALOG10(WORK(I+1)))
         IF(SMAT.EQ.WORK(I)) TERPS=(ALOG10(SN(JJ))-ALOG10(WORK(I+1)))
     X    /(ALOG10(WORK(I))-ALOG10(WORK(I+1)))
         GO TO 100
   90    CONTINUE
      ENDIF
  100 TERP(JJ,ISUBM)=TERPT*TERPS
  110 CONTINUE
  115 CONTINUE
*----
*  CHECK FOR CONSISTENCY OF THE TERP FACTORS.
*----
      DO 130 JJ=1,NGRO
      CHECK=0.0D0
      DO 120 ISUBM=1,NSUBM
      CHECK=CHECK+TERP(JJ,ISUBM)
  120 CONTINUE
      IF(ABS(CHECK-1.0D0).GT.5.0D-3) THEN
         WRITE (IOUT,200) JJ,CHECK,(TERP(JJ,ISUBM),ISUBM=1,NSUBM)
         CALL XABORT('LIBTE2: INTERPOLATION FAILURE.')
      ENDIF
  130 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
*
  200 FORMAT (/51H LIBTE2: INCONSISTENT LAGRANGE INTERPOLATION FACTOR,
     1 9H IN GROUP,I4,8H. CHECK=,1P,E13.3/9H FACTORS=,1P,9E13.3/
     2 (9X,1P,9E13.3))
      END
