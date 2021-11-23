*DECK ALPRTB
      SUBROUTINE ALPRTB(NOR,IINI,DEMT,IER,WEIGHT,BASEPT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute a probability table preserving 2*NOR moments of a function
* using the modified Ribon approach.
*
*Copyright:
* Copyright (C) 1993 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NOR     the number of moments to preserve is 2*NOR.
* IINI    minimum order of the moment we want to preserve. we must
*         have 2-2*NOR <= IINI <= 0 (order 0 and 1 moments are always
*         preserved).
* DEMT    moments.
*
*Parameters: output
* IER     error flag (=0/=1 success/failure of the algorithm).
* WEIGHT  weights of the probability table.
* BASEPT  base points of the probability table.
*
*-----------------------------------------------------------------------
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NOR,IINI,IER
      DOUBLE PRECISION DEMT(IINI:2*NOR+IINI-1)
      REAL WEIGHT(NOR),BASEPT(NOR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXNOR=20)
      DOUBLE PRECISION DS(MAXNOR+1,MAXNOR+1),DDA(0:MAXNOR),DD,DSIGX
      COMPLEX*16 ROOTS(MAXNOR),CCC,DCC,XCC
      COMPLEX CGAR
      LOGICAL LFAIL
*
      IF(NOR.GT.MAXNOR) CALL XABORT('ALPRTB: STORAGE OVERFLOW.')
      IF(NOR.LE.0) CALL XABORT('ALPRTB: NEGATIVE OR ZERO VALUE OF NOR.')
      IF((2-2*NOR.GT.IINI).OR.(IINI.GT.0)) CALL XABORT('ALPRTB: INCONSI'
     1 //'STENT VALUE OF IINI.')
*
* BUILD THE MATRIX.
      DO 15 IOR=1,NOR
      DS(IOR,NOR+1)=-DEMT(NOR+IOR+IINI-1)
      DO 10 JOR=1,IOR
      DS(IOR,JOR)=DEMT(IOR+JOR+IINI-2)
      DS(JOR,IOR)=DEMT(IOR+JOR+IINI-2)
   10 CONTINUE
   15 CONTINUE
*
* L-D-L(T) FACTORIZATION OF THE MATRIX.
      DO 40 I=1,NOR
      DO 30 J=1,I-1
      DS(J,I)=DS(I,J)
      DO 20 K=1,J-1
      DS(J,I)=DS(J,I)-DS(K,I)*DS(J,K)
   20 CONTINUE
      DS(I,J)=DS(J,I)*DS(J,J)
      DS(I,I)=DS(I,I)-DS(J,I)*DS(I,J)
   30 CONTINUE
      IF(DS(I,I).EQ.0.D0) THEN
         IER=1
         RETURN
      ENDIF
      DS(I,I)=1.D0/DS(I,I)
   40 CONTINUE
*
* SOLUTION OF THE FACTORIZED SYSTEM TO OBTAIN THE DENOMINATOR OF THE
* PADE APPROXIMATION.
      DO 55 I=1,NOR
      DO 50 K=1,I-1
      DS(I,NOR+1)=DS(I,NOR+1)-DS(I,K)*DS(K,NOR+1)
   50 CONTINUE
   55 CONTINUE
      DO 60 I=1,NOR
      DS(I,NOR+1)=DS(I,NOR+1)*DS(I,I)
   60 CONTINUE
      DO 71 I=NOR,1,-1
      DO 70 K=I+1,NOR
      DS(I,NOR+1)=DS(I,NOR+1)-DS(K,I)*DS(K,NOR+1)
   70 CONTINUE
   71 CONTINUE
      DS(NOR+1,NOR+1)=1.0D0
*
* COMPUTE THE BASE POINTS AS THE ROOTS OF THE DENOMINATOR.
      CALL ALROOT(DS(1,NOR+1),NOR,ROOTS,LFAIL)
      IF(LFAIL) CALL XABORT('ALPRTB: POLYNOMIAL ROOT FINDING FAILURE.')
      DO 80 I=1,NOR
*
*     NEWTON IMPROVEMENT OF THE ROOTS.
      CCC=0.0D0
      XCC=1.0D0
      DO 74 J=0,NOR
      CCC=CCC+DS(J+1,NOR+1)*XCC
      XCC=XCC*ROOTS(I)
   74 CONTINUE
      DCC=0.0D0
      XCC=1.0D0
      DO 75 J=1,NOR
      DCC=DCC+DS(J+1,NOR+1)*XCC*REAL(J)
      XCC=XCC*ROOTS(I)
   75 CONTINUE
      ROOTS(I)=ROOTS(I)-CCC/DCC
*
      CGAR=CMPLX(ROOTS(I))
      IF(ABS(AIMAG(CGAR)).GT.1.0E-4*ABS(REAL(CGAR))) THEN
         IER=1
         RETURN
      ELSE
         BASEPT(I)=REAL(CMPLX(ROOTS(I)))
      ENDIF
   80 CONTINUE
*
* COMPUTE THE WEIGHTS.
      DO 130 I=1,NOR
      DSIGX=DBLE(ROOTS(I))
      DDA(0)=1.0D0
      J0=0
      DO 100 J=1,NOR
      IF(J.EQ.I) GO TO 100
      J0=J0+1
      DDA(J0)=DDA(J0-1)
      DO 90 K=1,J0-1
      DDA(J0-K)=DDA(J0-K-1)-DDA(J0-K)*DBLE(ROOTS(J))
   90 CONTINUE
      DDA(0)=-DDA(0)*DBLE(ROOTS(J))
  100 CONTINUE
      DD=0.0D0
      DO 110 J=0,NOR-1
      DD=DD+DDA(J)*DEMT((IINI-1)/2+J)
  110 CONTINUE
      DO 120 J=1,NOR
      IF(J.NE.I) DD=DD/(DBLE(ROOTS(J))-DSIGX)
  120 CONTINUE
      WEIGHT(I)=REAL(((-1.0D0)**(NOR-1))*DD*DSIGX**((1-IINI)/2))
  130 CONTINUE
*
* TEST THE CONSISTENCY OF THE SOLUTION.
      DO 140 I=1,NOR
      IF((WEIGHT(I).LE.0.0).OR.(BASEPT(I).LE.0.0)) THEN
         IER=1
         RETURN
      ENDIF
  140 CONTINUE
      IER=0
      RETURN
      END
