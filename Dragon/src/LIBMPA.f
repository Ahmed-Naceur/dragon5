*DECK LIBMPA
      SUBROUTINE LIBMPA(NOR,JINI,WEIGHT,BASEPT,DEMP,SP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute a set of partial base points preserving NOR partial moments
* of a function using the modified Ribon approach.
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
* NOR     number of partial moments to preserve.
* JINI    minimum order of the partial moment we want to preserve.
*         We must have 1-NOR <= JINI <= 0 (order 0 moment is always
*         preserved).
* WEIGHT  weights of the probability table.
* BASEPT  base points of the probability table.
* DEMP    partial moments.
*
*Parameters: output
* SP      base points for the partial cross section.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NOR,JINI
      DOUBLE PRECISION DEMP(JINI:NOR+JINI-1)
      REAL WEIGHT(NOR),BASEPT(NOR),SP(NOR)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXNOR=20)
      DOUBLE PRECISION DDA(0:MAXNOR),DD,DSIGX
*
      IF(NOR.GT.MAXNOR) CALL XABORT('LIBMPA: STORAGE OVERFLOW.')
      IF(NOR.LE.0) CALL XABORT('LIBMPA: NEGATIVE OR ZERO VALUE OF NOR.')
      IF((1-NOR.GT.JINI).OR.(JINI.GT.0)) CALL XABORT('LIBMPA: INCONSIST'
     1 //'ENT VALUE OF JINI.')
*
      DO 50 I=1,NOR
      DSIGX=DBLE(BASEPT(I))
      DDA(0)=1.0D0
      J0=0
      DO 20 J=1,NOR
      IF(J.EQ.I) GO TO 20
      J0=J0+1
      DDA(J0)=DDA(J0-1)
      DO 10 K=1,J0-1
      DDA(J0-K)=DDA(J0-K-1)-DDA(J0-K)*DBLE(BASEPT(J))
   10 CONTINUE
      DDA(0)=-DDA(0)*DBLE(BASEPT(J))
   20 CONTINUE
      DD=0.0D0
      DO 30 J=0,NOR-1
      DD=DD+DDA(J)*DEMP(J+JINI)
   30 CONTINUE
      DO 40 J=1,NOR
      IF(J.NE.I) DD=DD/(DBLE(BASEPT(J))-DSIGX)
   40 CONTINUE
      SP(I)=REAL(((-1.0D0)**(NOR-1))*DD*DSIGX**(-JINI))/WEIGHT(I)
   50 CONTINUE
      RETURN
      END
