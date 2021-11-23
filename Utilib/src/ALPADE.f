*DECK ALPADE
      SUBROUTINE ALPADE(NORIN,X,Y,EPSRID,NOR,A,B,PREC,IER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the polynomial coefficients of a Pade approximation using an
* inverse differences collocation.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NORIN   2*NORIN+1 is the number of collocation points.
* X       abscissa of the collocation points.
* Y       ordinates of the collocation points.
* EPSRID  epsilon used in polynomial simplification.
*
*Parameters: output
* NOR     order of the polynomials.
* A       polynomial coefficients of the numerator of the Pade
*         approximation. a(0) is the constant term.
* B       polynomial coefficients of the denominator of the Pade
*         approximation. b(0) is the constant term.
*         DOUBLE PRECISION A(0:NOR),B(0:NOR)
* PREC    accuracy of the fit.
* IER     error flag (=0: no error; =1: negative pole removing).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NORIN,IER
      REAL X(0:2*NORIN),Y(0:2*NORIN),PREC
      DOUBLE PRECISION EPSRID,A(0:NORIN),B(0:NORIN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXNOR=10)
      DOUBLE PRECISION DC(0:MAXNOR),DENOM,DY(0:2*MAXNOR,0:2*MAXNOR),
     1 ERR1,ERR2,GARINF,SCC,SDD
      COMPLEX*16 DDA(0:MAXNOR-2),DDB(0:MAXNOR-2),SIGX0(MAXNOR),
     1 SIGXW(MAXNOR)
      LOGICAL LINF,LFAIL
*
      IER=0
      IF(NORIN.GT.MAXNOR) CALL XABORT('ALPADE: INSUFFICIENT MAXNOR.')
      LINF=X(2*NORIN).GE.1.0E10
      JMAX=2*NORIN
      ERR2=DBLE(Y(1))
      GARINF=0.0D0
      DO 15 J=0,2*NORIN
      IF(X(J).LE.0.0) CALL XABORT('ALPADE: ZERO OR NEGATIVE ABSCISSA.')
      IF(Y(J).LE.0.0) CALL XABORT('ALPADE: ZERO OR NEGATIVE ORDINATE.')
      ERR1=0.0D0
      DO 10 I=J,2*NORIN
      IF(J.EQ.0) THEN
         DY(I,J)=DBLE(Y(I))
      ELSE IF(LINF.AND.(MOD(J,2).EQ.1).AND.(I.EQ.2*NORIN)) THEN
         DENOM=DY(I,J-1)-DY(J-1,J-1)
         IF(DENOM.EQ.0.0) CALL XABORT('ALPADE: ALGORITHM FAILURE(1).')
         DY(I,J)=1.0D0/DENOM
      ELSE IF(LINF.AND.(I.EQ.2*NORIN)) THEN
         DENOM=DY(I,J-1)
         IF(DENOM.EQ.0.0) CALL XABORT('ALPADE: ALGORITHM FAILURE(2).')
         DY(I,J)=1.0D0/DENOM
      ELSE
         DENOM=DY(I,J-1)-DY(J-1,J-1)
         IF(DENOM.EQ.0.0) CALL XABORT('ALPADE: ALGORITHM FAILURE(3).')
         DY(I,J)=(DBLE(X(I))-DBLE(X(J-1)))/DENOM
      ENDIF
      ERR1=MAX(ERR1,ABS(DY(I,J)-DY(J,J)))
   10 CONTINUE
      IF(MOD(J,2).EQ.0) GARINF=GARINF+DY(J,J)
      IF(LINF.AND.(ERR1.LE.1.0D-6*ERR2).AND.(ABS(GARINF-Y(2*NORIN)).LE.
     1 1.0D-5*ABS(GARINF))) THEN
         JMAX=J
         GO TO 20
      ENDIF
      ERR2=ERR1
   15 CONTINUE
*
   20 IF(MOD(JMAX,2).NE.0) CALL XABORT('ALPADE: ALGORITHM FAILURE(4).')
      N=0
      MM=JMAX-1
      A(0)=DY(JMAX,JMAX)
      B(0)=1.0D0
      NOR=JMAX/2
      DO 60 K=1,NOR
      DC(0)=0.0D0
      DO 30 I=0,N
      DC(I+1)=B(I)
      DC(I)=DC(I)-B(I)*X(MM)+A(I)*DY(MM,MM)
   30 CONTINUE
      MM=MM-1
      B(0)=0.0D0
      DO 40 I=0,N
      B(I+1)=A(I)
      B(I)=B(I)-A(I)*X(MM)+DC(I)*DY(MM,MM)
   40 CONTINUE
      B(N+1)=B(I)+DC(N+1)*DY(MM,MM)
      DO 50 I=0,N+1
      A(I)=B(I)
      B(I)=DC(I)
   50 CONTINUE
      MM=MM-1
      N=N+1
   60 CONTINUE
*
* POLYNOMIAL SIMPLIFICATION.
      DDA(0)=A(NOR)
      DDB(0)=B(NOR)
      IF(NOR.EQ.0) GO TO 120
      CALL ALROOT(A,NOR,SIGX0,LFAIL)
      IF(LFAIL) CALL XABORT('ALPADE: POLYNOMIAL ROOT FINDING FAILURE.')
      CALL ALROOT(B,NOR,SIGXW,LFAIL)
      IF(LFAIL) CALL XABORT('ALPADE: POLYNOMIAL ROOT FINDING FAILURE.')
      IJNOR=1
   70 XXX=ABS(REAL(CMPLX(SIGXW(IJNOR)))-REAL(CMPLX(SIGX0(IJNOR))))
      IF(XXX.LT.EPSRID*ABS(SIGXW(IJNOR))) THEN
         NOR=NOR-1
         DO 80 I=IJNOR,NOR
         SIGX0(I)=SIGX0(I+1)
         SIGXW(I)=SIGXW(I+1)
   80    CONTINUE
      ELSE IF((DBLE(SIGXW(IJNOR)).GT.0.).AND.(DIMAG(SIGXW(IJNOR)).EQ.0.)
     1 ) THEN
         IER=1
         NOR=NOR-1
         DO 90 I=IJNOR,NOR
         SIGX0(I)=SIGX0(I+1)
         SIGXW(I)=SIGXW(I+1)
   90    CONTINUE
      ELSE
         IJNOR=IJNOR+1
      ENDIF
      IF(IJNOR.LE.NOR) GO TO 70
      IF(NOR.LT.0) CALL XABORT('ALPADE: ALGORITHM FAILURE(5).')
      DO 110 I=1,NOR
      DDA(I)=DDA(I-1)
      DDB(I)=DDB(I-1)
      DO 100 J=I-1,1,-1
      DDA(J)=DDA(J-1)-DDA(J)*SIGX0(I)
      DDB(J)=DDB(J-1)-DDB(J)*SIGXW(I)
  100 CONTINUE
      DDA(0)=-DDA(0)*SIGX0(I)
      DDB(0)=-DDB(0)*SIGXW(I)
  110 CONTINUE
  120 DENOM=DBLE(DDB(NOR))
      DO 130 I=0,NOR
      A(I)=DBLE(DDA(I))/DENOM
      B(I)=DBLE(DDB(I))/DENOM
  130 CONTINUE
*
* TEST THE ACCURACY OF THE PADE APPROXIMATION.
      PREC=0.0
      PREC1=0.0
      DO 150 I=0,2*NORIN
      SCC=A(NOR)
      SDD=B(NOR)
      IF(X(I).LT.1.0E10) THEN
         DO 140 INOR=NOR-1,0,-1
         SCC=A(INOR)+SCC*X(I)
         SDD=B(INOR)+SDD*X(I)
  140    CONTINUE
      ENDIF
      PREC=MAX(PREC,ABS(REAL(SCC/SDD)/Y(I)-1.0))
      PREC1=MAX(PREC1,ABS(Y(2*NORIN)/Y(I)-1.0))
  150 CONTINUE
      IF((IER.NE.0).AND.(PREC.GT.0.99*PREC1)) THEN
*        USE A UNIFORM REPRESENTATION.
         NOR=0
         A(0)=DBLE(Y(2*NORIN))
         B(0)=1.0D0
         PREC=PREC1
      ENDIF
      RETURN
      END
