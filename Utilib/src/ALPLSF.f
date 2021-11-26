*DECK ALPLSF
      SUBROUTINE ALPLSF(IMETH,N,X,Y,EPSRID,LREAL,NOR,A,B,PREC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the polynomial coefficients of a Pade approximation using a
* direct least square procedure.
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
* IMETH   type of algorithm (=1: use QR; =2: use SVD; =3: use NNLS).
* N       number of collocation points.
* X       abscissa of the collocation points.
* Y       ordinates of the collocation points.
* EPSRID  epsilon used in polynomial simplification.
* LREAL   selection flag (=.true. to get rid of complex roots).
*
*Parameters: output
* NOR     order of the polynomials.
* A       polynomial coefficients of the numerator of the Pade
*         approximation. A(0) is the constant term.
* B       polynomial coefficients of the denominator of the Pade
*         approximation. B(0) is the constant term.
*         DOUBLE PRECISION A(0:NOR),B(0:NOR)
* PREC    accuracy of the fit.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMETH,N,NOR
      REAL X(N),Y(N),PREC
      DOUBLE PRECISION EPSRID,A(0:(N-1)/2),B(0:(N-1)/2)
      LOGICAL LREAL
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXNOR=10,MAXPTS=99)
      DOUBLE PRECISION BB(MAXPTS),GAR,PARAM(MAXPTS),AGAR(0:MAXNOR),
     1 BGAR(0:MAXNOR),CGAR(0:MAXNOR+1),W(MAXPTS),RV1(MAXPTS),SGN,
     2 GAROLD,YAPPR,RNORM,RMAX
      COMPLEX*16 SIGX0(MAXNOR+1),SIGXW(MAXNOR+1),DDAGAR(0:MAXNOR),
     1 DDBGAR(0:MAXNOR),WEIGH(MAXNOR+1),CC,DD,CCC,XCC,DCC
      LOGICAL LFAIL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: AA,V
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(AA(MAXPTS,MAXPTS),V(MAXPTS,MAXPTS))
*
      IF(N.GT.MAXPTS) CALL XABORT('ALPLSF: INSUFFICIENT MAXPTS.')
*
* NOR=0 CASE.
      GAROLD=0.0
      NOR=0
      PREC=0.0
      A(0)=DBLE(Y(N))
      B(0)=1.0D0
      DO 10 I=1,N-1
      PREC=MAX(PREC,ABS(Y(N)/Y(I)-1.0))
   10 CONTINUE
*
      RMAX=1.0D10
      DO 210 IINOR=1,MIN((N-1)/2,MAXNOR)
      INOR=IINOR
      IF(X(N).GE.0.99E10) THEN
         DO 41 I=1,N-1
         GAR=1.0D0
         IOF=0
         DO 20 J=0,INOR-1
         IOF=IOF+1
         AA(I,IOF)=GAR
         GAR=GAR*X(I)
   20    CONTINUE
         GAR=1.0D0
         DO 30 J=0,INOR-1
         IOF=IOF+1
         GAROLD=GAR
         AA(I,IOF)=-GAR*Y(I)
         GAR=GAR*X(I)
   30    CONTINUE
         BB(I)=(Y(I)-Y(N))*X(I)
         DO 40 J=1,2*INOR
         AA(I,J)=AA(I,J)/GAROLD
   40    CONTINUE
   41    CONTINUE
         IF(IMETH.EQ.1) THEN
            CALL ALST2F(MAXPTS,N-1,2*INOR,AA,RV1)
            CALL ALST2S(MAXPTS,N-1,2*INOR,AA,RV1,BB,PARAM)
         ELSE IF(IMETH.EQ.2) THEN
            CALL ALSVDF(AA,N-1,2*INOR,MAXPTS,MAXPTS,W,V,RV1)
            DO 45 J=1,2*INOR
            IF(W(J).EQ.0.0D0) CALL XABORT('ALPLSF: SVD FAILURE(1).')
   45       CONTINUE
            CALL ALSVDS(AA,W,V,N-1,2*INOR,MAXPTS,MAXPTS,BB,PARAM,RV1)
         ELSE IF(IMETH.EQ.3) THEN
            CALL ALNNLS(AA,N-1,2*INOR,MAXPTS,MAXPTS,BB,PARAM,RNORM,MODE)
            IF(MODE.NE.1) CALL XABORT('ALPLSF: NNLS FAILURE(1).')
            IF((INOR.GT.1).AND.(RNORM.GE.0.95D0*RMAX)) GO TO 210
            RMAX=RNORM
         ENDIF
         DO 50 I=0,INOR-1
         AGAR(I)=PARAM(I+1)
         BGAR(I)=PARAM(INOR+1+I)
   50    CONTINUE
         AGAR(INOR)=Y(N)
         BGAR(INOR)=1.0D0
      ELSE
         DO 81 I=1,N
         GAR=1.0D0
         IOF=0
         DO 60 J=0,INOR
         IOF=IOF+1
         AA(I,IOF)=GAR
         GAR=GAR*X(I)
   60    CONTINUE
         GAR=1.0D0
         DO 70 J=0,INOR-1
         IOF=IOF+1
         GAROLD=GAR
         AA(I,IOF)=-GAR*Y(I)
         GAR=GAR*X(I)
   70    CONTINUE
         BB(I)=Y(I)*X(I)
         DO 80 J=1,2*INOR+1
         AA(I,J)=AA(I,J)/GAROLD
   80    CONTINUE
   81    CONTINUE
         IF(IMETH.EQ.1) THEN
            CALL ALST2F(MAXPTS,N,2*INOR+1,AA,RV1)
            CALL ALST2S(MAXPTS,N,2*INOR+1,AA,RV1,BB,PARAM)
         ELSE IF(IMETH.EQ.2) THEN
            CALL ALSVDF(AA,N,2*INOR+1,MAXPTS,MAXPTS,W,V,RV1)
            DO 85 J=1,2*INOR
            IF(W(J).EQ.0.0D0) CALL XABORT('ALPLSF: SVD FAILURE(2).')
   85       CONTINUE
            CALL ALSVDS(AA,W,V,N,2*INOR+1,MAXPTS,MAXPTS,BB,PARAM,RV1)
         ELSE IF(IMETH.EQ.3) THEN
            CALL ALNNLS(AA,N,2*INOR+1,MAXPTS,MAXPTS,BB,PARAM,RNORM,MODE)
            IF(MODE.NE.1) CALL XABORT('ALPLSF: NNLS FAILURE(2).')
            IF((INOR.GT.1).AND.(RNORM.GE.0.95D0*RMAX)) GO TO 210
            RMAX=RNORM
         ENDIF
         DO 90 I=0,INOR
         AGAR(I)=PARAM(I+1)
         IF(I.EQ.INOR) THEN
            BGAR(I)=1.0D0
         ELSE
            BGAR(I)=PARAM(INOR+2+I)
         ENDIF
   90    CONTINUE
      ENDIF
*
* POLYNOMIAL SIMPLIFICATION.
      DDAGAR(0)=AGAR(INOR)
      DDBGAR(0)=BGAR(INOR)
      CALL ALROOT(AGAR,INOR,SIGX0,LFAIL)
      IF(LFAIL) GO TO 210
      CALL ALROOT(BGAR,INOR,SIGXW,LFAIL)
      IF(LFAIL) GO TO 210
      IJINOR=1
   95 XXX=REAL(ABS(DBLE(SIGXW(IJINOR))-DBLE(SIGX0(IJINOR))))
      IF(XXX.LT.EPSRID*ABS(DBLE(SIGXW(IJINOR)))) THEN
         INOR=INOR-1
         DO 100 I=IJINOR,INOR
         SIGX0(I)=SIGX0(I+1)
         SIGXW(I)=SIGXW(I+1)
  100    CONTINUE
      ELSE IF((DBLE(SIGXW(IJINOR)).GT.EPSRID).AND.
     >        (DIMAG(SIGXW(IJINOR)).EQ.0.0).AND.
     >        (IMETH.EQ.3)) THEN
         CALL XABORT('ALPLSF: NNLS FAILURE(3).')
      ELSE IF((DBLE(SIGXW(IJINOR)).GT.0.1*EPSRID).AND.
     >        (DIMAG(SIGXW(IJINOR)).EQ.0.0)) THEN
         GO TO 210
      ELSE
         IJINOR=IJINOR+1
      ENDIF
      IF(IJINOR.LE.INOR) GO TO 95
      IF(INOR.LT.0) CALL XABORT('ALPLSF: ALGORITHM FAILURE.')
      DO 120 I=1,INOR
      DDAGAR(I)=DDAGAR(I-1)
      DDBGAR(I)=DDBGAR(I-1)
      DO 110 J=I-1,1,-1
      DDAGAR(J)=DDAGAR(J-1)-DDAGAR(J)*SIGX0(I)
      DDBGAR(J)=DDBGAR(J-1)-DDBGAR(J)*SIGXW(I)
  110 CONTINUE
      DDAGAR(0)=-DDAGAR(0)*SIGX0(I)
      DDBGAR(0)=-DDBGAR(0)*SIGXW(I)
  120 CONTINUE
      DO 130 I=0,INOR
      AGAR(I)=DBLE(DDAGAR(I))/DBLE(DDBGAR(INOR))
      BGAR(I)=DBLE(DDBGAR(I))/DBLE(DDBGAR(INOR))
      IF(AGAR(I).LE.0.0D0) GO TO 210
      IF(BGAR(I).LE.0.0D0) GO TO 210
  130 CONTINUE
      SGN=1.0D0
      CGAR(0)=AGAR(0)
      DO 135 I=2,INOR+1
      SGN=-SGN
      CGAR(I-1)=SGN*(BGAR(I-2)+AGAR(I-1))
  135 CONTINUE
      CGAR(INOR+1)=-SGN
      CALL ALROOT(CGAR,INOR+1,SIGX0,LFAIL)
      IF(LFAIL) GO TO 210
*
* NEWTON IMPROVEMENT OF THE ROOTS.
      DO 138 I=1,INOR+1
      CCC=0.0D0
      XCC=1.0D0
      DO 136 J=0,INOR+1
      CCC=CCC+CGAR(J)*XCC
      XCC=XCC*SIGX0(I)
  136 CONTINUE
      DCC=0.0D0
      XCC=1.0D0
      DO 137 J=1,INOR+1
      DCC=DCC+CGAR(J)*XCC*REAL(J)
      XCC=XCC*SIGX0(I)
  137 CONTINUE
      SIGX0(I)=SIGX0(I)-CCC/DCC
  138 CONTINUE
*
      IF(LREAL) THEN
         DO 140 I=1,INOR+1
         IF(DBLE(SIGX0(I)).LT.1.0E-10) GO TO 210
         IF(DIMAG(SIGX0(I)).NE.0.0) GO TO 210
  140    CONTINUE
      ENDIF
*
* COMPUTE THE WEIGHTS.
      DO 170 I=1,INOR+1
      CC=(1.0D0,0.0D0)
      DD=0.0D0
      DO 150 JNOR=0,INOR
      DD=DD+BGAR(JNOR)*CC
      CC=-CC*SIGX0(I)
  150 CONTINUE
      DO 160 J=1,INOR+1
      IF(J.NE.I) DD=DD/(SIGX0(J)-SIGX0(I))
  160 CONTINUE
      WEIGH(I)=DD
  170 CONTINUE
*
* TEST THE ACCURACY OF THE PADE APPROXIMATION.
      PREC1=0.0
      DO 190 I=1,N
      CC=0.0D0
      DD=0.0D0
      DO 180 JNOR=1,INOR+1
      CC=CC+WEIGH(JNOR)/(SIGX0(JNOR)+X(I))
      DD=DD+WEIGH(JNOR)*SIGX0(JNOR)/(SIGX0(JNOR)+X(I))
  180 CONTINUE
      YAPPR=DBLE(DD/CC)
      PREC1=MAX(PREC1,ABS(REAL(YAPPR)/Y(I)-1.0))
  190 CONTINUE
*
      IF(PREC1.LT.0.95*PREC) THEN
         NOR=INOR
         PREC=PREC1
         DO 200 I=0,NOR
         A(I)=AGAR(I)
         B(I)=BGAR(I)
  200    CONTINUE
      ENDIF
  210 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(V,AA)
      RETURN
      END
