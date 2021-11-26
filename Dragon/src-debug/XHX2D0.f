*DECK XHX2D0
      SUBROUTINE XHX2D0 (NGPT,ZGAUS,WGAUS,COTE,SIGT,TRONC,PII,PIS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute DP0 collision, leakage and transmission probabilities for
* hexagonal 2D geometries.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): M. Ouisloumen
*
*Parameters: input
* NGPT    number of Gauss integration points.
* COTE    length of one of sides of the hexagon.
* SIGT    total cross section.
* TRONC   voided block cutoff criterion.
* ZGAUS   Gauss-Legendre integration points.
* WGAUS   Gauss-Legendre integration weights.
*
*Parameters: output
* PII     volume to volume reduced probability.
* PIS     leakage probability (PIS(i) volume to side i).
* PSS     transmission probability (PSS(i,j) side i to side j).
*
*Comments:
* Faces identification for hexagon
*                   side 4
*                  xxxxxxxx
*                 x        x
*       side 5   x          x side 3
*               x            x
*              x              x
*               x            x
*       side 6   x          x side 2
*                 x        x
*                  xxxxxxxx
*                   side 1
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER    NGPT
      REAL       ZGAUS(NGPT),WGAUS(NGPT),COTE,SIGT,TRONC,PII,PIS(6),
     +           PSS(6,6)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600,MKI5=600)
      PARAMETER (PI=3.141592653589793,SQRT3=1.732050807568877,
     +           SQRT2=1.414213562373095,ALOG2=.693147180559945,
     +           ALOG3=1.0986122886681097,ALOGX=.7676517525907618)
*
      REAL       TAU(3),FKI3(3),FKI4(3),FKI5(3),FKI6(3)
      INTEGER    IROT(6,6)
      DOUBLE PRECISION P(3),PIS10
*----
*  ASSUME THAT BICKLEY KI TABLES HAVE THE SAME TABULATION POINTS AND
*  THE SAME TRUNCATION LIMIT
*----
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5),PAS5,XLIM5,L5
*
      SAVE IROT
      DATA IROT /
     +          0, 1, 2, 3, 2, 1,
     +          1, 0, 1, 2, 3, 2,
     +          2, 1, 0, 1, 2, 3,
     +          3, 2, 1, 0, 1, 2,
     +          2, 3, 2, 1, 0, 1,
     +          1, 2, 3, 2, 1, 0/
*
      FUNC3(X,K)=BI3(K)+X*(BI31(K)+X*BI32(K))
      FUNC4(X,K)=BI4(K)+X*(BI41(K)+X*BI42(K))
      FUNC5(X,K)=BI5(K)+X*(BI51(K)+X*BI52(K))
*----
*  INITIALIZATION OF COLLISION PROBABILITIES
*----
      P(3)=0.
      P(2)=0.
      P(1)=0.
*----
*  COMPUTE CORDE =4*V*SIGMA/S=SQRT(3)*COTE*SIGMA  (AVERAGE CORDE)
*----
      CORDE=SQRT3*COTE*SIGT
      IF (CORDE.LE.TRONC) GO TO 300
*----
*  CONSIDER EXPLICIT INTEGRATION OF F PSS
*----
      PI12=PI/12.
*
*  CONSIDER TWO CASES  1) IF ZGAUS(I)<0 => 1/COSFI>1/SINA>1/COSA
*                      2) IF ZGAUS(I)>0 => 1/COSFI>1/COSA>1/SINA
*
      NGPT2=IFIX(FLOAT(NGPT)/2.)
      DO 50 I=1,NGPT
      FI=PI12*(1.+ZGAUS(I))
      COSFI=COS(FI)
      SINFI=SIN(FI)
      COSA=SQRT3*COSFI-SINFI
      AUX=SQRT3*SINFI
      SINA=AUX+COSFI
      SINB=COSFI-AUX
*----
*  OPTICAL LENGHTS
*----
      TAU(1)=CORDE/COSFI
      IAUX1=2
      IAUX2=3
      IF(I.GT.NGPT2) THEN
        IAUX1=3
        IAUX2=2
      ENDIF
      TAU(IAUX1)=CORDE/SINA
      TAU(IAUX2)=CORDE/COSA
*
      LB=4
      IF(TAU(1).LT.XLIM3) THEN
        LB=1
      ELSEIF(TAU(2).LT.XLIM3) THEN
        LB=2
      ELSEIF(TAU(3).LT.XLIM3) THEN
        LB=3
      ENDIF
*
      DO 10 J=LB,3
      K1=NINT(PAS3*TAU(J))
      FKI3(J)=FUNC3(TAU(J),K1)
      FKI4(J)=FUNC4(TAU(J),K1)
      FKI5(J)=FUNC5(TAU(J),K1)
      FKI6(J)=.8*FKI4(J)+.2*TAU(J)*(FKI3(J)-FKI5(J))
   10 CONTINUE
*
      WEIGHT=0.0
      GO TO (20,25,30,50),LB
*----
*  PSS CALCULATION
*----
   20 WEIGHT=WGAUS(I)
      P(3)=P(3)+WEIGHT*SINB*FKI3(1)
      P(2)=P(2)+WEIGHT*COSFI*SINA*(FKI4(IAUX1)-FKI4(1))
      GO TO 30
*
   25 WEIGHT=WGAUS(I)
      P(2)=P(2)+WEIGHT*COSFI*SINA*FKI4(IAUX1)
   30 P(1)=P(1)+WEIGHT*SINFI*COSA*(BI4(0)-FKI4(IAUX2))
   50 CONTINUE
*----
*  NORMALIZATION
*----
      X1=1./(3.*CORDE)
      P(1)=X1*P(1)
      P(2)=X1*P(2)
      P(3)=P(3)/3.
      PIS10=(1.-2.*(P(1)+P(2))-P(3))/(6.*CORDE)
      PII=(1.-6.*REAL(PIS10))/SIGT
*
      GO TO 350
*----
*  USE SERIES EXPANSION FOR CORDE->0: TAYLOR SERIES OF KI FUNCTIONS
*----
  300 TAU0=CORDE*.5
      TAU02=TAU0*TAU0
      AUX=SQRT3/PI
      AUX0=1.-.5*SQRT3
      AUX1=AUX*ALOG3-.33333333333333333
      AUX2=2./SQRT3-(2.+ALOGX)/3.
      P(1)=AUX0-.5*(TAU0*AUX1-TAU02*AUX2)
      AUX1=AUX*(2.5*ALOG3-4.*ALOG2)-.5
      AUX2=5.*SQRT3-9.-2.*ALOG3+.5*ALOGX
      P(2)=SQRT3-1.5+TAU0*AUX1-TAU02*AUX2/3.
      AUX0=2.*AUX0
      AUX1=AUX*(.5*ALOG3-ALOG2)
      P(3)=AUX0-8.*(1./6.+AUX1)*TAU0-4.*TAU02*(AUX0-.5*ALOG3)
      PII=COTE*SQRT3*(4.*SQRT3-8.-2.*ALOGX+10.*ALOG3)/12.
      PIS10=(1.-SIGT*PII)/6.
*
  350 CONTINUE
*----
*  TRANSMISSION MATRIX
*----
      DO  59 I=1,6
      DO  58 J=1,6
      IB=IROT(J,I)
      IF(IB.GT.0) THEN
         PSS(I,J)=REAL(P(IB))
      ELSE
         PSS(I,J)=0.
      ENDIF
   58 CONTINUE
   59 CONTINUE
*----
*  LEAKAGE PRABABILITIES MATRIX
*----
      DO 56 I=1,6
      PIS(I)=REAL(PIS10)
   56 CONTINUE
*
      RETURN
      END
