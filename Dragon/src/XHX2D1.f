*DECK XHX2D1
      SUBROUTINE XHX2D1 (NGPT,ZGAUS,WGAUS,COTE,SIGT,TRONC,PII,PIS,PSS,
     + P)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute DP1 collision, leakage and transmission probabilities for
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
*Parameters: scratch
* P       undefined.
*
*Comments:
*  Faces identification for hexagon
*                                         side a,b,c
*                  side 4,5,6             dir a -> isotropic
*                   xxxxxxxx              dir c -> tangent to surface
*                  x        x             dir b -> normal to surface
*     side 7,8,9  x          x side 1,2,3
*                x            x
*               x              x
*                x            x
*  side 10,11,12  x          x side 16,17,18
*                  x        x
*                   xxxxxxxx
*                side 13,14,15
*
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL FUNC3,FUNC4,FUNC5,X
      INTEGER    NGPT
      REAL       ZGAUS(NGPT),WGAUS(NGPT),COTE,SIGT,TRONC,PII,PIS(18),
     +           PSS(18,18)
      DOUBLE PRECISION P(16)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600,MKI5=600)
      PARAMETER (PI=3.141592653589793,SQRT3=1.732050807568877,
     +           SQRT2=1.414213562373095,ALOG2=.693147180559945,
     +           ALOG3=1.0986122886681097,ALOGX=.7676517525907618)
*
      REAL       TAU(3),FKI3(3),FKI4(3),FKI5(3),FKI6(3)
      INTEGER    IROT(18,18)
      DOUBLE PRECISION PIS10,PIS11
*----
*  ASSUME THAT BICKLEY KI TABLES HAVE THE SAME TABULATION POINTS AND
*  THE SAME TRUNCATION LIMIT.
*----
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5),PAS5,XLIM5,L5
*
      SAVE IROT
      DATA IROT /
     +          0, 0, 0, 1, 2,-3, 7, 8,-9,13,14, 0, 7, 8, 9, 1, 2, 3,
     +          0, 0, 0, 2, 4, 5, 8,10,11,14,15, 0, 8,10,-11, 2, 4,-5,
     +          0, 0, 0, 3,-5, 6, 9,-11,12,0, 0,16,-9,11,12,-3, 5,6,
     +          1, 2, 3, 0, 0, 0, 1, 2,-3, 7, 8,-9,13,14, 0, 7, 8, 9,
     +          2, 4,-5, 0, 0, 0, 2, 4, 5, 8,10,11,14,15, 0, 8,10,-11,
     +         -3, 5, 6, 0, 0, 0, 3,-5, 6, 9,-11,12,0, 0,16,-9,11,12,
     +          7, 8, 9, 1, 2, 3, 0, 0, 0, 1, 2,-3, 7, 8,-9,13,14, 0,
     +          8,10,-11, 2, 4,-5, 0, 0, 0, 2, 4, 5, 8,10,11,14,15, 0,
     +         -9,11,12,-3, 5, 6, 0, 0, 0, 3,-5, 6, 9,-11,12,0, 0,16,
     +          13,14, 0, 7, 8, 9, 1, 2, 3, 0, 0, 0, 1, 2,-3, 7, 8,-9,
     +          14,15, 0, 8,10,-11, 2, 4,-5, 0, 0, 0, 2, 4, 5, 8,10,11,
     +          0, 0,16,-9,11,12,-3, 5, 6, 0, 0, 0, 3,-5, 6, 9,-11,12,
     +          7, 8,-9, 13,14, 0, 7, 8, 9, 1, 2, 3, 0, 0, 0, 1, 2,-3,
     +          8,10,11,14,15, 0, 8,10,-11, 2, 4,-5, 0, 0, 0, 2, 4, 5,
     +          9,-11,12,  0, 0,16,-9,11,12,-3, 5, 6, 0, 0, 0, 3,-5, 6,
     +          1, 2,-3, 7, 8,-9, 13,14, 0, 7, 8, 9, 1, 2, 3, 0, 0, 0,
     +          2, 4, 5, 8,10,11,14,15, 0, 8,10,-11, 2, 4,-5, 0, 0, 0,
     +          3,-5, 6, 9,-11,12,  0, 0,16,-9,11,12,-3, 5, 6, 0, 0, 0/
*
      FUNC3(X,K)=BI3(K)+X*(BI31(K)+X*BI32(K))
      FUNC4(X,K)=BI4(K)+X*(BI41(K)+X*BI42(K))
      FUNC5(X,K)=BI5(K)+X*(BI51(K)+X*BI52(K))
*----
*  INITIALIZATION OF COLLISION PROBABILITIES
*----
      CALL XDDSET(P,16,0.0D0)
*----
*  COMPUTE CORDE =4*V*SIGMA/S=SQRT(3)*COTE*SIGMA  (AVERAGE CORDE)
*----
      S2S3=SQRT2*SQRT3
      S3DS2=SQRT3/SQRT2
      CAUX=SQRT3*COTE
      CORDE=CAUX*SIGT
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
      AUX=SQRT3*COSFI
      COSA=AUX-SINFI
      COSB=AUX+SINFI
      AUX=SQRT3*SINFI
      SINA=AUX+COSFI
      SINB=COSFI-AUX
*----
*  WEIGHTS TIMES BICKLEY NAYLOR FUNCTIONS
*----
      W006=SINFI*COSA
      W106=W006*COSB
      W005=COSFI*SINA
      W105=W005*COSB
      W204=COSFI*SINB
      W114=SINFI*SINFI*SINB
      W224=COSFI*W204
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
      WEIGHT=0.0
      GO TO (20,25,30,50),LB
*----
*  PSS CALCULATION
*----
   20 WEIGHT=WGAUS(I)
      P(13)=P(13)+WEIGHT*SINB*FKI3(1)
      P(14)=P(14)+WEIGHT*W204*FKI4(1)
      P(16)=P(16)+WEIGHT*W114*FKI5(1)
      P(15)=P(15)+WEIGHT*W224*FKI5(1)
      AUX5=WEIGHT*W005
      P(7)=P(7)+AUX5*(FKI4(IAUX1)-FKI4(1))
      P(9)=P(9)+WEIGHT*W105*(FKI5(IAUX1)-FKI5(1))
      AUX=AUX5*(FKI6(IAUX1)-FKI6(1))
      P(11)=P(11)+AUX
      P(10)=P(10)+W005*AUX
      GO TO 30
*
   25 WEIGHT=WGAUS(I)
      AUX5=WEIGHT*W005
      P(7)=P(7)+AUX5*FKI4(IAUX1)
      P(9)=P(9)+WEIGHT*W105*FKI5(IAUX1)
      AUX=AUX5*FKI6(IAUX1)
      P(11)=P(11)+AUX
      P(10)=P(10)+W005*AUX
   30 P(1)=P(1)+WEIGHT*W006*(BI4(0)-FKI4(IAUX2))
      P(3)=P(3)+WEIGHT*W106*(BI5(0)-FKI5(IAUX2))
      AUX=W006*WEIGHT*(.533333333333333-FKI6(IAUX2))
      P(5)=P(5)+AUX
      P(4)=P(4)+AUX*W006
   50 CONTINUE
*----
*  NORMALIZATION
*----
      CORDI=1./CORDE
      X1=CORDI/3.
      X2=CORDI*SQRT3
      X3=3.*CORDI
      X4=X2/SQRT2
      P(1)=X1*P(1)
      P(3)=X2*P(3)/6.
      P(5)=-X4*P(5)
      P(4)=X3*P(4)
      P(7)=X1*P(7)
      P(9)=X1*P(9)*.5
      P(11)=-X4*P(11)
      P(10)=X3*P(10)
      P(13)=P(13)/3.
      P(14)=SQRT2*P(14)
      P(16)=4.*P(16)/3.
      P(15)=6.*P(15)
      AUX=2.*SQRT2
      P(2)=S3DS2*P(3)-AUX*P(1)
      P(8)=3.*S3DS2*P(9)-AUX*P(7)
      P(14)=P(14)-AUX*P(13)
      COEF=1./(6.*CORDE)
      PIS10=(1.-2.*(P(1)+P(7))-P(13))*COEF
      PIS(1)=REAL(PIS10)
      PIS11=-(2.*(P(2)+P(8))+P(14))*COEF
      PIS(2)=REAL(PIS11)
      PII=(1.-6.*PIS(1))/SIGT
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
      P(5)=-S2S3*(1.125*AUX0-.5*TAU0*AUX1+.375*TAU02*AUX2)
      P(3)=(-(1.5*ALOGX+1.-SQRT3)*TAU0*.25+(1.-TAU02)/9+
     +             AUX*TAU02/3.)*SQRT3
      P(4)=2.25*(1.25*SQRT3-2.)-TAU0*(2.-3.*AUX)+TAU02*
     +           (2.25*ALOGX-1.5)
      XAUX=5.*SQRT3-9.
      XAUX0=SQRT3-1.5
      AUX1=AUX*(2.5*ALOG3-4.*ALOG2)-.5
      AUX2=XAUX-2.*ALOG3+.5*ALOGX
      P(7)=XAUX0+TAU0*AUX1-TAU02*AUX2/3.
      P(11)=-S2S3*(1.125*XAUX0+TAU0*AUX1-.25*AUX2*TAU02)
      P(9)=SQRT3/9.+.25*(XAUX-.5*SQRT3*(2.*ALOG3-ALOGX))*TAU0-
     +           (9.*ALOG3-16.*ALOG2-3.)*TAU02/(3.*PI)
      P(10)=2.25*(SQRT3-.75)+TAU0*(3.*AUX-6.)-1.5*TAU02*
     +           (9.*(2.-SQRT3)+1.5*ALOGX-6.*ALOG3)
      AUX0=2.*AUX0
      AUX1=AUX*(.5*ALOG3-ALOG2)
      P(13)=AUX0-8.*(1./6.+AUX1)*TAU0-4.*TAU02*(AUX0-.5*ALOG3)
      P(15)=4.5*(2.5-SQRT3)-8.*TAU0+36.*TAU02*AUX0
      AUX2=XAUX-2.*ALOG3+.5*ALOGX
      P(16)=3.5-2.*SQRT3-64.*TAU0*(AUX1+1./12.)/3.+8.*TAU02*
     +           (.5*ALOG3-2.*AUX0)
      PSQ3=SQRT3/PI
      P(14)=SQRT2*(2./3.-6.*TAU0*AUX0+TAU02*(4.+24.*PSQ3*
     + (.5*ALOG3-ALOG2)))
      PII=COTE*SQRT3*(4.*SQRT3-8.-2.*ALOGX+10.*ALOG3)/12.
      PIS(1)=(1.-SIGT*PII)/6.
      PIS(2)=-SQRT2*(2.5-2.25*ALOG3+TAU0*(9.+PSQ3*(4.-8.*ALOG2+
     +  3.*ALOG3)-8./SQRT3-(20.*ALOG3-4.*ALOGX)/3.))/12.
      AUX=2.*SQRT2
      P(2)=S3DS2*P(3)-AUX*P(1)
      P(8)=3.*S3DS2*P(9)-AUX*P(7)
      P(14)=P(14)-AUX*P(13)
*
  350 CONTINUE
*----
*  COMPUTE REMAINING PROBABILITIES
*----
      P(4)=P(4)-4.*SQRT3*P(3)+8.*P(1)
      P(5)=P(5)+AUX*P(3)
      P(6)=(S2S3*P(5)+8.*(-SQRT3*P(3)+P(1))-P(4))*2./9.
      P(10)=P(10)-4.*SQRT2*P(8)-8.*P(7)
      P(11)=P(11)+AUX*P(9)
      P(12)=(8.*(P(7)-SQRT3*P(9))-S2S3*P(11)-P(10))*2./9.
      P(15)=P(15)-4.*SQRT2*P(14)-8.*P(13)
      P(3)=-P(3)
      P(5)=-P(5)
      P(9)=-P(9)
      P(11)=-P(11)
*----
*  TRANSMISSION MATRIX
*----
      DO 59 I=1,18
      DO 58 J=1,18
      IB=IROT(J,I)
      IF(IB.LT.0) THEN
         PSS(I,J)=-REAL(P(-IB))
      ELSEIF(IB.GT.0) THEN
         PSS(I,J)=REAL(P(IB))
      ELSE
         PSS(I,J)=0.
      ENDIF
   58 CONTINUE
   59 CONTINUE
*----
*  LEAKAGE PRABABILITIES MATRIX
*----
      PIS(3)=0.
      K=3
      DO 56 I=1,5
      K=K+3
      PIS(K-2)=PIS(1)
      PIS(K-1)=PIS(2)
      PIS(K)=0.
   56 CONTINUE
*
      RETURN
      END
