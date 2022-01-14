*DECK LIBBAS
      SUBROUTINE LIBBAS(NISO,AT,AKT,AMT,T,IX,V,DV,NDTE,P,XS,X,E)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Scattering kernel based on the free gas model of Brown and St. John.
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
* NISO    number of terms in the model:
*         NISO=1 : pure free gas model;
*         NISO=2 : Brown and St. John model.
* AT      potential microscopic cross section.
* AKT     exponential constant in the model. Equal to zero for the
*         pure free gas model.
* AMT     isotope mass divided by neutron mass.
* T       absolute temperature divided by 293.6K.
* IX      number of thermal groups.
* V       neutron velocities.
* DV      used to transform velocity to energy.
* NDTE    first dimension of matrix P.
*
*Parameters: output
* P       scattering kernel. The first index is for secondary neutrons.
* XS      scattering microscopic cross section.
*
*Parameters: scratch
* X       temporary storage.
* E       temporary storage.
*
*Reference:
* H. C. Honeck, 'The distribution of thermal neutrons in space and
* energy in reactor lattices. Part 1: theory', Nucl. Sci. Eng., 8,
* 193 (1960).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER  NISO,IX,NDTE
      REAL     AT(NISO),AKT(NISO),AMT(NISO),T,V(IX),DV(IX),P(NDTE,IX),
     1         XS(IX),X(IX),E(IX)
*----
* LOCAL VARIABLES
*----
      INTEGER  I,N,J
      REAL     AM,BETS,BET,TAUS,TAU,ALPHA,THETA,ZETA,OMEG,CONST,
     1         CONST1,TEM,WA,WB,EA,EB,AFERF
*
      DO 110 I=1,IX
      XS(I)=0.0
      DO 100 J=1,IX
      P(I,J)=0.0
  100 CONTINUE
  110 CONTINUE
      DO 500 N=1,NISO
      IF(AT(N).EQ.0.0) GO TO 500
      AM=AMT(N)
      BETS=AM/T
      BET=SQRT(BETS)
      DO 120 I=1,IX
      X(I)=BET*V(I)
      E(I)=X(I)*X(I)
  120 CONTINUE
      TAUS=BETS/(BETS+AKT(N))
      TAU=SQRT(TAUS)
      ALPHA=AKT(N)*TAUS/BETS
      THETA=(AM+1.0)/(2.0*AM*TAU)
      ZETA=TAU-THETA
      OMEG=TAUS*(BETS+(AM+1.0)*AKT(N))/BETS
      CONST=(AT(N)*TAU*TAUS*(AM+1.0)*(AM+1.0))/(4.0*AM*OMEG)
      DO 400 I=1,IX
      DO 300 J=I,IX
      CONST1=CONST*X(I)/X(J)
      WA=ALPHA*E(J)
      TEM=0.0
      IF(WA.GE.50.0) GO TO 250
      EA=AFERF((THETA*X(I))+(ZETA*X(J)))+AFERF((THETA*X(I))-(ZETA*X(J)))
      TEM=CONST1*EA*EXP(-WA)
  250 WB=(OMEG*E(I)-E(J))/AM
      IF(WB.GE.50.0) GO TO 260
      EB=AFERF((THETA*X(J))-(ZETA*X(I)))-AFERF((ZETA*X(I))+(THETA*X(J)))
      TEM=TEM+CONST1*EB*EXP(-WB)
  260 IF(TEM.LE.1.E-15) GO TO 350
      P(I,J)=TEM
  300 CONTINUE
  350 EA=TAU*X(I)
      WA=ALPHA*E(I)
      IF(WA.LT.50.0) GO TO 352
      WA=0.0
      GOTO 353
  352 WA=EXP(-WA)
  353 EB=WA*AFERF(EA)*(EA+(0.5/EA))
      IF(E(I).LT.50.0) GO TO 355
      WB=0.0
      GOTO 356
  355 WB=EXP(-E(I))
  356 XS(I)=XS(I)+(AT(N)*TAUS*TAUS/BET)*(EB+0.5641896*WB)
  400 CONTINUE
  500 CONTINUE
      DO 610 I=1,IX
      E(I)=0.0
      WA=(V(I)*V(I))/T
      IF(WA.GE.50.0) GO TO 610
      E(I)=V(I)*V(I)*EXP(-WA)
  610 CONTINUE
      DO 630 I=1,IX
      DO 620 J=I,IX
      IF(E(I).LE.1.E-20) GO TO 620
      P(J,I)=P(I,J)*E(J)/E(I)
  620 CONTINUE
  630 CONTINUE
      DO 650 J=1,IX
      TEM=0.0
      DO 640 I=1,IX
      TEM=TEM+P(I,J)*DV(I)
  640 CONTINUE
      P(J,J)=((XS(J)-TEM)/DV(J))+P(J,J)
  650 CONTINUE
      DO 690 I=1,IX
      XS(I)=XS(I)/V(I)
      DO 680 J=1,IX
      P(I,J)=P(I,J)*DV(J)/V(I)
  680 CONTINUE
  690 CONTINUE
      RETURN
      END
