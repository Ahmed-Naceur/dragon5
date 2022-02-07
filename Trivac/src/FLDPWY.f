*DECK FLDPWY
      SUBROUTINE FLDPWY(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,
     > DIFF,F2Y,F3W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* compute the Piolat contribution to the current-current tranverse
* couplings for the Thomas-Raviart-Schneider method.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* LL4W    number of currents in direction W.
* LL4X    number of currents in direction X.
* LL4Y    number of currents in direction Y.
* NBLOS   number of lozenges in one ADI direction.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* CTRAN   tranverse coupling Piolat unit matrix.
* IPERT   mixture permutation index.
* KN      ADI permutation indices for the volumes and currents.
* DIFF    inverse diffusion coefficients.
* F2W     right-hand-side vector in direction W.
* F2X     right-hand-side vector in direction X.
* F2Y     right-hand-side vector in direction Y.
*
*Parameters: output
* F3W     result of matrix multiplication in direction W.
* F3X     result of matrix multiplication in direction X.
* F3Y     result of matrix multiplication in direction Y.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LL4W,LL4X,LL4Y,NBLOS,IELEM,IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL DIFF(NBLOS),F2Y(LL4Y),F3W(LL4W)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*
      NELEM=(IELEM+1)*IELEM
      NELEH=NELEM*IELEM
      NUM=0
      DO 30 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 30
      NUM=NUM+1
      ITRS=KN(NUM,3)
      DINV=DIFF(KEL)
      DO 25 I1=0,IELEM-1
      DO 20 I0=1,NELEM
      I=I1*NELEM+I0
      KNW1=KN(ITRS,3+I)
      IF(KNW1.EQ.0) GO TO 20
      INW1=ABS(KNW1)
      DO 10 J0=1,NELEM
      J=I1*NELEM+J0
      KNY2=KN(NUM,3+5*NELEH+J)
      IF(KNY2.EQ.0) GO TO 10
      INY2=ABS(KNY2)-LL4W-LL4X
      SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
      F3W(INW1)=F3W(INW1)-SG*DINV*REAL(CTRAN(I0,J0))*F2Y(INY2)
   10 CONTINUE
   20 CONTINUE
   25 CONTINUE
   30 CONTINUE
      RETURN
      END
*
      SUBROUTINE FLDPWX(LL4W,LL4X,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     > F2X,F3W)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LL4W,LL4X,NBLOS,IELEM,IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL DIFF(NBLOS),F2X(LL4X),F3W(LL4W)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*
      NELEM=(IELEM+1)*IELEM
      NELEH=NELEM*IELEM
      NUM=0
      DO 60 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 60
      NUM=NUM+1
      DINV=DIFF(KEL)
      DO 55 I1=0,IELEM-1
      DO 50 I0=1,NELEM
      I=I1*NELEM+I0
      KNX1=KN(NUM,3+2*NELEH+I)
      IF(KNX1.EQ.0) GO TO 50
      INX1=ABS(KNX1)-LL4W
      DO 40 J0=1,NELEM
      J=I1*NELEM+J0
      KNW2=KN(NUM,3+NELEH+J)
      IF(KNW2.EQ.0) GO TO 40
      INW2=ABS(KNW2)
      SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
      F3W(INW2)=F3W(INW2)-SG*DINV*REAL(CTRAN(I0,J0))*F2X(INX1)
   40 CONTINUE
   50 CONTINUE
   55 CONTINUE
   60 CONTINUE
      RETURN
      END
*
      SUBROUTINE FLDPXW(LL4W,LL4X,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     > F2W,F3X)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LL4W,LL4X,NBLOS,IELEM,IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL DIFF(NBLOS),F2W(LL4W),F3X(LL4X)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*
      NELEM=(IELEM+1)*IELEM
      NELEH=NELEM*IELEM
      NUM=0
      DO 90 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 90
      NUM=NUM+1
      DINV=DIFF(KEL)
      DO 85 I1=0,IELEM-1
      DO 80 I0=1,NELEM
      I=I1*NELEM+I0
      KNX1=KN(NUM,3+2*NELEH+I)
      IF(KNX1.EQ.0) GO TO 80
      INX1=ABS(KNX1)-LL4W
      DO 70 J0=1,NELEM
      J=I1*NELEM+J0
      KNW2=KN(NUM,3+NELEH+J)
      IF(KNW2.EQ.0) GO TO 70
      INW2=ABS(KNW2)
      SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
      F3X(INX1)=F3X(INX1)-SG*DINV*REAL(CTRAN(I0,J0))*F2W(INW2)
   70 CONTINUE
   80 CONTINUE
   85 CONTINUE
   90 CONTINUE
      RETURN
      END
*
      SUBROUTINE FLDPXY(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     > F2Y,F3X)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LL4W,LL4X,LL4Y,NBLOS,IELEM,IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL DIFF(NBLOS),F2Y(LL4Y),F3X(LL4X)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*
      NELEM=(IELEM+1)*IELEM
      NELEH=NELEM*IELEM
      NUM=0
      DO 120 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 120
      NUM=NUM+1
      DINV=DIFF(KEL)
      DO 115 I1=0,IELEM-1
      DO 110 I0=1,NELEM
      I=I1*NELEM+I0
      KNY1=KN(NUM,3+4*NELEH+I)
      IF(KNY1.EQ.0) GO TO 110
      INY1=ABS(KNY1)-LL4W-LL4X
      DO 100 J0=1,NELEM
      J=I1*NELEM+J0
      KNX2=KN(NUM,3+3*NELEH+J)
      IF(KNX2.EQ.0) GO TO 100
      INX2=ABS(KNX2)-LL4W
      SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
      F3X(INX2)=F3X(INX2)-SG*DINV*REAL(CTRAN(I0,J0))*F2Y(INY1)
  100 CONTINUE
  110 CONTINUE
  115 CONTINUE
  120 CONTINUE
      RETURN
      END
*
      SUBROUTINE FLDPYX(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     > F2X,F3Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LL4W,LL4X,LL4Y,NBLOS,IELEM,IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL DIFF(NBLOS),F2X(LL4X),F3Y(LL4Y)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*
      NELEM=(IELEM+1)*IELEM
      NELEH=NELEM*IELEM
      NUM=0
      DO 150 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 150
      NUM=NUM+1
      DINV=DIFF(KEL)
      DO 145 I1=0,IELEM-1
      DO 140 I0=1,NELEM
      I=I1*NELEM+I0
      KNY1=KN(NUM,3+4*NELEH+I)
      IF(KNY1.EQ.0) GO TO 140
      INY1=ABS(KNY1)-LL4W-LL4X
      DO 130 J0=1,NELEM
      J=I1*NELEM+J0
      KNX2=KN(NUM,3+3*NELEH+J)
      IF(KNX2.EQ.0) GO TO 130
      INX2=ABS(KNX2)-LL4W
      SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
      F3Y(INY1)=F3Y(INY1)-SG*DINV*REAL(CTRAN(I0,J0))*F2X(INX2)
  130 CONTINUE
  140 CONTINUE
  145 CONTINUE
  150 CONTINUE
      RETURN
      END
*
      SUBROUTINE FLDPYW(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     > F2W,F3Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LL4W,LL4X,LL4Y,NBLOS,IELEM,IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL DIFF(NBLOS),F2W(LL4W),F3Y(LL4Y)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*
      NELEM=(IELEM+1)*IELEM
      NELEH=NELEM*IELEM
      NUM=0
      DO 180 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 180
      NUM=NUM+1
      ITRS=KN(NUM,3)
      DINV=DIFF(KEL)
      DO 175 I1=0,IELEM-1
      DO 170 I0=1,NELEM
      I=I1*NELEM+I0
      KNW1=KN(ITRS,3+I)
      IF(KNW1.EQ.0) GO TO 170
      INW1=ABS(KNW1)
      DO 160 J0=1,NELEM
      J=I1*NELEM+J0
      KNY2=KN(NUM,3+5*NELEH+J)
      IF(KNY2.EQ.0) GO TO 160
      INY2=ABS(KNY2)-LL4W-LL4X
      SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
      F3Y(INY2)=F3Y(INY2)-SG*DINV*REAL(CTRAN(I0,J0))*F2W(INW1)
  160 CONTINUE
  170 CONTINUE
  175 CONTINUE
  180 CONTINUE
      RETURN
      END
