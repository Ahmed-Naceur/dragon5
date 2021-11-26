*DECK TRIRWW
      SUBROUTINE TRIRWW (IR,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,MUW,
     1 A11W,ISPLH,R,Q,RH,QH,RT,QT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a mesh corner finite difference
* discretization in hexagonal geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud
*
*Parameters: input
* IR      first dimension of matrix SGD.
* NEL     total number of finite elements.
* ll4     order of system matrices.
* VOL     volume of each element.
* MAT     mixture index assigned to each element.
* XSGD    nuclear properties, derivatives or first variations of
*         nuclear properties per material mixture:
*         XSGD(L,1): W-, X-, and Y-oriented diffusion coefficients;
*         XSGD(L,3): Z-oriented diffusion coefficients;
*         XSGD(L,4): removal macroscopic cross section.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* KN      element-ordered unknown list (dimensionned to KN(ICOF*NEL)
*         where ICOF=12 or 14).
* QFR     element-ordered boundary conditions.
* MUW     W-oriented compressed storage mode indices.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IPX     X-oriented permutation matrices.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
* ISPLH   hexagonal mesh-splitting flag:
*         =1 for complete hexagons; >1 for triangular elements.
* R       unit matrix.
* Q       unit matrix.
* RH      unit matrix.
* QH      unit matrix.
* RT      unit matrix.
* QT      unit matrix.
*
*Parameters: output
* A11W    W-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
* A11X    X-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
* A11Y    Y-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
* A11Z    Z-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(*),MUW(LL4),ISPLH
      REAL VOL(NEL),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),A11W(*),R(2,2),
     > Q(2,2),RH(6,6),QH(6,6),RT(3,3),QT(3,3)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(8,25)
      REAL R2DP(4)
      DOUBLE PRECISION QTHP(14,14),QTHZ(14,14),RTHG(14,14),
     >           HW(14,14),HX(14,14),HY(14,14),HZ(14,14)
      DOUBLE PRECISION RR,QQP,QQZ,VOL0,VOL1,DZ,VAR1
      DATA R2DP / 4*0.25 /
*----
*  ASSEMBLY OF MATRIX A11W
*----
      CALL TRIRMA(ISPLH,R,Q,RH,QH,RT,QT,LL,LC,ISR,QTHP,QTHZ,RTHG,HW,HX,
     > HY,HZ)
      NUM1=0
      NUM2=0
      VOL1=SIDE*SIDE
      DO 160 K=1,NEL
         L=MAT(K)
         IF(L.EQ.0) GO TO 160
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 150
         DZ=ZZ(K)
         DO 110 I=1,LL
            INW1=KN(NUM1+I)
            IF(INW1.EQ.0) GO TO 110
            KEY0=MUW(INW1)-INW1
            DO 100 J=1,LL
               INW2=KN(NUM1+J)
               IF(INW2.EQ.0) GO TO 100
               IF(INW2.EQ.INW1) THEN
                  QQP=QTHP(I,J)*DZ
                  QQZ=QTHZ(I,J)*VOL1/DZ
                  KEY=KEY0+INW2
                  VAR1=QQP*XSGD(L,1)+QQZ*XSGD(L,3)
                  A11W(KEY)=A11W(KEY)+REAL(VAR1)
               ELSE IF((INW2.LT.INW1).AND.(HW(I,J).NE.0.0)) THEN
                  QQP=QTHP(I,J)*HW(I,J)*DZ
                  KEY=KEY0+INW2
                  A11W(KEY)=A11W(KEY)+REAL(QQP)*XSGD(L,1)
               ENDIF
  100       CONTINUE
            RR=RTHG(I,I)*VOL1*DZ
            KEY=KEY0+INW1
            A11W(KEY)=A11W(KEY)+REAL(RR)*XSGD(L,4)
  110    CONTINUE
         DO 140 IC=1,8
            QFR1=QFR(NUM2+IC)
            IF(QFR1.EQ.0.0) GO TO 140
            IF(IC.LT.7) THEN
               DO 120 I1=1,4
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 120
                  KEY=MUW(INW1)
                  RR=R2DP(I1)
                  A11W(KEY)=A11W(KEY)+REAL(RR)*QFR1
  120          CONTINUE
            ELSE
               DO 130 I1=1,LC
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 130
                  KEY=MUW(INW1)
                  RR=RTHG(I1,I1)
                  A11W(KEY)=A11W(KEY)+REAL(RR)*QFR1
  130          CONTINUE
            ENDIF
  140    CONTINUE
  150    NUM1=NUM1+LL
         NUM2=NUM2+8
  160 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIRWX (IR,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,MUX,IPX,
     > A11X,ISPLH,R,Q,RH,QH,RT,QT)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(*),MUX(LL4),IPX(LL4),ISPLH
      REAL VOL(NEL),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),A11X(*),R(2,2),
     > Q(2,2),RH(6,6),QH(6,6),RT(3,3),QT(3,3)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(8,25)
      REAL R2DP(4)
      DOUBLE PRECISION QTHP(14,14),QTHZ(14,14),RTHG(14,14),
     >           HW(14,14),HX(14,14),HY(14,14),HZ(14,14)
      DOUBLE PRECISION RR,QQP,QQZ,VOL0,VOL1,DZ,VAR1
      DATA R2DP / 4*0.25 /
*----
*  ASSEMBLY OF MATRIX A11X
*----
      CALL TRIRMA(ISPLH,R,Q,RH,QH,RT,QT,LL,LC,ISR,QTHP,QTHZ,RTHG,HW,HX,
     > HY,HZ)
      NUM1=0
      NUM2=0
      VOL1=SIDE*SIDE
      DO 230 K=1,NEL
         L=MAT(K)
         IF(L.EQ.0) GO TO 230
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 220
         DZ=ZZ(K)
         DO 180 I=1,LL
            INW1=KN(NUM1+I)
            IF(INW1.EQ.0) GO TO 180
            INX1=IPX(INW1)
            KEY0=MUX(INX1)-INX1
            DO 170 J=1,LL
               INW2=KN(NUM1+J)
               IF(INW2.EQ.0) GO TO 170
               INX2=IPX(INW2)
               IF(INX2.EQ.INX1) THEN
                  QQP=QTHP(I,J)*DZ
                  QQZ=QTHZ(I,J)*VOL1/DZ
                  KEY=KEY0+INX2
                  VAR1=QQP*XSGD(L,1)+QQZ*XSGD(L,3)
                  A11X(KEY)=A11X(KEY)+REAL(VAR1)
               ELSE IF((INX2.LT.INX1).AND.(HX(I,J).NE.0.0)) THEN
                  QQP=QTHP(I,J)*HX(I,J)*DZ
                  KEY=KEY0+INX2
                  A11X(KEY)=A11X(KEY)+REAL(QQP)*XSGD(L,1)
               ENDIF
  170       CONTINUE
            RR=RTHG(I,I)*VOL1*DZ
            KEY=KEY0+INX1
            A11X(KEY)=A11X(KEY)+REAL(RR)*XSGD(L,4)
  180    CONTINUE
         DO 210 IC=1,8
            QFR1=QFR(NUM2+IC)
            IF(QFR1.EQ.0.0) GO TO 210
            IF(IC.LT.7) THEN
               DO 190 I1=1,4
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 190
                  INX1=IPX(INW1)
                  KEY=MUX(INX1)
                  RR=R2DP(I1)
                  A11X(KEY)=A11X(KEY)+REAL(RR)*QFR1
  190          CONTINUE
            ELSE
               DO 200 I1=1,LC
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 200
                  INX1=IPX(INW1)
                  KEY=MUX(INX1)
                  RR=RTHG(I1,I1)
                  A11X(KEY)=A11X(KEY)+REAL(RR)*QFR1
  200          CONTINUE
            ENDIF
  210    CONTINUE
  220    NUM1=NUM1+LL
         NUM2=NUM2+8
  230 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIRWY (IR,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,MUY,IPY,
     > A11Y,ISPLH,R,Q,RH,QH,RT,QT)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(*),MUY(LL4),IPY(LL4),ISPLH
      REAL VOL(NEL),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),A11Y(*),R(2,2),
     > Q(2,2),RH(6,6),QH(6,6),RT(3,3),QT(3,3)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(8,25)
      REAL R2DP(4)
      DOUBLE PRECISION QTHP(14,14),QTHZ(14,14),RTHG(14,14),
     >           HW(14,14),HX(14,14),HY(14,14),HZ(14,14)
      DOUBLE PRECISION RR,QQP,QQZ,VOL0,VOL1,DZ,VAR1
      DATA R2DP / 4*0.25 /
*----
*  ASSEMBLY OF MATRIX A11Y
*----
      CALL TRIRMA(ISPLH,R,Q,RH,QH,RT,QT,LL,LC,ISR,QTHP,QTHZ,RTHG,HW,HX,
     > HY,HZ)
      NUM1=0
      NUM2=0
      VOL1=SIDE*SIDE
      DO 300 K=1,NEL
         L=MAT(K)
         IF(L.EQ.0) GO TO 300
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 290
         DZ=ZZ(K)
         DO 250 I=1,LL
            INW1=KN(NUM1+I)
            IF(INW1.EQ.0) GO TO 250
            INY1=IPY(INW1)
            KEY0=MUY(INY1)-INY1
            DO 240 J=1,LL
               INW2=KN(NUM1+J)
               IF(INW2.EQ.0) GO TO 240
               INY2=IPY(INW2)
               IF(INY2.EQ.INY1) THEN
                  QQP=QTHP(I,J)*DZ
                  QQZ=QTHZ(I,J)*VOL1/DZ
                  KEY=KEY0+INY2
                  VAR1=QQP*XSGD(L,1)+QQZ*XSGD(L,3)
                  A11Y(KEY)=A11Y(KEY)+REAL(VAR1)
               ELSE IF((INY2.LT.INY1).AND.(HY(I,J).NE.0.0)) THEN
                  QQP=QTHP(I,J)*HY(I,J)*DZ
                  KEY=KEY0+INY2
                  A11Y(KEY)=A11Y(KEY)+REAL(QQP)*XSGD(L,1)
               ENDIF
  240       CONTINUE
            RR=RTHG(I,I)*VOL1*DZ
            KEY=KEY0+INY1
            A11Y(KEY)=A11Y(KEY)+REAL(RR)*XSGD(L,4)
  250    CONTINUE
         DO 280 IC=1,8
            QFR1=QFR(NUM2+IC)
            IF(QFR1.EQ.0.0) GO TO 280
            IF(IC.LT.7) THEN
               DO 260 I1=1,4
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 260
                  INY1=IPY(INW1)
                  KEY=MUY(INY1)
                  RR=R2DP(I1)
                  A11Y(KEY)=A11Y(KEY)+REAL(RR)*QFR1
  260          CONTINUE
            ELSE
               DO 270 I1=1,LC
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 270
                  INY1=IPY(INW1)
                  KEY=MUY(INY1)
                  RR=RTHG(I1,I1)
                  A11Y(KEY)=A11Y(KEY)+REAL(RR)*QFR1
  270          CONTINUE
            ENDIF
  280    CONTINUE
  290    NUM1=NUM1+LL
         NUM2=NUM2+8
  300 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIRWZ (IR,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,MUZ,IPZ,
     > A11Z,ISPLH,R,Q,RH,QH,RT,QT)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(*),MUZ(LL4),IPZ(LL4),ISPLH
      REAL VOL(NEL),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),A11Z(*),R(2,2),
     > Q(2,2),RH(6,6),QH(6,6),RT(3,3),QT(3,3)
*----
*  LOCAL VARIABLES
*----
      INTEGER ISR(8,25)
      REAL R2DP(4)
      DOUBLE PRECISION QTHP(14,14),QTHZ(14,14),RTHG(14,14),
     >           HW(14,14),HX(14,14),HY(14,14),HZ(14,14)
      DOUBLE PRECISION RR,QQP,QQZ,VOL0,VOL1,DZ,VAR1
      DATA R2DP / 4*0.25 /
*----
*  ASSEMBLY OF MATRIX A11Z
*----
      CALL TRIRMA(ISPLH,R,Q,RH,QH,RT,QT,LL,LC,ISR,QTHP,QTHZ,RTHG,HW,HX,
     > HY,HZ)
      NUM1=0
      NUM2=0
      VOL1=SIDE*SIDE
      DO 360 K=1,NEL
         L=MAT(K)
         IF(L.EQ.0) GO TO 360
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 350
         DZ=ZZ(K)
         DO 320 I=1,LL
            INW1=KN(NUM1+I)
            IF(INW1.EQ.0) GO TO 320
            INZ1=IPZ(INW1)
            KEY0=MUZ(INZ1)-INZ1
            DO 310 J=1,LL
               INW2=KN(NUM1+J)
               IF(INW2.EQ.0) GO TO 310
               INZ2=IPZ(INW2)
               IF(INZ2.EQ.INZ1) THEN
                  QQP=QTHP(I,J)*DZ
                  QQZ=QTHZ(I,J)*VOL1/DZ
                  KEY=KEY0+INZ2
                  VAR1=QQP*XSGD(L,1)+QQZ*XSGD(L,3)
                  A11Z(KEY)=A11Z(KEY)+REAL(VAR1)
               ELSE IF((INZ2.LT.INZ1).AND.(HZ(I,J).NE.0.0)) THEN
                  QQZ=QTHZ(I,J)*VOL1/DZ
                  KEY=KEY0+INZ2
                  A11Z(KEY)=A11Z(KEY)+REAL(QQZ)*XSGD(L,1)
               ENDIF
  310       CONTINUE
            RR=RTHG(I,I)*VOL1*DZ
            KEY=KEY0+INZ1
            A11Z(KEY)=A11Z(KEY)+REAL(RR)*XSGD(L,4)
  320    CONTINUE
         DO 340 IC=1,8
            QFR1=QFR(NUM2+IC)
            IF(QFR1.EQ.0.0) GO TO 340
            IF(IC.LT.7) THEN
               DO 330 I1=1,4
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 330
                  INZ1=IPZ(INW1)
                  KEY=MUZ(INZ1)
                  RR=R2DP(I1)
                  A11Z(KEY)=A11Z(KEY)+REAL(RR)*QFR1
  330          CONTINUE
            ELSE
               DO 335 I1=1,LC
                  I=ISR(IC,I1)
                  INW1=KN(NUM1+I)
                  IF(INW1.EQ.0) GO TO 335
                  INZ1=IPZ(INW1)
                  KEY=MUZ(INZ1)
                  RR=RTHG(I1,I1)
                  A11Z(KEY)=A11Z(KEY)+REAL(RR)*QFR1
  335          CONTINUE
            ENDIF
  340    CONTINUE
  350    NUM1=NUM1+LL
         NUM2=NUM2+8
  360 CONTINUE
      RETURN
      END
