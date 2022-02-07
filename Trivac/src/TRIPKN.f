*DECK TRIPKN
      SUBROUTINE TRIPKN (IELEM,LX,LY,LZ,L4,CYLIND,XXX,YYY,ZZZ,XX,YY,ZZ,
     1 DD,KN,QFR,IQFR,VOL,MAT,NCODE,ICODE,ZCODE,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a primal formulation of the finite element
* discretization in a 3-D geometry.
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
* IMPX    print parameter.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
* CYLIND  cylindrical geometry flag (set with CYLIND=.true.).
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NCODE   type of boundary condition applied on each side:
*         I=1: X-; I=2: X+; I=3: Y-; I=4: Y+; I=5: Z-; I=6: Z+;
*         NCODE(I)=1: VOID;  NCODE(I)=2: REFL;  NCODE(I)=4: TRAN;
*         NCODE(I)=5: SYME;  NCODE(I)=7: ZERO;  NCODE(I)=20: CYLI.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(i)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
*
*Parameters: output
* L4      total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      used with cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,LX,LY,LZ,L4,KN(LX*LY*LZ*(IELEM+1)**3),MAT(LX*LY*LZ),
     1 IQFR(6*LX*LY*LZ),NCODE(6),ICODE(6),IMPX
      REAL XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),XX(LX*LY*LZ),YY(LX*LY*LZ),
     1 ZZ(LX*LY*LZ),DD(LX*LY*LZ),QFR(6*LX*LY*LZ),VOL(LX*LY*LZ),ZCODE(6)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL LL1,LL2
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP,IWRK
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IP((1+IELEM*LX)*(1+IELEM*LY)*(1+IELEM*LZ)),
     1 IWRK((1+IELEM*LX)*(1+IELEM*LY)*(1+IELEM*LZ)))
*
      IF(IMPX.GT.0) WRITE(6,500) LX,LY,LZ
      MAXIP=(1+IELEM*LX)*(1+IELEM*LY)*(1+IELEM*LZ)
      LC=1+IELEM
      LL=LC*LC*LC
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE ELEMENTS
*----
      IX=LX*(LC-1)+1
      IXY=(LY*(LC-1)+1)*IX
      IXYZ=(LZ*(LC-1)+1)*IXY
      NUM1=0
      NUM2=0
      KEL=0
      DO 182 K0=1,LZ
      DO 181 K1=1,LY
      DO 180 K2=1,LX
      KEL=KEL+1
      XX(KEL)=0.0
      YY(KEL)=0.0
      ZZ(KEL)=0.0
      VOL(KEL)=0.0
      IF(MAT(KEL).LE.0) GO TO 180
      XX(KEL)=XXX(K2+1)-XXX(K2)
      YY(KEL)=YYY(K1+1)-YYY(K1)
      ZZ(KEL)=ZZZ(K0+1)-ZZZ(K0)
      IF(CYLIND) DD(KEL)=0.5*(XXX(K2)+XXX(K2+1))
      IND1=(LC-1)*((K0-1)*IXY+(K1-1)*IX+(K2-1))
      L=0
      DO 12 I=1,LC
      DO 11 J=1,LC
      DO 10 K=1,LC
      L=L+1
      KN(NUM1+L)=IND1+(I-1)*IXY+(J-1)*IX+K
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
      DO 20 IC=1,6
      QFR(NUM2+IC)=0.0
      IQFR(NUM2+IC)=0
   20 CONTINUE
      KK1=KEL-1
      KK2=KEL+1
      KK3=KEL-LX
      KK4=KEL+LX
      KK5=KEL-LX*LY
      KK6=KEL+LX*LY
      FRX=1.0
      FRY=1.0
      FRZ=1.0
*----
*  VOID, REFL OR ZERO BOUNDARY CONTITION
*----
      IF(K2.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK1).EQ.0)
      ENDIF
      IF(LL1) THEN
         IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM2+1)=ALB(ZCODE(1))
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM2+1)=1.0
            IQFR(NUM2+1)=ICODE(1)
         ELSE IF(NCODE(1).EQ.7) THEN
            L=0
            DO 32 I=1,LC
            DO 31 J=1,LC
            DO 30 K=1,LC
            L=L+1
            IF(K.EQ.1) KN(NUM1+L)=0
   30       CONTINUE
   31       CONTINUE
   32       CONTINUE
         ENDIF
      ENDIF
*
      IF(K2.EQ.LX) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK2).EQ.0)
      ENDIF
      IF(LL1) THEN
         IF((NCODE(2).EQ.1).AND.(ICODE(2).EQ.0)) THEN
            QFR(NUM2+2)=ALB(ZCODE(2))
         ELSE IF(NCODE(2).EQ.1) THEN
            QFR(NUM2+2)=1.0
            IQFR(NUM2+2)=ICODE(2)
         ELSE IF(NCODE(2).EQ.7) THEN
            L=0
            DO 42 I=1,LC
            DO 41 J=1,LC
            DO 40 K=1,LC
            L=L+1
            IF(K.EQ.LC) KN(NUM1+L)=0
   40       CONTINUE
   41       CONTINUE
   42       CONTINUE
         ENDIF
      ENDIF
*
      IF(K1.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK3).EQ.0)
      ENDIF
      IF(LL1) THEN
         IF((NCODE(3).EQ.1).AND.(ICODE(3).EQ.0)) THEN
            QFR(NUM2+3)=ALB(ZCODE(3))
         ELSE IF(NCODE(3).EQ.1) THEN
            QFR(NUM2+3)=1.0
            IQFR(NUM2+3)=ICODE(3)
         ELSE IF(NCODE(3).EQ.7) THEN
            L=0
            DO 52 I=1,LC
            DO 51 J=1,LC
            DO 50 K=1,LC
            L=L+1
            IF(J.EQ.1) KN(NUM1+L)=0
   50       CONTINUE
   51       CONTINUE
   52       CONTINUE
         ENDIF
      ENDIF
*
      IF(K1.EQ.LY) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK4).EQ.0)
      ENDIF
      IF(LL1) THEN
         IF((NCODE(4).EQ.1).AND.(ICODE(4).EQ.0)) THEN
            QFR(NUM2+4)=ALB(ZCODE(4))
         ELSE IF(NCODE(4).EQ.1) THEN
            QFR(NUM2+4)=1.0
            IQFR(NUM2+4)=ICODE(4)
         ELSE IF(NCODE(4).EQ.7) THEN
            L=0
            DO 62 I=1,LC
            DO 61 J=1,LC
            DO 60 K=1,LC
            L=L+1
            IF(J.EQ.LC) KN(NUM1+L)=0
   60       CONTINUE
   61       CONTINUE
   62       CONTINUE
         ENDIF
      ENDIF
*
      IF(K0.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK5).EQ.0)
      ENDIF
      IF(LL1) THEN
         IF((NCODE(5).EQ.1).AND.(ICODE(5).EQ.0)) THEN
            QFR(NUM2+5)=ALB(ZCODE(5))
         ELSE IF(NCODE(5).EQ.1) THEN
            QFR(NUM2+5)=1.0
            IQFR(NUM2+5)=ICODE(5)
         ELSE IF(NCODE(5).EQ.7) THEN
            L=0
            DO 72 I=1,LC
            DO 71 J=1,LC
            DO 70 K=1,LC
            L=L+1
            IF(I.EQ.1) KN(NUM1+L)=0
   70       CONTINUE
   71       CONTINUE
   72       CONTINUE
         ENDIF
      ENDIF
*
      IF(K0.EQ.LZ) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK6).EQ.0)
      ENDIF
      IF(LL1) THEN
         IF((NCODE(6).EQ.1).AND.(ICODE(6).EQ.0)) THEN
            QFR(NUM2+6)=ALB(ZCODE(6))
         ELSE IF(NCODE(6).EQ.1) THEN
            QFR(NUM2+6)=1.0
            IQFR(NUM2+6)=ICODE(6)
         ELSE IF(NCODE(6).EQ.7) THEN
            L=0
            DO 82 I=1,LC
            DO 81 J=1,LC
            DO 80 K=1,LC
            L=L+1
            IF(I.EQ.LC) KN(NUM1+L)=0
   80       CONTINUE
   81       CONTINUE
   82       CONTINUE
         ENDIF
      ENDIF
*----
*  TRAN BOUNDARY CONDITION
*----
      IF((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
         DO 91 I=1,LC
         DO 90 J=1,LC
         M=(I-1)*LC*LC+(J-1)*LC+LC
         KN(NUM1+M)=KN(NUM1+M)-IX+1
   90    CONTINUE
   91    CONTINUE
      ENDIF
      IF((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
         DO 101 I=1,LC
         DO 100 K=1,LC
         M=(I-1)*LC*LC+(LC-1)*LC+K
         KN(NUM1+M)=KN(NUM1+M)-IXY+IX
  100    CONTINUE
  101    CONTINUE
      ENDIF
      IF((K0.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
         DO 111 J=1,LC
         DO 110 K=1,LC
         M=(LC-1)*LC*LC+(J-1)*LC+K
         KN(NUM1+M)=KN(NUM1+M)-IXYZ+IXY
  110    CONTINUE
  111    CONTINUE
      ENDIF
*----
*  SYME BOUNDARY CONDITION
*----
      IF((NCODE(1).EQ.5).AND.(K2.EQ.1)) THEN
         QFR(NUM2+1)=QFR(NUM2+2)
         IQFR(NUM2+1)=IQFR(NUM2+2)
         FRX=0.5
         DO 122 I=1,LC
         DO 121 J=1,LC
         DO 120 K=1,(LC+1)/2
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(I-1)*LC*LC+(J-1)*LC+(LC-K+1)
         KN(NUM1+L)=KN(NUM1+M)
  120    CONTINUE
  121    CONTINUE
  122    CONTINUE
      ELSE IF((NCODE(2).EQ.5).AND.(K2.EQ.LX)) THEN
         QFR(NUM2+2)=QFR(NUM2+1)
         IQFR(NUM2+2)=IQFR(NUM2+1)
         FRX=0.5
         DO 132 I=1,LC
         DO 131 J=1,LC
         DO 130 K=(LC+2)/2,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(I-1)*LC*LC+(J-1)*LC+(LC-K+1)
         KN(NUM1+L)=KN(NUM1+M)
  130    CONTINUE
  131    CONTINUE
  132    CONTINUE
      ENDIF
      IF((NCODE(3).EQ.5).AND.(K1.EQ.1)) THEN
         QFR(NUM2+3)=QFR(NUM2+4)
         IQFR(NUM2+3)=IQFR(NUM2+4)
         FRY=0.5
         DO 142 I=1,LC
         DO 141 J=1,(LC+1)/2
         DO 140 K=1,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(I-1)*LC*LC+(LC-J)*LC+K
         KN(NUM1+L)=KN(NUM1+M)
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
      ELSE IF((NCODE(4).EQ.5).AND.(K1.EQ.LY)) THEN
         QFR(NUM2+4)=QFR(NUM2+3)
         IQFR(NUM2+4)=IQFR(NUM2+3)
         FRY=0.5
         DO 152 I=1,LC
         DO 151 J=(LC+2)/2,LC
         DO 150 K=1,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(I-1)*LC*LC+(LC-J)*LC+K
         KN(NUM1+L)=KN(NUM1+M)
  150    CONTINUE
  151    CONTINUE
  152    CONTINUE
      ENDIF
      IF((NCODE(5).EQ.5).AND.(K0.EQ.1)) THEN
         QFR(NUM2+5)=QFR(NUM2+6)
         IQFR(NUM2+5)=IQFR(NUM2+6)
         FRZ=0.5
         DO 162 I=1,(LC+1)/2
         DO 161 J=1,LC
         DO 160 K=1,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(LC-I)*LC*LC+(J-1)*LC+K
         KN(NUM1+L)=KN(NUM1+M)
  160    CONTINUE
  161    CONTINUE
  162    CONTINUE
      ELSE IF((NCODE(6).EQ.5).AND.(K0.EQ.LZ)) THEN
         QFR(NUM2+6)=QFR(NUM2+5)
         IQFR(NUM2+6)=IQFR(NUM2+5)
         FRZ=0.5
         DO 172 I=(LC+2)/2,LC
         DO 171 J=1,LC
         DO 170 K=1,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(LC-I)*LC*LC+(J-1)*LC+K
         KN(NUM1+L)=KN(NUM1+M)
  170    CONTINUE
  171    CONTINUE
  172    CONTINUE
      ENDIF
*
      VOL0=XX(KEL)*YY(KEL)*ZZ(KEL)*FRX*FRY*FRZ
      IF(CYLIND) VOL0=6.2831853072*DD(KEL)*VOL0
      VOL(KEL)=VOL0
      QFR(NUM2+1)=QFR(NUM2+1)*VOL0/XX(KEL)
      QFR(NUM2+2)=QFR(NUM2+2)*VOL0/XX(KEL)
      QFR(NUM2+3)=QFR(NUM2+3)*VOL0/YY(KEL)
      QFR(NUM2+4)=QFR(NUM2+4)*VOL0/YY(KEL)
      QFR(NUM2+5)=QFR(NUM2+5)*VOL0/ZZ(KEL)
      QFR(NUM2+6)=QFR(NUM2+6)*VOL0/ZZ(KEL)
      NUM1=NUM1+LL
      NUM2=NUM2+6
  180 CONTINUE
  181 CONTINUE
  182 CONTINUE
* END OF THE MAIN LOOP OVER ELEMENTS.
*
*----
*  PROCESSING OF 1-D AND 2-D CASES
*----
      LL1=(LX.EQ.1).AND.(NCODE(1).EQ.2).AND.(NCODE(2).EQ.5)
     1 .AND.(IELEM.GT.1)
      LL2=(LX.EQ.1).AND.(NCODE(1).EQ.5).AND.(NCODE(2).EQ.2)
     1 .AND.(IELEM.GT.1)
      IF(LL1.OR.LL2) THEN
         NUM1=0
         DO 200 KEL=1,LX*LY*LZ
         IF(MAT(KEL).EQ.0) GO TO 200
         DO 192 I=1,LC
         DO 191 J=1,LC
         DO 190 K=2,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(I-1)*LC*LC+(J-1)*LC+1
         KN(NUM1+L)=KN(NUM1+M)
  190    CONTINUE
  191    CONTINUE
  192    CONTINUE
         NUM1=NUM1+LL
  200    CONTINUE
      ENDIF
      LL1=(LY.EQ.1).AND.(NCODE(3).EQ.2).AND.(NCODE(4).EQ.5)
     1 .AND.(IELEM.GT.1)
      LL2=(LY.EQ.1).AND.(NCODE(3).EQ.5).AND.(NCODE(4).EQ.2)
     1 .AND.(IELEM.GT.1)
      IF(LL1.OR.LL2) THEN
         NUM1=0
         DO 220 KEL=1,LX*LY*LZ
         IF(MAT(KEL).EQ.0) GO TO 220
         DO 212 I=1,LC
         DO 211 J=2,LC
         DO 210 K=1,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(I-1)*LC*LC+K
         KN(NUM1+L)=KN(NUM1+M)
  210    CONTINUE
  211    CONTINUE
  212    CONTINUE
         NUM1=NUM1+LL
  220    CONTINUE
      ENDIF
      LL1=(LZ.EQ.1).AND.(NCODE(5).EQ.2).AND.(NCODE(6).EQ.5)
     1 .AND.(IELEM.GT.1)
      LL2=(LZ.EQ.1).AND.(NCODE(5).EQ.5).AND.(NCODE(6).EQ.2)
     1 .AND.(IELEM.GT.1)
      IF(LL1.OR.LL2) THEN
         NUM1=0
         DO 240 KEL=1,LX*LY*LZ
         IF(MAT(KEL).EQ.0) GO TO 240
         DO 232 I=2,LC
         DO 231 J=1,LC
         DO 230 K=1,LC
         L=(I-1)*LC*LC+(J-1)*LC+K
         M=(J-1)*LC+K
         KN(NUM1+L)=KN(NUM1+M)
  230    CONTINUE
  231    CONTINUE
  232    CONTINUE
         NUM1=NUM1+LL
  240    CONTINUE
      ENDIF
*----
*  JUXTAPOSITION OF A CHECKERBOARD OVER THE REACTOR DOMAIN
*----
      LZTOT=LZ*(LC-1)+1
      LYTOT=LY*(LC-1)+1
      LXTOT=LX*(LC-1)+1
      DO 250 I=1,LXTOT*LYTOT*LZTOT
      IWRK(I)=-1
  250 CONTINUE
      NUM1=0
      KEL=0
      DO 272 K0=1,LZ
      LK0=(K0-1)*(LC-1)
      DO 271 K1=1,LY
      LK1=(K1-1)*(LC-1)
      DO 270 K2=1,LX
      KEL=KEL+1
      IF(MAT(KEL).EQ.0) GO TO 270
      LK2=(K2-1)*(LC-1)
      L=0
      DO 262 IK0=LK0+1,LK0+LC
      I0=(IK0-1)*LXTOT*LYTOT
      DO 261 IK1=LK1+1,LK1+LC
      I1=I0+(IK1-1)*LXTOT
      DO 260 IK2=LK2+1,LK2+LC
      I2=I1+IK2
      L=L+1
      IND1=KN(NUM1+L)
      IF(IND1.EQ.0) THEN
         IWRK(I2)=0
         GO TO 260
      ENDIF
      IF(IWRK(I2).EQ.-1) THEN
         IWRK(I2)=IND1
      ELSE IF(IWRK(I2).EQ.0) THEN
         KN(NUM1+L)=0
      ELSE IF(IWRK(I2).NE.IND1) THEN
         CALL XABORT('TRIPKN: FAILURE OF THE RENUMBERING ALGORITHM(1).')
      ENDIF
  260 CONTINUE
  261 CONTINUE
  262 CONTINUE
      NUM1=NUM1+LL
  270 CONTINUE
  271 CONTINUE
  272 CONTINUE
*----
*  CALCULATION OF PERMUTATION VECTOR IP AND RENUMBERING OF UNKNOWNS
*----
      DO 280 I=1,MAXIP
      IP(I)=0
  280 CONTINUE
      L4=0
      IF(NCODE(1).EQ.5) THEN
         K2MIN=1+LC/2
      ELSE
         K2MIN=1
      ENDIF
      DO 292 K0=1,LZTOT
      IK0=(K0-1)*LXTOT*LYTOT
      DO 291 K1=1,LYTOT
      IK1=IK0+(K1-1)*LXTOT
      DO 290 K2=K2MIN,LXTOT
      I=IWRK(IK1+K2)
      IF(I.LE.0) GO TO 290
      IF(I.GT.MAXIP) THEN
         CALL XABORT('TRIPKN: FAILURE OF THE RENUMBERING ALGORITHM(2).')
      ENDIF
      IF(IP(I).EQ.0) THEN
         L4=L4+1
         IP(I)=L4
      ENDIF
  290 CONTINUE
  291 CONTINUE
  292 CONTINUE
      DO 300 K=1,NUM1
      KNK=KN(K)
      IF(KNK.NE.0) KN(K)=IP(KNK)
  300 CONTINUE
      IF(IMPX.GT.0) WRITE (6,510) L4
      IF(IMPX.GT.2) WRITE (6,520) (VOL(I),I=1,LX*LY*LZ)
      IF(L4.EQ.0) THEN
         CALL XABORT('TRIPKN: FAILURE OF THE RENUMBERING ALGORITHM(3).')
      ENDIF
*
      IF(IMPX.LT.2) RETURN
      IF(IELEM.EQ.1) THEN
         WRITE (6,530)
         NUM1=0
         NUM2=0
         DO 310 KEL=1,LX*LY*LZ
         IF(MAT(KEL).LE.0) GO TO 310
         WRITE (6,540) KEL,(KN(NUM1+I),I=1,LL),(QFR(NUM2+I),I=1,6)
         NUM1=NUM1+LL
         NUM2=NUM2+6
  310    CONTINUE
      ELSE
         WRITE (6,590)
         NUM1=0
         DO 320 KEL=1,LX*LY*LZ
         IF(MAT(KEL).LE.0) GO TO 320
         WRITE (6,600) KEL,(KN(NUM1+I),I=1,LL)
         NUM1=NUM1+LL
  320    CONTINUE
         WRITE (6,610)
         NUM2=0
         DO 330 KEL=1,LX*LY*LZ
         IF(MAT(KEL).LE.0) GO TO 330
         WRITE (6,620) KEL,(QFR(NUM2+I),I=1,6)
         NUM2=NUM2+6
  330    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IP,IWRK)
      RETURN
*
  500 FORMAT(/38H TRIPKN: PRIMAL FINITE ELEMENT METHOD.//7H NUMBER,
     1 27H OF ELEMENTS ALONG X AXIS =,I3/20X,14HALONG Y AXIS =,I3/
     2 20X,14HALONG Z AXIS =,I3)
  510 FORMAT(31H NUMBER OF UNKNOWNS PER GROUP =,I8)
  520 FORMAT(/20H VOLUMES PER ELEMENT/(1X,1P,10E13.4))
  530 FORMAT(/22H NUMBERING OF UNKNOWNS//8H ELEMENT,5X,7HNUMBERS,
     1 41X,23HVOID BOUNDARY CONDITION)
  540 FORMAT(1X,I6,2X,8I6,2X,1P,6E11.2)
  590 FORMAT(/22H NUMBERING OF UNKNOWNS//5H ELE-/5H MENT,3X,
     1 7HNUMBERS)
  600 FORMAT(1X,I6,2X,20I6/(9X,20I6))
  610 FORMAT(///24H VOID BOUNDARY CONDITION//8H ELEMENT,5X,3HQFR)
  620 FORMAT(1X,I6,4X,1P,6E11.2)
      END
