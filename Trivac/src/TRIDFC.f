*DECK TRIDFC
      SUBROUTINE TRIDFC(IMPX,LX,LY,LZ,CYLIND,NCODE,ICODE,ZCODE,MAT,XXX,
     1 YYY,ZZZ,LL4,VOL,XX,YY,ZZ,DD,KN,QFR,IQFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh centered finite difference (CHEBY
* type) or nodal collocation discretization in a 3-D geometry.
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
* LL4     total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      used with cylindrical geometry.
* KN      element-ordered unknown list:
*         .GT.0; neighbour index;
*         =-1;   void/albedo boundary condition;
*         =-2;   reflection boundary condition;
*         =-3;   ZERO flux boundary condition;
*         =-4;   SYME boundary condition (axial symmetry).
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
*
*-----------------------------------------------------------------------
*
      INTEGER IMPX,LX,LY,LZ,NCODE(6),ICODE(6),MAT(LX*LY*LZ),LL4,
     1 KN(6*LX*LY*LZ),IQFR(6*LX*LY*LZ)
      REAL ZCODE(6),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),VOL(LX*LY*LZ),
     1 XX(LX*LY*LZ),YY(LX*LY*LZ),ZZ(LX*LY*LZ),DD(LX*LY*LZ),
     2 QFR(6*LX*LY*LZ)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL LL1,LALB
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE ELEMENTS
*----
      IF(IMPX.GT.0) WRITE(6,700) LX,LY,LZ
      LXY=LX*LY
      NUM1=0
      KEL=0
      DO 22 K0=1,LZ
      DO 21 K1=1,LY
      DO 20 K2=1,LX
      KEL=KEL+1
      XX(KEL)=0.0
      YY(KEL)=0.0
      ZZ(KEL)=0.0
      VOL(KEL)=0.0
      IF(MAT(KEL).LE.0) GO TO 20
      XX(KEL)=XXX(K2+1)-XXX(K2)
      YY(KEL)=YYY(K1+1)-YYY(K1)
      ZZ(KEL)=ZZZ(K0+1)-ZZZ(K0)
      IF(CYLIND) DD(KEL)=0.5*(XXX(K2)+XXX(K2+1))
      DO 10 IC=1,6
      QFR(NUM1+IC)=0.0
      IQFR(NUM1+IC)=0
   10 CONTINUE
      FRX=1.0
      FRY=1.0
      FRZ=1.0
      KK1=KEL-1
      KK2=KEL+1
      KK3=KEL-LX
      KK4=KEL+LX
      KK5=KEL-LXY
      KK6=KEL+LXY
*----
*  VOID, REFL, ZERO OR CYLI BOUNDARY CONTITION
*----
      IF(K2.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK1).EQ.0)
      ENDIF
      IF(LL1) THEN
         LALB=(NCODE(1).EQ.1).OR.(NCODE(1).EQ.6)
         IF(LALB.AND.(ICODE(1).EQ.0)) THEN
            KK1=-1
            QFR(NUM1+1)=ALB(ZCODE(1))
         ELSE IF(LALB) THEN
            KK1=-1
            QFR(NUM1+1)=1.0
            IQFR(NUM1+1)=ICODE(1)
         ELSE IF(NCODE(1).EQ.2) THEN
            KK1=-2
         ELSE IF(NCODE(1).EQ.7) THEN
            KK1=-3
         ELSE IF(NCODE(1).EQ.20) THEN
            KK1=-1
         ENDIF
      ENDIF
*
      IF(K2.EQ.LX) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK2).EQ.0)
      ENDIF
      IF(LL1) THEN
         LALB=(NCODE(2).EQ.1).OR.(NCODE(2).EQ.6)
         IF(LALB.AND.(ICODE(2).EQ.0)) THEN
            KK2=-1
            QFR(NUM1+2)=ALB(ZCODE(2))
         ELSE IF(LALB) THEN
            KK2=-1
            QFR(NUM1+2)=1.0
            IQFR(NUM1+2)=ICODE(2)
         ELSE IF(NCODE(2).EQ.2) THEN
            KK2=-2
         ELSE IF(NCODE(2).EQ.7) THEN
            KK2=-3
         ELSE IF(NCODE(2).EQ.20) THEN
            KK2=-1
         ENDIF
      ENDIF
*
      IF(K1.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK3).EQ.0)
      ENDIF
      IF(LL1) THEN
         LALB=(NCODE(3).EQ.1).OR.(NCODE(3).EQ.6)
         IF(LALB.AND.(ICODE(3).EQ.0)) THEN
            KK3=-1
            QFR(NUM1+3)=ALB(ZCODE(3))
         ELSE IF(LALB) THEN
            KK3=-1
            QFR(NUM1+3)=1.0
            IQFR(NUM1+3)=ICODE(3)
         ELSE IF(NCODE(3).EQ.2) THEN
            KK3=-2
         ELSE IF(NCODE(3).EQ.7) THEN
            KK3=-3
         ELSE IF(NCODE(3).EQ.20) THEN
            KK3=-1
         ENDIF
      ENDIF
*
      IF(K1.EQ.LY) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK4).EQ.0)
      ENDIF
      IF(LL1) THEN
         LALB=(NCODE(4).EQ.1).OR.(NCODE(4).EQ.6)
         IF(LALB.AND.(ICODE(4).EQ.0)) THEN
            KK4=-1
            QFR(NUM1+4)=ALB(ZCODE(4))
         ELSE IF(LALB) THEN
            KK4=-1
            QFR(NUM1+4)=1.0
            IQFR(NUM1+4)=ICODE(4)
         ELSE IF(NCODE(4).EQ.2) THEN
            KK4=-2
         ELSE IF(NCODE(4).EQ.7) THEN
            KK4=-3
         ELSE IF(NCODE(4).EQ.20) THEN
            KK4=-1
         ENDIF
      ENDIF
*
      IF(K0.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK5).EQ.0)
      ENDIF
      IF(LL1) THEN
         LALB=(NCODE(5).EQ.1).OR.(NCODE(5).EQ.6)
         IF(LALB.AND.(ICODE(5).EQ.0)) THEN
            KK5=-1
            QFR(NUM1+5)=ALB(ZCODE(5))
         ELSE IF(LALB) THEN
            KK5=-1
            QFR(NUM1+5)=1.0
            IQFR(NUM1+5)=ICODE(5)
         ELSE IF(NCODE(5).EQ.2) THEN
            KK5=-2
         ELSE IF(NCODE(5).EQ.7) THEN
            KK5=-3
         ELSE IF(NCODE(5).EQ.20) THEN
            KK5=-1
         ENDIF
      ENDIF
*
      IF(K0.EQ.LZ) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KK6).EQ.0)
      ENDIF
      IF(LL1) THEN
         LALB=(NCODE(6).EQ.1).OR.(NCODE(6).EQ.6)
         IF(LALB.AND.(ICODE(6).EQ.0)) THEN
            KK6=-1
            QFR(NUM1+6)=ALB(ZCODE(6))
         ELSE IF(LALB) THEN
            KK6=-1
            QFR(NUM1+6)=1.0
            IQFR(NUM1+6)=ICODE(6)
         ELSE IF(NCODE(6).EQ.2) THEN
            KK6=-2
         ELSE IF(NCODE(6).EQ.7) THEN
            KK6=-3
         ELSE IF(NCODE(6).EQ.20) THEN
            KK6=-1
         ENDIF
      ENDIF
*----
*  TRAN BOUNDARY CONDITION
*----
      IF((K2.EQ.1).AND.(NCODE(1).EQ.4)) THEN
         KK1=KEL+LX-1
      ENDIF
      IF((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
         KK2=KEL+1-LX
      ENDIF
      IF((K1.EQ.1).AND.(NCODE(3).EQ.4)) THEN
         KK3=KEL+(LY-1)*LX
      ENDIF
      IF((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
         KK4=KEL-(LY-1)*LX
      ENDIF
      IF((K0.EQ.1).AND.(NCODE(5).EQ.4)) THEN
         KK5=KEL+(LZ-1)*LXY
      ENDIF
      IF((K0.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
         KK6=KEL-(LZ-1)*LXY
      ENDIF
*----
*  SYME BOUNDARY CONDITION
*----
      IF((NCODE(1).EQ.5).AND.(K2.EQ.1)) THEN
         KK1=-4
         FRX=0.5
      ELSE IF((NCODE(2).EQ.5).AND.(K2.EQ.LX)) THEN
         KK2=-4
         FRX=0.5
      ENDIF
      IF((NCODE(3).EQ.5).AND.(K1.EQ.1)) THEN
         KK3=-4
         FRY=0.5
      ELSE IF((NCODE(4).EQ.5).AND.(K1.EQ.LY)) THEN
         KK4=-4
         FRY=0.5
      ENDIF
      IF((NCODE(5).EQ.5).AND.(K0.EQ.1)) THEN
         KK5=-4
         FRZ=0.5
      ELSE IF((NCODE(6).EQ.5).AND.(K0.EQ.LZ)) THEN
         KK6=-4
         FRZ=0.5
      ENDIF
*
      VOL0=XX(KEL)*YY(KEL)*ZZ(KEL)*FRX*FRY*FRZ
      IF(CYLIND) VOL0=6.2831853072*DD(KEL)*VOL0
      VOL(KEL)=VOL0
      KN(NUM1+1)=KK1
      KN(NUM1+2)=KK2
      KN(NUM1+3)=KK3
      KN(NUM1+4)=KK4
      KN(NUM1+5)=KK5
      KN(NUM1+6)=KK6
      NUM1=NUM1+6
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
* END OF THE MAIN LOOP OVER ELEMENTS.
*
      LL4=0
      DO 40 KEL=1,LXY*LZ
      IF(MAT(KEL).NE.0) LL4=LL4+1
   40 CONTINUE
*
      IF(IMPX.GE.2) THEN
         WRITE(6,720) (VOL(I),I=1,LXY*LZ)
         WRITE(6,750)
         NUM1=0
         DO 50 KEL=1,LXY*LZ
         IF(MAT(KEL).LE.0) GO TO 50
         WRITE (6,760) KEL,(KN(NUM1+I),I=1,6),(QFR(NUM1+I),I=1,6)
         NUM1=NUM1+6
   50    CONTINUE
      ENDIF
      RETURN
*
  700 FORMAT(/53H TRIDFC: MESH CENTERED FINITE DIFFERENCE OR NODAL COL,
     1 16HLOCATION METHOD.//34H NUMBER OF ELEMENTS ALONG X AXIS =,I3/
     2 20X,14HALONG Y AXIS =,I3/20X,14HALONG Z AXIS =,I3)
  720 FORMAT(/20H VOLUMES PER ELEMENT/(1X,1P,10E13.4))
  750 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//
     1 8H ELEMENT,5X,7HNUMBERS,50X,23HVOID BOUNDARY CONDITION)
  760 FORMAT(1X,I6,7X,6I8,6X,6F9.2)
      END
