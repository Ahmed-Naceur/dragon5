*DECK BIVPKN
      SUBROUTINE BIVPKN (MAXEV,IMPX,LX,LY,CYLIND,IELEM,L4,NCODE,ICODE,
     1 ZCODE,MAT,VOL,XXX,YYY,XX,YY,DD,KN,QFR,IQFR,BFR,MU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh corner finite difference or primal
* finite element discretization in a 2-D geometry. This version does
* not support diagonal symmetries.
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
* MAXEV   allocated storage for vector MU.
* IMPX    print parameter.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry)
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NCODE   type of boundary condition applied on each side
*         (i=1: X-  i=2: X+  i=3: Y-  i=4: Y+):
*         NCODE(I)=1: VOID;   NCODE(I)=2: REFL;   NCODE(I)=4: TRAN;
*         NCODE(I)=5: SYME;   NCODE(I)=7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on
*         each side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
*
*Parameters: output
* L4      total number of unknown (variational coefficients) per
*         energy group (order of system matrices).
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* DD      values used with a cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* BFR     element-ordered surface fractions.
* MU      compressed storage mode indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXEV,IMPX,LX,LY,IELEM,L4,NCODE(4),ICODE(4),MAT(LX*LY),
     1 KN(LX*LY*IELEM*IELEM),IQFR(4*LX*LY),MU(MAXEV)
      REAL ZCODE(4),VOL(LX*LY),XXX(LX+1),YYY(LY+1),XX(LX*LY),YY(LX*LY),
     1 DD(LX*LY),QFR(4*LX*LY),BFR(4*LX*LY)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL LOG1,LOG2,LOG3,LOG4
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP,IWRK
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IP((IELEM*LX+1)*(IELEM*LY+1)))
      ALLOCATE(IWRK((IELEM*LX+1)*(IELEM*LY+1)))
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE ELEMENTS.
*----
      IF(IMPX.GT.0) WRITE(6,700) LX,LY
      LC=1+IELEM
      LL=LC*LC
      IX=LX*(LC-1)+1
      IXY=(LY*(LC-1)+1)*IX
      SURFTOT=0.0
      NUM1=0
      NUM2=0
      KEL=0
      DO 151 K1=1,LY
      DO 150 K2=1,LX
      KEL=KEL+1
      XX(KEL)=0.0
      YY(KEL)=0.0
      VOL(KEL)=0.0
      IF(MAT(KEL).LE.0) GO TO 150
      XX(KEL)=XXX(K2+1)-XXX(K2)
      YY(KEL)=YYY(K1+1)-YYY(K1)
      IF(CYLIND) DD(KEL)=0.5*(XXX(K2)+XXX(K2+1))
      IND1=(LC-1)*((K1-1)*IX+(K2-1))
      L=0
      DO 15 I=1,LC
      DO 10 J=1,LC
      L=L+1
      KN(NUM1+L)=IND1+(I-1)*IX+J
   10 CONTINUE
   15 CONTINUE
      DO 20 IC=1,4
      QFR(NUM2+IC)=0.0
      IQFR(NUM2+IC)=0
      BFR(NUM2+IC)=0.0
   20 CONTINUE
      KK1=KEL-1
      KK2=KEL+1
      KK3=KEL-LX
      KK4=KEL+LX
      FRX=1.0
      FRY=1.0
*----
*  VOID, REFL OR ZERO BOUNDARY CONDITION.
*----
      IF(K2.EQ.1) THEN
         LOG1=.TRUE.
      ELSE
         LOG1=(MAT(KK1).EQ.0)
      ENDIF
      IF(LOG1) THEN
         IF(NCODE(1).EQ.1) THEN
            IF(ICODE(1).EQ.0) THEN
              QFR(NUM2+1)=ALB(ZCODE(1))
            ELSE
              QFR(NUM2+1)=1.0
              IQFR(NUM2+1)=ICODE(1)
            ENDIF
         ELSE IF(NCODE(1).EQ.7) THEN
            L=0
            DO 35 I=1,LC
            DO 30 J=1,LC
            L=L+1
            IF(J.EQ.1) KN(NUM1+L)=0
   30       CONTINUE
   35       CONTINUE
         ENDIF
      ENDIF
*
      IF(K2.EQ.LX) THEN
         LOG2=.TRUE.
      ELSE
         LOG2=(MAT(KK2).EQ.0)
      ENDIF
      IF(LOG2) THEN
         IF(NCODE(2).EQ.1) THEN
            IF(ICODE(2).EQ.0) THEN
              QFR(NUM2+2)=ALB(ZCODE(2))
            ELSE
              QFR(NUM2+2)=1.0
              IQFR(NUM2+2)=ICODE(2)
            ENDIF
         ELSE IF(NCODE(2).EQ.7) THEN
            L=0
            DO 45 I=1,LC
            DO 40 J=1,LC
            L=L+1
            IF(J.EQ.LC) KN(NUM1+L)=0
   40       CONTINUE
   45       CONTINUE
         ENDIF
      ENDIF
*
      IF(K1.EQ.1) THEN
         LOG3=.TRUE.
      ELSE
         LOG3=(MAT(KK3).EQ.0)
      ENDIF
      IF(LOG3) THEN
         IF(NCODE(3).EQ.1) THEN
            IF(ICODE(3).EQ.0) THEN
              QFR(NUM2+3)=ALB(ZCODE(3))
            ELSE
              QFR(NUM2+3)=1.0
              IQFR(NUM2+3)=ICODE(3)
            ENDIF
         ELSE IF(NCODE(3).EQ.7) THEN
            L=0
            DO 55 I=1,LC
            DO 50 J=1,LC
            L=L+1
            IF(I.EQ.1) KN(NUM1+L)=0
   50       CONTINUE
   55       CONTINUE
         ENDIF
      ENDIF
*
      IF(K1.EQ.LY) THEN
         LOG4=.TRUE.
      ELSE
         LOG4=(MAT(KK4).EQ.0)
      ENDIF
      IF(LOG4) THEN
         IF(NCODE(4).EQ.1) THEN
            IF(ICODE(4).EQ.0) THEN
              QFR(NUM2+4)=ALB(ZCODE(4))
            ELSE
              QFR(NUM2+4)=1.0
              IQFR(NUM2+4)=ICODE(4)
            ENDIF
         ELSE IF(NCODE(4).EQ.7) THEN
            L=0
            DO 65 I=1,LC
            DO 60 J=1,LC
            L=L+1
            IF(I.EQ.LC) KN(NUM1+L)=0
   60       CONTINUE
   65       CONTINUE
         ENDIF
      ENDIF
*----
*  TRAN BOUNDARY CONDITION.
*----
      IF((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
         DO 70 I=1,LC
         M=(I-1)*LC+LC
         KN(NUM1+M)=KN(NUM1+M)-IX+1
   70    CONTINUE
      ENDIF
      IF((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
         DO 80 I=1,LC
         M=(LC-1)*LC+I
         KN(NUM1+M)=KN(NUM1+M)-IXY+IX
   80    CONTINUE
      ENDIF
*----
*  SYME BOUNDARY CONDITION.
*----
      IF((NCODE(1).EQ.5).AND.(K2.EQ.1)) THEN
         QFR(NUM2+1)=QFR(NUM2+2)
         IQFR(NUM2+1)=IQFR(NUM2+2)
         FRX=0.5
         DO 95 I=1,LC
         DO 90 J=1,(LC+1)/2
         L=(I-1)*LC+J
         M=(I-1)*LC+(LC-J+1)
         KN(NUM1+L)=KN(NUM1+M)
   90    CONTINUE
   95    CONTINUE
      ELSE IF((NCODE(2).EQ.5).AND.(K2.EQ.LX)) THEN
         QFR(NUM2+2)=QFR(NUM2+1)
         IQFR(NUM2+2)=IQFR(NUM2+1)
         FRX=0.5
         DO 105 I=1,LC
         DO 100 J=(LC+2)/2,LC
         L=(I-1)*LC+J
         M=(I-1)*LC+(LC-J+1)
         KN(NUM1+L)=KN(NUM1+M)
  100    CONTINUE
  105    CONTINUE
      ENDIF
      IF((NCODE(3).EQ.5).AND.(K1.EQ.1)) THEN
         QFR(NUM2+3)=QFR(NUM2+4)
         IQFR(NUM2+3)=IQFR(NUM2+4)
         FRY=0.5
         DO 115 I=1,(LC+1)/2
         DO 110 J=1,LC
         L=(I-1)*LC+J
         M=(LC-I)*LC+J
         KN(NUM1+L)=KN(NUM1+M)
  110    CONTINUE
  115    CONTINUE
      ELSE IF((NCODE(4).EQ.5).AND.(K1.EQ.LY)) THEN
         QFR(NUM2+4)=QFR(NUM2+3)
         IQFR(NUM2+4)=IQFR(NUM2+3)
         FRY=0.5
         DO 125 I=(LC+2)/2,LC
         DO 120 J=1,LC
         L=(I-1)*LC+J
         M=(LC-I)*LC+J
         KN(NUM1+L)=KN(NUM1+M)
  120    CONTINUE
  125    CONTINUE
      ENDIF
*
      VOL0=XX(KEL)*YY(KEL)*FRX*FRY
      IF(CYLIND) THEN
         VOL0=6.2831853072*DD(KEL)*VOL0
      ENDIF
      VOL(KEL)=VOL0
      QFR(NUM2+1)=QFR(NUM2+1)*VOL0/XX(KEL)
      QFR(NUM2+2)=QFR(NUM2+2)*VOL0/XX(KEL)
      QFR(NUM2+3)=QFR(NUM2+3)*VOL0/YY(KEL)
      QFR(NUM2+4)=QFR(NUM2+4)*VOL0/YY(KEL)
*
      IF(((NCODE(1).EQ.1).OR.(NCODE(1).EQ.7)).AND.LOG1)
     1 BFR(NUM2+1)=VOL0/XX(KEL)
      IF(((NCODE(2).EQ.1).OR.(NCODE(2).EQ.7)).AND.LOG2)
     1 BFR(NUM2+2)=VOL0/XX(KEL)
      IF(((NCODE(3).EQ.1).OR.(NCODE(3).EQ.7)).AND.LOG3)
     1 BFR(NUM2+3)=VOL0/YY(KEL)
      IF(((NCODE(4).EQ.1).OR.(NCODE(4).EQ.7)).AND.LOG4)
     1 BFR(NUM2+4)=VOL0/YY(KEL)
      SURFTOT=SURFTOT+BFR(NUM2+1)+BFR(NUM2+2)+BFR(NUM2+3)+BFR(NUM2+4)
      NUM1=NUM1+LL
      NUM2=NUM2+4
  150 CONTINUE
  151 CONTINUE
* END OF THE MAIN LOOP OVER THE ELEMENTS.
*
*----
* COMPUTE THE SURFACE FRACTIONS.
*----
      IF(SURFTOT.GT.0.0) THEN
         DO 155 I=1,4*LX*LY
         BFR(I)=BFR(I)/SURFTOT
  155    CONTINUE
      ENDIF
*----
*  TREATMENT OF 1D CASES.
*----
      LOG1=(LX.EQ.1).AND.(NCODE(1).EQ.2).AND.(NCODE(2).EQ.5)
     1 .AND.(IELEM.GT.1)
      LOG2=(LX.EQ.1).AND.(NCODE(1).EQ.5).AND.(NCODE(2).EQ.2)
     1 .AND.(IELEM.GT.1)
      IF(LOG1.OR.LOG2) THEN
         NUM1=0
         DO 170 KEL=1,LX*LY
         IF(MAT(KEL).EQ.0) GO TO 170
         DO 165 I=1,LC
         DO 160 J=2,LC
         KN(NUM1+(I-1)*LC+J)=KN(NUM1+(I-1)*LC+1)
  160    CONTINUE
  165    CONTINUE
         NUM1=NUM1+LL
  170    CONTINUE
      ENDIF
      LOG1=(LY.EQ.1).AND.(NCODE(3).EQ.2).AND.(NCODE(4).EQ.5)
     1 .AND.(IELEM.GT.1)
      LOG2=(LY.EQ.1).AND.(NCODE(3).EQ.5).AND.(NCODE(4).EQ.2)
     1 .AND.(IELEM.GT.1)
      IF(LOG1.OR.LOG2) THEN
         NUM1=0
         DO 190 KEL=1,LX*LY
         IF(MAT(KEL).EQ.0) GO TO 190
         DO 185 I=2,LC
         DO 180 J=1,LC
         KN(NUM1+(I-1)*LC+J)=KN(NUM1+J)
  180    CONTINUE
  185    CONTINUE
         NUM1=NUM1+LL
  190    CONTINUE
      ENDIF
*----
*  JUXTAPOSITION OF A CHECKERBOARD OVER THE REACTOR DOMAIN.
*----
      LYTOT=LY*(LC-1)+1
      LXTOT=LX*(LC-1)+1
      DO 220 I=1,LXTOT*LYTOT
      IWRK(I)=-1
  220 CONTINUE
      NUM1=0
      KEL=0
      DO 245 K1=1,LY
      LK1=(K1-1)*(LC-1)
      DO 240 K2=1,LX
      KEL=KEL+1
      IF(MAT(KEL).EQ.0) GO TO 240
      LK2=(K2-1)*(LC-1)
      L=0
      DO 235 IK1=LK1+1,LK1+LC
      I1=(IK1-1)*LXTOT
      DO 230 IK2=LK2+1,LK2+LC
      I2=I1+IK2
      L=L+1
      IND1=KN(NUM1+L)
      IF(IND1.EQ.0) THEN
         IWRK(I2)=0
         GO TO 230
      ENDIF
      IF(IWRK(I2).EQ.-1) THEN
         IWRK(I2)=IND1
      ELSE IF(IWRK(I2).EQ.0) THEN
         KN(NUM1+L)=0
      ELSE IF(IWRK(I2).NE.IND1) THEN
         CALL XABORT('BIVPKN: FAILURE OF THE RENUMBERING ALGORITHM(1).')
      ENDIF
  230 CONTINUE
  235 CONTINUE
      NUM1=NUM1+LL
  240 CONTINUE
  245 CONTINUE
*----
*  COMPUTE THE PERMUTATION VECTOR IP AND RENUMBER THE UNKNOWNS.
*----
      DO 250 I=1,MAXEV
      IP(I)=0
  250 CONTINUE
      L4=0
      IF(NCODE(1).EQ.5) THEN
         K2MIN=1+LC/2
      ELSE
         K2MIN=1
      ENDIF
      DO 265 K1=1,LYTOT
      IK1=(K1-1)*LXTOT
      DO 260 K2=K2MIN,LXTOT
      I=IWRK(IK1+K2)
      IF(I.LE.0) GO TO 260
      IF(I.GT.MAXEV) THEN
         CALL XABORT('BIVPKN: FAILURE OF THE RENUMBERING ALGORITHM(2).')
      ENDIF
      IF(IP(I).EQ.0) THEN
         L4=L4+1
         IP(I)=L4
      ENDIF
  260 CONTINUE
  265 CONTINUE
      DO 270 K=1,NUM1
      KNK=KN(K)
      IF(KNK.NE.0) KN(K)=IP(KNK)
  270 CONTINUE
      IF(L4.EQ.0) THEN
         CALL XABORT('BIVPKN: FAILURE OF THE RENUMBERING ALGORITHM(3).')
      ELSE IF(L4.GT.MAXEV) THEN
         CALL XABORT('BIVPKN: INSUFFICIENT MAXEV.')
      ENDIF
      IF(IMPX.GT.2) WRITE (6,745) (VOL(I),I=1,LX*LY)
*----
*  COMPUTE THE SYSTEM MATRIX BANDWIDTH.
*----
      DO 450 I=1,L4
      MU(I)=0
  450 CONTINUE
      NUM1=0
      DO 480 K=1,LX*LY
      IF(MAT(K).LE.0) GO TO 480
      DO 470 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 470
      DO 460 J=1,LL
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 460
      MU(IND1)=MAX0(MU(IND1),IND1-IND2+1)
  460 CONTINUE
  470 CONTINUE
      NUM1=NUM1+LL
  480 CONTINUE
      IIMAX=0
      DO 490 I=1,L4
      IIMAX=IIMAX+MU(I)
      MU(I)=IIMAX
  490 CONTINUE
*
      IF(IMPX.GT.2) THEN
         WRITE (6,720) IIMAX
         NUM1=0
         NUM2=0
         IF(IELEM.EQ.1) THEN
            WRITE (6,750)
            DO 500 K=1,LX*LY
            IF(MAT(K).LE.0) GO TO 500
            WRITE (6,755) K,(KN(NUM1+I),I=1,LL),(QFR(NUM2+I),I=1,4),
     1      (BFR(NUM2+I),I=1,4)
            NUM1=NUM1+LL
            NUM2=NUM2+4
500         CONTINUE
         ELSE IF(IELEM.EQ.2) THEN
            WRITE (6,760)
            DO 510 K=1,LX*LY
            IF(MAT(K).LE.0) GO TO 510
            WRITE (6,765) K,(KN(NUM1+I),I=1,LL),(QFR(NUM2+I),I=1,4)
            NUM1=NUM1+LL
            NUM2=NUM2+4
510         CONTINUE
            NUM2=0
            WRITE (6,830)
            DO 515 K=1,LX*LY
            IF(MAT(K).LE.0) GO TO 515
            WRITE (6,820) K,(BFR(NUM2+I),I=1,4)
            NUM2=NUM2+4
515         CONTINUE
         ELSE IF((IELEM.EQ.3).OR.(IELEM.EQ.4)) THEN
            WRITE (6,790)
            DO 530 K=1,LX*LY
            IF(MAT(K).LE.0) GO TO 530
            WRITE (6,800) K,(KN(NUM1+I),I=1,LL)
            NUM1=NUM1+LL
530         CONTINUE
            WRITE (6,810)
            DO 540 K=1,LX*LY
            IF(MAT(K).LE.0) GO TO 540
            WRITE (6,820) K,(QFR(NUM2+I),I=1,4)
            NUM2=NUM2+4
540         CONTINUE
            NUM2=0
            WRITE (6,830)
            DO 550 K=1,LX*LY
            IF(MAT(K).LE.0) GO TO 550
            WRITE (6,820) K,(BFR(NUM2+I),I=1,4)
            NUM2=NUM2+4
550         CONTINUE
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IWRK,IP)
      RETURN
*
  700 FORMAT(/38H BIVPKN: PRIMAL FINITE ELEMENT METHOD.//7H NUMBER,
     1 27H OF ELEMENTS ALONG X AXIS =,I3/26H NUMBER OF ELEMENTS ALONG ,
     2 8HY AXIS =,I3)
  720 FORMAT(/52H NUMBER OF TERMS IN THE COMPRESSED SYSTEM MATRICES =,
     1 I7)
  745 FORMAT(/20H VOLUMES PER ELEMENT/(1X,1P,10E13.4))
  750 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//
     1 8H ELEMENT,5X,7HNUMBERS,23X,23HVOID BOUNDARY CONDITION,25X,
     2 17HSURFACE FRACTIONS)
  755 FORMAT (3X,I4,7X,4I5,6X,1P,4E11.2,5X,4E10.2)
  760 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//
     1 8H ELEMENT,5X,7HNUMBERS,47X,23HVOID BOUNDARY CONDITION)
  765 FORMAT (3X,I4,7X,9I5,6X,1P,4E11.2)
  790 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//5H ELE-/5H MENT,
     1 3X,7HNUMBERS)
  800 FORMAT (1X,I4,2X,25I5)
  810 FORMAT (///24H VOID BOUNDARY CONDITION//8H ELEMENT,5X,3HQFR)
  820 FORMAT (3X,I4,4X,1P,4E10.2)
  830 FORMAT (///17H SURFACE FRACTION//8H ELEMENT,5X,3HBFR)
      END
