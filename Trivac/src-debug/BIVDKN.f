*DECK BIVDKN
      SUBROUTINE BIVDKN (MAXEV,IMPX,LX,LY,CYLIND,IELEM,ICOL,L4,NCODE,
     1 ICODE,ZCODE,MAT,VOL,XXX,YYY,XX,YY,DD,KN,QFR,IQFR,BFR,IDL,MU)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mixed-dual formulation of the finite-
* element discretization in a 2-D geometry. This version does not
* support diagonal symmetries.
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
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NCODE   type of boundary condition applied on each side
*         (i=1: X-  i=2: X+  i=3: Y-  i=4: Y+):
*         NCODE(I)=1: VOID;   NCODE(I)=2: REFL;   NCODE(I)=4: TRAN;
*         NCODE(I)=5: SYME;   NCODE(I)=7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
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
* DD      value used with a cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* BFR     element-ordered surface fractions.
* IDL     position of integrated fluxes into unknown vector.
* MU      compressed storage mode indices.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXEV,IMPX,LX,LY,IELEM,ICOL,L4,NCODE(4),ICODE(4),
     1 MAT(LX*LY),KN(5*LX*LY),IQFR(4*LX*LY),IDL(LX*LY),MU(MAXEV)
      REAL ZCODE(4),VOL(LX*LY),XXX(LX+1),YYY(LY+1),XX(LX*LY),YY(LX*LY),
     1 DD(LX*LY),QFR(4*LX*LY),BFR(4*LX*LY)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL COND,LOG1,LOG2,LOG3,LOG4
      CHARACTER TEXT8*8
      REAL ZALB(4)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IP(MAXEV))
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE ELEMENTS.
*----
      DO 10 I=1,4
      IF(ZCODE(I).NE.1.0) THEN
         ZALB(I)=2.0*(1.0+ZCODE(I))/(1.0-ZCODE(I))
      ELSE
         ZALB(I)=1.0E20
      ENDIF
   10 CONTINUE
      IF(IMPX.GT.0) WRITE(6,700) LX,LY
      CALL XDISET(KN,5*LX*LY,0)
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
      IF(MAT(KEL).EQ.0) GO TO 150
      XX(KEL)=XXX(K2+1)-XXX(K2)
      YY(KEL)=YYY(K1+1)-YYY(K1)
      IF(CYLIND) DD(KEL)=0.5*(XXX(K2)+XXX(K2+1))
      IND1=(K1-1)*(3*LX+1)
      KN(NUM1+1)=IND1+LX+2*K2
      KN(NUM1+2)=IND1+LX+2*K2-1
      KN(NUM1+3)=IND1+LX+2*K2+1
      KN(NUM1+4)=IND1+K2
      KN(NUM1+5)=IND1+3*LX+K2+1
      CALL XDRSET(QFR(NUM2+1),4,0.0)
      CALL XDISET(IQFR(NUM2+1),4,0)
      CALL XDRSET(BFR(NUM2+1),4,0.0)
      FRX=1.0
      FRY=1.0
*----
*  VOID, REFL OR ZERO BOUNDARY CONTITION.
*----
      IF(K2.EQ.1) THEN
         LOG1=.TRUE.
      ELSE
         LOG1=(MAT(KEL-1).EQ.0)
      ENDIF
      IF(LOG1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            KN(NUM1+2)=0
         ELSE IF(NCODE(1).EQ.1) THEN
            IF(ICODE(1).EQ.0) THEN
              QFR(NUM2+1)=ZALB(1)
            ELSE
              QFR(NUM2+1)=1.0
              IQFR(NUM2+1)=ICODE(1)
            ENDIF
         ENDIF
      ENDIF
*
      IF(K2.EQ.LX) THEN
         LOG2=.TRUE.
      ELSE
         LOG2=(MAT(KEL+1).EQ.0)
      ENDIF
      IF(LOG2) THEN
         COND=(NCODE(2).EQ.2).OR.((NCODE(2).EQ.1).AND.(ZCODE(2).EQ.1.0))
         IF(COND) THEN
            KN(NUM1+3)=0
         ELSE IF(NCODE(2).EQ.1) THEN
            IF(ICODE(2).EQ.0) THEN
              QFR(NUM2+2)=ZALB(2)
            ELSE
              QFR(NUM2+2)=1.0
              IQFR(NUM2+2)=ICODE(2)
            ENDIF
         ENDIF
      ENDIF
*
      IF(K1.EQ.1) THEN
         LOG3=.TRUE.
      ELSE
         LOG3=(MAT(KEL-LX).EQ.0)
      ENDIF
      IF(LOG3) THEN
         COND=(NCODE(3).EQ.2).OR.((NCODE(3).EQ.1).AND.(ZCODE(3).EQ.1.0))
         IF(COND) THEN
            KN(NUM1+4)=0
         ELSE IF(NCODE(3).EQ.1) THEN
            IF(ICODE(3).EQ.0) THEN
              QFR(NUM2+3)=ZALB(3)
            ELSE
              QFR(NUM2+3)=1.0
              IQFR(NUM2+3)=ICODE(3)
            ENDIF
         ENDIF
      ENDIF
*
      IF(K1.EQ.LY) THEN
         LOG4=.TRUE.
      ELSE
         LOG4=(MAT(KEL+LX).EQ.0)
      ENDIF
      IF(LOG4) THEN
         COND=(NCODE(4).EQ.2).OR.((NCODE(4).EQ.1).AND.(ZCODE(4).EQ.1.0))
         IF(COND) THEN
            KN(NUM1+5)=0
         ELSE IF(NCODE(4).EQ.1) THEN
            IF(ICODE(4).EQ.0) THEN
              QFR(NUM2+4)=ZALB(4)
            ELSE
              QFR(NUM2+4)=1.0
              IQFR(NUM2+4)=ICODE(4)
            ENDIF
         ENDIF
      ENDIF
*----
*  TRAN BOUNDARY CONDITION.
*----
      IF((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
         KN(NUM1+3)=KN(NUM1+3)-2*LX
      ENDIF
      IF((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
         KN(NUM1+5)=K2
      ENDIF
*----
*  SYME BOUNDARY CONDITION.
*----
      IF((NCODE(1).EQ.5).AND.(K2.EQ.1)) THEN
         QFR(NUM2+1)=QFR(NUM2+2)
         IQFR(NUM2+1)=IQFR(NUM2+2)
         FRX=0.5
         KN(NUM1+2)=-KN(NUM1+3)
      ELSE IF((NCODE(2).EQ.5).AND.(K2.EQ.LX)) THEN
         QFR(NUM2+2)=QFR(NUM2+1)
         IQFR(NUM2+2)=IQFR(NUM2+1)
         FRX=0.5
         KN(NUM1+3)=-KN(NUM1+2)
      ENDIF
      IF((NCODE(3).EQ.5).AND.(K1.EQ.1)) THEN
         QFR(NUM2+3)=QFR(NUM2+4)
         FRY=0.5
         KN(NUM1+4)=-KN(NUM1+5)
      ELSE IF((NCODE(4).EQ.5).AND.(K1.EQ.LY)) THEN
         QFR(NUM2+4)=QFR(NUM2+3)
         IQFR(NUM2+4)=IQFR(NUM2+3)
         FRY=0.5
         KN(NUM1+5)=-KN(NUM1+4)
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
      NUM1=NUM1+5
      NUM2=NUM2+4
  150 CONTINUE
  151 CONTINUE
* END OF THE MAIN LOOP OVER ELEMENTS.
*
* COMPUTE THE SURFACE FRACTIONS.
      IF(SURFTOT.GT.0.0) THEN
         DO 155 I=1,4*LX*LY
         BFR(I)=BFR(I)/SURFTOT
  155    CONTINUE
      ENDIF
*----
*  REMOVING THE UNUSED UNKNOWNS INDICES FROM KN.
*----
      LL4=LY*(3*LX+1)+LX
      WRITE (TEXT8,'(I8)') LL4
      IF(LL4.GT.MAXEV) CALL XABORT('BIVDKN: MAXEV SHOULD BE INCREASED '
     1 //'TO'//TEXT8//'.')
      DO 160 IND=1,LL4
      IP(IND)=0
  160 CONTINUE
      DO 170 NUM1=1,5*LX*LY
      IF(KN(NUM1).NE.0) IP(ABS(KN(NUM1)))=1
  170 CONTINUE
      L4=0
      DO 180 IND=1,LL4
      IF(IP(IND).EQ.1) THEN
         L4=L4+1
         IP(IND)=L4
      ENDIF
  180 CONTINUE
      DO 190 NUM1=1,5*LX*LY
      IF(KN(NUM1).NE.0) KN(NUM1)=SIGN(IP(ABS(KN(NUM1))),KN(NUM1))
  190 CONTINUE
*----
*  PROCESS CASES WITH IELEM.GT.1.
*----
      IF(IELEM.GT.1) THEN
         LL4=0
         DO 220 IND=1,L4
         IP(IND)=LL4+1
         NUM1=0
         DO 210 KEL=1,LX*LY
         IF(MAT(KEL).EQ.0) GO TO 210
         IF(ABS(KN(NUM1+1)).EQ.IND) THEN
            LL4=LL4+IELEM**2
            GO TO 220
         ELSE
            DO 200 I=2,5
            IF(ABS(KN(NUM1+I)).EQ.IND) THEN
               LL4=LL4+IELEM
               GO TO 220
            ENDIF
  200       CONTINUE
         ENDIF
         NUM1=NUM1+5
  210    CONTINUE
         CALL XABORT('BIVDKN: FAILURE OF THE RENUMBERING ALGORITHM.')
  220    CONTINUE
         L4=LL4
         DO 230 NUM1=1,5*LX*LY
         IF(KN(NUM1).NE.0) KN(NUM1)=SIGN(IP(ABS(KN(NUM1))),KN(NUM1))
  230    CONTINUE
      ENDIF
      NUM1=0
      DO 235 KEL=1,LX*LY
      IDL(KEL)=0
      IF(MAT(KEL).EQ.0) GO TO 235
      IDL(KEL)=KN(NUM1+1)
      NUM1=NUM1+5
  235 CONTINUE
      WRITE (TEXT8,'(I8)') L4
      IF(L4.GT.MAXEV) CALL XABORT('BIVDKN: MAXEV SHOULD BE INCREASED TO'
     1 //TEXT8//'.')
      IF(IMPX.GT.2) WRITE(6,710) (VOL(I),I=1,LX*LY)
*----
*  COMPUTE THE SYSTEM MATRIX BANDWIDTH.
*----
      DO 240 I=1,L4
      MU(I)=1
  240 CONTINUE
      NUM1=0
      DO 270 KEL=1,LX*LY
      IF(MAT(KEL).EQ.0) GO TO 270
      DO 260 I0=1,IELEM
      INX1=ABS(KN(NUM1+2))+I0-1
      INX2=ABS(KN(NUM1+3))+I0-1
      INY1=ABS(KN(NUM1+4))+I0-1
      INY2=ABS(KN(NUM1+5))+I0-1
      DO 250 J0=1,IELEM
      JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
      IF(IELEM.GE.4) MU(JND1)=MAX(MU(JND1),J0)
      IF(KN(NUM1+2).NE.0) THEN
         MU(JND1)=MAX(MU(JND1),JND1-INX1+1)
         MU(INX1)=MAX(MU(INX1),INX1-JND1+1)
      ENDIF
      IF(KN(NUM1+3).NE.0) THEN
         MU(INX2)=MAX(MU(INX2),INX2-JND1+1)
         MU(JND1)=MAX(MU(JND1),JND1-INX2+1)
      ENDIF
      JND1=KN(NUM1+1)+(J0-1)*IELEM+I0-1
      IF(IELEM.GE.4) MU(JND1)=MAX(MU(JND1),(J0-1)*IELEM+1)
      IF(KN(NUM1+4).NE.0) THEN
         MU(JND1)=MAX(MU(JND1),JND1-INY1+1)
         MU(INY1)=MAX(MU(INY1),INY1-JND1+1)
      ENDIF
      IF(KN(NUM1+5).NE.0) THEN
         MU(INY2)=MAX(MU(INY2),INY2-JND1+1)
         MU(JND1)=MAX(MU(JND1),JND1-INY2+1)
      ENDIF
  250 CONTINUE
      IF(ICOL.NE.2) THEN
         IF((KN(NUM1+2).NE.0).AND.(KN(NUM1+3).NE.0)) THEN
            MU(INX2)=MAX(MU(INX2),INX2-INX1+1)
            MU(INX1)=MAX(MU(INX1),INX1-INX2+1)
         ENDIF
         IF((KN(NUM1+4).NE.0).AND.(KN(NUM1+5).NE.0)) THEN
            MU(INY2)=MAX(MU(INY2),INY2-INY1+1)
            MU(INY1)=MAX(MU(INY1),INY1-INY2+1)
         ENDIF
      ENDIF
  260 CONTINUE
      NUM1=NUM1+5
  270 CONTINUE
      IIMAX=0
      DO 280 I=1,L4
      IIMAX=IIMAX+MU(I)
      MU(I)=IIMAX
  280 CONTINUE
*
      IF(IMPX.GT.2) THEN
         WRITE (6,720) IIMAX
         NUM1=0
         NUM2=0
         WRITE (6,750)
         DO 500 K=1,LX*LY
         IF(MAT(K).EQ.0) GO TO 500
         WRITE (6,755) K,(KN(NUM1+I),I=1,5),(QFR(NUM2+I),I=1,4),
     1   (BFR(NUM2+I),I=1,4)
         NUM1=NUM1+5
         NUM2=NUM2+4
  500    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IP)
      RETURN
*
  700 FORMAT(/42H BIVDKN: MIXED-DUAL FINITE ELEMENT METHOD.//7H NUMBER,
     1 27H OF ELEMENTS ALONG X AXIS =,I3/26H NUMBER OF ELEMENTS ALONG ,
     2 8HY AXIS =,I3)
  710 FORMAT(/20H VOLUMES PER ELEMENT/(1X,1P,10E13.4))
  720 FORMAT(/52H NUMBER OF TERMS IN THE COMPRESSED SYSTEM MATRICES =,
     1 I7)
  750 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//
     1 8H ELEMENT,2X,7HNUMBERS,30X,23HVOID BOUNDARY CONDITION,23X,
     2 17HSURFACE FRACTIONS)
  755 FORMAT (1X,I4,2X,5I7,2X,1P,4E11.2,3X,4E10.2)
      END
