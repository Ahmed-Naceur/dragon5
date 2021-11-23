*DECK TRIDKN
      SUBROUTINE TRIDKN(IMPX,LX,LY,LZ,CYLIND,IELEM,L4,LL4F,LL4X,LL4Y,
     1 LL4Z,NCODE,ICODE,ZCODE,MAT,VOL,XXX,YYY,ZZZ,XX,YY,ZZ,DD,KN,QFR,
     2 IQFR,IDL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a Thomas-Raviart (dual) formulation of the
* finite element discretization in a 3-D Cartesian geometry.
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
* LL4F    number of flux unknowns.
* LL4X    number of X-directed currents
* LL4Y    number of Y-directed currents
* LL4Z    number of Z-directed currents
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      used with cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
* IDL     position of integrated fluxes into unknown vector.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,LX,LY,LZ,IELEM,L4,LL4F,LL4X,LL4Y,LL4Z,NCODE(6),
     1 ICODE(6),MAT(LX*LY*LZ),KN(LX*LY*LZ*(1+6*IELEM**2)),
     2 IQFR(6*LX*LY*LZ),IDL(LX*LY*LZ)
      REAL ZCODE(6),VOL(LX*LY*LZ),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),
     1 XX(LX*LY*LZ),YY(LX*LY*LZ),ZZ(LX*LY*LZ),DD(LX*LY*LZ),
     2 QFR(6*LX*LY*LZ)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL COND,LL1
      REAL ZALB(6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: IP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IP((LX+1)*LY*LZ*IELEM*IELEM + LX*(LY+1)*LZ*IELEM*IELEM
     1  + LX*LY*(LZ+1)*IELEM*IELEM + LX*LY*LZ*IELEM*IELEM*IELEM))
*----
*  IDENTIFICATION OF THE GEOMETRY. MAIN LOOP OVER THE ELEMENTS
*----
      DO 10 I=1,6
      IF(ZCODE(I).NE.1.0) THEN
         ZALB(I)=2.0*(1.0+ZCODE(I))/(1.0-ZCODE(I))
      ELSE
         ZALB(I)=1.0E20
      ENDIF
   10 CONTINUE
      IF(IMPX.GT.0) WRITE(6,700) LX,LY,LZ
      L2=LX*LY*LZ
      CALL XDISET(KN,L2*(1+6*IELEM**2),0)
      LL4F0=LX*LY*LZ*IELEM**3
      LL4X0=(LX+1)*LY*LZ*IELEM**2
      LL4Y0=LX*(LY+1)*LZ*IELEM**2
      LL4Z0=LX*LY*(LZ+1)*IELEM**2
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
      IF(MAT(KEL).EQ.0) GO TO 180
      XX(KEL)=XXX(K2+1)-XXX(K2)
      YY(KEL)=YYY(K1+1)-YYY(K1)
      ZZ(KEL)=ZZZ(K0+1)-ZZZ(K0)
      IF(CYLIND) DD(KEL)=0.5*(XXX(K2)+XXX(K2+1))
      KN(NUM1+1)=((K0-1)*LX*LY+(K1-1)*LX+K2-1)*IELEM**3 + 1
      DO 20 IEL=1,IELEM**2
      KN(NUM1+1+IEL)=LL4F0+((K0-1)*LY+K1-1)*(LX+1)*IELEM**2+(LX+1)*
     1 (IEL-1)+K2
      KN(NUM1+1+IELEM**2+IEL)=KN(NUM1+1+IEL)+1
      KN(NUM1+1+2*IELEM**2+IEL)=LL4F0+LL4X0+((K0-1)*LX+K2-1)*(LY+1)*
     1 IELEM**2+(LY+1)*(IEL-1)+K1
      KN(NUM1+1+3*IELEM**2+IEL)=KN(NUM1+1+2*IELEM**2+IEL)+1
      KN(NUM1+1+4*IELEM**2+IEL)=LL4F0+LL4X0+LL4Y0+((K1-1)*LX+K2-1)*
     1 (LZ+1)*IELEM**2+(LZ+1)*(IEL-1)+K0
      KN(NUM1+1+5*IELEM**2+IEL)=KN(NUM1+1+4*IELEM**2+IEL)+1
   20 CONTINUE
      CALL XDRSET(QFR(NUM2+1),6,0.0)
      CALL XDISET(IQFR(NUM2+1),6,0)
      FRX=1.0
      FRY=1.0
      FRZ=1.0
*----
*  VOID, REFL OR ZERO BOUNDARY CONTITION
*----
      IF(K2.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL-1).EQ.0)
      ENDIF
      IF(LL1) THEN
         COND=(NCODE(1).EQ.2).OR.((NCODE(1).EQ.1).AND.(ZCODE(1).EQ.1.0))
         IF(COND) THEN
            DO 30 IEL=1,IELEM**2
            KN(NUM1+1+IEL)=0
   30       CONTINUE
         ELSE IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
            QFR(NUM2+1)=ZALB(1)
         ELSE IF(NCODE(1).EQ.1) THEN
            QFR(NUM2+1)=1.0
            IQFR(NUM2+1)=ICODE(1)
         ENDIF
      ENDIF
*
      IF(K2.EQ.LX) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL+1).EQ.0)
      ENDIF
      IF(LL1) THEN
         COND=(NCODE(2).EQ.2).OR.((NCODE(2).EQ.1).AND.(ZCODE(2).EQ.1.0))
         IF(COND) THEN
            DO 40 IEL=1,IELEM**2
            KN(NUM1+1+IELEM**2+IEL)=0
   40       CONTINUE
         ELSE IF((NCODE(2).EQ.1).AND.(ICODE(2).EQ.0)) THEN
            QFR(NUM2+2)=ZALB(2)
         ELSE IF(NCODE(2).EQ.1) THEN
            QFR(NUM2+2)=1.0
            IQFR(NUM2+2)=ICODE(2)
         ENDIF
      ENDIF
*
      IF(K1.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL-LX).EQ.0)
      ENDIF
      IF(LL1) THEN
         COND=(NCODE(3).EQ.2).OR.((NCODE(3).EQ.1).AND.(ZCODE(3).EQ.1.0))
         IF(COND) THEN
            DO 50 IEL=1,IELEM**2
            KN(NUM1+1+2*IELEM**2+IEL)=0
   50       CONTINUE
         ELSE IF((NCODE(3).EQ.1).AND.(ICODE(3).EQ.0)) THEN
            QFR(NUM2+3)=ZALB(3)
         ELSE IF(NCODE(3).EQ.1) THEN
            QFR(NUM2+3)=1.0
            IQFR(NUM2+3)=ICODE(3)
         ENDIF
      ENDIF
*
      IF(K1.EQ.LY) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL+LX).EQ.0)
      ENDIF
      IF(LL1) THEN
         COND=(NCODE(4).EQ.2).OR.((NCODE(4).EQ.1).AND.(ZCODE(4).EQ.1.0))
         IF(COND) THEN
            DO 60 IEL=1,IELEM**2
            KN(NUM1+1+3*IELEM**2+IEL)=0
   60       CONTINUE
         ELSE IF((NCODE(4).EQ.1).AND.(ICODE(4).EQ.0)) THEN
            QFR(NUM2+4)=ZALB(4)
         ELSE IF(NCODE(4).EQ.1) THEN
            QFR(NUM2+4)=1.0
            IQFR(NUM2+4)=ICODE(4)
         ENDIF
      ENDIF
*
      IF(K0.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL-LX*LY).EQ.0)
      ENDIF
      IF(LL1) THEN
         COND=(NCODE(5).EQ.2).OR.((NCODE(5).EQ.1).AND.(ZCODE(5).EQ.1.0))
         IF(COND) THEN
            DO 70 IEL=1,IELEM**2
            KN(NUM1+1+4*IELEM**2+IEL)=0
   70       CONTINUE
         ELSE IF((NCODE(5).EQ.1).AND.(ICODE(5).EQ.0)) THEN
            QFR(NUM2+5)=ZALB(5)
         ELSE IF(NCODE(5).EQ.1) THEN
            QFR(NUM2+5)=1.0
            IQFR(NUM2+5)=ICODE(5)
         ENDIF
      ENDIF
*
      IF(K0.EQ.LZ) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL+LX*LY).EQ.0)
      ENDIF
      IF(LL1) THEN
         COND=(NCODE(6).EQ.2).OR.((NCODE(6).EQ.1).AND.(ZCODE(6).EQ.1.0))
         IF(COND) THEN
            DO 80 IEL=1,IELEM**2
            KN(NUM1+1+5*IELEM**2+IEL)=0
   80       CONTINUE
         ELSE IF((NCODE(6).EQ.1).AND.(ICODE(6).EQ.0)) THEN
            QFR(NUM2+6)=ZALB(6)
         ELSE IF(NCODE(6).EQ.1) THEN
            QFR(NUM2+6)=1.0
            IQFR(NUM2+6)=ICODE(6)
         ENDIF
      ENDIF
*----
*  TRAN BOUNDARY CONDITION
*----
      IF((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
         DO 90 IEL=1,IELEM**2
         KN(NUM1+1+IELEM**2+IEL)=KN(NUM1+1+IELEM**2+IEL)-LX
   90    CONTINUE
      ENDIF
      IF((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
         DO 100 IEL=1,IELEM**2
         KN(NUM1+1+3*IELEM**2+IEL)=KN(NUM1+1+3*IELEM**2+IEL)-LY
  100    CONTINUE
      ENDIF
      IF((K0.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
         DO 110 IEL=1,IELEM**2
         KN(NUM1+1+5*IELEM**2+IEL)=KN(NUM1+1+5*IELEM**2+IEL)-LZ
  110    CONTINUE
      ENDIF
*----
*  SYME BOUNDARY CONDITION
*----
      IF((NCODE(1).EQ.5).AND.(K2.EQ.1)) THEN
         QFR(NUM2+1)=QFR(NUM2+2)
         IQFR(NUM2+1)=IQFR(NUM2+2)
         FRX=0.5
         DO 120 IEL=1,IELEM**2
         KN(NUM1+1+IEL)=-KN(NUM1+1+IELEM**2+IEL)
  120    CONTINUE
      ELSE IF((NCODE(2).EQ.5).AND.(K2.EQ.LX)) THEN
         QFR(NUM2+2)=QFR(NUM2+1)
         IQFR(NUM2+2)=IQFR(NUM2+1)
         FRX=0.5
         DO 130 IEL=1,IELEM**2
         KN(NUM1+1+IELEM**2+IEL)=-KN(NUM1+1+IEL)
  130    CONTINUE
      ENDIF
      IF((NCODE(3).EQ.5).AND.(K1.EQ.1)) THEN
         QFR(NUM2+3)=QFR(NUM2+4)
         IQFR(NUM2+3)=IQFR(NUM2+4)
         FRY=0.5
         DO 140 IEL=1,IELEM**2
         KN(NUM1+1+2*IELEM**2+IEL)=-KN(NUM1+1+3*IELEM**2+IEL)
  140    CONTINUE
      ELSE IF((NCODE(4).EQ.5).AND.(K1.EQ.LY)) THEN
         QFR(NUM2+4)=QFR(NUM2+3)
         IQFR(NUM2+4)=IQFR(NUM2+3)
         FRY=0.5
         DO 150 IEL=1,IELEM**2
         KN(NUM1+1+3*IELEM**2+IEL)=-KN(NUM1+1+2*IELEM**2+IEL)
  150    CONTINUE
      ENDIF
      IF((NCODE(5).EQ.5).AND.(K0.EQ.1)) THEN
         QFR(NUM2+5)=QFR(NUM2+6)
         IQFR(NUM2+5)=IQFR(NUM2+6)
         FRZ=0.5
         DO 160 IEL=1,IELEM**2
         KN(NUM1+1+4*IELEM**2+IEL)=-KN(NUM1+1+5*IELEM**2+IEL)
  160    CONTINUE
      ELSE IF((NCODE(6).EQ.5).AND.(K0.EQ.LZ)) THEN
         QFR(NUM2+6)=QFR(NUM2+5)
         IQFR(NUM2+6)=IQFR(NUM2+5)
         FRZ=0.5
         DO 170 IEL=1,IELEM**2
         KN(NUM1+1+5*IELEM**2+IEL)=-KN(NUM1+1+4*IELEM**2+IEL)
  170    CONTINUE
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
      NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  180 CONTINUE
  181 CONTINUE
  182 CONTINUE
* END OF THE MAIN LOOP OVER ELEMENTS.
*
*----
*  REMOVING THE UNUSED UNKNOWNS INDICES FROM KN
*----
      CALL XDISET(IP,LL4F0+LL4X0+LL4Y0+LL4Z0,0)
      DO 190 NUM1=1,L2*(1+6*IELEM**2)
      IF(KN(NUM1).NE.0) IP(ABS(KN(NUM1)))=1
  190 CONTINUE
      LL4F=0
      IND=0
      DO 200 KEL=1,L2
      IF(IP(IND+1).EQ.1) THEN
         DO 195 IEL=1,IELEM**3
         LL4F=LL4F+1
         IP(IND+IEL)=LL4F
  195    CONTINUE
      ENDIF
      IND=IND+IELEM**3
  200 CONTINUE
      LL4X=0
      DO 210 IND=LL4F0+1,LL4F0+LL4X0
      IF(IP(IND).EQ.1) THEN
         LL4X=LL4X+1
         IP(IND)=LL4F+LL4X
      ENDIF
  210 CONTINUE
      LL4Y=0
      DO 220 IND=LL4F0+LL4X0+1,LL4F0+LL4X0+LL4Y0
      IF(IP(IND).EQ.1) THEN
         LL4Y=LL4Y+1
         IP(IND)=LL4F+LL4X+LL4Y
      ENDIF
  220 CONTINUE
      LL4Z=0
      DO 230 IND=LL4F0+LL4X0+LL4Y0+1,LL4F0+LL4X0+LL4Y0+LL4Z0
      IF(IP(IND).EQ.1) THEN
         LL4Z=LL4Z+1
         IP(IND)=LL4F+LL4X+LL4Y+LL4Z
      ENDIF
  230 CONTINUE
      DO 240 NUM1=1,L2*(1+6*IELEM**2)
      IF(KN(NUM1).NE.0) KN(NUM1)=SIGN(IP(ABS(KN(NUM1))),KN(NUM1))
  240 CONTINUE
      L4=LL4F+LL4X+LL4Y+LL4Z
      NUM1=0
      DO 250 KEL=1,L2
      IDL(KEL)=0
      IF(MAT(KEL).EQ.0) GO TO 250
      IDL(KEL)=KN(NUM1+1)
      NUM1=NUM1+1+6*IELEM**2
  250 CONTINUE
*
      IF(IMPX.GT.0) WRITE(6,710) L4
      IF(IMPX.GT.2) THEN
         WRITE(6,720) (VOL(I),I=1,L2)
         NUM1=0
         WRITE (6,730)
         DO 500 K=1,L2
         IF(MAT(K).EQ.0) GO TO 500
         WRITE (6,740) K,KN(NUM1+1),'X',(KN(NUM1+I),I=2,1+2*IELEM**2)
         WRITE (6,750) 'Y',(KN(NUM1+I),I=2+2*IELEM**2,1+4*IELEM**2)
         WRITE (6,750) 'Z',(KN(NUM1+I),I=2+4*IELEM**2,1+6*IELEM**2)
         NUM1=NUM1+1+6*IELEM**2
  500    CONTINUE
         WRITE (6,760)
         NUM2=0
         DO 510 K=1,L2
         IF(MAT(K).EQ.0) GO TO 510
         WRITE (6,770) K,(QFR(NUM2+I),I=1,6)
         NUM2=NUM2+6
  510    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IP)
      RETURN
*
  700 FORMAT(/42H TRIDKN: MIXED-DUAL FINITE ELEMENT METHOD.//7H NUMBER,
     1 27H OF ELEMENTS ALONG X AXIS =,I3/20X,14HALONG Y AXIS =,I3/
     2 20X,14HALONG Z AXIS =,I3)
  710 FORMAT(31H NUMBER OF UNKNOWNS PER GROUP =,I8)
  720 FORMAT(/20H VOLUMES PER ELEMENT/(1X,1P,10E13.4))
  730 FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-)//8H ELEMENT,8X,
     1 4HFLUX,6X,8HCURRENTS,89(1H.))
  740 FORMAT (1X,I6,5X,I8,6X,A1,12I8/(27X,12I8))
  750 FORMAT (26X,A1,12I8/(27X,12I8))
  760 FORMAT(/8H ELEMENT,3X,23HVOID BOUNDARY CONDITION)
  770 FORMAT (1X,I6,5X,1P,6E10.1)
      END
