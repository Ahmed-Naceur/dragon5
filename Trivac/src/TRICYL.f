*DECK TRICYL
      SUBROUTINE TRICYL(MAXMIX,IMPX,ICHX,IDIM,LX,LY,LZ,XX,YY,ZZ,VOL,
     1 MAT,NCODE,ZALB,NR0,RR0,XR0,ANG,SGD,QFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the albedo term corresponding to a cylinderized boundary
* in Cartesian geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* MAXMIX  first dimension of matrix SGD.
* IMPX    print parameter (equal to zero for no print).
* ICHX    type of finite element approximation:
*         =1 primal (Lagrangian) finite elements or mesh corner finit
*         differences;
*         =2 dual finite elements;
*         =3 or 4 nodal collocation method or mesh centered finite
*         differences.
* IDIM    number of dimensions.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* VOL     volume of each element.
* MAT     mixture index of each element.
* NCODE   type of boundary condition applied on each side
*         (i=1: X-  i=2: X+  i=3: Y-  i=4: Y+):
*         NCODE(I)=1: VOID;   NCODE(I)=2: REFL;   NCODE(I)=5: SYME;
*         NCODE(I)=7: ZERO;   NCODE(I)=20: VOID on cylindrical boundary.
* ZALB    albedo function corresponding to boundary condition 'VOID' on
*         each side (ZALB(I)=0.0 by default).
* NR0     number of radii.
* RR0     radii.
* XR0     coordinates on principal axis.
* ANG     angles for applying circular correction.
* SGD     directional diffusion coefficients per mixture.
*
*Parameters: output
* QFR     boundary transmission factor.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXMIX,IMPX,ICHX,IDIM,LX,LY,LZ,MAT(LX*LY*LZ),NCODE(6),NR0
      REAL XX(LX*LY*LZ),YY(LX*LY*LZ),ZZ(LX*LY*LZ),VOL(LX*LY*LZ),
     1 ZALB(6),RR0(NR0),XR0(NR0),ANG(NR0),SGD(MAXMIX,3),QFR(6*LX*LY*LZ)
*----
*  LOCAL VARIABLES
*----
      LOGICAL LL1
      CHARACTER*4 CAXE(3)
      REAL CENTER(3),CELEM(3)
      REAL, DIMENSION(:), ALLOCATABLE :: XXX,YYY,ZZZ
      DATA CAXE / '(X) ', '(Y) ', '(Z) ' /
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XXX(LX+1),YYY(LY+1),ZZZ(LZ+1))
*----
*  DETERMINE CARTESIAN COORDINATES
*----
      KEL=0
      ZZZ(1)=0.0
      DO 12 K0=1,LZ
      YYY(1)=0.0
      DO 11 K1=1,LY
      XXX(1)=0.0
      DO 10 K2=1,LX
      KEL=KEL+1
      IF(MAT(KEL).LE.0) GO TO 10
      XXX(K2+1)=XXX(K2)+XX(KEL)
      YYY(K1+1)=YYY(K1)+YY(KEL)
      ZZZ(K0+1)=ZZZ(K0)+ZZ(KEL)
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
*
      CALL TRIKAX (IDIM,NCODE,XXX,YYY,ZZZ,LX,LY,LZ,IAXIS,CENTER)
      IF((IAXIS.GT.0).AND.(IMPX.GT.0)) THEN
         WRITE(6,600) CAXE(IAXIS),
     1   CAXE(MOD(IAXIS  ,3)+1), CENTER(MOD(IAXIS  ,3)+1),
     2   CAXE(MOD(IAXIS+1,3)+1), CENTER(MOD(IAXIS+1,3)+1)
      ENDIF
      IF(NR0.LE.0) CALL XABORT('TRICYL: B.C. RADIUS NOT DEFINED.')
*
      NUM2=0
      KEL=0
      DO 152 K0=1,LZ
      DO 151 K1=1,LY
      DO 150 K2=1,LX
      KEL=KEL+1
      L=MAT(KEL)
      IF(L.LE.0) GO TO 150
*
      IF(K2.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL-1).EQ.0)
      ENDIF
      IF(LL1.AND.(NCODE(1).EQ.20)) THEN
         CELEM(1)=XXX(K2)
         CELEM(2)=0.5*(YYY(K1+1)+YYY(K1))
         CELEM(3)=0.5*(ZZZ(K0+1)+ZZZ(K0))
         CALL TRIZNR(IMPX,1,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,QFRI,
     1   QTRI)
         IF(ICHX.EQ.2) THEN
            QFR(NUM2+1)=(SGD(L,1)*ZALB(1)+QTRI)/(SGD(L,1)*QFRI)
         ELSE
            QFR(NUM2+1)=SGD(L,1)*QFRI*ZALB(1)/(SGD(L,1)+QTRI*ZALB(1))
         ENDIF
         IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
            QFR(NUM2+1)=QFR(NUM2+1)*VOL(KEL)/(XXX(K2+1)-XXX(K2))
         ENDIF
      ENDIF
*
      IF(K2.EQ.LX) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL+1).EQ.0)
      ENDIF
      IF(LL1.AND.(NCODE(2).EQ.20)) THEN
         CELEM(1)=XXX(K2+1)
         CELEM(2)=0.5*(YYY(K1+1)+YYY(K1))
         CELEM(3)=0.5*(ZZZ(K0+1)+ZZZ(K0))
         CALL TRIZNR(IMPX,2,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,QFRI,
     1   QTRI)
         IF(ICHX.EQ.2) THEN
            QFR(NUM2+2)=(SGD(L,1)*ZALB(2)+QTRI)/(SGD(L,1)*QFRI)
         ELSE
            QFR(NUM2+2)=SGD(L,1)*QFRI*ZALB(2)/(SGD(L,1)+QTRI*ZALB(2))
         ENDIF
         IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
            QFR(NUM2+2)=QFR(NUM2+2)*VOL(KEL)/(XXX(K2+1)-XXX(K2))
         ENDIF
      ENDIF
*
      IF(K1.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL-LX).EQ.0)
      ENDIF
      IF(LL1.AND.(NCODE(3).EQ.20)) THEN
         CELEM(1)=0.5*(XXX(K2+1)+XXX(K2))
         CELEM(2)=YYY(K1)
         CELEM(3)=0.5*(ZZZ(K0+1)+ZZZ(K0))
         CALL TRIZNR(IMPX,3,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,QFRI,
     1   QTRI)
         IF(ICHX.EQ.2) THEN
            QFR(NUM2+3)=(SGD(L,2)*ZALB(3)+QTRI)/(SGD(L,2)*QFRI)
         ELSE
            QFR(NUM2+3)=SGD(L,2)*QFRI*ZALB(3)/(SGD(L,2)+QTRI*ZALB(3))
         ENDIF
         IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
            QFR(NUM2+3)=QFR(NUM2+3)*VOL(KEL)/(YYY(K1+1)-YYY(K1))
         ENDIF
      ENDIF
*
      IF(K1.EQ.LY) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL+LX).EQ.0)
      ENDIF
      IF(LL1.AND.(NCODE(4).EQ.20)) THEN
         CELEM(1)=0.5*(XXX(K2+1)+XXX(K2))
         CELEM(2)=YYY(K1+1)
         CELEM(3)=0.5*(ZZZ(K0+1)+ZZZ(K0))
         CALL TRIZNR(IMPX,4,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,QFRI,
     1   QTRI)
         IF(ICHX.EQ.2) THEN
            QFR(NUM2+4)=(SGD(L,2)*ZALB(4)+QTRI)/(SGD(L,2)*QFRI)
         ELSE
            QFR(NUM2+4)=SGD(L,2)*QFRI*ZALB(4)/(SGD(L,2)+QTRI*ZALB(4))
         ENDIF
         IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
            QFR(NUM2+4)=QFR(NUM2+4)*VOL(KEL)/(YYY(K1+1)-YYY(K1))
         ENDIF
      ENDIF
*
      IF(K0.EQ.1) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL-LX*LY).EQ.0)
      ENDIF
      IF(LL1.AND.(NCODE(5).EQ.20)) THEN
         CELEM(1)=0.5*(XXX(K2+1)+XXX(K2))
         CELEM(2)=0.5*(YYY(K1+1)+YYY(K1))
         CELEM(3)=ZZZ(K0)
         CALL TRIZNR(IMPX,5,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,QFRI,
     1   QTRI)
         IF(ICHX.EQ.2) THEN
            QFR(NUM2+5)=(SGD(L,3)*ZALB(5)+QTRI)/(SGD(L,3)*QFRI)
         ELSE
            QFR(NUM2+5)=SGD(L,3)*QFRI*ZALB(5)/(SGD(L,3)+QTRI*ZALB(5))
         ENDIF
         IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
            QFR(NUM2+5)=QFR(NUM2+5)*VOL(KEL)/(ZZZ(K0+1)-ZZZ(K0))
         ENDIF
      ENDIF
*
      IF(K0.EQ.LZ) THEN
         LL1=.TRUE.
      ELSE
         LL1=(MAT(KEL+LX*LY).EQ.0)
      ENDIF
      IF(LL1.AND.(NCODE(6).EQ.20)) THEN
         CELEM(1)=0.5*(XXX(K2+1)+XXX(K2))
         CELEM(2)=0.5*(YYY(K1+1)+YYY(K1))
         CELEM(3)=ZZZ(K0+1)
         CALL TRIZNR(IMPX,6,CENTER,CELEM,IAXIS,NR0,RR0,XR0,ANG,QFRI,
     1   QTRI)
         IF(ICHX.EQ.2) THEN
            QFR(NUM2+6)=(SGD(L,3)*ZALB(6)+QTRI)/(SGD(L,3)*QFRI)
         ELSE
            QFR(NUM2+6)=SGD(L,3)*QFRI*ZALB(6)/(SGD(L,3)+QTRI*ZALB(6))
         ENDIF
         IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
            QFR(NUM2+6)=QFR(NUM2+6)*VOL(KEL)/(ZZZ(K0+1)-ZZZ(K0))
         ENDIF
      ENDIF
*
      IF((NCODE(1).EQ.5).AND.(NCODE(2).EQ.20).AND.(LX.EQ.1)) THEN
         QFR(NUM2+1)=QFR(NUM2+2)
      ELSE IF((NCODE(1).EQ.20).AND.(NCODE(2).EQ.5).AND.(LX.EQ.1)) THEN
         QFR(NUM2+2)=QFR(NUM2+1)
      ENDIF
      IF((NCODE(3).EQ.5).AND.(NCODE(4).EQ.20).AND.(LY.EQ.1)) THEN
         QFR(NUM2+3)=QFR(NUM2+4)
      ELSE IF((NCODE(3).EQ.20).AND.(NCODE(4).EQ.5).AND.(LY.EQ.1)) THEN
         QFR(NUM2+4)=QFR(NUM2+3)
      ENDIF
      IF((NCODE(5).EQ.5).AND.(NCODE(6).EQ.20).AND.(LZ.EQ.1)) THEN
         QFR(NUM2+5)=QFR(NUM2+6)
      ELSE IF((NCODE(5).EQ.20).AND.(NCODE(6).EQ.5).AND.(LZ.EQ.1)) THEN
         QFR(NUM2+6)=QFR(NUM2+5)
      ENDIF
*
      NUM2=NUM2+6
  150 CONTINUE
  151 CONTINUE
  152 CONTINUE
*
      IF(IMPX.GE.2) THEN
         WRITE (6,610)
         NUM2=0
         DO 160 KEL=1,LX*LY*LZ
         IF(MAT(KEL).LE.0) GO TO 160
         WRITE (6,620) KEL,(QFR(NUM2+I),I=1,6)
         NUM2=NUM2+6
  160    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XXX,YYY,ZZZ)
      RETURN
*
  600 FORMAT (/52H TRICYL: CYLINDRICAL ALBEDO BOUNDARY CONDITION ON A ,
     1 17HCYLINDER OF AXIS ,A4/
     2 9X,12HCENTER IS ( ,A4,1H=,1P,E15.7,3H , ,A4,1H=,E15.7 ,1H) )
  610 FORMAT(///53H VOID BOUNDARY CONDITION WITH CYLINDRICAL CORRECTION:
     1 //8H ELEMENT,5X,3HQFR)
  620 FORMAT(1X,I6,4X,1P,6E11.2)
      END
