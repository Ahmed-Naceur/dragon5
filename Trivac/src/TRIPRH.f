*DECK TRIPRH
      SUBROUTINE TRIPRH(ISPLH,IPTRK,LX,LZ,LL4,SIDE,ZZZ,ZZ,KN,QFR,IQFR,
     1 VOL,MAT,NCODE,ICODE,ZCODE,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a mesh corner finite difference or
* Lagrangian finite element discretization of a 3-D hexagonal geometry.
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
* ISPLH   type of hexagonal finite element:
*         =1 for hexagonal element with 6 points;
*         =2 for hexagonal element with 7 points;
*         =3 for triangular element.
* IPTRK   L_TRACK pointer to the tracking information.
* IMPX    print parameter.
* LX      number of elements.
* LZ      number of axial planes.
* NCODE   type of boundary condition applied on each side (I=1: hbc):
*         NCODE(I)=1: VOID;          =2: REFL;        =5: SYME;
*                 =7: ZERO.
* ICODE   physical albedo index on each side of the domain.
* ZCODE   albedo corresponding to boundary condition 'VOID' on each
*         side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each element.
* SIDE    side of the hexagon.
* ZZZ     Z-coordinates of the axial planes.
*
*Parameters: output
* LL4     order of system matrices.
* ZZ      axial width of each element.
* VOL     volume of each element.
* KN      element-ordered unknown list. Dimensionned to LC*LX*LZ
*         where LC= 14 for triangle and 12 for hexagon.
* QFR     element-ordered boundary conditions.
* IQFR    element-ordered physical albedo indices.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER ISPLH,LX,LZ,LL4,KN(*),IQFR(8*LX*LZ),MAT(LX*LZ),NCODE(6),
     1 ICODE(6),IMPX
      REAL SIDE,ZZZ(LZ+1),ZZ(LX*LZ),QFR(8*LX*LZ),VOL(LX*LZ),ZCODE(6)
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*
      IPAR=ISPLH
      IK=0
      IF(ISPLH.EQ.1) THEN
         IK=12
      ELSE IF(ISPLH.EQ.2) THEN
         IK=14
      ELSE
         CALL XABORT('TRIPRH: DISCRETIZATION NOT AVAILABLE.')
      ENDIF
      CALL TRIHEX(IPAR+2,LX,LZ,LL4,MAT,KN,NCODE,IPTRK)
*----
*  COMPUTE BOUNDARY CONDITIONS
*----
      FRZ=1.
      KEL=0
      NUM1=0
      DO 15 KZ=1,LZ
      DO 10 KX=1,LX
         KEL=KEL + 1
         ZZ(KEL)=0.0
         VOL(KEL)=0.0
         IF(MAT(KEL).LE.0) GO TO 10
         ZZ(KEL)=ZZZ(KZ+1) - ZZZ(KZ)
         DO 20 IC=1,6
            QFR(NUM1+IC)=0.0
            IQFR(NUM1+IC)=0
            NV=NEIGHB (KX,IC,9,LX,POIDS)
            IF(NV.GT.LX) THEN
               IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
                  QFR(NUM1+IC)=ALB(ZCODE(1))
               ELSE IF(NCODE(1).EQ.1) THEN
                  QFR(NUM1+IC)=1.0
                  IQFR(NUM1+IC)=ICODE(1)
               ENDIF
            ELSE IF(MAT(NV).LE.0) THEN
               IF((NCODE(1).EQ.1).AND.(ICODE(1).EQ.0)) THEN
                  QFR(NUM1+IC)=ALB(ZCODE(1))
               ELSE IF(NCODE(1).EQ.1) THEN
                  QFR(NUM1+IC)=1.0
                  IQFR(NUM1+IC)=ICODE(1)
               ENDIF
            ENDIF
 20      CONTINUE
         QFR(NUM1+7)=0.0
         QFR(NUM1+8)=0.0
         IF((NCODE(5).EQ.1).AND.(KZ.EQ.1).AND.(ICODE(5).EQ.0)) THEN
            QFR(NUM1+7)=ALB(ZCODE(5))
         ELSE IF((NCODE(5).EQ.1).AND.(KZ.EQ.1)) THEN
            QFR(NUM1+7)=1.0
            IQFR(NUM1+7)=ICODE(5)
         ENDIF
         IF((NCODE(6).EQ.1).AND.(KZ.EQ.LZ).AND.(ICODE(6).EQ.0)) THEN
            QFR(NUM1+8)=ALB(ZCODE(6))
         ELSE IF((NCODE(6).EQ.1).AND.(KZ.EQ.LZ)) THEN
            QFR(NUM1+8)=1.0
            IQFR(NUM1+8)=ICODE(6)
         ENDIF
         IF((NCODE(5).EQ.5).AND.(KZ.EQ.1)) THEN
            QFR(NUM1+7)=QFR(NUM1+8)
            IQFR(NUM1+7)=IQFR(NUM1+8)
            FRZ=0.5
         ELSE IF((NCODE(6).EQ.5).AND.(KZ.EQ.LZ)) THEN
            QFR(NUM1+7)=QFR(NUM1+8)
            IQFR(NUM1+7)=IQFR(NUM1+8)
            FRZ=0.5
         ENDIF
         ZZ(KEL)=ZZ(KEL)*FRZ
*
*        COMPUTE VOLUMES.
         VOL(KEL)=2.59807587*SIDE*SIDE*ZZ(KEL)
*
         DO 30 IC=1,6
         QFR(NUM1+IC)=QFR(NUM1+IC)*SIDE*ZZ(KEL)
 30      CONTINUE
         QFR(NUM1+7)=QFR(NUM1+7)*SIDE*SIDE
         QFR(NUM1+8)=QFR(NUM1+8)*SIDE*SIDE
         NUM1=NUM1+8
 10   CONTINUE
 15   CONTINUE
      IF(IMPX.GT.2) WRITE(6,720) (VOL(I),I=1,LX*LZ)
*
      IF(IMPX.GT.2) THEN
         NUM1=0
         NUM2=0
         WRITE(6,730)
         DO 510 KZ=1,LZ
         WRITE(6,'(/13H PLANE NUMBER,I6)') KZ
         IF(IK.EQ.12) WRITE(6,740)
         IF(IK.EQ.14) WRITE(6,745)
         DO 500 KX=1,LX
         IF(MAT(KX+(KZ-1)*LX).LE.0) GO TO 500
         K=KX+(KZ-1)*LX
         IF(IK.EQ.12)
     >      WRITE(6,750) K,(KN(NUM1+I),I=1,12),(QFR(NUM2+I),I=1,8)
         IF(IK.EQ.14)
     >      WRITE(6,760) K,(KN(NUM1+I),I=1,14),(QFR(NUM2+I),I=1,8)
         NUM1=NUM1+IK
         NUM2=NUM2+8
 500     CONTINUE
 510     CONTINUE
      ENDIF
      RETURN
*
720   FORMAT(/20H VOLUMES PER ELEMENT/(1X,10(1X,E12.5)))
730   FORMAT(/22H NUMBERING OF UNKNOWNS/1X,21(1H-))
740   FORMAT(/8H ELEMENT,3X,7HNUMBERS,58X,23HVOID BOUNDARY CONDITION)
745   FORMAT(/8H ELEMENT,3X,7HNUMBERS,68X,23HVOID BOUNDARY CONDITION)
750   FORMAT (3X,I4,4X,12I5,3X,8F6.2)
760   FORMAT (3X,I4,4X,14I5,3X,8F6.2)
      END
