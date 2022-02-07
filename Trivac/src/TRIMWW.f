*DECK TRIMWW
      SUBROUTINE TRIMWW(IR,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,KN,QFR,MUW,
     1 IPW,IPR,A11W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a mesh centered finite difference
* discretization in hexagonal geometry (complete hexagons).
* Note: system matrices should be initialized by the calling program.
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
* SGD     nuclear properties per material mixtures:
*         SGD(L,1)= W-, X-, and Y-oriented diffusion coefficients;
*         SGD(L,3)= Z-oriented diffusion coefficients;
*         SGD(L,4)= removal macroscopic cross section.
* XSGD    nuclear properties (IPR=0), derivatives (IPR=1) or first
*         variations (IPR=2 or 3) of nuclear properties per material
*         mixture.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* MUW     W-oriented compressed storage mode indices.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IPW     W-oriented permutation matrices.
* IPX     X-oriented permutation matrices.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
* IPR     type of assembly:
*         =0: compute the system matrices;
*         =1: compute the derivative of system matrices;
*         =2 or =3: compute the variation of system matrices.
*
*Parameters: output
* A11W    W-oriented matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUW(LL4).
* A11X    X-oriented matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUX(LL4).
* A11Y    Y-oriented matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUY(LL4).
* A11Z    Z-oriented matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUZ(LL4).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(8*NEL),MUW(LL4),IPW(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11W(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(8),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  ASSEMBLY OF MATRIX A11W
*----
      ALLOCATE(IGAR(LL4))
      LL=0
      DO 10 K=1,NEL
      IF(MAT(K).LE.0) GO TO 10
      LL=LL+1
      IGAR(LL)=K
   10 CONTINUE
      NUM1=0
      KEL=0
      DO 70 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 70
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 60
      KEL=KEL+1
*
      CALL TRIHCO (IR,K,NEL,VOL0,MAT,SGD(1,1),XSGD(1,1),SIDE,ZZ,
     1 KN(NUM1+1),QFR(NUM1+1),IGAR,IPR,A1)
      KK1=KN(NUM1+6)
      KK2=KN(NUM1+3)
*
      INW1=IPW(KEL)
      KEY0=MUW(INW1)-INW1
      IF(KK1.GT.0) THEN
         INW2=IPW(KK1)
         IF(INW2.LT.INW1) THEN
            KEY=KEY0+INW2
            A11W(KEY)=A11W(KEY)-REAL(A1(6))
         ENDIF
      ENDIF
      IF(KK2.GT.0) THEN
         INW2=IPW(KK2)
         IF(INW2.LT.INW1) THEN
            KEY=KEY0+INW2
            A11W(KEY)=A11W(KEY)-REAL(A1(3))
         ENDIF
      ENDIF
      KEY=KEY0+INW1
      VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)+A1(7)+A1(8)
      A11W(KEY)=A11W(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
   60 NUM1=NUM1+8
   70 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
*
      SUBROUTINE TRIMWX (IR,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,KN,QFR,MUX,
     1              IPX,IPR,A11X)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(8*NEL),MUX(LL4),IPX(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11X(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(8),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  ASSEMBLY OF MATRIX A11X
*----
      ALLOCATE(IGAR(LL4))
      LL=0
      DO 80 K=1,NEL
      IF(MAT(K).LE.0) GO TO 80
      LL=LL+1
      IGAR(LL)=K
   80 CONTINUE
      NUM1=0
      KEL=0
      DO 140 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 140
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 130
      KEL=KEL+1
*
      CALL TRIHCO (IR,K,NEL,VOL0,MAT,SGD(1,1),XSGD(1,1),SIDE,ZZ,
     1 KN(NUM1+1),QFR(NUM1+1),IGAR,IPR,A1)
      KK3=KN(NUM1+1)
      KK4=KN(NUM1+4)
*
      INX1=IPX(KEL)
      KEY0=MUX(INX1)-INX1
      IF(KK3.GT.0) THEN
         INX2=IPX(KK3)
         IF(INX2.LT.INX1) THEN
            KEY=KEY0+INX2
            A11X(KEY)=A11X(KEY)-REAL(A1(1))
         ENDIF
      ENDIF
      IF(KK4.GT.0) THEN
         INX2=IPX(KK4)
         IF(INX2.LT.INX1) THEN
            KEY=KEY0+INX2
            A11X(KEY)=A11X(KEY)-REAL(A1(4))
         ENDIF
      ENDIF
      KEY=KEY0+INX1
      VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)+A1(7)+A1(8)
      A11X(KEY)=A11X(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
  130 NUM1=NUM1+8
  140 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
*
      SUBROUTINE TRIMWY (IR,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,KN,QFR,
     1              MUY,IPY,IPR,A11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(8*NEL),MUY(LL4),IPY(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11Y(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(8),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  ASSEMBLY OF MATRIX A11Y
*----
      ALLOCATE(IGAR(LL4))
      LL=0
      DO 85 K=1,NEL
      IF(MAT(K).LE.0) GO TO 85
      LL=LL+1
      IGAR(LL)=K
   85 CONTINUE
      NUM1=0
      KEL=0
      DO 145 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 145
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 135
      KEL=KEL+1
*
      CALL TRIHCO (IR,K,NEL,VOL0,MAT,SGD(1,1),XSGD(1,1),SIDE,ZZ,
     1 KN(NUM1+1),QFR(NUM1+1),IGAR,IPR,A1)
      KK5=KN(NUM1+2)
      KK6=KN(NUM1+5)
*
      INY1=IPY(KEL)
      KEY0=MUY(INY1)-INY1
      IF(KK5.GT.0) THEN
         INY2=IPY(KK5)
         IF(INY2.LT.INY1) THEN
            KEY=KEY0+INY2
            A11Y(KEY)=A11Y(KEY)-REAL(A1(2))
         ENDIF
      ENDIF
      IF(KK6.GT.0) THEN
         INY2=IPY(KK6)
         IF(INY2.LT.INY1) THEN
            KEY=KEY0+INY2
            A11Y(KEY)=A11Y(KEY)-REAL(A1(5))
         ENDIF
      ENDIF
      KEY=KEY0+INY1
      VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)+A1(7)+A1(8)
      A11Y(KEY)=A11Y(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
  135 NUM1=NUM1+8
  145 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
*
      SUBROUTINE TRIMWZ (IR,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,KN,QFR,
     1              MUZ,IPZ,IPR,A11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,NEL,LL4,MAT(NEL),KN(8*NEL),MUZ(LL4),IPZ(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11Z(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(8),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  ASSEMBLY OF MATRIX A11Z
*----
      ALLOCATE(IGAR(LL4))
      LL=0
      DO 150 K=1,NEL
      IF(MAT(K).LE.0) GO TO 150
      LL=LL+1
      IGAR(LL)=K
  150 CONTINUE
      NUM1=0
      KEL=0
      DO 210 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 210
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 200
      KEL=KEL+1
*
      CALL TRIHCO (IR,K,NEL,VOL0,MAT,SGD(1,1),XSGD(1,1),SIDE,ZZ,
     1 KN(NUM1+1),QFR(NUM1+1),IGAR,IPR,A1)
      KK7=KN(NUM1+7)
      KK8=KN(NUM1+8)
*
      INZ1=IPZ(KEL)
      KEY0=MUZ(INZ1)-INZ1
      IF(KK7.GT.0) THEN
         INZ2=IPZ(KK7)
         IF(INZ2.LT.INZ1) THEN
            KEY=KEY0+INZ2
            A11Z(KEY)=A11Z(KEY)-REAL(A1(7))
         ENDIF
      ENDIF
      IF(KK8.GT.0) THEN
        INZ2=IPZ(KK8)
        IF(INZ2.LT.INZ1) THEN
            KEY=KEY0+INZ2
            A11Z(KEY)=A11Z(KEY)-REAL(A1(8))
         ENDIF
      ENDIF
      KEY=KEY0+INZ1
      VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)+A1(7)+A1(8)
      A11Z(KEY)=A11Z(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
  200 NUM1=NUM1+8
  210 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
