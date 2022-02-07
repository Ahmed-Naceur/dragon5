*DECK TRIMTW
      SUBROUTINE TRIMTW(ISPLH,IR,NEL,LL4,VOL,MAT,MATN,SGD,XSGD,SIDE,
     1 ZZ,KN,QFR,MUW,IPW,IPR,A11W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a mesh centered finite difference
* discretization in hexagonal geometry (triangular sub meshs).
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
* ISPLH   used to compute the number of triangles as 6*(ISPLH-1)**2.
* IR      first dimension of matrix SGD.
* NEL     total number of finite elements.
* LL4     order of system matrices.
* VOL     volume of each element.
* MAT     mixture index assigned to each hexagon.
* MATN    mixture index assigned to each triangle.
* SGD     nuclear properties per material mixtures:
*         SGD(L,1): W-, X-, and Y-oriented diffusion coefficients;
*         SGD(L,3): Z-oriented diffusion coefficients;
*         SGD(L,4): removal macroscopic cross section.
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
* IPR     type of calculation:
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
      INTEGER ISPLH,IR,NEL,LL4,MAT(NEL),MATN(LL4),
     1 KN((18*(ISPLH-1)**2+3)*NEL),MUW(LL4),IPW(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11W(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(5),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWRK
*----
*  ASSEMBLY OF MATRIX A11W
*----
      NUM1 = 0
      NUM2 = 0
      NTPH = 6 * (ISPLH-1)**2
      NTPL = 1 + 2 * (ISPLH-1)
      NVT1 = NTPL + 2 * (ISPLH-2) + NTPH / 2
      NVT2 = NTPH - NTPL - (ISPLH-4) * (NTPL+2)
      NVT3 = NTPH - (ISPLH-4) * NTPL
      IVAL = 3*NTPH+8
      IF(ISPLH.EQ.3) NVT2 = NTPH
      IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
      IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
      ICR = ISAU*(1+2*(ISPLH-2))
      ALLOCATE(IWRK(NEL))
      MEL = 0
      DO 10 M=1,NEL
         IF(MAT(M).LE.0) GO TO 10
         MEL = MEL + 1
         IWRK(MEL) = M
 10   CONTINUE
      DO 40 K=1,NEL
         L = MAT(K)
         IF(L.EQ.0) GO TO 40
         VOL0 = VOL(K)/NTPH
         IF(VOL0.EQ.0.0) GO TO 30
         KK4=KN(NUM1+3*NTPH+7)
         KK5=KN(NUM1+3*NTPH+8)
         IF(KK4.GT.0) KK4 = IWRK(KK4)
         IF(KK5.GT.0) KK5 = IWRK(KK5)
         DO 20 I = 1,NTPH
*
            CALL TRINEI (3,1,1,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,
     >                   NTPH,NTPL,NVT1,NVT2,NVT3,IVAL,KN)
*
            CALL TRITCO (NEL,LL4,ISPLH,IR,IQF,K,KK1,KK2,KK3,KK4,KK5,
     >      VOL0,MAT,MATN,SGD(1,1),XSGD(1,1),SIDE,ZZ,QFR(NUM2+1),IPR,A1)
*
            INW1=IPW(KEL)
            KEY0=MUW(INW1)-INW1
            IF(KK1.GT.0) THEN
               INW2=IPW(KK1)
               IF(INW2.LT.INW1) THEN
                  KEY=KEY0+INW2
                  A11W(KEY)=A11W(KEY)-REAL(A1(1))/2.
               ENDIF
            ENDIF
            IF(KK2.GT.0) THEN
               INW2=IPW(KK2)
               IF(INW2.LT.INW1) THEN
                  KEY=KEY0+INW2
                  A11W(KEY)=A11W(KEY)-REAL(A1(2))/2.
               ENDIF
            ENDIF
            KEY=KEY0+INW1
            VAR1 = A1(1)+A1(2)+A1(3)+A1(4)+A1(5)
            A11W(KEY)=A11W(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
   20    CONTINUE
   30    NUM1=NUM1+IVAL
         NUM2=NUM2+8
   40 CONTINUE
      DEALLOCATE(IWRK)
      RETURN
      END
*
      SUBROUTINE TRIMTX (ISPLH,IR,NEL,LL4,VOL,MAT,MATN,SGD,XSGD,SIDE,
     1              ZZ,KN,QFR,MUX,IPX,IPR,A11X)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ISPLH,IR,NEL,LL4,MAT(NEL),MATN(LL4),
     1 KN((18*(ISPLH-1)**2+3)*NEL),MUX(LL4),IPX(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11X(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(5),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWRK
*----
*  ASSEMBLY OF MATRIX A11X
*----
      NUM1=0
      NUM2=0
      NTPH = 6*(ISPLH-1)**2
      NTPL = 1+2*(ISPLH-1)
      NVT1 = NTPL + 2 * (ISPLH-2) + NTPH / 2
      NVT2 = NTPH - NTPL - (ISPLH-4) * (NTPL+2)
      NVT3 = NTPH - (ISPLH-4) * NTPL
      IVAL = 3*NTPH+8
      IF(ISPLH.EQ.3) NVT2 = NTPH
      IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
      IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
      ICR = ISAU*(1+2*(ISPLH-2))
      ALLOCATE(IWRK(NEL))
      MEL = 0
      DO 105 M=1,NEL
         IF(MAT(M).LE.0) GO TO 105
         MEL = MEL + 1
         IWRK(MEL) = M
105   CONTINUE
      DO 130 K=1,NEL
         L = MAT(K)
         IF(L.EQ.0) GO TO 130
         VOL0 = VOL(K)/NTPH
         IF(VOL0.EQ.0.0) GO TO 120
         KK4=KN(NUM1+3*NTPH+7)
         KK5=KN(NUM1+3*NTPH+8)
         IF(KK4.GT.0) KK4 = IWRK(KK4)
         IF(KK5.GT.0) KK5 = IWRK(KK5)
         DO 110 I = 1,NTPH
*
            CALL TRINEI (3,2,1,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,
     >                   NTPH,NTPL,NVT1,NVT2,NVT3,IVAL,KN)
*
            CALL TRITCO (NEL,LL4,ISPLH,IR,IQF,K,KK1,KK2,KK3,KK4,KK5,
     >      VOL0,MAT,MATN,SGD(1,1),XSGD(1,1),SIDE,ZZ,QFR(NUM2+1),IPR,A1)
*
            INX1=IPX(KEL)
            KEY0=MUX(INX1)-INX1
            IF(KK1.GT.0) THEN
               INX2=IPX(KK1)
               IF(INX2.LT.INX1) THEN
                  KEY=KEY0+INX2
                  A11X(KEY)=A11X(KEY)-REAL(A1(1))/2.
               ENDIF
            ENDIF
            IF(KK2.GT.0) THEN
               INX2=IPX(KK2)
               IF(INX2.LT.INX1) THEN
                  KEY=KEY0+INX2
                  A11X(KEY)=A11X(KEY)-REAL(A1(2))/2.
               ENDIF
            ENDIF
            KEY=KEY0+INX1
            VAR1 = A1(1)+A1(2)+A1(3)+A1(4)+A1(5)
            A11X(KEY)=A11X(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
  110    CONTINUE
  120    NUM1=NUM1+IVAL
         NUM2=NUM2+8
  130 CONTINUE
      DEALLOCATE(IWRK)
      RETURN
      END
*
      SUBROUTINE TRIMTY (ISPLH,IR,NEL,LL4,VOL,MAT,MATN,SGD,XSGD,SIDE,
     1              ZZ,KN,QFR,MUY,IPY,IPR,A11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ISPLH,IR,NEL,LL4,MAT(NEL),MATN(LL4),
     1 KN((18*(ISPLH-1)**2+3)*NEL),MUY(LL4),IPY(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11Y(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(5),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWRK
*----
*  ASSEMBLY OF MATRIX A11Y
*----
      NUM1=0
      NUM2=0
      NTPH = 6*(ISPLH-1)**2
      NTPL = 1+2*(ISPLH-1)
      NVT1 = NTPL + 2 * (ISPLH-2) + NTPH / 2
      NVT2 = NTPH - NTPL - (ISPLH-4) * (NTPL+2)
      NVT3 = NTPH - (ISPLH-4) * NTPL
      IVAL = 3*NTPH+8
      IF(ISPLH.EQ.3) NVT2 = NTPH
      IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
      IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
      ICR = ISAU*(1+2*(ISPLH-2))
      ALLOCATE(IWRK(NEL))
      MEL = 0
      DO 205 M=1,NEL
         IF(MAT(M).LE.0) GO TO 205
         MEL = MEL + 1
         IWRK(MEL) = M
205   CONTINUE
      DO 230 K=1,NEL
         L = MAT(K)
         IF(L.EQ.0) GO TO 230
         VOL0 = VOL(K)/NTPH
         IF(VOL0.EQ.0.0) GO TO 220
         KK4=KN(NUM1+3*NTPH+7)
         KK5=KN(NUM1+3*NTPH+8)
         IF(KK4.GT.0) KK4 = IWRK(KK4)
         IF(KK5.GT.0) KK5 = IWRK(KK5)
         DO 210 I = 1,NTPH
*
            CALL TRINEI (3,3,1,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,
     >                   NTPH,NTPL,NVT1,NVT2,NVT3,IVAL,KN)
*
            CALL TRITCO (NEL,LL4,ISPLH,IR,IQF,K,KK1,KK2,KK3,KK4,KK5,
     >      VOL0,MAT,MATN,SGD(1,1),XSGD(1,1),SIDE,ZZ,QFR(NUM2+1),IPR,A1)
*
            INY1=IPY(KEL)
            KEY0=MUY(INY1)-INY1
            IF(KK1.GT.0) THEN
               INY2=IPY(KK1)
               IF(INY2.LT.INY1) THEN
                  KEY=KEY0+INY2
                  A11Y(KEY)=A11Y(KEY)-REAL(A1(1))/2.
               ENDIF
            ENDIF
            IF(KK2.GT.0) THEN
               INY2=IPY(KK2)
               IF(INY2.LT.INY1) THEN
                  KEY=KEY0+INY2
                  A11Y(KEY)=A11Y(KEY)-REAL(A1(2))/2.
               ENDIF
            ENDIF
            KEY=KEY0+INY1
            VAR1 = A1(1)+A1(2)+A1(3)+A1(4)+A1(5)
            A11Y(KEY)=A11Y(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
  210    CONTINUE
  220    NUM1=NUM1+IVAL
         NUM2=NUM2+8
  230 CONTINUE
      DEALLOCATE(IWRK)
      RETURN
      END
*
      SUBROUTINE TRIMTZ (ISPLH,IR,NEL,LL4,VOL,MAT,MATN,SGD,XSGD,SIDE,
     1              ZZ,KN,QFR,MUZ,IPZ,IPR,A11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ISPLH,IR,NEL,LL4,MAT(NEL),MATN(LL4),
     1 KN((18*(ISPLH-1)**2+3)*NEL),MUZ(LL4),IPZ(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),SIDE,ZZ(NEL),QFR(8*NEL),
     1 A11Z(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION A1(5),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWRK
*----
*  ASSEMBLY OF MATRIX A11Z
*----
      NUM1=0
      NUM2=0
      NTPH = 6*(ISPLH-1)**2
      NTPL = 1+2*(ISPLH-1)
      NVT1 = NTPL + 2 * (ISPLH-2) + NTPH / 2
      NVT2 = NTPH - NTPL - (ISPLH-4) * (NTPL+2)
      NVT3 = NTPH - (ISPLH-4) * NTPL
      IVAL = 3*NTPH+8
      IF(ISPLH.EQ.3) NVT2 = NTPH
      IF(ISPLH.LE.3) ISAU =   2*(ISPLH-2)
      IF(ISPLH.GE.4) ISAU =   6*(ISPLH-3)
      ICR = ISAU*(1+2*(ISPLH-2))
      ALLOCATE(IWRK(NEL))
      MEL = 0
      DO 305 M=1,NEL
         IF(MAT(M).LE.0) GO TO 305
         MEL = MEL + 1
         IWRK(MEL) = M
305   CONTINUE
      DO 330 K=1,NEL
         L = MAT(K)
         IF(L.EQ.0) GO TO 330
         VOL0 = VOL(K)/NTPH
         IF(VOL0.EQ.0.0) GO TO 320
         DO 310 I = 1,NTPH
*
            CALL TRINEI (3,1,1,ISPLH,ICR,I,KK1,KK2,KK3,KEL,IQF,NUM1,
     >                   NTPH,NTPL,NVT1,NVT2,NVT3,IVAL,KN)
            KK4 = KN(NUM1+NTPH+I)
            KK5 = KN(NUM1+2*NTPH+I)
            LK4 = KK4
            LK5 = KK5
            IF(LK4.GT.0) LK4 = IWRK(KN(NUM1+3*NTPH+7))
            IF(LK5.GT.0) LK5 = IWRK(KN(NUM1+3*NTPH+8))
*
            CALL TRITCO (NEL,LL4,ISPLH,IR,IQF,K,KK1,KK2,KK3,LK4,LK5,
     >      VOL0,MAT,MATN,SGD(1,1),XSGD(1,1),SIDE,ZZ,QFR(NUM2+1),IPR,A1)
*
            INZ1=IPZ(KEL)
            KEY0=MUZ(INZ1)-INZ1
            IF(KK4.GT.0) THEN
               INZ2=IPZ(KK4)
               IF(INZ2.LT.INZ1) THEN
                  KEY=KEY0+INZ2
                  A11Z(KEY)=A11Z(KEY)-REAL(A1(4))
               ENDIF
            ENDIF
            IF(KK5.GT.0) THEN
               INZ2=IPZ(KK5)
               IF(INZ2.LT.INZ1) THEN
                  KEY=KEY0+INZ2
                  A11Z(KEY)=A11Z(KEY)-REAL(A1(5))
               ENDIF
            ENDIF
            KEY=KEY0+INZ1
            VAR1 = A1(1)+A1(2)+A1(3)+A1(4)+A1(5)
            A11Z(KEY)=A11Z(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
  310    CONTINUE
  320    NUM1=NUM1+IVAL
         NUM2=NUM2+8
  330 CONTINUE
      DEALLOCATE(IWRK)
      RETURN
      END
