*DECK ALVDLM
      SUBROUTINE ALVDLM (LTSW,ASS,VEC,Z,MU1,ITY,ISEG,LON,NBL,LBL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* multiplication of a symmetric matrix in compressed diagonal storage
* mode by a vector. Supervectorial version.
*
*Copyright:
* Copyright (C) 1992 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* LTSW    maximum bandwidth. =2 for tridiagonal systems.
* ASS     LDL(T) factors of the coefficient matrix in compressed
*         diagonal storage mode. DIMENSION ASS(ISEG,MU1(L4))
* VEC     vector to multiply.
* Z       vector that will be added to the result if ITY=2.
*         DIMENSION VEC(ISEG,L4),Z(ISEG,L4) WHERE L4=SUM(LBL(I))
* MU1     position of each diagonal element in vector ASS.
*         DIMENSION MU1(L4)
* ITY     type of multiplication (ITY=1: Z=ASS*VEC;
*         ITY=2: Z=Z+(ASS-DIAG(ASS))*VEC).
* ISEG    number of elements in a vector register.
* LON     number of groups of linear systems.
* NBL     number of linear systems in each group. DIMENSION NBL(LON)
* LBL     number of unknowns in each group. DIMENSION LBL(LON)
*
*Parameters: output
* Z       solution of the multiplication. DIMENSION Z(ISEG,L4)
*
*-----------------------------------------------------------------------
*
      DIMENSION ASS(ISEG,*),VEC(ISEG,*),Z(ISEG,*),MU1(*),NBL(*),LBL(*)
      LBL0=0
      IMAX=0
      DO 300 ILON=1,LON
      L4=LBL(ILON)
      NBANC=NBL(ILON)
      IF((LTSW.GT.2).AND.(ITY.EQ.1)) THEN
*        CALCULATION OF Z=ASS*VEC.
CDIR$ SHORTLOOP
         DO 30 IB=1,NBANC
         Z(IB,LBL0+1)=ASS(IB,MU1(LBL0+1))*VEC(IB,LBL0+1)
   30    CONTINUE
         I1=MU1(LBL0+1)+1
         DO 90 K=LBL0+2,LBL0+L4
         I2=MU1(K)
         KEY1=I2-K
CDIR$ SHORTLOOP
         DO 40 IB=1,NBANC
         Z(IB,K)=0.0
   40    CONTINUE
         DO 60 L=K+I1-I2,K-1
         KEYL=KEY1+L
CDIR$ SHORTLOOP
         DO 50 IB=1,NBANC
         Z(IB,K)=Z(IB,K)+ASS(IB,KEYL)*VEC(IB,L)
         Z(IB,L)=Z(IB,L)+ASS(IB,KEYL)*VEC(IB,K)
   50    CONTINUE
   60    CONTINUE
CDIR$ SHORTLOOP
         DO 80 IB=1,NBANC
         Z(IB,K)=Z(IB,K)+ASS(IB,KEY1+K)*VEC(IB,K)
   80    CONTINUE
         I1=I2+1
   90    CONTINUE
      ELSE IF((LTSW.GT.2).AND.(ITY.EQ.2)) THEN
*        CALCULATION OF Z=Z+(ASS-DIAG(ASS))*VEC.
         I1=MU1(LBL0+1)+1
         DO 150 K=LBL0+2,LBL0+L4
         I2=MU1(K)
         KEY1=I2-K
         IF(I1.EQ.I2) GO TO 150
         DO 130 L=K+I1-I2,K-1
         KEYL=KEY1+L
CDIR$ SHORTLOOP
         DO 120 IB=1,NBANC
         Z(IB,K)=Z(IB,K)+ASS(IB,KEYL)*VEC(IB,L)
         Z(IB,L)=Z(IB,L)+ASS(IB,KEYL)*VEC(IB,K)
  120    CONTINUE
  130    CONTINUE
         I1=I2+1
  150    CONTINUE
      ELSE IF((LTSW.EQ.2).AND.(ITY.EQ.1)) THEN
*        CALCULATION OF Z=ASS*VEC FOR A 3-DIAGONAL SYSTEM.
CDIR$ SHORTLOOP
         DO 180 IB=1,NBANC
         Z(IB,LBL0+1)=ASS(IB,IMAX+1)*VEC(IB,LBL0+1)
  180    CONTINUE
         I1=2
         DO 230 K=LBL0+2,LBL0+L4
         KEYL=IMAX+I1
CDIR$ SHORTLOOP
         DO 210 IB=1,NBANC
         Z(IB,K)=ASS(IB,KEYL)*VEC(IB,K-1)+ASS(IB,KEYL+1)*VEC(IB,K)
         Z(IB,K-1)=Z(IB,K-1)+ASS(IB,KEYL)*VEC(IB,K)
  210    CONTINUE
         I1=I1+2
  230    CONTINUE
         IMAX=IMAX+I1-1
      ELSE IF((LTSW.EQ.2).AND.(ITY.EQ.2)) THEN
*        CALCULATION OF Z=Z+(ASS-DIAG(ASS))*VEC FOR A 3-DIAGONAL SYSTEM.
         I1=2
         DO 280 K=LBL0+2,LBL0+L4
         KEYL=IMAX+I1
CDIR$ SHORTLOOP
         DO 260 IB=1,NBANC
         Z(IB,K)=Z(IB,K)+ASS(IB,KEYL)*VEC(IB,K-1)
         Z(IB,K-1)=Z(IB,K-1)+ASS(IB,KEYL)*VEC(IB,K)
  260    CONTINUE
         I1=I1+2
  280    CONTINUE
         IMAX=IMAX+I1-1
      ENDIF
      LBL0=LBL0+L4
  300 CONTINUE
      RETURN
      END
