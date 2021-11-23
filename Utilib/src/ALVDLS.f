*DECK ALVDLS
      SUBROUTINE ALVDLS (LTSW,MU1,ASS,F,ISEG,LON,NBL,LBL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* solution of a symmetric linear system where the coefficient matrix
* have been previously factorized as LDL(T). Supervectorial version.
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
* MU1     position of each diagonal element in vector ass.
*         DIMENSION MU1(L4) where L4=SUM(LBL(I))
* ASS     LDL(T) factors of the coefficient matrix in compressed
*         diagonal storage mode. DIMENSION ASS(ISEG,MU1(L4))
* F       right-hand side of the linear system. DIMENSION F(ISEG,L4)
* ISEG    number of elements in a vector register.
* LON     number of groups of linear systems.
* NBL     number of linear systems in each group. DIMENSION NBL(LON)
* LBL     number of unknowns in each group. DIMENSION LBL(LON)
*
*Parameters: output
* F       solution of the linear system. DIMENSION F(ISEG,L4)
*
*-----------------------------------------------------------------------
*
      INTEGER ISEG,LON,MU1(*),NBL(LON),LBL(LON)
      REAL ASS(ISEG,*),F(ISEG,*)
      REAL, DIMENSION(:), ALLOCATABLE :: T
*
      ALLOCATE(T(ISEG))
      LBL0=0
      IMAX=0
      DO 200 ILON=1,LON
      L4=LBL(ILON)
      NBANC=NBL(ILON)
      IF(L4.EQ.1) THEN
         IMAX=IMAX+1
CDIR$ SHORTLOOP
         DO 10 IB=1,NBANC
         F(IB,LBL0+1)=F(IB,LBL0+1)*ASS(IB,IMAX)
   10    CONTINUE
      ELSE IF(LTSW.GT.2) THEN
         IMAX=MU1(LBL0+L4)
         K1=MU1(LBL0+1)+1
         DO 55 I=LBL0+2,LBL0+L4
         K2=MU1(I)
         KJ=I-K2+K1
CDIR$ SHORTLOOP
         DO 20 IB=1,NBANC
         T(IB)=-F(IB,I)
   20    CONTINUE
         DO 40 K=K1,K2-1
CDIR$ SHORTLOOP
         DO 30 IB=1,NBANC
         T(IB)=T(IB)+F(IB,KJ)*ASS(IB,K)
   30    CONTINUE
         KJ=KJ+1
   40    CONTINUE
         K1=K2+1
CDIR$ SHORTLOOP
         DO 50 IB=1,NBANC
         F(IB,I)=-T(IB)
   50    CONTINUE
   55    CONTINUE
*
         DO 65 I=LBL0+1,LBL0+L4
         K1=MU1(I)
CDIR$ SHORTLOOP
         DO 60 IB=1,NBANC
         F(IB,I)=F(IB,I)*ASS(IB,K1)
   60    CONTINUE
   65    CONTINUE
*
         K2=IMAX
         DO 100 I=LBL0+L4,LBL0+2,-1
CDIR$ SHORTLOOP
         DO 70 IB=1,NBANC
         T(IB)=-F(IB,I)
   70    CONTINUE
         K1=MU1(I-1)+1
         KJ=I-K2+K1
         DO 90 K=K1,K2-1
CDIR$ SHORTLOOP
         DO 80 IB=1,NBANC
         F(IB,KJ)=F(IB,KJ)+ASS(IB,K)*T(IB)
   80    CONTINUE
         KJ=KJ+1
   90    CONTINUE
         K2=K1-1
  100    CONTINUE
      ELSE IF(LTSW.EQ.2) THEN
         K1=IMAX+2
         DO 130 I=LBL0+2,LBL0+L4
         KJ=I-1
CDIR$ SHORTLOOP
         DO 110 IB=1,NBANC
         T(IB)=-F(IB,I)+F(IB,KJ)*ASS(IB,K1)
  110    CONTINUE
CDIR$ SHORTLOOP
         DO 120 IB=1,NBANC
         F(IB,I)=-T(IB)
  120    CONTINUE
         K1=K1+2
  130    CONTINUE
*
         DO 145 I=LBL0+1,LBL0+L4
         K1=IMAX+2*(I-LBL0)-1
CDIR$ SHORTLOOP
         DO 140 IB=1,NBANC
         F(IB,I)=F(IB,I)*ASS(IB,K1)
  140    CONTINUE
  145    CONTINUE
*
         K1=IMAX+2*L4-2
         DO 170 I=LBL0+L4,LBL0+2,-1
         KJ=I-1
CDIR$ SHORTLOOP
         DO 150 IB=1,NBANC
         T(IB)=-F(IB,I)
  150    CONTINUE
CDIR$ SHORTLOOP
         DO 160 IB=1,NBANC
         F(IB,KJ)=F(IB,KJ)+ASS(IB,K1)*T(IB)
  160    CONTINUE
         K1=K1-2
  170    CONTINUE
         IMAX=IMAX+2*L4-1
      ENDIF
      LBL0=LBL0+L4
  200 CONTINUE
      DEALLOCATE(T)
      RETURN
      END
