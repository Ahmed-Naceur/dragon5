*DECK KINT02
      SUBROUTINE KINT02(MAXKN,SGD,IELEM,ICHX,IDIM,NREG,LL4,NBMIX,MAT,
     1 KN,VOL,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mixed-dual finite element
* diffusion approximation (Cartesian geometry). Special version for
* Trivac.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXKN   dimension of array KN.
* SGD     mixture-ordered cross sections.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICHX    type of discretization method:
*         =2: dual finite element approximations;
*         =3: nodal collocation method with full tensorial products;
*         =4: nodal collocation method with serendipity approximation.
* IDIM    number of dimensions.
* NREG    number of elements in Trivac.
* LL4     number of unknowns per group in Trivac.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* VOL     volume of regions.
* F2      vector to multiply.
*
*Parameters: output
* F3      result of the multiplication.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXKN,IELEM,ICHX,IDIM,NREG,LL4,NBMIX,MAT(NREG),KN(MAXKN)
      REAL SGD(NBMIX),VOL(NREG),F2(LL4),F3(LL4)
*----
*  LOCAL VARIABLES
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*
      IORD(J,K,L,LL,IEL,IW)=(IEL*L+K)*LL*IEL+(1+IEL*(IW-1))+J
*
      IORL(J,K,L,LL,IEL,IW)=
     1 1+LL*(L*(IEL*(IEL+1))/2-(L*(L-1)*(3*IEL-L+2))/6
     2 +K*(IEL-L)-(K*(K-1))/2)+(IEL-K-L)*(IW-1)+J
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IGAR(NREG))
*
      IF(ICHX.EQ.2) THEN
*        DUAL FINITE ELEMENT METHOD.
         NUM1=0
         DO 30 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 30
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 20
         DO 12 K3=0,IELEM-1
         DO 11 K2=0,IELEM-1
         DO 10 K1=0,IELEM-1
         IND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
         F3(IND1)=F3(IND1)+VOL0*SGD(L)*F2(IND1)
   10    CONTINUE
   11    CONTINUE
   12    CONTINUE
   20    NUM1=NUM1+1+6*IELEM**2
   30    CONTINUE
      ELSE IF(ICHX.EQ.3) THEN
*        NODAL COLLOCATION METHOD WITH FULL TENSORIAL PRODUCTS.
         LNUN=0
         DO 40 K=1,NREG
         IF(MAT(K).EQ.0) GO TO 40
         LNUN=LNUN+1
         IGAR(K)=LNUN
   40    CONTINUE
*
         DO 70 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 70
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 70
         DO 65 I3=0,IELEM-1
         DO 60 I2=0,IELEM-1
         DO 50 I1=0,IELEM-1
         INX1=IORD(I1,I2,I3,LNUN,IELEM,IGAR(K))
         F3(INX1)=F3(INX1)+VOL0*SGD(L)*F2(INX1)
   50    CONTINUE
         IF((IDIM.EQ.1).AND.(I2.EQ.0)) GO TO 70
         IF((IDIM.EQ.2).AND.(I2.EQ.IELEM-1)) GO TO 70
   60    CONTINUE
   65    CONTINUE
   70    CONTINUE
      ELSE IF(ICHX.EQ.4) THEN
*        NODAL COLLOCATION METHOD WITH SERENDIPITY APPROXIMATION.
         LNUN=0
         DO 80 K=1,NREG
         IF(MAT(K).EQ.0) GO TO 80
         LNUN=LNUN+1
         IGAR(K)=LNUN
   80    CONTINUE
*
         DO 110 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 110
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 110
         DO 105 I3=0,IELEM-1
         DO 100 I2=0,IELEM-1-I3
         DO 90 I1=0,IELEM-1-I2-I3
         INX1=IORL(I1,I2,I3,LNUN,IELEM,IGAR(K))
         F3(INX1)=F3(INX1)+VOL0*SGD(L)*F2(INX1)
   90    CONTINUE
         IF((IDIM.EQ.1).AND.(I2.EQ.0)) GO TO 110
         IF((IDIM.EQ.2).AND.(I2.EQ.IELEM-1)) GO TO 110
  100    CONTINUE
  105    CONTINUE
  110    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IGAR)
      RETURN
      END
