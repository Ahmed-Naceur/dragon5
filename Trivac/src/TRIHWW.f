*DECK TRIHWW
      SUBROUTINE TRIHWW(NBMIX,NBLOS,IELEM,LL4F,LL4W,MAT,SIDE,ZZ,FRZ,
     1 QFR,IPERT,KN,SGD,XSGD,MUW,IPBBW,LC,R,V,BBW,TTF,AW,C11W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a Thomas-Raviart-Schneider (dual)
* finite element method in hexagonal 3-D diffusion approximation.
* Note: system matrices should be initialized by the calling program.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NBMIX   number of mixtures.
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* ISPLH   mesh-splitting index. Each hexagon is splitted into 3*ISPLH**2
*         lozenges.
* LL4F    number of flux components.
* LL4W    number of W-directed currents.
* LL4X    number of X-directed currents.
* LL4Y    number of Y-directed currents.
* LL4Z    number of Z-directed currents.
* MAT     mixture index assigned to each element.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* FRZ     volume fractions for the axial SYME boundary condition.
* QFR     element-ordered boundary conditions.
* IPERT   mixture permutation index.
* KN      ADI permutation indices for the volumes and currents.
* SGD     nuclear properties by material mixture:
*         SGD(L,1)= X-oriented diffusion coefficients;
*         SGD(L,2)= Y-oriented diffusion coefficients;
*         SGD(L,3)= Z-oriented diffusion coefficients;
*         SGD(L,4)= removal macroscopic cross section.
* XSGD    one over nuclear properties.
* MUW     W-directed compressed storage mode indices.
* MUX     X-directed compressed storage mode indices.
* MUY     Y-directed compressed storage mode indices.
* MUZ     Z-directed compressed storage mode indices.
* IPBBW   W-directed perdue storage indices.
* IPBBX   X-directed perdue storage indices.
* IPBBY   Y-directed perdue storage indices.
* IPBBZ   Z-directed perdue storage indices.
* LC      order of the unit matrices.
* R       unit matrix.
* V       unit matrix.
* BBW     W-directed flux-current matrices.
* BBX     X-directed flux-current matrices.
* BBY     Y-directed flux-current matrices.
* BBZ     Z-directed flux-current matrices.
*
*Parameters: output
* TTF     flux-flux matrices.
* AW      W-directed main current-current matrices. Dimensionned to
*         MUW(LL4W).
* AX      X-directed main current-current matrices. Dimensionned to
*         MUX(LL4X).
* AY      Y-directed main current-current matrices. Dimensionned to
*         MUY(LL4Y).
* AZ      Z-directed main current-current matrices. Dimensionned to
*         MUZ(LL4Z).
* C11W    W-directed main current-current matrices to be factorized.
*         Dimensionned to MUW(LL4W).
* C11X    X-directed main current-current matrices to be factorized.
*         Dimensionned to MUX(LL4X).
* C11Y    Y-directed main current-current matrices to be factorized.
*         Dimensionned to MUY(LL4Y).
* C11Z    Z-directed main current-current matrices to be factorized.
*         Dimensionned to MUZ(LL4Z).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,LL4F,LL4W,MAT(3,NBLOS),IPERT(NBLOS),
     1 KN(NBLOS,3+6*(IELEM+2)*IELEM**2),MUW(LL4W),IPBBW(2*IELEM,LL4W),LC
      REAL SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),SGD(NBMIX,4),
     1 XSGD(NBMIX,4),R(LC,LC),V(LC,LC-1),TTF(LL4F),BBW(2*IELEM,LL4W),
     2 AW(*),C11W(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION FFF,TTTT
      REAL QQ(5,5)
*----
*  W-ORIENTED COUPLINGS
*----
      DO 25 I0=1,IELEM
      DO 20 J0=1,IELEM
      FFF=0.0D0
      DO 10 K0=2,IELEM
      FFF=FFF+V(K0,I0)*V(K0,J0)/R(K0,K0)
   10 CONTINUE
      IF(ABS(FFF).LE.1.0E-6) FFF=0.0D0
      QQ(I0,J0)=REAL(FFF)
   20 CONTINUE
   25 CONTINUE
*
      NELEH=(IELEM+1)*IELEM**2
      IIMAW=MUW(LL4W)
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 50 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 50
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 50
      DZ=ZZ(1,IPERT(KEL))
      VOL0=REAL(TTTT*DZ*FRZ(KEL))
      DINV=XSGD(IBM,1)
      SIG3=SGD(IBM,3)/(DZ*DZ)
      SIG4=SGD(IBM,4)
      DO 34 K5=0,1
      DO 33 K4=0,IELEM-1
      DO 32 K3=0,IELEM-1
      DO 31 K2=1,IELEM+1
      KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INW1=ABS(KNW1)
      DO 30 K1=1,IELEM+1
      KNW2=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      INW2=ABS(KNW2)
      IF((KNW2.NE.0).AND.(KNW1.NE.0)) THEN
         L=MUW(INW1)-INW1+INW2
         SG=REAL(SIGN(1,KNW1)*SIGN(1,KNW2))
         IF(K1.LE.K2) AW(L)=AW(L)-(4./3.)*SG*VOL0*DINV*R(K2,K1)
         IF(K1.EQ.K2) THEN
            IF((K1.EQ.1).AND.(K5.EQ.0)) AW(L)=AW(L)-QFR(NUM,1)
            IF((K1.EQ.IELEM+1).AND.(K5.EQ.1)) AW(L)=AW(L)-QFR(NUM,2)
         ENDIF
      ENDIF
   30 CONTINUE
   31 CONTINUE
   32 CONTINUE
   33 CONTINUE
   34 CONTINUE
      DO 42 K3=0,IELEM-1
      DO 41 K2=0,IELEM-1
      DO 40 K1=0,IELEM-1
      JND1=(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      JND2=(KN(NUM,1)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      JND3=(KN(NUM,2)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
      TTF(JND1)=TTF(JND1)+VOL0*SIG4+VOL0*QQ(K3+1,K3+1)*SIG3
      TTF(JND2)=TTF(JND2)+VOL0*SIG4+VOL0*QQ(K3+1,K3+1)*SIG3
      TTF(JND3)=TTF(JND3)+VOL0*SIG4+VOL0*QQ(K3+1,K3+1)*SIG3
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
   50 CONTINUE
*----
*  COMPUTE THE W-ORIENTED SYSTEM MATRIX AFTER FLUX ELIMINATION
*----
      DO 60 I0=1,IIMAW
      C11W(I0)=-AW(I0)
   60 CONTINUE
      MUIM1=0
      DO 90 I=1,LL4W
      MUI=MUW(I)
      DO 80 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 75 I0=1,2*IELEM
      II=IPBBW(I0,I)
      IF(II.EQ.0) GO TO 80
      DO 70 J0=1,2*IELEM
      JJ=IPBBW(J0,J)
      IF(II.EQ.JJ) C11W(KEY)=C11W(KEY)+BBW(I0,I)*BBW(J0,J)/TTF(II)
   70 CONTINUE
   75 CONTINUE
   80 CONTINUE
      MUIM1=MUI
   90 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIHWX(NBMIX,NBLOS,IELEM,LL4F,LL4W,LL4X,MAT,SIDE,ZZ,
     1 FRZ,QFR,IPERT,KN,XSGD,MUX,IPBBX,LC,R,BBX,TTF,AX,C11X)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,LL4F,LL4W,LL4X,MAT(3,NBLOS),
     1 MUX(LL4X),IPBBX(2*IELEM,LL4X),LC,IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),XSGD(NBMIX,4),
     1 R(LC,LC),TTF(LL4F),BBX(2*IELEM,LL4X),AX(*),C11X(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT
*----
*  X-ORIENTED COUPLINGS
*----
      NELEH=(IELEM+1)*IELEM**2
      IIMAX=MUX(LL4X)
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 120 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 120
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 120
      VOL0=REAL(TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL))
      DINV=XSGD(IBM,1)
      DO 114 K5=0,1
      DO 113 K4=0,IELEM-1
      DO 112 K3=0,IELEM-1
      DO 111 K2=1,IELEM+1
      KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INX1=ABS(KNX1)-LL4W
      DO 110 K1=1,IELEM+1
      KNX2=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      INX2=ABS(KNX2)-LL4W
      IF((KNX2.NE.0).AND.(KNX1.NE.0)) THEN
         L=MUX(INX1)-INX1+INX2
         SG=REAL(SIGN(1,KNX1)*SIGN(1,KNX2))
         IF(K1.LE.K2) AX(L)=AX(L)-(4./3.)*SG*VOL0*DINV*R(K2,K1)
         IF(K1.EQ.K2) THEN
            IF((K1.EQ.1).AND.(K5.EQ.0)) AX(L)=AX(L)-QFR(NUM,3)
            IF((K1.EQ.IELEM+1).AND.(K5.EQ.1)) AX(L)=AX(L)-QFR(NUM,4)
         ENDIF
      ENDIF
  110 CONTINUE
  111 CONTINUE
  112 CONTINUE
  113 CONTINUE
  114 CONTINUE
  120 CONTINUE
*----
*  COMPUTE THE X-ORIENTED SYSTEM MATRIX AFTER FLUX ELIMINATION
*----
      DO 130 I0=1,IIMAX
      C11X(I0)=-AX(I0)
  130 CONTINUE
      MUIM1=0
      DO 160 I=1,LL4X
      MUI=MUX(I)
      DO 150 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 145 I0=1,2*IELEM
      II=IPBBX(I0,I)
      IF(II.EQ.0) GO TO 150
      DO 140 J0=1,2*IELEM
      JJ=IPBBX(J0,J)
      IF(II.EQ.JJ) C11X(KEY)=C11X(KEY)+BBX(I0,I)*BBX(J0,J)/TTF(II)
  140 CONTINUE
  145 CONTINUE
  150 CONTINUE
      MUIM1=MUI
  160 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIHWY(NBMIX,NBLOS,IELEM,LL4F,LL4W,LL4X,LL4Y,MAT,
     1 SIDE,ZZ,FRZ,QFR,IPERT,KN,XSGD,MUY,IPBBY,LC,R,BBY,TTF,AY,C11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,LL4F,LL4W,LL4X,LL4Y,MAT(3,NBLOS),
     1 MUY(LL4Y),IPBBY(2*IELEM,LL4Y),LC,IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),XSGD(NBMIX,4),
     1 R(LC,LC),TTF(LL4F),BBY(2*IELEM,LL4Y),AY(*),C11Y(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT
*----
*  Y-ORIENTED COUPLINGS
*----
      NELEH=(IELEM+1)*IELEM**2
      IIMAY=MUY(LL4Y)
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 220 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 220
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 220
      VOL0=REAL(TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL))
      DINV=XSGD(IBM,1)
      DO 214 K5=0,1
      DO 213 K4=0,IELEM-1
      DO 212 K3=0,IELEM-1
      DO 211 K2=1,IELEM+1
      KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INY1=ABS(KNY1)-LL4W-LL4X
      DO 210 K1=1,IELEM+1
      KNY2=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      INY2=ABS(KNY2)-LL4W-LL4X
      IF((KNY2.NE.0).AND.(KNY1.NE.0)) THEN
         L=MUY(INY1)-INY1+INY2
         SG=REAL(SIGN(1,KNY1)*SIGN(1,KNY2))
         IF(K1.LE.K2) AY(L)=AY(L)-(4./3.)*SG*VOL0*DINV*R(K2,K1)
         IF(K1.EQ.K2) THEN
            IF((K1.EQ.1).AND.(K5.EQ.0)) AY(L)=AY(L)-QFR(NUM,5)
            IF((K1.EQ.IELEM+1).AND.(K5.EQ.1)) AY(L)=AY(L)-QFR(NUM,6)
         ENDIF
      ENDIF
  210 CONTINUE
  211 CONTINUE
  212 CONTINUE
  213 CONTINUE
  214 CONTINUE
  220 CONTINUE
*----
*  COMPUTE THE Y-ORIENTED SYSTEM MATRIX AFTER FLUX ELIMINATION
*----
      DO 230 I0=1,IIMAY
      C11Y(I0)=-AY(I0)
  230 CONTINUE
      MUIM1=0
      DO 260 I=1,LL4Y
      MUI=MUY(I)
      DO 250 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 245 I0=1,2*IELEM
      II=IPBBY(I0,I)
      IF(II.EQ.0) GO TO 250
      DO 240 J0=1,2*IELEM
      JJ=IPBBY(J0,J)
      IF(II.EQ.JJ) C11Y(KEY)=C11Y(KEY)+BBY(I0,I)*BBY(J0,J)/TTF(II)
  240 CONTINUE
  245 CONTINUE
  250 CONTINUE
      MUIM1=MUI
  260 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIHWZ(NBMIX,NBLOS,IELEM,ICOL,LL4F,LL4W,LL4X,LL4Y,
     1 LL4Z,MAT,SIDE,ZZ,FRZ,QFR,IPERT,KN,XSGD,MUZ,IPBBZ,LC,R,BBZ,TTF,
     2 AZ,C11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,ICOL,LL4F,LL4W,LL4X,LL4Y,LL4Z,
     1 MAT(3,NBLOS),MUZ(LL4Z),IPBBZ(2*IELEM,LL4Z),LC,IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),XSGD(NBMIX,4),
     1 R(LC,LC),TTF(LL4F),BBZ(2*IELEM,LL4Z),AZ(*),C11Z(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT
*----
*  Z-ORIENTED COUPLINGS
*----
      NELEH=(IELEM+1)*IELEM**2
      IIMAZ=MUZ(LL4Z)
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 340 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 340
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 340
      VOL0=REAL(TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL))
      DINV=XSGD(IBM,3)
      DO 292 K5=0,2 ! THREE LOZENGES PER HEXAGON
      DO 291 K2=0,IELEM-1
      DO 290 K1=0,IELEM-1
      KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
      KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
      INZ1=ABS(KNZ1)-LL4W-LL4X-LL4Y
      INZ2=ABS(KNZ2)-LL4W-LL4X-LL4Y
      IF(KNZ1.NE.0) THEN
         KEY=MUZ(INZ1)
         AZ(KEY)=AZ(KEY)-VOL0*R(1,1)*DINV-QFR(NUM,7)
      ENDIF
      IF(KNZ2.NE.0) THEN
         KEY=MUZ(INZ2)
         AZ(KEY)=AZ(KEY)-VOL0*R(IELEM+1,IELEM+1)*DINV-QFR(NUM,8)
      ENDIF
      IF((ICOL.NE.2).AND.(KNZ1.NE.0).AND.(KNZ2.NE.0)) THEN
         IF(INZ2.GT.INZ1) KEY=MUZ(INZ2)-INZ2+INZ1
         IF(INZ2.LE.INZ1) KEY=MUZ(INZ1)-INZ1+INZ2
         SG=REAL(SIGN(1,KNZ1)*SIGN(1,KNZ2))
         IF(INZ1.EQ.INZ2) SG=2.0*SG
         AZ(KEY)=AZ(KEY)-SG*VOL0*R(IELEM+1,1)*DINV
      ENDIF
  290 CONTINUE
  291 CONTINUE
  292 CONTINUE
  340 CONTINUE
*
      DO 350 I0=1,IIMAZ
      C11Z(I0)=-AZ(I0)
  350 CONTINUE
      MUIM1=0
      DO 380 I=1,LL4Z
      MUI=MUZ(I)
      DO 370 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 365 I0=1,2*IELEM
      II=IPBBZ(I0,I)
      IF(II.EQ.0) GO TO 370
      DO 360 J0=1,2*IELEM
      JJ=IPBBZ(J0,J)
      IF(II.EQ.JJ) C11Z(KEY)=C11Z(KEY)+BBZ(I0,I)*BBZ(J0,J)/TTF(II)
  360 CONTINUE
  365 CONTINUE
  370 CONTINUE
      MUIM1=MUI
  380 CONTINUE
      RETURN
      END
