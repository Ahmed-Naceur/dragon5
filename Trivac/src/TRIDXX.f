*DECK TRIDXX
      SUBROUTINE TRIDXX(NBMIX,CYLIND,IELEM,ICOL,NEL,LL4F,LL4X,MAT,VOL,
     1 XX,YY,ZZ,DD,KN,QFR,SGD,XSGD,MUX,IPBBX,LC,R,V,BBX,TTF,AX,C11X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a Thomas-Raviart (dual) finite element
* method in Cartesian 3-D diffusion approximation.
* Note: system matrices should be initialized by the calling program.
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
* NBMIX   number of mixtures.
* CYLIND  cylindrical geometry flag (set with CYLIND=.true.).
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NEL     total number of finite elements.
* LL4F    number of flux components.
* LL4X    number of X-directed currents.
* LL4Y    number of Y-directed currents.
* LL4Z    number of Z-directed currents.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      used with cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* SGD     nuclear properties by material mixture:
*         SGD(L,1)= X-oriented diffusion coefficients;
*         SGD(L,2)= Y-oriented diffusion coefficients;
*         SGD(L,3)= Z-oriented diffusion coefficients;
*         SGD(L,4)= removal macroscopic cross section.
* XSGD    one over nuclear properties.
* MUX     X-directed compressed storage mode indices.
* MUY     Y-directed compressed storage mode indices.
* MUZ     Z-directed compressed storage mode indices.
* IPBBX   X-directed perdue storage indices.
* IPBBY   Y-directed perdue storage indices.
* IPBBZ   Z-directed perdue storage indices.
* LC      order of the unit matrices.
* R       unit matrix.
* V       unit matrix.
* BBX     X-directed flux-current matrices.
* BBY     Y-directed flux-current matrices.
* BBZ     Z-directed flux-current matrices.
*
*Parameters: output
* TTF     flux-flux matrices.
* AX      X-directed main current-current matrices. Dimensionned to
*         MUX(LL4X).
* AY      Y-directed main current-current matrices. Dimensionned to
*         MUY(LL4Y).
* AZ      Z-directed main current-current matrices. Dimensionned to
*         MUZ(LL4Z).
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
      INTEGER NBMIX,IELEM,ICOL,NEL,LL4F,LL4X,MAT(NEL),
     1 KN(NEL*(1+6*IELEM**2)),MUX(LL4X),IPBBX(2*IELEM,LL4X),LC
      REAL VOL(NEL),XX(NEL),YY(NEL),ZZ(NEL),DD(NEL),QFR(6*NEL),
     1 SGD(NBMIX,4),XSGD(NBMIX,4),R(LC,LC),V(LC,LC-1),TTF(LL4F),
     2 BBX(2*IELEM,LL4X),AX(*),C11X(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION FFF
      REAL QQ(5,5)
*----
*  X-ORIENTED COUPLINGS
*----
      IF((CYLIND).AND.((IELEM.GT.1).OR.(ICOL.NE.2)))
     1 CALL XABORT('TRIDXX: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
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
      NUM1=0
      NUM2=0
      DO 60 IE=1,NEL
      L=MAT(IE)
      IF(L.EQ.0) GO TO 60
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 50
      DX=XX(IE)
      DY=YY(IE)
      DZ=ZZ(IE)
      IF(CYLIND) THEN
         DIN=1.0-0.5*DX/DD(IE)
         DOT=1.0+0.5*DX/DD(IE)
      ELSE
         DIN=1.0
         DOT=1.0
      ENDIF
*
      DO 45 K3=0,IELEM-1
      DO 40 K2=0,IELEM-1
      KN1=KN(NUM1+2+K3*IELEM+K2)
      KN2=KN(NUM1+2+IELEM**2+K3*IELEM+K2)
      INX1=ABS(KN1)-LL4F
      INX2=ABS(KN2)-LL4F
      DO 30 K1=0,IELEM-1
      JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
      TTF(JND1)=TTF(JND1)+VOL0*SGD(L,4)+VOL0*QQ(K1+1,K1+1)*SGD(L,1)/
     1 (DX*DX)
      TTF(JND1)=TTF(JND1)+VOL0*QQ(K2+1,K2+1)*SGD(L,2)/(DY*DY)
      TTF(JND1)=TTF(JND1)+VOL0*QQ(K3+1,K3+1)*SGD(L,3)/(DZ*DZ)
   30 CONTINUE
      IF(KN1.NE.0) THEN
         KEY=MUX(INX1)
         AX(KEY)=AX(KEY)-DIN*(VOL0*R(1,1)*XSGD(L,1)+QFR(NUM2+1))
      ENDIF
      IF(KN2.NE.0) THEN
         KEY=MUX(INX2)
         AX(KEY)=AX(KEY)-DOT*(VOL0*R(IELEM+1,IELEM+1)*XSGD(L,1)
     1   +QFR(NUM2+2))
      ENDIF
      IF((ICOL.NE.2).AND.(KN1.NE.0).AND.(KN2.NE.0)) THEN
         IF(INX2.GT.INX1) KEY=MUX(INX2)-INX2+INX1
         IF(INX2.LE.INX1) KEY=MUX(INX1)-INX1+INX2
         SG=REAL(SIGN(1,KN1)*SIGN(1,KN2))
         IF(INX1.EQ.INX2) SG=2.0*SG
         AX(KEY)=AX(KEY)-SG*VOL0*R(IELEM+1,1)*XSGD(L,1)
      ENDIF
   40 CONTINUE
   45 CONTINUE
   50 NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
   60 CONTINUE
*
      DO 121 I0=1,MUX(LL4X)
      C11X(I0)=-AX(I0)
  121 CONTINUE
      MUIM1=0
      DO 716 I=1,LL4X
      MUI=MUX(I)
      DO 715 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 714 I0=1,2*IELEM
      II=IPBBX(I0,I)
      IF(II.EQ.0) GO TO 715
      DO 713 J0=1,2*IELEM
      JJ=IPBBX(J0,J)
      IF(II.EQ.JJ) C11X(KEY)=C11X(KEY)+BBX(I0,I)*BBX(J0,J)/TTF(II)
  713 CONTINUE
  714 CONTINUE
  715 CONTINUE
      MUIM1=MUI
  716 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIDXY(NBMIX,IELEM,ICOL,NEL,LL4F,LL4X,LL4Y,MAT,VOL,YY,
     1 KN,QFR,XSGD,MUY,IPBBY,LC,R,BBY,TTF,AY,C11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,IELEM,ICOL,NEL,LL4F,LL4X,LL4Y,MAT(NEL),
     1 KN(NEL*(1+6*IELEM**2)),MUY(LL4Y),IPBBY(2*IELEM,LL4Y),LC
      REAL VOL(NEL),YY(NEL),QFR(6*NEL),XSGD(NBMIX,4),R(LC,LC),TTF(LL4F),
     1 BBY(2*IELEM,LL4Y),AY(*),C11Y(*)
*----
*  Y-ORIENTED COUPLINGS
*----
      NUM1=0
      NUM2=0
      DO 240 IE=1,NEL
      L=MAT(IE)
      IF(L.EQ.0) GO TO 240
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 230
      DY=YY(IE)
*
      DO 195 K3=0,IELEM-1
      DO 190 K1=0,IELEM-1
      KN1=KN(NUM1+2+2*IELEM**2+K3*IELEM+K1)
      KN2=KN(NUM1+2+3*IELEM**2+K3*IELEM+K1)
      INY1=ABS(KN1)-LL4F-LL4X
      INY2=ABS(KN2)-LL4F-LL4X
      IF(KN1.NE.0) THEN
         KEY=MUY(INY1)
         AY(KEY)=AY(KEY)-VOL0*R(1,1)*XSGD(L,2)-QFR(NUM2+3)
      ENDIF
      IF(KN2.NE.0) THEN
         KEY=MUY(INY2)
         AY(KEY)=AY(KEY)-VOL0*R(IELEM+1,IELEM+1)*XSGD(L,2)
     1   -QFR(NUM2+4)
      ENDIF
      IF((ICOL.NE.2).AND.(KN1.NE.0).AND.(KN2.NE.0)) THEN
         IF(INY2.GT.INY1) KEY=MUY(INY2)-INY2+INY1
         IF(INY2.LE.INY1) KEY=MUY(INY1)-INY1+INY2
         SG=REAL(SIGN(1,KN1)*SIGN(1,KN2))
         IF(INY1.EQ.INY2) SG=2.0*SG
         AY(KEY)=AY(KEY)-SG*VOL0*R(IELEM+1,1)*XSGD(L,2)
      ENDIF
  190 CONTINUE
  195 CONTINUE
  230 NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  240 CONTINUE
*
      DO 212 I0=1,MUY(LL4Y)
      C11Y(I0)=-AY(I0)
  212 CONTINUE
      MUIM1=0
      DO 216 I=1,LL4Y
      MUI=MUY(I)
      DO 215 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 214 I0=1,2*IELEM
      II=IPBBY(I0,I)
      IF(II.EQ.0) GO TO 215
      DO 213 J0=1,2*IELEM
      JJ=IPBBY(J0,J)
      IF(II.EQ.JJ) C11Y(KEY)=C11Y(KEY)+BBY(I0,I)*BBY(J0,J)/TTF(II)
  213 CONTINUE
  214 CONTINUE
  215 CONTINUE
      MUIM1=MUI
  216 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIDXZ(NBMIX,IELEM,ICOL,NEL,LL4F,LL4X,LL4Y,LL4Z,MAT,
     1 VOL,ZZ,KN,QFR,XSGD,MUZ,IPBBZ,LC,R,BBZ,TTF,AZ,C11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,IELEM,ICOL,NEL,LL4F,LL4X,LL4Y,LL4Z,MAT(NEL),
     1 KN(NEL*(1+6*IELEM**2)),MUZ(LL4Z),IPBBZ(2*IELEM,LL4Z),LC
      REAL VOL(NEL),ZZ(NEL),QFR(6*NEL),XSGD(NBMIX,4),R(LC,LC),TTF(LL4F),
     1 BBZ(2*IELEM,LL4Z),AZ(*),C11Z(*)
*----
*  Z-ORIENTED COUPLINGS
*----
      NUM1=0
      NUM2=0
      DO 340 IE=1,NEL
      L=MAT(IE)
      IF(L.EQ.0) GO TO 340
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 330
      DZ=ZZ(IE)
*
      DO 295 K2=0,IELEM-1
      DO 290 K1=0,IELEM-1
      KN1=KN(NUM1+2+4*IELEM**2+K2*IELEM+K1)
      KN2=KN(NUM1+2+5*IELEM**2+K2*IELEM+K1)
      INZ1=ABS(KN1)-LL4F-LL4X-LL4Y
      INZ2=ABS(KN2)-LL4F-LL4X-LL4Y
      IF(KN1.NE.0) THEN
         KEY=MUZ(INZ1)
         AZ(KEY)=AZ(KEY)-VOL0*R(1,1)*XSGD(L,3)-QFR(NUM2+5)
      ENDIF
      IF(KN2.NE.0) THEN
         KEY=MUZ(INZ2)
         AZ(KEY)=AZ(KEY)-VOL0*R(IELEM+1,IELEM+1)*XSGD(L,3)
     1   -QFR(NUM2+6)
      ENDIF
      IF((ICOL.NE.2).AND.(KN1.NE.0).AND.(KN2.NE.0)) THEN
         IF(INZ2.GT.INZ1) KEY=MUZ(INZ2)-INZ2+INZ1
         IF(INZ2.LE.INZ1) KEY=MUZ(INZ1)-INZ1+INZ2
         SG=REAL(SIGN(1,KN1)*SIGN(1,KN2))
         IF(INZ1.EQ.INZ2) SG=2.0*SG
         AZ(KEY)=AZ(KEY)-SG*VOL0*R(IELEM+1,1)*XSGD(L,3)
      ENDIF
  290 CONTINUE
  295 CONTINUE
  330 NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  340 CONTINUE
*
      DO 312 I0=1,MUZ(LL4Z)
      C11Z(I0)=-AZ(I0)
  312 CONTINUE
      MUIM1=0
      DO 316 I=1,LL4Z
      MUI=MUZ(I)
      DO 315 J=I-(MUI-MUIM1)+1,I
      KEY=MUI-I+J
      DO 314 I0=1,2*IELEM
      II=IPBBZ(I0,I)
      IF(II.EQ.0) GO TO 315
      DO 313 J0=1,2*IELEM
      JJ=IPBBZ(J0,J)
      IF(II.EQ.JJ) C11Z(KEY)=C11Z(KEY)+BBZ(I0,I)*BBZ(J0,J)/TTF(II)
  313 CONTINUE
  314 CONTINUE
  315 CONTINUE
      MUIM1=MUI
  316 CONTINUE
      RETURN
      END
