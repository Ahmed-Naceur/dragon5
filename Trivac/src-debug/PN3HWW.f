*DECK PN3HWW
      SUBROUTINE PN3HWW(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,
     1 MAT,SIGT,SIGTI,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUW,IPBBW,LC,R,V,BBW,
     2 TTF,AW,C11W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a Thomas-Raviart-Schneider (dual)
* finite element method in hexagonal 3-D simplified PN approximation.
* Note: system matrices should be initialized by the calling program.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
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
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* LL4F    number of flux components.
* LL4W    number of W-directed currents.
* LL4X    number of X-directed currents.
* LL4Y    number of Y-directed currents.
* LL4Z    number of Z-directed currents.
* MAT     mixture index assigned to each lozenge.
* SIGT    total minus self-scattering macroscopic cross sections.
*         SIGT(:,NAN) generally contains the total cross section only.
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* FRZ     volume fractions for the axial SYME boundary condition.
* QFR     element-ordered boundary conditions.
* IPERT   mixture permutation index.
* KN      ADI permutation indices for the volumes and currents.
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
*         MUW(LL4W)*NLF/2.
* AX      X-directed main current-current matrices. Dimensionned to
*         MUX(LL4X)*NLF/2.
* AY      Y-directed main current-current matrices. Dimensionned to
*         MUY(LL4Y)*NLF/2.
* AZ      Z-directed main current-current matrices. Dimensionned to
*         MUZ(LL4Z)*NLF/2.
* C11W    W-directed main current-current matrices to be factorized.
*         Dimensionned to MUW(LL4W)*NLF/2.
* C11X    X-directed main current-current matrices to be factorized.
*         Dimensionned to MUX(LL4X)*NLF/2.
* C11Y    Y-directed main current-current matrices to be factorized.
*         Dimensionned to MUY(LL4Y)*NLF/2.
* C11Z    Z-directed main current-current matrices to be factorized.
*         Dimensionned to MUZ(LL4Z)*NLF/2.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,
     1 MAT(3,NBLOS),MUW(LL4W),IPBBW(2*IELEM,LL4W),LC,IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIGT(NBMIX,NAN),SIGTI(NBMIX,NAN),SIDE,ZZ(3,NBLOS),FRZ(NBLOS),
     1 QFR(NBLOS,8),R(LC,LC),V(LC,LC-1),BBW(2*IELEM,LL4W),
     2 TTF(LL4F*NLF/2),AW(*),C11W(*)
*----
*  LOCAL VARIABLES
*----
      REAL QQ(5,5)
      DOUBLE PRECISION FFF,TTTT,VOL0,GARS,GARSI,FACT,VAR1
*----
*  W-ORIENTED COUPLINGS
*----
      ZMARS=0.0
      IF(ICOL.EQ.3) THEN
         IF(NVD.EQ.0) THEN
            NZMAR=NLF+1
         ELSE IF(NVD.EQ.1) THEN
            NZMAR=NLF
         ELSE IF(NVD.EQ.2) THEN
            NZMAR=65
         ENDIF
      ELSE
         NZMAR=65
      ENDIF
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
      MUMAX=MUW(LL4W)
      DO 120 IL=0,NLF-1
      IF(MOD(IL,2).EQ.1) ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE W-ORIENTED COEFFICIENT MATRICES AT ORDER IL.
*----
      NELEH=(IELEM+1)*IELEM**2
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 70 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 70
      IBM=MAT(1,IPERT(KEL))
      NUM=NUM+1
      IF(IBM.EQ.0) GO TO 70
      DZ=ZZ(1,IPERT(KEL))
      VOL0=TTTT*DZ*FRZ(KEL)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION.
         VAR1=FACT*VOL0*GARS
         DO 32 K3=0,IELEM-1
         DO 31 K2=0,IELEM-1
         DO 30 K1=0,IELEM-1
         IOF=(IL/2)*LL4F
         JND1=IOF+(((NUM-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         JND2=IOF+(((KN(NUM,1)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         JND3=IOF+(((KN(NUM,2)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         TTF(JND1)=TTF(JND1)+REAL(VAR1)
         TTF(JND2)=TTF(JND2)+REAL(VAR1)
         TTF(JND3)=TTF(JND3)+REAL(VAR1)
   30    CONTINUE
   31    CONTINUE
   32    CONTINUE
      ELSE
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION OF
*        THE EVEN PARITY EQUATION.
         GARSI=SIGTI(IBM,MIN(IL+1,NAN))
         IF(IELEM.GT.1) THEN
            KOFF=((IL-1)/2)*LL4F
            DO 42 K3=0,IELEM-1
            DO 41 K2=0,IELEM-1
            DO 40 K1=0,IELEM-1
            JND1=KOFF+(((NUM-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
            JND2=KOFF+(((KN(NUM,1)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
            JND3=KOFF+(((KN(NUM,2)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
            VAR1=(REAL(IL)**2)*VOL0*QQ(K3+1,K3+1)*GARSI/(FACT*DZ*DZ)
            TTF(JND1)=TTF(JND1)+REAL(VAR1)
            TTF(JND2)=TTF(JND2)+REAL(VAR1)
            TTF(JND3)=TTF(JND3)+REAL(VAR1)
            IF(IL.LE.NLF-3) THEN
              JND1=JND1+LL4F
              JND2=JND2+LL4F
              JND3=JND3+LL4F
              VAR1=(REAL(IL+1)**2)*VOL0*QQ(K3+1,K3+1)*GARSI/(FACT*DZ*DZ)
              TTF(JND1)=TTF(JND1)+REAL(VAR1)
              TTF(JND2)=TTF(JND2)+REAL(VAR1)
              TTF(JND3)=TTF(JND3)+REAL(VAR1)
            ENDIF
   40       CONTINUE
   41       CONTINUE
   42       CONTINUE
         ENDIF
*
*        ODD PARITY EQUATION.
         DO 63 K5=0,1 ! TWO LOZENGES PER HEXAGON
         DO 62 K4=0,IELEM-1
         DO 61 K3=0,IELEM-1
         DO 60 K2=1,IELEM+1
         KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         INW1=ABS(KNW1)
         DO 50 K1=1,IELEM+1
         KNW2=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
         INW2=ABS(KNW2)
         IF((KNW2.NE.0).AND.(KNW1.NE.0).AND.(INW1.GE.INW2)) THEN
            KEY=((IL-1)/2)*MUMAX+MUW(INW1)-INW1+INW2
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNW2))
            VAR1=(4./3.)*SG*FACT*VOL0*GARS*R(K2,K1)
            AW(KEY)=AW(KEY)-REAL(VAR1)
         ENDIF
   50    CONTINUE
         IF(KNW1.NE.0) THEN
            KEY=((IL-1)/2)*MUMAX+MUW(INW1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,1)*ZMARS
               AW(KEY)=AW(KEY)-REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,2)*ZMARS
               AW(KEY)=AW(KEY)-REAL(VAR1)
            ENDIF
         ENDIF
   60    CONTINUE
   61    CONTINUE
   62    CONTINUE
   63    CONTINUE
      ENDIF
   70 CONTINUE
*
      IF(MOD(IL,2).EQ.1) THEN
         DO 80 I0=1,MUMAX
         C11W(((IL-1)/2)*MUMAX+I0)=-AW(((IL-1)/2)*MUMAX+I0)
   80    CONTINUE
         MUIM1=0
         DO 110 I=1,LL4W
         MUI=MUW(I)
         DO 100 J=I-(MUI-MUIM1)+1,I
         KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
         DO 95 I0=1,2*IELEM
         II=IPBBW(I0,I)
         IF(II.EQ.0) GO TO 100
         DO 90 J0=1,2*IELEM
         JJ=IPBBW(J0,J)
         IF(II.EQ.JJ) C11W(KEY)=C11W(KEY)+REAL(IL**2)*BBW(I0,I)*
     1   BBW(J0,J)/TTF(((IL-1)/2)*LL4F+II)
   90    CONTINUE
   95    CONTINUE
  100    CONTINUE
         MUIM1=MUI
  110    CONTINUE
      ENDIF
  120 CONTINUE
      RETURN
      END
*
      SUBROUTINE PN3HWX(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,
     1 LL4X,MAT,SIGT,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUX,IPBBX,LC,R,BBX,TTF,
     2 AX,C11X)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,LL4X,
     1 MAT(3,NBLOS),MUX(LL4X),IPBBX(2*IELEM,LL4X),LC,IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIGT(NBMIX,NAN),SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),
     1 R(LC,LC),BBX(2*IELEM,LL4X),TTF(LL4F*NLF/2),AX(*),C11X(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT,VOL0,GARS,FACT,VAR1
*----
*  X-ORIENTED COUPLINGS
*----
      ZMARS=0.0
      IF(ICOL.EQ.3) THEN
         IF(NVD.EQ.0) THEN
            NZMAR=NLF+1
         ELSE IF(NVD.EQ.1) THEN
            NZMAR=NLF
         ELSE IF(NVD.EQ.2) THEN
            NZMAR=65
         ENDIF
      ELSE
         NZMAR=65
      ENDIF
      MUMAX=MUX(LL4X)
      DO 200 IL=1,NLF-1,2
      ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE X-ORIENTED COEFFICIENT MATRICES AT ODD ORDER IL.
*----
      NELEH=(IELEM+1)*IELEM**2
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 150 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 150
      IBM=MAT(1,IPERT(KEL))
      NUM=NUM+1
      IF(IBM.EQ.0) GO TO 150
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
*
      DO 143 K5=0,1 ! TWO LOZENGES PER HEXAGON
      DO 142 K4=0,IELEM-1
      DO 141 K3=0,IELEM-1
      DO 140 K2=1,IELEM+1
      KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INX1=ABS(KNX1)-LL4W
      DO 130 K1=1,IELEM+1
      KNX2=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      INX2=ABS(KNX2)-LL4W
      IF((KNX2.NE.0).AND.(KNX1.NE.0).AND.(INX1.GE.INX2)) THEN
         KEY=((IL-1)/2)*MUMAX+MUX(INX1)-INX1+INX2
         SG=REAL(SIGN(1,KNX1)*SIGN(1,KNX2))
         VAR1=(4./3.)*SG*FACT*VOL0*GARS*R(K2,K1)
         AX(KEY)=AX(KEY)-REAL(VAR1)
      ENDIF
  130 CONTINUE
      IF(KNX1.NE.0) THEN
         KEY=((IL-1)/2)*MUMAX+MUX(INX1)
         IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
            VAR1=0.5*FACT*QFR(NUM,3)*ZMARS
            AX(KEY)=AX(KEY)-REAL(VAR1)
         ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
            VAR1=0.5*FACT*QFR(NUM,4)*ZMARS
            AX(KEY)=AX(KEY)-REAL(VAR1)
         ENDIF
      ENDIF
  140 CONTINUE
  141 CONTINUE
  142 CONTINUE
  143 CONTINUE
  150 CONTINUE
*
      DO 160 I0=1,MUMAX
      C11X(((IL-1)/2)*MUMAX+I0)=-AX(((IL-1)/2)*MUMAX+I0)
  160 CONTINUE
      MUIM1=0
      DO 190 I=1,LL4X
      MUI=MUX(I)
      DO 180 J=I-(MUI-MUIM1)+1,I
      KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
      DO 175 I0=1,2*IELEM
      II=IPBBX(I0,I)
      IF(II.EQ.0) GO TO 180
      DO 170 J0=1,2*IELEM
      JJ=IPBBX(J0,J)
      IF(II.EQ.JJ) THEN
         VAR1=REAL(IL**2)*BBX(I0,I)*BBX(J0,J)/TTF(((IL-1)/2)*LL4F+II)
         C11X(KEY)=C11X(KEY)+REAL(VAR1)
      ENDIF
  170 CONTINUE
  175 CONTINUE
  180 CONTINUE
      MUIM1=MUI
  190 CONTINUE
  200 CONTINUE
      RETURN
      END
*
      SUBROUTINE PN3HWY(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,
     1 LL4X,LL4Y,MAT,SIGT,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUY,IPBBY,LC,R,BBY,
     2 TTF,AY,C11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,LL4X,LL4Y,
     1 MAT(3,NBLOS),MUY(LL4Y),IPBBY(2*IELEM,LL4Y),LC,IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIGT(NBMIX,NAN),SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),
     1 R(LC,LC),BBY(2*IELEM,LL4Y),TTF(LL4F*NLF/2),AY(*),C11Y(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT,VOL0,GARS,FACT,VAR1
*----
*  Y-ORIENTED COUPLINGS
*----
      ZMARS=0.0
      IF(ICOL.EQ.3) THEN
         IF(NVD.EQ.0) THEN
            NZMAR=NLF+1
         ELSE IF(NVD.EQ.1) THEN
            NZMAR=NLF
         ELSE IF(NVD.EQ.2) THEN
            NZMAR=65
         ENDIF
      ELSE
         NZMAR=65
      ENDIF
      MUMAX=MUY(LL4Y)
      DO 280 IL=1,NLF-1,2
      ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE Y-ORIENTED COEFFICIENT MATRICES AT ODD ORDER IL.
*----
      NELEH=(IELEM+1)*IELEM**2
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 230 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 230
      IBM=MAT(1,IPERT(KEL))
      NUM=NUM+1
      IF(IBM.EQ.0) GO TO 230
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
*
      DO 223 K5=0,1 ! TWO LOZENGES PER HEXAGON
      DO 222 K4=0,IELEM-1
      DO 221 K3=0,IELEM-1
      DO 220 K2=1,IELEM+1
      KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
      INY1=ABS(KNY1)-LL4W-LL4X
      DO 210 K1=1,IELEM+1
      KNY2=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
      INY2=ABS(KNY2)-LL4W-LL4X
      IF((KNY2.NE.0).AND.(KNY1.NE.0).AND.(INY1.GE.INY2)) THEN
         KEY=((IL-1)/2)*MUMAX+MUY(INY1)-INY1+INY2
         SG=REAL(SIGN(1,KNY1)*SIGN(1,KNY2))
         VAR1=(4./3.)*SG*FACT*VOL0*GARS*R(K2,K1)
         AY(KEY)=AY(KEY)-REAL(VAR1)
      ENDIF
  210 CONTINUE
      IF(KNY1.NE.0) THEN
         KEY=((IL-1)/2)*MUMAX+MUY(INY1)
         IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
            VAR1=0.5*FACT*QFR(NUM,5)*ZMARS
            AY(KEY)=AY(KEY)-REAL(VAR1)
         ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
            VAR1=0.5*FACT*QFR(NUM,6)*ZMARS
            AY(KEY)=AY(KEY)-REAL(VAR1)
         ENDIF
      ENDIF
  220 CONTINUE
  221 CONTINUE
  222 CONTINUE
  223 CONTINUE
  230 CONTINUE
*
      DO 240 I0=1,MUMAX
      C11Y(((IL-1)/2)*MUMAX+I0)=-AY(((IL-1)/2)*MUMAX+I0)
  240 CONTINUE
      MUIM1=0
      DO 270 I=1,LL4Y
      MUI=MUY(I)
      DO 260 J=I-(MUI-MUIM1)+1,I
      KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
      DO 255 I0=1,2*IELEM
      II=IPBBY(I0,I)
      IF(II.EQ.0) GO TO 260
      DO 250 J0=1,2*IELEM
      JJ=IPBBY(J0,J)
      IF(II.EQ.JJ) C11Y(KEY)=C11Y(KEY)+REAL(IL**2)*BBY(I0,I)*
     1 BBY(J0,J)/TTF(((IL-1)/2)*LL4F+II)
  250 CONTINUE
  255 CONTINUE
  260 CONTINUE
      MUIM1=MUI
  270 CONTINUE
  280 CONTINUE
      RETURN
      END
*
      SUBROUTINE PN3HWZ(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,
     1 LL4X,LL4Y,LL4Z,MAT,SIGT,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUZ,IPBBZ,LC,
     2 R,BBZ,TTF,AZ,C11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,LL4X,
     1 LL4Y,LL4Z,MAT(3,NBLOS),MUZ(LL4Z),IPBBZ(2*IELEM,LL4Z),LC,
     2 IPERT(NBLOS),KN(NBLOS,3+6*(IELEM+2)*IELEM**2)
      REAL SIGT(NBMIX,NAN),SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),
     1 R(LC,LC),BBZ(2*IELEM,LL4Z),TTF(LL4F*NLF/2),AZ(*),C11Z(*)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION TTTT,VOL0,GARS,FACT,VAR1
*----
*  Z-ORIENTED COUPLINGS
*----
      ZMARS=0.0
      IF(ICOL.EQ.3) THEN
         IF(NVD.EQ.0) THEN
            NZMAR=NLF+1
         ELSE IF(NVD.EQ.1) THEN
            NZMAR=NLF
         ELSE IF(NVD.EQ.2) THEN
            NZMAR=65
         ENDIF
      ELSE
         NZMAR=65
      ENDIF
      MUMAX=MUZ(LL4Z)
      DO 360 IL=1,NLF-1,2
      ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE Z-ORIENTED COEFFICIENT MATRICES AT ODD ORDER IL.
*----
      NELEH=(IELEM+1)*IELEM**2
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 310 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 310
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 310
      NUM=NUM+1
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
*
      DO 302 K5=0,2 ! THREE LOZENGES PER HEXAGON
      DO 301 K2=0,IELEM-1
      DO 300 K1=0,IELEM-1
      KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
      INZ1=ABS(KNZ1)-LL4W-LL4X-LL4Y
      KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
      INZ2=ABS(KNZ2)-LL4W-LL4X-LL4Y
      IF(KNZ1.NE.0) THEN
        KEY=((IL-1)/2)*MUMAX+MUZ(INZ1)
        VAR1=FACT*VOL0*GARS*R(1,1)+0.5*FACT*QFR(NUM,7)*ZMARS
        AZ(KEY)=AZ(KEY)-REAL(VAR1)
      ENDIF
      IF(KNZ2.NE.0) THEN
        KEY=((IL-1)/2)*MUMAX+MUZ(INZ2)
        VAR1=FACT*VOL0*GARS*R(IELEM+1,IELEM+1)+0.5*FACT*QFR(NUM,8)*ZMARS
        AZ(KEY)=AZ(KEY)-REAL(VAR1)
      ENDIF
      IF((ICOL.NE.2).AND.(KNZ1.NE.0).AND.(KNZ2.NE.0)) THEN
        IF(INZ2.GT.INZ1) KEY=((IL-1)/2)*MUMAX+MUZ(INZ2)-INZ2+INZ1
        IF(INZ2.LE.INZ1) KEY=((IL-1)/2)*MUMAX+MUZ(INZ1)-INZ1+INZ2
        SG=REAL(SIGN(1,KNZ1)*SIGN(1,KNZ2))
        IF(INZ1.EQ.INZ2) SG=2.0*SG
        VAR1=SG*FACT*VOL0*GARS*R(IELEM+1,1)
        AZ(KEY)=AZ(KEY)-REAL(VAR1)
      ENDIF
  300 CONTINUE
  301 CONTINUE
  302 CONTINUE
  310 CONTINUE
*
      DO 320 I0=1,MUMAX
      C11Z(((IL-1)/2)*MUMAX+I0)=-AZ(((IL-1)/2)*MUMAX+I0)
  320 CONTINUE
      MUIM1=0
      DO 350 I=1,LL4Z
      MUI=MUZ(I)
      DO 340 J=I-(MUI-MUIM1)+1,I
      KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
      DO 335 I0=1,2*IELEM
      II=IPBBZ(I0,I)
      IF(II.EQ.0) GO TO 340
      DO 330 J0=1,2*IELEM
      JJ=IPBBZ(J0,J)
      IF(II.EQ.JJ) C11Z(KEY)=C11Z(KEY)+REAL(IL**2)*BBZ(I0,I)*
     1 BBZ(J0,J)/TTF(((IL-1)/2)*LL4F+II)
  330 CONTINUE
  335 CONTINUE
  340 CONTINUE
      MUIM1=MUI
  350 CONTINUE
  360 CONTINUE
      RETURN
      END
