*DECK PN3DXX
      SUBROUTINE PN3DXX(NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,
     1 SIGT,SIGTI,MAT,VOL,XX,YY,ZZ,KN,QFR,MUX,IPBBX,LC,R,V,BBX,TTF,
     2 AX,C11X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a Thomas-Raviart (dual) finite element
* method in 3-D simplified PN approximation. Note: system matrices
* should be initialized by the calling program.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NBMIX   number of mixtures.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NEL     total number of finite elements.
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* LL4F    number of flux components.
* LL4X    number of X-directed currents.
* LL4Y    number of Y-directed currents.
* LL4Z    number of Z-directed currents.
* SIGT    total minus self-scattering macroscopic cross sections.
*         SIGT(:,NAN) generally contains the total cross section only.
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
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
*         MUX(LL4X)*NLF/2.
* AY      Y-directed main current-current matrices. Dimensionned to
*         MUY(LL4Y)*NLF/2.
* AZ      Z-directed main current-current matrices. Dimensionned to
*         MUZ(LL4Z)*NLF/2.
* C11X    X-directed main current-current matrices to be factorized.
*         Dimensionned to MUX(LL4X)*NLF/2.
* C11Y    Y-directed main current-current matrices to be factorized.
*         Dimensionned to MUY(LL4Y)*NLF/2.
* C11Z    Z-directed main current-current matrices to be factorized.
*         Dimensionned to MUZ(LL4Z)*NLF/2.
*
*Reference(s):
* J.J. Lautard, D. Schneider, A.M. Baudron, "Mixed Dual Methods for
* Neutronic Reactor Core Calculations in the CRONOS System," Proc.
* Int. Conf. on Mathematics and Computation, Reactor Physics and
* Environmental Analysis in Nuclear Applications, Madrid, Spain,
* September 27-30, 1999.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,MAT(NEL),
     1 KN(NEL*(1+6*IELEM**2)),MUX(LL4X),IPBBX(2*IELEM,LL4X),LC
      REAL SIGT(NBMIX,NAN),SIGTI(NBMIX,NAN),VOL(NEL),XX(NEL),YY(NEL),
     1 ZZ(NEL),QFR(6*NEL),R(LC,LC),V(LC,LC-1),BBX(2*IELEM,LL4X),
     2 TTF(LL4F*NLF/2),AX(*),C11X(*)
*----
*  LOCAL VARIABLES
*----
      REAL QQ(5,5)
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
      DO 25 I0=1,IELEM
      DO 20 J0=1,IELEM
      FFF=0.0
      DO 10 K0=2,IELEM
      FFF=FFF+V(K0,I0)*V(K0,J0)/R(K0,K0)
   10 CONTINUE
      IF(ABS(FFF).LE.1.0E-6) FFF=0.0
      QQ(I0,J0)=FFF
   20 CONTINUE
   25 CONTINUE
      MUMAX=MUX(LL4X)
      DO 170 IL=0,NLF-1
      IF(MOD(IL,2).EQ.1) ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE X-ORIENTED COEFFICIENT MATRICES AT ORDER IL.
*----
      NUM1=0
      NUM2=0
      DO 120 IE=1,NEL
      IBM=MAT(IE)
      IF(IBM.EQ.0) GO TO 120
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 110
      DX=XX(IE)
      DY=YY(IE)
      DZ=ZZ(IE)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION.
         DO 32 K3=0,IELEM-1
         DO 31 K2=0,IELEM-1
         DO 30 K1=0,IELEM-1
         KEY=(IL/2)*LL4F+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
         TTF(KEY)=TTF(KEY)+FACT*VOL0*GARS
   30    CONTINUE
   31    CONTINUE
   32    CONTINUE
      ELSE
         GARSI=SIGTI(IBM,MIN(IL+1,NAN))
         DO 105 K3=0,IELEM-1
         DO 100 K2=0,IELEM-1
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION OF
*        THE EVEN PARITY EQUATION.
         DO 40 K1=0,IELEM-1
         JND1=KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
         KEY=((IL-1)/2)*LL4F+JND1
         TTF(KEY)=TTF(KEY)+(REAL(IL)**2)*VOL0*QQ(K1+1,K1+1)*GARSI/(FACT*
     1   DX*DX)
         IF(IL.LE.NLF-3) THEN
            KEY=((IL+2)/2)*LL4F+JND1
            TTF(KEY)=TTF(KEY)+(REAL(IL+1)**2)*VOL0*QQ(K1+1,K1+1)*GARSI/
     1      (FACT*DX*DX)
         ENDIF
         KEY=((IL-1)/2)*LL4F+JND1
         TTF(KEY)=TTF(KEY)+(REAL(IL)**2)*VOL0*QQ(K2+1,K2+1)*GARSI/
     1   (FACT*DY*DY)
         IF(IL.LE.NLF-3) THEN
            KEY=((IL+2)/2)*LL4F+JND1
            TTF(KEY)=TTF(KEY)+(REAL(IL+1)**2)*VOL0*QQ(K2+1,K2+1)*GARSI/
     1      (FACT*DY*DY)
         ENDIF
         KEY=((IL-1)/2)*LL4F+JND1
         TTF(KEY)=TTF(KEY)+(REAL(IL)**2)*VOL0*QQ(K3+1,K3+1)*GARSI/
     1   (FACT*DZ*DZ)
         IF(IL.LE.NLF-3) THEN
            KEY=((IL+2)/2)*LL4F+JND1
            TTF(KEY)=TTF(KEY)+(REAL(IL+1)**2)*VOL0*QQ(K3+1,K3+1)*GARSI/
     1      (FACT*DZ*DZ)
         ENDIF
   40    CONTINUE
*
*        ODD PARITY EQUATION.
         DO 55 IC=1,2
         IF(IC.EQ.1) IIC=1
         IF(IC.EQ.2) IIC=IELEM+1
         KN1=KN(NUM1+2+(IC-1)*IELEM**2+K3*IELEM+K2)
         IND1=ABS(KN1)-LL4F
         S1=REAL(SIGN(1,KN1))
         DO 50 JC=1,2
         IF(JC.EQ.1) JJC=1
         IF(JC.EQ.2) JJC=IELEM+1
         KN2=KN(NUM1+2+(JC-1)*IELEM**2+K3*IELEM+K2)
         IND2=ABS(KN2)-LL4F
         IF((KN1.NE.0).AND.(KN2.NE.0).AND.(IND1.GE.IND2)) THEN
            S2=REAL(SIGN(1,KN2))
            KEY=((IL-1)/2)*MUMAX+MUX(IND1)-IND1+IND2
            AX(KEY)=AX(KEY)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS
         ENDIF
   50    CONTINUE
   55    CONTINUE
*
         KN1=KN(NUM1+2+K3*IELEM+K2)
         KN2=KN(NUM1+2+IELEM**2+K3*IELEM+K2)
         IND1=ABS(KN1)-LL4F
         IND2=ABS(KN2)-LL4F
         IF((QFR(NUM2+1).NE.0.0).AND.(KN1.NE.0)) THEN
            KEY=((IL-1)/2)*MUMAX+MUX(IND1)
            AX(KEY)=AX(KEY)-0.5*FACT*QFR(NUM2+1)*ZMARS
         ENDIF
         IF((QFR(NUM2+2).NE.0.0).AND.(KN2.NE.0)) THEN
            KEY=((IL-1)/2)*MUMAX+MUX(IND2)
            AX(KEY)=AX(KEY)-0.5*FACT*QFR(NUM2+2)*ZMARS
         ENDIF
  100    CONTINUE
  105    CONTINUE
      ENDIF
  110 NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  120 CONTINUE
*
      IF(MOD(IL,2).EQ.1) THEN
         DO 130 I0=1,MUMAX
         C11X(((IL-1)/2)*MUMAX+I0)=-AX(((IL-1)/2)*MUMAX+I0)
  130    CONTINUE
         MUIM1=0
         DO 160 I=1,LL4X
         MUI=MUX(I)
         DO 150 J=I-(MUI-MUIM1)+1,I
         KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
         DO 145 I0=1,2*IELEM
         II=IPBBX(I0,I)
         IF(II.EQ.0) GO TO 150
         DO 140 J0=1,2*IELEM
         JJ=IPBBX(J0,J)
         IF(II.EQ.JJ) C11X(KEY)=C11X(KEY)+REAL(IL**2)*BBX(I0,I)*
     1   BBX(J0,J)/TTF(((IL-1)/2)*LL4F+II)
  140    CONTINUE
  145    CONTINUE
  150    CONTINUE
         MUIM1=MUI
  160    CONTINUE
      ENDIF
  170 CONTINUE
      RETURN
      END
*
      SUBROUTINE PN3DXY(NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,LL4Y,
     1 SIGT,MAT,VOL,YY,KN,QFR,MUY,IPBBY,LC,R,BBY,TTF,AY,C11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,LL4Y,MAT(NEL),
     1 KN(NEL*(1+6*IELEM**2)),MUY(LL4Y),IPBBY(2*IELEM,LL4Y),LC
      REAL SIGT(NBMIX,NAN),VOL(NEL),YY(NEL),QFR(6*NEL),R(LC,LC),
     1 BBY(2*IELEM,LL4Y),TTF(LL4F*NLF/2),AY(*),C11Y(*)
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
      DO 320 IL=1,NLF-1,2
      ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE Y-ORIENTED COEFFICIENT MATRICES AT ODD ORDER IL.
*----
      NUM1=0
      NUM2=0
      DO 270 IE=1,NEL
      IBM=MAT(IE)
      IF(IBM.EQ.0) GO TO 270
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 260
      DY=YY(IE)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
*
      DO 255 K3=0,IELEM-1
      DO 250 K1=0,IELEM-1
      DO 205 IC=3,4
      IF(IC.EQ.3) IIC=1
      IF(IC.EQ.4) IIC=IELEM+1
      KN1=KN(NUM1+2+(IC-1)*IELEM**2+K3*IELEM+K1)
      IND1=ABS(KN1)-LL4F-LL4X
      S1=REAL(SIGN(1,KN1))
      DO 200 JC=3,4
      IF(JC.EQ.3) JJC=1
      IF(JC.EQ.4) JJC=IELEM+1
      KN2=KN(NUM1+2+(JC-1)*IELEM**2+K3*IELEM+K1)
      IND2=ABS(KN2)-LL4F-LL4X
      IF((KN1.NE.0).AND.(KN2.NE.0).AND.(IND1.GE.IND2)) THEN
         S2=REAL(SIGN(1,KN2))
         KEY=((IL-1)/2)*MUMAX+MUY(IND1)-IND1+IND2
         AY(KEY)=AY(KEY)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS
      ENDIF
  200 CONTINUE
  205 CONTINUE
*
      KN1=KN(NUM1+2+2*IELEM**2+K3*IELEM+K1)
      KN2=KN(NUM1+2+3*IELEM**2+K3*IELEM+K1)
      IND1=ABS(KN1)-LL4F-LL4X
      IND2=ABS(KN2)-LL4F-LL4X
      IF((QFR(NUM2+3).NE.0.0).AND.(KN1.NE.0)) THEN
         KEY=((IL-1)/2)*MUMAX+MUY(IND1)
         AY(KEY)=AY(KEY)-0.5*FACT*QFR(NUM2+3)*ZMARS
      ENDIF
      IF((QFR(NUM2+4).NE.0.0).AND.(KN2.NE.0)) THEN
         KEY=((IL-1)/2)*MUMAX+MUY(IND2)
         AY(KEY)=AY(KEY)-0.5*FACT*QFR(NUM2+4)*ZMARS
      ENDIF
  250 CONTINUE
  255 CONTINUE
  260 NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  270 CONTINUE
*
      DO 280 I0=1,MUMAX
      C11Y(((IL-1)/2)*MUMAX+I0)=-AY(((IL-1)/2)*MUMAX+I0)
  280 CONTINUE
      MUIM1=0
      DO 310 I=1,LL4Y
      MUI=MUY(I)
      DO 300 J=I-(MUI-MUIM1)+1,I
      KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
      DO 295 I0=1,2*IELEM
      II=IPBBY(I0,I)
      IF(II.EQ.0) GO TO 300
      DO 290 J0=1,2*IELEM
      JJ=IPBBY(J0,J)
      IF(II.EQ.JJ) C11Y(KEY)=C11Y(KEY)+REAL(IL**2)*BBY(I0,I)*BBY(J0,J)/
     1 TTF(((IL-1)/2)*LL4F+II)
  290 CONTINUE
  295 CONTINUE
  300 CONTINUE
      MUIM1=MUI
  310 CONTINUE
  320 CONTINUE
      RETURN
      END
*
      SUBROUTINE PN3DXZ(NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,LL4Y,
     1 LL4Z,SIGT,MAT,VOL,ZZ,KN,QFR,MUZ,IPBBZ,LC,R,BBZ,TTF,AZ,C11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,LL4Y,LL4Z,
     1 MAT(NEL),KN(NEL*(1+6*IELEM**2)),MUZ(LL4Z),IPBBZ(2*IELEM,LL4Z),LC
      REAL SIGT(NBMIX,NAN),VOL(NEL),ZZ(NEL),QFR(6*NEL),R(LC,LC),
     1 BBZ(2*IELEM,LL4Z),TTF(LL4F*NLF/2),AZ(*),C11Z(*)
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
      DO 470 IL=1,NLF-1,2
      ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE Z-ORIENTED COEFFICIENT MATRICES AT ORDER IL.
*----
      NUM1=0
      NUM2=0
      DO 420 IE=1,NEL
      IBM=MAT(IE)
      IF(IBM.EQ.0) GO TO 420
      VOL0=VOL(IE)
      IF(VOL0.EQ.0.0) GO TO 410
      DZ=ZZ(IE)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
*
      DO 405 K2=0,IELEM-1
      DO 400 K1=0,IELEM-1
      DO 355 IC=5,6
      IF(IC.EQ.5) IIC=1
      IF(IC.EQ.6) IIC=IELEM+1
      KN1=KN(NUM1+2+(IC-1)*IELEM**2+K2*IELEM+K1)
      IND1=ABS(KN1)-LL4F-LL4X-LL4Y
      S1=REAL(SIGN(1,KN1))
      DO 350 JC=5,6
      IF(JC.EQ.5) JJC=1
      IF(JC.EQ.6) JJC=IELEM+1
      KN2=KN(NUM1+2+(JC-1)*IELEM**2+K2*IELEM+K1)
      IND2=ABS(KN2)-LL4F-LL4X-LL4Y
      IF((KN1.NE.0).AND.(KN2.NE.0).AND.(IND1.GE.IND2)) THEN
         S2=REAL(SIGN(1,KN2))
         KEY=((IL-1)/2)*MUMAX+MUZ(IND1)-IND1+IND2
         AZ(KEY)=AZ(KEY)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS
      ENDIF
  350 CONTINUE
  355 CONTINUE
*
      KN1=KN(NUM1+2+4*IELEM**2+K2*IELEM+K1)
      KN2=KN(NUM1+2+5*IELEM**2+K2*IELEM+K1)
      IND1=ABS(KN1)-LL4F-LL4X-LL4Y
      IND2=ABS(KN2)-LL4F-LL4X-LL4Y
      IF((QFR(NUM2+5).NE.0.0).AND.(KN1.NE.0)) THEN
         KEY=((IL-1)/2)*MUMAX+MUZ(IND1)
         AZ(KEY)=AZ(KEY)-0.5*FACT*QFR(NUM2+5)*ZMARS
      ENDIF
      IF((QFR(NUM2+6).NE.0.0).AND.(KN2.NE.0)) THEN
         KEY=((IL-1)/2)*MUMAX+MUZ(IND2)
         AZ(KEY)=AZ(KEY)-0.5*FACT*QFR(NUM2+6)*ZMARS
      ENDIF
  400 CONTINUE
  405 CONTINUE
  410 NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  420 CONTINUE
*
      DO 430 I0=1,MUMAX
      C11Z(((IL-1)/2)*MUMAX+I0)=-AZ(((IL-1)/2)*MUMAX+I0)
  430 CONTINUE
      MUIM1=0
      DO 460 I=1,LL4Z
      MUI=MUZ(I)
      DO 450 J=I-(MUI-MUIM1)+1,I
      KEY=((IL-1)/2)*MUMAX+(MUI-I+J)
      DO 445 I0=1,2*IELEM
      II=IPBBZ(I0,I)
      IF(II.EQ.0) GO TO 450
      DO 440 J0=1,2*IELEM
      JJ=IPBBZ(J0,J)
      IF(II.EQ.JJ) C11Z(KEY)=C11Z(KEY)+REAL(IL**2)*BBZ(I0,I)*BBZ(J0,J)/
     1 TTF(((IL-1)/2)*LL4F+II)
  440 CONTINUE
  445 CONTINUE
  450 CONTINUE
      MUIM1=MUI
  460 CONTINUE
  470 CONTINUE
      RETURN
      END
