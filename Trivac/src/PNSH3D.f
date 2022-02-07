*DECK PNSH3D
      SUBROUTINE PNSH3D (ITY,IPR,NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,L4,
     1 LL4F,LL4W,LL4X,LL4Y,MAT,SIGT,SIGTI,SIDE,ZZ,FRZ,QFR,IPERT,KN,LC,
     2 R,V,CTRAN,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source calculation for a SPN approximation in TRIVAC, including
* neighbour Legendre and out-of-group contributions.
* Raviart-Thomas-Schneider method in hexagonal geometry.
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
* ITY     type of assembly:
*         =0: leakage-removal matrix assembly; =1: cross section matrix
*         assembly.
* IPR     type of assembly:
*         =0: contains system matrices;
*         =1: contains derivative of these matrices;
*         =2: contains first variation of these matrices;
*         =3: contains addition of first vatiation to unperturbed
*         system matrices.
* NBMIX   number of mixtures.
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* LL4F    number of flux components.
* LL4W    number of W-directed currents.
* LL4X    number of X-directed currents.
* LL4Y    number of Y-directed currents.
* MAT     index-number of the mixture type assigned to each volume.
* SIGT    macroscopic cross sections ordered by mixture.
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
* LC      order of the unit matrices.
* R       unit Cartesian mass matrix.
* V       unit nodal coupling matrix.
* CTRAN   tranverse coupling Piolat unit matrix.
* FUNKNO  initial fluxes.
* SUNKNO  initial sources.
*
*Parameters: output
* SUNKNO  modified sources.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,IPR,NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,L4,LL4F,LL4W,
     1 LL4X,LL4Y,MAT(3,NBLOS),IPERT(NBLOS),
     2 KN(NBLOS,3+6*(IELEM+2)*IELEM**2),LC
      REAL SIGT(NBMIX,NAN),SIGTI(NBMIX,NAN),SIDE,ZZ(3,NBLOS),
     1 FRZ(NBLOS),QFR(NBLOS,8),R(LC,LC),V(LC,LC-1),SUNKNO(L4*NLF/2),
     2 FUNKNO(L4*NLF/2)
      DOUBLE PRECISION CTRAN((IELEM+1)*IELEM,(IELEM+1)*IELEM)
*----
*  LOCAL VARIABLES
*----
      REAL QQ(5,5)
      DOUBLE PRECISION FFF,TTTT,UUUU,VOL0,GARS,GARSI,FACT,VAR1
      REAL, DIMENSION(:), ALLOCATABLE :: DIFF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DIFF(NBLOS))
*
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
      DO 16 I0=1,IELEM
      DO 15 J0=1,IELEM
      FFF=0.0D0
      DO 10 K0=2,IELEM
      FFF=FFF+V(K0,I0)*V(K0,J0)/R(K0,K0)
   10 CONTINUE
      IF(ABS(FFF).LE.1.0E-6) FFF=0.0D0
      QQ(I0,J0)=REAL(FFF)
   15 CONTINUE
   16 CONTINUE
*----
*  MAIN LOOP OVER LEGENDRE ORDERS FOR THE FLUX.
*----
      DO 200 IL=0,NLF-1
      IF((ITY.EQ.1).AND.(IL.GE.NAN)) GO TO 200
      FACT=REAL(2*IL+1)
*----
*  RECOVER CROSS SECTIONS FOR THE PIOLAT TERMS.
*----
      IF(MOD(IL,2).EQ.1) THEN
         DO 20 KEL=1,NBLOS
         DIFF(KEL)=0.0
         IF(IPERT(KEL).GT.0) THEN
            IBM=MAT(1,IPERT(KEL))
            IF(IBM.GT.0) THEN
               GARS=SIGT(IBM,MIN(IL+1,NAN))
               VAR1=FACT*ZZ(1,IPERT(KEL))*FRZ(KEL)*GARS
               DIFF(KEL)=REAL(VAR1)
            ENDIF
         ENDIF
   20    CONTINUE
      ENDIF
*----
*  COMPUTE THE SOURCE AT ORDER IL.
*----
      NELEH=(IELEM+1)*IELEM**2
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      JOFF=(IL/2)*L4
      NUM=0
      DO 180 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 180
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 180
      DZ=ZZ(1,IPERT(KEL))
      VOL0=TTTT*DZ*FRZ(KEL)
      UUUU=SIDE*DZ*FRZ(KEL)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION
         DO 27 K3=0,IELEM-1
         DO 26 K2=0,IELEM-1
         DO 25 K1=0,IELEM-1
         JND1=JOFF+(((NUM-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         JND2=JOFF+(((KN(NUM,1)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         JND3=JOFF+(((KN(NUM,2)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         SUNKNO(JND1)=SUNKNO(JND1)+REAL(FACT*VOL0*GARS*FUNKNO(JND1))
         SUNKNO(JND2)=SUNKNO(JND2)+REAL(FACT*VOL0*GARS*FUNKNO(JND2))
         SUNKNO(JND3)=SUNKNO(JND3)+REAL(FACT*VOL0*GARS*FUNKNO(JND3))
   25    CONTINUE
   26    CONTINUE
   27    CONTINUE
         IF(ITY.EQ.1) GO TO 180
         IF((IPR.EQ.1).OR.(IPR.EQ.2)) GO TO 180
*
         DO 34 K5=0,1 ! TWO LOZENGES PER HEXAGON
         DO 33 K4=0,IELEM-1
         DO 32 K3=0,IELEM-1
         DO 31 K2=1,IELEM+1
         KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         INW1=JOFF+LL4F+ABS(KNW1)
         INX1=JOFF+LL4F+ABS(KNX1)
         INY1=JOFF+LL4F+ABS(KNY1)
         DO 30 K1=0,IELEM-1
         IF(V(K2,K1+1).EQ.0.0) GO TO 30
         IF(K5.EQ.0) THEN
            SSS=(-1.0)**K1
            JND1=JOFF+(((NUM-1)*IELEM+K4)*IELEM+K3)*IELEM+K1+1
            JND2=JOFF+(((KN(NUM,1)-1)*IELEM+K4)*IELEM+K3)*IELEM+K1+1
            JND3=JOFF+(((KN(NUM,2)-1)*IELEM+K4)*IELEM+K3)*IELEM+K1+1
         ELSE
            SSS=1.0
            JND1=JOFF+(((KN(NUM,1)-1)*IELEM+K4)*IELEM+K1)*IELEM+K3+1
            JND2=JOFF+(((KN(NUM,2)-1)*IELEM+K4)*IELEM+K1)*IELEM+K3+1
            JND3=JOFF+(((KN(NUM,3)-1)*IELEM+K4)*IELEM+K1)*IELEM+K3+1
         ENDIF
         IF(KNW1.NE.0) THEN
            SG=REAL(SIGN(1,KNW1))
            VAR1=SG*SSS*REAL(IL+1)*UUUU*V(K2,K1+1)*FUNKNO(INW1)
            SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)
         ENDIF
         IF(KNX1.NE.0) THEN
            SG=REAL(SIGN(1,KNX1))
            VAR1=SG*SSS*REAL(IL+1)*UUUU*V(K2,K1+1)*FUNKNO(INX1)
            SUNKNO(JND2)=SUNKNO(JND2)+REAL(VAR1)
         ENDIF
         IF(KNY1.NE.0) THEN
            SG=REAL(SIGN(1,KNY1))
            VAR1=SG*SSS*REAL(IL+1)*UUUU*V(K2,K1+1)*FUNKNO(INY1)
            SUNKNO(JND3)=SUNKNO(JND3)+REAL(VAR1)
         ENDIF
         IF(IL.GE.2) THEN
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               VAR1=SG*SSS*REAL(IL)*UUUU*V(K2,K1+1)*FUNKNO(INW1-L4)
               SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               VAR1=SG*SSS*REAL(IL)*UUUU*V(K2,K1+1)*FUNKNO(INX1-L4)
               SUNKNO(JND2)=SUNKNO(JND2)+REAL(VAR1)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               VAR1=SG*SSS*REAL(IL)*UUUU*V(K2,K1+1)*FUNKNO(INY1-L4)
               SUNKNO(JND3)=SUNKNO(JND3)+REAL(VAR1)
            ENDIF
         ENDIF
   30    CONTINUE
   31    CONTINUE
   32    CONTINUE
   33    CONTINUE
   34    CONTINUE
         DO 43 K5=0,2 ! THREE LOZENGES PER HEXAGON
         DO 42 K2=0,IELEM-1
         DO 41 K1=0,IELEM-1
         KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
         KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
         INZ1=JOFF+LL4F+ABS(KNZ1)
         INZ2=JOFF+LL4F+ABS(KNZ2)
         DO 40 K3=0,IELEM-1
         IF(K5.EQ.0) THEN
            JND1=JOFF+((((NUM-1)*IELEM)+K3)*IELEM+K2)*IELEM+K1+1
         ELSE
            JND1=JOFF+(((KN(NUM,K5)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         ENDIF
         IF(KNZ1.NE.0) THEN
            SG=REAL(SIGN(1,KNZ1))
            VAR1=SG*(VOL0/DZ)*REAL(IL+1)*V(1,K3+1)*FUNKNO(INZ1)
            SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)
         ENDIF
         IF(KNZ2.NE.0) THEN
            SG=REAL(SIGN(1,KNZ2))
            VAR1=SG*(VOL0/DZ)*REAL(IL+1)*V(IELEM+1,K3+1)*FUNKNO(INZ2)
            SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)
         ENDIF
         IF(IL.GE.2) THEN
            IF(KNZ1.NE.0) THEN
               SG=REAL(SIGN(1,KNZ1))
               VAR1=SG*(VOL0/DZ)*REAL(IL)*V(1,K3+1)*FUNKNO(INZ1-L4)
               SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)
            ENDIF
            IF(KNZ2.NE.0) THEN
               SG=REAL(SIGN(1,KNZ2))
               VAR1=SG*(VOL0/DZ)*REAL(IL)*V(IELEM+1,K3+1)*
     1         FUNKNO(INZ2-L4)
               SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)
            ENDIF
         ENDIF
   40    CONTINUE
   41    CONTINUE
   42    CONTINUE
   43    CONTINUE
      ELSE IF(MOD(IL,2).EQ.1) THEN
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION OF
*        THE EVEN PARITY EQUATION.
         GARSI=SIGTI(IBM,MIN(IL+1,NAN))
         IF(IELEM.GT.1) THEN
            DO 72 K3=0,IELEM-1
            DO 71 K2=0,IELEM-1
            DO 70 K1=0,IELEM-1
            IF(QQ(K3+1,K3+1).EQ.0.0) GO TO 70
            JND1=JOFF+(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            JND2=JOFF+(KN(NUM,1)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            JND3=JOFF+(KN(NUM,2)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            VAR1=(REAL(IL)**2)*VOL0*QQ(K3+1,K3+1)*GARSI/(FACT*DZ*DZ)
            SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)*FUNKNO(JND1)
            SUNKNO(JND2)=SUNKNO(JND2)+REAL(VAR1)*FUNKNO(JND2)
            SUNKNO(JND3)=SUNKNO(JND3)+REAL(VAR1)*FUNKNO(JND3)
            IF(IL.LE.NLF-3) THEN
               KND1=JND1+L4
               KND2=JND2+L4
               KND3=JND3+L4
               VAR1=(REAL(IL)*REAL(IL+1))*VOL0*QQ(K3+1,K3+1)*GARSI/
     1         (FACT*DZ*DZ)
               SUNKNO(KND1)=SUNKNO(KND1)+REAL(VAR1)*FUNKNO(JND1)
               SUNKNO(KND2)=SUNKNO(KND2)+REAL(VAR1)*FUNKNO(JND2)
               SUNKNO(KND3)=SUNKNO(KND3)+REAL(VAR1)*FUNKNO(JND3)
*
               SUNKNO(JND1)=SUNKNO(JND1)+REAL(VAR1)*FUNKNO(KND1)
               SUNKNO(JND2)=SUNKNO(JND2)+REAL(VAR1)*FUNKNO(KND2)
               SUNKNO(JND3)=SUNKNO(JND3)+REAL(VAR1)*FUNKNO(KND3)
*
               VAR1=(REAL(IL+1)**2)*VOL0*QQ(K3+1,K3+1)*GARSI/
     1         (FACT*DZ*DZ)
               SUNKNO(KND1)=SUNKNO(KND1)+REAL(VAR1)*FUNKNO(KND1)
               SUNKNO(KND2)=SUNKNO(KND2)+REAL(VAR1)*FUNKNO(KND2)
               SUNKNO(KND3)=SUNKNO(KND3)+REAL(VAR1)*FUNKNO(KND3)
            ENDIF
   70       CONTINUE
   71       CONTINUE
   72       CONTINUE
         ENDIF
*        ODD PARITY EQUATION.
         DO 84 K5=0,1 ! TWO LOZENGES PER HEXAGON
         DO 83 K4=0,IELEM-1
         DO 82 K3=0,IELEM-1
         DO 81 K2=1,IELEM+1
         KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         INW1=JOFF+LL4F+ABS(KNW1)
         INX1=JOFF+LL4F+ABS(KNX1)
         INY1=JOFF+LL4F+ABS(KNY1)
         DO 80 K1=1,IELEM+1
         KNW2=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
         KNX2=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
         KNY2=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K1)
         INW2=JOFF+LL4F+ABS(KNW2)
         INX2=JOFF+LL4F+ABS(KNX2)
         INY2=JOFF+LL4F+ABS(KNY2)
         IF((KNW2.NE.0).AND.(KNW1.NE.0)) THEN
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNW2))
            VAR1=(4./3.)*SG*FACT*VOL0*GARS*R(K2,K1)*FUNKNO(INW2)
            SUNKNO(INW1)=SUNKNO(INW1)-REAL(VAR1)
         ENDIF
         IF((KNX2.NE.0).AND.(KNX1.NE.0)) THEN
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNX2))
            VAR1=(4./3.)*SG*FACT*VOL0*GARS*R(K2,K1)*FUNKNO(INX2)
            SUNKNO(INX1)=SUNKNO(INX1)-REAL(VAR1)
         ENDIF
         IF((KNY2.NE.0).AND.(KNY1.NE.0)) THEN
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNY2))
            VAR1=(4./3.)*SG*FACT*VOL0*GARS*R(K2,K1)*FUNKNO(INY2)
            SUNKNO(INY1)=SUNKNO(INY1)-REAL(VAR1)
         ENDIF
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
   83    CONTINUE
   84    CONTINUE
         DO 94 K5=0,2 ! THREE LOZENGES PER HEXAGON
         DO 93 K2=0,IELEM-1
         DO 92 K1=0,IELEM-1
         DO 91 IC=1,2
         IF(IC.EQ.1) IIC=1
         IF(IC.EQ.2) IIC=IELEM+1
         KNZ1=KN(NUM,3+6*NELEH+((2*K5+IC-1)*IELEM+K2)*IELEM+K1+1)
         INZ1=JOFF+LL4F+ABS(KNZ1)
         DO 90 JC=1,2
         IF(JC.EQ.1) JJC=1
         IF(JC.EQ.2) JJC=IELEM+1
         KNZ2=KN(NUM,3+6*NELEH+((2*K5+JC-1)*IELEM+K2)*IELEM+K1+1)
         INZ2=JOFF+LL4F+ABS(KNZ2)
         IF((KNZ1.NE.0).AND.(KNZ2.NE.0)) THEN
            SG=REAL(SIGN(1,KNZ1)*SIGN(1,KNZ2))
            VAR1=SG*FACT*VOL0*GARS*R(IIC,JJC)*FUNKNO(INZ2)
            SUNKNO(INZ1)=SUNKNO(INZ1)-REAL(VAR1)
         ENDIF
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
   93    CONTINUE
   94    CONTINUE
         IF(ITY.EQ.1) GO TO 180
*----
*  BOUNDARY CONDITIONS.
*----
         DO 133 K5=0,1 ! TWO LOZENGES PER HEXAGON
         DO 132 K4=0,IELEM-1
         DO 131 K3=0,IELEM-1
         DO 130 K2=1,IELEM+1,IELEM
         KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         INW1=JOFF+LL4F+ABS(KNW1)
         INX1=JOFF+LL4F+ABS(KNX1)
         INY1=JOFF+LL4F+ABS(KNY1)
         IF(KNW1.NE.0) THEN
            DO 100 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INW2=(IL2/2)*L4+LL4F+ABS(KNW1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,1)*ZMARS*FUNKNO(INW2)
               SUNKNO(INW1)=SUNKNO(INW1)-REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,2)*ZMARS*FUNKNO(INW2)
               SUNKNO(INW1)=SUNKNO(INW1)-REAL(VAR1)
            ENDIF
  100       CONTINUE
         ENDIF
         IF(KNX1.NE.0) THEN
            DO 110 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INX2=(IL2/2)*L4+LL4F+ABS(KNX1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,3)*ZMARS*FUNKNO(INX2)
               SUNKNO(INX1)=SUNKNO(INX1)-REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,4)*ZMARS*FUNKNO(INX2)
               SUNKNO(INX1)=SUNKNO(INX1)-REAL(VAR1)
            ENDIF
  110       CONTINUE
         ENDIF
         IF(KNY1.NE.0) THEN
            DO 120 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INY2=(IL2/2)*L4+LL4F+ABS(KNY1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,5)*ZMARS*FUNKNO(INY2)
               SUNKNO(INY1)=SUNKNO(INY1)-REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,6)*ZMARS*FUNKNO(INY2)
               SUNKNO(INY1)=SUNKNO(INY1)-REAL(VAR1)
            ENDIF
  120       CONTINUE
         ENDIF
  130    CONTINUE
  131    CONTINUE
  132    CONTINUE
  133    CONTINUE
         IF((IPR.EQ.1).OR.(IPR.EQ.2)) GO TO 180
         DO 153 K5=0,2 ! THREE LOZENGES PER HEXAGON
         DO 152 K2=0,IELEM-1
         DO 151 K1=0,IELEM-1
         DO 150 IC=1,2
         KNZ1=KN(NUM,3+6*NELEH+((2*K5+IC-1)*IELEM+K2)*IELEM+K1+1)
         INZ1=JOFF+LL4F+ABS(KNZ1)
         IF(KNZ1.NE.0) THEN
            DO 140 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INZ2=(IL2/2)*L4+LL4F+ABS(KNZ1)
            IF(IC.EQ.1) THEN
               VAR1=0.5*FACT*QFR(NUM,7)*ZMARS*FUNKNO(INZ2)
               SUNKNO(INZ1)=SUNKNO(INZ1)-REAL(VAR1)
            ELSE IF(IC.EQ.2) THEN
               VAR1=0.5*FACT*QFR(NUM,8)*ZMARS*FUNKNO(INZ2)
               SUNKNO(INZ1)=SUNKNO(INZ1)-REAL(VAR1)
            ENDIF
  140       CONTINUE
         ENDIF
  150    CONTINUE
  151    CONTINUE
  152    CONTINUE
  153    CONTINUE
*
         DO 164 K5=0,1 ! TWO LOZENGES PER HEXAGON
         DO 163 K4=0,IELEM-1
         DO 162 K3=0,IELEM-1
         DO 161 K2=1,IELEM+1
         KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         INW1=JOFF+LL4F+ABS(KNW1)
         INX1=JOFF+LL4F+ABS(KNX1)
         INY1=JOFF+LL4F+ABS(KNY1)
         DO 160 K1=0,IELEM-1
         IF(V(K2,K1+1).EQ.0.0) GO TO 160
         IF(K5.EQ.0) THEN
            SSS=(-1.0)**K1
            JND1=JOFF+(((NUM-1)*IELEM+K4)*IELEM+K3)*IELEM+K1+1
            JND2=JOFF+(((KN(NUM,1)-1)*IELEM+K4)*IELEM+K3)*IELEM+K1+1
            JND3=JOFF+(((KN(NUM,2)-1)*IELEM+K4)*IELEM+K3)*IELEM+K1+1
         ELSE
            SSS=1.0
            JND1=JOFF+(((KN(NUM,1)-1)*IELEM+K4)*IELEM+K1)*IELEM+K3+1
            JND2=JOFF+(((KN(NUM,2)-1)*IELEM+K4)*IELEM+K1)*IELEM+K3+1
            JND3=JOFF+(((KN(NUM,3)-1)*IELEM+K4)*IELEM+K1)*IELEM+K3+1
         ENDIF
         IF(KNW1.NE.0) THEN
            SG=REAL(SIGN(1,KNW1))
            VAR1=SG*SSS*REAL(IL)*UUUU*V(K2,K1+1)*FUNKNO(JND1)
            SUNKNO(INW1)=SUNKNO(INW1)+REAL(VAR1)
         ENDIF
         IF(KNX1.NE.0) THEN
            SG=REAL(SIGN(1,KNX1))
            VAR1=SG*SSS*REAL(IL)*UUUU*V(K2,K1+1)*FUNKNO(JND2)
            SUNKNO(INX1)=SUNKNO(INX1)+REAL(VAR1)
         ENDIF
         IF(KNY1.NE.0) THEN
            SG=REAL(SIGN(1,KNY1))
            VAR1=SG*SSS*REAL(IL)*UUUU*V(K2,K1+1)*FUNKNO(JND3)
            SUNKNO(INY1)=SUNKNO(INY1)+REAL(VAR1)
         ENDIF
         IF(IL.LE.NLF-3) THEN
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               VAR1=SG*SSS*REAL(IL+1)*UUUU*V(K2,K1+1)*FUNKNO(JND1+L4)
               SUNKNO(INW1)=SUNKNO(INW1)+REAL(VAR1)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               VAR1=SG*SSS*REAL(IL+1)*UUUU*V(K2,K1+1)*FUNKNO(JND2+L4)
               SUNKNO(INX1)=SUNKNO(INX1)+REAL(VAR1)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               VAR1=SG*SSS*REAL(IL+1)*UUUU*V(K2,K1+1)*FUNKNO(JND3+L4)
               SUNKNO(INY1)=SUNKNO(INY1)+REAL(VAR1)
            ENDIF
         ENDIF
  160    CONTINUE
  161    CONTINUE
  162    CONTINUE
  163    CONTINUE
  164    CONTINUE
         DO 173 K5=0,2 ! THREE LOZENGES PER HEXAGON
         DO 172 K2=0,IELEM-1
         DO 171 K1=0,IELEM-1
         KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
         KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
         INZ1=JOFF+LL4F+ABS(KNZ1)
         INZ2=JOFF+LL4F+ABS(KNZ2)
         DO 170 K3=0,IELEM-1
         IF(K5.EQ.0) THEN
            JND1=JOFF+((((NUM-1)*IELEM)+K3)*IELEM+K2)*IELEM+K1+1
         ELSE
            JND1=JOFF+(((KN(NUM,K5)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
         ENDIF
         IF(KNZ1.NE.0) THEN
            SG=REAL(SIGN(1,KNZ1))
            VAR1=SG*(VOL0/DZ)*REAL(IL)*V(1,K3+1)*FUNKNO(JND1)
            SUNKNO(INZ1)=SUNKNO(INZ1)+REAL(VAR1)
         ENDIF
         IF(KNZ2.NE.0) THEN
            SG=REAL(SIGN(1,KNZ2))
            VAR1=SG*(VOL0/DZ)*REAL(IL)*V(IELEM+1,K3+1)*FUNKNO(JND1)
            SUNKNO(INZ2)=SUNKNO(INZ2)+REAL(VAR1)
         ENDIF
         IF(IL.LE.NLF-3) THEN
            IF(KNZ1.NE.0) THEN
               SG=REAL(SIGN(1,KNZ1))
               VAR1=SG*(VOL0/DZ)*REAL(IL+1)*V(1,K3+1)*FUNKNO(JND1+L4)
               SUNKNO(INZ1)=SUNKNO(INZ1)+REAL(VAR1)
            ENDIF
            IF(KNZ2.NE.0) THEN
               SG=REAL(SIGN(1,KNZ2))
               VAR1=SG*(VOL0/DZ)*REAL(IL+1)*V(IELEM+1,K3+1)*
     1         FUNKNO(JND1+L4)
               SUNKNO(INZ2)=SUNKNO(INZ2)+REAL(VAR1)
            ENDIF
         ENDIF
  170    CONTINUE
  171    CONTINUE
  172    CONTINUE
  173    CONTINUE
      ENDIF
  180 CONTINUE
      IF(MOD(IL,2).EQ.1) THEN
         IOFW=JOFF+LL4F
         IOFX=JOFF+LL4F+LL4W
         IOFY=JOFF+LL4F+LL4W+LL4X
         CALL FLDPWY(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     1   FUNKNO(IOFY+1),SUNKNO(IOFW+1))
         CALL FLDPWX(LL4W,LL4X,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     1   FUNKNO(IOFX+1),SUNKNO(IOFW+1))
         CALL FLDPXW(LL4W,LL4X,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     1   FUNKNO(IOFW+1),SUNKNO(IOFX+1))
         CALL FLDPXY(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     1   FUNKNO(IOFY+1),SUNKNO(IOFX+1))
         CALL FLDPYX(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     1   FUNKNO(IOFX+1),SUNKNO(IOFY+1))
         CALL FLDPYW(LL4W,LL4X,LL4Y,NBLOS,IELEM,CTRAN,IPERT,KN,DIFF,
     1   FUNKNO(IOFW+1),SUNKNO(IOFY+1))
      ENDIF
  200 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIFF)
      RETURN
      END
