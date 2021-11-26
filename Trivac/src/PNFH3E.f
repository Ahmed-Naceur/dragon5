*DECK PNFH3E
      SUBROUTINE PNFH3E (IL,NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,L4,LL4F,
     1 MAT,SIGTI,SIDE,ZZ,FRZ,QFR,IPERT,KN,LC,R,V,SUNKNO,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a one-group SPN flux iteration in hexagonal 3D geometry.
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
* IL      current Legendre order.
* NBMIX   number of mixtures.
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*          =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* LL4F    number of flux components.
* MAT     index-number of the mixture type assigned to each volume.
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
* SUNKNO  sources.
* FUNKNO  initial fluxes.
*
*Parameters: output
* FUNKNO  right-hand-side of the linear system.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IL,NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,L4,LL4F,
     1 MAT(3,NBLOS),IPERT(NBLOS),KN(NBLOS,3+6*(IELEM+2)*IELEM**2),LC
      REAL SIGTI(NBMIX,NAN),SIDE,ZZ(3,NBLOS),FRZ(NBLOS),QFR(NBLOS,8),
     1 R(LC,LC),V(LC,LC-1),SUNKNO(L4*NLF/2),FUNKNO(L4*NLF/2)
*----
*  LOCAL VARIABLES
*----
      REAL QQ(5,5)
      DOUBLE PRECISION FFF,TTTT,UUUU,VOL0,GARSI,FACT,VAR1
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
      JOFF=(IL/2)*L4
      FACT=REAL(2*IL+1)
      IF(MOD(IL,2).EQ.0) THEN
         DO 20 I=1,L4
         FUNKNO(JOFF+I)=SUNKNO(JOFF+I)
   20    CONTINUE
      ENDIF
*----
*  COMPUTE THE SOLUTION AT ORDER IL.
*----
      NELEH=(IELEM+1)*IELEM**2
      TTTT=0.5D0*SQRT(3.D00)*SIDE*SIDE
      NUM=0
      DO 150 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 150
      NUM=NUM+1
      DZ=ZZ(1,IPERT(KEL))
      VOL0=TTTT*DZ*FRZ(KEL)
      UUUU=SIDE*DZ*FRZ(KEL)
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION
         IF(IL.GE.2) THEN
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
            VAR1=SSS*REAL(IL)*UUUU*V(K2,K1+1)
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(VAR1)*FUNKNO(INW1-L4)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               FUNKNO(JND2)=FUNKNO(JND2)-SG*REAL(VAR1)*FUNKNO(INX1-L4)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               FUNKNO(JND3)=FUNKNO(JND3)-SG*REAL(VAR1)*FUNKNO(INY1-L4)
            ENDIF
   30       CONTINUE
   31       CONTINUE
   32       CONTINUE
   33       CONTINUE
   34       CONTINUE
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
               VAR1=SG*(VOL0/DZ)*REAL(IL)*V(1,K3+1)*FUNKNO(INZ1-L4)
               FUNKNO(JND1)=FUNKNO(JND1)-REAL(VAR1)
            ENDIF
            IF(KNZ2.NE.0) THEN
               SG=REAL(SIGN(1,KNZ2))
               VAR1=SG*(VOL0/DZ)*REAL(IL)*V(IELEM+1,K3+1)*
     1         FUNKNO(INZ2-L4)
               FUNKNO(JND1)=FUNKNO(JND1)-REAL(VAR1)
            ENDIF
   40       CONTINUE
   41       CONTINUE
   42       CONTINUE
   43       CONTINUE
         ENDIF
      ELSE
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION
*        OF THE EVEN PARITY EQUATION.
         IBM=MAT(1,IPERT(KEL))
         IF(IBM.EQ.0) GO TO 150
         IF(IELEM.GT.1) THEN
            DO 52 K3=0,IELEM-1
            DO 51 K2=0,IELEM-1
            DO 50 K1=0,IELEM-1
            IF(QQ(K3+1,K3+1).EQ.0.0) GO TO 50
            JND1=JOFF+(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            JND2=JOFF+(KN(NUM,1)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            JND3=JOFF+(KN(NUM,2)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            IF(IL.GE.3) THEN
               GARSI=SIGTI(IBM,MIN(IL-1,NAN))
               KND1=JND1-L4
               KND2=JND2-L4
               KND3=JND3-L4
               VAR1=(REAL(IL-1)*REAL(IL-2))*VOL0*QQ(K3+1,K3+1)*GARSI
     1         /(REAL(2*IL-3)*DZ*DZ)
               FUNKNO(JND1)=FUNKNO(JND1)-REAL(VAR1)*FUNKNO(KND1)
               FUNKNO(JND2)=FUNKNO(JND2)-REAL(VAR1)*FUNKNO(KND2)
               FUNKNO(JND3)=FUNKNO(JND3)-REAL(VAR1)*FUNKNO(KND3)
            ENDIF
            IF(IL.LE.NLF-3) THEN
               GARSI=SIGTI(IBM,MIN(IL+1,NAN))
               KND1=JND1+L4
               KND2=JND2+L4
               KND3=JND3+L4
               VAR1=(REAL(IL)*REAL(IL+1))*VOL0*QQ(K3+1,K3+1)*GARSI
     1         /(FACT*DZ*DZ)
               FUNKNO(JND1)=FUNKNO(JND1)-REAL(VAR1)*FUNKNO(KND1)
               FUNKNO(JND2)=FUNKNO(JND2)-REAL(VAR1)*FUNKNO(KND2)
               FUNKNO(JND3)=FUNKNO(JND3)-REAL(VAR1)*FUNKNO(KND3)
            ENDIF
   50       CONTINUE
   51       CONTINUE
   52       CONTINUE
         ENDIF
*
*        ODD PARITY EQUATION
         DO 93 K5=0,1 ! TWO LOZENGES PER HEXAGON
         DO 92 K4=0,IELEM-1
         DO 91 K3=0,IELEM-1
         DO 90 K2=1,IELEM+1
         KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
         INW1=JOFF+LL4F+ABS(KNW1)
         INX1=JOFF+LL4F+ABS(KNX1)
         INY1=JOFF+LL4F+ABS(KNY1)
         IF(KNW1.NE.0) THEN
            DO 60 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 60
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INW2=(IL2/2)*L4+LL4F+ABS(KNW1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,1)*ZMARS*FUNKNO(INW2)
               FUNKNO(INW1)=FUNKNO(INW1)+REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,2)*ZMARS*FUNKNO(INW2)
               FUNKNO(INW1)=FUNKNO(INW1)+REAL(VAR1)
            ENDIF
   60       CONTINUE
         ENDIF
         IF(KNX1.NE.0) THEN
            DO 70 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 70
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INX2=(IL2/2)*L4+LL4F+ABS(KNX1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,3)*ZMARS*FUNKNO(INX2)
               FUNKNO(INX1)=FUNKNO(INX1)+REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,4)*ZMARS*FUNKNO(INX2)
               FUNKNO(INX1)=FUNKNO(INX1)+REAL(VAR1)
            ENDIF
   70       CONTINUE
         ENDIF
         IF(KNY1.NE.0) THEN
            DO 80 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 80
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INY2=(IL2/2)*L4+LL4F+ABS(KNY1)
            IF((K2.EQ.1).AND.(K5.EQ.0)) THEN
               VAR1=0.5*FACT*QFR(NUM,5)*ZMARS*FUNKNO(INY2)
               FUNKNO(INY1)=FUNKNO(INY1)+REAL(VAR1)
            ELSE IF((K2.EQ.IELEM+1).AND.(K5.EQ.1)) THEN
               VAR1=0.5*FACT*QFR(NUM,6)*ZMARS*FUNKNO(INY2)
               FUNKNO(INY1)=FUNKNO(INY1)+REAL(VAR1)
            ENDIF
   80       CONTINUE
         ENDIF
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
   93    CONTINUE
         DO 122 K5=0,2 ! THREE LOZENGES PER HEXAGON
         DO 121 K2=0,IELEM-1
         DO 120 K1=0,IELEM-1
         KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
         KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
         INZ1=JOFF+LL4F+ABS(KNZ1)
         INZ2=JOFF+LL4F+ABS(KNZ2)
         IF((QFR(NUM,7).NE.0.0).AND.(KNZ1.NE.0)) THEN
*           ZINF SIDE.
            DO 100 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 100
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+LL4F+ABS(KNZ1)
            VAR1=0.5*FACT*QFR(NUM,7)*ZMARS*FUNKNO(INDL)
            FUNKNO(INZ1)=FUNKNO(INZ1)+REAL(VAR1)
  100       CONTINUE
         ENDIF
         IF((QFR(NUM,8).NE.0.0).AND.(KNZ2.NE.0)) THEN
*           ZSUP SIDE.
            DO 110 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 110
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+LL4F+ABS(KNZ2)
            VAR1=0.5*FACT*QFR(NUM,8)*ZMARS*FUNKNO(INDL)
            FUNKNO(INZ2)=FUNKNO(INZ2)+REAL(VAR1)
  110       CONTINUE
         ENDIF
  120    CONTINUE
  121    CONTINUE
  122    CONTINUE
*
         IF(IL.LE.NLF-3) THEN
            DO 134 K5=0,1 ! TWO LOZENGES PER HEXAGON
            DO 133 K4=0,IELEM-1
            DO 132 K3=0,IELEM-1
            DO 131 K2=1,IELEM+1
            KNW1=KN(NUM,3+K5*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
            KNX1=KN(NUM,3+(K5+2)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
            KNY1=KN(NUM,3+(K5+4)*NELEH+(K4*IELEM+K3)*(IELEM+1)+K2)
            INW1=JOFF+LL4F+ABS(KNW1)
            INX1=JOFF+LL4F+ABS(KNX1)
            INY1=JOFF+LL4F+ABS(KNY1)
            DO 130 K1=0,IELEM-1
            IF(V(K2,K1+1).EQ.0.0) GO TO 130
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
            VAR1=SSS*REAL(IL+1)*UUUU*V(K2,K1+1)
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               FUNKNO(INW1)=FUNKNO(INW1)-SG*REAL(VAR1)*FUNKNO(JND1+L4)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               FUNKNO(INX1)=FUNKNO(INX1)-SG*REAL(VAR1)*FUNKNO(JND2+L4)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               FUNKNO(INY1)=FUNKNO(INY1)-SG*REAL(VAR1)*FUNKNO(JND3+L4)
            ENDIF
  130       CONTINUE
  131       CONTINUE
  132       CONTINUE
  133       CONTINUE
  134       CONTINUE
            DO 143 K5=0,2 ! THREE LOZENGES PER HEXAGON
            DO 142 K2=0,IELEM-1
            DO 141 K1=0,IELEM-1
            KNZ1=KN(NUM,3+6*NELEH+2*K5*IELEM**2+K2*IELEM+K1+1)
            KNZ2=KN(NUM,3+6*NELEH+(2*K5+1)*IELEM**2+K2*IELEM+K1+1)
            INZ1=JOFF+LL4F+ABS(KNZ1)
            INZ2=JOFF+LL4F+ABS(KNZ2)
            DO 140 K3=0,IELEM-1
            IF(K5.EQ.0) THEN
               JND1=JOFF+((((NUM-1)*IELEM)+K3)*IELEM+K2)*IELEM+K1+1
            ELSE
               JND1=JOFF+(((KN(NUM,K5)-1)*IELEM+K3)*IELEM+K2)*IELEM+K1+1
            ENDIF
            IF(KNZ1.NE.0) THEN
               SG=REAL(SIGN(1,KNZ1))
               VAR1=SG*(VOL0/DZ)*REAL(IL+1)*V(1,K3+1)*FUNKNO(JND1+L4)
               FUNKNO(INZ1)=FUNKNO(INZ1)-REAL(VAR1)
            ENDIF
            IF(KNZ2.NE.0) THEN
               SG=REAL(SIGN(1,KNZ2))
               VAR1=SG*(VOL0/DZ)*REAL(IL+1)*V(IELEM+1,K3+1)*
     1         FUNKNO(JND1+L4)
               FUNKNO(INZ2)=FUNKNO(INZ2)-REAL(VAR1)
            ENDIF
  140       CONTINUE
  141       CONTINUE
  142       CONTINUE
  143       CONTINUE
         ENDIF
      ENDIF
  150 CONTINUE
      RETURN
      END
