*DECK PNSZ3D
      SUBROUTINE PNSZ3D(ITY,IPR,NREG,IELEM,ICOL,XX,YY,ZZ,MAT,VOL,NBMIX,
     1 NLF,NVD,NAN,SIGT,SIGTI,L4,KN,QFR,LC,R,V,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source calculation for a SPN approximation in TRIVAC, including
* neighbour Legendre and out-of-group contributions.
* Raviart-Thomas method in Cartesian geometry.
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
* ITY     type of assembly:
*         =0: leakage-removal matrix assembly; =1: cross section matrix
*         assembly.
* IPR     type of assembly:
*         =0: contains system matrices;
*         =1: contains derivative of these matrices;
*         =2: contains first variation of these matrices;
*         =3: contains addition of first vatiation to unperturbed
*         system matrices.
* NREG    total number of regions.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* NBMIX   number of mixtures.
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* SIGT    macroscopic cross sections ordered by mixture.
*         SIGT(:,NAN) generally contains the total cross section only.
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only.
* L4      order of the profiled system matrices.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* LC      order of the unit matrices.
* R       unit Cartesian mass matrix.
* V       unit nodal coupling matrix.
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
      INTEGER ITY,IPR,NREG,IELEM,ICOL,MAT(NREG),NBMIX,NLF,NVD,NAN,L4,
     1 KN(NREG*(1+6*IELEM**2)),LC
      REAL XX(NREG),YY(NREG),ZZ(NREG),VOL(NREG),SIGT(NBMIX,NAN),
     1 SIGTI(NBMIX,NAN),QFR(6*NREG),R(LC,LC),V(LC,LC-1),
     2 SUNKNO(L4*NLF/2),FUNKNO(L4*NLF/2)
*----
*  LOCAL VARIABLES
*----
      REAL QQ(5,5)
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
      DO 12 I0=1,IELEM
      DO 11 J0=1,IELEM
      QQ(I0,J0)=0.0
      DO 10 K0=2,IELEM
      QQ(I0,J0)=QQ(I0,J0)+V(K0,I0)*V(K0,J0)/R(K0,K0)
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
      DO 200 IL=0,NLF-1
      IF((ITY.EQ.1).AND.(IL.GE.NAN)) GO TO 200
      FACT=REAL(2*IL+1)
*----
*  COMPUTE THE SOURCE AT ORDER IL.
*----
      NUM1=0
      NUM2=0
      DO 190 K=1,NREG
      IBM=MAT(K)
      IF(IBM.EQ.0) GO TO 190
      VOL0=VOL(K)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION
         DO 55 K3=0,IELEM-1
         DO 50 K2=0,IELEM-1
         DO 20 K1=0,IELEM-1
         JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
         SUNKNO(JND1)=SUNKNO(JND1)+FACT*VOL0*GARS*FUNKNO(JND1)
   20    CONTINUE
         IF(ITY.EQ.1) GO TO 50
         IF((IPR.EQ.1).OR.(IPR.EQ.2)) GO TO 50
*
         KN1=KN(NUM1+2+K3*IELEM+K2)
         KN2=KN(NUM1+2+IELEM**2+K3*IELEM+K2)
         KN3=KN(NUM1+2+2*IELEM**2+K3*IELEM+K2)
         KN4=KN(NUM1+2+3*IELEM**2+K3*IELEM+K2)
         KN5=KN(NUM1+2+4*IELEM**2+K3*IELEM+K2)
         KN6=KN(NUM1+2+5*IELEM**2+K3*IELEM+K2)
         IND1=(IL/2)*L4+ABS(KN1)
         IND2=(IL/2)*L4+ABS(KN2)
         IND3=(IL/2)*L4+ABS(KN3)
         IND4=(IL/2)*L4+ABS(KN4)
         IND5=(IL/2)*L4+ABS(KN5)
         IND6=(IL/2)*L4+ABS(KN6)
         DO 30 K1=0,IELEM-1
         JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
         IF(KN1.NE.0) THEN
            SG=REAL(SIGN(1,KN1))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1      FUNKNO(IND1)/XX(K)
         ENDIF
         IF(KN2.NE.0) THEN
            SG=REAL(SIGN(1,KN2))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(IELEM+1,K1+1)
     1      *FUNKNO(IND2)/XX(K)
         ENDIF
         JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
         IF(KN3.NE.0) THEN
            SG=REAL(SIGN(1,KN3))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1      FUNKNO(IND3)/YY(K)
         ENDIF
         IF(KN4.NE.0) THEN
            SG=REAL(SIGN(1,KN4))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(IELEM+1,K1+1)
     1      *FUNKNO(IND4)/YY(K)
         ENDIF
         JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
         IF(KN5.NE.0) THEN
            SG=REAL(SIGN(1,KN5))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1      FUNKNO(IND5)/ZZ(K)
         ENDIF
         IF(KN6.NE.0) THEN
            SG=REAL(SIGN(1,KN6))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(IELEM+1,K1+1)
     1      *FUNKNO(IND6)/ZZ(K)
         ENDIF
   30    CONTINUE
         IF(IL.GE.2) THEN
            KN1=KN(NUM1+2+K3*IELEM+K2)
            KN2=KN(NUM1+2+IELEM**2+K3*IELEM+K2)
            KN3=KN(NUM1+2+2*IELEM**2+K3*IELEM+K2)
            KN4=KN(NUM1+2+3*IELEM**2+K3*IELEM+K2)
            KN5=KN(NUM1+2+4*IELEM**2+K3*IELEM+K2)
            KN6=KN(NUM1+2+5*IELEM**2+K3*IELEM+K2)
            IND1=((IL-2)/2)*L4+ABS(KN1)
            IND2=((IL-2)/2)*L4+ABS(KN2)
            IND3=((IL-2)/2)*L4+ABS(KN3)
            IND4=((IL-2)/2)*L4+ABS(KN4)
            IND5=((IL-2)/2)*L4+ABS(KN5)
            IND6=((IL-2)/2)*L4+ABS(KN6)
            DO 40 K1=0,IELEM-1
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            IF(KN1.NE.0) THEN
               SG=REAL(SIGN(1,KN1))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(1,K1+1)*
     1         FUNKNO(IND1)/XX(K)
            ENDIF
            IF(KN2.NE.0) THEN
               SG=REAL(SIGN(1,KN2))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(IND2)/XX(K)
            ENDIF
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            IF(KN3.NE.0) THEN
               SG=REAL(SIGN(1,KN3))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(1,K1+1)*
     1         FUNKNO(IND3)/YY(K)
            ENDIF
            IF(KN4.NE.0) THEN
               SG=REAL(SIGN(1,KN4))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(IND4)/YY(K)
            ENDIF
            JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            IF(KN5.NE.0) THEN
               SG=REAL(SIGN(1,KN5))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(1,K1+1)*
     1         FUNKNO(IND5)/ZZ(K)
            ENDIF
            IF(KN6.NE.0) THEN
               SG=REAL(SIGN(1,KN6))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(IND6)/ZZ(K)
            ENDIF
   40       CONTINUE
         ENDIF
   50    CONTINUE
   55    CONTINUE
      ELSE IF(MOD(IL,2).EQ.1) THEN
         GARSI=SIGTI(IBM,MIN(IL+1,NAN))
         DO 185 K3=0,IELEM-1
         DO 180 K2=0,IELEM-1
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION OF
*        THE EVEN PARITY EQUATION.
         IF(IELEM.GT.1) THEN
            DO 60 K1=0,IELEM-1
            IF(QQ(K1+1,K1+1).EQ.0.0) GO TO 60
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)**2)*VOL0*QQ(K1+1,K1+1)*
     1      GARSI*FUNKNO(JND1)/(FACT*XX(K)*XX(K))
*
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)**2)*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(JND1)/(FACT*YY(K)*YY(K))
*
            JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)**2)*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(JND1)/(FACT*ZZ(K)*ZZ(K))
            IF(IL.LE.NLF-3) THEN
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
               KND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
               JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL+1)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(JND1)/(FACT*XX(K)*XX(K))
*
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
               KND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
               JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL+1)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(JND1)/(FACT*YY(K)*YY(K))
*
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
               KND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*ZZ(K)*ZZ(K))
               JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*ZZ(K)*ZZ(K))
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL+1)*REAL(IL+1))*VOL0*
     1         QQ(K1+1,K1+1)*GARSI*FUNKNO(JND1)/(FACT*ZZ(K)*ZZ(K))
            ENDIF
   60       CONTINUE
         ENDIF
*        ODD PARITY EQUATION.
         DO 75 IC=1,2
         IF(IC.EQ.1) IIC=1
         IF(IC.EQ.2) IIC=IELEM+1
         KN1=KN(NUM1+2+(IC-1)*IELEM**2+K3*IELEM+K2)
         IND1=(IL/2)*L4+ABS(KN1)
         S1=REAL(SIGN(1,KN1))
         DO 70 JC=1,2
         KN2=KN(NUM1+2+(JC-1)*IELEM**2+K3*IELEM+K2)
         IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
            JJC=1
            IF(JC.EQ.2) JJC=IELEM+1
            IND2=(IL/2)*L4+ABS(KN2)
            S2=REAL(SIGN(1,KN2))
            SUNKNO(IND1)=SUNKNO(IND1)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS*
     1      FUNKNO(IND2)
         ENDIF
   70    CONTINUE
   75    CONTINUE
         DO 85 IC=3,4
         IF(IC.EQ.3) IIC=1
         IF(IC.EQ.4) IIC=IELEM+1
         KN1=KN(NUM1+2+(IC-1)*IELEM**2+K3*IELEM+K2)
         IND1=(IL/2)*L4+ABS(KN1)
         S1=REAL(SIGN(1,KN1))
         DO 80 JC=3,4
         KN2=KN(NUM1+2+(JC-1)*IELEM**2+K3*IELEM+K2)
         IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
            JJC=1
            IF(JC.EQ.4) JJC=IELEM+1
            IND2=(IL/2)*L4+ABS(KN2)
            S2=REAL(SIGN(1,KN2))
            SUNKNO(IND1)=SUNKNO(IND1)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS*
     1      FUNKNO(IND2)
         ENDIF
   80    CONTINUE
   85    CONTINUE
         DO 95 IC=5,6
         IF(IC.EQ.5) IIC=1
         IF(IC.EQ.6) IIC=IELEM+1
         KN1=KN(NUM1+2+(IC-1)*IELEM**2+K3*IELEM+K2)
         IND1=(IL/2)*L4+ABS(KN1)
         S1=REAL(SIGN(1,KN1))
         DO 90 JC=5,6
         KN2=KN(NUM1+2+(JC-1)*IELEM**2+K3*IELEM+K2)
         IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
            JJC=1
            IF(JC.EQ.6) JJC=IELEM+1
            IND2=(IL/2)*L4+ABS(KN2)
            S2=REAL(SIGN(1,KN2))
            SUNKNO(IND1)=SUNKNO(IND1)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS*
     1      FUNKNO(IND2)
         ENDIF
   90    CONTINUE
   95    CONTINUE
         IF(ITY.EQ.1) GO TO 180
*
         KN1=KN(NUM1+2+K3*IELEM+K2)
         KN2=KN(NUM1+2+IELEM**2+K3*IELEM+K2)
         KN3=KN(NUM1+2+2*IELEM**2+K3*IELEM+K2)
         KN4=KN(NUM1+2+3*IELEM**2+K3*IELEM+K2)
         KN5=KN(NUM1+2+4*IELEM**2+K3*IELEM+K2)
         KN6=KN(NUM1+2+5*IELEM**2+K3*IELEM+K2)
         IND1=(IL/2)*L4+ABS(KN1)
         IND2=(IL/2)*L4+ABS(KN2)
         IND3=(IL/2)*L4+ABS(KN3)
         IND4=(IL/2)*L4+ABS(KN4)
         IND5=(IL/2)*L4+ABS(KN5)
         IND6=(IL/2)*L4+ABS(KN6)
         IF((QFR(NUM2+1).NE.0.0).AND.(KN1.NE.0)) THEN
*           XINF SIDE.
            DO 100 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN1)
            SUNKNO(IND1)=SUNKNO(IND1)-0.5*FACT*QFR(NUM2+1)*ZMARS*
     1      FUNKNO(INDL)
  100       CONTINUE
         ENDIF
         IF((QFR(NUM2+2).NE.0.0).AND.(KN2.NE.0)) THEN
*           XSUP SIDE.
            DO 110 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN2)
            SUNKNO(IND2)=SUNKNO(IND2)-0.5*FACT*QFR(NUM2+2)*ZMARS*
     1      FUNKNO(INDL)
  110       CONTINUE
         ENDIF
         IF((QFR(NUM2+3).NE.0.0).AND.(KN3.NE.0)) THEN
*           YINF SIDE.
            DO 120 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN3)
            SUNKNO(IND3)=SUNKNO(IND3)-0.5*FACT*QFR(NUM2+3)*ZMARS*
     1      FUNKNO(INDL)
  120       CONTINUE
         ENDIF
         IF((QFR(NUM2+4).NE.0.0).AND.(KN4.NE.0)) THEN
*           YSUP SIDE.
            DO 130 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN4)
            SUNKNO(IND4)=SUNKNO(IND4)-0.5*FACT*QFR(NUM2+4)*ZMARS*
     1      FUNKNO(INDL)
  130       CONTINUE
         ENDIF
         IF((QFR(NUM2+5).NE.0.0).AND.(KN5.NE.0)) THEN
*           ZINF SIDE.
            DO 140 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN5)
            SUNKNO(IND5)=SUNKNO(IND5)-0.5*FACT*QFR(NUM2+5)*ZMARS*
     1      FUNKNO(INDL)
  140       CONTINUE
         ENDIF
         IF((QFR(NUM2+6).NE.0.0).AND.(KN6.NE.0)) THEN
*           ZSUP SIDE.
            DO 150 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN6)
            SUNKNO(IND6)=SUNKNO(IND6)-0.5*FACT*QFR(NUM2+6)*ZMARS*
     1      FUNKNO(INDL)
  150       CONTINUE
         ENDIF
*
         IF((IPR.EQ.1).OR.(IPR.EQ.2)) GO TO 180
         DO 160 K1=0,IELEM-1
         JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
         IF(KN1.NE.0) THEN
            SG=REAL(SIGN(1,KN1))
            SUNKNO(IND1)=SUNKNO(IND1)+SG*REAL(IL)*VOL0*V(1,K1+1)*
     1      FUNKNO(JND1)/XX(K)
         ENDIF
         IF(KN2.NE.0) THEN
            SG=REAL(SIGN(1,KN2))
            SUNKNO(IND2)=SUNKNO(IND2)+SG*REAL(IL)*VOL0*V(IELEM+1,K1+1)*
     1      FUNKNO(JND1)/XX(K)
         ENDIF
         JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
         IF(KN3.NE.0) THEN
            SG=REAL(SIGN(1,KN3))
            SUNKNO(IND3)=SUNKNO(IND3)+SG*REAL(IL)*VOL0*V(1,K1+1)*
     1      FUNKNO(JND1)/YY(K)
         ENDIF
         IF(KN4.NE.0) THEN
            SG=REAL(SIGN(1,KN4))
            SUNKNO(IND4)=SUNKNO(IND4)+SG*REAL(IL)*VOL0*V(IELEM+1,K1+1)*
     1      FUNKNO(JND1)/YY(K)
         ENDIF
         JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
         IF(KN5.NE.0) THEN
            SG=REAL(SIGN(1,KN5))
            SUNKNO(IND5)=SUNKNO(IND5)+SG*REAL(IL)*VOL0*V(1,K1+1)*
     1      FUNKNO(JND1)/ZZ(K)
         ENDIF
         IF(KN6.NE.0) THEN
            SG=REAL(SIGN(1,KN6))
            SUNKNO(IND6)=SUNKNO(IND6)+SG*REAL(IL)*VOL0*V(IELEM+1,K1+1)*
     1      FUNKNO(JND1)/ZZ(K)
         ENDIF
  160    CONTINUE
         IF(IL.LE.NLF-3) THEN
            DO 170 K1=0,IELEM-1
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            IF(KN1.NE.0) THEN
               SG=REAL(SIGN(1,KN1))
               SUNKNO(IND1)=SUNKNO(IND1)+SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1         FUNKNO(JND1)/XX(K)
            ENDIF
            IF(KN2.NE.0) THEN
               SG=REAL(SIGN(1,KN2))
               SUNKNO(IND2)=SUNKNO(IND2)+SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(JND1)/XX(K)
            ENDIF
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            IF(KN3.NE.0) THEN
               SG=REAL(SIGN(1,KN3))
               SUNKNO(IND3)=SUNKNO(IND3)+SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1         FUNKNO(JND1)/YY(K)
            ENDIF
            IF(KN4.NE.0) THEN
               SG=REAL(SIGN(1,KN4))
               SUNKNO(IND4)=SUNKNO(IND4)+SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(JND1)/YY(K)
            ENDIF
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            IF(KN5.NE.0) THEN
               SG=REAL(SIGN(1,KN5))
               SUNKNO(IND5)=SUNKNO(IND5)+SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1         FUNKNO(JND1)/ZZ(K)
            ENDIF
            IF(KN6.NE.0) THEN
               SG=REAL(SIGN(1,KN6))
               SUNKNO(IND6)=SUNKNO(IND6)+SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(JND1)/ZZ(K)
            ENDIF
  170       CONTINUE
         ENDIF
  180    CONTINUE
  185    CONTINUE
      ENDIF
      NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  190 CONTINUE
  200 CONTINUE
      RETURN
      END
