*DECK PNFL3E
      SUBROUTINE PNFL3E (IL,NREG,IELEM,ICOL,XX,YY,ZZ,MAT,VOL,NBMIX,NLF,
     1 NVD,NAN,SIGTI,L4,KN,QFR,LC,R,V,SUNKNO,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a one-group SPN flux iteration in Cartesian 3D geometry.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IL      current Legendre order.
* NREG    total number of regions.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*          =2 (parabolic); =3 (cubic); =4 (quartic).
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
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only.
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
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
      INTEGER IL,NREG,IELEM,ICOL,MAT(NREG),NBMIX,NLF,NVD,NAN,L4,
     1 KN(NREG*(1+6*IELEM**2)),LC
      REAL XX(NREG),YY(NREG),ZZ(NREG),VOL(NREG),SIGTI(NBMIX,NAN),
     1 QFR(6*NREG),R(LC,LC),V(LC,LC-1),SUNKNO(L4*NLF/2),FUNKNO(L4*NLF/2)
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
      FACT=REAL(2*IL+1)
      IF(MOD(IL,2).EQ.0) THEN
         DO 20 I=1,L4
         FUNKNO((IL/2)*L4+I)=SUNKNO((IL/2)*L4+I)
   20    CONTINUE
      ENDIF
*----
*  COMPUTE THE SOLUTION AT ORDER IL.
*----
      NUM1=0
      NUM2=0
      DO 150 K=1,NREG
      IBM=MAT(K)
      IF(IBM.EQ.0) GO TO 150
      VOL0=VOL(K)
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION
         IF(IL.GE.2) THEN
            DO 32 K3=0,IELEM-1
            DO 31 K2=0,IELEM-1
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
            DO 30 K1=0,IELEM-1
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            IF(KN1.NE.0) THEN
               SG=REAL(SIGN(1,KN1))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(1,K1+1)*
     1         FUNKNO(IND1)/XX(K)
            ENDIF
            IF(KN2.NE.0) THEN
               SG=REAL(SIGN(1,KN2))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(IND2)/XX(K)
            ENDIF
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            IF(KN3.NE.0) THEN
               SG=REAL(SIGN(1,KN3))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(1,K1+1)*
     1         FUNKNO(IND3)/YY(K)
            ENDIF
            IF(KN4.NE.0) THEN
               SG=REAL(SIGN(1,KN4))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(IND4)/YY(K)
            ENDIF
            JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            IF(KN5.NE.0) THEN
               SG=REAL(SIGN(1,KN5))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(1,K1+1)*
     1         FUNKNO(IND5)/ZZ(K)
            ENDIF
            IF(KN6.NE.0) THEN
               SG=REAL(SIGN(1,KN6))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(IND6)/ZZ(K)
            ENDIF
   30       CONTINUE
   31       CONTINUE
   32       CONTINUE
         ENDIF
      ELSE
         DO 145 K3=0,IELEM-1
         DO 140 K2=0,IELEM-1
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION
*        OF THE EVEN PARITY EQUATION.
         IF((IL.GE.3).AND.(IELEM.GT.1)) THEN
            GARSI=SIGTI(IBM,MIN(IL-1,NAN))
            DO 40 K1=0,IELEM-1
            IF(QQ(K1+1,K1+1).EQ.0.0) GO TO 40
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            KND1=((IL-2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL-1)*REAL(IL-2))*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(REAL(2*IL-3)*XX(K)*XX(K))
*
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            KND1=((IL-2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL-1)*REAL(IL-2))*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(REAL(2*IL-3)*YY(K)*YY(K))
*
            JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            KND1=((IL-2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL-1)*REAL(IL-2))*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(REAL(2*IL-3)*ZZ(K)*ZZ(K))
   40       CONTINUE
         ENDIF
         IF((IL.LE.NLF-3).AND.(IELEM.GT.1)) THEN
            GARSI=SIGTI(IBM,MIN(IL+1,NAN))
            DO 50 K1=0,IELEM-1
            IF(QQ(K1+1,K1+1).EQ.0.0) GO TO 50
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL)*REAL(IL+1))*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
*
            JND1=(IL/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL)*REAL(IL+1))*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
*
            JND1=(IL/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL)*REAL(IL+1))*VOL0*
     1      QQ(K1+1,K1+1)*GARSI*FUNKNO(KND1)/(FACT*ZZ(K)*ZZ(K))
   50       CONTINUE
         ENDIF
*
*        ODD PARITY EQUATION
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
            DO 60 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 60
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN1)
            FUNKNO(IND1)=FUNKNO(IND1)+0.5*FACT*QFR(NUM2+1)*ZMARS*
     1      FUNKNO(INDL)
   60       CONTINUE
         ENDIF
         IF((QFR(NUM2+2).NE.0.0).AND.(KN2.NE.0)) THEN
*           XSUP SIDE.
            DO 70 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 70
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN2)
            FUNKNO(IND2)=FUNKNO(IND2)+0.5*FACT*QFR(NUM2+2)*ZMARS*
     1      FUNKNO(INDL)
   70       CONTINUE
         ENDIF
         IF((QFR(NUM2+3).NE.0.0).AND.(KN3.NE.0)) THEN
*           YINF SIDE.
            DO 80 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 80
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN3)
            FUNKNO(IND3)=FUNKNO(IND3)+0.5*FACT*QFR(NUM2+3)*ZMARS*
     1      FUNKNO(INDL)
   80       CONTINUE
         ENDIF
         IF((QFR(NUM2+4).NE.0.0).AND.(KN4.NE.0)) THEN
*           YSUP SIDE.
            DO 90 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 90
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN4)
            FUNKNO(IND4)=FUNKNO(IND4)+0.5*FACT*QFR(NUM2+4)*ZMARS*
     1      FUNKNO(INDL)
   90       CONTINUE
         ENDIF
         IF((QFR(NUM2+5).NE.0.0).AND.(KN5.NE.0)) THEN
*           ZINF SIDE.
            DO 100 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 100
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN5)
            FUNKNO(IND5)=FUNKNO(IND5)+0.5*FACT*QFR(NUM2+5)*ZMARS*
     1      FUNKNO(INDL)
  100       CONTINUE
         ENDIF
         IF((QFR(NUM2+6).NE.0.0).AND.(KN6.NE.0)) THEN
*           ZSUP SIDE.
            DO 110 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 110
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INDL=(IL2/2)*L4+ABS(KN6)
            FUNKNO(IND6)=FUNKNO(IND6)+0.5*FACT*QFR(NUM2+6)*ZMARS*
     1      FUNKNO(INDL)
  110       CONTINUE
         ENDIF
         IF(IL.LE.NLF-3) THEN
            DO 130 K1=0,IELEM-1
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K2)*IELEM+K1
            IF(KN1.NE.0) THEN
               SG=REAL(SIGN(1,KN1))
               FUNKNO(IND1)=FUNKNO(IND1)-SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1         FUNKNO(JND1)/XX(K)
            ENDIF
            IF(KN2.NE.0) THEN
               SG=REAL(SIGN(1,KN2))
               FUNKNO(IND2)=FUNKNO(IND2)-SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(JND1)/XX(K)
            ENDIF
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K3*IELEM+K1)*IELEM+K2
            IF(KN3.NE.0) THEN
               SG=REAL(SIGN(1,KN3))
               FUNKNO(IND3)=FUNKNO(IND3)-SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1         FUNKNO(JND1)/YY(K)
            ENDIF
            IF(KN4.NE.0) THEN
               SG=REAL(SIGN(1,KN4))
               FUNKNO(IND4)=FUNKNO(IND4)-SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(JND1)/YY(K)
            ENDIF
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(K1*IELEM+K3)*IELEM+K2
            IF(KN5.NE.0) THEN
               SG=REAL(SIGN(1,KN5))
               FUNKNO(IND5)=FUNKNO(IND5)-SG*REAL(IL+1)*VOL0*V(1,K1+1)*
     1         FUNKNO(JND1)/ZZ(K)
            ENDIF
            IF(KN6.NE.0) THEN
               SG=REAL(SIGN(1,KN6))
               FUNKNO(IND6)=FUNKNO(IND6)-SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,K1+1)*FUNKNO(JND1)/ZZ(K)
            ENDIF
  130       CONTINUE
         ENDIF
  140    CONTINUE
  145    CONTINUE
      ENDIF
      NUM1=NUM1+1+6*IELEM**2
      NUM2=NUM2+6
  150 CONTINUE
      RETURN
      END
