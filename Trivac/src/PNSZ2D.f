*DECK PNSZ2D
      SUBROUTINE PNSZ2D(ITY,NREG,IELEM,ICOL,XX,YY,MAT,VOL,NBMIX,NLF,
     1 NVD,NAN,SIGT,SIGTI,L4,KN,QFR,LC,R,V,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source calculation for a SPN approximation in BIVAC, including
* neighbour Legendre and out-of-group contributions.
* Raviart-Thomas method in Cartesian geometry.
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
* ITY     type of assembly:
*         =0: leakage-removal matrix assembly; =1: cross section matrix
*         assembly.
* NREG    total number of regions.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
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
      INTEGER ITY,NREG,IELEM,ICOL,MAT(NREG),NBMIX,NLF,NVD,NAN,L4,
     1 KN(5*NREG),LC
      REAL XX(NREG),YY(NREG),VOL(NREG),SIGT(NBMIX,NAN),SIGTI(NBMIX,NAN),
     1 QFR(4*NREG),R(LC,LC),V(LC,LC-1),SUNKNO(L4*NLF/2),FUNKNO(L4*NLF/2)
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
      DO 170 IL=0,NLF-1
      IF((ITY.EQ.1).AND.(IL.GE.NAN)) GO TO 170
      FACT=REAL(2*IL+1)
*----
*  COMPUTE THE SOURCE AT ORDER IL.
*----
      NUM1=0
      NUM2=0
      DO 160 K=1,NREG
      IBM=MAT(K)
      IF(IBM.EQ.0) GO TO 160
      VOL0=VOL(K)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
         DO 50 I0=1,IELEM
         DO 20 J0=1,IELEM
         JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
         SUNKNO(JND1)=SUNKNO(JND1)+FACT*VOL0*GARS*FUNKNO(JND1)
   20    CONTINUE
         IF(ITY.EQ.1) GO TO 50
*
         IND1=(IL/2)*L4+ABS(KN(NUM1+2))+I0-1
         IND2=(IL/2)*L4+ABS(KN(NUM1+3))+I0-1
         IND3=(IL/2)*L4+ABS(KN(NUM1+4))+I0-1
         IND4=(IL/2)*L4+ABS(KN(NUM1+5))+I0-1
         DO 30 J0=1,IELEM
         JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
         IF(KN(NUM1+2).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+2)))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(1,J0)*
     1      FUNKNO(IND1)/XX(K)
         ENDIF
         IF(KN(NUM1+3).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+3)))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(IELEM+1,J0)*
     1      FUNKNO(IND2)/XX(K)
         ENDIF
         JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
         IF(KN(NUM1+4).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+4)))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(1,J0)*
     1      FUNKNO(IND3)/YY(K)
         ENDIF
         IF(KN(NUM1+5).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+5)))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL+1)*VOL0*V(IELEM+1,J0)*
     1      FUNKNO(IND4)/YY(K)
         ENDIF
   30    CONTINUE
         IF(IL.GE.2) THEN
            IND1=((IL-2)/2)*L4+ABS(KN(NUM1+2))+I0-1
            IND2=((IL-2)/2)*L4+ABS(KN(NUM1+3))+I0-1
            IND3=((IL-2)/2)*L4+ABS(KN(NUM1+4))+I0-1
            IND4=((IL-2)/2)*L4+ABS(KN(NUM1+5))+I0-1
            DO 40 J0=1,IELEM
            JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            IF(KN(NUM1+2).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+2)))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(1,J0)*
     1         FUNKNO(IND1)/XX(K)
            ENDIF
            IF(KN(NUM1+3).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+3)))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(IELEM+1,J0)*
     1         FUNKNO(IND2)/XX(K)
            ENDIF
            JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            IF(KN(NUM1+4).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+4)))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(1,J0)*
     1         FUNKNO(IND3)/YY(K)
            ENDIF
            IF(KN(NUM1+5).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+5)))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(IL)*VOL0*V(IELEM+1,J0)*
     1         FUNKNO(IND4)/YY(K)
            ENDIF
   40       CONTINUE
         ENDIF
   50    CONTINUE
      ELSE IF(MOD(IL,2).EQ.1) THEN
         GARSI=SIGTI(IBM,MIN(IL+1,NAN))
         DO 150 I0=1,IELEM
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION OF
*        THE EVEN PARITY EQUATION.
         IF(IELEM.GT.1) THEN
            DO 65 J0=1,IELEM
            DO 60 K0=1,IELEM
            IF(QQ(J0,K0).EQ.0.0) GO TO 60
            JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            KND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+K0-1
            SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)**2)*VOL0*QQ(J0,K0)*
     1      GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
            JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            KND1=(IL/2)*L4+KN(NUM1+1)+(K0-1)*IELEM+I0-1
            SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)**2)*VOL0*QQ(J0,K0)*
     1      GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
            IF(IL.LE.NLF-3) THEN
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
               KND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+K0-1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
               KND1=(IL/2)*L4+KN(NUM1+1)+(K0-1)*IELEM+I0-1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
               JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+K0-1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
               JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K0-1)*IELEM+I0-1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL)*REAL(IL+1))*VOL0*
     1         QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+K0-1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL+1)*REAL(IL+1))*VOL0*
     1         QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
               JND1=((IL+2)/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
               KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K0-1)*IELEM+I0-1
               SUNKNO(JND1)=SUNKNO(JND1)+(REAL(IL+1)*REAL(IL+1))*VOL0*
     1         QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
            ENDIF
   60       CONTINUE
   65       CONTINUE
         ENDIF
*        ODD PARITY EQUATION.
         DO 75 IC=1,2
         IIC=1
         IF(IC.EQ.2) IIC=IELEM+1
         IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))+I0-1
         S1=REAL(SIGN(1,KN(NUM1+1+IC)))
         DO 70 JC=1,2
         JJC=1
         IF(JC.EQ.2) JJC=IELEM+1
         IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))+I0-1
         IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0)) THEN
            S2=REAL(SIGN(1,KN(NUM1+1+JC)))
            SUNKNO(IND1)=SUNKNO(IND1)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS*
     1      FUNKNO(IND2)
         ENDIF
   70    CONTINUE
   75    CONTINUE
         DO 85 IC=3,4
         IF(IC.EQ.3) IIC=1
         IF(IC.EQ.4) IIC=IELEM+1
         IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))+I0-1
         S1=REAL(SIGN(1,KN(NUM1+1+IC)))
         DO 80 JC=3,4
         IF(JC.EQ.3) JJC=1
         IF(JC.EQ.4) JJC=IELEM+1
         IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))+I0-1
         IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0)) THEN
            S2=REAL(SIGN(1,KN(NUM1+1+JC)))
            SUNKNO(IND1)=SUNKNO(IND1)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS*
     1      FUNKNO(IND2)
         ENDIF
   80    CONTINUE
   85    CONTINUE
         IF(ITY.EQ.1) GO TO 150
*
         IND1=(IL/2)*L4+ABS(KN(NUM1+2))+I0-1
         IND2=(IL/2)*L4+ABS(KN(NUM1+3))+I0-1
         IND3=(IL/2)*L4+ABS(KN(NUM1+4))+I0-1
         IND4=(IL/2)*L4+ABS(KN(NUM1+5))+I0-1
         IF((QFR(NUM2+1).NE.0.0).AND.(KN(NUM1+2).NE.0)) THEN
*           XINF SIDE.
            DO 90 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+2))+I0-1
            SUNKNO(IND1)=SUNKNO(IND1)-0.5*FACT*QFR(NUM2+1)*ZMARS*
     1      FUNKNO(IND5)
   90       CONTINUE
         ENDIF
         IF((QFR(NUM2+2).NE.0.0).AND.(KN(NUM1+3).NE.0)) THEN
*           XSUP SIDE.
            DO 100 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+3))+I0-1
            SUNKNO(IND2)=SUNKNO(IND2)-0.5*FACT*QFR(NUM2+2)*ZMARS*
     1      FUNKNO(IND5)
  100       CONTINUE
         ENDIF
         IF((QFR(NUM2+3).NE.0.0).AND.(KN(NUM1+4).NE.0)) THEN
*           YINF SIDE.
            DO 110 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+4))+I0-1
            SUNKNO(IND3)=SUNKNO(IND3)-0.5*FACT*QFR(NUM2+3)*ZMARS*
     1      FUNKNO(IND5)
  110       CONTINUE
         ENDIF
         IF((QFR(NUM2+4).NE.0.0).AND.(KN(NUM1+5).NE.0)) THEN
*           YSUP SIDE.
            DO 120 IL2=1,NLF-1,2
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+5))+I0-1
            SUNKNO(IND4)=SUNKNO(IND4)-0.5*FACT*QFR(NUM2+4)*ZMARS*
     1      FUNKNO(IND5)
  120       CONTINUE
         ENDIF
*
         DO 130 J0=1,IELEM
         JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
         IF(KN(NUM1+2).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+2)))
            SUNKNO(IND1)=SUNKNO(IND1)+SG*REAL(IL)*VOL0*V(1,J0)*
     1      FUNKNO(JND1)/XX(K)
         ENDIF
         IF(KN(NUM1+3).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+3)))
            SUNKNO(IND2)=SUNKNO(IND2)+SG*REAL(IL)*VOL0*V(IELEM+1,J0)*
     1      FUNKNO(JND1)/XX(K)
         ENDIF
         JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
         IF(KN(NUM1+4).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+4)))
            SUNKNO(IND3)=SUNKNO(IND3)+SG*REAL(IL)*VOL0*V(1,J0)*
     1      FUNKNO(JND1)/YY(K)
         ENDIF
         IF(KN(NUM1+5).NE.0) THEN
            SG=REAL(SIGN(1,KN(NUM1+5)))
            SUNKNO(IND4)=SUNKNO(IND4)+SG*REAL(IL)*VOL0*V(IELEM+1,J0)*
     1      FUNKNO(JND1)/YY(K)
         ENDIF
  130    CONTINUE
         IF(IL.LE.NLF-3) THEN
            DO 140 J0=1,IELEM
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            IF(KN(NUM1+2).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+2)))
               SUNKNO(IND1)=SUNKNO(IND1)+SG*REAL(IL+1)*VOL0*V(1,J0)*
     1         FUNKNO(JND1)/XX(K)
            ENDIF
            IF(KN(NUM1+3).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+3)))
               SUNKNO(IND2)=SUNKNO(IND2)+SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,J0)*FUNKNO(JND1)/XX(K)
            ENDIF
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            IF(KN(NUM1+4).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+4)))
               SUNKNO(IND3)=SUNKNO(IND3)+SG*REAL(IL+1)*VOL0*V(1,J0)*
     1         FUNKNO(JND1)/YY(K)
            ENDIF
            IF(KN(NUM1+5).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+5)))
               SUNKNO(IND4)=SUNKNO(IND4)+SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,J0)*FUNKNO(JND1)/YY(K)
            ENDIF
  140       CONTINUE
         ENDIF
  150    CONTINUE
      ENDIF
      NUM1=NUM1+5
      NUM2=NUM2+4
  160 CONTINUE
  170 CONTINUE
      RETURN
      END
