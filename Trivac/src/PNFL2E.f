*DECK PNFL2E
      SUBROUTINE PNFL2E (NREG,IELEM,ICOL,XX,YY,MAT,VOL,NBMIX,NLF,NVD,
     1 NAN,SIGTI,L4,KN,QFR,MU,IIMAX,LC,R,V,SYS,SUNKNO,FUNKNO,NADI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a one-group SPN flux iteration in Cartesian 2D geometry.
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
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only.
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* MU      profiled storage indices for matrix SYS.
* IIMAX   dimension of SYS.
* LC      order of the unit matrices.
* R       unit Cartesian mass matrix.
* V       unit nodal coupling matrix.
* SYS     LU factors of the system matrix.
* SUNKNO  sources.
* FUNKNO  initial fluxes.
* NADI    number of inner ADI iterations.
*
*Parameters: output
* FUNKNO  fluxes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,IELEM,ICOL,MAT(NREG),NBMIX,NLF,NVD,NAN,L4,KN(5*NREG),
     1 MU(L4),IIMAX,LC,NADI
      REAL XX(NREG),YY(NREG),VOL(NREG),SIGTI(NBMIX,NAN),QFR(4*NREG),
     1 R(LC,LC),V(LC,LC-1),SYS(IIMAX),SUNKNO(L4*NLF/2),FUNKNO(L4*NLF/2)
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
      MUMAX=MU(L4)
      DO 170 IADI=1,MAX(1,NADI)
      DO 160 IL=0,NLF-1
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
            DO 35 I0=1,IELEM
            IND1=((IL-2)/2)*L4+ABS(KN(NUM1+2))+I0-1
            IND2=((IL-2)/2)*L4+ABS(KN(NUM1+3))+I0-1
            IND3=((IL-2)/2)*L4+ABS(KN(NUM1+4))+I0-1
            IND4=((IL-2)/2)*L4+ABS(KN(NUM1+5))+I0-1
            DO 30 J0=1,IELEM
            JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            IF(KN(NUM1+2).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+2)))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(1,J0)*
     1         FUNKNO(IND1)/XX(K)
            ENDIF
            IF(KN(NUM1+3).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+3)))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(IELEM+1,J0)*
     1         FUNKNO(IND2)/XX(K)
            ENDIF
            JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            IF(KN(NUM1+4).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+4)))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(1,J0)*
     1         FUNKNO(IND3)/YY(K)
            ENDIF
            IF(KN(NUM1+5).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+5)))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*REAL(IL)*VOL0*V(IELEM+1,J0)*
     1         FUNKNO(IND4)/YY(K)
            ENDIF
   30       CONTINUE
   35       CONTINUE
         ENDIF
      ELSE
         DO 140 I0=1,IELEM
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION
*        OF THE EVEN PARITY EQUATION.
         IF((IL.GE.3).AND.(IELEM.GT.1)) THEN
            GARSI=SIGTI(IBM,MIN(IL-1,NAN))
            DO 65 J0=1,IELEM
            DO 50 K0=1,IELEM
            IF(QQ(J0,K0).EQ.0.0) GO TO 50
            JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            KND1=((IL-2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+K0-1
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL-1)*REAL(IL-2))*VOL0*
     1      QQ(J0,K0)*GARSI*FUNKNO(KND1)/(REAL(2*IL-3)*XX(K)*XX(K))
   50       CONTINUE
            DO 60 K0=1,IELEM
            IF(QQ(J0,K0).EQ.0.0) GO TO 60
            JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            KND1=((IL-2)/2)*L4+KN(NUM1+1)+(K0-1)*IELEM+I0-1
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL-1)*REAL(IL-2))*VOL0*
     1      QQ(J0,K0)*GARSI*FUNKNO(KND1)/(REAL(2*IL-3)*YY(K)*YY(K))
   60       CONTINUE
   65       CONTINUE
         ENDIF
         IF((IL.LE.NLF-3).AND.(IELEM.GT.1)) THEN
            GARSI=SIGTI(IBM,MIN(IL+1,NAN))
            DO 85 J0=1,IELEM
            DO 70 K0=1,IELEM
            IF(QQ(J0,K0).EQ.0.0) GO TO 70
            JND1=(IL/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            KND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+K0-1
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL)*REAL(IL+1))*VOL0*
     1      QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*XX(K)*XX(K))
   70       CONTINUE
            DO 80 K0=1,IELEM
            IF(QQ(J0,K0).EQ.0.0) GO TO 80
            JND1=(IL/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            KND1=((IL+2)/2)*L4+KN(NUM1+1)+(K0-1)*IELEM+I0-1
            FUNKNO(JND1)=FUNKNO(JND1)-(REAL(IL)*REAL(IL+1))*VOL0*
     1      QQ(J0,K0)*GARSI*FUNKNO(KND1)/(FACT*YY(K)*YY(K))
   80       CONTINUE
   85       CONTINUE
         ENDIF
*
*        ODD PARITY EQUATION
         IND1=(IL/2)*L4+ABS(KN(NUM1+2))+I0-1
         IND2=(IL/2)*L4+ABS(KN(NUM1+3))+I0-1
         IND3=(IL/2)*L4+ABS(KN(NUM1+4))+I0-1
         IND4=(IL/2)*L4+ABS(KN(NUM1+5))+I0-1
         IF((QFR(NUM2+1).NE.0.0).AND.(KN(NUM1+2).NE.0)) THEN
*           XINF SIDE.
            DO 90 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 90
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+2))+I0-1
            FUNKNO(IND1)=FUNKNO(IND1)+0.5*FACT*QFR(NUM2+1)*ZMARS*
     1      FUNKNO(IND5)
   90       CONTINUE
         ENDIF
         IF((QFR(NUM2+2).NE.0.0).AND.(KN(NUM1+3).NE.0)) THEN
*           XSUP SIDE.
            DO 100 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 100
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+3))+I0-1
            FUNKNO(IND2)=FUNKNO(IND2)+0.5*FACT*QFR(NUM2+2)*ZMARS*
     1      FUNKNO(IND5)
  100       CONTINUE
         ENDIF
         IF((QFR(NUM2+3).NE.0.0).AND.(KN(NUM1+4).NE.0)) THEN
*           YINF SIDE.
            DO 110 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 110
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+4))+I0-1
            FUNKNO(IND3)=FUNKNO(IND3)+0.5*FACT*QFR(NUM2+3)*ZMARS*
     1      FUNKNO(IND5)
  110       CONTINUE
         ENDIF
         IF((QFR(NUM2+4).NE.0.0).AND.(KN(NUM1+5).NE.0)) THEN
*           YSUP SIDE.
            DO 120 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 120
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            IND5=(IL2/2)*L4+ABS(KN(NUM1+5))+I0-1
            FUNKNO(IND4)=FUNKNO(IND4)+0.5*FACT*QFR(NUM2+4)*ZMARS*
     1      FUNKNO(IND5)
  120       CONTINUE
         ENDIF
         IF(IL.LE.NLF-3) THEN
            DO 130 J0=1,IELEM
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(I0-1)*IELEM+J0-1
            IF(KN(NUM1+2).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+2)))
               FUNKNO(IND1)=FUNKNO(IND1)-SG*REAL(IL+1)*VOL0*V(1,J0)*
     1         FUNKNO(JND1)/XX(K)
            ENDIF
            IF(KN(NUM1+3).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+3)))
               FUNKNO(IND2)=FUNKNO(IND2)-SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,J0)*FUNKNO(JND1)/XX(K)
            ENDIF
            JND1=((IL+2)/2)*L4+KN(NUM1+1)+(J0-1)*IELEM+I0-1
            IF(KN(NUM1+4).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+4)))
               FUNKNO(IND3)=FUNKNO(IND3)-SG*REAL(IL+1)*VOL0*V(1,J0)*
     1         FUNKNO(JND1)/YY(K)
            ENDIF
            IF(KN(NUM1+5).NE.0) THEN
               SG=REAL(SIGN(1,KN(NUM1+5)))
               FUNKNO(IND4)=FUNKNO(IND4)-SG*REAL(IL+1)*VOL0*
     1         V(IELEM+1,J0)*FUNKNO(JND1)/YY(K)
            ENDIF
  130       CONTINUE
         ENDIF
  140    CONTINUE
      ENDIF
      NUM1=NUM1+5
      NUM2=NUM2+4
  150 CONTINUE
      IF(MOD(IL,2).EQ.1) THEN
         CALL ALLDLS(L4,MU,SYS((IL/2)*MUMAX+1),FUNKNO((IL/2)*L4+1))
      ENDIF
  160 CONTINUE
  170 CONTINUE
      RETURN
      END
