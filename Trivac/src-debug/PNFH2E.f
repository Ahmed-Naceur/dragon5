*DECK PNFH2E
      SUBROUTINE PNFH2E (IELEM,ICOL,NBLOS,SIDE,NLF,NVD,L4,IPERT,KN,
     1 QFR,MU,IIMAX,LC,V,SYS,SUNKNO,FUNKNO,NADI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a one-group SPN flux iteration in hexagonal 2D geometry.
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
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* SIDE    side of the hexagons.
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* MU      profiled storage indices for matrix SYS.
* IIMAX   dimension of SYS.
* LC      order of the unit matrices.
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
      INTEGER IELEM,ICOL,NBLOS,NLF,NVD,L4,IPERT(NBLOS),
     1 KN(NBLOS,4+6*IELEM*(IELEM+1)),MU(L4),IIMAX,LC,NADI
      REAL SIDE,QFR(NBLOS,6),V(LC,LC-1),SYS(IIMAX),SUNKNO(L4*NLF/2),
     1 FUNKNO(L4*NLF/2)
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
      MUMAX=MU(L4)
      NELEM=IELEM*(IELEM+1)
      DO 170 IADI=1,MAX(1,NADI)
      DO 160 IL=0,NLF-1
      FACT=REAL(2*IL+1)
      IF(MOD(IL,2).EQ.0) THEN
         DO 10 I=1,L4
         FUNKNO((IL/2)*L4+I)=SUNKNO((IL/2)*L4+I)
   10    CONTINUE
      ENDIF
*----
*  COMPUTE THE SOLUTION AT ORDER IL.
*----
      NUM=0
      DO 150 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 150
      NUM=NUM+1
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION
         IF(IL.GE.2) THEN
            DO 33 K4=0,1
            DO 32 K3=0,IELEM-1
            DO 31 K2=1,IELEM+1
            KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
            KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
            KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
            INW1=((IL-2)/2)*L4+ABS(KNW1)
            INX1=((IL-2)/2)*L4+ABS(KNX1)
            INY1=((IL-2)/2)*L4+ABS(KNY1)
            DO 30 K1=0,IELEM-1
            IF(V(K2,K1+1).EQ.0.0) GO TO 30
            IF(K4.EQ.0) THEN
               SSS=(-1.0)**K1
               JND1=(IL/2)*L4+KN(NUM,1)+K3*IELEM+K1
               JND2=(IL/2)*L4+KN(NUM,2)+K3*IELEM+K1
               JND3=(IL/2)*L4+KN(NUM,3)+K3*IELEM+K1
            ELSE
               SSS=1.0
               JND1=(IL/2)*L4+KN(NUM,2)+K1*IELEM+K3
               JND2=(IL/2)*L4+KN(NUM,3)+K1*IELEM+K3
               JND3=(IL/2)*L4+KN(NUM,4)+K1*IELEM+K3
            ENDIF
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               FUNKNO(JND1)=FUNKNO(JND1)-SG*SSS*REAL(IL)*SIDE*
     1         V(K2,K1+1)*FUNKNO(INW1)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               FUNKNO(JND2)=FUNKNO(JND2)-SG*SSS*REAL(IL)*SIDE*
     1         V(K2,K1+1)*FUNKNO(INX1)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               FUNKNO(JND3)=FUNKNO(JND3)-SG*SSS*REAL(IL)*SIDE*
     1         V(K2,K1+1)*FUNKNO(INY1)
            ENDIF
   30       CONTINUE
   31       CONTINUE
   32       CONTINUE
   33       CONTINUE
         ENDIF
      ELSE
*        ODD PARITY EQUATION
         DO 142 K4=0,1
         DO 141 K3=0,IELEM-1
         DO 140 K2=1,IELEM+1
         KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
         KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
         INW1=(IL/2)*L4+ABS(KNW1)
         INX1=(IL/2)*L4+ABS(KNX1)
         INY1=(IL/2)*L4+ABS(KNY1)
         IF(KNW1.NE.0) THEN
            DO 90 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 90
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INW2=(IL2/2)*L4+ABS(KNW1)
            IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
               FUNKNO(INW1)=FUNKNO(INW1)+0.5*FACT*QFR(NUM,1)*ZMARS*
     1         FUNKNO(INW2)
            ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
               FUNKNO(INW1)=FUNKNO(INW1)+0.5*FACT*QFR(NUM,2)*ZMARS*
     1         FUNKNO(INW2)
            ENDIF
   90       CONTINUE
         ENDIF
         IF(KNX1.NE.0) THEN
            DO 100 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 100
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INX2=(IL2/2)*L4+ABS(KNX1)
            IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
               FUNKNO(INX1)=FUNKNO(INX1)+0.5*FACT*QFR(NUM,3)*ZMARS*
     1         FUNKNO(INX2)
            ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
               FUNKNO(INX1)=FUNKNO(INX1)+0.5*FACT*QFR(NUM,4)*ZMARS*
     1         FUNKNO(INX2)
            ENDIF
  100       CONTINUE
         ENDIF
         IF(KNY1.NE.0) THEN
            DO 110 IL2=1,NLF-1,2
            IF(IL2.EQ.IL) GO TO 110
            ZMARS=PNMAR2(NZMAR,IL2,IL)
            INY2=(IL2/2)*L4+ABS(KNY1)
            IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
               FUNKNO(INY1)=FUNKNO(INY1)+0.5*FACT*QFR(NUM,5)*ZMARS*
     1         FUNKNO(INY2)
            ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
               FUNKNO(INY1)=FUNKNO(INY1)+0.5*FACT*QFR(NUM,6)*ZMARS*
     1         FUNKNO(INY2)
            ENDIF
  110       CONTINUE
         ENDIF
         IF(IL.LE.NLF-3) THEN
            DO 130 K1=0,IELEM-1
            IF(V(K2,K1+1).EQ.0.0) GO TO 130
            IF(K4.EQ.0) THEN
               SSS=(-1.0)**K1
               JND1=((IL+2)/2)*L4+KN(NUM,1)+K3*IELEM+K1
               JND2=((IL+2)/2)*L4+KN(NUM,2)+K3*IELEM+K1
               JND3=((IL+2)/2)*L4+KN(NUM,3)+K3*IELEM+K1
            ELSE
               SSS=1.0
               JND1=((IL+2)/2)*L4+KN(NUM,2)+K1*IELEM+K3
               JND2=((IL+2)/2)*L4+KN(NUM,3)+K1*IELEM+K3
               JND3=((IL+2)/2)*L4+KN(NUM,4)+K1*IELEM+K3
            ENDIF
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               FUNKNO(INW1)=FUNKNO(INW1)-SG*SSS*REAL(IL+1)*SIDE*
     1         V(K2,K1+1)*FUNKNO(JND1)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               FUNKNO(INX1)=FUNKNO(INX1)-SG*SSS*REAL(IL+1)*SIDE*
     1         V(K2,K1+1)*FUNKNO(JND2)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               FUNKNO(INY1)=FUNKNO(INY1)-SG*SSS*REAL(IL+1)*SIDE*
     1         V(K2,K1+1)*FUNKNO(JND3)
            ENDIF
  130       CONTINUE
         ENDIF
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
      ENDIF
  150 CONTINUE
      IF(MOD(IL,2).EQ.1) THEN
         CALL ALLDLS(L4,MU,SYS((IL/2)*MUMAX+1),FUNKNO((IL/2)*L4+1))
      ENDIF
  160 CONTINUE
  170 CONTINUE
      RETURN
      END
