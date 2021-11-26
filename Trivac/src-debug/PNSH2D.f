*DECK PNSH2D
      SUBROUTINE PNSH2D(ITY,IELEM,ICOL,NBLOS,SIDE,MAT,NBMIX,NLF,NVD,
     1 NAN,SIGT,L4,IPERT,KN,QFR,LC,R,V,H,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Source calculation for a SPN approximation in BIVAC, including
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
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* SIDE    side of the hexagons.
* MAT     index-number of the mixture type assigned to each volume.
* NBMIX   number of mixtures.
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* SIGT    macroscopic cross sections ordered by mixture.
*         SIGT(:,NAN) generally contains the total cross section only.
* L4      order of the profiled system matrices.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* LC      order of the unit matrices.
* R       unit Cartesian mass matrix.
* V       unit nodal coupling matrix.
* H       Piolat (hexagonal) coupling matrix.
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
      INTEGER ITY,IELEM,ICOL,NBLOS,MAT(3,NBLOS),NBMIX,NLF,NVD,NAN,L4,
     1 IPERT(NBLOS),KN(NBLOS,4+6*IELEM*(IELEM+1)),LC
      REAL SIDE,SIGT(NBMIX,NAN),QFR(NBLOS,6),R(LC,LC),V(LC,LC-1),
     1 H(LC,LC-1),SUNKNO(L4*NLF/2),FUNKNO(L4*NLF/2)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXIEL=3)
      DOUBLE PRECISION CTRAN(MAXIEL*(MAXIEL+1),MAXIEL*(MAXIEL+1)),VAR1
*
      TTTT=REAL(0.5D0*SQRT(3.D00)*SIDE*SIDE)
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
      NELEM=IELEM*(IELEM+1)
      COEF=REAL(2.0D0*SIDE*SIDE/SQRT(3.D00))
*----
*  COMPUTE THE TRANVERSE COUPLING PIOLAT UNIT MATRIX
*----
      CALL XDDSET(CTRAN,MAXIEL**2,0.0D0)
      CNORM=REAL(SIDE*SIDE/SQRT(3.D00))
      I=0
      DO 22 JS=1,IELEM
      DO 21 JT=1,IELEM+1
      J=0
      I=I+1
      SSS=1.0
      DO 20 IT=1,IELEM
      DO 10 IS=1,IELEM+1
      J=J+1
      CTRAN(I,J)=SSS*CNORM*H(IS,JS)*H(JT,IT)
   10 CONTINUE
      SSS=-SSS
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
*
      DO 160 IL=0,NLF-1
      IF((ITY.EQ.1).AND.(IL.GE.NAN)) GO TO 160
      FACT=REAL(2*IL+1)
*----
*  COMPUTE THE SOURCE AT ORDER IL.
*----
      NUM=0
      DO 150 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 150
      NUM=NUM+1
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 150
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION.
         DO 35 K2=0,IELEM-1
         DO 30 K1=0,IELEM-1
         JND1=(IL/2)*L4+KN(NUM,1)+K2*IELEM+K1 ! w-oriented flux
         JND2=(IL/2)*L4+KN(NUM,2)+K2*IELEM+K1
         JND3=(IL/2)*L4+KN(NUM,3)+K2*IELEM+K1
         SUNKNO(JND1)=SUNKNO(JND1)+FACT*TTTT*GARS*FUNKNO(JND1)
         SUNKNO(JND2)=SUNKNO(JND2)+FACT*TTTT*GARS*FUNKNO(JND2)
         SUNKNO(JND3)=SUNKNO(JND3)+FACT*TTTT*GARS*FUNKNO(JND3)
   30    CONTINUE
   35    CONTINUE
         IF(ITY.EQ.1) GO TO 150
*
         DO 43 K4=0,1
         DO 42 K3=0,IELEM-1
         DO 41 K2=1,IELEM+1
         KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
         KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
         INW1=(IL/2)*L4+ABS(KNW1) ! w-oriented current
         INX1=(IL/2)*L4+ABS(KNX1)
         INY1=(IL/2)*L4+ABS(KNY1)
         DO 40 K1=0,IELEM-1
         IF(V(K2,K1+1).EQ.0.0) GO TO 40
         IF(K4.EQ.0) THEN
            SSS=(-1.0)**K1
            JND1=(IL/2)*L4+KN(NUM,1)+K3*IELEM+K1 ! w-oriented flux
            JND2=(IL/2)*L4+KN(NUM,2)+K3*IELEM+K1
            JND3=(IL/2)*L4+KN(NUM,3)+K3*IELEM+K1
         ELSE
            SSS=1.0
            JND1=(IL/2)*L4+KN(NUM,2)+K1*IELEM+K3
            JND2=(IL/2)*L4+KN(NUM,3)+K1*IELEM+K3
            JND3=(IL/2)*L4+KN(NUM,4)+K1*IELEM+K3
         ENDIF
         VAR1=SSS*REAL(IL+1)*SIDE*V(K2,K1+1)
         IF(KNW1.NE.0) THEN
            SG=REAL(SIGN(1,KNW1))
            SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(VAR1)*FUNKNO(INW1)
         ENDIF
         IF(KNX1.NE.0) THEN
            SG=REAL(SIGN(1,KNX1))
            SUNKNO(JND2)=SUNKNO(JND2)+SG*REAL(VAR1)*FUNKNO(INX1)
         ENDIF
         IF(KNY1.NE.0) THEN
            SG=REAL(SIGN(1,KNY1))
            SUNKNO(JND3)=SUNKNO(JND3)+SG*REAL(VAR1)*FUNKNO(INY1)
         ENDIF
         IF(IL.GE.2) THEN
            VAR1=SSS*REAL(IL)*SIDE*V(K2,K1+1)
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               SUNKNO(JND1)=SUNKNO(JND1)+SG*REAL(VAR1)*FUNKNO(INW1-L4)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               SUNKNO(JND2)=SUNKNO(JND2)+SG*REAL(VAR1)*FUNKNO(INX1-L4)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               SUNKNO(JND3)=SUNKNO(JND3)+SG*REAL(VAR1)*FUNKNO(INY1-L4)
            ENDIF
         ENDIF
   40    CONTINUE
   41    CONTINUE
   42    CONTINUE
   43    CONTINUE
      ELSE IF(MOD(IL,2).EQ.1) THEN
*        ODD PARITY EQUATION.
         DO 112 K4=0,1 ! TWO LOZENGES PER HEXAGON
         DO 111 K3=0,IELEM-1
         DO 110 K2=1,IELEM+1
         KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2) ! w-oriented current
         KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
         INW1=(IL/2)*L4+ABS(KNW1)
         INX1=(IL/2)*L4+ABS(KNX1)
         INY1=(IL/2)*L4+ABS(KNY1)
         DO 70 K1=1,IELEM+1
         KNW2=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K1) ! w-oriented current
         KNX2=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K1)
         KNY2=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K1)
         INW2=(IL/2)*L4+ABS(KNW2)
         INX2=(IL/2)*L4+ABS(KNX2)
         INY2=(IL/2)*L4+ABS(KNY2)
         VAR1=FACT*COEF*GARS*R(K2,K1)
         IF((KNW2.NE.0).AND.(KNW1.NE.0)) THEN
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNW2))
            SUNKNO(INW1)=SUNKNO(INW1)-SG*REAL(VAR1)*FUNKNO(INW2)
         ENDIF
         IF((KNX2.NE.0).AND.(KNX1.NE.0)) THEN
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNX2))
            SUNKNO(INX1)=SUNKNO(INX1)-SG*REAL(VAR1)*FUNKNO(INX2)
         ENDIF
         IF((KNY2.NE.0).AND.(KNY1.NE.0)) THEN
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNY2))
            SUNKNO(INY1)=SUNKNO(INY1)-SG*REAL(VAR1)*FUNKNO(INY2)
         ENDIF
   70    CONTINUE
         IF(ITY.EQ.0) THEN
*           BOUNDARY CONDITIONS.
            IF(KNW1.NE.0) THEN
               DO 80 IL2=1,NLF-1,2
               ZMARS=PNMAR2(NZMAR,IL2,IL)
               INW2=(IL2/2)*L4+ABS(KNW1)
               IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
                  VAR1=0.5*FACT*QFR(NUM,1)*ZMARS*FUNKNO(INW2)
                  SUNKNO(INW1)=SUNKNO(INW1)-REAL(VAR1)
               ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
                  VAR1=0.5*FACT*QFR(NUM,2)*ZMARS*FUNKNO(INW2)
                  SUNKNO(INW1)=SUNKNO(INW1)-REAL(VAR1)
               ENDIF
   80          CONTINUE
            ENDIF
            IF(KNX1.NE.0) THEN
               DO 90 IL2=1,NLF-1,2
               ZMARS=PNMAR2(NZMAR,IL2,IL)
               INX2=(IL2/2)*L4+ABS(KNX1)
               IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
                  VAR1=0.5*FACT*QFR(NUM,3)*ZMARS*FUNKNO(INX2)
                  SUNKNO(INX1)=SUNKNO(INX1)-REAL(VAR1)
               ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
                  VAR1=0.5*FACT*QFR(NUM,4)*ZMARS*FUNKNO(INX2)
                  SUNKNO(INX1)=SUNKNO(INX1)-REAL(VAR1)
               ENDIF
   90          CONTINUE
            ENDIF
            IF(KNY1.NE.0) THEN
               DO 100 IL2=1,NLF-1,2
               ZMARS=PNMAR2(NZMAR,IL2,IL)
               INY2=(IL2/2)*L4+ABS(KNY1)
               IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
                  VAR1=0.5*FACT*QFR(NUM,5)*ZMARS*FUNKNO(INY2)
                  SUNKNO(INY1)=SUNKNO(INY1)-REAL(VAR1)
               ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
                  VAR1=0.5*FACT*QFR(NUM,6)*ZMARS*FUNKNO(INY2)
                  SUNKNO(INY1)=SUNKNO(INY1)-REAL(VAR1)
               ENDIF
  100          CONTINUE
            ENDIF
         ENDIF
  110    CONTINUE
  111    CONTINUE
  112    CONTINUE
*
         ITRS=0
         DO I=1,NBLOS
            IF(KN(I,1).EQ.KN(NUM,4)) THEN
               ITRS=I
               GO TO 120
            ENDIF
         ENDDO
         CALL XABORT('PNDH2E: ITRS FAILURE.')
  120    DO 135 I=1,NELEM
         KNW1=KN(ITRS,4+I)
         KNX1=KN(NUM,4+2*NELEM+I)
         KNY1=KN(NUM,4+4*NELEM+I)
         INW1=(IL/2)*L4+ABS(KNW1)
         INX1=(IL/2)*L4+ABS(KNX1)
         INY1=(IL/2)*L4+ABS(KNY1)
         DO 130 J=1,NELEM
         KNW2=KN(NUM,4+NELEM+J)
         KNX2=KN(NUM,4+3*NELEM+J)
         KNY2=KN(NUM,4+5*NELEM+J)
         INW2=(IL/2)*L4+ABS(KNW2)
         INX2=(IL/2)*L4+ABS(KNX2)
         INY2=(IL/2)*L4+ABS(KNY2)
         VAR1=FACT*GARS*CTRAN(I,J)
         IF((KNY2.NE.0).AND.(KNW1.NE.0)) THEN
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
            SUNKNO(INY2)=SUNKNO(INY2)-SG*REAL(VAR1)*FUNKNO(INW1) ! y w
            SUNKNO(INW1)=SUNKNO(INW1)-SG*REAL(VAR1)*FUNKNO(INY2) ! w y
         ENDIF
         IF((KNW2.NE.0).AND.(KNX1.NE.0)) THEN
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
            SUNKNO(INX1)=SUNKNO(INX1)-SG*REAL(VAR1)*FUNKNO(INW2) ! x w
            SUNKNO(INW2)=SUNKNO(INW2)-SG*REAL(VAR1)*FUNKNO(INX1) ! w x
         ENDIF
         IF((KNX2.NE.0).AND.(KNY1.NE.0)) THEN
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
            SUNKNO(INY1)=SUNKNO(INY1)-SG*REAL(VAR1)*FUNKNO(INX2) ! y x
            SUNKNO(INX2)=SUNKNO(INX2)-SG*REAL(VAR1)*FUNKNO(INY1) ! x y
         ENDIF
  130    CONTINUE
  135    CONTINUE
         IF(ITY.EQ.1) GO TO 150
*
         DO 143 K4=0,1
         DO 142 K3=0,IELEM-1
         DO 141 K2=1,IELEM+1
         KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
         KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
         INW1=(IL/2)*L4+ABS(KNW1) ! w-oriented current
         INX1=(IL/2)*L4+ABS(KNX1)
         INY1=(IL/2)*L4+ABS(KNY1)
         DO 140 K1=0,IELEM-1
         IF(V(K2,K1+1).EQ.0.0) GO TO 140
         IF(K4.EQ.0) THEN
            SSS=(-1.0)**K1
            JND1=(IL/2)*L4+KN(NUM,1)+K3*IELEM+K1 ! w-oriented flux
            JND2=(IL/2)*L4+KN(NUM,2)+K3*IELEM+K1
            JND3=(IL/2)*L4+KN(NUM,3)+K3*IELEM+K1
         ELSE
            SSS=1.0
            JND1=(IL/2)*L4+KN(NUM,2)+K1*IELEM+K3
            JND2=(IL/2)*L4+KN(NUM,3)+K1*IELEM+K3
            JND3=(IL/2)*L4+KN(NUM,4)+K1*IELEM+K3
         ENDIF
         VAR1=SSS*REAL(IL)*SIDE*V(K2,K1+1)
         IF(KNW1.NE.0) THEN
            SG=REAL(SIGN(1,KNW1))
            SUNKNO(INW1)=SUNKNO(INW1)+SG*REAL(VAR1)*FUNKNO(JND1)
         ENDIF
         IF(KNX1.NE.0) THEN
            SG=REAL(SIGN(1,KNX1))
            SUNKNO(INX1)=SUNKNO(INX1)+SG*REAL(VAR1)*FUNKNO(JND2)
         ENDIF
         IF(KNY1.NE.0) THEN
            SG=REAL(SIGN(1,KNY1))
            SUNKNO(INY1)=SUNKNO(INY1)+SG*REAL(VAR1)*FUNKNO(JND3)
         ENDIF
         IF(IL.LE.NLF-3) THEN
            VAR1=SSS*REAL(IL+1)*SIDE*V(K2,K1+1)
            IF(KNW1.NE.0) THEN
               SG=REAL(SIGN(1,KNW1))
               SUNKNO(INW1)=SUNKNO(INW1)+SG*REAL(VAR1)*FUNKNO(JND1+L4)
            ENDIF
            IF(KNX1.NE.0) THEN
               SG=REAL(SIGN(1,KNX1))
               SUNKNO(INX1)=SUNKNO(INX1)+SG*REAL(VAR1)*FUNKNO(JND2+L4)
            ENDIF
            IF(KNY1.NE.0) THEN
               SG=REAL(SIGN(1,KNY1))
               SUNKNO(INY1)=SUNKNO(INY1)+SG*REAL(VAR1)*FUNKNO(JND3+L4)
            ENDIF
         ENDIF
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
  143    CONTINUE
      ENDIF
  150 CONTINUE
  160 CONTINUE
      RETURN
      END
