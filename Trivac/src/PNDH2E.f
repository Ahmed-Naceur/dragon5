*DECK PNDH2E
      SUBROUTINE PNDH2E(ITY,IELEM,ICOL,NBLOS,L4,NBMIX,IIMAX,SIDE,MAT,
     1 IPERT,SIGT,KN,QFR,NLF,NVD,NAN,MU,LC,R,V,H,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a within-group (leakage and removal) or out-of-group
* system matrix in a Thomas-Raviart-Schneider (dual) finite element
* simplified PN method approximation (2D hexagonal geometry).
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
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* NBMIX   number of mixtures.
* IIMAX   allocated dimension of array SYS.
* SIDE    side of the hexagons.
* MAT     mixture index assigned to each element.
* SIGT    total minus self-scattering macroscopic cross sections.
*         SIGT(:,NAN) generally contains the total cross section only.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* MU      indices used with compressed diagonal storage mode matrix SYS.
* LC      order of the unit matrices.
* R       Cartesian mass matrix.
* V       nodal coupling matrix.
* H       Piolat (hexagonal) coupling matrix.
*
*Parameters: output
* SYS     system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,IELEM,ICOL,NBLOS,L4,NBMIX,IIMAX,MAT(3,NBLOS),
     1 IPERT(NBLOS),KN(NBLOS,4+6*IELEM*(IELEM+1)),NLF,NVD,NAN,MU(L4),LC
      REAL SIDE,SIGT(NBMIX,NAN),QFR(NBLOS,6),R(LC,LC),V(LC,LC-1),
     1 H(LC,LC-1),SYS(IIMAX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXIEL=3)
      DOUBLE PRECISION CTRAN(MAXIEL*(MAXIEL+1),MAXIEL*(MAXIEL+1)),VAR1
*
      TTTT=REAL(0.5D0*SQRT(3.D00)*SIDE*SIDE)
      IF(IELEM.GT.MAXIEL) CALL XABORT('PNDH2E: MAXIEL OVERFLOW.')
      NZMAR=65
      IF(ICOL.EQ.3) THEN
         IF(NVD.EQ.0) THEN
            NZMAR=NLF+1
         ELSE IF(NVD.EQ.1) THEN
            NZMAR=NLF
         ELSE IF(NVD.EQ.2) THEN
            NZMAR=65
         ENDIF
      ENDIF
      MUMAX=MU(L4)
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
*----
*  ASSEMBLY OF THE MAIN COEFFICIENT MATRIX AT ORDER IL.
*----
      DO 100 IL=0,NLF-1
      ZMARS=0.0
      IF(MOD(IL,2).EQ.1) ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
      NUM=0
      KEY=0
      DO 90 KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) GO TO 90
      IBM=MAT(1,IPERT(KEL))
      IF(IBM.EQ.0) GO TO 90
      NUM=NUM+1
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION.
         DO 35 K2=0,IELEM-1
         DO 30 K1=0,IELEM-1
         JND1=KN(NUM,1)+K2*IELEM+K1
         JND2=KN(NUM,2)+K2*IELEM+K1
         JND3=KN(NUM,3)+K2*IELEM+K1
         KEY=(IL/2)*MUMAX+MU(JND1)
         SYS(KEY)=SYS(KEY)+FACT*TTTT*GARS
         KEY=(IL/2)*MUMAX+MU(JND2)
         SYS(KEY)=SYS(KEY)+FACT*TTTT*GARS
         KEY=(IL/2)*MUMAX+MU(JND3)
         SYS(KEY)=SYS(KEY)+FACT*TTTT*GARS
   30    CONTINUE
   35    CONTINUE
      ELSE
*        ODD PARITY EQUATION.
         DO 52 K4=0,1
         DO 51 K3=0,IELEM-1
         DO 50 K2=1,IELEM+1
         KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
         KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
         INW1=ABS(KNW1)
         INX1=ABS(KNX1)
         INY1=ABS(KNY1)
         DO 40 K1=1,IELEM+1
         KNW2=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K1)
         KNX2=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K1)
         KNY2=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K1)
         INW2=ABS(KNW2)
         INX2=ABS(KNX2)
         INY2=ABS(KNY2)
         IF((KNW2.NE.0).AND.(KNW1.NE.0).AND.(INW1.GE.INW2)) THEN
            KEY=(IL/2)*MUMAX+MU(INW1)-INW1+INW2
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNW2))
            SYS(KEY)=SYS(KEY)-SG*FACT*COEF*GARS*R(K2,K1)
         ENDIF
         IF((KNX2.NE.0).AND.(KNX1.NE.0).AND.(INX1.GE.INX2)) THEN
            KEY=(IL/2)*MUMAX+MU(INX1)-INX1+INX2
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNX2))
            SYS(KEY)=SYS(KEY)-SG*FACT*COEF*GARS*R(K2,K1)
         ENDIF
         IF((KNY2.NE.0).AND.(KNY1.NE.0).AND.(INY1.GE.INY2)) THEN
            KEY=(IL/2)*MUMAX+MU(INY1)-INY1+INY2
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNY2))
            SYS(KEY)=SYS(KEY)-SG*FACT*COEF*GARS*R(K2,K1)
         ENDIF
   40    CONTINUE
         IF(ITY.EQ.0) THEN
            IF(KNW1.NE.0) THEN
               KEY=(IL/2)*MUMAX+MU(INW1)
               IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
                  SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM,1)*ZMARS
               ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
                  SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM,2)*ZMARS
               ENDIF
            ENDIF
            IF(KNX1.NE.0) THEN
               KEY=(IL/2)*MUMAX+MU(INX1)
               IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
                  SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM,3)*ZMARS
               ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
                  SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM,4)*ZMARS
               ENDIF
            ENDIF
            IF(KNY1.NE.0) THEN
               KEY=(IL/2)*MUMAX+MU(INY1)
               IF((K2.EQ.1).AND.(K4.EQ.0)) THEN
                  SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM,5)*ZMARS
               ELSE IF((K2.EQ.IELEM+1).AND.(K4.EQ.1)) THEN
                  SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM,6)*ZMARS
               ENDIF
            ENDIF
         ENDIF
   50    CONTINUE
   51    CONTINUE
   52    CONTINUE
*
         ITRS=0
         DO I=1,NBLOS
            IF(KN(I,1).EQ.KN(NUM,4)) THEN
               ITRS=I
               GO TO 60
            ENDIF
         ENDDO
         CALL XABORT('PNDH2E: ITRS FAILURE.')
   60    DO 75 I=1,NELEM
         KNW1=KN(ITRS,4+I)
         KNX1=KN(NUM,4+2*NELEM+I)
         KNY1=KN(NUM,4+4*NELEM+I)
         INW1=ABS(KNW1)
         INX1=ABS(KNX1)
         INY1=ABS(KNY1)
         DO 70 J=1,NELEM
         KNW2=KN(NUM,4+NELEM+J)
         KNX2=KN(NUM,4+3*NELEM+J)
         KNY2=KN(NUM,4+5*NELEM+J)
         INW2=ABS(KNW2)
         INX2=ABS(KNX2)
         INY2=ABS(KNY2)
         VAR1=FACT*GARS*CTRAN(I,J)
         IF((KNY2.NE.0).AND.(KNW1.NE.0).AND.(INW1.LT.INY2)) THEN
            KEY=(IL/2)*MUMAX+MU(INY2)-INY2+INW1
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
            SYS(KEY)=SYS(KEY)-SG*REAL(VAR1) ! y w
         ELSE IF((KNY2.NE.0).AND.(KNW1.NE.0).AND.(INW1.GT.INY2)) THEN
            KEY=(IL/2)*MUMAX+MU(INW1)-INW1+INY2
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
            SYS(KEY)=SYS(KEY)-SG*REAL(VAR1) ! w y
         ENDIF
         IF((KNW2.NE.0).AND.(KNX1.NE.0).AND.(INW2.LT.INX1)) THEN
            KEY=(IL/2)*MUMAX+MU(INX1)-INX1+INW2
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
            SYS(KEY)=SYS(KEY)-SG*REAL(VAR1) ! x w
         ELSE IF((KNW2.NE.0).AND.(KNX1.NE.0).AND.(INW2.GT.INX1)) THEN
            KEY=(IL/2)*MUMAX+MU(INW2)-INW2+INX1
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
            SYS(KEY)=SYS(KEY)-SG*REAL(VAR1) ! w x
         ENDIF
         IF((KNX2.NE.0).AND.(KNY1.NE.0).AND.(INX2.LT.INY1)) THEN
            KEY=(IL/2)*MUMAX+MU(INY1)-INY1+INX2
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
            SYS(KEY)=SYS(KEY)-SG*REAL(VAR1) ! y x
         ELSE IF((KNX2.NE.0).AND.(KNY1.NE.0).AND.(INX2.GT.INY1)) THEN
            KEY=(IL/2)*MUMAX+MU(INX2)-INX2+INY1
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
            SYS(KEY)=SYS(KEY)-SG*REAL(VAR1) ! x y
         ENDIF
   70    CONTINUE
   75    CONTINUE
*
         IF(ITY.EQ.0) THEN
            DO 83 K4=0,1
            DO 82 K3=0,IELEM-1
            DO 81 K2=1,IELEM+1
            KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
            KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
            KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
            INW1=ABS(KNW1)
            INX1=ABS(KNX1)
            INY1=ABS(KNY1)
            DO 80 K1=0,IELEM-1
            IF(V(K2,K1+1).EQ.0.0) GO TO 80
            IF(K4.EQ.0) THEN
               SSS=(-1.0)**K1
               JND1=KN(NUM,1)+K3*IELEM+K1
               JND2=KN(NUM,2)+K3*IELEM+K1
               JND3=KN(NUM,3)+K3*IELEM+K1
            ELSE
               SSS=1.0
               JND1=KN(NUM,2)+K1*IELEM+K3
               JND2=KN(NUM,3)+K1*IELEM+K3
               JND3=KN(NUM,4)+K1*IELEM+K3
            ENDIF
            IF(KNW1.NE.0) THEN
               IF(JND1.GT.INW1) KEY=(IL/2)*MUMAX+MU(JND1)-JND1+INW1
               IF(JND1.LT.INW1) KEY=(IL/2)*MUMAX+MU(INW1)-INW1+JND1
               SG=REAL(SIGN(1,KNW1))
               SYS(KEY)=SYS(KEY)+SG*SSS*REAL(IL)*SIDE*V(K2,K1+1)
            ENDIF
            IF(KNX1.NE.0) THEN
               IF(JND2.GT.INX1) KEY=(IL/2)*MUMAX+MU(JND2)-JND2+INX1
               IF(JND2.LT.INX1) KEY=(IL/2)*MUMAX+MU(INX1)-INX1+JND2
               SG=REAL(SIGN(1,KNX1))
               SYS(KEY)=SYS(KEY)+SG*SSS*REAL(IL)*SIDE*V(K2,K1+1)
            ENDIF
            IF(KNY1.NE.0) THEN
               IF(JND3.GT.INY1) KEY=(IL/2)*MUMAX+MU(JND3)-JND3+INY1
               IF(JND3.LT.INY1) KEY=(IL/2)*MUMAX+MU(INY1)-INY1+JND3
               SG=REAL(SIGN(1,KNY1))
               SYS(KEY)=SYS(KEY)+SG*SSS*REAL(IL)*SIDE*V(K2,K1+1)
            ENDIF
   80       CONTINUE
   81       CONTINUE
   82       CONTINUE
   83       CONTINUE
         ENDIF
      ENDIF
   90 CONTINUE
  100 CONTINUE
      RETURN
      END
