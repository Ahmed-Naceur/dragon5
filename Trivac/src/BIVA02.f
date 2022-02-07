*DECK BIVA02
      SUBROUTINE BIVA02(ITY,SGD,CYLIND,IELEM,ICOL,NREG,LL4,NBMIX,IIMAX,
     1 XX,YY,DD,MAT,KN,QFR,VOL,MU,LC,R,V,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a within-group (leakage and removal) or out-of-group
* system matrix in mixed-dual finite element diffusion approximation
* (Cartesian geometry).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ITY     type of assembly: =0: leakage-removal matrix assembly;
*         =1: cross section matrix assembly.
* SGD     nuclear properties. SGD(:,1) and SGD(:,2) are diffusion
*         coefficients. SGD(:,3) are removal macroscopic cross sections.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry).
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* NREG    number of elements in BIVAC.
* LL4     number of unknowns per group in BIVAC.
* NBMIX   number of macro-mixtures.
* IIMAX   allocated dimension of array SYS.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* DD      value used used with a cylindrical geometry.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* VOL     volume of regions.
* MU      indices used with compressed diagonal storage mode matrix SYS.
* LC      number of polynomials in a complete 1-D basis.
* R       cartesian mass matrix.
* V       nodal coupling matrix.
*
*Parameters: output
* SYS     system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,IELEM,ICOL,NREG,LL4,NBMIX,IIMAX,MAT(NREG),KN(5*NREG),
     1 MU(LL4),LC
      REAL SGD(NBMIX,3),XX(NREG),YY(NREG),DD(NREG),QFR(4*NREG),
     1 VOL(NREG),R(LC,LC),V(LC,LC-1),SYS(IIMAX)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      REAL QQ(5,5)
*
      IF((CYLIND).AND.((IELEM.GT.1).OR.(ICOL.NE.2)))
     1 CALL XABORT('BIVA02: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
*----
*  ASSEMBLY OF A SYSTEM MATRIX.
*----
      IF(ITY.EQ.0) THEN
*        LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY.
         DO 12 I0=1,IELEM
         DO 11 J0=1,IELEM
         QQ(I0,J0)=0.0
         DO 10 K0=2,IELEM
         QQ(I0,J0)=QQ(I0,J0)+V(K0,I0)*V(K0,J0)/R(K0,K0)
   10    CONTINUE
   11    CONTINUE
   12    CONTINUE
         NUM1=0
         NUM2=0
         DO 80 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 80
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 70
         DX=XX(K)
         DY=YY(K)
         IF(CYLIND) THEN
            DIN=1.0-0.5*DX/DD(K)
            DOT=1.0+0.5*DX/DD(K)
         ELSE
            DIN=1.0
            DOT=1.0
         ENDIF
*
         DO 60 I0=1,IELEM
         INX1=ABS(KN(NUM1+2))+I0-1
         INX2=ABS(KN(NUM1+3))+I0-1
         INY1=ABS(KN(NUM1+4))+I0-1
         INY2=ABS(KN(NUM1+5))+I0-1
         DO 50 J0=1,IELEM
         JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
         KEY=MU(JND1)
         SYS(KEY)=SYS(KEY)+VOL0*SGD(L,3)
         DO 20 K0=1,J0
         IF(QQ(J0,K0).EQ.0.0) GO TO 20
         KND1=KN(NUM1+1)+(I0-1)*IELEM+K0-1
         KEY=MU(JND1)-JND1+KND1
         SYS(KEY)=SYS(KEY)+VOL0*QQ(J0,K0)*SGD(L,1)/(DX*DX)
   20    CONTINUE
         IF(KN(NUM1+2).NE.0) THEN
            IF(JND1.GT.INX1) KEY=MU(JND1)-JND1+INX1
            IF(JND1.LT.INX1) KEY=MU(INX1)-INX1+JND1
            SG=REAL(SIGN(1,KN(NUM1+2)))
            SYS(KEY)=SYS(KEY)+SG*(VOL0/DX)*DIN*V(1,J0)
         ENDIF
         IF(KN(NUM1+3).NE.0) THEN
            IF(INX2.GT.JND1) KEY=MU(INX2)-INX2+JND1
            IF(INX2.LT.JND1) KEY=MU(JND1)-JND1+INX2
            SG=REAL(SIGN(1,KN(NUM1+3)))
            SYS(KEY)=SYS(KEY)+SG*(VOL0/DX)*DOT*V(IELEM+1,J0)
         ENDIF
         JND1=KN(NUM1+1)+(J0-1)*IELEM+I0-1
         DO 30 K0=1,J0
         IF(QQ(J0,K0).EQ.0.0) GO TO 30
         KND1=KN(NUM1+1)+(K0-1)*IELEM+I0-1
         KEY=MU(JND1)-JND1+KND1
         SYS(KEY)=SYS(KEY)+VOL0*QQ(J0,K0)*SGD(L,2)/(DY*DY)
   30    CONTINUE
         IF(KN(NUM1+4).NE.0) THEN
            IF(JND1.GT.INY1) KEY=MU(JND1)-JND1+INY1
            IF(JND1.LT.INY1) KEY=MU(INY1)-INY1+JND1
            SG=REAL(SIGN(1,KN(NUM1+4)))
            SYS(KEY)=SYS(KEY)+SG*(VOL0/DY)*V(1,J0)
         ENDIF
         IF(KN(NUM1+5).NE.0) THEN
            IF(INY2.GT.JND1) KEY=MU(INY2)-INY2+JND1
            IF(INY2.LT.JND1) KEY=MU(JND1)-JND1+INY2
            SG=REAL(SIGN(1,KN(NUM1+5)))
            SYS(KEY)=SYS(KEY)+SG*(VOL0/DY)*V(IELEM+1,J0)
         ENDIF
   50    CONTINUE
         IF(KN(NUM1+2).NE.0) THEN
            KEY=MU(INX1)
            SYS(KEY)=SYS(KEY)-DIN*(VOL0*R(1,1)/SGD(L,1)+QFR(NUM2+1))
         ENDIF
         IF(KN(NUM1+3).NE.0) THEN
            KEY=MU(INX2)
            SYS(KEY)=SYS(KEY)-DOT*(VOL0*R(IELEM+1,IELEM+1)/SGD(L,1)
     1      +QFR(NUM2+2))
         ENDIF
         IF(KN(NUM1+4).NE.0) THEN
            KEY=MU(INY1)
            SYS(KEY)=SYS(KEY)-VOL0*R(1,1)/SGD(L,2)-QFR(NUM2+3)
         ENDIF
         IF(KN(NUM1+5).NE.0) THEN
            KEY=MU(INY2)
            SYS(KEY)=SYS(KEY)-VOL0*R(IELEM+1,IELEM+1)/SGD(L,2)
     1      -QFR(NUM2+4)
         ENDIF
         IF(ICOL.NE.2) THEN
            IF((KN(NUM1+2).NE.0).AND.(KN(NUM1+3).NE.0)) THEN
               IF(INX2.GT.INX1) KEY=MU(INX2)-INX2+INX1
               IF(INX2.LE.INX1) KEY=MU(INX1)-INX1+INX2
               SG=REAL(SIGN(1,KN(NUM1+2))*SIGN(1,KN(NUM1+3)))
               IF(INX1.EQ.INX2) SG=2.0*SG
               SYS(KEY)=SYS(KEY)-SG*VOL0*R(IELEM+1,1)/SGD(L,1)
            ENDIF
            IF((KN(NUM1+4).NE.0).AND.(KN(NUM1+5).NE.0)) THEN
               IF(INY2.GT.INY1) KEY=MU(INY2)-INY2+INY1
               IF(INY2.LE.INY1) KEY=MU(INY1)-INY1+INY2
               SG=REAL(SIGN(1,KN(NUM1+4))*SIGN(1,KN(NUM1+5)))
               IF(INY1.EQ.INY2) SG=2.0*SG
               SYS(KEY)=SYS(KEY)-SG*VOL0*R(IELEM+1,1)/SGD(L,2)
            ENDIF
         ENDIF
   60   CONTINUE
   70   NUM1=NUM1+5
        NUM2=NUM2+4
   80   CONTINUE
      ELSE
*        CROSS SECTION SYSTEM MATRIX ASSEMBLY.COMPONENTS WITH 1E-10
*        FACTORS ARE INTRODUCED TO MAKE THE MATRIX INVERTIBLE.
         NUM1=0
         DO 110 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 110
         VOL0=VOL(K)
         IF(VOL0.EQ.0.0) GO TO 100
         DO 95 I0=1,IELEM
         INX1=ABS(KN(NUM1+2))+I0-1
         INX2=ABS(KN(NUM1+3))+I0-1
         INY1=ABS(KN(NUM1+4))+I0-1
         INY2=ABS(KN(NUM1+5))+I0-1
         IF(KN(NUM1+2).NE.0) SYS(MU(INX1))=SYS(MU(INX1))+1.0E-30
         IF(KN(NUM1+3).NE.0) SYS(MU(INX2))=SYS(MU(INX2))+1.0E-30
         IF(KN(NUM1+4).NE.0) SYS(MU(INY1))=SYS(MU(INY1))+1.0E-30
         IF(KN(NUM1+5).NE.0) SYS(MU(INY2))=SYS(MU(INY2))+1.0E-30
         DO 90 J0=1,IELEM
         JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
         KEY=MU(JND1)
         SYS(KEY)=SYS(KEY)+VOL0*SGD(L,1)
   90    CONTINUE
   95    CONTINUE
  100    NUM1=NUM1+5
  110    CONTINUE
      ENDIF
      RETURN
      END
