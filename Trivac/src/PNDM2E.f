*DECK PNDM2E
      SUBROUTINE PNDM2E(ITY,NEL,L4,IELEM,ICOL,MAT,VOL,NBMIX,NLF,NVD,
     1 NAN,SIGT,SIGTI,XX,YY,KN,QFR,MU,IIMAX,LC,R,V,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a mixed-dual formulation of the
* simplified PN method in 2D Cartesian geometry.
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
* NEL     number of finite elements.
* L4      number of unknowns per energy group and per set of two
*         Legendre orders.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* ICOL    type of quadrature: =1 (analytical integration);
*         =2 (Gauss-Lobatto); =3 (Gauss-Legendre).
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* NBMIX   number of mixtures.
* NLF     number of Legendre orders for the flux (even number).
* NVD     type of void boundary condition if NLF>0 and ICOL=3.
* NAN     number of Legendre orders for the cross sections.
* SIGT    total minus self-scattering macroscopic cross sections.
*         SIGT(:,NAN) generally contains the total cross section only.
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* MU      compressed storage mode indices.
* IIMAX   dimension of vector SYS.
* LC      order of the unit matrices.
* R       unit Cartesian mass matrix.
* V       unit nodal coupling matrix.
*
*Parameters: output
* SYS     system matrix.
*
*Reference:
* J.J. Lautard, D. Schneider, A.M. Baudron, "Mixed Dual Methods for
* Neutronic Reactor Core Calculations in the CRONOS System,"
* Proc. Int. Conf. on Mathematics and Computation, Reactor
* Physics and Environmental Analysis in Nuclear Applications,
* Madrid, Spain, September 27-30, 1999.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,NEL,L4,IELEM,ICOL,MAT(NEL),NBMIX,NLF,NAN,KN(5*NEL),
     1 MU(L4),IIMAX,LC
      REAL VOL(NEL),SIGT(NBMIX,NAN),SIGTI(NBMIX,NAN),XX(NEL),YY(NEL),
     1 QFR(4*NEL),R(LC,LC),V(LC,LC-1),SYS(IIMAX)
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
      DO 100 IL=0,NLF-1
      ZMARS=0.0
      IF(MOD(IL,2).EQ.1) ZMARS=PNMAR2(NZMAR,IL,IL)
      FACT=REAL(2*IL+1)
*----
*  ASSEMBLY OF THE MAIN COEFFICIENT MATRIX AT ORDER IL.
*----
      NUM1=0
      NUM2=0
      KEY=0
      DO 90 K=1,NEL
      IBM=MAT(K)
      IF(IBM.EQ.0) GO TO 90
      VOL0=VOL(K)
      GARS=SIGT(IBM,MIN(IL+1,NAN))
      IF(MOD(IL,2).EQ.0) THEN
*        EVEN PARITY EQUATION.
         DO 25 I0=1,IELEM
         DO 20 J0=1,IELEM
         JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
         KEY=(IL/2)*MUMAX+MU(JND1)
         SYS(KEY)=SYS(KEY)+FACT*VOL0*GARS
   20    CONTINUE
   25    CONTINUE
      ELSE
         GARSI=SIGTI(IBM,MIN(IL+1,NAN))
         DO 80 I0=1,IELEM
*        PARTIAL INVERSION OF THE ODD PARITY EQUATION. MODIFICATION OF
*        THE EVEN PARITY EQUATION.
         DO 45 J0=1,IELEM
         JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
         DO 30 K0=1,J0
         IF(QQ(J0,K0).EQ.0.0) GO TO 30
         KND1=KN(NUM1+1)+(I0-1)*IELEM+K0-1
         KEY=(IL/2)*MUMAX+MU(JND1)-JND1+KND1
         SYS(KEY)=SYS(KEY)+(REAL(IL)**2)*VOL0*QQ(J0,K0)*GARSI/(FACT*
     1   XX(K)*XX(K))
         IF(IL.LE.NLF-3) THEN
            KEY=((IL+2)/2)*MUMAX+MU(JND1)-JND1+KND1
            SYS(KEY)=SYS(KEY)+(REAL(IL+1)**2)*VOL0*QQ(J0,K0)*GARSI/
     1      (FACT*XX(K)*XX(K))
         ENDIF
   30    CONTINUE
         JND1=KN(NUM1+1)+(J0-1)*IELEM+I0-1
         DO 40 K0=1,J0
         IF(QQ(J0,K0).EQ.0.0) GO TO 40
         KND1=KN(NUM1+1)+(K0-1)*IELEM+I0-1
         KEY=(IL/2)*MUMAX+MU(JND1)-JND1+KND1
         SYS(KEY)=SYS(KEY)+(REAL(IL)**2)*VOL0*QQ(J0,K0)*GARSI/(FACT*
     1   YY(K)*YY(K))
         IF(IL.LE.NLF-3) THEN
            KEY=((IL+2)/2)*MUMAX+MU(JND1)-JND1+KND1
            SYS(KEY)=SYS(KEY)+(REAL(IL+1)**2)*VOL0*QQ(J0,K0)*GARSI/
     1      (FACT*YY(K)*YY(K))
         ENDIF
   40    CONTINUE
   45    CONTINUE
*
*        ODD PARITY EQUATION.
         DO 55 IC=1,2
         IIC=1
         IF(IC.EQ.2) IIC=IELEM+1
         IND1=ABS(KN(NUM1+1+IC))+I0-1
         S1=REAL(SIGN(1,KN(NUM1+1+IC)))
         DO 50 JC=1,2
         JJC=1
         IF(JC.EQ.2) JJC=IELEM+1
         IND2=ABS(KN(NUM1+1+JC))+I0-1
         IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0).AND.
     1      (IND1.GE.IND2)) THEN
            S2=REAL(SIGN(1,KN(NUM1+1+JC)))
            KEY=(IL/2)*MUMAX+MU(IND1)-IND1+IND2
            SYS(KEY)=SYS(KEY)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS
         ENDIF
   50    CONTINUE
   55    CONTINUE
         DO 65 IC=3,4
         IIC=1
         IF(IC.EQ.4) IIC=IELEM+1
         IND1=ABS(KN(NUM1+1+IC))+I0-1
         S1=REAL(SIGN(1,KN(NUM1+1+IC)))
         DO 60 JC=3,4
         JJC=1
         IF(JC.EQ.4) JJC=IELEM+1
         IND2=ABS(KN(NUM1+1+JC))+I0-1
         IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0).AND.
     1      (IND1.GE.IND2)) THEN
            S2=REAL(SIGN(1,KN(NUM1+1+JC)))
            KEY=(IL/2)*MUMAX+MU(IND1)-IND1+IND2
            SYS(KEY)=SYS(KEY)-S1*S2*FACT*R(IIC,JJC)*VOL0*GARS
         ENDIF
   60    CONTINUE
   65    CONTINUE
         IF(ITY.EQ.1) GO TO 80
*
         IND1=ABS(KN(NUM1+2))+I0-1
         IND2=ABS(KN(NUM1+3))+I0-1
         IND3=ABS(KN(NUM1+4))+I0-1
         IND4=ABS(KN(NUM1+5))+I0-1
         IF((QFR(NUM2+1).NE.0.0).AND.(KN(NUM1+2).NE.0)) THEN
            KEY=(IL/2)*MUMAX+MU(IND1)
            SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM2+1)*ZMARS
         ENDIF
         IF((QFR(NUM2+2).NE.0.0).AND.(KN(NUM1+3).NE.0)) THEN
            KEY=(IL/2)*MUMAX+MU(IND2)
            SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM2+2)*ZMARS
         ENDIF
         IF((QFR(NUM2+3).NE.0.0).AND.(KN(NUM1+4).NE.0)) THEN
            KEY=(IL/2)*MUMAX+MU(IND3)
            SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM2+3)*ZMARS
         ENDIF
         IF((QFR(NUM2+4).NE.0.0).AND.(KN(NUM1+5).NE.0)) THEN
            KEY=(IL/2)*MUMAX+MU(IND4)
            SYS(KEY)=SYS(KEY)-0.5*FACT*QFR(NUM2+4)*ZMARS
         ENDIF
*
         DO 70 J0=1,IELEM
         JND1=KN(NUM1+1)+(I0-1)*IELEM+J0-1
         IF(KN(NUM1+2).NE.0) THEN
            S1=REAL(SIGN(1,KN(NUM1+2)))
            IF(JND1.GT.IND1) KEY=(IL/2)*MUMAX+MU(JND1)-JND1+IND1
            IF(JND1.LT.IND1) KEY=(IL/2)*MUMAX+MU(IND1)-IND1+JND1
            SYS(KEY)=SYS(KEY)+S1*REAL(IL)*VOL0*V(1,J0)/XX(K)
         ENDIF
         IF(KN(NUM1+3).NE.0) THEN
            S2=REAL(SIGN(1,KN(NUM1+3)))
            IF(JND1.GT.IND2) KEY=(IL/2)*MUMAX+MU(JND1)-JND1+IND2
            IF(JND1.LT.IND2) KEY=(IL/2)*MUMAX+MU(IND2)-IND2+JND1
            SYS(KEY)=SYS(KEY)+S2*REAL(IL)*VOL0*V(IELEM+1,J0)/XX(K)
         ENDIF
         JND1=KN(NUM1+1)+(J0-1)*IELEM+I0-1
         IF(KN(NUM1+4).NE.0) THEN
            S3=REAL(SIGN(1,KN(NUM1+4)))
            IF(JND1.GT.IND3) KEY=(IL/2)*MUMAX+MU(JND1)-JND1+IND3
            IF(JND1.LT.IND3) KEY=(IL/2)*MUMAX+MU(IND3)-IND3+JND1
            SYS(KEY)=SYS(KEY)+S3*REAL(IL)*VOL0*V(1,J0)/YY(K)
         ENDIF
         IF(KN(NUM1+5).NE.0) THEN
            S4=REAL(SIGN(1,KN(NUM1+5)))
            IF(JND1.GT.IND4) KEY=(IL/2)*MUMAX+MU(JND1)-JND1+IND4
            IF(JND1.LT.IND4) KEY=(IL/2)*MUMAX+MU(IND4)-IND4+JND1
            SYS(KEY)=SYS(KEY)+S4*REAL(IL)*VOL0*V(IELEM+1,J0)/YY(K)
         ENDIF
   70    CONTINUE
   80    CONTINUE
      ENDIF
      NUM1=NUM1+5
      NUM2=NUM2+4
   90 CONTINUE
  100 CONTINUE
      RETURN
      END
