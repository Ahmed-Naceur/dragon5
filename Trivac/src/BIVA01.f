*DECK BIVA01
      SUBROUTINE BIVA01(ITY,MAXKN,SGD,CYLIND,NREG,LL4,NBMIX,IIMAX,XX,
     1 YY,DD,MAT,KN,QFR,VOL,MU,LC,R,RS,Q,QS,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a within-group (leakage and removal) or out-of-group
* system matrix in mesh corner finite difference or finite element
* diffusion approximation (Cartesian geometry).
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
* MAXKN   dimension of array KN.
* SGD     nuclear properties. SGD(:,1) and SGD(:,2) are diffusion
*         coefficients. SGD(:,3) are removal macroscopic cross sections.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry).
* NREG    number of elements in BIVAC.
* LL4     order of matrix SYS.
* NBMIX   number of macro-mixtures.
* IIMAX   allocated dimension of array SYS.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* DD      value used with a cylindrical geometry.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* VOL     volume of regions.
* MU      indices used with compressed diagonal storage mode matrix SYS.
* LC      number of polynomials in a complete 1-D basis.
* R       Cartesian mass matrix.
* RS      cylindrical mass matrix.
* Q       Cartesian stiffness matrix.
* QS      cylindrical stiffness matrix.
*
*Parameters: output
* SYS     system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,MAXKN,NREG,LL4,NBMIX,IIMAX,MAT(NREG),KN(MAXKN),
     1 MU(LL4),LC
      REAL SGD(NBMIX,3),XX(NREG),YY(NREG),DD(NREG),QFR(4*NREG),
     1 VOL(NREG),R(LC,LC),RS(LC,LC),Q(LC,LC),QS(LC,LC),SYS(IIMAX)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(25),IJ2(25),ISR(4,5)
      REAL Q2DP1(25,25),Q2DP2(25,25),R2DP(25,25),Q2DC1(25,25),
     1 Q2DC2(25,25),R2DC(25,25)
*----
*  COMPUTE VECTORS IJ1, IJ2 AND MATRIX ISR.
*----
      LL=LC*LC
      DO 10 I=1,LL
      IJ1(I)=1+MOD(I-1,LC)
      IJ2(I)=1+(I-IJ1(I))/LC
   10 CONTINUE
      DO 20 I=1,LC
      ISR(1,I)=(I-1)*LC+1
      ISR(2,I)=I*LC
      ISR(3,I)=I
      ISR(4,I)=LL-LC+I
   20 CONTINUE
*----
*  COMPUTE THE CARTESIAN 2-D MASS AND STIFFNESS MATRICES FROM TENSORIAL
*  PRODUCTS OF 1-D MATRICES.
*----
      DO 40 I=1,LL
      I1=IJ1(I)
      I2=IJ2(I)
      DO 30 J=1,LL
      J1=IJ1(J)
      J2=IJ2(J)
      Q2DP1(I,J)=Q(I1,J1)*R(I2,J2)
      Q2DP2(I,J)=R(I1,J1)*Q(I2,J2)
      R2DP(I,J)=R(I1,J1)*R(I2,J2)
      Q2DC1(I,J)=QS(I1,J1)*R(I2,J2)
      Q2DC2(I,J)=RS(I1,J1)*Q(I2,J2)
      R2DC(I,J)=RS(I1,J1)*R(I2,J2)
   30 CONTINUE
   40 CONTINUE
*----
*  ASSEMBLY OF A SYSTEM MATRIX.
*----
      IF(ITY.EQ.0) THEN
*        LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY.
         NUM1=0
         NUM2=0
         DO 110 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 110
         IF(VOL(K).EQ.0.0) GO TO 100
         DX=XX(K)
         DY=YY(K)
         DO 60 I=1,LL
         IND1=KN(NUM1+I)
         IF(IND1.EQ.0) GO TO 60
         KEY1=MU(IND1)-IND1
         DO 50 J=1,LL
         IND2=KN(NUM1+J)
         IF((IND2.EQ.0).OR.(IND2.GT.IND1)) GO TO 50
         IF(CYLIND) THEN
            QQX=(Q2DP1(I,J)+Q2DC1(I,J)*DX/DD(K))/(DX*DX)
            QQY=(Q2DP2(I,J)+Q2DC2(I,J)*DX/DD(K))/(DY*DY)
            RR=R2DP(I,J)+R2DC(I,J)*DX/DD(K)
         ELSE
            QQX=Q2DP1(I,J)/(DX*DX)
            QQY=Q2DP2(I,J)/(DY*DY)
            RR=R2DP(I,J)
         ENDIF
         IF((QQX.EQ.0.0).AND.(QQY.EQ.0.0).AND.(RR.EQ.0.0)) GO TO 50
         KEY=KEY1+IND2
         SYS(KEY)=SYS(KEY)+(QQX*SGD(L,1)+QQY*SGD(L,2)+RR*SGD(L,3))
     1   *VOL(K)
   50    CONTINUE
   60    CONTINUE
         DO 90 IC=1,4
         QFR1=QFR(NUM2+IC)
         IF(QFR1.EQ.0.0) GO TO 90
         DO 80 I1=1,LC
         IND1=KN(NUM1+ISR(IC,I1))
         IF(IND1.EQ.0) GO TO 80
         KEY1=MU(IND1)-IND1
         DO 70 J1=1,LC
         IND2=KN(NUM1+ISR(IC,J1))
         IF((IND2.EQ.0).OR.(IND2.GT.IND1)) GO TO 70
         IF(CYLIND) THEN
            CRZ=0.0
            IF(IC.EQ.1) THEN
               CRZ=-0.5*R(I1,J1)
            ELSE IF(IC.EQ.2) THEN
               CRZ=0.5*R(I1,J1)
            ELSE IF(IC.EQ.3) THEN
               CRZ=RS(I1,J1)
            ELSE IF(IC.EQ.4) THEN
               CRZ=RS(I1,J1)
            ENDIF
            RR=R(I1,J1)+CRZ*DX/DD(K)
         ELSE
            RR=R(I1,J1)
         ENDIF
         IF(RR.EQ.0.0) GO TO 70
         KEY=KEY1+IND2
         SYS(KEY)=SYS(KEY)+RR*QFR1
   70    CONTINUE
   80    CONTINUE
   90    CONTINUE
  100    NUM1=NUM1+LL
         NUM2=NUM2+4
  110    CONTINUE
      ELSE
*        CROSS SECTION SYSTEM MATRIX ASSEMBLY.
         NUM1=0
         DO 150 K=1,NREG
         L=MAT(K)
         IF(L.EQ.0) GO TO 150
         IF(VOL(K).EQ.0.0) GO TO 140
         DX=XX(K)
         DO 130 I=1,LL
         IND1=KN(NUM1+I)
         IF(IND1.EQ.0) GO TO 130
         KEY1=MU(IND1)-IND1
         DO 120 J=1,LL
         IND2=KN(NUM1+J)
         IF((IND2.EQ.0).OR.(IND2.GT.IND1)) GO TO 120
         IF(CYLIND) THEN
            RR=R2DP(I,J)+R2DC(I,J)*DX/DD(K)
         ELSE
            RR=R2DP(I,J)
         ENDIF
         IF(RR.EQ.0.0) GO TO 120
         KEY=KEY1+IND2
         SYS(KEY)=SYS(KEY)+RR*SGD(L,1)*VOL(K)
  120    CONTINUE
  130    CONTINUE
  140    NUM1=NUM1+LL
  150    CONTINUE
      ENDIF
      RETURN
      END
