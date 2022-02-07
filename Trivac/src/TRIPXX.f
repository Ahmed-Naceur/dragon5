*DECK TRIPXX
      SUBROUTINE TRIPXX (IR,MAXKN,NEL,LL4,VOL,MAT,XSGD,XX,YY,ZZ,DD,KN,
     1 QFR,MUX,IPX,CYLIND,LC,T,TS,Q,QS,A11X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for a primal finite element method in 3-D.
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
* IR      first dimension for matrix SGD.
* MAXKN   first dimension for matrix KN.
* NEL     total number of finite elements.
* LL4     order of system matrices.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      values used with a cylindrical geometry.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* XSGD    nuclear properties, derivatives or first variations of
*         nuclear properties per material mixture:
*         XSGD(L,1): X-oriented diffusion coefficients;
*         XSGD(L,2): Y-oriented diffusion coefficients;
*         XSGD(L,3): Z-oriented diffusion coefficients;
*         XSGD(L,4): removal macroscopic cross section.
* MUX     X-oriented compressed storage mode indices.
* MUY     Y-oriented compressed storage mode indices.
* MUZ     Z-oriented compressed storage mode indices.
* IPX     X-oriented permutation matrices.
* IPY     Y-oriented permutation matrices.
* IPZ     Z-oriented permutation matrices.
* CYLIND  cylinderization flag (=.true. for cylindrical geometry).
* LC      order of the unit matrices.
* T       cartesian linear product vector.
* TS      cylindrical linear product vector.
* Q       cartesian stiffness matrix.
* QS      cylindrical stiffness matrix.
*
*Parameters: output
* A11X    X-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
* A11Y    Y-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
* A11Z    Z-oriented matrix corresponding to the divergence (i.e
*         leakage) and removal terms (should be initialized by the
*         calling program).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,MAXKN,NEL,LL4,MAT(NEL),KN(MAXKN),MUX(LL4),IPX(LL4),LC
      REAL VOL(NEL),XSGD(IR,4),XX(NEL),YY(NEL),ZZ(NEL),DD(NEL),
     1 QFR(6*NEL),T(LC),TS(LC),Q(LC,LC),QS(LC,LC),A11X(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION VAR1,VOL1,VOL2,VOL3,QQX,QQY,QQZ
      COMMON /ELEM2/LL,LCC,IJ1(125),IJ2(125),IJ3(125),ISR(6,25),
     1 Q3DP1(125,125),Q3DP2(125,125),Q3DP3(125,125),R3DP(125),
     2 Q3DC1(125,125),Q3DC2(125,125),Q3DC3(125,125),R3DC(125),
     3 R2DP(25),R2DC(25)
*----
*  X-DIRECTED COUPLINGS.
*
*  ASSEMBLY OF MATRIX A11X.
*----
      CALL TRIPMA(LC,T,TS,Q,QS)
      NUM1=0
      NUM2=0
      DO 90 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 90
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 80
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
      VOL1=VOL0/(DX*DX)
      VOL2=VOL0/(DY*DY)
      VOL3=VOL0/(DZ*DZ)
      DO 50 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 50
      INX1=IPX(IND1)
      KEY0=MUX(INX1)
      IF(CYLIND) THEN
         RR=(R3DP(I)+R3DC(I)*DX/DD(K))*VOL0
      ELSE
         RR=R3DP(I)*VOL0
      ENDIF
      A11X(KEY0)=A11X(KEY0)+RR*XSGD(L,4)
      KEY0=KEY0-INX1
      DO 40 J=1,LL
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 40
      INX2=IPX(IND2)
      IF(INX2.EQ.INX1) THEN
         IF(CYLIND) THEN
            QQX=(Q3DP1(I,J)+Q3DC1(I,J)*DX/DD(K))*VOL1
            QQY=(Q3DP2(I,J)+Q3DC2(I,J)*DX/DD(K))*VOL2
            QQZ=(Q3DP3(I,J)+Q3DC3(I,J)*DX/DD(K))*VOL3
         ELSE
            QQX=Q3DP1(I,J)*VOL1
            QQY=Q3DP2(I,J)*VOL2
            QQZ=Q3DP3(I,J)*VOL3
         ENDIF
         KEY=KEY0+INX2
         VAR1=QQX*XSGD(L,1)+QQY*XSGD(L,2)+QQZ*XSGD(L,3)
         A11X(KEY)=REAL(A11X(KEY)+VAR1)
      ELSE IF((INX2.LT.INX1).AND.(IJ2(I).EQ.IJ2(J)).AND.
     1 (IJ3(I).EQ.IJ3(J))) THEN
         IF(CYLIND) THEN
            QQX=(Q3DP1(I,J)+Q3DC1(I,J)*DX/DD(K))*VOL1
         ELSE
            QQX=Q3DP1(I,J)*VOL1
         ENDIF
         KEY=KEY0+INX2
         A11X(KEY)=REAL(A11X(KEY)+QQX*XSGD(L,1))
      ENDIF
   40 CONTINUE
   50 CONTINUE
      DO 70 IC=1,6
      QFR1=QFR(NUM2+IC)
      IF(QFR1.EQ.0.0) GO TO 70
      DO 60 I1=1,LCC
      I=ISR(IC,I1)
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 60
      INX1=IPX(IND1)
      KEY=MUX(INX1)
      IF(CYLIND) THEN
         IF(IC.EQ.1) THEN
            CRZ=-0.5*R2DP(I1)
         ELSE IF(IC.EQ.2) THEN
            CRZ=0.5*R2DP(I1)
         ELSE
            CRZ=R2DC(I1)
         ENDIF
         RR=(R2DP(I1)+CRZ*DX/DD(K))
      ELSE
         RR=R2DP(I1)
      ENDIF
      A11X(KEY)=A11X(KEY)+RR*QFR1
   60 CONTINUE
   70 CONTINUE
   80 NUM1=NUM1+LL
      NUM2=NUM2+6
   90 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIPXY (IR,MAXKN,NEL,LL4,VOL,MAT,XSGD,XX,YY,ZZ,DD,KN,
     1 QFR,MUY,IPY,CYLIND,LC,T,TS,Q,QS,A11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,MAXKN,NEL,LL4,MAT(NEL),KN(MAXKN),MUY(LL4),IPY(LL4),LC
      REAL VOL(NEL),XSGD(IR,4),XX(NEL),YY(NEL),ZZ(NEL),DD(NEL),
     1 QFR(6*NEL),T(LC),TS(LC),Q(LC,LC),QS(LC,LC),A11Y(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION VAR1,VOL1,VOL2,VOL3,QQX,QQY,QQZ
      COMMON /ELEM2/LL,LCC,IJ1(125),IJ2(125),IJ3(125),ISR(6,25),
     1 Q3DP1(125,125),Q3DP2(125,125),Q3DP3(125,125),R3DP(125),
     2 Q3DC1(125,125),Q3DC2(125,125),Q3DC3(125,125),R3DC(125),
     3 R2DP(25),R2DC(25)
*----
*  Y-DIRECTED COUPLINGS.
*
*  ASSEMBLY OF MATRIX A11Y.
*----
      CALL TRIPMA(LC,T,TS,Q,QS)
      NUM1=0
      NUM2=0
      DO 180 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 180
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 170
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
      VOL1=VOL0/(DX*DX)
      VOL2=VOL0/(DY*DY)
      VOL3=VOL0/(DZ*DZ)
      DO 140 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 140
      INY1=IPY(IND1)
      KEY0=MUY(INY1)
      IF(CYLIND) THEN
         RR=(R3DP(I)+R3DC(I)*DX/DD(K))*VOL0
      ELSE
         RR=R3DP(I)*VOL0
      ENDIF
      A11Y(KEY0)=A11Y(KEY0)+RR*XSGD(L,4)
      KEY0=KEY0-INY1
      DO 130 J=1,LL
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 130
      INY2=IPY(IND2)
      IF(INY2.EQ.INY1) THEN
         IF(CYLIND) THEN
            QQX=(Q3DP1(I,J)+Q3DC1(I,J)*DX/DD(K))*VOL1
            QQY=(Q3DP2(I,J)+Q3DC2(I,J)*DX/DD(K))*VOL2
            QQZ=(Q3DP3(I,J)+Q3DC3(I,J)*DX/DD(K))*VOL3
         ELSE
            QQX=Q3DP1(I,J)*VOL1
            QQY=Q3DP2(I,J)*VOL2
            QQZ=Q3DP3(I,J)*VOL3
         ENDIF
         KEY=KEY0+INY2
         VAR1=QQX*XSGD(L,1)+QQY*XSGD(L,2)+QQZ*XSGD(L,3)
         A11Y(KEY)=REAL(A11Y(KEY)+VAR1)
      ELSE IF((INY2.LT.INY1).AND.(IJ1(I).EQ.IJ1(J)).AND.
     1 (IJ3(I).EQ.IJ3(J))) THEN
         IF(CYLIND) THEN
            QQY=(Q3DP2(I,J)+Q3DC2(I,J)*DX/DD(K))*VOL2
         ELSE
            QQY=Q3DP2(I,J)*VOL2
         ENDIF
         KEY=KEY0+INY2
         A11Y(KEY)=REAL(A11Y(KEY)+QQY*XSGD(L,2))
      ENDIF
  130 CONTINUE
  140 CONTINUE
      DO 160 IC=1,6
      QFR1=QFR(NUM2+IC)
      IF(QFR1.EQ.0.0) GO TO 160
      DO 150 I1=1,LCC
      I=ISR(IC,I1)
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 150
      INY1=IPY(IND1)
      KEY=MUY(INY1)
      IF(CYLIND) THEN
         IF(IC.EQ.1) THEN
            CRZ=-0.5*R2DP(I1)
         ELSE IF(IC.EQ.2) THEN
            CRZ=0.5*R2DP(I1)
         ELSE
            CRZ=R2DC(I1)
         ENDIF
         RR=(R2DP(I1)+DX*CRZ/DD(K))
      ELSE
         RR=R2DP(I1)
      ENDIF
      A11Y(KEY)=A11Y(KEY)+RR*QFR1
  150 CONTINUE
  160 CONTINUE
  170 NUM1=NUM1+LL
      NUM2=NUM2+6
  180 CONTINUE
      RETURN
      END
*
      SUBROUTINE TRIPXZ (IR,MAXKN,NEL,LL4,VOL,MAT,XSGD,XX,YY,ZZ,DD,KN,
     1 QFR,MUZ,IPZ,CYLIND,LC,T,TS,Q,QS,A11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,MAXKN,NEL,LL4,MAT(NEL),KN(MAXKN),MUZ(LL4),IPZ(LL4),LC
      REAL VOL(NEL),XSGD(IR,4),XX(NEL),YY(NEL),ZZ(NEL),DD(NEL),
     1 QFR(6*NEL),T(LC),TS(LC),Q(LC,LC),QS(LC,LC),A11Z(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION VAR1,VOL1,VOL2,VOL3,QQX,QQY,QQZ
      COMMON /ELEM2/LL,LCC,IJ1(125),IJ2(125),IJ3(125),ISR(6,25),
     1 Q3DP1(125,125),Q3DP2(125,125),Q3DP3(125,125),R3DP(125),
     2 Q3DC1(125,125),Q3DC2(125,125),Q3DC3(125,125),R3DC(125),
     3 R2DP(25),R2DC(25)
*----
*  Z-DIRECTED COUPLINGS.
*
*  ASSEMBLY OF MATRIX A11Z.
*----
      CALL TRIPMA(LC,T,TS,Q,QS)
      NUM1=0
      NUM2=0
      DO 270 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 270
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 260
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
      VOL1=VOL0/(DX*DX)
      VOL2=VOL0/(DY*DY)
      VOL3=VOL0/(DZ*DZ)
      DO 230 I=1,LL
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 230
      INZ1=IPZ(IND1)
      KEY0=MUZ(INZ1)
      IF(CYLIND) THEN
         RR=(R3DP(I)+R3DC(I)*DX/DD(K))*VOL0
      ELSE
         RR=R3DP(I)*VOL0
      ENDIF
      A11Z(KEY0)=A11Z(KEY0)+RR*XSGD(L,4)
      KEY0=KEY0-INZ1
      DO 220 J=1,LL
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 220
      INZ2=IPZ(IND2)
      IF(INZ2.EQ.INZ1) THEN
         IF(CYLIND) THEN
            QQX=(Q3DP1(I,J)+Q3DC1(I,J)*DX/DD(K))*VOL1
            QQY=(Q3DP2(I,J)+Q3DC2(I,J)*DX/DD(K))*VOL2
            QQZ=(Q3DP3(I,J)+Q3DC3(I,J)*DX/DD(K))*VOL3
         ELSE
            QQX=Q3DP1(I,J)*VOL1
            QQY=Q3DP2(I,J)*VOL2
            QQZ=Q3DP3(I,J)*VOL3
         ENDIF
         KEY=KEY0+INZ2
         VAR1=QQX*XSGD(L,1)+QQY*XSGD(L,2)+QQZ*XSGD(L,3)
         A11Z(KEY)=REAL(A11Z(KEY)+VAR1)
      ELSE IF((INZ2.LT.INZ1).AND.(IJ1(I).EQ.IJ1(J)).AND.
     1 (IJ2(I).EQ.IJ2(J))) THEN
         IF(CYLIND) THEN
            QQZ=(Q3DP3(I,J)+Q3DC3(I,J)*DX/DD(K))*VOL3
         ELSE
            QQZ=Q3DP3(I,J)*VOL3
         ENDIF
         KEY=KEY0+INZ2
         A11Z(KEY)=REAL(A11Z(KEY)+QQZ*XSGD(L,3))
      ENDIF
  220 CONTINUE
  230 CONTINUE
      DO 250 IC=1,6
      QFR1=QFR(NUM2+IC)
      IF(QFR1.EQ.0.0) GO TO 250
      DO 240 I1=1,LCC
      I=ISR(IC,I1)
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 240
      INZ1=IPZ(IND1)
      KEY=MUZ(INZ1)
      IF(CYLIND) THEN
         IF(IC.EQ.1) THEN
            CRZ=-0.5*R2DP(I1)
         ELSE IF(IC.EQ.2) THEN
            CRZ=0.5*R2DP(I1)
         ELSE
            CRZ=R2DC(I1)
         ENDIF
         RR=(R2DP(I1)+DX*CRZ/DD(K))
      ELSE
         RR=R2DP(I1)
      ENDIF
      A11Z(KEY)=A11Z(KEY)+RR*QFR1
  240 CONTINUE
  250 CONTINUE
  260 NUM1=NUM1+LL
      NUM2=NUM2+6
  270 CONTINUE
      RETURN
      END
