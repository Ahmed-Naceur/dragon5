*DECK TRIMXX
      SUBROUTINE TRIMXX(IR,CYLIND,IELEM,IDIM,NEL,LL4,VOL,MAT,SGD,XSGD,
     1 XX,YY,ZZ,DD,KN,QFR,MUX,IPX,IPR,A11X)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices for mesh centered finite differences or
* nodal collocation method. Note: system matrices should be initialized
* by the calling program.
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
* IR      first dimension of matrices SGD and XSGD.
* CYLIND  cylindrical geometry flag (set with CYLIND  =.true.).
* IELEM   degree of the polynomial basis: =1 (linear/finite
*         differences); =2 (parabolic); =3 (cubic); =4 (quartic).
* IDIM    number of dimensions (1, 2 or 3).
* NEL     total number of finite elements.
* ll4     order of system matrices.
* VOL     volume of each element.
* MAT     mixture index assigned to each element.
* SGD     nuclear properties by material mixture:
*         SGD(L,1)  X-oriented diffusion coefficients;
*         SGD(L,2)  Y-oriented diffusion coefficients;
*         SGD(L,3)  Z-oriented diffusion coefficients;
*         SGD(L,4)  removal macroscopic cross section.
* XSGD    derivative of nuclear properties if IPR=1;
*         variation of nuclear properties if IPR=2 or IPR=3.
*         Note that XSGD=SGD if IPR=0.
* XX      X-directed mesh spacings.
* YY      Y-directed mesh spacings.
* ZZ      Z-directed mesh spacings.
* DD      used with cylindrical geometry.
* KN      element-ordered unknown list:
*         .GT.0: neighbour index;
*         =-1:   void/albedo boundary condition;
*         =-2:   reflection boundary condition;
*         =-3:   ZERO flux boundary condition;
*         =-4:   SYME boundary condition (axial symmetry).
* QFR     element-ordered boundary conditions.
* MUX     X-directed compressed storage mode indices.
* MUY     Y-directed compressed storage mode indices.
* MUZ     Z-directed compressed storage mode indices.
* IPX     permutation matrices.
* IPY     Y-directed permutation matrices.
* IPZ     Z-directed permutation matrices.
* IPR     type of assembly matrix calculation:
*         =0: compute the system matrices;
*         =1: compute the derivative of system matrices;
*         =2 or =3: compute the variation of system matrices.
*
*Parameters: output
* A11X    X-directed  matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUX(LL4).
* A11Y    Y-directed  matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUY(LL4).
* A11Z    Z-directed  matrices corresponding to the divergence (i.e
*         leakage) and removal terms. Dimensionned to MUZ(LL4).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,IELEM,IDIM,NEL,LL4,MAT(NEL),KN(6*NEL),MUX(LL4),
     1 IPX(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),XX(NEL),YY(NEL),ZZ(NEL),
     1 DD(NEL),QFR(6*NEL),A11X(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL LOGIC
      DOUBLE PRECISION RLL,R,S,QQ,PAIR,A1(6),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  STATEMENT FUNCTION
*----
      IORD(J,K,L,LL,IEL,IW)=(IEL*L+K)*LL*IEL+(1+IEL*(IW-1))+J
*----
*  X-ORIENTED COUPLINGS. ASSEMBLY OF MATRIX A11X
*----
      ALLOCATE(IGAR(NEL))
      LL=0
      DO 10 K=1,NEL
      IF(MAT(K).EQ.0) GO TO 10
      LL=LL+1
      IGAR(K)=LL
   10 CONTINUE
      RLL=REAL(IELEM*(IELEM+1))
      NUM1=0
      DO 70 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 70
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 60
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
*
      IF(IPR.EQ.0) THEN
         CALL TRICO (IELEM,IR,NEL,K,VOL0,MAT,XSGD(1,1),XX,YY,ZZ,DD,
     1   KN(NUM1+1),QFR(NUM1+1),CYLIND,A1)
      ELSE IF(IPR.GE.1) THEN
         CALL TRIDCO (IELEM,IR,NEL,K,VOL0,MAT,SGD(1,1),XSGD(1,1),XX,YY,
     1   ZZ,DD,KN(NUM1+1),QFR(NUM1+1),CYLIND,IPR,A1)
      ENDIF
      KK1=KN(NUM1+1)
      KK2=KN(NUM1+2)
      IF(KK1.EQ.-4) KK1=KK2
      IF(KK2.EQ.-4) KK2=KK1
*
      IF(IELEM.EQ.1) THEN
         IND1=IGAR(K)
         INX1=IPX(IND1)
         KEY0=MUX(INX1)-INX1
         IF(KK1.GT.0) THEN
            INX2=IPX(IGAR(KK1))
            IF(INX2.LT.INX1) THEN
               KEY=KEY0+INX2
               A11X(KEY)=A11X(KEY)-REAL(A1(1))
            ENDIF
         ENDIF
         IF(KK2.GT.0) THEN
            INX2=IPX(IGAR(KK2))
            IF(INX2.LT.INX1) THEN
               KEY=KEY0+INX2
               A11X(KEY)=A11X(KEY)-REAL(A1(2))
            ENDIF
         ENDIF
         KEY=KEY0+INX1
         VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)
         A11X(KEY)=A11X(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
      ELSE
         DO 55 I3=0,IELEM-1
         DO 50 I2=0,IELEM-1
         DO 40 I1=0,IELEM-1
         IND1=IORD(I1,I2,I3,LL,IELEM,IGAR(K))
         INX1=IPX(IND1)
         KEY0=MUX(INX1)-INX1
         QQ=SQRT(REAL(2*I1+1))*(RLL-REAL(I1*(I1+1)))/RLL
         IF(KK1.GT.0) THEN
            PAIR=(-1.0D0)**I1
            DO 20 I0=0,IELEM-1
            LOGIC=(KN((IGAR(KK1)-1)*6+1).NE.-4).OR.(MOD(I0+1,2).NE.0)
            INX2=IPX(IORD(I0,I2,I3,LL,IELEM,IGAR(KK1)))
            IF((INX2.LT.INX1).AND.LOGIC) THEN
               KEY=KEY0+INX2
               R=REAL(I0*(I0+1))
               S=SQRT(REAL(2*I0+1))
               VAR1=0.5D0*QQ*PAIR*S*(RLL-R)*A1(1)
               A11X(KEY)=A11X(KEY)-REAL(VAR1)
            ENDIF
   20       CONTINUE
         ENDIF
         IF(KK2.GT.0) THEN
            DO 25 I0=0,IELEM-1
            INX2=IPX(IORD(I0,I2,I3,LL,IELEM,IGAR(KK2)))
            IF(INX2.LT.INX1) THEN
               PAIR=(-1.0D0)**I0
               IF(KN(NUM1+2).EQ.-4) PAIR=1.0D0
               KEY=KEY0+INX2
               R=REAL(I0*(I0+1))
               S=SQRT(REAL(2*I0+1))
               VAR1=0.5D0*QQ*PAIR*S*(RLL-R)*A1(2)
               A11X(KEY)=A11X(KEY)-REAL(VAR1)
            ENDIF
   25       CONTINUE
         ENDIF
         KEY=KEY0+INX1-I1
         DO 30 I0=0,I1
         R=REAL(I0*(I0+1))
         S=SQRT(REAL(2*I0+1))
         PAIR=1.0D0+(-1.0D0)**(I0+I1)
         VAR1=QQ*(PAIR*S*R*XSGD(L,1)*VOL0/(DX*DX)+0.5D0*S*(RLL-R)*
     1   ((-1.0D0)**(I0+I1)*A1(1)+A1(2)))
         A11X(KEY+I0)=A11X(KEY+I0)+REAL(VAR1)
   30    CONTINUE
*
         KEY=KEY0+INX1
         R=REAL(I2*(I2+1))
         QQ=REAL(2*I2+1)*(RLL-R)/RLL
         VAR1=QQ*(2.0D0*R*XSGD(L,2)*VOL0/(DY*DY)+0.5D0*(RLL-R)*
     1   (A1(3)+A1(4)))
         A11X(KEY)=A11X(KEY)+REAL(VAR1)
*
         R=REAL(I3*(I3+1))
         QQ=REAL(2*I3+1)*(RLL-R)/RLL
         VAR1=QQ*(2.0D0*R*XSGD(L,3)*VOL0/(DZ*DZ)+0.5D0*(RLL-R)*
     1   (A1(5)+A1(6)))+XSGD(L,4)*VOL0
         A11X(KEY)=A11X(KEY)+REAL(VAR1)
*
   40    CONTINUE
         IF((IDIM.EQ.1).AND.(I2.EQ.0)) GO TO 60
         IF((IDIM.EQ.2).AND.(I2.EQ.IELEM-1)) GO TO 60
   50    CONTINUE
   55    CONTINUE
      ENDIF
   60 NUM1=NUM1+6
   70 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
*
      SUBROUTINE TRIMXY(IR,CYLIND,IELEM,IDIM,NEL,LL4,VOL,MAT,SGD,XSGD,
     1 XX,YY,ZZ,DD,KN,QFR,MUY,IPY,IPR,A11Y)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,IELEM,IDIM,NEL,LL4,MAT(NEL),KN(6*NEL),MUY(LL4),
     1 IPY(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),XX(NEL),YY(NEL),ZZ(NEL),
     1 DD(NEL),QFR(6*NEL),A11Y(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL LOGIC
      DOUBLE PRECISION RLL,R,S,QQ,PAIR,A1(6),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  STATEMENT FUNCTION
*----
      IORD(J,K,L,LL,IEL,IW)=(IEL*L+K)*LL*IEL+(1+IEL*(IW-1))+J
*----
*  Y-ORIENTED COUPLINGS. ASSEMBLY OF MATRIX A11Y
*----
      ALLOCATE(IGAR(NEL))
      LL=0
      DO 80 K=1,NEL
      IF(MAT(K).EQ.0) GO TO 80
      LL=LL+1
      IGAR(K)=LL
   80 CONTINUE
      RLL=REAL(IELEM*(IELEM+1))
      NUM1=0
      DO 140 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 140
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 130
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
*
      IF(IPR.EQ.0) THEN
         CALL TRICO (IELEM,IR,NEL,K,VOL0,MAT,XSGD(1,1),XX,YY,ZZ,DD,
     1   KN(NUM1+1),QFR(NUM1+1),CYLIND,A1)
      ELSE IF(IPR.GE.1) THEN
         CALL TRIDCO (IELEM,IR,NEL,K,VOL0,MAT,SGD(1,1),XSGD(1,1),XX,YY,
     1   ZZ,DD,KN(NUM1+1),QFR(NUM1+1),CYLIND,IPR,A1)
      ENDIF
      KK3=KN(NUM1+3)
      KK4=KN(NUM1+4)
      IF(KK3.EQ.-4) KK3=KK4
      IF(KK4.EQ.-4) KK4=KK3
*
      IF(IELEM.EQ.1) THEN
         INY1=IPY(IGAR(K))
         KEY0=MUY(INY1)-INY1
         IF(KK3.GT.0) THEN
            INY2=IPY(IGAR(KK3))
            IF(INY2.LT.INY1) THEN
               KEY=KEY0+INY2
               A11Y(KEY)=A11Y(KEY)-REAL(A1(3))
            ENDIF
         ENDIF
         IF(KK4.GT.0) THEN
            INY2=IPY(IGAR(KK4))
            IF(INY2.LT.INY1) THEN
               KEY=KEY0+INY2
               A11Y(KEY)=A11Y(KEY)-REAL(A1(4))
            ENDIF
         ENDIF
         KEY=KEY0+INY1
         VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)
         A11Y(KEY)=A11Y(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
      ELSE
         DO 125 I3=0,IELEM-1
         DO 120 I2=0,IELEM-1
         DO 110 I1=0,IELEM-1
         INY1=IPY(IORD(I2,I1,I3,LL,IELEM,IGAR(K)))
         KEY0=MUY(INY1)-INY1
         QQ=SQRT(REAL(2*I1+1))*(RLL-REAL(I1*(I1+1)))/RLL
         IF(KK3.GT.0) THEN
            PAIR=(-1.0D0)**I1
            DO 90 I0=0,IELEM-1
            LOGIC=(KN((IGAR(KK3)-1)*6+3).NE.-4).OR.(MOD(I0+1,2).NE.0)
            INY2=IPY(IORD(I2,I0,I3,LL,IELEM,IGAR(KK3)))
            IF((INY2.LT.INY1).AND.LOGIC) THEN
               KEY=KEY0+INY2
               R=REAL(I0*(I0+1))
               S=SQRT(REAL(2*I0+1))
               VAR1=0.5D0*QQ*PAIR*S*(RLL-R)*A1(3)
               A11Y(KEY)=A11Y(KEY)-REAL(VAR1)
            ENDIF
   90       CONTINUE
         ENDIF
         IF(KK4.GT.0) THEN
            DO 95 I0=0,IELEM-1
            INY2=IPY(IORD(I2,I0,I3,LL,IELEM,IGAR(KK4)))
            IF(INY2.LT.INY1) THEN
               PAIR=(-1.0D0)**I0
               IF(KN(NUM1+4).EQ.-4) PAIR=1.0D0
               KEY=KEY0+INY2
               R=REAL(I0*(I0+1))
               S=SQRT(REAL(2*I0+1))
               VAR1=0.5D0*QQ*PAIR*S*(RLL-R)*A1(4)
               A11Y(KEY)=A11Y(KEY)-REAL(VAR1)
            ENDIF
   95       CONTINUE
         ENDIF
         KEY=KEY0+INY1-I1
         DO 100 I0=0,I1
         R=REAL(I0*(I0+1))
         S=SQRT(REAL(2*I0+1))
         PAIR=1.0D0+(-1.0D0)**(I0+I1)
         VAR1=QQ*(PAIR*S*R*XSGD(L,2)*VOL0/(DY*DY)+0.5D0*S*(RLL-R)*
     1   ((-1.0D0)**(I0+I1)*A1(3)+A1(4)))
         A11Y(KEY+I0)=A11Y(KEY+I0)+REAL(VAR1)
  100    CONTINUE
*
         KEY=KEY0+INY1
         R=REAL(I2*(I2+1))
         QQ=REAL(2*I2+1)*(RLL-R)/RLL
         VAR1=QQ*(2.0D0*R*XSGD(L,1)*VOL0/(DX*DX)+0.5D0*(RLL-R)*
     1   (A1(1)+A1(2)))
         A11Y(KEY)=A11Y(KEY)+REAL(VAR1)
*
         R=REAL(I3*(I3+1))
         QQ=REAL(2*I3+1)*(RLL-R)/RLL
         VAR1=QQ*(2.0D0*R*XSGD(L,3)*VOL0/(DZ*DZ)+0.5D0*(RLL-R)*
     1   (A1(5)+A1(6)))+XSGD(L,4)*VOL0
         A11Y(KEY)=A11Y(KEY)+REAL(VAR1)
*
  110    CONTINUE
         IF((IDIM.EQ.2).AND.(I2.EQ.IELEM-1)) GO TO 130
  120    CONTINUE
  125    CONTINUE
      ENDIF
  130 NUM1=NUM1+6
  140 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
*
      SUBROUTINE TRIMXZ(IR,CYLIND,IELEM,NEL,LL4,VOL,MAT,SGD,XSGD,XX,YY,
     1 ZZ,DD,KN,QFR,MUZ,IPZ,IPR,A11Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,IELEM,NEL,LL4,MAT(NEL),KN(6*NEL),MUZ(LL4),IPZ(LL4),IPR
      REAL VOL(NEL),SGD(IR,4),XSGD(IR,4),XX(NEL),YY(NEL),ZZ(NEL),
     1 DD(NEL),QFR(6*NEL),A11Z(*)
      LOGICAL CYLIND
*----
*  LOCAL VARIABLES
*----
      LOGICAL LOGIC
      DOUBLE PRECISION RLL,R,S,QQ,PAIR,A1(6),VAR1
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGAR
*----
*  STATEMENT FUNCTION
*----
      IORD(J,K,L,LL,IEL,IW)=(IEL*L+K)*LL*IEL+(1+IEL*(IW-1))+J
*----
*  Z-ORIENTED COUPLINGS. ASSEMBLY OF MATRIX A11Z
*----
      ALLOCATE(IGAR(NEL))
      LL=0
      DO 150 K=1,NEL
      IF(MAT(K).EQ.0) GO TO 150
      LL=LL+1
      IGAR(K)=LL
  150 CONTINUE
      RLL=REAL(IELEM*(IELEM+1))
      NUM1=0
      DO 210 K=1,NEL
      L=MAT(K)
      IF(L.EQ.0) GO TO 210
      VOL0=VOL(K)
      IF(VOL0.EQ.0.0) GO TO 200
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
*
      IF(IPR.EQ.0) THEN
         CALL TRICO (IELEM,IR,NEL,K,VOL0,MAT,XSGD(1,1),XX,YY,ZZ,DD,
     1   KN(NUM1+1),QFR(NUM1+1),CYLIND,A1)
      ELSE IF(IPR.GE.1) THEN
         CALL TRIDCO (IELEM,IR,NEL,K,VOL0,MAT,SGD(1,1),XSGD(1,1),XX,YY,
     1   ZZ,DD,KN(NUM1+1),QFR(NUM1+1),CYLIND,IPR,A1)
      ENDIF
      KK5=KN(NUM1+5)
      KK6=KN(NUM1+6)
      IF(KK5.EQ.-4) KK5=KK6
      IF(KK6.EQ.-4) KK6=KK5
*
      IF(IELEM.EQ.1) THEN
         INZ1=IPZ(IGAR(K))
         KEY0=MUZ(INZ1)-INZ1
         IF(KK5.GT.0) THEN
            INZ2=IPZ(IGAR(KK5))
            IF(INZ2.LT.INZ1) THEN
               KEY=KEY0+INZ2
               A11Z(KEY)=A11Z(KEY)-REAL(A1(5))
            ENDIF
         ENDIF
         IF(KK6.GT.0) THEN
            INZ2=IPZ(IGAR(KK6))
            IF(INZ2.LT.INZ1) THEN
               KEY=KEY0+INZ2
               A11Z(KEY)=A11Z(KEY)-REAL(A1(6))
            ENDIF
         ENDIF
         KEY=KEY0+INZ1
         VAR1=A1(1)+A1(2)+A1(3)+A1(4)+A1(5)+A1(6)
         A11Z(KEY)=A11Z(KEY)+REAL(VAR1)+XSGD(L,4)*VOL0
      ELSE
         DO 192 I3=0,IELEM-1
         DO 191 I2=0,IELEM-1
         DO 190 I1=0,IELEM-1
         INZ1=IPZ(IORD(I2,I3,I1,LL,IELEM,IGAR(K)))
         KEY0=MUZ(INZ1)-INZ1
         QQ=SQRT(REAL(2*I1+1))*(RLL-REAL(I1*(I1+1)))/RLL
         IF(KK5.GT.0) THEN
            PAIR=(-1.0D0)**I1
            DO 160 I0=0,IELEM-1
            LOGIC=(KN((IGAR(KK5)-1)*6+5).NE.-4).OR.(MOD(I0+1,2).NE.0)
            INZ2=IPZ(IORD(I2,I3,I0,LL,IELEM,IGAR(KK5)))
            IF((INZ2.LT.INZ1).AND.LOGIC) THEN
               KEY=KEY0+INZ2
               R=REAL(I0*(I0+1))
               S=SQRT(REAL(2*I0+1))
               VAR1=0.5D0*QQ*PAIR*S*(RLL-R)*A1(5)
               A11Z(KEY)=A11Z(KEY)-REAL(VAR1)
            ENDIF
  160       CONTINUE
         ENDIF
         IF(KK6.GT.0) THEN
            DO 165 I0=0,IELEM-1
            INZ2=IPZ(IORD(I2,I3,I0,LL,IELEM,IGAR(KK6)))
            IF(INZ2.LT.INZ1) THEN
               PAIR=(-1.0D0)**I0
               IF(KN(NUM1+6).EQ.-4) PAIR=1.0D0
               KEY=KEY0+INZ2
               R=REAL(I0*(I0+1))
               S=SQRT(REAL(2*I0+1))
               VAR1=0.5D0*QQ*PAIR*S*(RLL-R)*A1(6)
               A11Z(KEY)=A11Z(KEY)-REAL(VAR1)
            ENDIF
  165       CONTINUE
         ENDIF
         KEY=KEY0+INZ1-I1
         DO 170 I0=0,I1
         R=REAL(I0*(I0+1))
         S=SQRT(REAL(2*I0+1))
         PAIR=1.0D0+(-1.0D0)**(I0+I1)
         VAR1=QQ*(PAIR*S*R*XSGD(L,3)*VOL0/(DZ*DZ)+0.5D0*S*(RLL-R)*
     1   ((-1.0D0)**(I0+I1)*A1(5)+A1(6)))
         A11Z(KEY+I0)=A11Z(KEY+I0)+REAL(VAR1)
  170    CONTINUE
*
         KEY=KEY0+INZ1
         R=REAL(I2*(I2+1))
         QQ=REAL(2*I2+1)*(RLL-R)/RLL
         VAR1=QQ*(2.0D0*R*XSGD(L,1)*VOL0/(DX*DX)+0.5D0*(RLL-R)*
     1   (A1(1)+A1(2)))
         A11Z(KEY)=A11Z(KEY)+REAL(VAR1)
*
         R=REAL(I3*(I3+1))
         QQ=REAL(2*I3+1)*(RLL-R)/RLL
         VAR1=QQ*(2.0D0*R*XSGD(L,2)*VOL0/(DY*DY)+0.5D0*(RLL-R)*
     1   (A1(3)+A1(4)))+XSGD(L,4)*VOL0
         A11Z(KEY)=A11Z(KEY)+REAL(VAR1)
*
  190    CONTINUE
  191    CONTINUE
  192    CONTINUE
      ENDIF
  200 NUM1=NUM1+6
  210 CONTINUE
      DEALLOCATE(IGAR)
      RETURN
      END
