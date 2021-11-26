*DECK TRIDCO
      SUBROUTINE TRIDCO (IELEM,IR,NEL,K,VOL0,MAT,DIF,DDF,XX,YY,ZZ,DD,KN,
     1 QFR,CYLIND,IPR,A)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the derivative or variation of mesh centered finite difference
* coefficients in element K.
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
* IELEM   degree of the polynomial basis: =1 (linear/finite
*         differences); =2 (parabolic); =3 (cubic); =4 (quartic).
* IR      first dimension of matrix DIF.
* NEL     total number of finite elements.
* K       index of finite element under consideration.
* VOL0    volume of finite element under consideration.
* MAT     mixture index assigned to each element.
* DIF     directional diffusion coefficients
* DDF     derivative or variation of directional diffusion coefficients.
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
* CYLIND  cylindrical geometry flag (set with CYLIND  =.true.).
* IPR     type of coefficient calculation:
*         .eq.1 take derivative of MCFD coefficients;
*         .ge.2 take variation of MCFD coefficients.
*
*Parameters: output
* A       derivative or variation of mesh centered finite difference
*         coefficients.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,IR,NEL,K,MAT(NEL),KN(6),IPR
      REAL VOL0,DIF(IR,3),DDF(IR,3),XX(NEL),YY(NEL),ZZ(NEL),DD(NEL),
     1 QFR(6)
      LOGICAL CYLIND
      DOUBLE PRECISION A(6)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION VHARM,DHARM,DIN,DOT
*
*     VARIATION FORMULA:
      VHARM(X1,X2,DIF1,DIF2,DDF1,DDF2)=2.0D0*((DIF1+DDF1)*(DIF2+DDF2)
     1 /(X1*(DIF2+DDF2)+X2*(DIF1+DDF1))-DIF1*DIF2/(X1*DIF2+X2*DIF1))
*     DERIVATIVE FORMULA:
      DHARM(X1,X2,DIF1,DIF2,DDF1,DDF2)=2.0D0*(X1*DIF2*DIF2*DDF1+
     1 X2*DIF1*DIF1*DDF2)/(X1*DIF2+X2*DIF1)**2
*
      DENOM=REAL((IELEM+1)*IELEM)
      L=MAT(K)
      DX=XX(K)
      DY=YY(K)
      DZ=ZZ(K)
      IF(CYLIND) THEN
         DIN=1.0D0-0.5D0*DX/DD(K)
         DOT=1.0D0+0.5D0*DX/DD(K)
      ELSE
         DIN=1.0D0
         DOT=1.0D0
      ENDIF
      KK1=KN(1)
      KK2=KN(2)
      KK3=KN(3)
      KK4=KN(4)
      KK5=KN(5)
      KK6=KN(6)
      IF(IPR.EQ.1) THEN
*        DERIVATIVE FORMULA.
*        X- SIDE:
         IF(KK1.GT.0) THEN
            A(1)=DHARM(DX,XX(KK1),DIF(L,1),DIF(MAT(KK1),1),DDF(L,1),
     1      DDF(MAT(KK1),1))*DIN*VOL0/DX
         ELSE IF(KK1.EQ.-1) THEN
            A(1)=DHARM(DX,DX,DIF(L,1),DX*QFR(1)/DENOM,DDF(L,1),0.0)
     1      *DIN*VOL0/DX
         ELSE IF(KK1.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(1)=2.0D0*DDF(L,1)*DIN*VOL0/(DX*DX)
         ENDIF
*        X+ SIDE:
         IF(KK2.GT.0) THEN
            A(2)=DHARM(DX,XX(KK2),DIF(L,1),DIF(MAT(KK2),1),DDF(L,1),
     1      DDF(MAT(KK2),1))*DOT*VOL0/DX
         ELSE IF(KK2.EQ.-1) THEN
            A(2)=DHARM(DX,DX,DIF(L,1),DX*QFR(2)/DENOM,DDF(L,1),0.0)
     1      *DOT*VOL0/DX
         ELSE IF(KK2.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(2)=2.0D0*DDF(L,1)*DOT*VOL0/(DX*DX)
         ELSE IF(KK2.EQ.-4) THEN
            IF(KK1.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (1).')
            A(2)=A(1)
         ENDIF
         IF(KK1.EQ.-4) THEN
            IF(KK2.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (2).')
            A(1)=A(2)
         ENDIF
*        Y- SIDE:
         IF(KK3.GT.0) THEN
            A(3)=DHARM(DY,YY(KK3),DIF(L,2),DIF(MAT(KK3),2),DDF(L,2),
     1      DDF(MAT(KK3),2))*VOL0/DY
         ELSE IF(KK3.EQ.-1) THEN
            A(3)=DHARM(DY,DY,DIF(L,2),DY*QFR(3)/DENOM,DDF(L,2),0.0)
     1      *VOL0/DY
         ELSE IF(KK3.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,2)*VOL0/(DY*DY)
         ENDIF
*        Y+ SIDE:
         IF(KK4.GT.0) THEN
            A(4)=DHARM(DY,YY(KK4),DIF(L,2),DIF(MAT(KK4),2),DDF(L,2),
     1      DDF(MAT(KK4),2))*VOL0/DY
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=DHARM(DY,DY,DIF(L,2),DY*QFR(4)/DENOM,DDF(L,2),0.0)
     1      *VOL0/DY
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*DDF(L,2)*VOL0/(DY*DY)
         ELSE IF(KK4.EQ.-4) THEN
            IF(KK3.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (3).')
            A(4)=A(3)
         ENDIF
         IF(KK3.EQ.-4) THEN
            IF(KK4.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (4).')
            A(3)=A(4)
         ENDIF
*        Z- SIDE:
         IF(KK5.GT.0) THEN
            A(5)=DHARM(DZ,ZZ(KK5),DIF(L,3),DIF(MAT(KK5),3),DDF(L,3),
     1      DDF(MAT(KK5),3))*VOL0/DZ
         ELSE IF(KK5.EQ.-1) THEN
            A(5)=DHARM(DZ,DZ,DIF(L,3),DZ*QFR(5)/DENOM,DDF(L,3),0.0)
     1      *VOL0/DZ
         ELSE IF(KK5.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(5)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*        Z+ SIDE:
         IF(KK6.GT.0) THEN
            A(6)=DHARM(DZ,ZZ(KK6),DIF(L,3),DIF(MAT(KK6),3),DDF(L,3),
     1      DDF(MAT(KK6),3))*VOL0/DZ
         ELSE IF(KK6.EQ.-1) THEN
            A(6)=DHARM(DZ,DZ,DIF(L,3),DZ*QFR(6)/DENOM,DDF(L,3),0.0)
     1      *VOL0/DZ
         ELSE IF(KK6.EQ.-2) THEN
            A(6)=0.0D0
         ELSE IF(KK6.EQ.-3) THEN
            A(6)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ELSE IF(KK6.EQ.-4) THEN
            IF(KK5.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (5).')
            A(6)=A(5)
         ENDIF
         IF(KK5.EQ.-4) THEN
            IF(KK6.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (6).')
            A(5)=A(6)
         ENDIF
      ELSE
*        VARIATION FORMULA.
*        X- SIDE:
         IF(KK1.GT.0) THEN
            A(1)=VHARM(DX,XX(KK1),DIF(L,1),DIF(MAT(KK1),1),DDF(L,1),
     1      DDF(MAT(KK1),1))*DIN*VOL0/DX
         ELSE IF(KK1.EQ.-1) THEN
            A(1)=VHARM(DX,DX,DIF(L,1),DX*QFR(1)/DENOM,DDF(L,1),0.0)
     1      *DIN*VOL0/DX
         ELSE IF(KK1.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(1)=2.0D0*DDF(L,1)*DIN*VOL0/(DX*DX)
         ENDIF
*        X+ SIDE:
         IF(KK2.GT.0) THEN
            A(2)=VHARM(DX,XX(KK2),DIF(L,1),DIF(MAT(KK2),1),DDF(L,1),
     1      DDF(MAT(KK2),1))*DOT*VOL0/DX
         ELSE IF(KK2.EQ.-1) THEN
            A(2)=VHARM(DX,DX,DIF(L,1),DX*QFR(2)/DENOM,DDF(L,1),0.0)
     1      *DOT*VOL0/DX
         ELSE IF(KK2.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(2)=2.0D0*DDF(L,1)*DOT*VOL0/(DX*DX)
         ELSE IF(KK2.EQ.-4) THEN
            IF(KK1.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (7).')
            A(2)=A(1)
         ENDIF
         IF(KK1.EQ.-4) THEN
            IF(KK2.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (8).')
            A(1)=A(2)
         ENDIF
*        Y- SIDE:
         IF(KK3.GT.0) THEN
            A(3)=VHARM(DY,YY(KK3),DIF(L,2),DIF(MAT(KK3),2),DDF(L,2),
     1      DDF(MAT(KK3),2))*VOL0/DY
         ELSE IF(KK3.EQ.-1) THEN
            A(3)=VHARM(DY,DY,DIF(L,2),DY*QFR(3)/DENOM,DDF(L,2),0.0)
     1      *VOL0/DY
         ELSE IF(KK3.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,2)*VOL0/(DY*DY)
         ENDIF
*        Y+ SIDE:
         IF(KK4.GT.0) THEN
            A(4)=VHARM(DY,YY(KK4),DIF(L,2),DIF(MAT(KK4),2),DDF(L,2),
     1      DDF(MAT(KK4),2))*VOL0/DY
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=VHARM(DY,DY,DIF(L,2),DY*QFR(4)/DENOM,DDF(L,2),0.0)
     1      *VOL0/DY
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*DDF(L,2)*VOL0/(DY*DY)
         ELSE IF(KK4.EQ.-4) THEN
            IF(KK3.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (9).')
            A(4)=A(3)
         ENDIF
         IF(KK3.EQ.-4) THEN
            IF(KK4.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (10).')
            A(3)=A(4)
         ENDIF
*        Z- SIDE:
         IF(KK5.GT.0) THEN
            A(5)=VHARM(DZ,ZZ(KK5),DIF(L,3),DIF(MAT(KK5),3),DDF(L,3),
     1      DDF(MAT(KK5),3))*VOL0/DZ
         ELSE IF(KK5.EQ.-1) THEN
            A(5)=VHARM(DZ,DZ,DIF(L,3),DZ*QFR(5)/DENOM,DDF(L,3),0.0)
     1      *VOL0/DZ
         ELSE IF(KK5.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(5)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*        Z+ SIDE:
         IF(KK6.GT.0) THEN
            A(6)=VHARM(DZ,ZZ(KK6),DIF(L,3),DIF(MAT(KK6),3),DDF(L,3),
     1      DDF(MAT(KK6),3))*VOL0/DZ
         ELSE IF(KK6.EQ.-1) THEN
            A(6)=VHARM(DZ,DZ,DIF(L,3),DZ*QFR(6)/DENOM,DDF(L,3),0.0)
     1      *VOL0/DZ
         ELSE IF(KK6.EQ.-2) THEN
            A(6)=0.0D0
         ELSE IF(KK6.EQ.-3) THEN
            A(6)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ELSE IF(KK6.EQ.-4) THEN
            IF(KK5.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (11).')
            A(6)=A(5)
         ENDIF
         IF(KK5.EQ.-4) THEN
            IF(KK6.EQ.-4) CALL XABORT('TRIDCO: INCONSISTENT SYME (12).')
            A(5)=A(6)
         ENDIF
      ENDIF
      RETURN
      END
