*DECK TRIHCO
      SUBROUTINE TRIHCO (IR,K,NEL,VOL0,MAT,DIF,DDF,SIDE,ZZ,KN,QFR,IWRK,
     1 IPR,A)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the value or the derivative or variation of mesh centered
* finite difference coefficients in element K for hexagonal geometry.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Benaboud
*
*Parameters: input
* IR      first dimension of matrix DIF.
* K       index of finite element under consideration.
* NEL     total number of finite elements.
* VOL0    volume of finite element under consideration.
* MAT     mixture index assigned to each element.
* DIF     directional diffusion coefficients.
* DDF     derivative or variation of directional diffusion coefficients.
* SIDE    side of the hexagons.
* ZZ      Z-directed mesh spacings.
* KN      element-ordered unknown list:
*         .GT.0 neighbour index;
*         =-1   void/albedo boundary condition;
*         =-2   reflection boundary condition;
*         =-3   ZERO flux boundary condition;
*         =-4   SYME boundary condition (axial symmetry).
* QFR     element-ordered boundary conditions.
* IWRK    non-void indices.
* IPR     type of MCFD coefficients:
*         .eq.0 direct MCFD coefficients calculation;
*         .eq.1 take derivative of MCFD coefficients;
*         .ge.2 take variation of MCFD coefficients.
*
*Parameters: output
* A       value or derivative or variation of mesh centered finite
*         difference coefficients.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IR,K,MAT(NEL),KN(8),IWRK(NEL),IPR
      REAL VOL0,DIF(IR,3),DDF(IR,3),SIDE,ZZ(NEL),QFR(8)
      DOUBLE PRECISION A(8)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION SHARM,VHARM,DHARM
*
*     FORMULE DIRECTE:
      SHARM(X1,X2,DIF1,DIF2)=2.0D0*DIF1*DIF2/(X1*DIF2+X2*DIF1)
*     FORMULE DE VARIATION:
      VHARM(X1,X2,DIF1,DIF2,DDF1,DDF2)=2.0D0*((DIF1+DDF1)*(DIF2+DDF2)
     1 /(X1*(DIF2+DDF2)+X2*(DIF1+DDF1))-DIF1*DIF2/(X1*DIF2+X2*DIF1))
*     FORMULE DE DERIVEE:
      DHARM(X1,X2,DIF1,DIF2,DDF1,DDF2)=2.0D0*(X1*DIF2*DIF2*DDF1+
     1 X2*DIF1*DIF1*DDF2)/(X1*DIF2+X2*DIF1)**2
*
      DENOM=2.0
      L=MAT(K)
      DZ=ZZ(K)
      DS=SQRT(3.0)*SIDE
      IF(IPR.EQ.0) THEN
*        COTE W NEGATIF:
         KK1=KN(6)
         IF(KK1.GT.0) THEN
            A(6)=SHARM(DS,DS,DIF(L,1),DIF(MAT(IWRK(KK1)),1))*DZ*SIDE
         ELSE IF(KK1.EQ.-1) THEN
            A(6)=SHARM(DS,DS,DIF(L,1),DS*QFR(6)/DENOM)*DZ*SIDE
         ELSE IF(KK1.EQ.-2) THEN
            A(6)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(6)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DZ*SIDE
         ENDIF
*        COTE W POSITIF:
         KK2=KN(3)
         IF(KK2.GT.0) THEN
            A(3)=SHARM(DS,DS,DIF(L,1),DIF(MAT(IWRK(KK2)),1))*DZ*SIDE
         ELSE IF(KK2.EQ.-1) THEN
            A(3)=SHARM(DS,DS,DIF(L,1),DS*QFR(3)/DENOM)*DZ*SIDE
         ELSE IF(KK2.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(3)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DZ*SIDE
         ENDIF
*        COTE X NEGATIF:
         KK3=KN(1)
         IF(KK3.GT.0) THEN
            A(1)=SHARM(DS,DS,DIF(L,1),DIF(MAT(IWRK(KK3)),1))*DZ*SIDE
         ELSE IF(KK3.EQ.-1) THEN
            A(1)=SHARM(DS,DS,DIF(L,1),DS*QFR(1)/DENOM)*DZ*SIDE
         ELSE IF(KK3.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(1)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DZ*SIDE
         ENDIF
*        COTE X POSITIF:
         KK4=KN(4)
         IF(KK4.GT.0) THEN
            A(4)=SHARM(DS,DS,DIF(L,1),DIF(MAT(IWRK(KK4)),1))*DZ*SIDE
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=SHARM(DS,DS,DIF(L,1),DS*QFR(4)/DENOM)*DZ*SIDE
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DZ*SIDE
         ENDIF
*        COTE Y NEGATIF:
         KK5=KN(2)
         IF(KK5.GT.0) THEN
            A(2)=SHARM(DS,DS,DIF(L,1),DIF(MAT(IWRK(KK5)),1))*DZ*SIDE
         ELSE IF(KK5.EQ.-1) THEN
            A(2)=SHARM(DS,DS,DIF(L,1),DS*QFR(2)/DENOM)*DZ*SIDE
         ELSE IF(KK5.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(2)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DZ*SIDE
         ENDIF
*        COTE Y POSITIF:
         KK6=KN(5)
         IF(KK6.GT.0) THEN
            A(5)=SHARM(DS,DS,DIF(L,1),DIF(MAT(IWRK(KK6)),1))*DZ*SIDE
         ELSE IF(KK6.EQ.-1) THEN
            A(5)=SHARM(DS,DS,DIF(L,1),DS*QFR(5)/DENOM)*DZ*SIDE
         ELSE IF(KK6.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK6.EQ.-3) THEN
            A(5)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DZ*SIDE
         ENDIF
*        COTE Z NEGATIF:
         KK7=KN(7)
         IF(KK7.GT.0) THEN
            A(7)=SHARM(DZ,ZZ(IWRK(KK7)),DIF(L,3),DIF(MAT(IWRK(KK7)),3))
     *      *VOL0/DZ
         ELSE IF(KK7.EQ.-1) THEN
            A(7)=SHARM(DZ,DZ,DIF(L,3),DZ*QFR(7)/DENOM)*VOL0/DZ
         ELSE IF(KK7.EQ.-2) THEN
            A(7)=0.0D0
         ELSE IF(KK7.EQ.-3) THEN
            A(7)=2.0D0*SHARM(DZ,DZ,DIF(L,3),DIF(L,3))*VOL0/DZ
         ENDIF
*        COTE Z POSITIF:
         KK8=KN(8)
         IF(KK8.GT.0) THEN
            A(8)=SHARM(DZ,ZZ(IWRK(KK8)),DIF(L,3),DIF(MAT(IWRK(KK8)),3))
     *      *VOL0/DZ
         ELSE IF(KK8.EQ.-1) THEN
            A(8)=SHARM(DZ,DZ,DIF(L,3),DZ*QFR(8)/DENOM)*VOL0/DZ
         ELSE IF(KK8.EQ.-2) THEN
            A(8)=0.0D0
         ELSE IF(KK8.EQ.-3) THEN
            A(8)=2.0D0*SHARM(DZ,DZ,DIF(L,3),DIF(L,3))*VOL0/DZ
         ENDIF
      ELSE IF(IPR.EQ.1) THEN
*        FORMULE DE DERIVEE.
*        COTE W NEGATIF:
         KK1=KN(6)
         IF(KK1.GT.0) THEN
            KK1=IWRK(KK1)
            A(6)=DHARM(DS,DS,DIF(L,1),DIF(MAT(KK1),1),DDF(L,1),
     1           DDF(MAT(KK1),1))*DZ*SIDE
         ELSE IF(KK1.EQ.-1) THEN
            A(6)=DHARM(DS,DS,DIF(L,1),DS*QFR(6)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK1.EQ.-2) THEN
            A(6)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(6)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE W POSITIF:
         KK2=KN(3)
         IF(KK2.GT.0) THEN
            KK2=IWRK(KK2)
            A(3)=DHARM(DS,DS,DIF(L,1),DIF(MAT(KK2),1),DDF(L,1),
     1           DDF(MAT(KK2),1))*DZ*SIDE
         ELSE IF(KK2.EQ.-1) THEN
            A(3)=DHARM(DS,DS,DIF(L,1),DS*QFR(3)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK2.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE X NEGATIF:
         KK3=KN(1)
         IF(KK3.GT.0) THEN
            KK3=IWRK(KK3)
            A(1)=DHARM(DS,DS,DIF(L,1),DIF(MAT(KK3),1),DDF(L,1),
     1           DDF(MAT(KK3),1))*DZ*SIDE
         ELSE IF(KK3.EQ.-1) THEN
            A(1)=DHARM(DS,DS,DIF(L,1),DS*QFR(1)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK3.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE X POSITIF:
         KK4=KN(4)
         IF(KK4.GT.0) THEN
            KK4=IWRK(KK4)
            A(4)=DHARM(DS,DS,DIF(L,1),DIF(MAT(KK4),1),DDF(L,1),
     1           DDF(MAT(KK4),1))*DZ*SIDE
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=DHARM(DS,DS,DIF(L,1),DS*QFR(4)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE Y NEGATIF:
         KK5=KN(2)
         IF(KK5.GT.0) THEN
            KK5=IWRK(KK5)
            A(2)=DHARM(DS,DS,DIF(L,1),DIF(MAT(KK5),1),DDF(L,1),
     1           DDF(MAT(KK5),1))*DZ*SIDE
         ELSE IF(KK5.EQ.-1) THEN
            A(2)=DHARM(DS,DS,DIF(L,1),DZ*QFR(2)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK5.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(2)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE Y POSITIF:
         KK6=KN(5)
         IF(KK6.GT.0) THEN
            KK6=IWRK(KK6)
            A(5)=DHARM(DS,DS,DIF(L,1),DIF(MAT(KK6),1),DDF(L,1),
     1           DDF(MAT(KK6),1))*DZ*SIDE
         ELSE IF(KK6.EQ.-1) THEN
            A(5)=DHARM(DS,DS,DIF(L,1),DS*QFR(5)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK6.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK6.EQ.-3) THEN
            A(5)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE Z NEGATIF:
         KK7=KN(7)
         IF(KK7.GT.0) THEN
            KK7=IWRK(KK7)
            A(7)=DHARM(DZ,ZZ(KK7),DIF(L,3),DIF(MAT(KK7),3),DDF(L,3),
     1           DDF(MAT(KK7),3))*VOL0/DZ
         ELSE IF(KK7.EQ.-1) THEN
            A(7)=DHARM(DZ,DZ,DIF(L,3),DZ*QFR(7)/DENOM,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK7.EQ.-2) THEN
            A(7)=0.0D0
         ELSE IF(KK7.EQ.-3) THEN
            A(7)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*        COTE Z POSITIF:
         KK8=KN(8)
         IF(KK8.GT.0) THEN
            KK8=IWRK(KK8)
            A(8)=DHARM(DZ,ZZ(KK8),DIF(L,3),DIF(MAT(KK8),3),DDF(L,3),
     1           DDF(MAT(KK8),3))*VOL0/DZ
         ELSE IF(KK8.EQ.-1) THEN
            A(8)=DHARM(DZ,DZ,DIF(L,3),DZ*QFR(8)/DENOM,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK8.EQ.-2) THEN
            A(8)=0.0D0
         ELSE IF(KK8.EQ.-3) THEN
            A(8)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
      ELSE
*        FORMULE DE VARIATION.
*        COTE W NEGATIF:
         KK1=KN(6)
         IF(KK1.GT.0) THEN
            KK1=IWRK(KK1)
            A(6)=VHARM(DS,DS,DIF(L,1),DIF(MAT(KK1),1),DDF(L,1),
     1           DDF(MAT(KK1),1))*DZ*SIDE
         ELSE IF(KK1.EQ.-1) THEN
            A(6)=VHARM(DS,DS,DIF(L,1),DS*QFR(6)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK1.EQ.-2) THEN
            A(6)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(6)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE W POSITIF:
         KK2=KN(3)
         IF(KK2.GT.0) THEN
            KK2=IWRK(KK2)
            A(3)=VHARM(DS,DS,DIF(L,1),DIF(MAT(KK2),1),DDF(L,1),
     1           DDF(MAT(KK2),1))*DZ*SIDE
         ELSE IF(KK2.EQ.-1) THEN
            A(3)=VHARM(DS,DS,DIF(L,1),DS*QFR(3)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK2.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE X NEGATIF:
         KK3=KN(1)
         IF(KK3.GT.0) THEN
            KK3=IWRK(KK3)
            A(1)=VHARM(DS,DS,DIF(L,1),DIF(MAT(KK3),1),DDF(L,1),
     1           DDF(MAT(KK3),1))*DZ*SIDE
         ELSE IF(KK3.EQ.-1) THEN
            A(1)=VHARM(DS,DS,DIF(L,1),DS*QFR(1)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK3.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(1)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE X POSITIF:
         KK4=KN(4)
         IF(KK4.GT.0) THEN
            KK4=IWRK(KK4)
            A(4)=VHARM(DS,DS,DIF(L,1),DIF(MAT(KK4),1),DDF(L,1),
     1           DDF(MAT(KK4),1))*DZ*SIDE
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=VHARM(DS,DS,DIF(L,1),DS*QFR(4)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE Y NEGATIF:
         KK5=KN(2)
         IF(KK5.GT.0) THEN
            KK5=IWRK(KK5)
            A(2)=VHARM(DS,DS,DIF(L,1),DIF(MAT(KK5),1),DDF(L,1),
     1           DDF(MAT(KK5),1))*DZ*SIDE
         ELSE IF(KK5.EQ.-1) THEN
            A(2)=VHARM(DS,DS,DIF(L,1),DS*QFR(2)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK5.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(2)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE Y POSITIF:
         KK6=KN(5)
         IF(KK6.GT.0) THEN
            KK6=IWRK(KK6)
            A(5)=VHARM(DS,DS,DIF(L,1),DIF(MAT(KK6),1),DDF(L,1),
     1           DDF(MAT(KK6),1))*DZ*SIDE
         ELSE IF(KK6.EQ.-1) THEN
            A(5)=VHARM(DS,DS,DIF(L,1),DS*QFR(5)/DENOM,DDF(L,1),0.0)
     1           *DZ*SIDE
         ELSE IF(KK6.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK6.EQ.-3) THEN
            A(5)=2.0D0*DDF(L,1)*DZ*SIDE/DS
         ENDIF
*        COTE Z NEGATIF:
         KK7=KN(7)
         IF(KK7.GT.0) THEN
            KK7=IWRK(KK7)
            A(7)=VHARM(DZ,ZZ(KK7),DIF(L,3),DIF(MAT(KK7),3),DDF(L,3),
     1           DDF(MAT(KK7),3))*VOL0/DZ
         ELSE IF(KK7.EQ.-1) THEN
            A(7)=VHARM(DZ,DZ,DIF(L,3),DZ*QFR(7)/DENOM,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK7.EQ.-2) THEN
            A(7)=0.0D0
         ELSE IF(KK7.EQ.-3) THEN
            A(7)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*        COTE Z POSITIF:
         KK8=KN(8)
         IF(KK8.GT.0) THEN
            KK8=IWRK(KK8)
            A(8)=VHARM(DZ,ZZ(KK8),DIF(L,3),DIF(MAT(KK8),3),DDF(L,3),
     1           DDF(MAT(KK8),3))*VOL0/DZ
         ELSE IF(KK8.EQ.-1) THEN
            A(8)=VHARM(DZ,DZ,DIF(L,3),DZ*QFR(8)/DENOM,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK8.EQ.-2) THEN
            A(8)=0.0D0
         ELSE IF(KK8.EQ.-3) THEN
            A(8)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
      ENDIF
      RETURN
      END
