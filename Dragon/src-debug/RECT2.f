*DECK RECT2
      SUBROUTINE RECT2 (NA,A,B,SIGT,TRONC,PVV,PVS,PSS,ALPA,PWA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Integration of the collision probabilities of an homogeneous 2-D
* rectangle (DP-0 surface angular flux approximation).
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
* NA      number of Gauss-Legendre base points.
* A       X-thickness of the rectangle.
* B       Y-thickness of the rectangle.
* SIGT    cross section.
* TRONC   voided block criterion.
* ALPA    Gauss-Legendre base points.
* PWA     Gauss-Legendre weights.
*
*Parameters: output
* PVV     volume to volume reduced probability.
* PVS     volume to surface probabilities:
*         XINF surfaces 1, 2 and 3; XSUP surfaces 4, 5 and 6;
*         YINF surfaces 7, 8 and 9; YSUP surfaces 10, 11 and 12.
* PSS     surface to surface probabilities in the following order:
*         PSS(i,j) is the probability from surface i to surface j.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NA
      REAL A,B,SIGT,TRONC,PVV,PVS(12),PSS(12,12),ALPA(NA),PWA(NA)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600,MKI5=600)
      PARAMETER (PI=3.141592654,ZI40=0.6666666667,ZI50=0.589048623,
     1 ZI60=0.533333333,COEF1=4.24264069,COEF2=2.82842712)
      LOGICAL VOID
      DOUBLE PRECISION Y1,Y2,AOB,BOA,R,DPREC2,DPREC5,PHI,PHJ,RM,CO,SI,
     1 Z1,Z2,POP,COSI,CO2,ZI3,ZI4,ZI5,ZI6,DEN1,DEN2,GAR(20),QUOT(2),
     2 PHIB(2),ABSC(2),PBB(28)
      INTEGER ISN(12,12)
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5),PAS5,XLIM5,L5
      SAVE ISN
      DATA ISN/  0,  0,  0,  4, 12,  0,  1,  9,-19,  1,  9, 19,
     1           0,  0,  0,  8, 16,  0,  5, 13,-23,  5, 13, 23,
     2           0,  0,  0,  0,  0, 28, 17, 21,-25,-17,-21,-25,
     3           4, 12,  0,  0,  0,  0,  1,  9, 19,  1,  9,-19,
     4           8, 16,  0,  0,  0,  0,  5, 13, 23,  5, 13,-23,
     5           0,  0, 28,  0,  0,  0,-17,-21,-25, 17, 21,-25,
     6           3, 11, 20,  3, 11,-20,  0,  0,  0,  2, 10,  0,
     7           7, 15, 24,  7, 15,-24,  0,  0,  0,  6, 14,  0,
     8         -18,-22,-27, 18, 22,-27,  0,  0,  0,  0,  0, 26,
     9           3, 11,-20,  3, 11, 20,  2, 10,  0,  0,  0,  0,
     1           7, 15,-24,  7, 15, 24,  6, 14,  0,  0,  0,  0,
     2          18, 22,-27,-18,-22,-27,  0,  0, 26,  0,  0,  0/
*
      AOB=A/B
      BOA=B/A
      X=2.0*A*B/(A+B)
      VOID=(X*SIGT.LE.TRONC)
      SIG=1.0E10
      IF(VOID) GO TO 105
      SIG=1.0/SIGT
      IF(A.EQ.B) GO TO 60
*----
*  RECTANGULAR CELL
*----
      QUOT(1)=BOA
      QUOT(2)=AOB
      PHIB(1)=ATAN(AOB)
      PHIB(2)=ATAN(BOA)
      ABSC(1)=B
      ABSC(2)=A
      DO 10 I=1,20
      GAR(I)=0.0
10    CONTINUE
      DO 55 IL=1,NA
      X=0.5*(ALPA(IL)+1.0)
      WA=PWA(IL)
      L=0
      DO 50 IZ=1,2
      PHI=PHIB(IZ)*X
      CO=COS(PHI)
      SI=SIN(PHI)
      COSI=CO*SI
      CO2=CO*CO
      Z1=PHIB(IZ)*WA*COSI/(PI*A)
      Z2=PHIB(IZ)*WA*(CO-QUOT(IZ)*SI)*2.0/PI
      POP=SIGT*ABSC(IZ)/CO
      ZI3=0.0
      ZI4=0.0
      ZI5=0.0
      IF(POP.GE.XLIM3) GO TO 40
      K=NINT(POP*PAS3)
      ZI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
      ZI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
      ZI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
40    ZI6=0.8*ZI4+0.2*POP*(ZI3-ZI5)
      GAR(L+1)=GAR(L+1)+(ZI40-ZI4)*Z1
      GAR(L+2)=GAR(L+2)+(ZI50-ZI5)*Z1*CO
      GAR(L+3)=GAR(L+3)+(ZI50-ZI5)*Z1*SI
      GAR(L+4)=GAR(L+4)+(ZI60-ZI6)*Z1*COSI
      GAR(L+5)=GAR(L+5)+(ZI60-ZI6)*Z1*CO2
      GAR(L+6)=GAR(L+6)+(ZI60-ZI6)*Z1
      GAR(L+7)=GAR(L+7)+ZI3*Z2
      GAR(L+8)=GAR(L+8)+ZI4*Z2*CO
      GAR(L+9)=GAR(L+9)+ZI5*Z2*CO2
      GAR(L+10)=GAR(L+10)+ZI5*Z2
      L=L+10
50    CONTINUE
55    CONTINUE
      PBB(1)=(GAR(1)+GAR(11))*SIG
      PBB(2)=GAR(7)
      PBB(4)=GAR(17)
      PBB(5)=(GAR(3)+GAR(12))*SIG
      PBB(6)=GAR(8)
      PBB(8)=GAR(18)
      PBB(9)=(GAR(2)+GAR(13))*SIG
      PBB(10)=GAR(8)
      PBB(12)=GAR(18)
      PBB(13)=(GAR(4)+GAR(14))*SIG
      PBB(14)=GAR(9)
      PBB(16)=GAR(19)
      PBB(17)=PBB(9)
      PBB(19)=PBB(5)
      PBB(21)=(GAR(5)+GAR(16)-GAR(15))*SIG
      PBB(23)=(GAR(6)-GAR(5)+GAR(15))*SIG
      PBB(25)=PBB(13)
      PBB(26)=GAR(10)-GAR(9)
      PBB(28)=GAR(20)-GAR(19)
      GO TO 110
*----
*  SQUARE CELL
*----
60    PHI1=0.25*PI
      DO 65 I=1,8
      GAR(I)=0.0
65    CONTINUE
      DO 100 IL=1,NA
      X=0.5*(ALPA(IL)+1.0)
      WA=PWA(IL)
      PHI=PHI1*X
      CO=COS(PHI)
      SI=SIN(PHI)
      COSI=CO*SI
      Z1=0.25*WA*COSI/A
      Z2=0.5*WA*(CO-SI)
      POP=SIGT*A/CO
      ZI3=0.0
      ZI4=0.0
      ZI5=0.0
      IF(POP.GE.XLIM3) GO TO 90
      K=NINT(POP*PAS3)
      ZI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
      ZI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
      ZI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
90    ZI6=0.8*ZI4+0.2*POP*(ZI3-ZI5)
      GAR(1)=GAR(1)+(ZI40-ZI4)*Z1
      GAR(2)=GAR(2)+(ZI50-ZI5)*Z1*(CO+SI)
      GAR(3)=GAR(3)+(ZI60-ZI6)*Z1*COSI
      GAR(4)=GAR(4)+(ZI60-ZI6)*Z1
      GAR(5)=GAR(5)+ZI3*Z2
      GAR(6)=GAR(6)+ZI4*Z2*CO
      GAR(7)=GAR(7)+ZI5*Z2*CO*CO
      GAR(8)=GAR(8)+ZI5*Z2
100   CONTINUE
      PBB(1)=2.0*GAR(1)*SIG
      PBB(2)=GAR(5)
      PBB(4)=PBB(2)
      PBB(5)=GAR(2)*SIG
      PBB(6)=GAR(6)
      PBB(8)=PBB(6)
      PBB(9)=PBB(5)
      PBB(10)=PBB(6)
      PBB(12)=PBB(6)
      PBB(13)=2.0*GAR(3)*SIG
      PBB(14)=GAR(7)
      PBB(16)=PBB(14)
      PBB(17)=GAR(2)*SIG
      PBB(19)=PBB(17)
      PBB(21)=GAR(4)*SIG
      PBB(23)=PBB(21)
      PBB(25)=2.0*GAR(3)*SIG
      PBB(26)=GAR(8)-GAR(7)
      PBB(28)=PBB(26)
      GO TO 110
*----
*  VOIDED CELL
*----
105   Z1=1.0/(8.0*A*B)
      IF(A.GT.1.0E2*B) THEN
         PHI=0.5D0*PI
         R=A*(1.0D0+BOA*BOA*(0.5D0-0.125D0*BOA*BOA))
         RM=BOA*(0.5D0-0.125D0*BOA*BOA)
         Y1=LOG(DBLE((R+A)/B))
         PBB(1)=(0.5D0*(1.0D0-RM)+SIGT*B*LOG(B/R)/PI)*BOA
         PBB(2)=(R-B)/A-4.0D0*SIGT*(B*PHI+B*BOA*LOG(B/R))/PI
     1   +SIGT*SIGT*B*B*(Y1-(R-B)/A)
         PBB(4)=RM
         PBB(5)=((1.0-RM)/3.0-0.25*SIGT*B*Y1)*BOA
         PBB(6)=4.0*PHI/(3.0*PI)-SIGT*BOA*(R-B)
     1   +2.0*SIGT*SIGT*(B*B*PHI+B*B*BOA*LOG(B/R))/PI
         PBB(8)=RM*2.0/3.0
         PBB(9)=(-0.25*SIGT*B*(1.0-BOA*2.0/3.0))*BOA
         PBB(13)=(0.125*A/R-SIGT*B/6.0)*BOA
         PBB(14)=(R*(R-B)+A*A)/(4.0*A*R)-4.0*SIGT*B*PHI/(3.0*PI)
     1   +0.5*SIGT*SIGT*B*BOA*(R-B)
         PBB(16)=RM*4.0/9.0
         PBB(21)=0.125*(1.0-RM*RM*B/R)*BOA
         PBB(23)=0.125*((B+B-B*B/R)/B-RM)*BOA
         PBB(26)=0.25*PBB(2)*(R-B)/R
         PBB(28)=0.25*RM*RM*B/R
         PVS(1)=REAL(-0.5D0*BOA*LOG(B/R)/PI)
         PVS(7)=REAL((PHI+0.5D0*BOA*LOG(B/R))/PI-0.25D0*SIGT*(B*Y1-
     1   BOA*(R-B)))
         DPREC2=Z1*B*(2.0*(R-B)+B*(1.0-BOA*2.0/3.0))
     1   -0.5*SIGT*B*(PHI+BOA*LOG(B/R))/PI
         DPREC5=Z1*B*B*Y1
         PVV=REAL(0.5D0*B*Y1-0.5D0*BOA*(R-B))
      ELSE IF(B.GT.1.0E2*A) THEN
         PHJ=0.5D0*PI
         R=B*(1.0D0+AOB*AOB*(0.5D0-0.125D0*AOB*AOB))
         RM=AOB*(0.5D0-0.125D0*AOB*AOB)
         Y2=LOG(DBLE((R+B)/A))
         PBB(1)=0.5D0*(1.0D0-RM)+SIGT*A*LOG(A/R)/PI
         PBB(2)=RM
         PBB(4)=(R-A)/B-4.0D0*SIGT*(A*PHJ+A*AOB*LOG(A/R))/PI
     1   +SIGT*SIGT*A*A*(Y2-(R-A)/B)
         PBB(5)=-0.25*SIGT*A*(1.0-AOB*2.0/3.0)
         PBB(6)=RM*2.0/3.0
         PBB(8)=4.0*PHJ/(3.0*PI)-SIGT*AOB*(R-A)
     1   +2.0*SIGT*SIGT*(A*A*PHJ+A*A*AOB*LOG(A/R))/PI
         PBB(9)=(1.0-RM)/3.0-0.25*SIGT*A*Y2
         PBB(13)=0.125*B/R-SIGT*A/6.0
         PBB(14)=RM*4.0/9.0
         PBB(16)=(R*(R-A)+B*B)/(4.0*B*R)-4.0*SIGT*A*PHJ/(3.0*PI)
     1   +0.5*SIGT*SIGT*A*AOB*(R-A)
         PBB(21)=0.125*((A+A-A*A/R)/A-RM)
         PBB(23)=0.125*(1.0-RM*RM*A/R)
         PBB(26)=0.25*RM*RM*A/R
         PBB(28)=0.25*PBB(4)*(R-A)/R
         PVS(1)=REAL((PHJ+0.5D0*AOB*LOG(A/R))/PI-0.25D0*SIGT*(A*Y2-
     1   AOB*(R-A)))
         PVS(7)=REAL(-0.5D0*AOB*LOG(A/R)/PI)
         DPREC2=Z1*A*A*Y2
         DPREC5=Z1*A*(2.0*(R-A)+A*(1.0-AOB*2.0/3.0))
     1   -0.5*SIGT*A*(PHJ+AOB*LOG(A/R))/PI
         PVV=REAL(0.5D0*A*Y2-0.5D0*AOB*(R-A))
      ELSE
         PHI=ATAN(AOB)
         PHJ=ATAN(BOA)
         R=SQRT(A*A+B*B)
         Y1=LOG(DBLE((R+A)/B))
         Y2=LOG(DBLE((R+B)/A))
         PBB(1)=(A+B-R)/(A+A)+SIGT*(A*LOG(A/R)+B*BOA*LOG(B/R))/PI
     1   +SIGT*SIGT*(R*R*R-A*A*A-B*B*B)/(6.0*A)
         PBB(2)=(R-B)/A-4.0*SIGT*(B*PHI+B*BOA*LOG(B/R))/PI
     1   +SIGT*SIGT*B*B*(Y1-(R-B)/A)
         PBB(4)=(R-A)/B-4.0*SIGT*A*(PHJ+AOB*LOG(A/R))/PI
     1   +SIGT*SIGT*A*A*(Y2-(R-A)/B)
         PBB(5)=(2.0*PHI*B)/(3.0*PI*A)+0.25*SIGT*(-B*BOA*Y1+R-A)
     1   +SIGT*SIGT*(B*B-B*B*BOA*PHI-A*A*LOG(A/R))/(3.0*PI)
         PBB(6)=4.0*PHI/(3.0*PI)-SIGT*BOA*(R-B)
     1   +2.0*SIGT*SIGT*(B*B*PHI+B*B*BOA*LOG(B/R))/PI
         PBB(8)=4.0*PHJ/(3.0*PI)-SIGT*AOB*(R-A)
     1   +2.0*SIGT*SIGT*A*A*(PHJ+AOB*LOG(A/R))/PI
         PBB(9)=2.0*PHJ/(3.0*PI)+0.25*SIGT*(-A*Y2+BOA*(R-B))
     1   +SIGT*SIGT*(-B*B*BOA*LOG(B/R)+A*B-A*A*PHJ)/(3.0*PI)
         PBB(13)=B/(8.0*R)-SIGT*(B*BOA*PHI+A*PHJ-B)/(3.0*PI)
     1   +SIGT*SIGT*(B*B*BOA*Y1+A*A*Y2-B*R)/12.0
         PBB(14)=(R*(R-B)+A*A)/(4.0*A*R)-4.0*SIGT*B*PHI/(3.0*PI)
     1   +0.5*SIGT*SIGT*B*BOA*(R-B)
         PBB(16)=(R*(R-A)+B*B)/(4.0*B*R)-4.0*SIGT*A*PHJ/(3.0*PI)
     1   +0.5*SIGT*SIGT*A*AOB*(R-A)
         PBB(21)=(A+A+B-A*A/R-R)/(8.0*A)+2.0*SIGT*A*LOG(A/R)/(3.0*PI)
     1   +SIGT*SIGT*(2.0*A*(R-A)+B*BOA*(R-B))/12.0
         PBB(23)=(A+B+B-B*B/R-R)/(8.0*A)+2.0*SIGT*B*BOA*LOG(B/R)/
     1   (3.0*PI)+SIGT*SIGT*(2.0*B*BOA*(R-B)-A*(R-A))/12.0
         PBB(26)=0.25*PBB(2)*(R-B)/R
         PBB(28)=0.25*PBB(4)*(R-A)/R
         PVS(1)=REAL((0.5*AOB*LOG(A/R)-0.5*BOA*LOG(B/R)+PHJ)/PI
     1   -0.25*SIGT*((R*R*R-A*A*A-B*B*B)/(3.0*A*B)-(R-A)*AOB+A*Y2))
         PVS(7)=REAL((0.5*BOA*LOG(B/R)-0.5*AOB*LOG(A/R)+PHI)/PI
     1   -0.25*SIGT*((R*R*R-A*A*A-B*B*B)/(3.0*A*B)-(R-B)*BOA+B*Y1))
         DPREC2=Z1*(B*(R-B)+A*A*Y2)-SIGT*(A+2.0*B*BOA*LOG(B/R)
     1   +3.0*B*PHI-A*AOB*PHJ)/(6.0*PI)
         DPREC5=Z1*(A*(R-A)+B*B*Y1)-SIGT*(B+2.0*A*AOB*LOG(A/R)
     1   +3.0*A*PHJ-B*BOA*PHI)/(6.0*PI)
         PVV=REAL((A*A*A+B*B*B-R*R*R)/(6.0*A*B)+0.5*B*Y1+0.5*A*Y2)
      ENDIF
      PBB(10)=PBB(6)
      PBB(12)=PBB(8)
      PBB(17)=PBB(9)
      PBB(19)=PBB(5)
      PBB(25)=PBB(13)
      PVS(2)=REAL(COEF1*DPREC5-COEF2*PVS(1))
      PVS(3)=0.0
      PVS(8)=REAL(COEF1*DPREC2-COEF2*PVS(7))
      PVS(9)=0.0
*
110   PBB(3)=PBB(1)*AOB
      PBB(7)=PBB(9)*AOB
      PBB(11)=PBB(5)*AOB
      PBB(15)=PBB(13)*AOB
      PBB(18)=PBB(19)*AOB
      PBB(20)=PBB(17)*AOB
      PBB(22)=PBB(23)*AOB
      PBB(24)=PBB(21)*AOB
      PBB(27)=PBB(25)*AOB
*----
*  ORTHONORMALIZATION
*----
      DO 120 I=1,4
      DEN1=PBB(4+I)
      DEN2=PBB(8+I)
      PBB(4+I)=COEF1*DEN1-COEF2*PBB(I)
      PBB(8+I)=COEF1*DEN2-COEF2*PBB(I)
      PBB(12+I)=18.0*PBB(12+I)-12.0*(DEN1+DEN2)+8.0*PBB(I)
      PBB(24+I)=4.0*PBB(24+I)
120   CONTINUE
      DO 130 I=1,2
      DEN1=PBB(16+I)
      DEN2=PBB(18+I)
      PBB(16+I)=2.0*DEN1
      PBB(18+I)=2.0*DEN2
      PBB(20+I)=2.0*(COEF1*PBB(20+I)-COEF2*DEN1)
      PBB(22+I)=2.0*(COEF1*PBB(22+I)-COEF2*DEN2)
130   CONTINUE
*
      IF(.NOT.VOID) THEN
         PVS(7)=REAL(0.25*SIG*(1.0-2.0*PBB(1)-PBB(2))/B)
         PVS(8)=REAL(-0.25*SIG*(2.0*PBB(9)+PBB(10))/B)
         PVS(9)=0.0
         PVS(1)=REAL(0.25*SIG*(1.0-2.0*PBB(3)-PBB(4))/A)
         PVS(2)=REAL(-0.25*SIG*(2.0*PBB(11)+PBB(12))/A)
         PVS(3)=0.0
         PVV=REAL(SIG*(1.0-2.0*(PVS(7)+PVS(1))))
      ENDIF
      PVS(4)=PVS(1)
      PVS(5)=PVS(2)
      PVS(6)=PVS(3)
      PVS(10)=PVS(7)
      PVS(11)=PVS(8)
      PVS(12)=PVS(9)
      DO 150 JC=1,12
      DO 140 IC=1,12
      PSS(IC,JC)=0.0
      IB=ISN(IC,JC)
      IF(IB.LT.0) THEN
         PSS(IC,JC)=REAL(-PBB(-IB))
      ELSE IF(IB.GT.0) THEN
         PSS(IC,JC)=REAL(PBB(IB))
      ENDIF
  140 CONTINUE
  150 CONTINUE
      RETURN
      END
