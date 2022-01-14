*DECK RECT1
      SUBROUTINE RECT1 (NA,A,B,SIGT,TRONC,PII,PVS,PSS,ALPA,PWA)
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
* PII     volume to volume reduced probability.
* PVS     volume to surface probabilities:
*         XINF surface 1;   XSUP surface 2;
*         YINF surface 3;   YSUP surface 4.
* PSS     surface to surface probabilities in the following order:
*         PSS(i,j) is the probability from surface i to surface j.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NA
      REAL A,B,SIGT,TRONC,PII,PVS(4),PSS(4,4),ALPA(NA),PWA(NA)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600)
      PARAMETER (PI=3.141592654,ZI40=0.6666666667)
      DOUBLE PRECISION Y1,Y2,AOB,BOA,R,PHI,PHJ,RM,CO,SI,Z1,Z2,POP,ZI3,
     1 ZI4,GAR(4),QUOT(2),PHIB(2),ABSC(2),PBB(5)
      INTEGER ISL(4,4)
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      SAVE ISL
      DATA ISL/5,4,1,1,
     1         4,5,1,1,
     2         3,3,5,2,
     3         3,3,2,5/
*
      AOB=A/B
      BOA=B/A
      X=2.0*A*B/(A+B)
      IF(X*SIGT.LE.TRONC) GO TO 120
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
      DO 10 I=1,4
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
      Z1=PHIB(IZ)*WA*CO*SI/(PI*A)
      Z2=PHIB(IZ)*WA*(CO-QUOT(IZ)*SI)*2.0/PI
      POP=SIGT*ABSC(IZ)/CO
      ZI3=0.0
      ZI4=0.0
      IF(POP.GE.XLIM3) GO TO 40
      K=NINT(POP*PAS3)
      ZI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
      ZI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
40    GAR(L+1)=GAR(L+1)+(ZI40-ZI4)*Z1
      GAR(L+2)=GAR(L+2)+ZI3*Z2
      L=L+2
50    CONTINUE
55    CONTINUE
      PBB(1)=(GAR(1)+GAR(3))*SIG
      PBB(2)=GAR(2)
      PBB(4)=GAR(4)
      GO TO 110
*----
*  SQUARE CELL
*----
60    PHI1=0.25*PI
      GAR(1)=0.0
      GAR(2)=0.0
      DO 100 IL=1,NA
      X=0.5*(ALPA(IL)+1.0)
      WA=PWA(IL)
      PHI=PHI1*X
      CO=COS(PHI)
      SI=SIN(PHI)
      Z1=0.25*WA*CO*SI/A
      Z2=0.5*WA*(CO-SI)
      POP=SIGT*A/CO
      ZI3=0.0
      ZI4=0.0
      IF(POP.GE.XLIM3) GO TO 90
      K=NINT(POP*PAS3)
      ZI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
      ZI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
90    GAR(1)=GAR(1)+(ZI40-ZI4)*Z1
      GAR(2)=GAR(2)+ZI3*Z2
100   CONTINUE
      PBB(1)=2.0*GAR(1)*SIG
      PBB(2)=GAR(2)
      PBB(4)=PBB(2)
*
110   PBB(3)=PBB(1)*AOB
      PVS(3)=REAL(0.25*SIG*(1.0-2.0*PBB(1)-PBB(2))/B)
      PVS(1)=REAL(0.25*SIG*(1.0-2.0*PBB(3)-PBB(4))/A)
      PII=REAL(SIG*(1.0-2.0*(PVS(3)+PVS(1))))
      GO TO 130
*----
*  VOIDED CELL
*----
120   IF(A.GT.1.0E2*B) THEN
         PHI=0.5D0*PI
         R=A*(1.0D0+BOA*BOA*(0.5D0-0.125D0*BOA*BOA))
         RM=BOA*(0.5D0-0.125D0*BOA*BOA)
         Y1=LOG(DBLE((R+A)/B))
         PBB(3)=0.5D0*(1.0D0-RM)+SIGT*B*LOG(B/R)/PI
         PBB(1)=PBB(3)*BOA
         PBB(2)=(R-B)/A-4.0D0*SIGT*(B*PHI+B*BOA*LOG(B/R))/PI
     1   +SIGT*SIGT*B*B*(Y1-(R-B)/A)
         PBB(4)=RM
         PVS(3)=REAL((PHI+0.5D0*BOA*LOG(B/R))/PI-0.25D0*SIGT*(B*Y1-
     1   BOA*(R-B)))
         PVS(1)=REAL(-0.5D0*BOA*LOG(B/R)/PI)
         PII=REAL(0.5D0*B*Y1-0.5D0*BOA*(R-B))
      ELSE IF(B.GT.1.0E2*A) THEN
         PHJ=0.5D0*PI
         R=B*(1.0D0+AOB*AOB*(0.5D0-0.125D0*AOB*AOB))
         RM=AOB*(0.5D0-0.125D0*AOB*AOB)
         Y2=LOG(DBLE((R+B)/A))
         PBB(1)=0.5D0*(1.0D0-RM)+SIGT*A*LOG(A/R)/PI
         PBB(3)=PBB(1)*AOB
         PBB(2)=RM
         PBB(4)=(R-A)/B-4.0D0*SIGT*(A*PHJ+A*AOB*LOG(A/R))/PI
     1   +SIGT*SIGT*A*A*(Y2-(R-A)/B)
         PVS(3)=REAL(-0.5D0*AOB*LOG(A/R)/PI)
         PVS(1)=REAL((PHJ+0.5D0*AOB*LOG(A/R))/PI-0.25D0*SIGT*(A*Y2-
     1   AOB*(R-A)))
         PII=REAL(0.5D0*A*Y2-0.5D0*AOB*(R-A))
      ELSE
         PHI=ATAN(AOB)
         PHJ=ATAN(BOA)
         R=SQRT(A*A+B*B)
         Y1=LOG(DBLE((R+A)/B))
         Y2=LOG(DBLE((R+B)/A))
         PBB(1)=(A+B-R)/(A+A)+SIGT*(A*LOG(A/R)+B*BOA*LOG(B/R))/PI
     1   +SIGT*SIGT*(R*R*R-A*A*A-B*B*B)/(6.0D0*A)
         PBB(3)=PBB(1)*AOB
         PBB(2)=(R-B)/A-4.0D0*SIGT*(B*PHI+B*BOA*LOG(B/R))/PI
     1   +SIGT*SIGT*B*B*(Y1-(R-B)/A)
         PBB(4)=(R-A)/B-4.0D0*SIGT*(A*PHJ+A*AOB*LOG(A/R))/PI
     1   +SIGT*SIGT*A*A*(Y2-(R-A)/B)
         PVS(3)=REAL((0.5D0*BOA*LOG(B/R)-0.5D0*AOB*LOG(A/R)+PHI)/PI
     1   -0.25D0*SIGT*((R*R*R-A*A*A-B*B*B)/(3.0D0*A*B)-(R-B)*BOA+B*Y1))
         PVS(1)=REAL((0.5D0*AOB*LOG(A/R)-0.5D0*BOA*LOG(B/R)+PHJ)/PI
     1   -0.25D0*SIGT*((R*R*R-A*A*A-B*B*B)/(3.0D0*A*B)-(R-A)*AOB+A*Y2))
         PII=REAL((A*A*A+B*B*B-R*R*R)/6.0D0*A*B+0.5D0*B*Y1+0.5D0*A*Y2)
      ENDIF
  130 PVS(2)=PVS(1)
      PVS(4)=PVS(3)
      PBB(5)=0.0
      DO 150 JC=1,4
      DO 140 IC=1,4
      PSS(IC,JC)=REAL(PBB(ISL(IC,JC)))
  140 CONTINUE
  150 CONTINUE
      RETURN
      END
