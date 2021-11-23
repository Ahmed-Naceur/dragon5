*DECK TRITCO
      SUBROUTINE TRITCO (NEL,LL4,ISPLH,IR,IQF,K,KK1,KK2,KK3,KK4,KK5,
     1 VOL0,MAT,MATN,DIF,DDF,SIDE,ZZ,QFR,IPR,A)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Mesh centered finite difference coefficients in hexagonal geometry
* with triangular sub meshing.
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
* NEL     total number of finite elements.
* LL4     order of the system matrices.
* ISPLH   number of triangles (equal to 6*(ISPLH-1)**2).
* IR      first dimension for matrices DIF and DDF.
* IQF     index in array QFR.
* K       index of finite element.
* KK1     first neighbour of the triangular finite element.
* KK2     second neighbour of the triangular finite element.
* KK3     third neighbour of the triangular finite element.
* KK4     fourth neighbour of the triangular finite element.
* KK5     fifth neighbour of the triangular finite element.
* VOL0    volume of the finite element.
* MAT     mixture index assigned to each hexagon.
* MATN    mixture index assigned to each triangle.
* DIF     directional diffusion coefficients.
* DDF     variation of directional diffusion coefficients.
* SIDE    side of an hexagon.
* ZZ      Z-directed mesh spacings.
* QFR     element-ordered boundary conditions.
* IPR     type of matrix assembly:
*         =0: compute the system matrices;
*         =1: compute the derivative of system matrices;
*         =2 or =3: compute the variation of system matrices.
*
*Parameters: output
* A       mesh centered finite difference coefficients.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NEL,LL4,ISPLH,IR,IQF,K,KK1,KK2,KK3,KK4,KK5,MAT(NEL),
     1 MATN(LL4),IPR
      REAL VOL0,DIF(IR,3),DDF(IR,3),SIDE,ZZ(NEL),QFR(8)
      DOUBLE PRECISION A(5)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION SHARM,DHARM,VHARM
*     FORMULA WIHOUT VARIATION OR DERIVATIVE.
      SHARM(X1,X2,DIF1,DIF2)=2.0D0*DIF1*DIF2/(X1*DIF2+X2*DIF1)
*     FORMULA WITH DERIVATIVE.
      DHARM(X1,X2,DIF1,DIF2,DDF1,DDF2)=2.0D0*(X1*DIF2*DIF2*DDF1+
     1 X2*DIF1*DIF1*DDF2)/(X1*DIF2+X2*DIF1)**2
*     FORMULA WITH VARIATION.
      VHARM(X1,X2,DIF1,DIF2,DDF1,DDF2)=2.0D0*((DIF1+DDF1)*(DIF2+DDF2)
     1 /(X1*(DIF2+DDF2)+X2*(DIF1+DDF1))-DIF1*DIF2/(X1*DIF2+X2*DIF1))
*
      L=MAT(K)
      DZ=ZZ(K)
      DS=SIDE/(SQRT(3.0)*(ISPLH-1))
      DT=SIDE/(ISPLH-1)
      IF(IPR.EQ.0) THEN
*        FORMULE DIRECTE.
         IF(KK1.GT.0) THEN
            A(1)=SHARM(DS,DS,DIF(L,1),DIF(MATN(KK1),1))*DT*DZ
         ELSE IF(KK1.EQ.-1) THEN
            A(1)=SHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0)*DT*DZ
         ELSE IF(KK1.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(1)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DT*DZ
         ENDIF
*
         IF(KK2.GT.0) THEN
            A(2)=SHARM(DS,DS,DIF(L,1),DIF(MATN(KK2),1))*DT*DZ
         ELSE IF(KK2.EQ.-1) THEN
            A(2)=SHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0)*DT*DZ
         ELSE IF(KK2.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(2)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DT*DZ
         ENDIF
*
         IF(KK3.GT.0) THEN
            A(3)=SHARM(DS,DS,DIF(L,1),DIF(MATN(KK3),1))*DT*DZ
         ELSE IF(KK3.EQ.-1) THEN
            A(3)=SHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0)*DT*DZ
         ELSE IF(KK3.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(3)=2.0D0*SHARM(DS,DS,DIF(L,1),DIF(L,1))*DT*DZ
         ENDIF
*
         IF(KK4.GT.0) THEN
            A(4)=SHARM(DZ,ZZ(KK4),DIF(L,1),DIF(MAT(KK4),1))*VOL0/DZ
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=SHARM(DZ,DZ,DIF(L,1),DZ*QFR(7)/2.0)*VOL0/DZ
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*SHARM(DZ,DZ,DIF(L,1),DIF(L,1))*VOL0/DZ
         ENDIF
*
         IF(KK5.GT.0) THEN
            A(5)=SHARM(DZ,ZZ(KK5),DIF(L,1),DIF(MAT(KK5),1))*VOL0/DZ
         ELSE IF(KK5.EQ.-1) THEN
            A(5)=SHARM(DZ,DZ,DIF(L,1),DZ*QFR(8)/2.0)*VOL0/DZ
         ELSE IF(KK5.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(5)=2.0D0*SHARM(DZ,DZ,DIF(L,1),DIF(L,1))*VOL0/DZ
         ENDIF
*
      ELSE IF(IPR.EQ.1) THEN
*        FORMULE DE DERIVEE.
         IF(KK1.GT.0) THEN
            A(1)=DHARM(DS,DS,DIF(L,1),DIF(MATN(KK1),1),DDF(L,1),
     1           DDF(MATN(KK1),1))*DZ*DT
         ELSE IF(KK1.EQ.-1) THEN
            A(1)=DHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0,DDF(L,1),0.0)
     1           *DZ*DT
         ELSE IF(KK1.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(1)=2.0D0*DDF(L,1)*DZ*DT/DS
         ENDIF
*
         IF(KK2.GT.0) THEN
            A(2)=DHARM(DS,DS,DIF(L,1),DIF(MATN(KK2),1),DDF(L,1),
     1           DDF(MATN(KK2),1))*DZ*DT
         ELSE IF(KK2.EQ.-1) THEN
            A(2)=DHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0,DDF(L,1),0.0)
     1           *DZ*DT
         ELSE IF(KK2.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(2)=2.0D0*DDF(L,1)*DZ*DT/DS
         ENDIF
*
         IF(KK3.GT.0) THEN
            A(3)=DHARM(DS,DS,DIF(L,1),DIF(MATN(KK3),1),DDF(L,1),
     1           DDF(MATN(KK3),1))*DZ*DT
         ELSE IF(KK3.EQ.-1) THEN
            A(3)=DHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0,DDF(L,1),0.0)
     1           *DZ*DT
         ELSE IF(KK3.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,1)*DZ*DT/DS
         ENDIF
*
         IF(KK4.GT.0) THEN
            A(4)=DHARM(DZ,ZZ(KK4),DIF(L,3),DIF(MAT(KK4),3),DDF(L,3),
     1           DDF(MAT(KK4),3))*VOL0/DZ
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=DHARM(DZ,DZ,DIF(L,3),DZ*QFR(7)/2.0,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*
         IF(KK5.GT.0) THEN
            A(5)=DHARM(DZ,ZZ(KK5),DIF(L,3),DIF(MAT(KK5),3),DDF(L,3),
     1           DDF(MAT(KK5),3))*VOL0/DZ
         ELSE IF(KK5.EQ.-1) THEN
            A(5)=DHARM(DZ,DZ,DIF(L,3),DZ*QFR(8)/2.0,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK5.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(5)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*
      ELSE IF(IPR.GE.2) THEN
*        FORMULE DE VARIATION.
         IF(KK1.GT.0) THEN
            A(1)=VHARM(DS,DS,DIF(L,1),DIF(MATN(KK1),1),DDF(L,1),
     1           DDF(MATN(KK1),1))*DZ*DT
         ELSE IF(KK1.EQ.-1) THEN
            A(1)=VHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0,DDF(L,1),0.0)
     1           *DZ*DT
         ELSE IF(KK1.EQ.-2) THEN
            A(1)=0.0D0
         ELSE IF(KK1.EQ.-3) THEN
            A(1)=2.0D0*DDF(L,1)*DZ*DT/DS
         ENDIF
*
         IF(KK2.GT.0) THEN
            A(2)=VHARM(DS,DS,DIF(L,1),DIF(MATN(KK2),1),DDF(L,1),
     1           DDF(MATN(KK2),1))*DZ*DT
         ELSE IF(KK2.EQ.-1) THEN
            A(2)=VHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0,DDF(L,1),0.0)
     1           *DZ*DT
         ELSE IF(KK2.EQ.-2) THEN
            A(2)=0.0D0
         ELSE IF(KK2.EQ.-3) THEN
            A(2)=2.0D0*DDF(L,1)*DZ*DT/DS
         ENDIF
*
         IF(KK3.GT.0) THEN
            A(3)=VHARM(DS,DS,DIF(L,1),DIF(MATN(KK3),1),DDF(L,1),
     1           DDF(MATN(KK3),1))*DZ*DT
         ELSE IF(KK3.EQ.-1) THEN
            A(3)=VHARM(DS,DS,DIF(L,1),DS*QFR(IQF)/2.0,DDF(L,1),0.0)
     1           *DZ*DT
         ELSE IF(KK3.EQ.-2) THEN
            A(3)=0.0D0
         ELSE IF(KK3.EQ.-3) THEN
            A(3)=2.0D0*DDF(L,1)*DZ*DT/DS
         ENDIF
*
         IF(KK4.GT.0) THEN
            A(4)=VHARM(DZ,ZZ(KK4),DIF(L,3),DIF(MAT(KK4),3),DDF(L,3),
     1           DDF(MAT(KK4),3))*VOL0/DZ
         ELSE IF(KK4.EQ.-1) THEN
            A(4)=VHARM(DZ,DZ,DIF(L,3),DZ*QFR(7)/2.0,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK4.EQ.-2) THEN
            A(4)=0.0D0
         ELSE IF(KK4.EQ.-3) THEN
            A(4)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*
         IF(KK5.GT.0) THEN
            A(5)=VHARM(DZ,ZZ(KK5),DIF(L,3),DIF(MAT(KK5),3),DDF(L,3),
     1           DDF(MAT(KK5),3))*VOL0/DZ
         ELSE IF(KK5.EQ.-1) THEN
            A(5)=VHARM(DZ,DZ,DIF(L,3),DZ*QFR(8)/2.0,DDF(L,3),0.0)
     1           *VOL0/DZ
         ELSE IF(KK5.EQ.-2) THEN
            A(5)=0.0D0
         ELSE IF(KK5.EQ.-3) THEN
            A(5)=2.0D0*DDF(L,3)*VOL0/(DZ*DZ)
         ENDIF
*
      ENDIF
      RETURN
      END
