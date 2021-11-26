*DECK SYBRN2
      SUBROUTINE SYBRN2 (NREG,NSURF,A,B,Z,IZ,VOL,SIGT,TRONC,PVS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the DP-1 leakage and transmission probabilities for an
* heterogeneous non-sectorized square or rectangular cell. The tracks
* are computed by SYBRTK.
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
* NREG    number of regions in the cell.
* NSURF   number of surfaces.
* A       dimension of the external X sides.
* B       dimension of the external Y sides.
* Z       real tracking information.
* IZ      integer tracking information.
* VOL     volumes.
* SIGT    total macroscopic cross section.
* TRONC   voided block criterion.
*
*Parameters: output
* PVS     volume to surface probabilities:
*         XINF surfaces 1, 2 and 3;   XSUP surfaces 4, 5 and 6;
*         YINF surfaces 7, 8 and 9;   YSUP surfaces 10, 11 and 12.
* PSS     surface to surface probabilities in the following order:
*         PSS(i,j) is the probability from surface i to surface j.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,NSURF,IZ(*)
      REAL A,B,Z(*),VOL(NREG),SIGT(NREG),TRONC,PVS(3*NSURF,NREG),
     1 PSS(3*NSURF,3*NSURF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600,MKI5=600)
      PARAMETER (PI=3.141592654,ZI30=0.785398164,ZI40=0.666666667)
      LOGICAL ICARE
      REAL KI3,KI4,KI5
      INTEGER ISN(12,12)
      REAL PBB(28)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: COSINU
*----
*  BICKLEY TABLES
*----
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5),PAS5,XLIM5,L5
*
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
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(COSINU(2,NSURF))
*----
*  INTEGRATION USING THE TRACKING
*----
      ICARE=IZ(1).EQ.2
      ZERO=TRONC*(A+B)/(2.0*A*B)
      AOB=A/B
      CALL XDRSET(PBB,28,0.0)
      DO 15 I=1,12
      DO 10 J=1,NREG
      PVS(I,J)=0.0
   10 CONTINUE
   15 CONTINUE
      IZ0=2
      IZR=4
      DO 205 IA=1,IZ(2)
      DO 20 I=1,NSURF
      COSINU(1,I)=Z(IZR+1)
      COSINU(2,I)=Z(IZR+2)
      IZR=IZR+2
   20 CONTINUE
      MNT=IZ(IZ0+1)
      IZ0=IZ0+2
      IZR=IZR+1
      DO 200 IMNT=1,MNT
      NH=IZ(IZ0+1)
      NX=IZ(IZ0+2)
      ISURF=IZ(IZ0+3)+1
      JSURF=IZ(IZ0+NH+4)+1
      DO 190 INX=1,NX
      Z1=Z(IZR+1)
      IZR=IZR+1
      IF((ISURF.EQ.3).AND.(JSURF.EQ.2)) THEN
         Z1=0.5*Z1/A
         Z2=Z1*COSINU(1,ISURF)
         Z3=Z1*COSINU(2,ISURF)
         Z4=Z1*COSINU(1,ISURF)*COSINU(1,JSURF)
         Z5=Z1*COSINU(2,ISURF)*COSINU(1,JSURF)
         KI3=ZI30
         KI4=ZI40
         POP=0.0
         DO 40 I=1,NH
         III=IZ(IZ0+3+I)-NSURF+1
         SIGTI=SIGT(III)
         POP0=POP
         POP=POP+SIGTI*Z(IZR+I)
         IF(POP.LT.XLIM3) GO TO 30
         IF(SIGTI.LE.ZERO) GO TO 50
         PVS(7,III)=PVS(7,III)+2.0*KI3*Z1
         PVS(8,III)=PVS(8,III)+2.0*KI4*Z3
         GO TO 50
   30    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(7,III)=PVS(7,III)+2.0*TABKI(2,POP0)*Z(IZR+I)*Z1
            PVS(8,III)=PVS(8,III)+2.0*KI3*Z(IZR+I)*Z3
         ELSE
            PVS(7,III)=PVS(7,III)+2.0*(KI3-WI3)*Z1
            PVS(8,III)=PVS(8,III)+2.0*(KI4-WI4)*Z3
         ENDIF
         KI3=WI3
         KI4=WI4
   40    CONTINUE
   50    KI3=ZI30
         KI4=ZI40
         POP=0.0
         K=0
         DO 80 I=1,NH
         III=IZ(IZ0+3+I)-NSURF+1
         SIGTI=SIGT(III)
         POP0=POP
         POP=POP+SIGTI*Z(IZR+NH+1-I)
         IF(POP.LT.XLIM3) GO TO 70
         IF(SIGTI.LE.ZERO) GO TO 185
         PVS(1,III)=PVS(1,III)+2.0*KI3*Z1*AOB
         PVS(2,III)=PVS(2,III)+2.0*KI4*Z2*AOB
         GO TO 185
   70    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+2.0*TABKI(2,POP0)*Z(IZR+NH+1-I)*Z1*AOB
            PVS(2,III)=PVS(2,III)+2.0*KI3*Z(IZR+NH+1-I)*Z2*AOB
         ELSE
            PVS(1,III)=PVS(1,III)+2.0*(KI3-WI3)*Z1*AOB
            PVS(2,III)=PVS(2,III)+2.0*(KI4-WI4)*Z2*AOB
         ENDIF
         KI3=WI3
         KI4=WI4
   80    CONTINUE
         KI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
         PBB(1)=PBB(1)+KI3*Z1
         PBB(5)=PBB(5)+KI4*Z2
         PBB(9)=PBB(9)+KI4*Z3
         PBB(13)=PBB(13)+KI5*Z4
         PBB(21)=PBB(21)+KI5*Z5
         PBB(23)=PBB(23)+KI5*(Z1-Z5)
      ELSE IF((ISURF.EQ.3).AND.(JSURF.EQ.4)) THEN
         Z1=Z1/A
         Z2=Z1*COSINU(2,ISURF)
         Z3=Z1*COSINU(2,ISURF)*COSINU(2,JSURF)
         KI3=ZI30
         KI4=ZI40
         POP=0.0
         K=0
         DO 100 I=1,NH
         III=IZ(IZ0+3+I)-NSURF+1
         SIGTI=SIGT(III)
         POP0=POP
         POP=POP+SIGTI*Z(IZR+I)
         IF(POP.LT.XLIM3) GO TO 90
         IF(SIGTI.LE.ZERO) GO TO 110
         PVS(7,III)=PVS(7,III)+KI3*Z1
         PVS(8,III)=PVS(8,III)+KI4*Z2
         GO TO 110
   90    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(7,III)=PVS(7,III)+TABKI(2,POP0)*Z(IZR+I)*Z1
            PVS(8,III)=PVS(8,III)+KI3*Z(IZR+I)*Z2
         ELSE
            PVS(7,III)=PVS(7,III)+(KI3-WI3)*Z1
            PVS(8,III)=PVS(8,III)+(KI4-WI4)*Z2
         ENDIF
         KI3=WI3
         KI4=WI4
  100    CONTINUE
  110    KI3=ZI30
         KI4=ZI40
         POP=0.0
         DO 130 I=1,NH
         III=IZ(IZ0+3+I)-NSURF+1
         SIGTI=SIGT(III)
         POP0=POP
         POP=POP+SIGTI*Z(IZR+NH+1-I)
         IF(POP.LT.XLIM3) GO TO 120
         IF(SIGTI.LE.ZERO) GO TO 185
         PVS(7,III)=PVS(7,III)+KI3*Z1
         PVS(8,III)=PVS(8,III)+KI4*Z2
         GO TO 185
  120    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(7,III)=PVS(7,III)+TABKI(2,POP0)*Z(IZR+NH+1-I)*Z1
            PVS(8,III)=PVS(8,III)+KI3*Z(IZR+NH+1-I)*Z2
         ELSE
            PVS(7,III)=PVS(7,III)+(KI3-WI3)*Z1
            PVS(8,III)=PVS(8,III)+(KI4-WI4)*Z2
         ENDIF
         KI3=WI3
         KI4=WI4
  130    CONTINUE
         KI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
         PBB(2)=PBB(2)+2.0*KI3*Z1
         PBB(6)=PBB(6)+2.0*KI4*Z2
         PBB(14)=PBB(14)+2.0*KI5*Z3
         PBB(26)=PBB(26)+2.0*KI5*(Z1-Z3)
      ELSE IF((ISURF.EQ.1).AND.(JSURF.EQ.2)) THEN
         Z1=Z1/B
         Z2=Z1*COSINU(1,ISURF)
         Z3=Z1*COSINU(1,ISURF)*COSINU(1,JSURF)
         KI3=ZI30
         KI4=ZI40
         POP=0.0
         DO 150 I=1,NH
         III=IZ(IZ0+3+I)-NSURF+1
         SIGTI=SIGT(III)
         POP0=POP
         POP=POP+SIGTI*Z(IZR+I)
         IF(POP.LT.XLIM3) GO TO 140
         IF(SIGTI.LE.ZERO) GO TO 160
         PVS(1,III)=PVS(1,III)+KI3*Z1
         PVS(2,III)=PVS(2,III)+KI4*Z2
         GO TO 160
  140    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+TABKI(2,POP0)*Z(IZR+I)*Z1
            PVS(2,III)=PVS(2,III)+KI3*Z(IZR+I)*Z2
         ELSE
            PVS(1,III)=PVS(1,III)+(KI3-WI3)*Z1
            PVS(2,III)=PVS(2,III)+(KI4-WI4)*Z2
         ENDIF
         KI3=WI3
         KI4=WI4
  150    CONTINUE
  160    KI3=ZI30
         KI4=ZI40
         POP=0.0
         K=0
         DO 180 I=1,NH
         III=IZ(IZ0+3+I)-NSURF+1
         SIGTI=SIGT(III)
         POP0=POP
         POP=POP+SIGTI*Z(IZR+NH+1-I)
         IF(POP.LT.XLIM3) GO TO 170
         IF(SIGTI.LE.ZERO) GO TO 185
         PVS(1,III)=PVS(1,III)+KI3*Z1
         PVS(2,III)=PVS(2,III)+KI4*Z2
         GO TO 185
  170    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+TABKI(2,POP0)*Z(IZR+NH+1-I)*Z1
            PVS(2,III)=PVS(2,III)+KI3*Z(IZR+NH+1-I)*Z2
         ELSE
            PVS(1,III)=PVS(1,III)+(KI3-WI3)*Z1
            PVS(2,III)=PVS(2,III)+(KI4-WI4)*Z2
         ENDIF
         KI3=WI3
         KI4=WI4
  180    CONTINUE
         KI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
         PBB(4)=PBB(4)+2.0*KI3*Z1
         PBB(8)=PBB(8)+2.0*KI4*Z2
         PBB(16)=PBB(16)+2.0*KI5*Z3
         PBB(28)=PBB(28)+2.0*KI5*(Z1-Z3)
      ENDIF
  185 IZR=IZR+NH
  190 CONTINUE
      IZ0=IZ0+NH+4
  200 CONTINUE
  205 CONTINUE
*----
*  APPLY SYMMETRIES      
*----
      IF(ICARE) THEN
         PBB(1)=2.0*PBB(1)
         PBB(5)=PBB(5)+PBB(9)
         PBB(9)=PBB(5)
         PBB(13)=2.0*PBB(13)
         PBB(21)=PBB(21)+PBB(23)
         PBB(23)=PBB(21)
         PBB(4)=PBB(2)
         PBB(8)=PBB(6)
         PBB(16)=PBB(14)
         PBB(28)=PBB(26)
         DO 210 I=1,NREG
         PVS(7,I)=PVS(7,I)+PVS(1,I)
         PVS(8,I)=PVS(8,I)+PVS(2,I)
         PVS(1,I)=PVS(7,I)
         PVS(2,I)=PVS(8,I)
  210    CONTINUE
      ENDIF
      PBB(10)=PBB(6)
      PBB(12)=PBB(8)
      PBB(17)=PBB(9)
      PBB(19)=PBB(5)
      PBB(25)=PBB(13)
      PBB(3)=PBB(1)*AOB
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
      Z1=Z(1)
      Z2=Z(2)
      Z3=Z(3)
      Z4=Z(4)
      DO 280 I=1,4
      DEN0=PBB(I)
      DEN1=PBB(4+I)
      DEN2=PBB(8+I)
      PBB(I)=Z1*Z1*DEN0
      PBB(4+I)=Z1*(Z2*DEN1-Z3*DEN0)
      PBB(8+I)=Z1*(Z2*DEN2-Z3*DEN0)
      PBB(12+I)=Z2*Z2*PBB(12+I)-Z2*Z3*(DEN1+DEN2)+Z3*Z3*DEN0
      PBB(24+I)=Z4*Z4*PBB(24+I)
  280 CONTINUE
      DO 290 I=1,2
      DEN1=PBB(16+I)
      DEN2=PBB(18+I)
      PBB(16+I)=Z1*Z4*DEN1
      PBB(18+I)=Z1*Z4*DEN2
      PBB(20+I)=(Z2*PBB(20+I)-Z3*DEN1)*Z4
      PBB(22+I)=(Z2*PBB(22+I)-Z3*DEN2)*Z4
  290 CONTINUE
      DO 300 J=1,NREG
      DEN1=PVS(7,J)
      DEN2=PVS(1,J)
      PVS(7,J)=Z1*Z1*DEN1
      PVS(1,J)=Z1*Z1*DEN2
      PVS(8,J)=Z1*(Z2*PVS(8,J)-Z3*DEN1)
      PVS(2,J)=Z1*(Z2*PVS(2,J)-Z3*DEN2)
  300 CONTINUE
*
      Z1=0.0
      Z2=0.0
      Z4=0.0
      Z5=0.0
      DO 510 J=1,NREG
      X=4.0*VOL(J)
      SIGTI=SIGT(J)
      IF(SIGTI.GT.ZERO) THEN
         Z1=Z1+PVS(7,J)
         Z2=Z2+PVS(8,J)
         Z4=Z4+PVS(1,J)
         Z5=Z5+PVS(2,J)
         X=X*SIGTI
      ENDIF
      PVS(7,J)=PVS(7,J)*A/X
      PVS(8,J)=PVS(8,J)*A/X
      PVS(1,J)=PVS(1,J)*B/X
      PVS(2,J)=PVS(2,J)*B/X
  510 CONTINUE
      IF((Z1.GT.TRONC).AND.(Z2.GT.TRONC).AND.(Z4.GT.TRONC).AND.(Z5.GT.
     1 TRONC)) THEN
         Z1=(1.0-2.0*PBB(1)-PBB(2))/Z1
         Z2=(-2.0*PBB(9)-PBB(10))/Z2
         Z4=(1.0-2.0*PBB(3)-PBB(4))/Z4
         Z5=(-2.0*PBB(11)-PBB(12))/Z5
         DO 520 J=1,NREG
         PVS(7,J)=PVS(7,J)*Z1
         PVS(8,J)=PVS(8,J)*Z2
         PVS(1,J)=PVS(1,J)*Z4
         PVS(2,J)=PVS(2,J)*Z5
  520    CONTINUE
      ENDIF
*
      DO 540 I=1,NREG
      PVS(4,I)=PVS(1,I)
      PVS(5,I)=PVS(2,I)
      PVS(10,I)=PVS(7,I)
      PVS(11,I)=PVS(8,I)
  540 CONTINUE
      DO 560 JC=1,12
      DO 550 IC=1,12
      IB=ISN(IC,JC)
      IF(IB.LT.0) THEN
         PSS(IC,JC)=-PBB(-IB)
      ELSE IF(IB.GT.0) THEN
         PSS(IC,JC)=PBB(IB)
      ELSE
         PSS(IC,JC)=0.0
      ENDIF
      CONTINUE
  550 CONTINUE
  560 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(COSINU)
      RETURN
      END
