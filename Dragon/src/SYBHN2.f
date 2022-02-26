*DECK SYBHN2
      SUBROUTINE SYBHN2(NREG,NSURF,SIDE,Z,IZ,VOL,SIGT,TRONC,PVS,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the DP-1 leakage and transmission probabilities for an
* heterogeneous non-sectorized hexagonal cell. The tracks are computed
* by SYBHTK.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NREG    number regions in the cell.
* NSURF   number of surfaces.
* SIDE    length of one of sides of the hexagon.
* Z       real integration mesh.
* IZ      integer integration mesh.
* VOL     volumes.
* SIGT    total cross sections.
* TRONC   voided block cutoff criterion.
*
*Parameters: output
* PVS     leakage probability:
*         PVS(j,i) for volume i to side j with i=1,nr and j=1,18.
* PSS     transmission probability:
*         PSS(i,j) for side i to side j with i=1,18 and j=1,18.
*
*Reference:
* M. Ouisloumen, Resolution par la methode des probabilites de
* collision de l'equation integrale du transport a deux et trois
* dimensions en geometrie hexagonale, Ph. D. thesis, Ecole
* Polytechnique de Montreal, Montreal, October 1993.
*
*Comments:
* hexagone surface identification. 
*                                        side a,b,c
*                 side 4,5,6             dir a -> isotropic
*                  xxxxxxxx              dir c -> tangent to surface
*                 x        x             dir b -> normal  to surface
*   side 7,8,9   x          x side 1,2,3
*               x            x
*              x              x
*               x            x
* side 10,11,12  x          x side 16,17,18
*                 x        x
*                  xxxxxxxx
*                side 13,14,15
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,NSURF,IZ(*)
      REAL Z(*),SIDE,SIGT(NREG),TRONC,VOL(NREG),PVS(3*NSURF,NREG),
     1 PSS(3*NSURF,3*NSURF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600,MKI5=600)
      PARAMETER (SQ3=1.732050807568877,SQ2=1.414213562373095,
     1 PI=3.141592653589793,ZI30=0.785398164,ZI40=0.666666667)
      REAL KI3,KI4,KI5
      REAL PBB(16)
      INTEGER IROT(18,18)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: COSINU
*----
*  BICKLEY TABLES
*----
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5),PAS5,XLIM5,L5
*
      SAVE IROT
      DATA IROT /
     +          0, 0, 0, 1, 2,-3, 7, 8,-9,13,14, 0, 7, 8, 9, 1, 2, 3,
     +          0, 0, 0, 2, 4, 5, 8,10,11,14,15, 0, 8,10,-11, 2, 4,-5,
     +          0, 0, 0, 3,-5, 6, 9,-11,12,0, 0,16,-9,11,12,-3, 5,6,
     +          1, 2, 3, 0, 0, 0, 1, 2,-3, 7, 8,-9,13,14, 0, 7, 8, 9,
     +          2, 4,-5, 0, 0, 0, 2, 4, 5, 8,10,11,14,15, 0, 8,10,-11,
     +         -3, 5, 6, 0, 0, 0, 3,-5, 6, 9,-11,12,0, 0,16,-9,11,12,
     +          7, 8, 9, 1, 2, 3, 0, 0, 0, 1, 2,-3, 7, 8,-9,13,14, 0,
     +          8,10,-11, 2, 4,-5, 0, 0, 0, 2, 4, 5, 8,10,11,14,15, 0,
     +         -9,11,12,-3, 5, 6, 0, 0, 0, 3,-5, 6, 9,-11,12,0, 0,16,
     +          13,14, 0, 7, 8, 9, 1, 2, 3, 0, 0, 0, 1, 2,-3, 7, 8,-9,
     +          14,15, 0, 8,10,-11, 2, 4,-5, 0, 0, 0, 2, 4, 5, 8,10,11,
     +          0, 0,16,-9,11,12,-3, 5, 6, 0, 0, 0, 3,-5, 6, 9,-11,12,
     +          7, 8,-9, 13,14, 0, 7, 8, 9, 1, 2, 3, 0, 0, 0, 1, 2,-3,
     +          8,10,11,14,15, 0, 8,10,-11, 2, 4,-5, 0, 0, 0, 2, 4, 5,
     +          9,-11,12,  0, 0,16,-9,11,12,-3, 5, 6, 0, 0, 0, 3,-5, 6,
     +          1, 2,-3, 7, 8,-9, 13,14, 0, 7, 8, 9, 1, 2, 3, 0, 0, 0,
     +          2, 4, 5, 8,10,11,14,15, 0, 8,10,-11, 2, 4,-5, 0, 0, 0,
     +          3,-5, 6, 9,-11,12,  0, 0,16,-9,11,12,-3, 5, 6, 0, 0, 0/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(COSINU(2,NSURF))
*----
*  INTEGRATION USING THE TRACKING
*----
      ZERO=TRONC/(SQ3*SIDE)
      CALL XDRSET(PBB,16,0.0)
      DO 15 I=1,18
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
      Z1=Z(IZR+1)/SIDE
      IZR=IZR+1
      DW1=Z1*COSINU(1,ISURF)
      DW2=Z1*COSINU(1,JSURF)
      IF((ISURF.EQ.5).AND.(JSURF.EQ.6)) THEN
         W610=SQ3*COSINU(2,ISURF)+COSINU(1,ISURF)
         W611=2.0*COSINU(2,JSURF)*COSINU(2,ISURF)
         Z2=Z1*W611
         Z3=Z1*W610
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
         PVS(1,III)=PVS(1,III)+KI3*Z1
         PVS(2,III)=PVS(2,III)+KI4*DW1
         GO TO 50
   30    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+TABKI(2,POP0)*Z(IZR+I)*Z1
            PVS(2,III)=PVS(2,III)+KI3*Z(IZR+I)*DW1
         ELSE
            PVS(1,III)=PVS(1,III)+(KI3-WI3)*Z1
            PVS(2,III)=PVS(2,III)+(KI4-WI4)*DW1
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
         PVS(1,III)=PVS(1,III)+KI3*Z1
         PVS(2,III)=PVS(2,III)+KI4*DW2
         GO TO 185
   70    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+TABKI(2,POP0)*Z(IZR+NH+1-I)*Z1
            PVS(2,III)=PVS(2,III)+KI3*Z(IZR+NH+1-I)*DW2
         ELSE
            PVS(1,III)=PVS(1,III)+(KI3-WI3)*Z1
            PVS(2,III)=PVS(2,III)+(KI4-WI4)*DW2
         ENDIF
         KI3=WI3
         KI4=WI4
   80    CONTINUE
         KI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
         PBB(1)=PBB(1)+KI3*Z1
         PBB(3)=PBB(3)+KI4*Z3
         PBB(5)=PBB(5)+KI5*Z1
         PBB(6)=PBB(6)+KI5*Z2
      ELSE IF((ISURF.EQ.5).AND.(JSURF.EQ.1)) THEN
         W610=SQ3*COSINU(1,JSURF)+COSINU(2,JSURF)
         W511=2.0*COSINU(2,ISURF)*COSINU(2,JSURF)
         Z2=Z1*W511
         Z3=Z1*W610
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
         PVS(1,III)=PVS(1,III)+KI3*Z1
         PVS(2,III)=PVS(2,III)+KI4*DW1
         GO TO 110
   90    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+TABKI(2,POP0)*Z(IZR+I)*Z1
            PVS(2,III)=PVS(2,III)+KI3*Z(IZR+I)*DW1
         ELSE
            PVS(1,III)=PVS(1,III)+(KI3-WI3)*Z1
            PVS(2,III)=PVS(2,III)+(KI4-WI4)*DW1
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
         PVS(1,III)=PVS(1,III)+KI3*Z1
         PVS(2,III)=PVS(2,III)+KI4*DW2
         GO TO 185
  120    K=NINT(POP*PAS3)
         WI3=BI3(K)+POP*(BI31(K)+POP*BI32(K))
         WI4=BI4(K)+POP*(BI41(K)+POP*BI42(K))
         IF(SIGTI.LE.ZERO) THEN
            PVS(1,III)=PVS(1,III)+TABKI(2,POP0)*Z(IZR+NH+1-I)*Z1
            PVS(2,III)=PVS(2,III)+KI3*Z(IZR+NH+1-I)*DW2
         ELSE
            PVS(1,III)=PVS(1,III)+(KI3-WI3)*Z1
            PVS(2,III)=PVS(2,III)+(KI4-WI4)*DW2
         ENDIF
         KI3=WI3
         KI4=WI4
  130    CONTINUE
         KI5=BI5(K)+POP*(BI51(K)+POP*BI52(K))
         PBB(7)=PBB(7)+KI3*Z1
         PBB(9)=PBB(9)+KI4*Z3
         PBB(11)=PBB(11)+KI5*Z1
         PBB(12)=PBB(12)+KI5*Z2
      ELSE IF((ISURF.EQ.5).AND.(JSURF.EQ.2)) THEN
         Z2=Z1*COSINU(1,ISURF)
         Z3=Z1*COSINU(1,ISURF)*COSINU(1,JSURF)
         KI3=ZI30
         KI4=ZI40
         POP=0.0
         K=0
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
         PBB(13)=PBB(13)+2.0*KI3*Z1
         PBB(14)=PBB(14)+2.0*KI4*Z2
         PBB(15)=PBB(15)+2.0*KI5*Z3
         PBB(16)=PBB(16)+2.0*KI5*(Z1-Z3)
      ENDIF
  185 IZR=IZR+NH
  190 CONTINUE
      IZ0=IZ0+NH+4
  200 CONTINUE
  205 CONTINUE
*----
*  PSS ORTHONORMALIZATION
*----
      E1=Z(1)
      E2=Z(2)
      E3=Z(3)
      E4=Z(4)
      P13=PBB(13)
      PBB(1)=PBB(1)*E1*E1
      PBB(7)=PBB(7)*E1*E1
      PBB(13)=P13*E1*E1
      PBB(3)=-0.25*PBB(3)*E1*E4*SQ3
      PBB(9)=-0.25*PBB(9)*E1*E4
      PBB(14)=E1*(PBB(14)*E2-E3*P13)
      SQRT6=.25*SQ3*E4*E2
      PBB(5)=SQRT6*PBB(5)+E3*PBB(3)/E1
      PBB(11)=SQRT6*PBB(11)+E3*PBB(9)/E1
      PBB(6)=-0.5*E4*E4*PBB(6)
      PBB(12)=-0.5*E4*E4*PBB(12)
      PBB(16)=E4*E4*PBB(16)
      PBB(15)=E2*E2*PBB(15)-E3*(2.*PBB(14)/E1+E3*P13)
*----
*  PIS NORMALIZATION
*----
      DO 210 I=1,NREG
      COEF=0.25*SIDE/VOL(I)
      IF(SIGT(I).LE.ZERO) THEN
         PVS(2,I)=COEF*E1*(E2*PVS(2,I)-E3*PVS(1,I))
         PVS(1,I)=COEF*PVS(1,I)*E1*E1
      ELSE
         SIGTI=COEF/SIGT(I)
         PVS(2,I)=SIGTI*E1*(E2*PVS(2,I)-E3*PVS(1,I))
         PVS(1,I)=SIGTI*PVS(1,I)*E1*E1
      ENDIF
  210 CONTINUE
*----
*  OTHER PROBABILITIES COMPUTATION
*----
      PBB(2)=(-SQ3*PBB(3)-4.*PBB(1))/SQ2
      PBB(4)=-4.5*PBB(6)+SQ3*(-SQ2*PBB(5)+8.*PBB(3))+8.*PBB(1)
      PBB(8)=(-3.*SQ3*PBB(9)-4.*PBB(7))/SQ2
      PBB(10)=-4.5*PBB(12)+SQ3*(SQ2*PBB(11)+8.*PBB(9))+8.*PBB(7)
*----
*  TRANSMISSION MATRIX
*----
      DO 225 I=1,18
      DO 220 J=1,18
      IB=IROT(J,I)
      IF(IB.LT.0) THEN
         PSS(I,J)=-PBB(-IB)
      ELSEIF(IB.GT.0) THEN
         PSS(I,J)=PBB(IB)
      ELSE
         PSS(I,J)=0.
      ENDIF
  220 CONTINUE
  225 CONTINUE
*----
*  LEAKAGE MARTIX
*----
      DO 235 J=1,NREG
      K=3
      DO 230 I=1,5
      K=K+3
      PVS(K-1,J)=PVS(2,J)
      PVS(K-2,J)=PVS(1,J)
  230 CONTINUE
  235 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(COSINU)
      RETURN
      END
