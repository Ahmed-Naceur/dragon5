*DECK BIVA05
      SUBROUTINE BIVA05(ITY,SGD,IELEM,NBLOS,LL4,NBMIX,IIMAX,SIDE,MAT,
     1 IPERT,KN,QFR,MU,LC,R,V,H,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a within-group (leakage and removal) or out-of-group
* system matrix in a Thomas-Raviart-Schneider (dual) finite element
* diffusion approximation (hexagonal geometry).
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
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
* SGD     nuclear properties. SGD(:,1) and SGD(:,2) are diffusion
*         coefficients. SGD(:,3) are removal macroscopic cross sections.
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*         =2 (parabolic); =3 (cubic); =4 (quartic).
* NBLOS   number of lozenges per direction, taking into account
*         mesh-splitting.
* LL4     number of unknowns per group in BIVAC.
* NBMIX   number of macro-mixtures.
* IIMAX   allocated dimension of array SYS.
* SIDE    side of the hexagons.
* MAT     mixture index per lozenge.
* IPERT   mixture permutation index.
* KN      element-ordered unknown list.
* QFR     element-ordered boundary conditions.
* MU      indices used with compressed diagonal storage mode matrix SYS.
* LC      order of the unit matrices.
* R       Cartesian mass matrix.
* V       nodal coupling matrix.
* H       Piolat (hexagonal) coupling matrix.
*
*Parameters: output
* SYS     system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,IELEM,NBLOS,LL4,NBMIX,IIMAX,MAT(3,NBLOS),IPERT(NBLOS),
     1 KN(NBLOS,4+6*IELEM*(IELEM+1)),MU(LL4),LC
      REAL SGD(NBMIX,3),SIDE,QFR(NBLOS,6),R(LC,LC),V(LC,LC-1),
     1 H(LC,LC-1),SYS(IIMAX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXIEL=3)
      DOUBLE PRECISION CTRAN(MAXIEL*(MAXIEL+1),MAXIEL*(MAXIEL+1))
*----
*  ASSEMBLY OF A SYSTEM MATRIX.
*----
      TTTT=0.5*SQRT(3.0)*SIDE*SIDE
      IF(IELEM.GT.MAXIEL) CALL XABORT('BIVA05: MAXIEL OVERFLOW.')
      IF(ITY.EQ.0) THEN
*        COMPUTE THE TRANVERSE COUPLING PIOLAT UNIT MATRIX
         CALL XDDSET(CTRAN,MAXIEL**2,0.0D0)
         CNORM=SIDE*SIDE/SQRT(3.0)
         I=0
         DO 22 JS=1,IELEM
         DO 21 JT=1,IELEM+1
         J=0
         I=I+1
         SSS=1.0
         DO 20 IT=1,IELEM
         DO 10 IS=1,IELEM+1
         J=J+1
         CTRAN(I,J)=SSS*CNORM*H(IS,JS)*H(JT,IT)
   10    CONTINUE
         SSS=-SSS
   20    CONTINUE
   21    CONTINUE
   22    CONTINUE
*
*        LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY
         NELEM=IELEM*(IELEM+1)
         COEF=2.0*SIDE*SIDE/SQRT(3.0)
         NUM=0
         DO 70 KEL=1,NBLOS
         IF(IPERT(KEL).EQ.0) GO TO 70
         IBM=MAT(1,IPERT(KEL))
         IF(IBM.EQ.0) GO TO 70
         NUM=NUM+1
         DINV=1.0/SGD(IBM,1)
         SIG=SGD(IBM,3)
         DO 43 K4=0,1
         DO 42 K3=0,IELEM-1
         DO 41 K2=1,IELEM+1
         KNW1=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K2)
         KNX1=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K2)
         KNY1=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K2)
         INW1=ABS(KNW1)
         INX1=ABS(KNX1)
         INY1=ABS(KNY1)
         DO 30 K1=1,IELEM+1
         KNW2=KN(NUM,4+K4*NELEM+K3*(IELEM+1)+K1)
         KNX2=KN(NUM,4+(K4+2)*NELEM+K3*(IELEM+1)+K1)
         KNY2=KN(NUM,4+(K4+4)*NELEM+K3*(IELEM+1)+K1)
         INW2=ABS(KNW2)
         INX2=ABS(KNX2)
         INY2=ABS(KNY2)
         IF((KNW2.NE.0).AND.(KNW1.NE.0).AND.(INW1.GE.INW2)) THEN
            L=MU(INW1)-INW1+INW2
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNW2))
            SYS(L)=SYS(L)-SG*COEF*DINV*R(K2,K1)
            IF(INW1.EQ.INW2) THEN
              IF((K1.EQ.1).AND.(K4.EQ.0)) SYS(L)=SYS(L)-QFR(NUM,1)
              IF((K1.EQ.IELEM+1).AND.(K4.EQ.1)) SYS(L)=SYS(L)-QFR(NUM,2)
            ENDIF
         ENDIF
         IF((KNX2.NE.0).AND.(KNX1.NE.0).AND.(INX1.GE.INX2)) THEN
            L=MU(INX1)-INX1+INX2
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNX2))
            SYS(L)=SYS(L)-SG*COEF*DINV*R(K2,K1)
            IF(INX1.EQ.INX2) THEN
              IF((K1.EQ.1).AND.(K4.EQ.0)) SYS(L)=SYS(L)-QFR(NUM,3)
              IF((K1.EQ.IELEM+1).AND.(K4.EQ.1)) SYS(L)=SYS(L)-QFR(NUM,4)
            ENDIF
         ENDIF
         IF((KNY2.NE.0).AND.(KNY1.NE.0).AND.(INY1.GE.INY2)) THEN
            L=MU(INY1)-INY1+INY2
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNY2))
            SYS(L)=SYS(L)-SG*COEF*DINV*R(K2,K1)
            IF(INY1.EQ.INY2) THEN
              IF((K1.EQ.1).AND.(K4.EQ.0)) SYS(L)=SYS(L)-QFR(NUM,5)
              IF((K1.EQ.IELEM+1).AND.(K4.EQ.1)) SYS(L)=SYS(L)-QFR(NUM,6)
            ENDIF
         ENDIF
   30    CONTINUE
         DO 40 K1=0,IELEM-1
         IF(V(K2,K1+1).EQ.0.0) GO TO 40
         IF(K4.EQ.0) THEN
            SSS=(-1.0)**K1
            JND1=KN(NUM,1)+K3*IELEM+K1
            JND2=KN(NUM,2)+K3*IELEM+K1
            JND3=KN(NUM,3)+K3*IELEM+K1
         ELSE
            SSS=1.0
            JND1=KN(NUM,2)+K1*IELEM+K3
            JND2=KN(NUM,3)+K1*IELEM+K3
            JND3=KN(NUM,4)+K1*IELEM+K3
         ENDIF
         IF(KNW1.NE.0) THEN
            L=MU(JND1)-JND1+INW1
            IF(JND1.LT.INW1) L=MU(INW1)-INW1+JND1
            SG=REAL(SIGN(1,KNW1))
            SYS(L)=SYS(L)+SG*SSS*SIDE*V(K2,K1+1)
         ENDIF
         IF(KNX1.NE.0) THEN
            L=MU(JND2)-JND2+INX1
            IF(JND2.LT.INX1) L=MU(INX1)-INX1+JND2
            SG=REAL(SIGN(1,KNX1))
            SYS(L)=SYS(L)+SG*SSS*SIDE*V(K2,K1+1)
         ENDIF
         IF(KNY1.NE.0) THEN
            L=MU(JND3)-JND3+INY1
            IF(JND3.LT.INY1) L=MU(INY1)-INY1+JND3
            SG=REAL(SIGN(1,KNY1))
            SYS(L)=SYS(L)+SG*SSS*SIDE*V(K2,K1+1)
         ENDIF
   40    CONTINUE
   41    CONTINUE
   42    CONTINUE
   43    CONTINUE
         ITRS=0
         DO I=1,NBLOS
            IF(KN(I,1).EQ.KN(NUM,4)) THEN
               ITRS=I
               GO TO 45
            ENDIF
         ENDDO
         CALL XABORT('BIVA05: ITRS FAILURE.')
   45    DO 55 I=1,NELEM
         KNW1=KN(ITRS,4+I)
         KNX1=KN(NUM,4+2*NELEM+I)
         KNY1=KN(NUM,4+4*NELEM+I)
         INW1=ABS(KNW1)
         INX1=ABS(KNX1)
         INY1=ABS(KNY1)
         DO 50 J=1,NELEM
         KNW2=KN(NUM,4+NELEM+J)
         KNX2=KN(NUM,4+3*NELEM+J)
         KNY2=KN(NUM,4+5*NELEM+J)
         INW2=ABS(KNW2)
         INX2=ABS(KNX2)
         INY2=ABS(KNY2)
         IF((KNY2.NE.0).AND.(KNW1.NE.0).AND.(INW1.LT.INY2)) THEN
            L=MU(INY2)-INY2+INW1
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
            SYS(L)=SYS(L)-SG*DINV*REAL(CTRAN(I,J)) ! y w
         ELSE IF((KNY2.NE.0).AND.(KNW1.NE.0).AND.(INW1.GT.INY2)) THEN
            L=MU(INW1)-INW1+INY2
            SG=REAL(SIGN(1,KNW1)*SIGN(1,KNY2))
            SYS(L)=SYS(L)-SG*DINV*REAL(CTRAN(I,J)) ! w y
         ENDIF
         IF((KNW2.NE.0).AND.(KNX1.NE.0).AND.(INW2.LT.INX1)) THEN
            L=MU(INX1)-INX1+INW2
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
            SYS(L)=SYS(L)-SG*DINV*REAL(CTRAN(I,J)) ! x w
         ELSE IF((KNW2.NE.0).AND.(KNX1.NE.0).AND.(INW2.GT.INX1)) THEN
            L=MU(INW2)-INW2+INX1
            SG=REAL(SIGN(1,KNX1)*SIGN(1,KNW2))
            SYS(L)=SYS(L)-SG*DINV*REAL(CTRAN(I,J)) ! w x
         ENDIF
         IF((KNX2.NE.0).AND.(KNY1.NE.0).AND.(INX2.LT.INY1)) THEN
            L=MU(INY1)-INY1+INX2
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
            SYS(L)=SYS(L)-SG*DINV*REAL(CTRAN(I,J)) ! y x
         ELSE IF((KNX2.NE.0).AND.(KNY1.NE.0).AND.(INX2.GT.INY1)) THEN
            L=MU(INX2)-INX2+INY1
            SG=REAL(SIGN(1,KNY1)*SIGN(1,KNX2))
            SYS(L)=SYS(L)-SG*DINV*REAL(CTRAN(I,J)) ! x y
         ENDIF
   50    CONTINUE
   55    CONTINUE
         DO 65 K2=0,IELEM-1
         DO 60 K1=0,IELEM-1
         JND1=KN(NUM,1)+K2*IELEM+K1
         JND2=KN(NUM,2)+K2*IELEM+K1
         JND3=KN(NUM,3)+K2*IELEM+K1
         L=MU(JND1)
         SYS(L)=SYS(L)+TTTT*SIG
         L=MU(JND2)
         SYS(L)=SYS(L)+TTTT*SIG
         L=MU(JND3)
         SYS(L)=SYS(L)+TTTT*SIG
   60    CONTINUE
   65    CONTINUE
   70    CONTINUE
      ELSE
*        CROSS SECTION SYSTEM MATRIX ASSEMBLY
         NUM=0
         DO 90 KEL=1,NBLOS
         IF(IPERT(KEL).EQ.0) GO TO 90
         IBM=MAT(1,IPERT(KEL))
         IF(IBM.EQ.0) GO TO 90
         NUM=NUM+1
         SIG=SGD(IBM,1)
         DO 85 K2=0,IELEM-1
         DO 80 K1=0,IELEM-1
         JND1=KN(NUM,1)+K2*IELEM+K1
         JND2=KN(NUM,2)+K2*IELEM+K1
         JND3=KN(NUM,3)+K2*IELEM+K1
         L=MU(JND1)
         SYS(L)=SYS(L)+TTTT*SIG
         L=MU(JND2)
         SYS(L)=SYS(L)+TTTT*SIG
         L=MU(JND3)
         SYS(L)=SYS(L)+TTTT*SIG
   80    CONTINUE
   85    CONTINUE
   90    CONTINUE
      ENDIF
      RETURN
      END
