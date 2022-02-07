*DECK BIVA03
      SUBROUTINE BIVA03(ITY,MAXKN,MAXQF,SGD,NREG,LL4,ISPLH,NELEM,NBMIX,
     1 IIMAX,SIDE,MAT,KN,QFR,VOL,MU,R,RH,QH,RT,QT,SYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a within-group (leakage and removal) or out-of-group
* system matrix in mesh-corner finite-difference diffusion
* approximation (hexagonal geometry).
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
* ITY     type of assembly: =0: leakage-removal matrix assembly;
*         =1: cross section matrix assembly.
* MAXKN   dimension of array KN.
* MAXQF   dimension of array QFR.
* SGD     nuclear properties. SGD(:,1) and SGD(:,2) are diffusion
*         coefficients. SGD(:,3) are removal macroscopic cross sections.
* NREG    number of hexagons in BIVAC.
* LL4     order of the matrix SYS.
* ISPLH   hexagonal geometry flag:
*         =1: hexagonal elements; >1: triangular elements.
* NELEM   number of finite elements (hexagons or triangles) excluding
*         the virtual elements.
* NBMIX   number of macro-mixtures.
* IIMAX   allocated dimension of array SYS.
* SIDE    side of the hexagons.
* MAT     mixture index per hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* VOL     volume of the hexagons.
* MU      indices used with the compressed diagonal storage mode matrix
*         SYS.
* R       unit matrix.
* RH      unit matrix.
* QH      unit matrix.
* RT      unit matrix.
* QT      unit matrix.
*
*Parameters: output
* SYS     system matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITY,MAXKN,MAXQF,NREG,LL4,ISPLH,NELEM,NBMIX,IIMAX,
     1 MAT(NREG),KN(MAXKN),MU(LL4)
      REAL SGD(NBMIX,3),SIDE,QFR(MAXQF),VOL(NREG),R(2,2),RH(6,6),
     1 QH(6,6),RT(3,3),QT(3,3),SYS(IIMAX)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION RR,RRH,QQH
      INTEGER ISR(6,2),ISRH(6,2),ISRT(3,2)
      REAL RH2(6,6),QH2(6,6)
      DATA ISRH/2,1,4,5,6,3,1,4,5,6,3,2/
      DATA ISRT/1,2,3,2,3,1/
*----
*  RECOVER THE HEXAGONAL MASS (RH2) AND STIFFNESS (QH2) MATRICES.
*----
      IF(ISPLH.EQ.1) THEN
*        HEXAGONAL BASIS.
         LH=6
         DO 15 I=1,6
         DO 10 J=1,2
         ISR(I,J)=ISRH(I,J)
   10    CONTINUE
   15    CONTINUE
         DO 25 I=1,6
         DO 20 J=1,6
         RH2(I,J)=RH(I,J)
         QH2(I,J)=QH(I,J)
   20    CONTINUE
   25    CONTINUE
         CONST=1.5*SQRT(3.0)
         CONSB=2.0*SQRT(3.0)/3.0
         AA=SIDE
      ELSE
*        TRIANGULAR BASIS.
         LH=3
         DO 35 I=1,3
         DO 30 J=1,2
         ISR(I,J)=ISRT(I,J)
   30    CONTINUE
   35    CONTINUE
         DO 45 I=1,3
         DO 40 J=1,3
         RH2(I,J)=RT(I,J)
         QH2(I,J)=QT(I,J)
   40    CONTINUE
   45    CONTINUE
         CONST=0.25*SQRT(3.0)
         CONSB=2.0*SQRT(3.0)
         AA=SIDE/REAL(ISPLH-1)
      ENDIF
*----
*  ASSEMBLY OF A SYSTEM MATRIX.
*----
      IF(ITY.EQ.0) THEN
*        LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY.
         NUM1=0
         DO 105 K=1,NELEM
         KHEX=KN(NUM1+LH+1)
         IF(VOL(KHEX).EQ.0.0) GO TO 100
         L=MAT(KHEX)
         VOL0=QFR(NUM1+LH+1)
         DO 60 I=1,LH
         IND1=KN(NUM1+I)
         IF(IND1.EQ.0) GO TO 60
         KEY1=MU(IND1)-IND1
         DO 50 J=1,LH
         IND2=KN(NUM1+J)
         IF((IND2.EQ.0).OR.(IND2.GT.IND1)) GO TO 50
         QQH=QH2(I,J)/(CONST*AA*AA)
         RRH=RH2(I,J)/CONST
         IF((QQH.EQ.0.0).AND.(RRH.EQ.0.0)) GO TO 50
         KEY=KEY1+IND2
         SYS(KEY)=SYS(KEY)+REAL(QQH*SGD(L,1)+RRH*SGD(L,3))*VOL0
   50    CONTINUE
   60    CONTINUE
         DO 90 IC=1,LH
         QFR1=QFR(NUM1+IC)
         IF(QFR1.EQ.0.0) GO TO 90
         DO 80 I1=1,2
         IND1=KN(NUM1+ISR(IC,I1))
         IF(IND1.EQ.0) GO TO 80
         KEY1=MU(IND1)-IND1
         DO 70 J1=1,2
         IND2=KN(NUM1+ISR(IC,J1))
         IF((IND2.EQ.0).OR.(IND2.GT.IND1)) GO TO 70
         RR=R(I1,J1)
         IF(RR.EQ.0.0) GO TO 70
         KEY=KEY1+IND2
         SYS(KEY)=SYS(KEY)+REAL(RR)*QFR1
   70    CONTINUE
   80    CONTINUE
   90    CONTINUE
  100    NUM1=NUM1+LH+1
  105    CONTINUE
      ELSE
*        CROSS SECTION SYSTEM MATRIX ASSEMBLY
         NUM1=0
         DO 135 K=1,NELEM
         KHEX=KN(NUM1+LH+1)
         IF(VOL(KHEX).EQ.0.0) GO TO 130
         L=MAT(KHEX)
         VOL0=QFR(NUM1+LH+1)
         DO 120 I=1,LH
         IND1=KN(NUM1+I)
         IF(IND1.EQ.0) GO TO 120
         KEY1=MU(IND1)-IND1
         DO 110 J=1,LH
         IND2=KN(NUM1+J)
         IF((IND2.EQ.0).OR.(IND2.GT.IND1)) GO TO 110
         RRH=RH2(I,J)/CONST
         IF(RRH.EQ.0.0) GO TO 110
         KEY=KEY1+IND2
         SYS(KEY)=SYS(KEY)+REAL(RRH)*SGD(L,1)*VOL0
  110    CONTINUE
  120    CONTINUE
  130    NUM1=NUM1+LH+1
  135    CONTINUE
      ENDIF
      RETURN
      END
