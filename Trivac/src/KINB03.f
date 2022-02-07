*DECK KINB03
      SUBROUTINE KINB03(MAXKN,MAXQF,SGD,NREG,LL4,ISPLH,NELEM,NBMIX,
     1 MAT,KN,QFR,VOL,RH,RT,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mesh-corner finite-
* difference diffusion approximation (hexagonal geometry). Special
* version for Bivac.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXKN   dimension of array KN.
* MAXQF   dimension of array QFR.
* SGD     mixture-ordered cross sections.
* NREG    number of hexagons in Bivac.
* LL4     number of unknowns (order of the system matrices).
* ISPLH   hexagonal geometry flag:
*         =1: hexagonal elements; >1: triangular elements.
* NELEM   number of finite elements (hexagons or triangles) excluding
*         the virtual elements.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* VOL     volume of the hexagons.
* RH      unit matrix.
* RT      unit matrix.
* F2      vector to multiply.
*
*Parameters: output
* F3      result of the multiplication.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXKN,MAXQF,NREG,LL4,ISPLH,NELEM,NBMIX,MAT(NREG),KN(MAXKN)
      REAL SGD(NBMIX),QFR(MAXQF),VOL(NREG),RH(6,6),RT(3,3),F2(LL4),
     1 F3(LL4)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION RRH
      INTEGER ISR(6,2),ISRH(6,2),ISRT(3,2)
      REAL RH2(6,6)
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
   20    CONTINUE
   25    CONTINUE
         CONST=1.5*SQRT(3.0)
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
   40    CONTINUE
   45    CONTINUE
         CONST=0.25*SQRT(3.0)
      ENDIF
*----
*  MULTIPLICATION
*----
      NUM1=0
      DO 80 K=1,NELEM
      KHEX=KN(NUM1+LH+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 70
      L=MAT(KHEX)
      VOL0=QFR(NUM1+LH+1)
      DO 60 I=1,LH
      IND1=KN(NUM1+I)
      IF(IND1.EQ.0) GO TO 60
      DO 50 J=1,LH
      IND2=KN(NUM1+J)
      IF(IND2.EQ.0) GO TO 50
      RRH=RH2(I,J)/CONST
      IF(RRH.EQ.0.0) GO TO 50
      F3(IND1)=F3(IND1)+REAL(RRH)*SGD(L)*VOL0*F2(IND2)
   50 CONTINUE
   60 CONTINUE
   70 NUM1=NUM1+LH+1
   80 CONTINUE
      RETURN
      END
