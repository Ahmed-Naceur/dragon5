*DECK KINB04
      SUBROUTINE KINB04(MAXKN,MAXQF,SGD,NREG,LL4,ISPLH,NBMIX,MAT,KN,
     1 QFR,VOL,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multiplication of a matrix by a vector in mesh-centered finite-
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
* LL4     number of unknowns per group in Bivac. Equal to the number
*         of finite elements (hexagons or triangles) excluding the
*         virtual elements.
* ISPLH   type of hexagonal mesh-splitting:
*         =1: hexagonal elements; >1: triangular elements.
* NBMIX   number of macro-mixtures.
* MAT     mixture index per hexagon.
* KN      element-ordered unknown list.
* QFR     element-ordered information.
* VOL     volume of hexagons.
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
      INTEGER MAXKN,MAXQF,NREG,LL4,ISPLH,NBMIX,MAT(NREG),KN(MAXKN)
      REAL SGD(NBMIX),QFR(MAXQF),VOL(NREG),F2(LL4),F3(LL4)
*
      IF(ISPLH.EQ.1) THEN
         NSURF=6
      ELSE
         NSURF=3
      ENDIF
*----
*  MULTIPLICATION.
*----
      NUM1=0
      DO 20 IND1=1,LL4
      KHEX=KN(NUM1+NSURF+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 10
      L=MAT(KHEX)
      F3(IND1)=F3(IND1)+SGD(L)*QFR(NUM1+NSURF+1)*F2(IND1)
   10 NUM1=NUM1+NSURF+1
   20 CONTINUE
      RETURN
      END
