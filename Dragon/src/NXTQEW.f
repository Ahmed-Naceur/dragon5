*DECK NXTQEW
      SUBROUTINE NXTQEW(NDIM  ,NANGL ,NQUAD ,NBANGL,DQUAD,
     >                  DANGLT,DDENWT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define quadrature angles for a given tracking option.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
*  G. Marleau, R. Roy, M. Hampartzounian
*
*Parameters: input
* NDIM    number of dimensions for geometry.
* NANGL   quadrature order.
* NQUAD   number of quadrant (in 3-D) and quarter (in 2-D).
* NBANGL  number of angles.
* DQUAD   relative density of each quadrant.
*
*Parameters: output
* DANGLT  director cosines of angles.
* DDENWT  angular density for each angle.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  Extracted from the subroutine XELTS2 of EXCELL.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NDIM,NANGL,NQUAD,NBANGL
      DOUBLE PRECISION DQUAD(NQUAD)
      DOUBLE PRECISION DANGLT(NDIM,NQUAD,NBANGL),DDENWT(NQUAD,NBANGL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTQEW')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IANG,NO2,IPOS,ICUR,IEND
      DOUBLE PRECISION DDA,X,Y,Z
*----
*  Data
*----
      INTEGER          INSN( 9),JNMU( 9),MUT(120), ETT(120),XHT(120)
      REAL             SNT (63)
      SAVE             INSN,JNMU,MUT,ETT,XHT,SNT
      DATA             INSN/  0,  1,  3,  7, 13, 21, 32, 46, 63/
      DATA             JNMU/  0,  1,  4, 10, 20, 35, 56, 84,120/
      DATA  MUT /  1,
     >             1,  1,  2,
     >             1,  3,  1,  4,  4,  2,
     >             1,  3,  3,  1,  4,  6,  4,  5,  5,  2,
     >             1,  3,  3,  3,  1,  5,  7,  7,  5,  4,
     >             8,  4,  6,  6,  2,
     >             1,  3,  3,  3,  3,  1,  4,  8, 10,  8,
     >             4,  6, 11, 11,  6,  7,  9,  7,  5,  5,
     >             2,
     >             1,  3,  3,  3,  3,  3,  1,  4,  9, 11,
     >            11,  9,  4,  6, 12, 14, 12,  6,  8, 13,
     >            13,  8,  7, 10,  7,  5,  5,  2,
     >             1,  3,  3,  3,  3,  3,  3,  1,  4, 10,
     >            12, 12, 12, 10,  4,  6, 13, 16, 16, 13,
     >             6,  8, 15, 17, 15,  8,  9, 14, 14,  9,
     >             7, 11,  7,  5,  5,  2/
      DATA  ETT /  1,
     >             1,  2,  1,
     >             1,  4,  2,  3,  4,  1,
     >             1,  4,  5,  2,  3,  6,  5,  3,  4,  1,
     >             1,  5,  4,  6,  2,  3,  7,  8,  6,  3,
     >             7,  4,  3,  5,  1,
     >             1,  4,  6,  7,  5,  2,  3,  8, 11,  9,
     >             5,  3, 10, 11,  7,  3,  8,  6,  3,  4,
     >             1,
     >             1,  4,  6,  8,  7,  5,  2,  3,  9, 12,
     >            13, 10,  5,  3, 11, 14, 13,  7,  3, 11,
     >            12,  8,  3,  9,  6,  3,  4,  1,
     >             1,  4,  6,  8,  9,  7,  5,  2,  3, 10,
     >            13, 15, 14, 11,  5,  3, 12, 16, 17, 14,
     >             7,  3, 12, 16, 15,  9,  3, 12, 13,  8,
     >             3, 10,  6,  3,  4,  1/
      DATA  XHT /  1,
     >             2,  1,  1,
     >             2,  4,  1,  4,  3,  1,
     >             2,  5,  4,  1,  5,  6,  3,  4,  3,  1,
     >             2,  6,  4,  5,  1,  6,  8,  7,  3,  4,
     >             7,  3,  5,  3,  1,
     >             2,  5,  7,  6,  4,  1,  5,  9, 11,  8,
     >             3,  7, 11, 10,  3,  6,  8,  3,  4,  3,
     >             1,
     >             2,  5,  7,  8,  6,  4,  1,  5, 10, 13,
     >            12,  9,  3,  7, 13, 14, 11,  3,  8, 12,
     >            11,  3,  6,  9,  3,  4,  3,  1,
     >             2,  5,  7,  9,  8,  6,  4,  1,  5, 11,
     >            14, 15, 13, 10,  3,  7, 14, 17, 16, 12,
     >             3,  9, 15, 16, 12,  3,  8, 13, 12,  3,
     >             6, 10,  3,  4,  3,  1/
      DATA             SNT /  .577350269,
     > .350021174,  .868890300,
     > .2561429,  .9320846, .2663443,  .6815646,
     > .1971380,  .9603506, .2133981,  .5512958,  .8065570,  .5773503,
     > .1631408,  .9730212, .1755273,  .6961286,  .4567576,  .8721024,
     > .4897749,  .7212773, .1370611,  .9810344,  .1497456,  .3911744,
     > .9080522,  .6040252, .7827706,  .4213515,  .8030727,  .4249785,
     > .6400755,
     > .1196230,  .9855865, .1301510,  .3399238,  .9314035,  .5326134,
     > .8362916,  .7010923, .3700559,  .8521252,  .3736108,  .5691823,
     > .7324250,  .577350269,
     > .1050159,  .9889102, .1152880,  .3016701,  .9464163,  .4743525,
     > .8727534,  .6327389, .7657351,  .3284315,  .8855877,  .3332906,
     > .5107319,  .7925089, .6666774,  .5215431,  .6752671/
*----
*  Start processing
*----
      DDA=DBLE(8*NBANGL)
      NO2   = NANGL/2
      IPOS  = INSN( NO2 )
      ICUR  = JNMU( NO2 )
      IEND  = JNMU( NO2 + 1)
      DO IANG=1,NBANGL
        ICUR=ICUR+1
        X  = DBLE(SNT( MUT(ICUR) + IPOS ))
        Y  = DBLE(SNT( ETT(ICUR) + IPOS ))
        Z  = DBLE(SNT( XHT(ICUR) + IPOS ))
        DANGLT(1,1,IANG)=X
        DANGLT(2,1,IANG)=Y
        DANGLT(3,1,IANG)=Z
        DDENWT(1,IANG)=DQUAD(1)*DDA
        DANGLT(1,2,IANG)=-X
        DANGLT(2,2,IANG)=Y
        DANGLT(3,2,IANG)=Z
        DDENWT(2,IANG)=DQUAD(2)*DDA
        DANGLT(1,3,IANG)=X
        DANGLT(2,3,IANG)=-Y
        DANGLT(3,3,IANG)=Z
        DDENWT(3,IANG)=DQUAD(3)*DDA
        DANGLT(1,4,IANG)=-X
        DANGLT(2,4,IANG)=-Y
        DANGLT(3,4,IANG)=Z
        DDENWT(4,IANG)=DQUAD(4)*DDA
      ENDDO
*----
*  Processing finished: return
*----
      RETURN
*----
*  Output formats
*----
      END
