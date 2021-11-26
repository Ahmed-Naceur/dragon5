*DECK XELEQN
      SUBROUTINE XELEQN( NDIM, NANGLE, ANGEQN )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Copy generated  angles according to the EQN standard.
*
*Copyright:
* Copyright (C) 1989 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NDIM    number of dimensions (2 or 3).                    
* NANGLE  number of angles.                                 
*
*Parameters: output
* ANGEQN  basis for angles in 2D or 3D.                      
*
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
*
      INTEGER    NDIM, NANGLE
*
      REAL       SN2 ( 1), SN4 ( 2), SN6 ( 4), SN8 ( 6), SN10( 8),
     >           SN12(11), SN14(14), SN16(17), SNT (63),
     >           ANGEQN( NDIM, NDIM), THETA, DTHETA
      INTEGER    MU2 ( 1), MU4 ( 3), MU6 ( 6), MU8 (10), MU10(15),
     >           MU12(21), MU14(28), MU16(36), MUT(120),
     >           ET2 ( 1), ET4 ( 3), ET6 ( 6), ET8 (10), ET10(15),
     >           ET12(21), ET14(28), ET16(36), ETT(120),
     >           XH2 ( 1), XH4 ( 3), XH6 ( 6), XH8 (10), XH10(15),
     >           XH12(21), XH14(28), XH16(36), XHT(120),
     >           INSN( 9), JNMU( 9)
      EQUIVALENCE (SNT( 1), SN2), (SNT( 2), SN4), (SNT( 4), SN6),
     >            (SNT( 8), SN8), (SNT(14),SN10), (SNT(22),SN12),
     >            (SNT(33),SN14), (SNT(47),SN16)
      EQUIVALENCE (MUT( 1), MU2), (MUT( 2), MU4), (MUT( 5), MU6),
     >            (MUT(11), MU8), (MUT(21),MU10), (MUT(36),MU12),
     >            (MUT(57),MU14), (MUT(85),MU16)
      EQUIVALENCE (ETT( 1), ET2), (ETT( 2), ET4), (ETT( 5), ET6),
     >            (ETT(11), ET8), (ETT(21),ET10), (ETT(36),ET12),
     >            (ETT(57),ET14), (ETT(85),ET16)
      EQUIVALENCE (XHT( 1), XH2), (XHT( 2), XH4), (XHT( 5), XH6),
     >            (XHT(11), XH8), (XHT(21),XH10), (XHT(36),XH12),
     >            (XHT(57),XH14), (XHT(85),XH16)
      INTEGER     NANG, NO2LIM, INDEL, NO2, ICUR, IPOS, IEND
      REAL        XPOS, YPOS, ZPOS, X, Y, Z, SUPX, SUPY, SUPZ,
     >            OOSUPX, OOSUPY, OOSUPZ, XOSUPX, YOSUPY, ZOSUPZ
      REAL        PI
      PARAMETER ( PI = 3.1415926535 )
      SAVE
*
      DATA  NANG,  NO2LIM / -1,  8  /
      DATA  INSN/  0,  1,  3,  7, 13, 21, 32, 46, 63/
      DATA  JNMU/  0,  1,  4, 10, 20, 35, 56, 84,120/
*
      DATA  SN2 /  .577350269/
      DATA  SN4 /  .350021174,  .868890300/
      DATA  SN6 /  .2561429  ,  .9320846  ,
     >             .2663443  ,  .6815646  /
      DATA  SN8 /  .1971380  ,  .9603506  ,
     >             .2133981  ,  .5512958  ,  .8065570  ,
     >             .5773503  /
      DATA  SN10/  .1631408  ,  .9730212  ,
     >             .1755273  ,  .6961286  ,
     >                          .4567576  ,  .8721024  ,
     >             .4897749  ,  .7212773  /
      DATA  SN12/  .1370611  ,  .9810344  ,
     >             .1497456  ,  .3911744  ,  .9080522  ,
     >                          .6040252  ,  .7827706  ,
     >             .4213515  ,  .8030727  ,
     >             .4249785  ,  .6400755  /
      DATA  SN14/  .1196230  ,  .9855865  ,
     >             .1301510  ,  .3399238  ,  .9314035  ,
     >                          .5326134  ,  .8362916  ,
     >                          .7010923  ,
     >             .3700559  ,  .8521252  ,
     >             .3736108  ,  .5691823  ,  .7324250  ,
     >             .577350269/
      DATA  SN16/  .1050159  ,  .9889102  ,
     >             .1152880  ,  .3016701  ,  .9464163  ,
     >                          .4743525  ,  .8727534  ,
     >                          .6327389  ,  .7657351  ,
     >             .3284315  ,  .8855877  ,
     >             .3332906  ,  .5107319  ,  .7925089  ,
     >                          .6666774  ,
     >             .5215431  ,  .6752671  /
*
      DATA  MU2 /  1/
      DATA  MU4 /  1,  1,  2/
      DATA  MU6 /  1,  3,  1,  4,  4,  2/
      DATA  MU8 /  1,  3,  3,  1,  4,  6,  4,  5,  5,  2/
      DATA  MU10/  1,  3,  3,  3,  1,  5,  7,  7,  5,  4,
     >             8,  4,  6,  6,  2/
      DATA  MU12/  1,  3,  3,  3,  3,  1,  4,  8, 10,  8,
     >             4,  6, 11, 11,  6,  7,  9,  7,  5,  5,
     >             2/
      DATA  MU14/  1,  3,  3,  3,  3,  3,  1,  4,  9, 11,
     >            11,  9,  4,  6, 12, 14, 12,  6,  8, 13,
     >            13,  8,  7, 10,  7,  5,  5,  2/
      DATA  MU16/  1,  3,  3,  3,  3,  3,  3,  1,  4, 10,
     >            12, 12, 12, 10,  4,  6, 13, 16, 16, 13,
     >             6,  8, 15, 17, 15,  8,  9, 14, 14,  9,
     >             7, 11,  7,  5,  5,  2/
*
      DATA  ET2 /  1/
      DATA  ET4 /  1,  2,  1/
      DATA  ET6 /  1,  4,  2,  3,  4,  1/
      DATA  ET8 /  1,  4,  5,  2,  3,  6,  5,  3,  4,  1/
      DATA  ET10/  1,  5,  4,  6,  2,  3,  7,  8,  6,  3,
     >             7,  4,  3,  5,  1/
      DATA  ET12/  1,  4,  6,  7,  5,  2,  3,  8, 11,  9,
     >             5,  3, 10, 11,  7,  3,  8,  6,  3,  4,
     >             1/
      DATA  ET14/  1,  4,  6,  8,  7,  5,  2,  3,  9, 12,
     >            13, 10,  5,  3, 11, 14, 13,  7,  3, 11,
     >            12,  8,  3,  9,  6,  3,  4,  1/
      DATA  ET16/  1,  4,  6,  8,  9,  7,  5,  2,  3, 10,
     >            13, 15, 14, 11,  5,  3, 12, 16, 17, 14,
     >             7,  3, 12, 16, 15,  9,  3, 12, 13,  8,
     >             3, 10,  6,  3,  4,  1/
*
      DATA  XH2 /  1/
      DATA  XH4 /  2,  1,  1/
      DATA  XH6 /  2,  4,  1,  4,  3,  1/
      DATA  XH8 /  2,  5,  4,  1,  5,  6,  3,  4,  3,  1/
      DATA  XH10/  2,  6,  4,  5,  1,  6,  8,  7,  3,  4,
     >             7,  3,  5,  3,  1/
      DATA  XH12/  2,  5,  7,  6,  4,  1,  5,  9, 11,  8,
     >             3,  7, 11, 10,  3,  6,  8,  3,  4,  3,
     >             1/
      DATA  XH14/  2,  5,  7,  8,  6,  4,  1,  5, 10, 13,
     >            12,  9,  3,  7, 13, 14, 11,  3,  8, 12,
     >            11,  3,  6,  9,  3,  4,  3,  1/
      DATA  XH16/  2,  5,  7,  9,  8,  6,  4,  1,  5, 11,
     >            14, 15, 13, 10,  3,  7, 14, 17, 16, 12,
     >             3,  9, 15, 16, 12,  3,  8, 13, 12,  3,
     >             6, 10,  3,  4,  3,  1/
*
      IF( NDIM.EQ.3 )THEN
         IF( NANGLE.NE.NANG )THEN
            NANG  = NANGLE
            INDEL = 0
            NO2   = NANGLE/2
            IF( NO2.EQ.0 )RETURN
            IF( NO2.LT.1 .OR. NO2.GT.NO2LIM )
     >         CALL XABORT('XELEQN: TOO MANY ANGLES ')
            IPOS  = INSN( NO2 )
            ICUR  = JNMU( NO2 )
            IEND  = JNMU( NO2 + 1)
         ENDIF
         INDEL = INDEL + 1
         IF    ( MOD(INDEL, 3).EQ.1 )THEN
            IF    ( MOD(INDEL, 4).EQ.1 )THEN
               ICUR  = ICUR + 1
               IF( ICUR.GT.IEND )
     >            CALL XABORT('XELEQN: NO MORE ANGLES ')
               XPOS  = SNT( MUT(ICUR) + IPOS )
               YPOS  = SNT( ETT(ICUR) + IPOS )
               ZPOS  = SNT( XHT(ICUR) + IPOS )
               X     = XPOS
               Y     = YPOS
               Z     = ZPOS
               SUPX  = SQRT( 1.0 - X * X )
               SUPY  = SQRT( 1.0 - Y * Y )
               SUPZ  = SQRT( 1.0 - Z * Z )
               OOSUPX= 1.0 / SUPX
               OOSUPY= 1.0 / SUPY
               OOSUPZ= 1.0 / SUPZ
            ELSEIF( MOD(INDEL, 4).EQ.2 )THEN
               X     = -XPOS
               Y     =  YPOS
            ELSEIF( MOD(INDEL, 4).EQ.3 )THEN
               X     =  XPOS
               Y     = -YPOS
            ELSE
               X     = -XPOS
               Y     = -YPOS
            ENDIF
            XOSUPX=  X  / SUPX
            YOSUPY=  Y  / SUPY
            ZOSUPZ=  Z  / SUPZ
*
*           SOLID ANGLE DIRECTION
            ANGEQN( 1, 1 )= X
            ANGEQN( 2, 1 )= Y
            ANGEQN( 3, 1 )= Z
*
*           DIRECTIONS PERPENDICULAR TO THIS SOLID ANGLE
            ANGEQN( 1, 2 )= -Y * OOSUPZ
            ANGEQN( 2, 2 )=  X * OOSUPZ
            ANGEQN( 3, 2 )=         0.0
*
            ANGEQN( 1, 3 )=  X * ZOSUPZ
            ANGEQN( 2, 3 )=  Y * ZOSUPZ
            ANGEQN( 3, 3 )=      - SUPZ
         ELSEIF( MOD(INDEL, 3).EQ.2 )THEN
*
*           SOLID ANGLE DIRECTION
            ANGEQN( 1, 1 )= X
            ANGEQN( 2, 1 )= Y
            ANGEQN( 3, 1 )= Z
*
*           DIRECTIONS PERPENDICULAR TO THIS SOLID ANGLE
            ANGEQN( 1, 2 )= -Z * OOSUPY
            ANGEQN( 2, 2 )=         0.0
            ANGEQN( 3, 2 )=  X * OOSUPY
*
            ANGEQN( 1, 3 )=  X * YOSUPY
            ANGEQN( 2, 3 )=      - SUPY
            ANGEQN( 3, 3 )=  Z * YOSUPY
         ELSE
*
*           SOLID ANGLE DIRECTION
            ANGEQN( 1, 1 )= X
            ANGEQN( 2, 1 )= Y
            ANGEQN( 3, 1 )= Z
*
*           DIRECTIONS PERPENDICULAR TO THIS SOLID ANGLE
            ANGEQN( 1, 2 )=         0.0
            ANGEQN( 2, 2 )= -Z * OOSUPX
            ANGEQN( 3, 2 )=  Y * OOSUPX
*
            ANGEQN( 1, 3 )=      - SUPX
            ANGEQN( 2, 3 )=  Y * XOSUPX
            ANGEQN( 3, 3 )=  Z * XOSUPX
         ENDIF
      ELSEIF( NDIM.EQ.2 )THEN
         IF( NANGLE.NE.NANG )THEN
            NANG  = NANGLE
            IF( NANG.EQ.0 )RETURN
            DTHETA =   PI / NANG
            IF( NANG.GT.0 )THEN
               THETA  = -0.5 * DTHETA
            ELSE
               THETA  =  0.5 * DTHETA
            ENDIF
            INDEL  = 0
         ENDIF
         INDEL = INDEL + 1
         IF( INDEL.GT.NANG ) CALL XABORT( 'XELEQN: NO MORE ANGLES ' )
         THETA = THETA + DTHETA
*
*        SOLID ANGLE DIRECTION
         ANGEQN( 1, 1 )=  COS(THETA)
         ANGEQN( 2, 1 )=  SIN(THETA)
*
*        DIRECTIONS PERPENDICULAR TO THIS SOLID ANGLE
         ANGEQN( 1, 2 )= -SIN(THETA)
         ANGEQN( 2, 2 )=  COS(THETA)
      ELSE
         CALL XABORT( 'XELEQN: *** FALSE NDIM VALUE')
      ENDIF
      RETURN
      END
