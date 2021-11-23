*DECK ALGJP
      SUBROUTINE ALGJP(NGPT,ZJKSI,WJKSI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* returns Gauss-Jacobi integration points and weights for integration
* between 0.0 and 1.0 to the order specified
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Reference:
* integral of X F(X) DX , from 0.0 to 1.0 is given by summ from i=1 to
* NGPT of F(ZJKSI(I))*WJKSI(I)
* the points ZGKSI and weights WJKSI are generated according
* to the technique described in Handbook of mathematical functions
* M. Abramowitz and I. Stegun, Dover Publication Inc. (1972).
*
*Parameters: input
* NGPT    number of gauss-jacobi points
*
*Parameters: output
* ZGKSI   integration points
* WGKSI   integration weights
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGPT
      REAL ZJKSI(NGPT),WJKSI(NGPT)
*----
*  LOCAL VARIABLES
*----
      PARAMETER ( I1= 1, I2= 2, I3= 4, I4= 7, I5=11, I6=16,
     >            I7=22, I8=29,IFN=36)
      REAL       RZKSI(IFN),RWKSI(IFN)
C  N=1
      DATA (RZKSI(I),I=I1,I2-1)   /0.6666666667/
      DATA (RWKSI(I),I=I1,I2-1)   /0.5000000000/
C  N=2
      DATA (RZKSI(I),I=I2,I3-1)   /0.3550510257,0.8449489743/
      DATA (RWKSI(I),I=I2,I3-1)   /0.1819586183,0.3180413817/
C  N=3
      DATA (RZKSI(I),I=I3,I4-1)   /.2123405382,.5905331356,.9114120405/
      DATA (RWKSI(I),I=I3,I4-1)   /.0698269799,.2292411064,.2009319137/
C  N=4
      DATA (RZKSI(I),I=I4,I5-1)
     >            /.1397598643,.4164095676,.7231569864,.9428958039/
      DATA (RWKSI(I),I=I4,I5-1)
     >            /.0311809710,.1298475476,.2034645680,.1355069134/
C  N=5
      DATA (RZKSI(I),I=I5,I6-1)   /.0985350858,.3045357266,
     >                             .5620251898,.8019865821,.9601901429/
      DATA (RWKSI(I),I=I5,I6-1)   /.0157479145,.0739088701,
     >                             .1463869871,.1671746381,.0967815902/
C  N=6
      DATA (RZKSI(I),I=I6,I7-1)   /.0730543287,.2307661380,.4413284812,
     >                             .6630153097,.8519214003,.9706835728/
      DATA (RWKSI(I),I=I6,I7-1)   /.0087383018,.0439551656,.0986611509,
     >                             .1407925538,.1355424972,.0723103307/
C  N=7
      DATA (RZKSI(I),I=I7,I8-1)   /.0562625605,.1802406917,.3526247171,
     >                 .5471536263,.7342101772,.8853209468,.9775206136/
      DATA (RWKSI(I),I=I7,I8-1)   /.0052143622,.0274083567,.0663846965,
     >                 .1071250657,.1273908973,.1105092582,.0559673634/
C  N=8
      DATA (RZKSI(I),I=I8,IFN)    /.0446339553,.1443662570,.2868247571,
     >     .4548133152,.6280678354,.7856915206,.9086763921,.9822200949/
      DATA (RWKSI(I),I=I8,IFN)    /.0032951914,.0178429027,.0454393195,
     >     .0791995995,.1060473594,.1125057995,.0911190236,.0445508044/
*
      IDEP=0
      IFIN=0
      IF( NGPT.EQ. 1 ) THEN
        IDEP=I1
        IFIN=I2-1
      ELSE IF( NGPT.EQ. 2 ) THEN
        IDEP=I2
        IFIN=I3-1
      ELSE IF( NGPT.EQ. 3 ) THEN
        IDEP=I3
        IFIN=I4-1
      ELSE IF( NGPT.EQ. 4 ) THEN
        IDEP=I4
        IFIN=I5-1
      ELSE IF( NGPT.EQ. 5 ) THEN
        IDEP=I5
        IFIN=I6-1
      ELSE IF( NGPT.EQ. 6 ) THEN
        IDEP=I6
        IFIN=I7-1
      ELSE IF( NGPT.EQ. 7 ) THEN
        IDEP=I7
        IFIN=I8-1
      ELSE IF( NGPT.EQ. 8 ) THEN
        IDEP=I8
        IFIN=IFN
      ELSE
        XINF=0.0
        XSUP=1.0
        CALL ALGPT(NGPT,XINF,XSUP,ZJKSI,WJKSI)
      ENDIF
      IF(NGPT.LE.8) THEN
C------
C  INITIALIZE ZJKSI AND WJKSI FROM DATA BASE
C------
        IUP=1
        DO 100 I=IDEP,IFIN
          ZJKSI(IUP)=RZKSI(I)
          WJKSI(IUP)=RWKSI(I)
          IUP=IUP+1
 100    CONTINUE
      ELSE
C------
C  USE GAUSS-LEGENDRE INTEGRATION POINTS INSTEAD OF GAUSS-JACOBI
C------
        DO 110 I=1,NGPT
          WJKSI(I)=WJKSI(I)*ZJKSI(I)
 110    CONTINUE
      ENDIF
      RETURN
      END
