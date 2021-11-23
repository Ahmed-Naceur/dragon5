*DECK ALCACT
      SUBROUTINE ALCACT(LCACT,NG,XG,WG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* set the weights and base points for the different polar quadrature:
* - "Cactus" (Halsall)
* - "optimized" (Leonard and extension by Le Tellier)
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and R. Le Tellier
*
*Reference:
* A. Leonard and C.T. Mc Daniel, "Optimal polar angles and weights for
* the characteristic method", Trans. Am. Nucl. Soc., 73, 172 (1995).
*
*Parameters: input
* LCACT   type of quadrature (=1,2: values used by Halsall in
*         Cactus (1980); =3: values proposed by Mc Daniel ng=2 only,
*         extended to 3 and 4 with the same approach by R. Le Tellier
*         in 06/2005)
* NG      number of weights/base points.
*
*Parameters: output
* XG      base points.
* WG      weights.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LCACT,NG
      REAL XG(NG),WG(NG)
*----
*  LOCAL VARIABLES
*----
      INTEGER  IG
      DOUBLE PRECISION PI, PHIM, PHIP, DELPHI, WTEST, XTEST, YTEST,
     >                 ZTEST, ZERO, ONE, HALF, QUART
      PARAMETER( PI=3.14159265358979D0, ZERO=0.D0, ONE=1.D0,
     >           HALF=0.5D0, QUART=0.25D0 )
* For "Optimized" Quadrature
      INTEGER I, I2, I3, I4, I5, IDEP, IFIN
      PARAMETER ( I2= 1, I3= 3, I4= 6, I5= 10)
      REAL     XOP0(I5-1), WOP0(I5-1) ! Quadrature which minimizes the error on Ki2 without constraints.
      REAL     XOP1(I5-1), WOP1(I5-1) ! Quadrature which minimizes the error on Ki2 with P1 constraints.
      REAL     XGOP(I5-1), WGOP(I5-1) ! Gauss optimized quadrature.
      SAVE XOP0, WOP0, XOP1, WOP1, XGOP, WGOP
* NG=2
      DATA (XOP0(I),I=I2,I3-1) / 0.273658, 0.865714/
      DATA (WOP0(I),I=I2,I3-1) / 0.139473, 0.860527/
      DATA (XOP1(I),I=I2,I3-1) / 0.340183, 0.894215/
      DATA (WOP1(I),I=I2,I3-1) / 0.194406, 0.805594/
      DATA (XGOP(I),I=I2,I3-1) / 0.399374, 0.914448/
      DATA (WGOP(I),I=I2,I3-1) / 0.250547, 0.749453/
* NG=3
      DATA (XOP0(I),I=I3,I4-1) / 0.891439, 0.395534, 0.099812/
      DATA (WOP0(I),I=I3,I4-1) / 0.793820, 0.188560, 0.017620/
      DATA (XOP1(I),I=I3,I4-1) / 0.131209, 0.478170, 0.920079/
      DATA (WOP1(I),I=I3,I4-1) / 0.029991, 0.250860, 0.719149/
      DATA (XGOP(I),I=I3,I4-1) / 0.231156, 0.639973, 0.954497/
      DATA (WGOP(I),I=I3,I4-1) / 0.085302, 0.341456, 0.573242/
* NG=4
      DATA (XOP0(I),I=I4,I5-1) / 0.464167, 0.908274, 0.166004, 0.042181/
      DATA (WOP0(I),I=I4,I5-1) / 0.218331, 0.746430, 0.032141, 0.003098/
      DATA (XOP1(I),I=I4,I5-1) / 0.054819, 0.212313, 0.546065, 0.932318/
      DATA (WOP1(I),I=I4,I5-1) / 0.005225, 0.051270, 0.272789, 0.670716/
      DATA (XGOP(I),I=I4,I5-1) / 0.152641, 0.450820, 0.769181, 0.972320/
      DATA (WGOP(I),I=I4,I5-1) / 0.037508, 0.167623, 0.338496, 0.456373/
!      DATA (XGOP(I),I=I5,I6-1) / 0.159153, 0.450941, 0.724750, 0.868405,
!     1                           0.980414/
!      DATA (WGOP(I),I=I5,I6-1) / 0.039977, 0.156595, 0.220669, 0.204037,
!     1                           0.378722/
*
      IF( LCACT.EQ.1 )THEN
*---
* CACTUS 1:
*---
*        Equal weight quadrature
         DELPHI= ONE/DBLE(NG)
         PHIM  = ZERO
         XTEST = ZERO
         DO 10 IG= 1, NG
            PHIP  = DBLE(IG) * DELPHI
            WTEST = PHIP - PHIM
            WG(IG)= REAL( WTEST )
            IF(IG.EQ.NG) THEN
               YTEST = PI/2.0
            ELSE
               YTEST = SQRT(ONE - PHIP * PHIP) * PHIP + ASIN(PHIP)
            ENDIF
            ZTEST = HALF * (YTEST - XTEST) / WTEST
            XG(IG)= REAL( SQRT(ONE - ZTEST*ZTEST) )
            PHIM  = PHIP
            XTEST = YTEST
 10      CONTINUE
      ELSEIF( LCACT.EQ.2 )THEN
*---
* CACTUS 2:
*---
*        Uniformly distributed angle quadrature
         DELPHI= PI/DBLE(2*NG)
         PHIM  = ZERO
         DO 20 IG= 1, NG
            PHIP  = DBLE(IG) * DELPHI
            WTEST = SIN(PHIP) - SIN(PHIM)
            WG(IG)= REAL( WTEST )
            XTEST = HALF * (PHIP - PHIM)
            YTEST = QUART * (SIN(PHIP+PHIP) - SIN(PHIM+PHIM))
            ZTEST = (XTEST + YTEST) / WTEST
            XG(IG)= REAL( SQRT(ONE - ZTEST*ZTEST) )
            PHIM= PHIP
 20      CONTINUE
      ELSEIF(( LCACT.GE.3 ).AND.( LCACT.LE.5 ))THEN
*---
* OPTIMIZED:
*---
         IDEP=0
         IFIN=0
         IF( NG.EQ.2 ) THEN
            IDEP=I2
            IFIN=I3-1
         ELSE IF( NG.EQ.3 ) THEN
            IDEP=I3
            IFIN=I4-1
         ELSE IF( NG.EQ.4 ) THEN
            IDEP=I4
            IFIN=I5-1
         ELSE
            CALL XABORT('ALCACT: LCACA=3 => NG= 2, 3 OR 4')
         ENDIF
         IF (LCACT.EQ.3) THEN
*        Quadrature which minimizes the error on Ki2 without constraints.   
            DO 30 IG= IDEP, IFIN
               XG(IG-IDEP+1)= REAL(SQRT(ONE - XOP0(IG)*XOP0(IG)))
               WG(IG-IDEP+1)= WOP0(IG)
 30         CONTINUE          
         ELSEIF (LCACT.EQ.4) THEN
*        Quadrature which minimizes the error on Ki2 with P1 constraints.
            DO 40 IG= IDEP, IFIN
               XG(IG-IDEP+1)= REAL(SQRT(ONE - XOP1(IG)*XOP1(IG)))
               WG(IG-IDEP+1)= WOP1(IG)
 40         CONTINUE 
         ELSEIF (LCACT.EQ.5) THEN
*        Gauss optimized quadrature.
            DO 50 IG= IDEP, IFIN
               XG(IG-IDEP+1)= REAL(SQRT(ONE - XGOP(IG)*XGOP(IG)))
               WG(IG-IDEP+1)= WGOP(IG)
 50         CONTINUE 
         ENDIF
      ELSE
         CALL XABORT('ALCACT: *LCACT* IN [1,5]') 
      ENDIF
      RETURN
      END
