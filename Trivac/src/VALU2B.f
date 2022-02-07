*DECK VALU2B
      SUBROUTINE VALU2B (LC,MKN,LX,LY,L4,X,Y,XXX,YYY,EVECT,ISS,KN,IXLG,
     + IYLG,E,AXY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution for PRIM method in 2D.
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
* LC      order of the unit matrices.
* MKN     second dimension for matrix KN.
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* L4      dimension of unknown array EVECT.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* Y       Cartesian coordinates along the Y axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* EVECT   variational coefficients of the flux.
* ISS     mixture index assigned to each element.
* KN      element-ordered unknown list.
* IELEM   MCFD polynomial order (IELEM=1 is the mesh centered finite
*         difference method).
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.
* E       Lagrange polynomial coefficients.
*                                                                      
*Parameters: output
* AXY     interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LC,MKN,LX,LY,L4,ISS(LX*LY),KN(LX*LY*MKN),IXLG,IYLG
      REAL X(IXLG),Y(IYLG),XXX(LX+1),YYY(LY+1),EVECT(L4),AXY(IXLG,IYLG),
     1 E(LC,LC)
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(125),IJ2(125)
      REAL FLX(5),FLY(5)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
      REAL, ALLOCATABLE, DIMENSION(:,:) ::COEF
*----
*  Scratch storage allocation
*----
      ALLOCATE(IWRK(LX*LY),COEF(LX*LY,MKN))
*----
*  Calculation of IJ integer arrays
*----
      LL=LC*LC
      DO 5 L=1,LL
      L1=1+MOD(L-1,LC)
      L2=1+(L-L1)/LC
      L3=1+MOD(L2-1,LC)
      IJ1(L)=L1
      IJ2(L)=L3
    5 CONTINUE
*
      NUM=0
      DO 10 I=1,LX*LY
      IWRK(I)=0
      IF (ISS(I).EQ.0) GO TO 10
      IWRK(I)=NUM
      NUM=NUM+1
  10  CONTINUE
*
      DO 110 J=1,IYLG
      ORDO=Y(J)
      DO 100 I=1,IXLG
      ABSC=X(I)
      AXY(I,J)=0.0
*                                                          
*     Find the finite element index containing the interpolation point
      IS=0
      JS=0
      DO 20 L=1,LX
      IS=L
      IF((ABSC.GE.XXX(L)).AND.(ABSC.LE.XXX(L+1))) GO TO 30
   20 CONTINUE
      CALL XABORT('VALU2B: WRONG INTERPOLATION(1).')
   30 DO 40 L=1,LY
      JS=L
      IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 70
   40 CONTINUE
      CALL XABORT('VALU2B: WRONG INTERPOLATION(2).')
   70 IEL=(JS-1)*LX+IS
*
      IF(ISS(IEL).EQ.0) GO TO 100
      NUM=IWRK(IEL)
      IF (NUM.NE.-1) THEN
        DO 85 M=1,LL
        I1=IJ1(M)
        I2=IJ2(M)
        COEF(IEL,M)=0.0
        DO 80 N=1,LL
        IND2=KN(LL*NUM+N)
        IF (IND2.EQ.0) GO TO 80
        J1=IJ1(N)
        J2=IJ2(N)
        COEF(IEL,M)=COEF(IEL,M)+E(I1,J1)*E(I2,J2)*EVECT(IND2)
   80   CONTINUE
   85   CONTINUE
        IWRK(IEL)=-1
      ENDIF
*
      U=(ABSC-0.5*(XXX(IS)+XXX(IS+1)))/(XXX(IS+1)-XXX(IS))
      FLX(1)=1.0
      FLX(2)=FLX(1)*U
      FLX(3)=FLX(2)*U
      FLX(4)=FLX(3)*U
      FLX(5)=FLX(4)*U
      V=(ORDO-0.5*(YYY(JS)+YYY(JS+1)))/(YYY(JS+1)-YYY(JS))
      FLY(1)=1.0
      FLY(2)=FLY(1)*V
      FLY(3)=FLY(2)*V
      FLY(4)=FLY(3)*V
      FLY(5)=FLY(4)*V
      DO 90 L=1,LL
      I1=IJ1(L)
      I2=IJ2(L)
      AXY(I,J)=AXY(I,J)+COEF(IEL,L)*FLX(I1)*FLY(I2)
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(COEF,IWRK)
      RETURN
      END
