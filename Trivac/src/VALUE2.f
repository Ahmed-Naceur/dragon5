*DECK VALUE2
      SUBROUTINE VALUE2 (LC,MKN,LX,LY,LZ,L4,X,Y,Z,XXX,YYY,ZZZ,EVECT,
     + ISS,KN,IXLG,IYLG,IZLG,E,AXYZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution for PRIM method in 3D.
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
* LZ      number of elements along the Z axis.
* L4      dimension of unknown array EVECT.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* Y       Cartesian coordinates along the Y axis where the flux is
*         interpolated.
* Z       Cartesian coordinates along the Z axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
* EVECT   variational coefficients of the flux.
* ISS     mixture index assigned to each element.
* KN      element-ordered unknown list.
* IELEM   MCFD polynomial order (IELEM=1 is the mesh centered finite
*         difference method).
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.
* IZLG    number of interpolated points according to Z.
* E       Lagrange polynomial coefficients.
*                                                                      
*Parameters: output
* AXYZ    interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER LC,MKN,LX,LY,LZ,L4,ISS(LX*LY*LZ),KN(LX*LY*LZ*MKN),IXLG,
     1 IYLG,IZLG
      REAL X(IXLG),Y(IYLG),Z(IZLG),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),
     1 EVECT(L4),AXYZ(IXLG,IYLG,IZLG),E(LC,LC)
*----
*  LOCAL VARIABLES
*----
      INTEGER IJ1(125),IJ2(125),IJ3(125)
      REAL FLX(5),FLY(5),FLZ(5)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
      REAL, ALLOCATABLE, DIMENSION(:,:) ::COEF
*----
*  Scratch storage allocation
*----
      ALLOCATE(IWRK(LX*LY*LZ),COEF(LX*LY*LZ,MKN))
*----
*  Calculation of IJ integer arrays
*----
      LL=LC*LC*LC
      DO 5 L=1,LL
      L1=1+MOD(L-1,LC)
      L2=1+(L-L1)/LC
      L3=1+MOD(L2-1,LC)
      IJ1(L)=L1
      IJ2(L)=L3
      IJ3(L)=1+(L2-L3)/LC
    5 CONTINUE
*
      NUM=0
      DO 10 I=1,LX*LY*LZ
      IWRK(I)=0
      IF (ISS(I).EQ.0) GO TO 10
      IWRK(I)=NUM
      NUM=NUM+1
  10  CONTINUE
*
      DO 120 K=1,IZLG
      COTE=Z(K)
      DO 110 J=1,IYLG
      ORDO=Y(J)
      DO 100 I=1,IXLG
      ABSC=X(I)
      AXYZ(I,J,K)=0.0
*                                                          
*     Find the finite element index containing the interpolation point
      IS=0
      JS=0
      KS=0
      DO 20 L=1,LX
      IS=L
      IF((ABSC.GE.XXX(L)).AND.(ABSC.LE.XXX(L+1))) GO TO 30
   20 CONTINUE
      CALL XABORT('VALUE2: WRONG INTERPOLATION(1).')
   30 DO 40 L=1,LY
      JS=L
      IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 50
   40 CONTINUE
      CALL XABORT('VALUE2: WRONG INTERPOLATION(2).')
   50 DO 60 L=1,LZ
      KS=L
      IF((COTE.GE.ZZZ(L)).AND.(COTE.LE.ZZZ(L+1))) GO TO 70
   60 CONTINUE
      CALL XABORT('VALUE2: WRONG INTERPOLATION(3).')
   70 IEL=(KS-1)*LX*LY+(JS-1)*LX+IS
*
      IF(ISS(IEL).EQ.0) GO TO 100
      NUM=IWRK(IEL)
      IF (NUM.NE.-1) THEN
        DO 85 M=1,LL
        I1=IJ1(M)
        I2=IJ2(M)
        I3=IJ3(M)
        COEF(IEL,M)=0.0
        DO 80 N=1,LL
        IND2=KN(LL*NUM+N)
        IF (IND2.EQ.0) GO TO 80
        J1=IJ1(N)
        J2=IJ2(N)
        J3=IJ3(N)
        COEF(IEL,M)=COEF(IEL,M)+E(I1,J1)*E(I2,J2)*E(I3,J3)*EVECT(IND2)
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
      W=(COTE-0.5*(ZZZ(KS)+ZZZ(KS+1)))/(ZZZ(KS+1)-ZZZ(KS))
      FLZ(1)=1.0
      FLZ(2)=FLZ(1)*W
      FLZ(3)=FLZ(2)*W
      FLZ(4)=FLZ(3)*W
      FLZ(5)=FLZ(4)*W
      DO 90 L=1,LL
      I1=IJ1(L)
      I2=IJ2(L)
      I3=IJ3(L)
      AXYZ(I,J,K)=AXYZ(I,J,K)+COEF(IEL,L)*FLX(I1)*FLY(I2)*FLZ(I3)
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(COEF,IWRK)
      RETURN
      END
