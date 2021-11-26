*DECK VALUE1
      SUBROUTINE VALUE1 (IDIM,LX,LY,LZ,L4,X,Y,Z,XXX,YYY,ZZZ,EVT,ISS,
     1 IELEM,IXLG,IYLG,IZLG,AXYZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution for MCFD method in 3D.
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
* IDIM    number of dimensions (1, 2 or 3).
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
* L4      dimension of unknown array EVT.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* Y       Cartesian coordinates along the Y axis where the flux is
*         interpolated.
* Z       Cartesian coordinates along the Z axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
* EVT     variational coefficients of the flux.
* ISS     mixture index assigned to each element.
* IELEM   MCFD polynomial order (IELEM=1 is the mesh centered finite
*         difference method).
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.                                    
* IZLG    number of interpolated points according to Z.                                    
*                                                                      
*Parameters: output
* AXYZ    interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IDIM,LX,LY,LZ,L4,ISS(LX*LY*LZ),IELEM,IXLG,IYLG,IZLG
      REAL X(IXLG),Y(IYLG),Z(IZLG),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),
     1 EVT(L4),AXYZ(IXLG,IYLG,IZLG)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
*----
*  Scratch storage allocation
*----
      ALLOCATE(IWRK(LX*LY*LZ))
*
      NUM=0
      DO 10 K=1,LX*LY*LZ
      IF (ISS(K).EQ.0) GO TO 10
      NUM=NUM+1
      IWRK(K)=NUM
  10  CONTINUE
*
      LL4=L4/IELEM**(IDIM-1)
      DO 130 K=1,IZLG
      COTE=Z(K)
      DO 120 J=1,IYLG
      ORDO=Y(J)
      DO 110 I=1,IXLG
      ABSC=X(I)
      GAR=0.0
*                                                          
*     Find the finite element index containing the interpolation point
      IS=0
      JS=0
      KS=0
      DO 20 L=1,LX
      IS=L
      IF((ABSC.GE.XXX(L)).AND.(ABSC.LE.XXX(L+1))) GO TO 30
   20 CONTINUE
      CALL XABORT('VALUE1: WRONG INTERPOLATION(1).')
   30 DO 40 L=1,LY
      JS=L
      IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 50
   40 CONTINUE
      CALL XABORT('VALUE1: WRONG INTERPOLATION(2).')
   50 DO 60 L=1,LZ
      KS=L
      IF((COTE.GE.ZZZ(L)).AND.(COTE.LE.ZZZ(L+1))) GO TO 70
   60 CONTINUE
      CALL XABORT('VALUE1: WRONG INTERPOLATION(3).')
   70 IEL=(KS-1)*LX*LY+(JS-1)*LX+IS
      IF(ISS(IEL).EQ.0) GO TO 100
      U=(ABSC-0.5*(XXX(IS)+XXX(IS+1)))/(XXX(IS+1)-XXX(IS))
      V=(ORDO-0.5*(YYY(JS)+YYY(JS+1)))/(YYY(JS+1)-YYY(JS))
      W=(COTE-0.5*(ZZZ(KS)+ZZZ(KS+1)))/(ZZZ(KS+1)-ZZZ(KS))
      L=1+IELEM*(IWRK(IEL)-1)
      DO 95 N3=0,IELEM-1
      DO 90 N2=0,IELEM-1
      DO 80 N1=0,IELEM-1
      GAR=GAR+VALPL(N1,U)*VALPL(N2,V)*VALPL(N3,W)*
     1    EVT(LL4*(IELEM*N3+N2)+N1+L)
   80 CONTINUE
      IF ((IDIM.EQ.1).AND.(N2.EQ.0)) GO TO 100
      IF ((IDIM.EQ.2).AND.(N2.EQ.IELEM-1)) GO TO 100
   90 CONTINUE
   95 CONTINUE
  100 AXYZ(I,J,K)=GAR
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(IWRK)
      RETURN
      END
