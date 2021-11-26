*DECK VALU1B
      SUBROUTINE VALU1B (IDIM,LX,LY,L4,X,Y,XXX,YYY,EVT,ISS,IELEM,IXLG,
     + IYLG,AXY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution for MCFD method in 2D.
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
* IDIM    number of dimensions (1 or 2).
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* L4      dimension of unknown array EVT.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* Y       Cartesian coordinates along the Y axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* EVT     variational coefficients of the flux.
* ISS     mixture index assigned to each element.
* IELEM   MCFD polynomial order (IELEM=1 is the mesh centered finite
*         difference method).
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.                                    
*                                                                      
*Parameters: output
* AXY     interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IDIM,LX,LY,L4,ISS(LX*LY),IELEM,IXLG,IYLG
      REAL X(IXLG),Y(IYLG),XXX(LX+1),YYY(LY+1),EVT(L4),AXY(IXLG,IYLG)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
*----
*  Scratch storage allocation
*----
      ALLOCATE(IWRK(LX*LY))
*
      NUM=0
      DO 10 K=1,LX*LY
      IF (ISS(K).EQ.0) GO TO 10
      NUM=NUM+1
      IWRK(K)=NUM
  10  CONTINUE
*
      LL4=L4/IELEM**(IDIM-1)
      DO 120 J=1,IYLG
      ORDO=Y(J)
      DO 110 I=1,IXLG
      ABSC=X(I)
      GAR=0.0
*                                                          
*     Find the finite element index containing the interpolation point
      IS=0
      JS=0
      DO 20 L=1,LX
      IS=L
      IF((ABSC.GE.XXX(L)).AND.(ABSC.LE.XXX(L+1))) GO TO 30
   20 CONTINUE
      CALL XABORT('VALU1B: WRONG INTERPOLATION(1).')
   30 DO 40 L=1,LY
      JS=L
      IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 70
   40 CONTINUE
      CALL XABORT('VALU1B: WRONG INTERPOLATION(2).')
   70 IEL=(JS-1)*LX+IS
      IF(ISS(IEL).EQ.0) GO TO 100
      U=(ABSC-0.5*(XXX(IS)+XXX(IS+1)))/(XXX(IS+1)-XXX(IS))
      V=(ORDO-0.5*(YYY(JS)+YYY(JS+1)))/(YYY(JS+1)-YYY(JS))
      L=1+IELEM*(IWRK(IEL)-1)
      DO 90 N2=0,IELEM-1
      DO 80 N1=0,IELEM-1
      GAR=GAR+VALPL(N1,U)*VALPL(N2,V)*EVT(LL4*N2+N1+L)
   80 CONTINUE
      IF ((IDIM.EQ.1).AND.(N2.EQ.0)) GO TO 100
   90 CONTINUE
  100 AXY(I,J)=GAR
  110 CONTINUE
  120 CONTINUE
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(IWRK)
      RETURN
      END
