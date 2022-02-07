*DECK VALUE4
      SUBROUTINE VALUE4(IELEM,NUN,LX,LY,LZ,X,Y,Z,XXX,YYY,ZZZ,EVECT,ISS,
     + KFLX,IXLG,IYLG,IZLG,AXYZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution for DUAL method in 3D.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input
* IELEM   finite element order
*         =1 : linear Raviart-Thomas
*         =2 : parabolic Raviart-Thomas
*         =3 : cubic Raviart-Thomas
*         =4 : quartic Raviart-Thomas
* NUN     number of unknowns
* LX      number of elements along the X axis.
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
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
* KFLX    correspondence between local and global numbering.
* IXLG    number of interpolated points according to X.
* IYLG    number of interpolated points according to Y.
* IZLG    number of interpolated points according to Z.
*
*Parameters: output
* AXYZ    interpolated fluxes.
*
*----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IELEM,NUN,LX,LY,LZ,IXLG,IYLG,IZLG,ISS(LX*LY*LZ),
     1   KFLX(LX*LY*LZ)
      REAL X(IXLG),Y(IYLG),Z(IZLG),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),
     1   EVECT(NUN),AXYZ(IXLG,IYLG,IZLG)
*----
*  LOCAL VARIABLES
*----
      INTEGER I,J,K,L,IS,JS,KS,IEL,I1,I2,I3,IE
      REAL COTE,ORDO,ABSC,COEF(2,5),FLX(5),FLY(5),FLZ(5)
      REAL U,V,W
*----
*  compute coefficient for legendre polynomials
*----
      CALL XDRSET(COEF,10,0.0)
      COEF(1,1)=1.0
      COEF(1,2)=2.*3.**0.5
      DO IE=1,3
        COEF(1,IE+2)=2.0*REAL(2*IE+1)/REAL(IE+1)
     1                     *(REAL(2*IE+3)/REAL(2*IE+1))**0.5
        COEF(2,IE+2)=REAL(IE)/REAL(IE+1)
     1                     *(REAL(2*IE+3)/REAL(2*IE-1))**0.5
      ENDDO
*----
*  perform interpolation
*----
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
      CALL XABORT('VALUE4: WRONG INTERPOLATION(1).')
   30 DO 40 L=1,LY
      JS=L
      IF((ORDO.GE.YYY(L)).AND.(ORDO.LE.YYY(L+1))) GO TO 50
   40 CONTINUE
      CALL XABORT('VALUE4: WRONG INTERPOLATION(2).')
   50 DO 60 L=1,LZ
      KS=L
      IF((COTE.GE.ZZZ(L)).AND.(COTE.LE.ZZZ(L+1))) GO TO 70
   60 CONTINUE
      CALL XABORT('VALUE4: WRONG INTERPOLATION(3).')
   70 IEL=(KS-1)*LX*LY+(JS-1)*LX+IS
C
      IF(ISS(IEL).EQ.0) GO TO 100
      U=(ABSC-0.5*(XXX(IS)+XXX(IS+1)))/(XXX(IS+1)-XXX(IS))
      FLX(1)=COEF(1,1)
      FLX(2)=COEF(1,2)*U
      V=(ORDO-0.5*(YYY(JS)+YYY(JS+1)))/(YYY(JS+1)-YYY(JS))
      FLY(1)=COEF(1,1)
      FLY(2)=COEF(1,2)*V
      W=(COTE-0.5*(ZZZ(KS)+ZZZ(KS+1)))/(ZZZ(KS+1)-ZZZ(KS))
      FLZ(1)=COEF(1,1)
      FLZ(2)=COEF(1,2)*W
      IF(IELEM.GE.2) THEN
        DO IE=2,IELEM
          FLX(IE+1)=FLX(IE)*U*COEF(1,IE+1)-FLX(IE-1)*COEF(2,IE+1)
          FLY(IE+1)=FLY(IE)*V*COEF(1,IE+1)-FLY(IE-1)*COEF(2,IE+1)
          FLZ(IE+1)=FLZ(IE)*W*COEF(1,IE+1)-FLZ(IE-1)*COEF(2,IE+1)
        ENDDO
      ENDIF
      DO 93 I3=1,IELEM
      DO 92 I2=1,IELEM
      DO 91 I1=1,IELEM
        L=(I3-1)*(IELEM)**2+(I2-1)*(IELEM)+I1
        AXYZ(I,J,K)=AXYZ(I,J,K)+EVECT(KFLX(IEL)+L-1)*FLX(I1)*FLY(I2)
     1                        *FLZ(I3)
   91 CONTINUE
   92 CONTINUE
   93 CONTINUE
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
      RETURN
      END
