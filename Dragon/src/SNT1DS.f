*DECK SNT1DS
      SUBROUTINE SNT1DS (IMPX,LX,NCODE,ZCODE,XXX,NLF,NSCT,U,W,ALPHA,
     1 PLZ,PL,VOL,IDL,SURF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a 1-D spherical geometry with discrete
* ordinates approximation of the flux.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IMPX    print parameter.
* LX      number of elements along the R axis.
* NCODE   type of boundary condition applied on each side
*         (i=1 R-;  i=2 R+):
*         =1 VOID;   =2 REFL;   =7 ZERO.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* XXX     coordinates along the R axis.
* NLF     order of the SN approximation (even number).
* NSCT    maximum number of spherical harmonics moments of the flux.
*
*Parameters: output
* U       base points in $\\xi$ of the axial quadrature.
* W       weights for the quadrature in $\\mu$.
* ALPHA   angular redistribution terms.
* PLZ     discrete values of the spherical harmonics corresponding
*         to the 1D quadrature.Used with zero-weight points.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* VOL     volume of each element.
* IDL     position of integrated fluxes into unknown vector.
* SURF    surfaces.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,LX,NCODE(2),NLF,NSCT,IDL(LX)
      REAL ZCODE(2),XXX(LX+1),U(NLF),W(NLF),PLZ(NSCT),PL(NSCT,NLF),
     1 ALPHA(NLF),VOL(LX),SURF(LX+1)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(PI=3.141592654)
*----
*  GENERATE QUADRATURE BASE POINTS AND CORRESPONDING WEIGHTS.
*----
      IF(MOD(NLF,2).EQ.1) CALL XABORT('SNT1DS: EVEN NLF EXPECTED.')
      IF(NLF.EQ.2) THEN
      	U(1)=-1.0/SQRT(3.0)
      	W(1)=1.0
      	U(2)=1.0/SQRT(3.0)
      	W(2)=1.0
      ELSE
        CALL ALGPT(NLF,-1.0,1.0,U,W)
      ENDIF
*----
*  COMPUTE ALPHA.
*----
      SUMETA=0.0
      DO IP=1,NLF
         SUMETA=SUMETA-2.0*W(IP)*U(IP)
         ALPHA(IP)=SUMETA
      ENDDO
*----
*  PRINT COSINES AND WEIGHTS.
*----
      IF(IMPX.GT.1) THEN
         WRITE(6,60) (U(M),M=1,NLF)
   60    FORMAT(//,1X,'THE POSITIVE QUADRATURE COSINES FOLLOW'//
     1   (1X,5E14.6))
         WRITE(6,70) (W(M),M=1,NLF)
   70    FORMAT(//,1X,'THE CORRESPONDING QUADRATURE WEIGHTS FOLLOW'//
     1   (1X,5E14.6))
         WRITE(6,76) (ALPHA(M),M=1,NLF)
   76    FORMAT(//,1X,'THE ALPHAS FOLLOW'//(1X,5E14.6))
      ENDIF
*----
*  GENERATE LEGENDRE POLYNOMIALS FOR SCATTERING SOURCE.
*----
      DO 150 M=1,NLF
      PL(1,M)=1.0
      IF(NSCT.GT.1) THEN
         PL(2,M)=U(M)
         DO 140 L=2,NSCT-1
         PL(L+1,M)=((2.0*L-1.0)*U(M)*PL(L,M)-(L-1)*PL(L-1,M))/REAL(L)
  140    CONTINUE
      ENDIF
  150 CONTINUE
      PLZ(1)=1.0
      IF(NSCT.GT.1) THEN
         PLZ(2)=-1.0
         DO 160 L=2,NSCT-1
         PLZ(L+1)=(-(2.0*L-1.0)*PLZ(L)-(L-1)*PLZ(L-1))/REAL(L)
  160    CONTINUE
      ENDIF
*----
*  COMPUTE SURFACES AND VOLUMES.
*----
      DO 200 I=1,LX+1
      SURF(I)=4.0*PI*XXX(I)*XXX(I)
  200 CONTINUE
      DO 210 I=1,LX
      IDL(I)=(I-1)*NSCT+1
      VOL(I)=4.0*PI*(XXX(I+1)**3-XXX(I)**3)/3.0
  210 CONTINUE
*----
*  SET BOUNDARY CONDITIONS.
*----
      IF(NCODE(2).EQ.2) ZCODE(2)=1.0
*
      RETURN
      END
