*DECK SNT1DC
      SUBROUTINE SNT1DC (IMPX,LX,NCODE,ZCODE,XXX,NLF,NPQ,NSCT,IQUAD,
     1 JOP,U,W,TPQ,UPQ,VPQ,WPQ,ALPHA,PLZ,PL,VOL,IDL,SURF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a 1-D cylindrical geometry with discrete
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
*         (i=1 R-; i=2 R+):
*         =1 VOID;   =2 REFL;  =7 ZERO.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* XXX     coordinates along the R axis.
* NLF     order of the SN approximation (even number).
* NPQ     number of SN directions in two octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* IQUAD   type of SN quadrature (1 Level symmetric, type IQUAD;
*         4 Gauss-Legendre and Gauss-Chebyshev; 10 product).
*
*Parameters: output
* JOP     number of base points per axial level in one octant.
* U       base points in $\\xi$ of the axial quadrature. Used with
*         zero-weight points.
* W       weights for the quadrature in $\\mu$.
* TPQ     base points in $\\xi$ of the 2D SN quadrature.
* UPQ     base points in $\\mu$ of the 2D SN quadrature.
* VPQ     base points in $\\eta$ of the 2D SN quadrature.
* WPQ     weights of the 2D SN quadrature.
* ALPHA   angular redistribution terms.
* PLZ     discrete values of the spherical harmonics corresponding
*         to the 1D quadrature. Used with zero-weight points.
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
      INTEGER IMPX,LX,NCODE(2),NLF,NPQ,NSCT,IQUAD,JOP(NLF/2),IDL(LX)
      REAL ZCODE(2),XXX(LX+1),U(NLF/2),W(NLF/2),PLZ(NSCT,NLF/2),
     1 PL(NSCT,NPQ),TPQ(NPQ),UPQ(NPQ),VPQ(NPQ),WPQ(NPQ),ALPHA(NPQ),
     2 VOL(LX),SURF(LX+1)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(PI=3.141592654)
      REAL, ALLOCATABLE, DIMENSION(:) :: TPQ2,UPQ2,VPQ2,WPQ2
*----
*  GENERATE QUADRATURE BASE POINTS AND CORRESPONDING WEIGHTS.
*----
      IF(MOD(NLF,2).EQ.1) CALL XABORT('SNT1DC: EVEN NLF EXPECTED.')
      M2=NLF/2
      ALLOCATE(TPQ2(NPQ/2),UPQ2(NPQ/2),VPQ2(NPQ/2),WPQ2(NPQ/2))
      IF(IQUAD.EQ.1) THEN
         CALL SNQU01(NLF,JOP,U,W,TPQ2,UPQ2,VPQ2,WPQ2)
      ELSE IF(IQUAD.EQ.2) THEN
         CALL SNQU02(NLF,JOP,U,W,TPQ2,UPQ2,VPQ2,WPQ2)
      ELSE IF(IQUAD.EQ.3) THEN
         CALL SNQU03(NLF,JOP,U,W,TPQ2,UPQ2,VPQ2,WPQ2)
      ELSE IF(IQUAD.EQ.4) THEN
         CALL SNQU04(NLF,JOP,U,W,TPQ2,UPQ2,VPQ2,WPQ2)
      ELSE IF(IQUAD.EQ.10) THEN
         CALL SNQU10(NLF,JOP,U,W,TPQ2,UPQ2,VPQ2,WPQ2)
      ELSE
         CALL XABORT('SNT1DC: UNKNOWN QUADRATURE TYPE.')
      ENDIF
      IPQ=0
      IPR=0
      WSUM=0.0
      DO IP=1,M2
         DO IQ=1,JOP(IP)
            IPQ=IPQ+1
            IPR=IPR+1
            TPQ(IPQ)=TPQ2(IPR)
            UPQ(IPQ)=UPQ2(IPR)
            VPQ(IPQ)=VPQ2(IPR)
            WPQ(IPQ)=WPQ2(IPR)
            WSUM=WSUM+WPQ(IPQ)
         ENDDO
         DO IQ=JOP(IP)+1,2*JOP(IP)
            IPQ=IPQ+1
            JQ=IQ-(JOP(IP)+1)+1
            JPQ=IPQ-2*JQ+1
            TPQ(IPQ)=TPQ(JPQ)
            UPQ(IPQ)=-UPQ(JPQ)
            VPQ(IPQ)=VPQ(JPQ)
            WPQ(IPQ)=WPQ(JPQ)
            WSUM=WSUM+WPQ(IPQ)
         ENDDO
      ENDDO
      DEALLOCATE(WPQ2,VPQ2,UPQ2,TPQ2)
      IF(IPQ.NE.NPQ) CALL XABORT('SNT1DC: BAD VALUE ON NPQ.')
*----
*  COMPUTE ALPHA.
*----
      IPQ=0
      DO IP=1,M2
         SUMETA=0.0
         DO IQ=1,2*JOP(IP)
            IPQ=IPQ+1
            ALPHA(IPQ)=SUMETA
            SUMETA=SUMETA+WPQ(IPQ)*UPQ(IPQ)
         ENDDO
      ENDDO
*----
*  PRINT COSINES AND WEIGHTS.
*----
      IF(IMPX.GT.1) THEN
         WRITE(6,'(/20H SNT1DC: WEIGHT SUM=,1P,E11.4)') WSUM
         WRITE(6,60) (U(M),M=1,M2)
   60    FORMAT(//,1X,'THE POSITIVE QUADRATURE COSINES FOLLOW'//
     1   (1X,5E14.6))
         WRITE(6,70) (W(M),M=1,M2)
   70    FORMAT(//,1X,'THE CORRESPONDING QUADRATURE WEIGHTS FOLLOW'//
     1   (1X,5E14.6))
         WRITE(6,74) (UPQ(M),M=1,NPQ)
   74    FORMAT(//,1X,'THE BASE POINTS (MU) FOLLOW'//(1X,5E14.6))
         WRITE(6,75) (WPQ(M),M=1,NPQ)
   75    FORMAT(//,1X,'THE WEIGHTS FOLLOW'//(1X,5E14.6))
         WRITE(6,76) (ALPHA(M),M=1,NPQ)
   76    FORMAT(//,1X,'THE ALPHAS FOLLOW'//(1X,5E14.6))
      ENDIF
*----
*  GENERATE SPHERICAL HARMONICS FOR SCATTERING SOURCE.
*----
      IOF=0
      DO 130 IL=0,NLF-1
      DO 120 IM=0,IL
      IF(MOD(IL+IM,2).EQ.1) GO TO 120
      IOF=IOF+1
      IF(IOF.GT.NSCT) GO TO 140
      DO 100 M=1,M2
      PLZ(IOF,M)=PNSH(IL,IM,U(M),-SQRT(1.0-U(M)*U(M)),0.0)
  100 CONTINUE
      DO 110 M=1,NPQ
      PL(IOF,M)=PNSH(IL,IM,TPQ(M),UPQ(M),VPQ(M))
  110 CONTINUE
  120 CONTINUE
  130 CONTINUE
*----
*  COMPUTE SURFACES AND VOLUMES.
*----
  140 DO 150 I=1,LX+1
      SURF(I)=2.0*PI*XXX(I)
  150 CONTINUE
      DO 160 I=1,LX
      IDL(I)=(I-1)*NSCT+1
      VOL(I)=PI*(XXX(I+1)*XXX(I+1)-XXX(I)*XXX(I))
  160 CONTINUE
*----
*  SET BOUNDARY CONDITIONS.
*----
      IF(NCODE(2).EQ.2) ZCODE(2)=1.0
*
      RETURN
      END
