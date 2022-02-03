*DECK SNFT1S
      SUBROUTINE SNFT1S(NREG,NMAT,NLF,NSCT,U,W,ALPHA,PLZ,PL,MAT,VOL,
     1 SURF,XXX,TOTAL,IGAV,QEXT,LFIXUP,CURR,FLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 1D spherical
* geometry.
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
* NREG    number of regions.
* NMAT    number of material mixtures.
* NLF     number of SN directions.
* NSCT    number of Legendre components in the flux.
*         =1 isotropic sources;
*         =2 linearly anisotropic sources.
* U       base points in $\\mu$ of the 1D quadrature.
* W       weights of the quadrature.
* ALPHA   angular redistribution terms.
* PLZ     discrete values of the spherical harmonics corresponding
*         to the 1D quadrature. Used with zero-weight points.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* SURF    surfaces surrounding each region.
* XXX     radii.
* TOTAL   macroscopic total cross sections.
* IGAV    type of condition at axial axis (=1 specular reflection;
*         =2 zero-weight reflection; =3 averaged reflection).
* QEXT    spherical harmonics components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
*
*Parameters: input/output
* CURR    entering current at input and leaving current at output.
*
*Parameters: output
* FLUX    spherical harmonics components of the flux.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,NMAT,NLF,NSCT,MAT(NREG),IGAV
      REAL U(NLF),W(NLF),ALPHA(NLF),PLZ(NSCT),PL(NSCT,NLF),VOL(NREG),
     1 SURF(NREG+1),XXX(NREG+1),TOTAL(0:NMAT),QEXT(NSCT,NREG),CURR,
     2 FLUX(NSCT,NREG)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION Q,E2,AFB,Q1,Q2,PSIA,WSIA,CURSUM
      PARAMETER(RLOG=1.0E-8)
      REAL, ALLOCATABLE, DIMENSION(:) :: AFGL,FLXB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(AFGL(NREG),FLXB(NLF/2))
*----
*  COMPUTE A NORMALIZATION CONSTANT.
*----
      DENOM=0.0
      DO 10 I=1,NLF
      IF(U(I).GT.0.0) DENOM=DENOM+W(I)*U(I)
   10 CONTINUE
*----
*  INITIALIZATION SWEEP (USING ZERO-WEIGHT POINTS). AFB IS THE EDGE
*  FLUX VALUE AND AFGL(I) IS THE CENTERED FLUX VALUE.
*----
      AFB=CURR/DENOM
      DO 30 I=NREG,1,-1 ! Spatial sweep
         Q=0.0D0
         IBM=MAT(I)
         IOF=0
         DO 20 IL=0,NSCT-1
         Q=Q+QEXT(IL+1,I)*PLZ(IL+1)/2.0
   20    CONTINUE
         Q1=(XXX(I+1)-XXX(I))*Q+2.0*AFB
         Q2=(XXX(I+1)-XXX(I))*TOTAL(IBM)+2.0
         E2=Q1/Q2
         IF(LFIXUP.AND.(E2.LE.RLOG)) E2=0.0
         AFB=2.0*E2-AFB
         IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
         AFGL(I)=REAL(E2)
         IF(LFIXUP.AND.(AFGL(I).LE.RLOG)) AFGL(I)=0.0
   30 CONTINUE
      IF(IGAV.EQ.2) THEN
         PSIA=0.0D0
         WSIA=0.0D0
      ELSE
         PSIA=AFB
      ENDIF
*----
*  BACKWARD SWEEP (FROM EXTERNAL SURFACE TOWARD CENTRAL AXIS).
*----
      CALL XDRSET(FLUX,NSCT*NREG,0.0)
      CURSUM=0.0D0
      ALPMIN=0.0
      DO 80 IP=1,NLF/2
      AFB=CURR/DENOM
      ALPMAX=ALPHA(IP)
      DO 70 I=NREG,1,-1 ! Spatial sweep
         Q=0.0
         IBM=MAT(I)
         IOF=0
         DO 40 IL=0,NSCT-1
         Q=Q+QEXT(IL+1,I)*PL(IL+1,IP)/2.0
   40    CONTINUE
         Q1=Q*VOL(I)-U(IP)*(SURF(I)+SURF(I+1))*AFB+0.5D0*(SURF(I+1)-
     1   SURF(I))*(ALPMIN+ALPMAX)*AFGL(I)/W(IP)
         Q2=TOTAL(IBM)*VOL(I)-2.0D0*U(IP)*SURF(I)+(SURF(I+1)-SURF(I))*
     1   ALPMAX/W(IP)
         E2=Q1/Q2
         IF(LFIXUP.AND.(E2.LE.RLOG)) E2=0.0
         AFB=2.0*E2-AFB ! IN SPACE
         IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
         AFGL(I)=2.0*REAL(E2)-AFGL(I) ! IN ANGLE
         IF(LFIXUP.AND.(AFGL(I).LE.RLOG)) AFGL(I)=0.0
         DO 60 K=1,NSCT
         FLUX(K,I)=FLUX(K,I)+W(IP)*REAL(E2)*PL(K,IP)
   60    CONTINUE
   70 CONTINUE
      IF(IGAV.EQ.1) THEN
         FLXB(IP)=REAL(AFB)
      ELSE IF(IGAV.EQ.2) THEN
         PSIA=PSIA+W(IP)*REAL(AFB)
         WSIA=WSIA+W(IP)
      ENDIF
      ALPMIN=ALPMAX
   80 CONTINUE
      IF(IGAV.EQ.2) PSIA=PSIA/WSIA
*----
*  FORWARD SWEEP (FROM CENTRAL AXIS TOWARD EXTERNAL SURFACE).
*----
      DO 130 IP=1+NLF/2,NLF
      ALPMAX=ALPHA(IP)
      IF(IGAV.EQ.1) THEN
         AFB=FLXB(NLF-IP+1)
      ELSE IF(IGAV.EQ.2) THEN
         AFB=PSIA
      ELSE IF(IGAV.EQ.3) THEN
         AFB=PSIA
      ENDIF
      DO 120 I=1,NREG ! Spatial sweep
         Q=0.0
         IBM=MAT(I)
         IOF=0
         DO 90 IL=0,NSCT-1
         Q=Q+QEXT(IL+1,I)*PL(IL+1,IP)/2.0
   90    CONTINUE
         Q1=Q*VOL(I)+U(IP)*(SURF(I)+SURF(I+1))*AFB+0.5D0*(SURF(I+1)-
     1   SURF(I))*(ALPMIN+ALPMAX)*AFGL(I)/W(IP)
         Q2=TOTAL(IBM)*VOL(I)+2.0D0*U(IP)*SURF(I+1)+(SURF(I+1)-SURF(I))*
     1   ALPMAX/W(IP)
         E2=Q1/Q2
         IF(LFIXUP.AND.(E2.LE.RLOG)) E2=0.0
         AFB=2.0*E2-AFB ! IN SPACE
         IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
         AFGL(I)=2.0*REAL(E2)-AFGL(I) ! IN ANGLE
         IF(LFIXUP.AND.(AFGL(I).LE.RLOG)) AFGL(I)=0.0
         DO 110 K=1,NSCT
         FLUX(K,I)=FLUX(K,I)+W(IP)*REAL(E2)*PL(K,IP)
  110    CONTINUE
  120 CONTINUE
      CURSUM=CURSUM+W(IP)*U(IP)*AFB
      ALPMIN=ALPMAX
  130 CONTINUE
      CURR=REAL(CURSUM)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLXB,AFGL)
      RETURN
      END
