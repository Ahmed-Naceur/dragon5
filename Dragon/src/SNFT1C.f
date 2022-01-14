*DECK SNFT1C
      SUBROUTINE SNFT1C(NREG,NMAT,M2,NPQ,ISCAT,NSCT,JOP,U,UPQ,WPQ,ALPHA,
     1 PLZ,PL,MAT,VOL,SURF,TOTAL,IGAV,QEXT,LFIXUP,CURR,FLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 1D cylindrical
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
* M2      number of axial $\\xi$ levels.
* NPQ     number of SN directions in one octant.
* ISCAT   anisotropy of one-speed sources:
*         =1 isotropic sources;
*         =2 linearly anisotropic sources.
* NSCT    number of spherical harmonics components in the flux.
* JOP     number of base points per axial level in one octant.
* U       base points in $\\mu$ of the 1D quadrature. Used with
*         zero-weight points.
* UPQ     base points in $\\mu$ of the 2D SN quadrature.
* WPQ     weights of the 2D SN quadrature.
* ALPHA   angular redistribution terms.
* PLZ     discrete values of the spherical harmonics corresponding
*         to the 1D quadrature. Used with zero-weight points.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* SURF    surfaces surrounding each region.
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
      INTEGER NREG,NMAT,M2,NPQ,ISCAT,NSCT,JOP(M2),MAT(NREG),IGAV
      REAL U(M2),UPQ(NPQ),WPQ(NPQ),ALPHA(NPQ),PLZ(NSCT,M2),PL(NSCT,NPQ),
     1 VOL(NREG),SURF(NREG+1),TOTAL(0:NMAT),QEXT(NSCT,NREG),CURR,
     2 FLUX(NSCT,NREG)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION Q,E2,AFB,Q1,Q2,PSIA,WSIA,CURSUM
      REAL, ALLOCATABLE, DIMENSION(:) :: FLXB,AFGL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FLXB(M2),AFGL(NREG))
*----
*  COMPUTE A NORMALIZATION CONSTANT.
*----
      DENOM=0.0
      DO 10 I=1,NPQ
      IF(UPQ(I).GT.0.0) DENOM=DENOM+4.0*WPQ(I)*UPQ(I)
   10 CONTINUE
*----
*  OUTER LOOP OVER AXIAL LEVELS.
*----
      CALL XDRSET(FLUX,NSCT*NREG,0.0)
      CURSUM=0.0D0
      IPQ=0
      DO 200 IP=1,M2
*----
*  INITIALIZATION SWEEP (USING ZERO-WEIGHT POINTS). AFB IS THE EDGE
*  FLUX VALUE AND AFGL(I) IS THE CENTERED FLUX VALUE.
*----
      AFB=CURR/DENOM
      DO 30 I=NREG,1,-1 ! Spatial sweep
         Q=0.0D0
         IBM=MAT(I)
         IOF=0
         DO 25 IL=0,ISCAT-1
         DO 20 IM=0,IL
         IF(MOD(IL+IM,2).EQ.1) GO TO 20
         IOF=IOF+1
         Q=Q+QEXT(IOF,I)*PLZ(IOF,IP)/(4.0*PI)
   20    CONTINUE
   25    CONTINUE
         Q1=-4.0D0*PI*SQRT(1.0-U(IP)*U(IP))/(SURF(I+1)-SURF(I))
         E2=(Q-Q1*AFB)/(TOTAL(IBM)-Q1)
         IF(LFIXUP.AND.(E2.LE.RLOG)) E2=0.0
         AFB=2.0*E2-AFB
         IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
         AFGL(I)=REAL(E2)
         IF(LFIXUP.AND.(AFGL(I).LE.RLOG)) AFGL(I)=0.0
   30 CONTINUE
      WSIA=0.0D0
      IF(IGAV.EQ.2) THEN
         PSIA=0.0D0
      ELSE
         PSIA=AFB
      ENDIF
*----
*  BACKWARD SWEEP (FROM EXTERNAL SURFACE TOWARD CENTRAL AXIS).
*----
      ALPMAX=0.0
      DO 80 M=2*JOP(IP),JOP(IP)+1,-1 ! Angular sweep
         AFB=CURR/DENOM
         ALPMIN=ALPHA(IPQ+M)
         DO 70 I=NREG,1,-1 ! Spatial sweep
            Q=0.0
            IBM=MAT(I)
            IOF=0
            DO 45 IL=0,ISCAT-1
            DO 40 IM=0,IL
            IF(MOD(IL+IM,2).EQ.1) GO TO 40
            IOF=IOF+1
            Q=Q+QEXT(IOF,I)*PL(IOF,IPQ+M)/(4.0*PI)
   40       CONTINUE
   45       CONTINUE
            Q1=Q*VOL(I)-UPQ(IPQ+M)*(SURF(I)+SURF(I+1))*AFB+(SURF(I+1)-
     1      SURF(I))*(ALPMIN+ALPMAX)*AFGL(I)/WPQ(IPQ+M)
            Q2=TOTAL(IBM)*VOL(I)-2.0D0*UPQ(IPQ+M)*SURF(I)+2.0D0*
     1      (SURF(I+1)-SURF(I))*ALPMIN/WPQ(IPQ+M)
            E2=Q1/Q2
            IF(LFIXUP.AND.(E2.LE.RLOG)) E2=0.0
            AFB=2.0*E2-AFB ! IN SPACE
            IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
            AFGL(I)=2.0*REAL(E2)-AFGL(I) ! IN ANGLE
            IF(LFIXUP.AND.(AFGL(I).LE.RLOG)) AFGL(I)=0.0
            DO 60 K=1,NSCT
            FLUX(K,I)=FLUX(K,I)+4.0*WPQ(IPQ+M)*REAL(E2)*PL(K,IPQ+M)
   60       CONTINUE
   70    CONTINUE
         IF(IGAV.EQ.1) THEN
            FLXB(2*JOP(IP)-M+1)=REAL(AFB)
         ELSE IF(IGAV.EQ.2) THEN
            PSIA=PSIA+WPQ(IPQ+M)*REAL(AFB)
            WSIA=WSIA+WPQ(IPQ+M)
         ENDIF
         ALPMAX=ALPMIN
   80 CONTINUE
      IF(IGAV.EQ.2) PSIA=PSIA/WSIA
*----
*  FORWARD SWEEP (FROM CENTRAL AXIS TOWARD EXTERNAL SURFACE).
*----
      DO 130 M=JOP(IP),1,-1 ! Angular sweep
         ALPMIN=ALPHA(IPQ+M)
         IF(IGAV.EQ.1) THEN
            AFB=FLXB(M)
         ELSE IF(IGAV.EQ.2) THEN
            AFB=PSIA
         ELSE IF(IGAV.EQ.3) THEN
            AFB=PSIA
         ENDIF
         DO 120 I=1,NREG ! Spatial sweep
            Q=0.0
            IBM=MAT(I)
            IOF=0
            DO 100 IL=0,ISCAT-1
            DO 90 IM=0,IL
            IF(MOD(IL+IM,2).EQ.1) GO TO 90
            IOF=IOF+1
            Q=Q+QEXT(IOF,I)*PL(IOF,IPQ+M)/(4.0*PI)
   90       CONTINUE
  100       CONTINUE
            Q1=Q*VOL(I)+UPQ(IPQ+M)*(SURF(I)+SURF(I+1))*AFB+(SURF(I+1)-
     1      SURF(I))*(ALPMIN+ALPMAX)*AFGL(I)/WPQ(IPQ+M)
            Q2=TOTAL(IBM)*VOL(I)+2.0D0*UPQ(IPQ+M)*SURF(I+1)+2.0D0*
     1      (SURF(I+1)-SURF(I))*ALPMIN/WPQ(IPQ+M)
            E2=Q1/Q2
            IF(LFIXUP.AND.(E2.LE.RLOG)) E2=0.0
            AFB=2.0*E2-AFB ! IN SPACE
            IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
            AFGL(I)=2.0*REAL(E2)-AFGL(I) ! IN ANGLE
            IF(LFIXUP.AND.(AFGL(I).LE.RLOG)) AFGL(I)=0.0
            DO 110 K=1,NSCT
            FLUX(K,I)=FLUX(K,I)+4.0*WPQ(IPQ+M)*REAL(E2)*PL(K,IPQ+M)
  110       CONTINUE
  120    CONTINUE
         CURSUM=CURSUM+4.0*WPQ(IPQ+M)*UPQ(IPQ+M)*AFB
         ALPMAX=ALPMIN
  130 CONTINUE
*----
*  END OF OUTER LOOP OVER AXIAL LEVELS
*----
      IPQ=IPQ+2*JOP(IP)
  200 CONTINUE
      IF(IPQ.NE.NPQ) CALL XABORT('SN1T1C: BUG.')
      CURR=REAL(CURSUM)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(AFGL,FLXB)
      RETURN
      END
