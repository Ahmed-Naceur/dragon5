*DECK SNFP1P
      SUBROUTINE SNFP1P(NREG,NMAT,ISCHM,IELEM,IBFP,NLF,NSCT,U,W,
     1 PL,MAT,VOL,TOTAL,ESTOPW,NCODE,ZCODE,DELTAE,QEXT,QEXT0,LFIXUP,
     2 FLUX,FLUX0,ALP,IG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 1D slab
* geometry. Albedo boundary conditions. Boltzmann-Fokker-Planck (BFP)
* discretization.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
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
* ISCHM   method of spatial discretisation:
*         =1 High-Order Diamond Differencing (HODD) - default;
*         =2 Discontinuous Galerkin finite element method (DG).
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - only for HODD, classical diamond scheme 
*         (default for HODD);
*         =2 linear (default for DG);
*         =3 parabolic;
*         =4 cubic - only for DG.
* IBFP    type of energy proparation relation:
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
* NLF     number of $\\mu$ levels.
* NSCT    number of Legendre components in the flux:
*         =1: isotropic sources;
*         =2: linearly anisotropic sources.
* U       base points in $\\mu$ of the SN quadrature.
* W       weights of the SN quadrature.
* PL      discrete values of the Legendre polynomials corresponding
*         to the SN quadrature.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* ESTOPW  stopping power.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* DELTAE  energy group width in MeV.
* QEXT    Legendre components of the fixed source.
* QEXT0   initial slowing-down angular fluxes.
* LFIXUP  flag to enable negative flux fixup.
*
*Parameters: output
* FLUX    Legendre components of the flux.
* FLUX0   final slowing-down angular fluxes.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,NMAT,ISCHM,IELEM,IBFP,NLF,NSCT,MAT(NREG),NCODE(2),IG
      REAL U(NLF),W(NLF),PL(NSCT,NLF),VOL(NREG),TOTAL(0:NMAT),
     1 ESTOPW(0:NMAT,2),ZCODE(2),DELTAE,QEXT(IELEM,NSCT,NREG),
     2 QEXT0(IELEM,NLF,NREG),FLUX(IELEM,NSCT,NREG),
     3 FLUX0(IELEM,NLF,NREG),ALP
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION AFB,FEP(3),CURR(2),FN,FD,G1,G2,
     1 CNP(NLF),CNM(NLF),AN(NLF+1),L0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Q
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLUXN,Q2
      PARAMETER(RLOG=1.0E-8)
      INTEGER REGI, REGF, REGSTEP,NLFI,NLFF
*----
*  ALLOCATABLE ARRAYS
*----
      ALLOCATE(Q(IELEM),Q2(IELEM,IELEM+1),FLUXN(NLF,NREG))
*----
*  PARAMETER VALIDATION.
*----
      IF(ISCHM.EQ.1) THEN 
         IF((IELEM.LT.1).OR.(IELEM.GT.3))
     1      CALL XABORT('SNFP1P: INVALID IELEM (DIAM) VALUE. ' 
     2      //'CHECK INPUT DATA FILE.')
      ELSEIF(ISCHM.EQ.2) THEN
         IF((IELEM.LT.1).OR.(IELEM.GT.4)) 
     1      CALL XABORT('SNFP1P: INVALID IELEM (DIAM) VALUE. '
     2      //'CHECK INPUT DATA FILE.(2)')
      ENDIF
*----
*  COMPUTE A NORMALIZATION CONSTANT.
*----
      DENOM=0.0
      DO 10 I=1+NLF/2,NLF
      DENOM=DENOM+W(I)*U(I)
   10 CONTINUE
*----
*  Compute the angular FP coefficients
*----
      AN(1)=0.0
      DO 11 I=2,NLF+1
      AN(I)=AN(I-1)-2*U(I-1)*W(I-1)
  11  CONTINUE

      FLUXN(:NLF,:NREG)=0.0
      L0=0.0

      DO 12 IP=1,NLF
      DO 13 I=1,NREG
      DO 14 IL=0,NSCT-1
      FLUXN(IP,I)=FLUXN(IP,I)+0.5*(2.0*IL+1)*FLUX(1,IL+1,I)*PL(IL+1,IP)
   14 CONTINUE
   13 CONTINUE

         IF(IP.EQ.1) THEN
         CNM(IP)=0.0
         CNP(IP)=AN(IP+1)/(W(IP)*(U(IP+1)-U(IP)))
         L0=CNM(IP)+CNP(IP)
         ELSE IF(IP.EQ.NLF) THEN
         CNM(IP)=AN(IP)/(W(IP)*(U(IP)-U(IP-1)))
         CNP(IP)=0.0
         IF(CNM(IP)+CNP(IP).GT.L0) L0=CNM(IP)+CNP(IP)
         ELSE
         CNP(IP)=AN(IP+1)/(W(IP)*(U(IP+1)-U(IP)))
         CNM(IP)=AN(IP)/(W(IP)*(U(IP)-U(IP-1)))
         IF(CNM(IP)+CNP(IP).GT.L0) L0=CNM(IP)+CNP(IP)
         ENDIF

   12 CONTINUE

*----
*  OUTER LOOP OVER MU LEVELS.
*----

      CALL XDRSET(FLUX,IELEM*NSCT*NREG,0.0)
      NSCT0=NSCT
      CURR(1)=0.0D0
      CURR(2)=0.0D0
*----
*  HIGH-ORDER DIAMOND DIFFERENCING METHOD
*----
      DO 400 M=1,2
      IF(M.EQ.2) THEN
      REGI=1
      REGF=NREG
      REGSTEP=1
      NLFI=NLF/2+1
      NLFF=NLF
      ELSE
      REGI=NREG
      REGF=1
      REGSTEP=-1
      NLFI=1
      NLFF=NLF/2
      ENDIF
      DO 200 IP=NLFI,NLFF

* SOURCE MONODIRECTIONNELLE  
      IF(M.EQ.2.AND.IP.EQ.NLF.AND.IG.EQ.1) THEN
      N_GENERATIONS=2
      AFB=1.5915494D9/(N_GENERATIONS*W(IP))
      ELSE
      AFB=0.0
      ENDIF

      DO 120 I=REGI,REGF,REGSTEP

         IBM=MAT(I)

         DO 105 IEL=1,IELEM
         Q(IEL)=0.0
         DO 100 IL=0,NSCT-1
         Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,IP)/2.0
  100    CONTINUE
  105    CONTINUE

         IF(IBFP.EQ.1) THEN
*          Galerkin type
           FN=2.0D0*ESTOPW(IBM,1)/DELTAE
           FD=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
           G1=ESTOPW(IBM,1)/ESTOPW(IBM,2)
           G2=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/ESTOPW(IBM,2)
         ELSE
*          Przybylski and Ligou type
           FN=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
           FD=2.0D0*ESTOPW(IBM,2)/DELTAE
           G1=1.0D0
           G2=2.0D0
         ENDIF

         FEP(:IELEM)=QEXT0(:IELEM,IP,I)

         IF(IELEM.EQ.1) THEN
           Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*ABS(U(IP))*AFB
     1      +0.5*ALP*(CNP(IP)*FLUXN(IP+1,I)+CNM(IP)*FLUXN(IP-1,I)
     2      +(L0-CNP(IP)-CNM(IP))*FLUXN(IP,I))*VOL(I)
           Q2(1,1)=(FD+TOTAL(IBM)+0.5*ALP*L0)*VOL(I)+2.0D0*ABS(U(IP))
         ELSE IF(IELEM.EQ.2) THEN
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
           Q2(1,2)=2.0D0*SQRT(3.0D0)*U(IP)
           Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
           Q2(2,1)=Q2(1,2)
           Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)-6.0D0*ABS(U(IP))
           Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1     U(IP)*AFB
         ELSE IF(IELEM.EQ.3) THEN
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*ABS(U(IP))
           Q2(1,2)=0.0
           Q2(1,3)=2.0D0*SQRT(5.0D0)*ABS(U(IP))
           Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*ABS(U(IP))*AFB
           Q2(2,1)=Q2(1,2)
           Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
           Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(IP)
           Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
           Q2(3,1)=Q2(1,3)
           Q2(3,2)=Q2(2,3)
           Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)+10.0D0*ABS(U(IP))
           Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)+2.0D0*SQRT(5.0D0)*
     1     ABS(U(IP))*AFB
         ENDIF

         CALL ALSBD(IELEM,1,Q2,IER,IELEM)
         IF(IER.NE.0) CALL XABORT('SNFP1P: SINGULAR MATRIX(5).')

         IF(IELEM.EQ.1) THEN
           IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
           AFB=2.0D0*Q2(1,2)-AFB
           IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
         ELSE IF(IELEM.EQ.2) THEN
           AFB=AFB-2.0D0*SQRT(3.0D0)*Q2(2,3)
         ELSE IF(IELEM.EQ.3) THEN
           AFB=2.0D0*Q2(1,4)+2.0D0*SQRT(5.0D0)*Q2(3,4)-AFB
         ENDIF

         FEP(:IELEM)=G2*Q2(:IELEM,IELEM+1)-G1*FEP(:IELEM)
         FLUX0(:IELEM,IP,I)=REAL(FEP(:IELEM))/DELTAE
         DO 115 K=1,NSCT
         DO 110 IEL=1,IELEM
         FLUX(IEL,K,I)=FLUX(IEL,K,I)+W(IP)*REAL(Q2(IEL,IELEM+1))*
     1    PL(K,IP)
  110    CONTINUE
  115    CONTINUE
  120 CONTINUE
  200 CONTINUE
  400 CONTINUE

*----
*  END OF OUTER LOOP OVER MU LEVELS
*---- 
      DEALLOCATE(Q,Q2,FLUXN)
      RETURN
      END
