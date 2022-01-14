*DECK SNFP1P
      SUBROUTINE SNFP1P(NREG,NMAT,ISCHM,IELEM,IBFP,NLF,NSCT,U,W,PL,MAT,
     1 VOL,TOTAL,ESTOPW,NCODE,ZCODE,DELTAE,QEXT,QEXT0,LFIXUP,FLUX,FLUX0)
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
      INTEGER NREG,NMAT,ISCHM,IELEM,IBFP,NLF,NSCT,MAT(NREG),NCODE(2)
      REAL U(NLF),W(NLF),PL(NSCT,NLF),VOL(NREG),TOTAL(0:NMAT),
     1 ESTOPW(0:NMAT,2),ZCODE(2),DELTAE,QEXT(IELEM,NSCT,NREG),
     2 QEXT0(IELEM,NLF,NREG),FLUX(IELEM,NSCT,NREG),FLUX0(IELEM,NLF,NREG)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION AFB,AFB1,AFB2,FEP(3),CURR(2),FN,FD,G1,G2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Q
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Q2
      PARAMETER(RLOG=1.0E-8)
*----
*  ALLOCATABLE ARRAYS
*----
      ALLOCATE(Q(IELEM),Q2(IELEM,IELEM+1))
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
*  OUTER LOOP OVER MU LEVELS.
*----
      CALL XDRSET(FLUX,IELEM*NSCT*NREG,0.0)
      NSCT0=NSCT
      CURR(1)=0.0D0
      CURR(2)=0.0D0
      DO 200 IP=1,NLF/2
      AFB=0.0D0
      AFB1=0.0D0
      AFB2=0.0D0
*----
*  HIGH-ORDER DIAMOND DIFFERENCING METHOD
*----
      IF(ISCHM.EQ.1) THEN 
*----
*  SHOOTING METHOD FOR NON-VACUUM RIGHT BOUNDARY CONDITION.
*----
      IF(ZCODE(2).NE.0.0) THEN
*----
*  BACKWARD SWEEP (FROM RIGHT TO LEFT).
*----
         DO 30 I=NREG,1,-1
            IBM=MAT(I)
            DO 25 IEL=1,IELEM
            Q(IEL)=0.0
            DO 20 IL=0,NSCT-1
            Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,IP)/2.0
   20       CONTINUE
   25       CONTINUE
            IF(IBFP.EQ.1) THEN
*             Galerkin type
              FN=2.0D0*ESTOPW(IBM,1)/DELTAE
              FD=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              G1=ESTOPW(IBM,1)/ESTOPW(IBM,2)
              G2=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/ESTOPW(IBM,2)
            ELSE
*             Przybylski and Ligou type
              FN=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              FD=2.0D0*ESTOPW(IBM,2)/DELTAE
              G1=1.0D0
              G2=2.0D0
            ENDIF
            FEP(:IELEM)=QEXT0(:IELEM,NLF/2+IP,I)
            IF(IELEM.EQ.1) THEN
              Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)-2.0D0*U(IP)*AFB
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)-2.0D0*U(IP)
            ELSE IF(IELEM.EQ.2) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
              Q2(1,2)=2.0D0*SQRT(3.0D0)*U(IP)
              Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)+6.0D0*U(IP)
              Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1        U(IP)*AFB
            ELSE IF(IELEM.EQ.3) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)-2.0D0*U(IP)
              Q2(1,2)=0.0
              Q2(1,3)=-2.0D0*SQRT(5.0D0)*U(IP)
              Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)-2.0D0*U(IP)*AFB
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
              Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(IP)
              Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
              Q2(3,1)=Q2(1,3)
              Q2(3,2)=Q2(2,3)
              Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)-10.0D0*U(IP)
              Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)-2.0D0*SQRT(5.0D0)*
     1        U(IP)*AFB
            ENDIF
            CALL ALSBD(IELEM,1,Q2,IER,IELEM)
            IF(IER.NE.0) CALL XABORT('SNFP1P: SINGULAR MATRIX(1).')
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
   30    CONTINUE
*----
*  FORWARD SWEEP (FROM LEFT TO RIGHT).
*----
         AFBA1=0.0
         IF(NCODE(1).EQ.4) THEN
            AFBA1=REAL(AFB)
            AFB=0.0
         ELSE
            AFB=ZCODE(1)*REAL(AFB)
         ENDIF
         DO 50 I=1,NREG
            IBM=MAT(I)
            DO 45 IEL=1,IELEM
            Q(IEL)=0.0
            DO 40 IL=0,NSCT-1
            Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,NLF-IP+1)/2.0
   40       CONTINUE
   45       CONTINUE
            IF(IBFP.EQ.1) THEN
*             Galerkin type
              FN=2.0D0*ESTOPW(IBM,1)/DELTAE
              FD=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              G1=ESTOPW(IBM,1)/ESTOPW(IBM,2)
              G2=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/ESTOPW(IBM,2)
            ELSE
*             Przybylski and Ligou type
              FN=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              FD=2.0D0*ESTOPW(IBM,2)/DELTAE
              G1=1.0D0
              G2=2.0D0
            ENDIF
            FEP(:IELEM)=QEXT0(:IELEM,NLF/2+IP,I)
            IF(IELEM.EQ.1) THEN
              Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*U(NLF-IP+1)*AFB
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*U(NLF-IP+1)
            ELSE IF(IELEM.EQ.2) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
              Q2(1,2)=2.0D0*SQRT(3.0D0)*U(NLF-IP+1)
              Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)-6.0D0*U(NLF-IP+1)
              Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1        U(NLF-IP+1)*AFB
            ELSE IF(IELEM.EQ.3) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*U(NLF-IP+1)
              Q2(1,2)=0.0
              Q2(1,3)=2.0D0*SQRT(5.0D0)*U(NLF-IP+1)
              Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*U(NLF-IP+1)*AFB
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
              Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(NLF-IP+1)
              Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
              Q2(3,1)=Q2(1,3)
              Q2(3,2)=Q2(2,3)
              Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)+10.0D0*U(NLF-IP+1)
              Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)+2.0D0*SQRT(5.0D0)*
     1        U(NLF-IP+1)*AFB
            ENDIF
            CALL ALSBD(IELEM,1,Q2,IER,IELEM)
            IF(IER.NE.0) CALL XABORT('SNFP1P: SINGULAR MATRIX(2).')
            IF(IELEM.EQ.1) THEN
              IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
              AFB=2.0D0*Q2(1,2)-AFB
              IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
            ELSE IF(IELEM.EQ.2) THEN
              AFB=AFB+2.0D0*SQRT(3.0D0)*Q2(2,3)
            ELSE IF(IELEM.EQ.3) THEN
              AFB=2.0D0*Q2(1,4)+2.0D0*SQRT(5.0D0)*Q2(3,4)-AFB
            ENDIF
            FEP(:IELEM)=G2*Q2(:IELEM,IELEM+1)-G1*FEP(:IELEM)
            FLUX0(:IELEM,NLF/2+IP,I)=REAL(FEP(:IELEM))/DELTAE
   50    CONTINUE
*----
*  BACKWARD SWEEP (FROM RIGHT TO LEFT).
*----
         AFBA2=0.0
         AFBA=0.0
         IF(NCODE(1).EQ.4) THEN
            AFBA2=REAL(AFB)
         ELSE
            AFBA=REAL(AFB)
         ENDIF
         AFB=1.0
         DO 70 I=NREG,1,-1
            IBM=MAT(I)
            DO 65 IEL=1,IELEM
            Q(IEL)=0.0
            DO 60 IL=0,NSCT-1
            Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,IP)/2.0
   60       CONTINUE
   65       CONTINUE
            IF(IBFP.EQ.1) THEN
*             Galerkin type
              FN=2.0D0*ESTOPW(IBM,1)/DELTAE
              FD=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              G1=ESTOPW(IBM,1)/ESTOPW(IBM,2)
              G2=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/ESTOPW(IBM,2)
            ELSE
*             Przybylski and Ligou type
              FN=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              FD=2.0D0*ESTOPW(IBM,2)/DELTAE
              G1=1.0D0
              G2=2.0D0
            ENDIF
            FEP(:IELEM)=QEXT0(:IELEM,NLF/2+IP,I)
            IF(IELEM.EQ.1) THEN
              Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)-2.0D0*U(IP)*AFB
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)-2.0D0*U(IP)
            ELSE IF(IELEM.EQ.2) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
              Q2(1,2)=2.0D0*SQRT(3.0D0)*U(IP)
              Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)+6.0D0*U(IP)
              Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1        U(IP)*AFB
            ELSE IF(IELEM.EQ.3) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)-2.0D0*U(IP)
              Q2(1,2)=0.0
              Q2(1,3)=-2.0D0*SQRT(5.0D0)*U(IP)
              Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)-2.0D0*U(IP)*AFB
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
              Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(IP)
              Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
              Q2(3,1)=Q2(1,3)
              Q2(3,2)=Q2(2,3)
              Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)-10.0D0*U(IP)
              Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)-2.0D0*SQRT(5.0D0)*
     1        U(IP)*AFB
            ENDIF
            CALL ALSBD(IELEM,1,Q2,IER,IELEM)
            IF(IER.NE.0) CALL XABORT('SNFP1P: SINGULAR MATRIX(3).')
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
   70    CONTINUE
*----
*  FORWARD SWEEP (FROM LEFT TO RIGHT).
*----
         IF(NCODE(1).EQ.4) THEN
            AFBB1=REAL(AFB)
            AFB1=AFBA1/(1.0+AFBA1-AFBB1)
            AFB=1.0
         ELSE
            AFB=ZCODE(1)*AFB
         ENDIF
         DO 90 I=1,NREG
            IBM=MAT(I)
            DO 85 IEL=1,IELEM
            Q(IEL)=0.0
            DO 80 IL=0,NSCT-1
            Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,NLF-IP+1)/2.0
   80       CONTINUE
   85       CONTINUE
            IF(IBFP.EQ.1) THEN
*             Galerkin type
              FN=2.0D0*ESTOPW(IBM,1)/DELTAE
              FD=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              G1=ESTOPW(IBM,1)/ESTOPW(IBM,2)
              G2=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/ESTOPW(IBM,2)
            ELSE
*             Przybylski and Ligou type
              FN=(ESTOPW(IBM,1)+ESTOPW(IBM,2))/DELTAE
              FD=2.0D0*ESTOPW(IBM,2)/DELTAE
              G1=1.0D0
              G2=2.0D0
            ENDIF
            FEP(:IELEM)=QEXT0(:IELEM,NLF/2+IP,I)
            IF(IELEM.EQ.1) THEN
              Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*U(NLF-IP+1)*AFB
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*U(NLF-IP+1)
            ELSE IF(IELEM.EQ.2) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
              Q2(1,2)=2.0D0*SQRT(3.0D0)*U(NLF-IP+1)
              Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)-6.0D0*U(NLF-IP+1)
              Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1        U(NLF-IP+1)*AFB
            ELSE IF(IELEM.EQ.3) THEN
              Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*U(NLF-IP+1)
              Q2(1,2)=0.0
              Q2(1,3)=2.0D0*SQRT(5.0D0)*U(NLF-IP+1)
              Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*U(NLF-IP+1)*AFB
              Q2(2,1)=Q2(1,2)
              Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
              Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(NLF-IP+1)
              Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
              Q2(3,1)=Q2(1,3)
              Q2(3,2)=Q2(2,3)
              Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)+10.0D0*U(NLF-IP+1)
              Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)+2.0D0*SQRT(5.0D0)*
     1        U(NLF-IP+1)*AFB
            ENDIF
            CALL ALSBD(IELEM,1,Q2,IER,IELEM)
            IF(IER.NE.0) CALL XABORT('SNFP1P: SINGULAR MATRIX(4).')
            IF(IELEM.EQ.1) THEN
              IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
              AFB=2.0D0*Q2(1,2)-AFB
              IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
            ELSE IF(IELEM.EQ.2) THEN
              AFB=AFB+2.0D0*SQRT(3.0D0)*Q2(2,3)
            ELSE IF(IELEM.EQ.3) THEN
              AFB=2.0D0*Q2(1,4)+2.0D0*SQRT(5.0D0)*Q2(3,4)-AFB
            ENDIF
            FEP(:IELEM)=G2*Q2(:IELEM,IELEM+1)-G1*FEP(:IELEM)
            FLUX0(:IELEM,NLF/2+IP,I)=REAL(FEP(:IELEM))/DELTAE
   90    CONTINUE
         IF(NCODE(1).EQ.4) THEN
            AFBB2=REAL(AFB)
            AFB2=AFBA2/(1.0+AFBA2-AFBB2)
         ELSE
            AFBB=REAL(AFB)
            AFB=ZCODE(2)*AFBA/(1.0+ZCODE(2)*(AFBA-AFBB))
         ENDIF
      ENDIF
*----
*  BACKWARD SWEEP (FROM RIGHT TO LEFT).
*----
      IF(NCODE(1).EQ.4) AFB=AFB1
      DO 120 I=NREG,1,-1
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
           Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)-2.0D0*U(IP)*AFB
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)-2.0D0*U(IP)
         ELSE IF(IELEM.EQ.2) THEN
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
           Q2(1,2)=2.0D0*SQRT(3.0D0)*U(IP)
           Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
           Q2(2,1)=Q2(1,2)
           Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)+6.0D0*U(IP)
           Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1     U(IP)*AFB
         ELSE IF(IELEM.EQ.3) THEN
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)-2.0D0*U(IP)
           Q2(1,2)=0.0
           Q2(1,3)=-2.0D0*SQRT(5.0D0)*U(IP)
           Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)-2.0D0*U(IP)*AFB
           Q2(2,1)=Q2(1,2)
           Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
           Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(IP)
           Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
           Q2(3,1)=Q2(1,3)
           Q2(3,2)=Q2(2,3)
           Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)-10.0D0*U(IP)
           Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)-2.0D0*SQRT(5.0D0)*
     1     U(IP)*AFB
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
     1   PL(K,IP)
  110    CONTINUE
  115    CONTINUE
  120 CONTINUE
      CURR(1)=CURR(1)+W(IP)*U(IP)*AFB
*----
*  FORWARD SWEEP (FROM LEFT TO RIGHT).
*----
      AFB=ZCODE(1)*AFB
      IF(NCODE(1).EQ.4) AFB=AFB2
      DO 150 I=1,NREG
         IBM=MAT(I)
         DO 135 IEL=1,IELEM
         Q(IEL)=0.0
         DO 130 IL=0,NSCT-1
         Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,NLF-IP+1)/2.0
  130    CONTINUE
  135    CONTINUE
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
         FEP(:IELEM)=QEXT0(:IELEM,NLF/2+IP,I)
         IF(IELEM.EQ.1) THEN
           Q2(1,2)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*U(NLF-IP+1)*AFB
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*U(NLF-IP+1)
         ELSE IF(IELEM.EQ.2) THEN
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)
           Q2(1,2)=2.0D0*SQRT(3.0D0)*U(NLF-IP+1)
           Q2(1,3)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)
           Q2(2,1)=Q2(1,2)
           Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)-6.0D0*U(NLF-IP+1)
           Q2(2,3)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*
     1     U(NLF-IP+1)*AFB
         ELSE IF(IELEM.EQ.3) THEN
           Q2(1,1)=(FD+TOTAL(IBM))*VOL(I)+2.0D0*U(NLF-IP+1)
           Q2(1,2)=0.0
           Q2(1,3)=2.0D0*SQRT(5.0D0)*U(NLF-IP+1)
           Q2(1,4)=FN*FEP(1)*VOL(I)+Q(1)*VOL(I)+2.0D0*U(NLF-IP+1)*AFB
           Q2(2,1)=Q2(1,2)
           Q2(2,2)=-(FD+TOTAL(IBM))*VOL(I)
           Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(NLF-IP+1)
           Q2(2,4)=-FN*FEP(2)*VOL(I)-Q(2)*VOL(I)
           Q2(3,1)=Q2(1,3)
           Q2(3,2)=Q2(2,3)
           Q2(3,3)=(FD+TOTAL(IBM))*VOL(I)+10.0D0*U(NLF-IP+1)
           Q2(3,4)=FN*FEP(3)*VOL(I)+Q(3)*VOL(I)+2.0D0*SQRT(5.0D0)*
     1     U(NLF-IP+1)*AFB
         ENDIF
         CALL ALSBD(IELEM,1,Q2,IER,IELEM)
         IF(IER.NE.0) CALL XABORT('SNFP1P: SINGULAR MATRIX(6).')
         IF(IELEM.EQ.1) THEN
           IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
           AFB=2.0D0*Q2(1,2)-AFB
           IF(LFIXUP.AND.(AFB.LE.RLOG)) AFB=0.0
         ELSE IF(IELEM.EQ.2) THEN
           AFB=AFB+2.0D0*SQRT(3.0D0)*Q2(2,3)
         ELSE IF(IELEM.EQ.3) THEN
           AFB=2.0D0*Q2(1,4)+2.0D0*SQRT(5.0D0)*Q2(3,4)-AFB
         ENDIF
         FEP(:IELEM)=G2*Q2(:IELEM,IELEM+1)-G1*FEP(:IELEM)
         FLUX0(:IELEM,NLF/2+IP,I)=REAL(FEP(:IELEM))/DELTAE
         DO 145 K=1,NSCT
         DO 140 IEL=1,IELEM
         FLUX(IEL,K,I)=FLUX(IEL,K,I)+W(NLF-IP+1)*REAL(Q2(IEL,IELEM+1))*
     1   PL(K,NLF-IP+1)
  140    CONTINUE
  145    CONTINUE
  150 CONTINUE
      CURR(2)=CURR(2)+W(NLF-IP+1)*U(NLF-IP+1)*AFB
*----
*  DISCONTINUOUS GALERKIN FINITE ELEMENT METHOD
*----
      ELSEIF(ISCHM.EQ.2) THEN
         CALL XABORT('SNFP1P: DISCONTINUOUS GALERKIN NOT IMPLEMENTED.')
      ENDIF
*----
*  END OF OUTER LOOP OVER MU LEVELS
*----
  200 CONTINUE
      DEALLOCATE(Q,Q2)
      RETURN
      END
