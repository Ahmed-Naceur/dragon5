*DECK SNFT1Q
      SUBROUTINE SNFT1Q(NREG,NMAT,ISCHM,IELEM,NLF,NSCT,U,W,PL,MAT,VOL,
     1 TOTAL,NCODE,ZCODE,QEXT,LFIXUP,FLUX,XNEI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 1D slab
* geometry. The shooting method is NOT used. LEGENDRE POLYNOMIALS 
* are used for the solution expansion. Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* NREG    number of regions.
* NMAT    number of material mixtures.
* ISCHM   method of spatial discretisation:
*         =1 High-Order Diamond Differencing (HODD) - default;
*         =2 Discontinuous Galerkin finite element method (DG).
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - classical diamond scheme, default for HODD;
*         =2 linear - default for DG;
*         =3 parabolic;
*         =4 cubic - only for DG.
* NLF     number of $\\mu$ levels.
* NSCT    number of Legendre components in the flux.
*         =1 isotropic sources;
*         =2 linearly anisotropic sources.
* U       base points in $\\mu$ of the SN quadrature.
* W       weights of the SN quadrature.
* PL      discrete values of the Legendre polynomials corresponding
*         to the SN quadrature.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* XNEI    not defined.
*
*Parameters: output
* FLUX    Legendre components of the flux.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NREG,NMAT,ISCHM,IELEM,NLF,NSCT,MAT(NREG),NCODE(2)
      REAL U(NLF),W(NLF),PL(NSCT,NLF),VOL(NREG),TOTAL(0:NMAT),ZCODE(2),
     1 QEXT(IELEM,NSCT,NREG),FLUX(IELEM,NSCT,NREG),XNEI(NLF)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION XNI
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
     1      CALL XABORT('SNFT1Q: INVALID IELEM (DIAM) VALUE. ' 
     2      //'CHECK INPUT DATA FILE.')
      ELSEIF(ISCHM.EQ.2) THEN
         IF((IELEM.LT.1).OR.(IELEM.GT.4)) 
     1      CALL XABORT('SNFT1Q: INVALID IELEM (DIAM) VALUE. '
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

      DO 200 M=1,NLF
*
      WEIGHT=W(M)
      IF(WEIGHT.EQ.0.0) GO TO 200
*
      IF(U(M).GT.0.0) THEN
         IF(NCODE(1).NE.4) XNEI(M)=XNEI(NLF-M+1)
      ELSEIF(U(M).LT.0.0) THEN
         IF(NCODE(2).NE.4) XNEI(M)=XNEI(NLF-M+1)
      ENDIF
*
      IND=2
      IF(U(M).GT.0.0) IND=1
      XNI=0.0D0
*
*----
*  HIGH-ORDER DIAMOND DIFFERENCING METHOD
*----
      IF(ISCHM.EQ.1) THEN 
      DO IZ=1,NREG
         I=IZ
         IF(IND.EQ.2) I=NREG+1-IZ
         IF(IZ.EQ.1) THEN
            IF(IND.EQ.1) THEN
               XNI=XNEI(M)*ZCODE(1)
            ELSE
               XNI=XNEI(M)*ZCODE(2)
            ENDIF
         ENDIF
*
         IBM=MAT(I)
         DO IEL=1,IELEM
            Q(IEL)=0.0
            DO IL=0,NSCT-1
               Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,M)/2.0
            ENDDO
         ENDDO
         CALL XDDSET(Q2,IELEM*(IELEM+1),0.0D0)
*
         IF(IELEM.EQ.1) THEN
            Q2(1,1)=TOTAL(IBM)*VOL(I)+2.0D0*ABS(U(M))
            Q2(1,2)=Q(1)*VOL(I)+2.0D0*ABS(U(M))*XNI
         ELSE IF(IELEM.EQ.2) THEN
            Q2(1,1)=TOTAL(IBM)*VOL(I)
            Q2(1,2)=2.0D0*SQRT(3.0D0)*U(M)
            Q2(1,3)=Q(1)*VOL(I)
            Q2(2,1)=Q2(1,2)
            Q2(2,2)=-TOTAL(IBM)*VOL(I)-6.0D0*ABS(U(M))
            Q2(2,3)=-Q(2)*VOL(I)+2.0D0*SQRT(3.0D0)*U(M)*XNI
         ELSE IF(IELEM.EQ.3) THEN
            Q2(1,1)=TOTAL(IBM)*VOL(I)+2.0D0*ABS(U(M))
            Q2(1,2)=0.0
            Q2(1,3)=2.0D0*SQRT(5.0D0)*ABS(U(M))
            Q2(1,4)=Q(1)*VOL(I)+2.0D0*ABS(U(M))*XNI
            Q2(2,1)=Q2(1,2)
            Q2(2,2)=-TOTAL(IBM)*VOL(I)
            Q2(2,3)=-2.0D0*SQRT(15.0D0)*U(M)
            Q2(2,4)=-Q(2)*VOL(I)
            Q2(3,1)=Q2(1,3)
            Q2(3,2)=Q2(2,3)
            Q2(3,3)=TOTAL(IBM)*VOL(I)+10.0D0*ABS(U(M))
            Q2(3,4)=Q(3)*VOL(I)+2.0D0*SQRT(5.0D0)*ABS(U(M))*XNI
         ENDIF
*
         CALL ALSBD(IELEM,1,Q2,IER,IELEM)
*
         IF(IER.NE.0) CALL XABORT('SNFT1D: SINGULAR MATRIX(1).')
         IF(IELEM.EQ.1) THEN
            IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
            XNI=2.0D0*Q2(1,2)-XNI
         ELSE IF(IELEM.EQ.2) THEN
            XNI=XNI+SIGN(1.0,U(M))*2.0D0*SQRT(3.0D0)*Q2(2,3)
         ELSE IF(IELEM.EQ.3) THEN
            XNI=2.0D0*Q2(1,4)+2.0D0*SQRT(5.0D0)*Q2(3,4)-XNI
         ENDIF
*
         DO K=1,NSCT
            DO IEL=1,IELEM
               FLUX(IEL,K,I)=FLUX(IEL,K,I)+W(M)*
     >          REAL(Q2(IEL,IELEM+1))*PL(K,M)
            ENDDO 
         ENDDO
      ENDDO
      XNEI(M)=REAL(XNI)
* -------------------------------------------------------------------- *
* -------------------------------------------------------------------- *
*----
*  DISCONTINUOUS GALERKIN FINITE ELEMENT METHOD
*----
      ELSEIF(ISCHM.EQ.2) THEN
      DO IZ=1,NREG
         I=IZ
         IF(IND.EQ.2) I=NREG+1-IZ
         IF(IZ.EQ.1) THEN
            IF(IND.EQ.1) THEN
               XNI=XNEI(M)*ZCODE(1)
            ELSE
               XNI=XNEI(M)*ZCODE(2)
            ENDIF
         ENDIF
*
         IBM=MAT(I)
         DO IEL=1,IELEM
            Q(IEL)=0.0
            DO IL=0,NSCT-1
               Q(IEL)=Q(IEL)+QEXT(IEL,IL+1,I)*PL(IL+1,M)/2.0
            ENDDO
         ENDDO
         CALL XDDSET(Q2,IELEM*(IELEM+1),0.0D0)
*
         IF(IELEM.EQ.1) THEN
            Q2(1,1) = ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,2) = Q(1)*VOL(I) 
            Q2(1,2) = Q2(1,2) + ABS(U(M))*XNI
         ELSE IF(IELEM.EQ.2) THEN
            Q2(1,1) = ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,2) = 3**(0.5D0)*U(M) 
            Q2(2,1) = -3**(0.5D0)*U(M) 
            Q2(2,2) = 3*ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,3) = Q(1)*VOL(I) 
            Q2(2,3) = Q(2)*VOL(I) 
            Q2(1,3) = Q2(1,3) + ABS(U(M))*XNI
            Q2(2,3) = Q2(2,3) - 3**(0.5D0)*U(M)*XNI
         ELSE IF(IELEM.EQ.3) THEN
            Q2(1,1) = ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,2) = 3**(0.5D0)*U(M) 
            Q2(1,3) = 5**(0.5D0)*ABS(U(M))
            Q2(2,1) = -3**(0.5D0)*U(M) 
            Q2(2,2) = 3*ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(2,3) = 2*15**(0.5D0)*U(M) - 3**(0.5D0)*5**(0.5D0)*U(M) 
            Q2(3,1) = 5**(0.5D0)*ABS(U(M))
            Q2(3,2) = -3**(0.5D0)*5**(0.5D0)*U(M) 
            Q2(3,3) = 5*ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,4) = Q(1)*VOL(I) 
            Q2(2,4) = Q(2)*VOL(I) 
            Q2(3,4) = Q(3)*VOL(I) 
            Q2(1,4) = Q2(1,4) + ABS(U(M))*XNI
            Q2(2,4) = Q2(2,4) - 3**(0.5D0)*U(M)*XNI
            Q2(3,4) = Q2(3,4) + 5**(0.5D0)*ABS(U(M))*XNI
         ELSE IF(IELEM.EQ.4) THEN
            Q2(1,1) = ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,2) = 3**(0.5D0)*U(M) 
            Q2(1,3) = 5**(0.5D0)*ABS(U(M))
            Q2(1,4) = 7**(0.5D0)*U(M) 
            Q2(2,1) = -3**(0.5D0)*U(M) 
            Q2(2,2) = 3*ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(2,3) = 2*15**(0.5D0)*U(M) - 3**(0.5D0)*5**(0.5D0)*U(M) 
            Q2(2,4) = 3**(0.5D0)*7**(0.5D0)*ABS(U(M))
            Q2(3,1) = 5**(0.5D0)*ABS(U(M))
            Q2(3,2) = -3**(0.5D0)*5**(0.5D0)*U(M) 
            Q2(3,3) = 5*ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(3,4) = 2*35**(0.5D0)*U(M) - 5**(0.5D0)*7**(0.5D0)*U(M) 
            Q2(4,1) = -7**(0.5D0)*U(M) 
            Q2(4,2) = 3**(0.5D0)*7**(0.5D0)*ABS(U(M))
            Q2(4,3) = -5**(0.5D0)*7**(0.5D0)*U(M) 
            Q2(4,4) = 7*ABS(U(M)) + TOTAL(IBM)*VOL(I) 
            Q2(1,5) = Q(1)*VOL(I) 
            Q2(2,5) = Q(2)*VOL(I) 
            Q2(3,5) = Q(3)*VOL(I) 
            Q2(4,5) = Q(4)*VOL(I) 
            Q2(1,5) = Q2(1,5) + ABS(U(M))*XNI
            Q2(2,5) = Q2(2,5) - 3**(0.5D0)*U(M)*XNI
            Q2(3,5) = Q2(3,5) + 5**(0.5D0)*ABS(U(M))*XNI
            Q2(4,5) = Q2(4,5) - 7**(0.5D0)*U(M)*XNI
         ENDIF
*
         CALL ALSBD(IELEM,1,Q2,IER,IELEM)
*
         IF(IER.NE.0) CALL XABORT('SNFT1D: SINGULAR MATRIX(1).')
         IF(IELEM.EQ.1) THEN
            IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
            XNI=2.0D0*Q2(1,2)-XNI
         ELSE IF(IELEM.EQ.2) THEN
            XNI=XNI+SIGN(1.0,U(M))*2.0D0*SQRT(3.0D0)*Q2(2,3)
         ELSE IF(IELEM.EQ.3) THEN
            XNI=2.0D0*Q2(1,4)+2.0D0*SQRT(5.0D0)*Q2(3,4)-XNI
         ENDIF
         IF(IELEM.EQ.1)THEN
            XNI = Q2(1,2)
         ELSEIF(IELEM.EQ.2)THEN
            XNI = Q2(1,3) + SIGN(1.0,U(M))*Q2(2,3)*3**(0.5D0)
         ELSEIF(IELEM.EQ.3)THEN
            XNI = Q2(1,4) + SIGN(1.0,U(M))*Q2(2,4)*3**(0.5D0) + 
     >         Q2(3,4)*5**(0.5D0)
         ELSEIF(IELEM.EQ.4)THEN
            XNI = Q2(1,5) + SIGN(1.0,U(M))*Q2(2,5)*3**(0.5D0) + 
     >         Q2(3,5)*5**(0.5D0) + SIGN(1.0,U(M))*Q2(4,5)*7**(0.5D0)
         ENDIF
*
         DO K=1,NSCT
            DO IEL=1,IELEM
               FLUX(IEL,K,I)=FLUX(IEL,K,I)+W(M)*
     >          REAL(Q2(IEL,IELEM+1))*PL(K,M)
            ENDDO 
         ENDDO
      ENDDO
      XNEI(M)=REAL(XNI)
      ENDIF
*----
*  END OF OUTER LOOP OVER MU LEVELS
*----
  200 CONTINUE
*
      DEALLOCATE(Q,Q2)
      RETURN
      END
