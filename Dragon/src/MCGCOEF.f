*DECK MCGCOEF
      SUBROUTINE MCGCOEF(NFUNL,NMU,ZMU,WZMU,NANGL,CAZ1,CAZ2,COEFI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find MOC coefficients for DD1 approximation up to P3 scattering order
* using a Gauss-Chebyshev or Leonard-McDaniel-Chebyshev quadrature.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* NFUNL   number of spherical harmonics components.
* NMU     number of polar angles.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* NANGL   number of azimuthal angles.
* CAZ1    first azimuthal cosines.
* CAZ2    second azimuthal cosines.
*
*Parameters: output
* COEFI   Gram-Schmidt MOC coefficients.
*
*Reference:
* A. Hebert, "High-Order Diamond Differencing Along Finite
* Characteristics," Nucl. Sci. Eng., 169, 81-97 (2011).
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFUNL,NMU
      REAL ZMU(NMU),WZMU(NMU)
      DOUBLE PRECISION CAZ1(NANGL),CAZ2(NANGL),COEFI(2*NFUNL,2*NFUNL)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: COEF,COEF2,COEF3
      DOUBLE PRECISION, DIMENSION(2,2) :: M
      INTEGER, DIMENSION(10), PARAMETER :: 
     >                        LL=(/ 0, 1, 1, 2, 2, 2, 3, 3, 3, 3 /)
      INTEGER, DIMENSION(10), PARAMETER ::
     >                        MM=(/ 0, -1, 1, -2, 0, 2, -3, -1, 1, 3 /)
*
      ALLOCATE(COEF(2*NFUNL,2*NFUNL))
      DO I=1,NFUNL
        IL=LL(I) ; IM=MM(I) ;
        DO J=1,NFUNL
          ILP=LL(J) ; IMP=MM(J) ;
          CALL MCGDYA(NMU,ZMU,WZMU,NANGL,CAZ1,CAZ2,IL,IM,ILP,IMP,M) ;
          COEF(I,J)=M(1,1) ; COEF(I,NFUNL+J)=M(1,2) ;
          COEF(NFUNL+I,J)=M(2,1) ; COEF(NFUNL+I,NFUNL+J)=M(2,2) ;
        ENDDO
      ENDDO
      IF(NFUNL == 1) THEN
        CALL ALPINVD(2*NFUNL,2*NFUNL,COEF,COEFI) ;
      ELSEIF(NFUNL == 3) THEN
        ALLOCATE(COEF2(2*NFUNL+1,2*NFUNL), COEF3(2*NFUNL,2*NFUNL+1))
        COEF2=0.0D0 ; COEF2(1:2*NFUNL,:)=COEF ;
        COEF2(2*NFUNL+1,2)=1.0D0 ; COEF2(2*NFUNL+1,NFUNL+3)=-1.0D0 ;
        CALL ALPINVD(2*NFUNL+1,2*NFUNL,COEF2,COEF3) ;
        COEFI=COEF3(:,1:2*NFUNL) ;
        DEALLOCATE(COEF2, COEF3)
      ELSEIF(NFUNL == 6) THEN
        ALLOCATE(COEF2(2*NFUNL+3,2*NFUNL), COEF3(2*NFUNL,2*NFUNL+3))
        COEF2=0.0D0 ; COEF2(1:2*NFUNL,:)=COEF ;
        COEF2(2*NFUNL+1,2)=1.0D0 ; COEF2(2*NFUNL+1,NFUNL+3)=-1.0D0 ;
        COEF2(2*NFUNL+2,5)=SQRT(3.0D0) ; COEF2(2*NFUNL+2,6)=1.0D0 ;
        COEF2(2*NFUNL+2,NFUNL+4)=1.0D0 ;
        COEF2(2*NFUNL+3,NFUNL+5)=SQRT(3.0D0) ;
        COEF2(2*NFUNL+3,NFUNL+6)=-1.0D0 ; COEF2(2*NFUNL+3,4)=1.0D0 ;
        CALL ALPINVD(2*NFUNL+3,2*NFUNL,COEF2,COEF3) ;
        COEFI=COEF3(:,1:2*NFUNL) ;
        DEALLOCATE(COEF2, COEF3)
      ELSEIF(NFUNL == 10) THEN
        ALLOCATE(COEF2(2*NFUNL+6,2*NFUNL), COEF3(2*NFUNL,2*NFUNL+6))
        COEF2=0.0D0 ; COEF2(1:2*NFUNL,:)=COEF ;
        COEF2(2*NFUNL+1,2)=1.0D0 ; COEF2(2*NFUNL+1,NFUNL+3)=-1.0D0 ;
        COEF2(2*NFUNL+2,5)=SQRT(3.0D0) ; COEF2(2*NFUNL+2,6)=1.0D0 ;
        COEF2(2*NFUNL+2,NFUNL+4)=1.0D0 ;
        COEF2(2*NFUNL+3,NFUNL+5)=SQRT(3.0D0) ;
        COEF2(2*NFUNL+3,NFUNL+6)=-1.0D0 ; COEF2(2*NFUNL+3,4)=1.0D0 ;
        COEF2(2*NFUNL+4,8)=1.0D0 ; COEF2(2*NFUNL+4,NFUNL+9)=-1.0D0 ;
        COEF2(2*NFUNL+5,7)=1.0D0 ; COEF2(2*NFUNL+5,NFUNL+10)=-1.0D0 ;
        COEF2(2*NFUNL+6,10)=1.0D0 ;
        COEF2(2*NFUNL+6,9)=SQRT(5.0D0/3.0D0) ;
        COEF2(2*NFUNL+6,NFUNL+8)=-SQRT(5.0D0/3.0D0) ;
        COEF2(2*NFUNL+6,NFUNL+7)=1.0D0 ;
        CALL ALPINVD(2*NFUNL+6,2*NFUNL,COEF2,COEF3) ;
        COEFI=COEF3(:,1:2*NFUNL) ;
        DEALLOCATE(COEF2, COEF3)
      ELSE
        CALL XABORT('MCGCOEF: NFUNL MUST BE = 1, 3, 6 OR 10')
      ENDIF
      DEALLOCATE(COEF)
      RETURN
      END
