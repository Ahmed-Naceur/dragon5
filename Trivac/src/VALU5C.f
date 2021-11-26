*DECK VALU5C
      SUBROUTINE VALU5C (KPMAC,LX,L4,NMIX,X,XXX,EVT,ISS,IXLG,ITRIAL,
     1 AXY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution for nodal expansion method in 1D.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* KPMAC   group directory in the macrolib.
* LX      number of elements along the X axis.
* L4      dimension of unknown array EVT.
* NMIX    number of mixtures.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* EVT     variational coefficients of the flux.
* ISS     mixture index assigned to each element.
* IXLG    number of interpolated points according to X.
* ITRIAL  type of expansion functions in the nodal calculation
*         (=1: polynomial; =2: hyperbolic).
*                                                                      
*Parameters: output
* AXY     interpolated fluxes.
*                                                                      
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPMAC
      INTEGER LX,L4,NMIX,ISS(LX),IXLG,ITRIAL
      REAL X(IXLG),XXX(LX+1),EVT(L4),AXY(IXLG)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: DIFF,SIGR,SIGW
*----
*  RECOVER REMOVAL CROSS SECTIONS AND DIFFUSION COEFFICIENTS
*----
      ALLOCATE(DIFF(NMIX),SIGR(NMIX),SIGW(NMIX))
      CALL LCMGET(KPMAC,'NTOT0',SIGR)
      CALL LCMGET(KPMAC,'SIGW00',SIGW)
      CALL LCMGET(KPMAC,'DIFF',DIFF)
      SIGR(:)=SIGR(:)-SIGW(:)
*----
*  PERFORM INTERPOLATION
*----
      DO I=1,IXLG
        ABSC=X(I)
        GAR=0.0
*                                                          
*       Find the finite element index containing the interpolation point
        IS=0
        DO 20 KEL=1,LX
          IS=KEL
          IF((ABSC.GE.XXX(KEL)).AND.(ABSC.LE.XXX(KEL+1))) GO TO 30
   20   CONTINUE
        CALL XABORT('VALU5C: WRONG INTERPOLATION.')
   30   IBM=ISS(IS)
        IF(IBM.EQ.0) GO TO 100
        ETA=(XXX(IS+1)-XXX(IS))*SQRT(SIGR(IBM)/DIFF(IBM))
        U=(ABSC-0.5*(XXX(IS)+XXX(IS+1)))/(XXX(IS+1)-XXX(IS))
        U=(ABSC-XXX(IS))/(XXX(IS+1)-XXX(IS))-0.5
        GAR=EVT((IS-1)*5+1)+EVT((IS-1)*5+2)*U+EVT((IS-1)*5+3)*
     1  (U**2-1.0/12.0)
        IF(ITRIAL.EQ.1) THEN
          GAR=GAR+EVT((IS-1)*5+4)*(U**2-0.25)*U+
     1    EVT((IS-1)*5+5)*(U**2-0.25)*(U**2-0.05)
        ELSE
          GAR=GAR+EVT((IS-1)*5+4)*SINH(ETA*U)+
     1    EVT((IS-1)*5+5)*(COSH(ETA*U)-2.0*SINH(ETA/2.0)/ETA)
        ENDIF
  100   AXY(I)=GAR
      ENDDO
      DEALLOCATE(SIGW,SIGR,DIFF)
      RETURN
      END
