*DECK VALU5C
      SUBROUTINE VALU5C (KPMAC,LX,NUN,NMIX,X,XXX,EVT,ISS,IXLG,ITRIAL,
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
* NUN     dimension of unknown array EVT.
* NMIX    number of mixtures.
* X       Cartesian coordinates along the X axis where the flux is
*         interpolated.
* XXX     Cartesian coordinates along the X axis.
* EVT     variational coefficients of the flux.
* ISS     mixture index assigned to each element.
* IXLG    number of interpolated points according to X.
* ITRIAL  type of expansion functions in the nodal calculation
*         (=0: CMFD; =1: polynomial NEM; =2: hyperbolic NEM).
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
      INTEGER LX,NUN,NMIX,ISS(LX),IXLG,ITRIAL
      REAL X(IXLG),XXX(LX+1),EVT(NUN),AXY(IXLG)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION  WORK1(4,5),WORK2(2,3)
      DOUBLE PRECISION GAR,ETA,ALP1,COEF,U
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
*       Find the node index containing the interpolation point
        IS=0
        DO 20 KEL=1,LX
          IS=KEL
          IF((ABSC.GE.XXX(KEL)).AND.(ABSC.LE.XXX(KEL+1))) GO TO 30
   20   CONTINUE
        CALL XABORT('VALU5C: WRONG INTERPOLATION.')
   30   IBM=ISS(IS)
        IF(IBM.EQ.0) GO TO 100
        ETA=(XXX(IS+1)-XXX(IS))*SQRT(SIGR(IBM)/DIFF(IBM))
        ALP1=ETA*COSH(ETA/2.0)-2.0*SINH(ETA/2.0)
        COEF=DIFF(IBM)/(XXX(IS+1)-XXX(IS))
        U=(ABSC-XXX(IS))/(XXX(IS+1)-XXX(IS))-0.5
        IF(ITRIAL.EQ.0) THEN
          WORK2(1,1)=COEF
          WORK2(1,2)=-3.0*COEF
          WORK2(1,3)=EVT(3*LX+IS)
          WORK2(2,1)=COEF
          WORK2(2,2)=3.0*COEF
          WORK2(2,3)=EVT(3*LX+IS+1)
          CALL ALSBD(3,1,WORK2,IER,3)
          IF(IER.NE.0) CALL XABORT('VALU5C: SINGULAR MATRIX(1).')
          GAR=EVT(IS)+WORK2(1,3)*U+WORK2(2,3)*(3.0*U**2-1.0/4.0)
        ELSE
          WORK1(:,:)=0.0
          WORK1(1,1)=-0.5
          WORK1(1,2)=0.5
          WORK1(1,5)=EVT(LX+IS)-EVT(IS)
          WORK1(2,1)=0.5
          WORK1(2,2)=0.5
          WORK1(2,5)=EVT(2*LX+IS)-EVT(IS)
          WORK1(3,1)=-COEF
          WORK1(3,2)=3.0*COEF
          WORK1(3,5)=EVT(3*LX+IS)
          WORK1(4,1)=-COEF
          WORK1(4,2)=-3.0*COEF
          WORK1(4,5)=EVT(3*LX+IS+1)
          IF(ITRIAL.EQ.1) THEN
            WORK1(3,3)=-0.5*COEF
            WORK1(3,4)=0.2*COEF
            WORK1(4,3)=-0.5*COEF
            WORK1(4,4)=-0.2*COEF
          ELSE
            WORK1(1,3)=-SINH(ETA/2.0)
            WORK1(1,4)=ALP1/ETA
            WORK1(2,3)=SINH(ETA/2.0)
            WORK1(2,4)=ALP1/ETA
            WORK1(3,3)=-COEF*ETA*COSH(ETA/2.0)
            WORK1(3,4)=COEF*ETA*SINH(ETA/2.0)
            WORK1(4,3)=-COEF*ETA*COSH(ETA/2.0)
            WORK1(4,4)=-COEF*ETA*SINH(ETA/2.0)
          ENDIF
          CALL ALSBD(4,1,WORK1,IER,4)
          IF(IER.NE.0) CALL XABORT('VALU5C: SINGULAR MATRIX(2).')
          GAR=EVT(IS)+WORK1(1,5)*U+WORK1(2,5)*(3.0*U**2-1.0/4.0)
          IF(ITRIAL.EQ.1) THEN
            GAR=GAR+WORK1(3,5)*(U**2-0.25)*U+
     1      WORK1(4,5)*(U**2-0.25)*(U**2-0.05)
          ELSE
            GAR=GAR+WORK1(3,5)*SINH(ETA*U)+
     1      WORK1(4,5)*(COSH(ETA*U)-2.0*SINH(ETA/2.0)/ETA)
          ENDIF
        ENDIF
  100   AXY(I)=REAL(GAR)
      ENDDO
      DEALLOCATE(SIGW,SIGR,DIFF)
      RETURN
      END
