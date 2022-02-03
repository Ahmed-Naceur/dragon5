*DECK EVOODE
      SUBROUTINE EVOODE(YSTART,NVAR,X1,X2,EPS,H1,NOK,NBAD,ITYPE,MU1,
     1 IMA,MAXA,NSUPF,NFISS,KFISS,YSF,ADPL,BDPL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Runge-Kutta or Kaps-Rentrop driver with adaptive stepsize control
* special version for isotopic depletion calculations.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* YSTART  dependent variable vector.
* NVAR    dimension of the dependent variable vector (number of
*         depleting isotopes).
* X1      initial value of the independent variable.
* X2      final value of the independent variable.
* EPS     required accuracy.
* H1      guessed first stepsize.
* NOK     number of good steps taken.
* NBAD    number of bad steps taken.
* ITYPE   type of ODE solution:
*         =1 fifth-order Runge-Kutta method;
*         =2 fourth-order Kaps-Rentrop method.
* MU1     position of each diagonal element in matrix ADPL.
* IMA     position of the first non-zero column element in matrix ADPL.
* MAXA    first dimension of matrix ADPL.
* NSUPF   number of depleting fission products.
* NFISS   number of fissile isotopes producing fission products.
* KFISS   position in chain of the fissile isotopes.
* YSF     components of the product of the fission yields and fission
*         rates.
* ADPL    depletion matrix components.
* BDPL    depletion source components.
*
*-----------------------------------------------------------------------
*
* REFERENCE:
*  W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY AND W.T. VETTERLING,
*  'NUMERICAL RECIPIES (FORTRAN VERSION)', CAMBRIDGE UNIVERSITY PRESS,
*  CHAPTER 15, CAMBRIDGE (1990).
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NVAR,NOK,NBAD,ITYPE,MU1(NVAR),IMA(NVAR),MAXA,NSUPF,NFISS,
     1 KFISS(NFISS)
      REAL YSTART(NVAR),X1,X2,EPS,H1,YSF(NFISS,NSUPF,2),ADPL(MAXA,2),
     1 BDPL(NVAR,2)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXSTP=50000,TINY=1.E-8)
      REAL, ALLOCATABLE, DIMENSION(:) :: YSCAL,Y
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(YSCAL(NVAR),Y(NVAR))
*
      NSUPL=NVAR-NSUPF
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0
      NBAD=0
      DO 10 I=1,NVAR
      Y(I)=YSTART(I)
   10 CONTINUE
      DO 50 NSTP=1,MAXSTP
        IF((X+H-X2)*(X+H-X1).GT.0.0) H=X2-X
        CALL ALLUM(NVAR,ADPL(1,1),Y,YSCAL,MU1,IMA)
        CALL ALLUM(NVAR,ADPL(1,2),Y,YSTART,MU1,IMA)
        DO 25 I=1,NSUPF
        DO 20 J=1,NFISS
        YSCAL(NSUPL+I)=YSCAL(NSUPL+I)+YSF(J,I,1)*Y(KFISS(J))
        YSTART(NSUPL+I)=YSTART(NSUPL+I)+YSF(J,I,2)*Y(KFISS(J))
   20   CONTINUE
   25   CONTINUE
        DO 30 I=1,NVAR
        YSCAL(I)=YSCAL(I)+BDPL(I,1)+X*(YSTART(I)+BDPL(I,2))
        YSCAL(I)=MAX(ABS(Y(I))+ABS(H*YSCAL(I)),TINY)
   30   CONTINUE
        IF(ITYPE.EQ.1) THEN
           CALL EVORK(Y,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,MU1,IMA,MAXA,
     1     NSUPF,NFISS,KFISS,YSF,ADPL,BDPL)
        ELSE IF(ITYPE.EQ.2) THEN
           CALL EVOKAP(Y,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,MU1,IMA,MAXA,
     1     NSUPF,NFISS,KFISS,YSF,ADPL,BDPL)
        ENDIF
        IF(HDID.EQ.H) THEN
          NOK=NOK+1
        ELSE
          NBAD=NBAD+1
        ENDIF
        IF((X-X2)*(X2-X1).GE.0.0) THEN
          DO 40 I=1,NVAR
          YSTART(I)=Y(I)
   40     CONTINUE
          GO TO 60
        ENDIF
        H=HNEXT
   50 CONTINUE
      CALL XABORT('EVOODE: TOO MANY STEPS.')
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   60 DEALLOCATE(Y,YSCAL)
      END
