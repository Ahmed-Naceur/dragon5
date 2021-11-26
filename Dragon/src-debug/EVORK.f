*DECK EVORK
      SUBROUTINE EVORK(Y,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,MU1,IMA,MAXA,
     1 NSUPF,NFISS,KFISS,YSF,ADPL,BDPL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Fifth-order Runge-Kutta Cash-Karp step with monitoring of local
* truncation error to ensure accuracy and adjust stepsize.
* Special version for isotopic depletion calculations.
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
* Y       dependent variable vector.
* N       size of the dependent variable vector.
* X       independent variable.
* HTRY    stepsize to be attempted.
* EPS     required accuracy.
* YSCAL   vector against which the error is scaled.
* HDID    stepsize that was actually accomplished.
* HNEXT   estimated next stepsize.
* MU1     position of each diagonal element in vectors ADPL and ASS.
* IMA     position of the first non-zero column element in vectors
*         ADPL and ASS.
* MAXA    first dimension of matrix ADPL.
* NSUPF   number of depleting fission products.
* NFISS   number of fissile isotopes producing fission products.
* KFISS   position in chain of the fissile isotopes.
* YSF     components of the product of the fission yields and fission
*         rates.
* ADPL    depletion matrix components.
* BDPL    depletion source components.
*
*Reference:
* W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
* "Numerical recipes in Fortran, Second edition, Chapter 16,
* Cambridge, 1992.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER N,MU1(N),IMA(N),MAXA,NSUPF,NFISS,KFISS(NFISS)
      REAL Y(N),X,HTRY,EPS,YSCAL(N),HDID,HNEXT,YSF(NFISS,NSUPF,2),
     1 ADPL(MAXA,2),BDPL(N,2)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      PARAMETER (SAFETY=0.85,PGROW=-.2,PSHRNK=-.25,GROW=1.5,SHRNK=0.5)
      CHARACTER HSMG*131
      REAL, ALLOCATABLE, DIMENSION(:) :: YTEMP,YGAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(YTEMP(N),YGAR(N),AK(N,6))
*
      NSUPL=N-NSUPF
      H=HTRY
   10 CALL ALLUM(N,ADPL(1,1),Y,AK(1,1),MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),Y,YGAR,MU1,IMA)
      DO 25 I=1,NSUPF
      DO 20 J=1,NFISS
      AK(NSUPL+I,1)=AK(NSUPL+I,1)+YSF(J,I,1)*Y(KFISS(J))
      YGAR(NSUPL+I)=YGAR(NSUPL+I)+YSF(J,I,2)*Y(KFISS(J))
   20 CONTINUE
   25 CONTINUE
      DO 30 I=1,N
      AK(I,1)=AK(I,1)+BDPL(I,1)+X*(YGAR(I)+BDPL(I,2))
      YTEMP(I)=Y(I)+H*B21*AK(I,1)
   30 CONTINUE
*
      CALL ALLUM(N,ADPL(1,1),YTEMP,AK(1,2),MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),YTEMP,YGAR,MU1,IMA)
      DO 45 I=1,NSUPF
      DO 40 J=1,NFISS
      AK(NSUPL+I,2)=AK(NSUPL+I,2)+YSF(J,I,1)*YTEMP(KFISS(J))
      YGAR(NSUPL+I)=YGAR(NSUPL+I)+YSF(J,I,2)*YTEMP(KFISS(J))
   40 CONTINUE
   45 CONTINUE
      DO 50 I=1,N
      AK(I,2)=AK(I,2)+BDPL(I,1)+(X+A2*H)*(YGAR(I)+BDPL(I,2))
      YTEMP(I)=Y(I)+H*(B31*AK(I,1)+B32*AK(I,2))
   50 CONTINUE
*
      CALL ALLUM(N,ADPL(1,1),YTEMP,AK(1,3),MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),YTEMP,YGAR,MU1,IMA)
      DO 65 I=1,NSUPF
      DO 60 J=1,NFISS
      AK(NSUPL+I,3)=AK(NSUPL+I,3)+YSF(J,I,1)*YTEMP(KFISS(J))
      YGAR(NSUPL+I)=YGAR(NSUPL+I)+YSF(J,I,2)*YTEMP(KFISS(J))
   60 CONTINUE
   65 CONTINUE
      DO 70 I=1,N
      AK(I,3)=AK(I,3)+BDPL(I,1)+(X+A3*H)*(YGAR(I)+BDPL(I,2))
      YTEMP(I)=Y(I)+H*(B41*AK(I,1)+B42*AK(I,2)+B43*AK(I,3))
   70 CONTINUE
*
      CALL ALLUM(N,ADPL(1,1),YTEMP,AK(1,4),MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),YTEMP,YGAR,MU1,IMA)
      DO 85 I=1,NSUPF
      DO 80 J=1,NFISS
      AK(NSUPL+I,4)=AK(NSUPL+I,4)+YSF(J,I,1)*YTEMP(KFISS(J))
      YGAR(NSUPL+I)=YGAR(NSUPL+I)+YSF(J,I,2)*YTEMP(KFISS(J))
   80 CONTINUE
   85 CONTINUE
      DO 90 I=1,N
      AK(I,4)=AK(I,4)+BDPL(I,1)+(X+A4*H)*(YGAR(I)+BDPL(I,2))
      YTEMP(I)=Y(I)+H*(B51*AK(I,1)+B52*AK(I,2)+B53*AK(I,3)+B54*AK(I,4))
   90 CONTINUE
*
      CALL ALLUM(N,ADPL(1,1),YTEMP,AK(1,5),MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),YTEMP,YGAR,MU1,IMA)
      DO 105 I=1,NSUPF
      DO 100 J=1,NFISS
      AK(NSUPL+I,5)=AK(NSUPL+I,5)+YSF(J,I,1)*YTEMP(KFISS(J))
      YGAR(NSUPL+I)=YGAR(NSUPL+I)+YSF(J,I,2)*YTEMP(KFISS(J))
  100 CONTINUE
  105 CONTINUE
      DO 110 I=1,N
      AK(I,5)=AK(I,5)+BDPL(I,1)+(X+A5*H)*(YGAR(I)+BDPL(I,2))
      YTEMP(I)=Y(I)+H*(B61*AK(I,1)+B62*AK(I,2)+B63*AK(I,3)+B64*AK(I,4)+
     1 B65*AK(I,5))
  110 CONTINUE
*
      CALL ALLUM(N,ADPL(1,1),YTEMP,AK(1,6),MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),YTEMP,YGAR,MU1,IMA)
      DO 125 I=1,NSUPF
      DO 120 J=1,NFISS
      AK(NSUPL+I,6)=AK(NSUPL+I,6)+YSF(J,I,1)*YTEMP(KFISS(J))
      YGAR(NSUPL+I)=YGAR(NSUPL+I)+YSF(J,I,2)*YTEMP(KFISS(J))
  120 CONTINUE
  125 CONTINUE
      DO 130 I=1,N
      AK(I,6)=AK(I,6)+BDPL(I,1)+(X+A6*H)*(YGAR(I)+BDPL(I,2))
      YTEMP(I)=Y(I)+H*(C1*AK(I,1)+C3*AK(I,3)+C4*AK(I,4)+C6*AK(I,6))
      YGAR(I)=H*(DC1*AK(I,1)+DC3*AK(I,3)+DC4*AK(I,4)+DC5*AK(I,5)+
     1 DC6*AK(I,6))
  130 CONTINUE
*
      ERRMAX=0.0
      DO 140 I=1,N
      ERRMAX=MAX(ERRMAX,ABS(YGAR(I)/YSCAL(I)))
  140 CONTINUE
      ERRMAX=ERRMAX/EPS
      IF (ERRMAX.EQ.0.0) THEN
        HDID=H
        HNEXT=GROW*H
        X=X+H
        DO 150 I=1,N
        Y(I)=YTEMP(I)
  150   CONTINUE
        GO TO 170
      ELSE IF (ERRMAX.LE.1.0) THEN
        HDID=H
        HNEXT=MIN(GROW,SAFETY*(ERRMAX**PGROW))*H
        X=X+H
        DO 160 I=1,N
        Y(I)=YTEMP(I)
  160   CONTINUE
        GO TO 170
      ELSE
        H=MAX(SHRNK,SAFETY*(ERRMAX**PSHRNK))*H
        XNEW=X+H
        IF (X.EQ.XNEW) THEN
           WRITE(HSMG,'(35HEVORK: STEPSIZE NOT SIGNIFICANT (H=,1P,E11.4,
     1     6H HTRY=,E11.4,2H).)') H,HTRY
           CALL XABORT(HSMG)
        ENDIF
        GO TO 10
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  170 DEALLOCATE(AK,YGAR,YTEMP)
      RETURN
      END
