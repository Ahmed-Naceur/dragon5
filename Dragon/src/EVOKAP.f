*DECK EVOKAP
      SUBROUTINE EVOKAP(Y,N,X,HTRY,EPS,YSCAL,HDID,HNEXT,MU1,IMA,MAXA,
     1 NSUPF,NFISS,KFISS,YSF,ADPL,BDPL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Fourth-order Kaps-Rentrop step for integrating stiff O.D.E.'s, with
* monitoring of local truncation error to adjust stepsize.
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
* W.H. Press and S.A. Teukolsky, 'Integrating stiff ordinary differen-
* tial equations', Computers in physics, 3 (3), 88 (May/June 1989).
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
      PARAMETER (MAXTRY=40,SAFETY=0.85,GROW=1.5,PGROW=-.25,SHRNK=0.5,
     1 PSHRNK=-1./3.)
      PARAMETER (GAM=.231,GAM21=-.270629667752/GAM,
     1 GAM31=.311254483294/GAM,GAM32=.852445628482E-2/GAM,
     2 GAM41=.282816832044/GAM,GAM42=-.457959483281/GAM,
     3 GAM43=-.111208333333/GAM,ALF21=.462,ALF31=-.815668168327E-1,
     4 ALF32=.961775150166,C1=.217487371653,C2=.486229037990,C3=0.,
     5 C4=.296283590357,CC1=-.717088504499,CC2=1.77617912176,
     6 CC3=-.590906172617E-1,GAM2X=GAM*(1.+GAM21),
     7 GAM3X=GAM*(1.+GAM31+GAM32),GAM4X=GAM*(1.+GAM41+GAM42+GAM43))
      CHARACTER HSMG*131
      REAL, ALLOCATABLE, DIMENSION(:) :: DYDX,TEMP,YSAV,DYSAV,DFDX,ASS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DYDX(N),TEMP(N),YSAV(N),DYSAV(N),DFDX(N),AK(N,4),
     1 ASS(IMA(N)))
*
      XSAV=X
      NSUPL=N-NSUPF
      CALL ALLUM(N,ADPL(1,1),Y,DYDX,MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),Y,TEMP,MU1,IMA)
      DO 15 I=1,NSUPF
      DO 10 J=1,NFISS
      DYDX(NSUPL+I)=DYDX(NSUPL+I)+YSF(J,I,1)*Y(KFISS(J))
      TEMP(NSUPL+I)=TEMP(NSUPL+I)+YSF(J,I,2)*Y(KFISS(J))
   10 CONTINUE
   15 CONTINUE
      DO 20 I=1,N
      DYDX(I)=DYDX(I)+X*TEMP(I)+BDPL(I,1)+X*BDPL(I,2)
      YSAV(I)=Y(I)
      DYSAV(I)=DYDX(I)
   20 CONTINUE
      H=HTRY
      DO 200 JTRY=1,MAXTRY
      HSQ=H*H
      CALL ALLUM(N,ADPL(1,2),YSAV,DFDX,MU1,IMA)
      DO 35 I=1,NSUPF
      DO 30 J=1,NFISS
      DFDX(NSUPL+I)=DFDX(NSUPL+I)+YSF(J,I,2)*YSAV(KFISS(J))
   30 CONTINUE
   35 CONTINUE
      DO 40 I=1,IMA(N)
      ASS(I)=-H*GAM*(ADPL(I,1)+XSAV*ADPL(I,2))
   40 CONTINUE
      DO 50 I=1,N
      DFDX(I)=DFDX(I)+BDPL(I,2)
      ASS(MU1(I))=1.+ASS(MU1(I))
   50 CONTINUE
      CALL ALLUF(N,ASS,MU1,IMA)
      DO 60 I=1,N
      AK(I,1)=H*DYSAV(I)+HSQ*GAM*DFDX(I)
   60 CONTINUE
      CALL ALLUS(NSUPL,MU1,IMA,ASS,AK(1,1))
      IF(NSUPF.GT.0) THEN
         DO 75 I=1,NSUPF
         DO 70 J=1,NFISS
         AK(NSUPL+I,1)=AK(NSUPL+I,1)+H*GAM*(YSF(J,I,1)+XSAV*YSF(J,I,2))
     1   *AK(KFISS(J),1)
   70    CONTINUE
   75    CONTINUE
         CALL ALLUS(NSUPF,MU1(NSUPL+1),IMA(NSUPL+1),ASS,AK(NSUPL+1,1))
      ENDIF
      DO 80 I=1,N
      Y(I)=YSAV(I)+ALF21*AK(I,1)
   80 CONTINUE
      X=XSAV+ALF21*H
      CALL ALLUM(N,ADPL(1,1),Y,DYDX,MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),Y,TEMP,MU1,IMA)
      DO 95 I=1,NSUPF
      DO 90 J=1,NFISS
      DYDX(NSUPL+I)=DYDX(NSUPL+I)+YSF(J,I,1)*Y(KFISS(J))
      TEMP(NSUPL+I)=TEMP(NSUPL+I)+YSF(J,I,2)*Y(KFISS(J))
   90 CONTINUE
   95 CONTINUE
      DO 100 I=1,N
      DYDX(I)=DYDX(I)+X*TEMP(I)+BDPL(I,1)+X*BDPL(I,2)
      AK(I,2)=H*DYDX(I)+HSQ*GAM2X*DFDX(I)+GAM21*AK(I,1)
  100 CONTINUE
      CALL ALLUS(NSUPL,MU1,IMA,ASS,AK(1,2))
      IF(NSUPF.GT.0) THEN
         DO 106 I=1,NSUPF
         DO 105 J=1,NFISS
         AK(NSUPL+I,2)=AK(NSUPL+I,2)+H*GAM*(YSF(J,I,1)+XSAV*YSF(J,I,2))
     1   *AK(KFISS(J),2)
  105    CONTINUE
  106    CONTINUE
         CALL ALLUS(NSUPF,MU1(NSUPL+1),IMA(NSUPL+1),ASS,AK(NSUPL+1,2))
      ENDIF
      DO 110 I=1,N
      AK(I,2)=AK(I,2)-GAM21*AK(I,1)
      Y(I)=YSAV(I)+ALF31*AK(I,1)+ALF32*AK(I,2)
  110 CONTINUE
      X=XSAV+(ALF31+ALF32)*H
      CALL ALLUM(N,ADPL(1,1),Y,DYDX,MU1,IMA)
      CALL ALLUM(N,ADPL(1,2),Y,TEMP,MU1,IMA)
      DO 125 I=1,NSUPF
      DO 120 J=1,NFISS
      DYDX(NSUPL+I)=DYDX(NSUPL+I)+YSF(J,I,1)*Y(KFISS(J))
      TEMP(NSUPL+I)=TEMP(NSUPL+I)+YSF(J,I,2)*Y(KFISS(J))
  120 CONTINUE
  125 CONTINUE
      DO 130 I=1,N
      DYDX(I)=DYDX(I)+X*TEMP(I)+BDPL(I,1)+X*BDPL(I,2)
      TEMP(I)=GAM31*AK(I,1)+GAM32*AK(I,2)
      AK(I,3)=H*DYDX(I)+GAM3X*HSQ*DFDX(I)+TEMP(I)
  130 CONTINUE
      CALL ALLUS(NSUPL,MU1,IMA,ASS,AK(1,3))
      IF(NSUPF.GT.0) THEN
         DO 136 I=1,NSUPF
         DO 135 J=1,NFISS
         AK(NSUPL+I,3)=AK(NSUPL+I,3)+H*GAM*(YSF(J,I,1)+XSAV*YSF(J,I,2))
     1   *AK(KFISS(J),3)
  135    CONTINUE
  136    CONTINUE
         CALL ALLUS(NSUPF,MU1(NSUPL+1),IMA(NSUPL+1),ASS,AK(NSUPL+1,3))
      ENDIF
      DO 140 I=1,N
      AK(I,3)=AK(I,3)-TEMP(I)
      TEMP(I)=GAM41*AK(I,1)+GAM42*AK(I,2)+GAM43*AK(I,3)
      AK(I,4)=H*DYDX(I)+HSQ*GAM4X*DFDX(I)+TEMP(I)
  140 CONTINUE
      CALL ALLUS(NSUPL,MU1,IMA,ASS,AK(1,4))
      IF(NSUPF.GT.0) THEN
         DO 146 I=1,NSUPF
         DO 145 J=1,NFISS
         AK(NSUPL+I,4)=AK(NSUPL+I,4)+H*GAM*(YSF(J,I,1)+XSAV*YSF(J,I,2))
     1   *AK(KFISS(J),4)
  145    CONTINUE
  146    CONTINUE
         CALL ALLUS(NSUPF,MU1(NSUPL+1),IMA(NSUPL+1),ASS,AK(NSUPL+1,4))
      ENDIF
      DO 150 I=1,N
      AK(I,4)=AK(I,4)-TEMP(I)
      Y(I)=YSAV(I)+C1*AK(I,1)+C2*AK(I,2)+C3*AK(I,3)+C4*AK(I,4)
      TEMP(I)=YSAV(I)+CC1*AK(I,1)+CC2*AK(I,2)+CC3*AK(I,3)
  150 CONTINUE
      X=XSAV+H
      IF (X.EQ.XSAV) THEN
         WRITE(HSMG,'(36HEVOKAP: STEPSIZE NOT SIGNIFICANT (H=,1P,E11.4,
     1   6H HTRY=,E11.4,2H).)') H,HTRY
         CALL XABORT(HSMG)
      ENDIF
      ERRMAX=0.
      DO 160 I=1,N
      ERRMAX=MAX(ERRMAX,ABS((Y(I)-TEMP(I))/YSCAL(I)))
  160 CONTINUE
      ERRMAX=ERRMAX/EPS
      IF (ERRMAX.EQ.0.) THEN
         HDID=H
         HNEXT=GROW*H
         GO TO 210
      ELSE IF (ERRMAX.LE.1.) THEN
         HDID=H
         HNEXT=MIN(GROW,SAFETY*(ERRMAX**PGROW))*H
         GO TO 210
      ELSE
         H=MAX(SHRNK,SAFETY*(ERRMAX**PSHRNK))*H
      ENDIF
  200 CONTINUE
      CALL XABORT('EVOKAP: EXCEEDED MAXTRY.')
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  210 DEALLOCATE(ASS,AK,DFDX,DYSAV,YSAV,TEMP,DYDX)
      RETURN
      END
