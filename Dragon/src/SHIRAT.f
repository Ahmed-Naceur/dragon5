*DECK SHIRAT
      SUBROUTINE SHIRAT(IMPX,NRAT,SIGX,DILUT,IGRP,SA,COEF,DENOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the NRAT-terms rational approximation coefficients for
* the SIGX-dependent fuel-to-fuel reduced collision probability in a
* closed cell.
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
*Parameters: input
* IMPX    print flag (no print if IMPX.lt.10).
* NRAT    number of terms in the pij rational approximation.
* SIGX    interpolation values for the resonant cross section of
*         the heavy nuclide.
* DILUT   interpolated macroscopic escape cross sections corresponding
*         to SIGX values.
* IGRP    group index.
*
*Parameters: output
* SA      asymptotic macroscopic escape cross section.
* COEF    numerator coefficients for the rational approximation
*         of fuel-to-fuel reduced collision probability.
* DENOM   denominator coefficients for the rational approximation
*         of fuel-to-fuel reduced collision probability.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,NRAT,IGRP
      REAL SIGX(2*NRAT-1),DILUT(2*NRAT-1),SA
      COMPLEX COEF(NRAT),DENOM(NRAT)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NMAX=11,EPSRID=1.0D-5,NORIN=(NMAX-1)/2)
      COMPLEX*16 E1,E2,E3,AAA,BBB,SQ1,TEST,D,E,SQRTM3,SIGXI,DD,
     1 CDENOM(NORIN+1)
      PARAMETER (SQRTM3=(0.0,1.73205080756888))
      DOUBLE PRECISION A(0:NORIN),B(0:NORIN),C(0:NORIN+1)
      CHARACTER HSMG*131
      REAL EPS,PREC
      COMPLEX EAV
      LOGICAL LFAIL
*
      IF(2*NRAT-1.GT.NMAX) CALL XABORT('SHIRAT: INCREASE NMAX.')
      DO 10 I=2,NRAT
      COEF(I)=0.0
      DENOM(I)=1.0
   10 CONTINUE
      EPS=0.0
      DO 15 I=1,2*NRAT-1,2
      EPS=MAX(EPS,ABS(DILUT(I)-DILUT(NRAT))/ABS(DILUT(I)))
   15 CONTINUE
      IF(EPS.LT.5.0E-4) THEN
         SA=DILUT(NRAT)
         COEF(1)=1.0
         DENOM(1)=DILUT(NRAT)
         IF(IMPX.GE.10) WRITE (6,110) SA
         RETURN
      ENDIF
*
      CALL ALPLSF(1,2*NRAT-1,SIGX,DILUT,EPSRID,.FALSE.,NOR,A,B,PREC)
      SA=REAL(A(NOR))
      QQ=A(1)+B(0)
      RR=A(0)
      IF(NOR.EQ.0) THEN
*        1-TERM RATIONAL APPROXIMATION.
         COEF(1)=1.0
         DENOM(1)=CMPLX(A(0),KIND=KIND(DENOM))
      ELSE IF(NOR.EQ.1) THEN
*        2-TERMS RATIONAL APPROXIMATION.
         AAA=QQ*QQ-4.0D0*RR
         AAA=SQRT(AAA)
         E1=0.5D0*(QQ+AAA)
         E2=0.5D0*(QQ-AAA)
         IF(ABS(DBLE(E1*E2)-RR).GT.1.0E-3*ABS(RR)) THEN
            WRITE (HSMG,'(42HSHIRAT: INTERPOLATION ALGORITHM FAILURE 1,,
     1      6H COEF=,1P,3E11.3)') QQ,RR,DBLE(E1*E2)
            CALL XABORT(HSMG)
         ENDIF
*
         COEF(1)=CMPLX((B(0)-E1)/(E2-E1))
         COEF(2)=CMPLX((B(0)-E2)/(E1-E2))
         DENOM(1)=CMPLX(E1)
         DENOM(2)=CMPLX(E2)
      ELSE IF(NOR.EQ.2) THEN
*        3-TERMS RATIONAL APPROXIMATION.
         PP=A(2)+B(1)
         AA=(3.0D0*QQ-PP**2)/3.0D0
         BB=(2.0D0*PP**3-9.0D0*PP*QQ+27.0D0*RR)/27.0D0
         SQ1=BB**2/4.0D0+AA**3/27.0D0
         TEST=BB/2.0D0-SQRT(SQ1)
         IF(DBLE(TEST).EQ.0.0) THEN
            AAA=0.0D0
         ELSE IF(DBLE(TEST).GT.0.0) THEN
            AAA=-(TEST)**(1.0D0/3.0D0)
         ELSE
            AAA=(-TEST)**(1.0D0/3.0D0)
         ENDIF
         TEST=BB/2.0D0+SQRT(SQ1)
         IF(DBLE(TEST).EQ.0.0) THEN
            BBB=0.0D0
         ELSE IF(DBLE(TEST).GT.0.0) THEN
            BBB=-(TEST)**(1.0D0/3.0D0)
         ELSE
            BBB=(-TEST)**(1.0D0/3.0D0)
         ENDIF
         E1=-(AAA+BBB-PP/3.0D0)
         E2=-(-(AAA+BBB)/2.0D0+(AAA-BBB)*SQRTM3/2.0D0-PP/3.0D0)
         E3=-(-(AAA+BBB)/2.0D0-(AAA-BBB)*SQRTM3/2.0D0-PP/3.0D0)
         IF(ABS(DBLE(E1*E2*E3)-RR).GT.1.0E-3*ABS(RR)) THEN
            WRITE (HSMG,'(42HSHIRAT: INTERPOLATION ALGORITHM FAILURE 2,,
     1      6H COEF=,1P,4E11.3)') PP,QQ,RR,DBLE(E1*E2*E3)
            CALL XABORT(HSMG)
         ENDIF
*
         SQ1=(0.5D0*B(1))**2-B(0)
         D=0.5D0*B(1)+SQRT(SQ1)
         E=0.5D0*B(1)-SQRT(SQ1)
         COEF(1)=CMPLX((D-E1)*(E-E1)/(E2-E1)/(E3-E1))
         COEF(2)=CMPLX((D-E2)*(E-E2)/(E1-E2)/(E3-E2))
         COEF(3)=CMPLX((D-E3)*(E-E3)/(E1-E3)/(E2-E3))
         DENOM(1)=CMPLX(E1)
         DENOM(2)=CMPLX(E2)
         DENOM(3)=CMPLX(E3)
      ELSE IF(NOR.GE.3) THEN
*        (NOR+1) TERMS RATIONAL APPROXIMATION.
         NORP1=NOR+1
         SGN=1.0D0
         C(0)=A(0)
         DO 25 I=2,NORP1
         SGN=-SGN
         C(I-1)=SGN*(B(I-2)+A(I-1))
   25    CONTINUE
         C(NORP1)=-SGN
         CALL ALROOT(C,NORP1,CDENOM,LFAIL)
         IF(LFAIL) CALL XABORT('SHIRAT: ROOT FINDING FAILURE.')
         DO 50 I=1,NORP1
         SIGXI=CDENOM(I)
         DENOM(I)=CMPLX(SIGXI)
         DD=SIGXI**(NORP1-1)
         SGN=1.0D0
         DO 30 J=NORP1-1,1,-1
         SGN=-SGN
         DD=DD+SGN*B(J-1)*SIGXI**(J-1)
   30    CONTINUE
         DO 40 J=1,NORP1
         IF(J.NE.I) DD=DD/(SIGXI-CDENOM(J))
   40    CONTINUE
         COEF(I)=CMPLX(DD)
   50    CONTINUE
      ELSE
         CALL XABORT('SHIRAT: PADE COLLOCATION FAILURE.')
      ENDIF
      IF(IMPX.GE.10) THEN
         WRITE (6,80) IGRP,(COEF(I),I=1,NOR+1)
         WRITE (6,90) (DENOM(I),I=1,NOR+1)
         WRITE (6,100)
         X=1.0D0
         DO 70 I=1,2*NRAT-1
         Z1=0.0D0
         Z2=0.0D0
         DO 60 J=0,NOR
         Z1=Z1+A(J)*X
         Z2=Z2+B(J)*X
         X=X*SIGX(I)
   60    CONTINUE
         WRITE (6,'(1X,I5,1P,3E13.5)') I,SIGX(I),DILUT(I),Z1/Z2
   70    CONTINUE
         WRITE (6,110) SA
      ENDIF
      EAV=0.0
      DO 75 I=1,NRAT
      EAV=EAV+COEF(I)*SQRT(DENOM(I))
   75 CONTINUE
      EAV=EAV*EAV
      IF(REAL(EAV).LT.0.0) THEN
         NALPHA=2*NRAT-1
         WRITE (6,120) (SIGX(I),I=1,NALPHA)
         WRITE (6,130) (DILUT(I),I=1,NALPHA)
         WRITE (6,80) IGRP,(COEF(I),I=1,NRAT)
         WRITE (6,90) (DENOM(I),I=1,NRAT)
         WRITE (HSMG,'(41HSHIRAT: RATIONAL EXPANSION FAILURE. EAV=(,
     1   1P,E10.3,1H,,E10.3,33H) HAS NEGATIVE REAL PART IN GROUP,I4,
     2   1H.)') EAV,IGRP
         CALL XABORT(HSMG)
      ELSE IF(ABS(AIMAG(EAV)).GT.5.0E-3*REAL(EAV)) THEN
         NALPHA=2*NRAT-1
         WRITE (6,120) (SIGX(I),I=1,NALPHA)
         WRITE (6,130) (DILUT(I),I=1,NALPHA)
         WRITE (6,80) IGRP,(COEF(I),I=1,NRAT)
         WRITE (6,90) (DENOM(I),I=1,NRAT)
         WRITE (6,'(/42H SHIRAT: RATIONAL EXPANSION WARNING. EAV=(,
     1   1P,E10.3,1H,,E10.3,35H) HAS LARGE IMAGINARY PART IN GROUP,
     2   I4,1H.)') EAV,IGRP
      ENDIF
      RETURN
*
   80 FORMAT(//52H RATIONAL APPROXIMATION COEFFICIENTS FOR FUEL-TO-FUE,
     1 59HL REDUCED COLLISION PROBABILITIES OF THE CLOSED CELL (GROUP,
     2 I5,2H):/1P,9X,10HNUMERATOR ,3(2H (,E11.4,1H,,E11.4,1H),:)/19X,
     3 3(2H (,E11.4,1H,,E11.4,1H),:))
   90 FORMAT(7X,12HDENOMINATOR ,3(2H (,E11.4,1H,,E11.4,1H),:)/19X,
     1 3(2H (,E11.4,1H,,E11.4,1H),:))
  100 FORMAT(/5X,1HI,9X,4HSIGX,8X,5HDILUT,10X,3HFIT)
  110 FORMAT(11X,8HINFINITE,1P,E13.5)
  120 FORMAT(//7H  SIGX:,1P,7E11.4)
  130 FORMAT(7H DILUT:,7E11.4/)
      END
