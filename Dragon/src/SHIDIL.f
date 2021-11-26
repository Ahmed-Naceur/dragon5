*DECK SHIDIL
      SUBROUTINE SHIDIL(NRAT,NALPHA,NBNRS,COEF,DENOM,DILUT,PICX,SIGX,
     1 DIST,VST,IMPX,LLL,XCOEF,XDENO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the zone-dependent weights and base points for a N-term
* rational approximation.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NRAT    number of terms in the pij rational approximation.
* NALPHA  number of available dilutions (NALPHA.ge.2*NRAT-1).
* NBNRS   number of totally correlated fuel regions.
* COEF    numerator for the fuel-to-fuel cp rational expansion.
* DENOM   base points for the fuel-to-fuel cp rational expansion.
* DILUT   average dilution.
* PICX    pic values.
* SIGX    resonant cross sections.
* DIST    number density ratio of the resonant isotope.
* VST     volumes of the resonant regions.
* IMPX    print flag (equal to zero for no print).
* LLL     energy group index.
*
*Parameters: output
* XCOEF   zone-dependent weights.
* XDENO   zone-dependent base points.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NRAT,NALPHA,NBNRS,IMPX,LLL
      REAL DILUT(NALPHA),SIGX(NALPHA),DIST(NBNRS),VST(NBNRS)
      DOUBLE PRECISION PICX(NALPHA,NBNRS)
      COMPLEX COEF(NRAT),DENOM(NRAT)
      COMPLEX*16 XCOEF(NRAT,NBNRS),XDENO(NRAT,NBNRS)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NMAX=11,NORIN=(NMAX-1)/2)
      DOUBLE PRECISION TOFIT(NMAX,2),SDDK(NMAX),SDDK2(NMAX),
     1 DA(0:NORIN),DB(0:NORIN),DC(0:NORIN),C(0:NORIN+1)
      COMPLEX*16 DD,E1,E2,AAA,SQRTM3,SIGXI,CDENOM(NORIN+1),DDC(0:NORIN)
      CHARACTER HSMG*131
      LOGICAL LFAIL
      PARAMETER (SQRTM3=(0.0,1.73205080756888))
*
      IF(NBNRS.EQ.1) THEN
         DO 10 I=1,NRAT
         XCOEF(I,1)=COEF(I)
         XDENO(I,1)=DENOM(I)
   10    CONTINUE
         RETURN
      ENDIF
      NORS=0
      DO 20 I=1,NRAT
      IF(COEF(I).NE.(0.0,0.0)) NORS=NORS+1
   20 CONTINUE
      IF(NORS.EQ.1) THEN
         DO 30 I=1,NBNRS
         XCOEF(1,I)=1.0D0
   30    CONTINUE
         GO TO 170
      ENDIF
      IF(NORS.GT.NORIN+1) CALL XABORT('SHIDIL: NORIN OVERFLOW.')
*----
*  TRANSFORM THE RATIONAL APPROXIMATION INTO A PADE REPRESENTATION
*----
      DO 40 I=0,NRAT-1
      DA(I)=0.0D0
      DB(I)=0.0D0
   40 CONTINUE
      DO 75 N=1,NRAT
      DDC(0)=(1.0D0,0.0D0)
      I0=0
      DO 60 I=1,NRAT
      IF((I.NE.N).AND.(COEF(I).NE.(0.0,0.0))) THEN
         I0=I0+1
         DDC(I0)=DDC(I0-1)
         DO 50 J=I0-1,1,-1
         DDC(J)=DDC(J-1)+DDC(J)*DENOM(I)
   50    CONTINUE
         DDC(0)=DDC(0)*DENOM(I)
      ENDIF
   60 CONTINUE
      DO 70 I=0,NRAT-1
      DA(I)=DA(I)+DBLE(COEF(N)*DENOM(N)*DDC(I))
      DB(I)=DB(I)+DBLE(COEF(N)*DDC(I))
   70 CONTINUE
   75 CONTINUE
      DO 80 I=0,NORS-1
      DA(I)=DA(I)/DB(NORS-1)
      DB(I)=DB(I)/DB(NORS-1)
   80 CONTINUE
*
      DO 100 IALP=1,2*NRAT-1
      GAR1=DA(NORS-1)
      DO 90 I=NORS-2,0,-1
      GAR1=DA(I)+GAR1*SIGX(IALP)
   90 CONTINUE
      SDDK(IALP)=DILUT(IALP)/GAR1
      SDDK2(IALP)=SDDK(IALP)*SDDK(IALP)
  100 CONTINUE
*----
*  PROCESS THE DISTRIBUTED DILUTIONS
*----
      DO 160 K=1,NBNRS
      DO 110 IALP=1,2*NRAT-1
      TOFIT(IALP,1)=SIGX(IALP)*DIST(K)
      DILUTM=1.0D0/PICX(IALP,K)-SIGX(IALP)*DIST(K)
      TOFIT(IALP,2)=DILUTM/SDDK(IALP)
  110 CONTINUE
      CALL ALDFIT(NALPHA,NORS-1,TOFIT(1,1),TOFIT(1,2),SDDK2,DC)
*
      QQ=DC(1)+DB(0)
      RR=DC(0)
      IF(NORS-1.EQ.0) THEN
*        1-TERM RATIONAL APPROXIMATION.
         XCOEF(1,K)=1.0D0
         XDENO(1,K)=DC(0)
      ELSE IF(NORS-1.EQ.1) THEN
*        2-TERMS RATIONAL APPROXIMATION.
         AAA=QQ*QQ-4.0D0*RR
         AAA=SQRT(AAA)
         E1=0.5D0*(QQ+AAA)
         E2=0.5D0*(QQ-AAA)
         IF(ABS(DBLE(E1*E2)-RR).GT.5.0E-3*ABS(RR)) THEN
            WRITE (HSMG,'(42HSHIDIL: INTERPOLATION ALGORITHM FAILURE 1,,
     1      6H COEF=,1P,3E11.3)') QQ,RR,DBLE(E1*E2)
            CALL XABORT(HSMG)
         ENDIF
*
         XCOEF(1,K)=(DB(0)-E1)/(E2-E1)
         XCOEF(2,K)=(DB(0)-E2)/(E1-E2)
         XDENO(1,K)=E1
         XDENO(2,K)=E2
      ELSE IF(NORS-1.GE.2) THEN
*        NORS-TERMS RATIONAL APPROXIMATION.
         SGN=1.0D0
         C(0)=DC(0)
         DO 120 I=2,NORS
         SGN=-SGN
         C(I-1)=SGN*(DB(I-2)/DIST(K)**(I-2)+DC(I-1))
  120    CONTINUE
         C(NORS)=-SGN
         CALL ALROOT(C,NORS,CDENOM,LFAIL)
         IF(LFAIL) CALL XABORT('SHIDIL: ROOT FINDING FAILURE.')
         DO 150 I=1,NORS
         SIGXI=CDENOM(I)
         XDENO(I,K)=CMPLX(SIGXI)
         DD=SIGXI**(NORS-1)
         SGN=1.0D0
         DO 130 J=NORS-1,1,-1
         SGN=-SGN
         DD=DD+SGN*DB(J-1)*SIGXI**(J-1)/DIST(K)**(J-1)
  130    CONTINUE
         DO 140 J=1,NORS
         IF(J.NE.I) DD=DD/(SIGXI-CDENOM(J))
  140    CONTINUE
         XCOEF(I,K)=CMPLX(DD)
  150    CONTINUE
      ELSE
         CALL XABORT('SHIDIL: PADE COLLOCATION FAILURE.')
      ENDIF
  160 CONTINUE
*
  170 DO 185 J=1,NBNRS
      DO 180 I=NORS+1,NRAT
      XCOEF(I,J)=(0.0,0.0)
  180 CONTINUE
  185 CONTINUE
      IF(IMPX.GE.10) THEN
         WRITE(6,'(/40H SHIDIL: ZONE-DEPENDENT WEIGHTS IN GROUP,I5)')
     1   LLL
         DO 190 I=1,NRAT
         WRITE(6,'(9H TERM NB.,I2,3X,1P,1H(,2E12.4,1H),:,2H (,2E12.4,
     1   1H),2H (,2E12.4,1H),:,2H (,2E12.4,1H),:/(14X,1H(,2E12.4,1H),
     2   :,2H (,2E12.4,1H),:,2H (,2E12.4,1H),:,2H (,2E12.4,1H)))')
     3   I,(XCOEF(I,J),J=1,NBNRS)
  190    CONTINUE
         WRITE(6,'(/36H SHIDIL: ZONE-DEPENDENT BASE POINTS:)')
         DO 200 I=1,NRAT
         WRITE(6,'(9H TERM NB.,I2,3X,1P,1H(,2E12.4,1H),:,2H (,2E12.4,
     1   1H),2H (,2E12.4,1H),:,2H (,2E12.4,1H),:/(14X,1H(,2E12.4,1H),
     2   :,2H (,2E12.4,1H),:,2H (,2E12.4,1H),:,2H (,2E12.4,1H)))')
     3   I,(XDENO(I,J),J=1,NBNRS)
  200    CONTINUE
      ENDIF
      IF(IMPX.GE.100) THEN
         DO 225 K=1,NBNRS
         WRITE(6,'(24H SHIDIL: RESONANT REGION,I4,1H:/14X,3HPIC,8X,
     1   3HFIT)') K
         DO 220 IALP=1,NALPHA
         E1=0.0
         DO 210 I=1,NRAT
         DD=SIGX(IALP)*DIST(K)+XDENO(I,K)
         E1=E1+XCOEF(I,K)/DD
  210    CONTINUE
         WRITE(6,'(3X,I3,1P,2E11.3,F10.2,1H%)') IALP,PICX(IALP,K),
     1   DBLE(E1),100.0*(DBLE(E1)-PICX(IALP,K))/ABS(PICX(IALP,K))
  220    CONTINUE
  225    CONTINUE
         WRITE(6,'(22H SHIDIL: OVERALL FUEL:/14X,3HPXX,8X,3HFIT)')
         GAR1=0.0
         DO 230 K=1,NBNRS
         GAR1=GAR1+VST(K)
  230    CONTINUE
         DO 260 IALP=1,NALPHA
         GAR2=0.0
         E1=0.0
         DO 250 K=1,NBNRS
         DO 240 I=1,NRAT
         DD=SIGX(IALP)*DIST(K)+XDENO(I,K)
         E1=E1+XCOEF(I,K)*VST(K)/(GAR1*DD)
  240    CONTINUE
         GAR2=GAR2+PICX(IALP,K)*VST(K)/GAR1
  250    CONTINUE
         WRITE(6,'(3X,I3,1P,2E11.3,F10.2,1H%)') IALP,GAR2,DBLE(E1),
     1   100.0*(DBLE(E1)-GAR2)/ABS(GAR2)
  260    CONTINUE
      ENDIF
      RETURN
      END
