*DECK LIBPTT
      SUBROUTINE LIBPTT(IGRP,NDIL,NPART,DILUT,XSDIL,GOLD,HNAMIS,IMPX,
     1 NOR,WEIGH,SIGX,SIGP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transform dilution dependent total and partial self-shielded cross
* section into probability tables.
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
* IGRP    energy group index of the probability table.
* NDIL    number of finite dilutions.
* NPART   number of partial cross sections.
* DILUT   dilutions with DILUT(NDIL+1)=1.e10.
* XSDIL   dilution dependent self-shielded cross sections:
*         XSDIL(I,1) total self-shielded cross sections;
*         XSDIL(I,2) nu*fission self-shielded cross sections;
*         XSDIL(I,3) P0 scattering cross sections;
*         etc.
*         XSDIL(NDIL+1,j) are the infinite dilution values.
* GOLD    Goldstein-Cohen parameter.
* HNAMIS  local name of the isotope:
*         HNAMIS(1:8)  is the local isotope name;
*         HNAMIS(9:12) is a suffix function of the mixture index.
* IMPX    print parameter (equal to zero for no print).
*
*Parameters: output
* NOR     order for the probability table.
* WEIGH   quadrature weights fot the probability table.
* SIGX    base points for the total cross sections.
* SIGP    base points for the partial cross sections.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER (MAXNOR=12)
      INTEGER IGRP,NDIL,NPART,IMPX,NOR
      REAL DILUT(NDIL+1),XSDIL(NDIL+1,NPART+1),GOLD,WEIGH(NOR),
     1 SIGX(NOR),SIGP(MAXNOR,NPART)
      CHARACTER HNAMIS*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (TARGET=0.2D-3,EPSRID=5.0D-4)
      REAL PRECA
      DOUBLE PRECISION SIGXI,CC,DD,EE,DENOM
      DOUBLE PRECISION DA(0:MAXNOR-1),DB(0:MAXNOR-1),DC(0:MAXNOR)
      COMPLEX*16 SIGX0(MAXNOR),CCC,DCC,XCC
      LOGICAL LCONV,LFAIL
      CHARACTER HSMG*131
      REAL, ALLOCATABLE, DIMENSION(:) :: WABS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: TEST,SDDK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: TOFIT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(TEST(NPART+2),SDDK(NDIL),TOFIT(NDIL,2))
*
      IF(IMPX.GT.5) THEN
         WRITE(6,'(/47H LIBPTT: DILUTION DEPENDANT CROSS SECTIONS OF I,
     1   8HSOTOPE '',A12,10H'' IN GROUP,I4,1H:/9X,10HDILUTIONS=,1P,
     2   9E12.4,:/(19X,9E12.4))') HNAMIS,IGRP,(DILUT(I),I=1,NDIL+1)
         WRITE(6,'(10X,9HTOTAL XS=,1P,9E12.4,:/(19X,9E12.4))')
     1   (XSDIL(I,1),I=1,NDIL+1)
         IF(GOLD.NE.1.0) THEN
            WRITE(6,'(6X,13HPRINCIPAL XS=,1P,9E12.4,:/(19X,9E12.4))')
     1      (XSDIL(I,1)-(1.0-GOLD)*XSDIL(I,3),I=1,NDIL+1)
         ENDIF
         IF(IMPX.GT.6) THEN
            DO 10 J=1,NPART
            WRITE(6,'(4X,10HPARTIAL XS,I4,1H=,1P,9E12.4,:/(19X,
     1      9E12.4))') J,(XSDIL(I,J+1),I=1,NDIL+1)
   10       CONTINUE
         ENDIF
      ENDIF
      IF(NPART.LT.2) CALL XABORT('LIBPTT: SCATTERING INFO MISSING.')
*----
*  CHECK IF THE ENERGY GROUP IS REALLY RESONANT.
*----
      IF(NDIL.EQ.0) THEN
         NOR=1
         WEIGH(1)=1.0
         SIGX(1)=XSDIL(1,1)
         DO 20 J=1,NPART
         SIGP(1,J)=XSDIL(1,J+1)
   20    CONTINUE
         PREC0=0.0D0
         PREC=0.0D0
         PREC1=0.0D0
         GO TO 480
      ENDIF
      PREC=0.0D0
      PREC0=0.0D0
      PREC1=0.0D0
      LCONV=.FALSE.
      DO 30 IDIL=1,NDIL
      ERR=ABS(DBLE((XSDIL(NDIL+1,1))/(XSDIL(IDIL,1)))-1.0D0)
      PREC=MAX(PREC,ERR)
      ERR=ABS(DBLE((XSDIL(NDIL+1,1)-(1.0D0-GOLD)*XSDIL(NDIL+1,3))/
     1 (XSDIL(IDIL,1)-(1.0-GOLD)*XSDIL(IDIL,3)))-1.0D0)
      PREC0=MAX(PREC0,ERR)
      ERR=ABS(DBLE((XSDIL(NDIL+1,1)-XSDIL(NDIL+1,3))/
     1 (XSDIL(IDIL,1)-XSDIL(IDIL,3)))-1.0D0)
      PREC1=MAX(PREC1,ERR)
      LCONV=LCONV.OR.(ABS(XSDIL(IDIL,1)).EQ.ABS(XSDIL(IDIL+1,1)))
   30 CONTINUE
      IF(IMPX.GT.3) WRITE(6,'(/36H LIBPTT: ORDER  1 PROBABILITY TABLE ,
     1 24HCALCULATION OF ISOTOPE '',A12,10H'' IN GROUP,I4,8H. ERROR=,1P,
     2 3D11.3,1H.)') HNAMIS,IGRP,PREC0,PREC,PREC1
      IF(PREC.LE.TARGET) THEN
         NOR=1
         WEIGH(1)=1.0
         SIGX(1)=XSDIL(NDIL+1,1)
         DO 40 J=1,NPART
         SIGP(1,J)=XSDIL(NDIL+1,J+1)
   40    CONTINUE
         GO TO 360
      ENDIF
      IF(LCONV) THEN
         WRITE(HSMG,'(45HLIBPTT: UNIFORM TOTAL XS IS NOT EXPECTED IN G,
     1   4HROUP,I4,1H.)') IGRP
         CALL XABORT(HSMG)
      ENDIF
*----
*  FIND THE PADE APPROXIMATION FOR ABS+XGOLD*SIGS CROSS SECTION USING
*  A PADE REGRESSION.
*----
      XGOLD=GOLD
   45 ALLOCATE(WABS(NDIL+1))
      DO 60 IDIL=1,NDIL+1
      WABS(IDIL)=REAL(XSDIL(IDIL,1)-(1.0D0-XGOLD)*XSDIL(IDIL,3))
   60 CONTINUE
      CALL ALPLSF(3,NDIL+1,DILUT,WABS,EPSRID,.TRUE.,NOR,DA,DB,PRECA)
      DEALLOCATE(WABS)
      NOR=NOR+1
      IF(NOR.GT.MAXNOR) CALL XABORT('LIBPTT: NOR IS TOO LARGE.')
      PREC0=DBLE(PRECA)
*
*     FIND THE BASE POINTS IN ABS+XGOLD*SIGS CROSS SECTION.
      SGN=1.0D0
      DC(0)=DA(0)
      DO 70 I=2,NOR
      SGN=-SGN
      DC(I-1)=SGN*(DB(I-2)+DA(I-1))
   70 CONTINUE
      DC(NOR)=-SGN
      CALL ALROOT(DC,NOR,SIGX0,LFAIL)
      IF(LFAIL) CALL XABORT('LIBPTT: POLYNOMIAL ROOT FINDING FAILURE.')
*
      DO 110 I=1,NOR
*
*     NEWTON IMPROVEMENT OF THE ROOTS.
      CCC=0.0D0
      XCC=1.0D0
      DO 80 J=0,NOR
      CCC=CCC+DC(J)*XCC
      XCC=XCC*SIGX0(I)
   80 CONTINUE
      DCC=0.0D0
      XCC=1.0D0
      DO 85 J=1,NOR
      DCC=DCC+DC(J)*XCC*REAL(J)
      XCC=XCC*SIGX0(I)
   85 CONTINUE
      SIGX0(I)=SIGX0(I)-CCC/DCC
*
*     COMPUTE THE WEIGHTS.
      IF(AIMAG(CMPLX(SIGX0(I))).NE.0.0) CALL XABORT('LIBPTT: COMPLEX '
     1 //'ROOT.')
      SIGXI=DBLE(SIGX0(I))
      CC=1.0D0
      DD=0.0D0
      DO 90 J=0,NOR-1
      DD=DD+DB(J)*CC
      CC=-CC*SIGXI
   90 CONTINUE
      DO 100 J=1,NOR
      IF(J.NE.I) DD=DD/(DBLE(SIGX0(J))-SIGXI)
  100 CONTINUE
      WEIGH(I)=REAL(DD)
  110 CONTINUE
*----
*  PROCESS THE TOTAL CROSS SECTIONS.
*----
      DO 210 IDIL=1,NDIL
      SCC=DA(NOR-1)
      DO 200 I=NOR-2,0,-1
      SCC=DA(I)+SCC*DILUT(IDIL)
  200 CONTINUE
      SDDK(IDIL)=(XSDIL(IDIL,1)-(1.0-XGOLD)*XSDIL(IDIL,3))/SCC
      TOFIT(IDIL,1)=DILUT(IDIL)
      TOFIT(IDIL,2)=XSDIL(IDIL,1)/SDDK(IDIL)
      SDDK(IDIL)=SDDK(IDIL)*SDDK(IDIL)
  210 CONTINUE
      IF(XGOLD.NE.1.0) THEN
         CALL ALDFIT(NDIL,NOR-1,TOFIT(1,1),TOFIT(1,2),SDDK,DA)
      ENDIF
      DO 220 I=0,NOR-1
      DA(I)=DA(I)*XSDIL(NDIL+1,1)/DA(NOR-1)
  220 CONTINUE
*----
*  COMPUTE THE BASE POINTS IN TOTAL CROSS SECTION.
*----
      DO 240 I=1,NOR
      SIGXI=DBLE(SIGX0(I))
      CC=1.0D0
      DD=0.0D0
      EE=0.0D0
      DO 230 J=0,NOR-1
      DD=DD+DA(J)*CC
      EE=EE+DB(J)*CC
      CC=-CC*SIGXI
  230 CONTINUE
      SIGX(I)=REAL(DD/EE)
      IF(SIGX(I).LT.0.0) THEN
         IF(XGOLD.EQ.1.0) CALL XABORT('LIBPTT: NEGATIVE BASE POINTS FO'
     1   //'R THE TOTAL CROSS SECTION.')
         XGOLD=MIN(1.0D0,XGOLD+0.1D0)
         GO TO 45
      ENDIF
  240 CONTINUE
*----
*  PROCESS THE PARTIAL CROSS SECTIONS.
*----
      DO 300 IPART=1,NPART
      IF(XSDIL(NDIL+1,IPART+1).EQ.0.0) THEN
         DO 250 I=1,NOR
         SIGP(I,IPART)=0.0
  250    CONTINUE
         GO TO 300
      ENDIF
      DO 260 IDIL=1,NDIL
      TOFIT(IDIL,1)=DILUT(IDIL)
      TOFIT(IDIL,2)=XSDIL(IDIL,IPART+1)/SQRT(SDDK(IDIL))
  260 CONTINUE
      CALL ALDFIT(NDIL,NOR-1,TOFIT(1,1),TOFIT(1,2),SDDK,DA)
      IF(DA(NOR-1).EQ.0.0) THEN
         DO 265 I=1,NOR
         SIGP(I,IPART)=XSDIL(NDIL+1,IPART+1)
  265    CONTINUE
         GO TO 300
      ENDIF
      DO 270 I=0,NOR-1
      DA(I)=DA(I)*XSDIL(NDIL+1,IPART+1)/DA(NOR-1)
  270 CONTINUE
*----
*  COMPUTE THE BASE POINTS IN PARTIAL CROSS SECTION.
*----
      DO 290 I=1,NOR
      SIGXI=DBLE(SIGX0(I))
      CC=1.0D0
      DD=0.0D0
      EE=0.0D0
      DO 280 J=0,NOR-1
      DD=DD+DA(J)*CC
      EE=EE+DB(J)*CC
      CC=-CC*SIGXI
  280 CONTINUE
      SIGP(I,IPART)=REAL(DD/EE)
  290 CONTINUE
  300 CONTINUE
*----
*  REMOVING SMALL PROBABILITIES.
*----
      INOR=0
  330 INOR=INOR+1
      IF(INOR.GT.NOR) GO TO 360
      IF(ABS(WEIGH(INOR)).LE.5.0E-7) THEN
         DO 355 JNOR=INOR+1,NOR
         WEIGH(JNOR-1)=WEIGH(JNOR)
         SIGX(JNOR-1)=SIGX(JNOR)
         DO 350 J=1,NPART
         SIGP(JNOR-1,J)=SIGP(JNOR,J)
  350    CONTINUE
  355    CONTINUE
         INOR=INOR-1
         NOR=NOR-1
      ENDIF
      GO TO 330
*----
*  NORMALIZE THE PROBABILITY TABLE TO INFINITE DILUTION X-S.
*----
  360 CC=0.0D0
      DO 390 I=1,NOR
      CC=CC+WEIGH(I)
  390 CONTINUE
      DO 400 I=1,NOR
      WEIGH(I)=WEIGH(I)/REAL(CC)
  400 CONTINUE
      CC=0.0D0
      DO 410 I=1,NOR
      CC=CC+WEIGH(I)*SIGX(I)
  410 CONTINUE
      IF(CC.NE.0.0) THEN
         DO 420 I=1,NOR
         SIGX(I)=SIGX(I)*XSDIL(NDIL+1,1)/REAL(CC)
  420    CONTINUE
      ENDIF
      DO 450 J=1,NPART
      CC=0.0D0
      DO 430 I=1,NOR
      CC=CC+WEIGH(I)*SIGP(I,J)
  430 CONTINUE
      IF(CC.NE.0.0) THEN
         DO 440 I=1,NOR
         SIGP(I,J)=SIGP(I,J)*XSDIL(NDIL+1,J+1)/REAL(CC)
  440    CONTINUE
      ENDIF
  450 CONTINUE
*----
*  TEST THE ACCURACY OF THE PROBABILITY TABLE.
*----
      PREC=0.0D0
      PREC1=0.0D0
      DO 470 IDIL=1,NDIL+1
         CC=0.0D0
         DD=0.0D0
         EE=0.0D0
         DO 460 I=1,NOR
         DENOM=SIGX(I)-(1.0-GOLD)*SIGP(I,2)+DILUT(IDIL)
         CC=CC+WEIGH(I)/DENOM
         DD=DD+WEIGH(I)*SIGX(I)/DENOM
         EE=EE+WEIGH(I)*(SIGX(I)-SIGP(I,2))/DENOM
  460    CONTINUE
         PREC=MAX(PREC,ABS((DD/CC)/DBLE(XSDIL(IDIL,1))-1.0D0))
         PREC1=MAX(PREC1,ABS((EE/CC)/DBLE(XSDIL(IDIL,1)-
     1   XSDIL(IDIL,3))-1.0D0))
  470 CONTINUE
  480 IF(IMPX.GT.2) WRITE(6,'(14H LIBPTT: ORDER,I3,15H PROBABILITY TA,
     1 28HBLE CALCULATION OF ISOTOPE '',A12,10H'' IN GROUP,I4,5H. ERR,
     2 3HOR=,1P,3D11.3,1H.)') NOR,HNAMIS,IGRP,PREC0,PREC,PREC1
*
      IF(((IMPX.GT.4).AND.(NOR.GT.1)).OR.(IMPX.GT.5)) THEN
         WRITE(6,'(/27H LIBPTT: PROBABILITY TABLE:/7X,11HPROBABILITY,
     1   2X,10HTOTAL-----,2X,22HPARTIAL CROSS SECTIONS)')
         DO 490 J=1,NPART+2
         TEST(J)=0.0D0
  490    CONTINUE
         DO 510 INOR=1,NOR
         TEST(1)=TEST(1)+WEIGH(INOR)
         TEST(2)=TEST(2)+WEIGH(INOR)*SIGX(INOR)
         DO 500 J=1,NPART
         TEST(J+2)=TEST(J+2)+WEIGH(INOR)*SIGP(INOR,J)
  500    CONTINUE
         WRITE(6,'(1X,I5,1P,7E12.4,:/(30X,5E12.4))') INOR,
     1   WEIGH(INOR),SIGX(INOR),(SIGP(INOR,J),J=1,NPART)
  510    CONTINUE
         WRITE(6,'(6H CHECK,1P,7E12.4,:/(30X,5E12.4))') (REAL(TEST(J)),
     1   J=1,NPART+2)
         TEST(1)=1.0D0
         DO 520 J=1,NPART+1
         TEST(J+1)=XSDIL(NDIL+1,J)
  520    CONTINUE
         WRITE(6,'(6H EXACT,1P,7E12.4,:/(30X,5E12.4))') (REAL(TEST(J)),
     1   J=1,NPART+2)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(TOFIT,SDDK,TEST)
      RETURN
      END
