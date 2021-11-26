*DECK CHAB01
      SUBROUTINE CHAB01(IPLIB,IMPX,IRHS,NGRP,NLEG,IMOD,TYPSEC,HISOT,
     1 VALUE,IGM,IGP,VAL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modify a specific isotope and reaction in a microlib.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   LCM pointer to the Microlib or Draglib.
* IMPX    print index.
* IRHS    type of IPLIB: =1: Microlib; =2: Draglib.
* NGRP    number of energy groups.
* NLEG    max Legendre order of scattering anisotropy (1=isotropic,
*         etc.).
* IMOD    type of modification: =1: complete replacement; =2: replace
*         specific values by VALUE; =3: increase by VALUE; =4: multiply
*         by VALUE.
* TYPSEC  name of reaction to modify.
* HISOT   name of isotope to modify.
* VALUE   value used in modification operation.
* IGM     first energy group to modify.
* IGP     last energy group to modify.
* VAL     array of values used if IMOD=1.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IMPX,IRHS,NGRP,NLEG,IMOD,IGM,IGP
      CHARACTER TYPSEC*8,HISOT*12
      REAL VALUE,VAL(NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6,NCAPT=5)
      CHARACTER AJUS(4)*4,HCAPT(NCAPT)*8,CM*2
      REAL, ALLOCATABLE, DIMENSION(:) :: XSECT,DELTA,FMULT,GAR1
*----
*  DATA STATEMENTS
*----
      DATA AJUS/'VALE','CONS','PLUS','MULT'/
      DATA HCAPT/'NG','NP','NA','ND','NT'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XSECT(NGRP),DELTA(NGRP),FMULT(NGRP),GAR1(NGRP))
*
      IF(IMPX.GT.0) WRITE(IOUT,'(/17H CHAB01: MODIFY (,A,11H) REACTION ,
     1 A,12H OF ISOTOPE ,A,1H.)') AJUS(IMOD),TYPSEC,HISOT
      CALL LCMLEN(IPLIB,TYPSEC,ILONG,ITYLCM)
      IF((ILONG.EQ.0.).AND.(TYPSEC(:4).NE.'CAPT')
     1 .AND.(TYPSEC(:2).NE.'NU')) THEN
         CALL XABORT('CHAB01: MISSING REACTION '//TYPSEC//'.')
      ENDIF
*----
*  MODIFY CROSS SECTION
*----
      CALL XDRSET(XSECT,NGRP,0.0)
      CALL XDRSET(GAR1,NGRP,0.0)
      IF(TYPSEC.EQ.'NTOT0') THEN
         CALL LCMGET(IPLIB,TYPSEC,XSECT)
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,XSECT,DELTA,FMULT)
         CALL LCMPUT(IPLIB,TYPSEC,NGRP,2,XSECT)
      ELSE IF(TYPSEC.EQ.'NG'.OR.TYPSEC.EQ.'NP'.OR.TYPSEC.EQ.'NA'.OR.
     1 TYPSEC.EQ.'ND'.OR.TYPSEC.EQ.'NT') THEN

*        application of the perturbation

         CALL LCMGET(IPLIB,TYPSEC,XSECT)
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,XSECT,DELTA,FMULT)
         CALL LCMPUT(IPLIB,TYPSEC,NGRP,2,XSECT)
         CALL XDRSET(XSECT,NGRP,0.0)
         CALL LCMGET(IPLIB,'NTOT0',XSECT)
         DO 10 IG1=1,NGRP
         XSECT(IG1)=XSECT(IG1)+DELTA(IG1)
   10    CONTINUE
         CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XSECT)
      ELSE IF(TYPSEC.EQ.'CAPT') THEN
         IF(IMOD.NE.4) CALL XABORT('CHAB01: ONLY MULT ALLOWED.')
         DO 320 ICAPT=1,NCAPT
         TYPSEC=HCAPT(ICAPT)
         CALL LCMLEN(IPLIB,TYPSEC,ILONG,ITYLCM)
         IF(ILONG.NE.0.0) THEN
*           application of the perturbation
            WRITE(IOUT,*) 'CHAB01: REACTION CAPTURE INCLUDES ',TYPSEC
            CALL XDRSET(XSECT,NGRP,0.0)
            CALL LCMGET(IPLIB,TYPSEC,XSECT)
            CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,XSECT,DELTA,FMULT)
            CALL LCMPUT(IPLIB,TYPSEC,NGRP,2,XSECT)
            CALL XDRSET(XSECT,NGRP,0.0)
            CALL LCMGET(IPLIB,'NTOT0',XSECT)
            DO 310 IG1=1,NGRP
            XSECT(IG1)=XSECT(IG1)+DELTA(IG1)
  310       CONTINUE
            CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XSECT)
         ENDIF
  320    CONTINUE
         TYPSEC='CAPT'
      ELSE IF(TYPSEC.EQ.'NELAS'.OR.TYPSEC.EQ.'NINEL') THEN
         CALL LCMGET(IPLIB,TYPSEC,XSECT)
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,XSECT,DELTA,FMULT)
         CALL LCMPUT(IPLIB,TYPSEC,NGRP,2,XSECT)
*
*        additive modification of P0 scattering information
         JMOD=3
         CALL CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,JMOD,0,IGM,IGP,DELTA,
     1   DELTA,FMULT)
*
*        multiplicative modification of transport correction
         CALL LCMLEN(IPLIB,'TRANC',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL XDRSET(XSECT,NGRP,0.0)
            CALL LCMGET(IPLIB,'TRANC',XSECT)
            DO 20 IG1=1,NGRP
            XSECT(IG1)=XSECT(IG1)*FMULT(IG1)
   20       CONTINUE
            CALL LCMPUT(IPLIB,'TRANC',NGRP,2,XSECT)
         ENDIF
*
*        multiplicative modification of Pn scattering information with
*        n>0
         DO 30 JL=1,NLEG-1
            WRITE(CM,'(I2.2)') JL
            CALL LCMLEN(IPLIB,'SCAT'//CM,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               JMOD=4
               CALL CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,JMOD,JL,IGM,IGP,
     1         FMULT,DELTA,FMULT)
            ENDIF
   30    CONTINUE
*
*        additive modification of total cross section
         CALL XDRSET(XSECT,NGRP,0.0)
         CALL LCMGET(IPLIB,'NTOT0',XSECT)
         DO 40 IG1=1,NGRP
         XSECT(IG1)=XSECT(IG1)+DELTA(IG1)
   40    CONTINUE
         CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XSECT)
      ELSE IF((TYPSEC.EQ.'N2N').OR.(TYPSEC.EQ.'N3N').OR.
     1        (TYPSEC.EQ.'N4N')) THEN

         CALL LCMGET(IPLIB,TYPSEC,GAR1)
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,GAR1,DELTA,FMULT)
         CALL LCMPUT(IPLIB,TYPSEC,NGRP,2,GAR1)
*
*        additive modification of total cross section
         CALL LCMGET(IPLIB,'NTOT0',XSECT)
         DO 50 IG1=1,NGRP
         XSECT(IG1)=XSECT(IG1)+DELTA(IG1)
   50    CONTINUE
         CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XSECT)
*
*        additive modification of P0 scattering information
         IF (TYPSEC.EQ.'N2N') THEN
            DO 60 IG1=1,NGRP
            DELTA(IG1)=2.0*DELTA(IG1)
   60      CONTINUE
         ELSE IF (TYPSEC.EQ.'N3N') THEN
            DO 70 IG1=1,NGRP
            DELTA(IG1)=3.0*DELTA(IG1)
   70       CONTINUE
         ELSE IF (TYPSEC.EQ.'N4N') THEN
            DO 80 IG1=1,NGRP
            DELTA(IG1)=4.0*DELTA(IG1)
   80       CONTINUE
         ENDIF
         JMOD=3
         CALL CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,JMOD,0,IGM,IGP,DELTA,
     1   DELTA,FMULT)
*
*        multiplicative modification of Pn scattering information with
*        n>0
         DO 90 JL=1,NLEG-1
            WRITE(CM,'(I2.2)') JL
            CALL LCMLEN(IPLIB,'SCAT'//CM,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               JMOD=4
               CALL CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,JMOD,JL,IGM,IGP,
     1         FMULT,DELTA,FMULT)
            ENDIF
   90    CONTINUE
      ELSE IF((TYPSEC(:4).EQ.'SIGS').OR.(TYPSEC(:4).EQ.'SCAT')) THEN
         READ(TYPSEC(5:6),'(I2)') IL
*        additive or multiplicative modification of Pn scattering
*        information
         CALL XDRSET(XSECT,NGRP,0.0)
         IF(IMOD.EQ.1) THEN
            DO 100 IG=IGM,IGP
            XSECT(IG)=VAL(IG)
  100       CONTINUE
         ELSE
            DO 110 IG=IGM,IGP
            XSECT(IG)=VALUE
  110       CONTINUE
         ENDIF
         CALL CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,IMOD,IL,IGM,IGP,XSECT,
     1   DELTA,FMULT)
*
*        multiplicative modification of transport correction
         CALL LCMLEN(IPLIB,'TRANC',ILONG,ITYLCM)
         IF((IL.LE.1).AND.(ILONG.GT.0)) THEN
            CALL XDRSET(XSECT,NGRP,0.0)
            CALL LCMGET(IPLIB,'TRANC',XSECT)
            DO 120 IG1=1,NGRP
            XSECT(IG1)=XSECT(IG1)*FMULT(IG1)
  120       CONTINUE
            CALL LCMPUT(IPLIB,'TRANC',NGRP,2,XSECT)
         ENDIF
*
*        additive modification of total cross-section
         IF(IL.EQ.0) THEN
            CALL XDRSET(XSECT,NGRP,0.0)
            CALL LCMGET(IPLIB,'NTOT0',XSECT)
            DO 130 IG1=1,NGRP
            XSECT(IG1)=XSECT(IG1)+DELTA(IG1)
  130    CONTINUE
            CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XSECT)
         ENDIF
*
*        multiplicative modification of Pn scattering information with
*        n>IL
         DO 140 JL=IL+1,NLEG-1
            WRITE(CM,'(I2.2)') JL
            CALL LCMLEN(IPLIB,'SCAT'//CM,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               JMOD=4
               CALL CHAB04(IPLIB,IMPX,IRHS,NGRP,NLEG,JMOD,JL,IGM,IGP,
     1         FMULT,DELTA,FMULT)
            ENDIF
  140    CONTINUE
      ELSE IF((TYPSEC.EQ.'NFTOT').OR.(TYPSEC.EQ.'NUSIGF')) THEN
         CALL LCMGET(IPLIB,'NFTOT',GAR1)
         CALL LCMGET(IPLIB,'NUSIGF',XSECT)
         DO 180 IG1=1,NGRP
         IF(GAR1(IG1).NE.0.0) THEN
            XSECT(IG1)=XSECT(IG1)/GAR1(IG1)
         ENDIF
  180    CONTINUE
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,GAR1,DELTA,FMULT)
         DO 190 IG1=1,NGRP
         XSECT(IG1)=GAR1(IG1)*XSECT(IG1)
  190    CONTINUE
         CALL LCMPUT(IPLIB,'NFTOT',NGRP,2,GAR1)
         CALL LCMPUT(IPLIB,'NUSIGF',NGRP,2,XSECT)
         CALL XDRSET(XSECT,NGRP,0.0)
         CALL LCMGET(IPLIB,'NTOT0',XSECT)
         DO 200 IG1=1,NGRP
         XSECT(IG1)=XSECT(IG1)+DELTA(IG1)
  200    CONTINUE
         CALL LCMPUT(IPLIB,'NTOT0',NGRP,2,XSECT)
      ELSE IF(TYPSEC.EQ.'NU') THEN
         CALL LCMGET(IPLIB,'NFTOT',GAR1)
         CALL LCMGET(IPLIB,'NUSIGF',XSECT)
         DO 210 IG1=1,NGRP
         IF(GAR1(IG1).NE.0.0) THEN
            XSECT(IG1)=XSECT(IG1)/GAR1(IG1)
         ENDIF
  210    CONTINUE
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,XSECT,DELTA,FMULT)
         DO 220 IG1=1,NGRP
         XSECT(IG1)=GAR1(IG1)*XSECT(IG1)
  220    CONTINUE
         CALL LCMPUT(IPLIB,'NUSIGF',NGRP,2,XSECT)
      ELSE IF(TYPSEC.EQ.'CHI') THEN
         CALL LCMGET(IPLIB,TYPSEC,GAR1)
         CALL CHAB02(NGRP,IMOD,VALUE,IGM,IGP,VAL,GAR1,DELTA,FMULT)
         SUM=0.0
         DO 230 IG1=1,NGRP
         SUM=SUM+GAR1(IG1)
  230    CONTINUE
         DO 240 IG1=1,NGRP
         GAR1(IG1)=GAR1(IG1)/SUM
  240    CONTINUE
         CALL LCMPUT(IPLIB,TYPSEC,NGRP,2,GAR1)
      ELSE
         CALL XABORT('CHAB01: UNKNOWN REACTION '//TYPSEC//'.')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR1,FMULT,DELTA,XSECT)
      RETURN
      END
