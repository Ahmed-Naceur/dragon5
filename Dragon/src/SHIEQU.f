*DECK SHIEQU
      SUBROUTINE SHIEQU(IPLIB,LEVEL,NGRO,NBISO,NBM,NBNRS,NRAT,MIX,
     1 ISONAM,NOCONV,ISONR,GC,COEF,DENOM,XCOEF,XDENO,DEN,IMPX,SN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the equivalent dilution.
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
* IPLIB   pointer to the internal microscopic cross section library
*         (L_LIBRARY signature).
* LEVEL   type of approximation (=1 use original Stamm'ler and
*         Nordheim approximations; =2 use Riemann integration method
*         with Nordheim approximation).
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBM     number of mixtures in the macrolib.
* NBNRS   number of totally correlated resonant regions.
* NRAT    number of terms in the rational pij expansion.
* MIX     mix number of each isotope.
* ISONAM  alias name of isotopes.
* NOCONV  mixture convergence flag (.TRUE. if mixture IBM
*         is not converged in group L).
* ISONR   resonant isotope indices.
* GC      Goldstein-Cohen parameters.
* COEF    zone-independent weights.
* DENOM   zone-independent base points.
* XCOEF   zone-dependent weights.
* XDENO   zone-dependent base points.
* DEN     isotopic number density.
* IMPX    print flag (equal to zero for no print).
*
*Parameters: output
* SN      equivalent dilution.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      LOGICAL NOCONV(NBM,NGRO)
      INTEGER LEVEL,NGRO,NBISO,NBM,NBNRS,NRAT,MIX(NBISO),
     1 ISONAM(3,NBISO),ISONR(NBNRS),IMPX
      REAL GC(NGRO,NBNRS),DEN(NBISO),SN(NGRO,NBISO)
      COMPLEX COEF(NRAT,NGRO),DENOM(NRAT,NGRO)
      COMPLEX*16 XCOEF(NRAT,NBNRS,NGRO),XDENO(NRAT,NBNRS,NGRO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXRAT=9)
      TYPE(C_PTR) JPLIB,KPLIB
      LOGICAL LD
      CHARACTER HNAMIS*12,HSMG*131
      COMPLEX*16 EAV,TTT,CBS(MAXRAT)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NFS
      REAL, ALLOCATABLE, DIMENSION(:) :: EBIN,TBIN,SBIN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NFS(NGRO))
*
      IF(NRAT.GT.MAXRAT) CALL XABORT('SHIEQU: MAXRAT OVERFLOW.')
      JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
      DO 200 IRS=1,NBNRS
      ISO=ISONR(IRS)
      KPLIB=LCMGIL(JPLIB,ISO) ! set ISO-th isotope
      WRITE(HNAMIS,'(3A4)') (ISONAM(I0,ISO),I0=1,3)
      CALL LCMLEN(KPLIB,'BIN-NFS',LENBIN,ITYLCM)
      IF((LENBIN.GT.0).AND.(LEVEL.EQ.2)) THEN
         CALL LCMGET(KPLIB,'BIN-NFS',NFS)
         LBIN=0
         DO 10 IGRP=1,NGRO
         LBIN=LBIN+NFS(IGRP)
   10    CONTINUE
         ALLOCATE(EBIN(LBIN+1),TBIN(LBIN),SBIN(LBIN))
         CALL LCMGET(KPLIB,'BIN-ENERGY',EBIN)
         CALL LCMGET(KPLIB,'BIN-NTOT0',TBIN)
         CALL LCMGET(KPLIB,'BIN-SIGS00',SBIN)
      ELSE
         CALL XDISET(NFS,NGRO,0)
      ENDIF
      LBIN=0
      DO 110 IGRP=1,NGRO
      IF(NOCONV(MIX(ISO),IGRP)) THEN
         EAV=0.0
         IF(NFS(IGRP).GT.0) THEN
*           RIEMANN INTEGRATION METHOD WITH NEWTON-RAPHSON ACCELERATION.
            IF(IMPX.GE.10) THEN
               WRITE(6,'(/17H SHIEQU: WEIGHTS:)')
               WRITE(6,290) IGRP,(COEF(I,IGRP),I=1,NRAT)
               WRITE(6,'(/21H SHIEQU: BASE POINTS:)')
               WRITE(6,290) IGRP,(DENOM(I,IGRP),I=1,NRAT)
            ENDIF
            DO 20 I=1,NRAT
            EAV=EAV+COEF(I,IGRP)*SQRT(DENOM(I,IGRP))
   20       CONTINUE
            SNI=(DBLE(EAV)**2)/DEN(ISO)
            UG=0.0
            BA=0.0
            BS=0.0
            DO 30 IGF=1,NFS(IGRP)
            DELM=LOG(EBIN(LBIN+IGF)/EBIN(LBIN+IGF+1))
            UG=UG+DELM
            SIGA=MAX(0.0,TBIN(LBIN+IGF)-SBIN(LBIN+IGF))
            SIGS=GC(IGRP,IRS)*MAX(0.0,SBIN(LBIN+IGF))
            BA=BA+DELM*SIGA/(SIGA+SIGS+SNI)
            BS=BS+DELM*SIGS/(SIGA+SIGS+SNI)
   30       CONTINUE
            BA=BA/UG
            BS=BS/UG
            ZCAL=SNI*BA/(1.0D0-BS)
            DO 50 I=1,NRAT
            CBS(I)=0.0
            IF(XCOEF(I,IRS,IGRP).EQ.0.0) GO TO 50
            EAV=XDENO(I,IRS,IGRP)/DEN(ISO)
            DO 40 IGF=1,NFS(IGRP)
            DELM=LOG(EBIN(LBIN+IGF)/EBIN(LBIN+IGF+1))
            SIGA=MAX(0.0,TBIN(LBIN+IGF)-SBIN(LBIN+IGF))
            SIGS=GC(IGRP,IRS)*MAX(0.0,SBIN(LBIN+IGF))
            CBS(I)=CBS(I)+DELM*SIGS/(SIGA+SIGS+EAV)
   40       CONTINUE
            CBS(I)=CBS(I)/UG
            CBS(I)=1.0D0/(1.0D0-CBS(I))
   50       CONTINUE
            TAUXA=0.0D0
            DO 60 IGF=1,NFS(IGRP)
            DELM=LOG(EBIN(LBIN+IGF)/EBIN(LBIN+IGF+1))
            SIGA=MAX(0.0,TBIN(LBIN+IGF)-SBIN(LBIN+IGF))
            SIGS=GC(IGRP,IRS)*MAX(0.0,SBIN(LBIN+IGF))
            TTT=0.0D0
            DO 55 I=1,NRAT
            IF(XCOEF(I,IRS,IGRP).EQ.0.0) GO TO 55
            EAV=XDENO(I,IRS,IGRP)/DEN(ISO)
            TTT=TTT+XCOEF(I,IRS,IGRP)*EAV*CBS(I)/(SIGA+SIGS+EAV)
   55       CONTINUE
            TAUXA=TAUXA+DELM*SIGA*MAX(0.0D0,DBLE(TTT))
   60       CONTINUE
            IF(TAUXA.EQ.0.0) THEN
               SNI=1.0E10
               GO TO 90
            ENDIF
            TAUXA=TAUXA/UG
            ITER=0
   70       ITER=ITER+1
            IF(IMPX.GE.10) THEN
               WRITE(6,'(15H SHIEQU: GROUP=,I4,11H ITERATION=,I3,
     1         10H DILUTION=,1P,E11.4,16H REFERENCE RATE=,E11.4,
     2         13H APPROXIMATE=,E11.4)') IGRP,ITER,SNI,TAUXA,ZCAL
            ENDIF
            IF(ABS(TAUXA-ZCAL).LE.1.0D-5*ABS(TAUXA)) GO TO 90
            IF(ITER.GT.20) THEN
               WRITE(6,'(15H SHIEQU: GROUP=,I4,10H DILUTION=,1P,E11.4,
     1         16H REFERENCE RATE=,E11.4,13H APPROXIMATE=,E11.4,
     2         9H ISOTOPE=,A12,1H.)') IGRP,SNI,TAUXA,ZCAL,HNAMIS
               IF(ABS(TAUXA-ZCAL).LE.5.0E-2*ABS(TAUXA)) THEN
                  WRITE(6,'(24H SHIEQU: *** WARNING ***)')
                  GO TO 90
               ENDIF
               CALL XABORT('SHIEQU: CONVERGENCE FAILURE.')
            ENDIF
            BA=0.0
            BS=0.0
            DBA=0.0
            DBS=0.0
            DO 80 IGF=1,NFS(IGRP)
            DELM=LOG(EBIN(LBIN+IGF)/EBIN(LBIN+IGF+1))
            SIGA=MAX(0.0,TBIN(LBIN+IGF)-SBIN(LBIN+IGF))
            SIGS=GC(IGRP,IRS)*MAX(0.0,SBIN(LBIN+IGF))
            BA=BA+DELM*SIGA/(SIGA+SIGS+SNI)
            BS=BS+DELM*SIGS/(SIGA+SIGS+SNI)
            DBA=DBA-DELM*SIGA/(SIGA+SIGS+SNI)**2
            DBS=DBS-DELM*SIGS/(SIGA+SIGS+SNI)**2
   80       CONTINUE
            BA=BA/UG
            BS=BS/UG
            DBA=DBA/UG
            DBS=DBS/UG
            ZCAL=SNI*BA/(1.0D0-BS)
            DZCAL=BA/(1.0D0-BS)+SNI*DBA/(1.0D0-BS)+SNI*BA*DBS/
     1      (1.0D0-BS)**2
            IF(DZCAL.LT.1.0D-15*ZCAL) GO TO 90
            SNI=MAX(1.0D0,SNI-REAL((ZCAL-TAUXA)/DZCAL))
            IF(SNI.GE.1.0E10) THEN
               SNI=1.0E10
               GO TO 90
            ENDIF
            GO TO 70
   90       IF(IMPX.GE.5) THEN
               WRITE(6,'(16H SHIEQU: REGION=,I3,7H GROUP=,I5,7H ISOTOP,
     1         3HE='',A12,19H'' NB.OF ITERATIONS=,I3,1H.)') IRS,IGRP,
     2         HNAMIS,ITER
            ENDIF
         ELSE
*           ORIGINAL STAMM'LER APPROXIMATION.
            DO 100 I=1,NRAT
            EAV=EAV+XCOEF(I,IRS,IGRP)*SQRT(XDENO(I,IRS,IGRP))
  100       CONTINUE
            SNI=MAX(1.0D0,(DBLE(EAV)**2)/DEN(ISO))
         ENDIF
         SN(IGRP,ISO)=REAL(SNI)
         IF(SN(IGRP,ISO).LE.0.0) THEN
            WRITE (HSMG,300) HNAMIS,SN(IGRP,ISO),IGRP
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
      LBIN=LBIN+NFS(IGRP)
  110 CONTINUE
      IF((LENBIN.GT.0).AND.(LEVEL.EQ.2)) THEN
         DEALLOCATE(SBIN,TBIN,EBIN)
      ENDIF
      IF(IMPX.GE.5) THEN
         LD=.FALSE.
         DO 120 IGRP=1,NGRO
         LD=LD.OR.NOCONV(MIX(ISO),IGRP)
  120    CONTINUE
         IF(LD) WRITE(6,310) HNAMIS,(SN(IGRP,ISO),IGRP=1,NGRO)
      ENDIF
  200 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NFS)
      RETURN
*
  290 FORMAT(10H GROUP NB.,I4,3X,1P,1H(,2E12.4,1H),:,2H (,2E12.4,1H),
     1 2H (,2E12.4,1H),:,2H (,2E12.4,1H),:/(15X,1H(,2E12.4,1H),:,2H (,
     2 2E12.4,1H),:,2H (,2E12.4,1H),:,2H (,2E12.4,1H)))
  300 FORMAT(30HSHIEQU: THE RESONANT ISOTOPE ',A12,16H' HAS A NEGATIVE,
     1 25H DILUTION CROSS-SECTION (,1P,E14.4,0P,10H) IN GROUP,I4,1H.)
  310 FORMAT(/31H SHIEQU: DILUTIONS OF ISOTOPE ',A12,2H':/(1P,10E12.4))
      END
