*DECK CHAB03
      SUBROUTINE CHAB03(IPLIB,IMPX,NGRP,NBIN,IMOD,TYPSEC,HISOT,VALUE,
     1 IGM,IGP,NFS,VAL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Modify a specific isotope and Autolib reaction in a Draglib.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   LCM pointer to the Draglib.
* IMPX    print index.
* NGRP    number of coarse energy groups.
* NBIN    number of fine energy groups.
* IMOD    type of modification: =1: complete replacement; =2: replace
*         specific values by VALUE; =3: increase by VALUE; =4: multiply
*         by VALUE.
* TYPSEC  name of reaction to modify.
* HISOT   name of isotope to modify.
* VALUE   value used in modification operation.
* IGM     first energy group to modify.
* IGP     last energy group to modify.
* NFS     number of fine groups per coarse group.
* VAL     array of values used if IMOD=1.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IMPX,NGRP,NBIN,IMOD,IGM,IGP,NFS(NGRP)
      CHARACTER TYPSEC*8,HISOT*12
      REAL VALUE,VAL(NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6)
      CHARACTER AJUS(4)*4
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,DELTA,FMULT,VALB
*----
*  DATA STATEMENTS
*----
      DATA AJUS/'VALE','CONS','PLUS','MULT'/
*----
*   CORRESPONDENCE BETWEEN BIN AND COARSE ENERGT GROUPS
*----
      IGMBIN=NBIN+1
      IGPBIN=0
      IBIN=0
      DO 10 IG=1,NGRP
      IF(IG.EQ.IGM) IGMBIN=IBIN+1
      IBIN=IBIN+NFS(IG)
      IF(IG.EQ.IGP) IGPBIN=IBIN
   10 CONTINUE
      IF(IGPBIN.LT.IGMBIN) RETURN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAR1(NBIN),DELTA(NBIN),FMULT(NBIN),VALB(NBIN))
*
      IF(IMOD.EQ.1) THEN
         IBIN=0
         DO 25 IG=1,NGRP
         DO 20 J=1,NFS(IG)
         IBIN=IBIN+1
         VALB(IBIN)=VAL(IG)
   20    CONTINUE
   25    CONTINUE
      ENDIF
*----
*   APPLY CORRECTION
*----
      IF(TYPSEC.EQ.'NTOT0') THEN
         CALL LCMLEN(IPLIB,'BIN-NTOT0',ILONG,ITYLCM)
         IF(ILONG.EQ.NBIN) THEN
            IF(IMPX.GT.0) WRITE(IOUT,'(/17H CHAB03: MODIFY (,A,5H) BIN,
     1      27H-NTOT0 REACTION OF ISOTOPE ,A,1H.)') AJUS(IMOD),HISOT
            CALL LCMGET(IPLIB,'BIN-NTOT0',GAR1)
            CALL CHAB02(NBIN,IMOD,VALUE,IGMBIN,IGPBIN,VALB,GAR1,DELTA,
     1      FMULT)
            CALL LCMPUT(IPLIB,'BIN-NTOT0',NBIN,2,GAR1)
         ENDIF
      ELSE IF((TYPSEC(:4).EQ.'SIGS').OR.(TYPSEC(:4).EQ.'SCAT')) THEN
         CALL LCMLEN(IPLIB,'BIN-SIGS00',ILONG,ITYLCM)
         IF(ILONG.EQ.NBIN) THEN
            IF(IMPX.GT.0) WRITE(IOUT,'(/17H CHAB03: MODIFY (,A,5H) BIN,
     1      28H-SIGS00 REACTION OF ISOTOPE ,A,1H.)') AJUS(IMOD),HISOT
            CALL LCMGET(IPLIB,'BIN-SIGS00',GAR1)
            CALL CHAB02(NBIN,IMOD,VALUE,IGMBIN,IGPBIN,VALB,GAR1,DELTA,
     1      FMULT)
            CALL LCMPUT(IPLIB,'BIN-SIGS00',NBIN,2,GAR1)
            CALL LCMGET(IPLIB,'BIN-NTOT0',GAR1)
            DO 30 IBIN=1,NBIN
            GAR1(IBIN)=GAR1(IBIN)+DELTA(IBIN)
   30       CONTINUE
            CALL LCMPUT(IPLIB,'BIN-NTOT0',NBIN,2,GAR1)
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VALB,FMULT,DELTA,GAR1)
      RETURN
      END
