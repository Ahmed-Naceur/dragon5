*DECK OUTAUX
      SUBROUTINE OUTAUX (IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,NEL,NUN,
     1 NALBP,NZS,NGCOND,MAT,VOL,IDL,EVECT,IHOM,IGCOND,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform homogenization into NZS regions and condensation into NGCOND
* macrogroups based on averaged fluxes contained in EVECT. Create an
* output extended macrolib containing homogenized volumes, integrated
* fluxes and cross sections.
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
* IPMAC1  L_MACROLIB pointer to the input macrolib.
* IPMAC2  L_MACROLIB pointer to the output extended macrolib.
* NBMIX   number of material mixtures.
* NL      scattering anisotropy.
* NBFIS   number of fissionable isotopes.
* NGRP    total number of energy groups.
* NEL     number of finite elements.
* NUN     number of unknowns per energy group.
* NALBP   number of physical albedos.
* NZS     number of homogenized regions so that NZS=max(IHOM(i)).
* NGCOND  number of macrogroups after energy condensation.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* IDL     position of the average flux component associated with
*         each volume.
* EVECT   unknowns.
* IHOM    homogenized index assigned to each element.
* IGCOND  limit of condensed groups.
* IMPX    print parameter (equal to zero for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1,IPMAC2
      PARAMETER(NREAC=11)
      INTEGER NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,NZS,NGCOND,MAT(NEL),
     1 IDL(NEL),IHOM(NEL),IGCOND(NGCOND),IMPX
      REAL VOL(NEL),EVECT(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,KPMAC1,JPMAC2,KPMAC2
      PARAMETER(NSTATE=40)
      CHARACTER HREAC(NREAC)*12,TEXT12*12,SUFF*2,TEXT6*6
      INTEGER IDATA(NSTATE)
      LOGICAL LNUSIG,LESTOP,LFIXE,LREAC(NREAC)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: IJJ,NJJ,IPOS
      REAL, DIMENSION(:), ALLOCATABLE :: VOLI,WORK,SCAT,RATE,GAR,DEN,
     1 DEN2
      REAL, DIMENSION(:,:), ALLOCATABLE :: FLINT,CHI,ZUFIS,ALBPGR,
     1 ALBP,OUTR,ESTOP,DEN3
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: OUTSC
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ACCUM
*----
*  DATA STATEMENT
*----
      DATA HREAC/'NTOT0','SIGW00','NUSIGF','NFTOT','H-FACTOR',
     1 'OVERV','DIFF','DIFFX','DIFFY','DIFFZ','C-FACTOR'/
*----
*  SCRATCH STORAGE ALLOCATION
*   OUTR(IBM,NREAC+1): volume
*   OUTR(IBM,NREAC+2): integrated direct flux
*   OUTR(IBM,NREAC+3): fission spectrum
*   OUTR(IBM,NREAC+4): fixed sources
*----
      ALLOCATE(VOLI(NZS),WORK(NZS),RATE(NZS),FLINT(NZS,NGRP),
     1 CHI(NBMIX,NBFIS),ZUFIS(NBMIX,NBFIS),OUTR(NZS+1,NREAC+4),
     2 OUTSC(NZS,NL+1,NGCOND),GAR(NGRP),ALBPGR(NALBP,NGRP),
     3 ALBP(NALBP,NGCOND),ESTOP(NZS,NGRP+1))
      ALLOCATE(ACCUM(NZS,NBFIS))
*
      CALL XDRSET(ALBP,NALBP*NGCOND,0.0)
      CALL XDRSET(ESTOP,NZS*(NGRP+1),0.0)
      LNUSIG=.FALSE.
      LESTOP=.FALSE.
      LFIXE=.FALSE.
      CALL XDLSET(LREAC,NREAC,.FALSE.)
*----
*  RECOVER PHYSICAL ALBEDOS.
*----
      IF(NALBP.GT.0) CALL LCMGET(IPMAC1,'ALBEDO',ALBPGR)
*----
*  DIRECT FLUX CALCULATION.
*----
      CALL XDRSET(VOLI,NZS,0.0)
      CALL XDRSET(FLINT,NZS*NGRP,0.0)
      DO 20 K=1,NEL
      IBM=IHOM(K)
      IPFL=IDL(K)
      IF((IBM.NE.0).AND.(MAT(K).NE.0).AND.(IPFL.NE.0)) THEN
         VOLI(IBM)=VOLI(IBM)+VOL(K)
         DO 10 IGR=1,NGRP
         FLINT(IBM,IGR)=FLINT(IBM,IGR)+EVECT(IPFL,IGR)*VOL(K)
   10    CONTINUE
      ENDIF
   20 CONTINUE
      CALL LCMPUT(IPMAC2,'VOLUME',NZS,2,VOLI)
*----
*  FISSION RATE CALCULATION.
*----
      IF(IMPX.GT.0) WRITE(6,'(/35H OUTAUX: REACTION RATE CALCULATION.)')
      JPMAC1=LCMGID(IPMAC1,'GROUP')
      JPMAC2=LCMLID(IPMAC2,'GROUP',NGCOND)
      IF(NBFIS.GT.0) THEN
         DO 70 IFISS=1,NBFIS
         CALL XDDSET(ACCUM(1,IFISS),NZS,0.0D0)
   70    CONTINUE
         DO 100 IGR=1,NGRP
         KPMAC1=LCMGIL(JPMAC1,IGR)
         CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
         DO 90 IFISS=1,NBFIS
         DO 80 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            ACCUM(IBM,IFISS)=ACCUM(IBM,IFISS)+EVECT(IPFL,IGR)*VOL(K)*
     1      ZUFIS(L,IFISS)
         ENDIF
   80    CONTINUE
   90    CONTINUE
  100    CONTINUE
      ENDIF
*----
*  LOOP OVER ENERGY GROUP LIST.
*----
      IGRFIN=0
      DO 500 IGRC=1,NGCOND
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRC)
      CALL XDRSET(OUTR,(NZS+1)*(NREAC+4),0.0)
      CALL XDRSET(OUTSC,NZS*(NL+1)*NGCOND,0.0)
      DO 350 IGR=IGRDEB,IGRFIN
      KPMAC1=LCMGIL(JPMAC1,IGR)
      DO 110 IBM=1,NZS
      OUTR(IBM,NREAC+2)=OUTR(IBM,NREAC+2)+FLINT(IBM,IGR)
  110 CONTINUE
*----
*  SET VOLUMES.
*----
      DO 120 IBM=1,NZS
      OUTR(IBM,NREAC+1)=VOLI(IBM)
  120 CONTINUE
*----
*  REACTION RATE CALCULATION.
*----
      DO 150 IREAC=1,NREAC
      CALL LCMLEN(KPMAC1,HREAC(IREAC),LENGT,ITYLCM)
      LREAC(IREAC)=LREAC(IREAC).OR.(LENGT.NE.0)
      IF((HREAC(IREAC).EQ.'H-FACTOR').AND.(LENGT.EQ.0)) THEN
         WRITE(6,'(/46H OUTAUX: *** WARNING *** NO H-FACTOR FOUND ON ,
     1   25HLCM. USE NU*SIGF INSTEAD.)')
         LNUSIG=.TRUE.
         GO TO 150
      ELSE IF(HREAC(IREAC).EQ.'NUSIGF') THEN
         GO TO 150
      ELSE IF(HREAC(IREAC).EQ.'SIGW00') THEN
         GO TO 150
      ELSE
         TEXT12=HREAC(IREAC)
      ENDIF
      IF(LENGT.GT.0) THEN
         IF(LENGT.GT.NBMIX) CALL XABORT('OUTAUX: INVALID LENGTH FOR '//
     1   HREAC(IREAC)//' CROSS SECTIONS.')
         CALL LCMGET(KPMAC1,TEXT12,WORK)
         CALL XDRSET(RATE,NZS,0.0)
         DO 130 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            RATE(IBM)=RATE(IBM)+EVECT(IPFL,IGR)*VOL(K)*WORK(L)
         ENDIF
  130    CONTINUE
         DO 140 IBM=1,NZS
         OUTR(IBM,IREAC)=OUTR(IBM,IREAC)+RATE(IBM)
  140    CONTINUE
      ENDIF
  150 CONTINUE
*----
*  FIXED SOURCES
*----
      CALL LCMLEN(KPMAC1,'FIXE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         LFIXE=.TRUE.
         IF(LENGT.GT.NBMIX) CALL XABORT('OUTAUX: INVALID LENGTH FOR '//
     1   'FIXE SOURCE.')
         CALL LCMGET(KPMAC1,'FIXE',WORK)
         CALL XDRSET(RATE,NZS,0.0)
         DO 160 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            RATE(IBM)=RATE(IBM)+VOL(K)*WORK(L)
         ENDIF
  160    CONTINUE
         DO 170 IBM=1,NZS
         OUTR(IBM,NREAC+4)=OUTR(IBM,NREAC+4)+RATE(IBM)
  170    CONTINUE
      ENDIF
*----
*  SCATTERING MATRIX INFORMATION IGR <-- JGR.
*----
      ALLOCATE(IJJ(NBMIX),NJJ(NBMIX),IPOS(NBMIX))
      ALLOCATE(SCAT(NBMIX*NGRP))
      DO 220 IL=1,NL
      WRITE(SUFF,'(I2.2)') IL-1
      CALL LCMLEN(KPMAC1,'NJJS'//SUFF,LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         IF(LENGT.GT.NBMIX) CALL XABORT('OUTAUX: INVALID LENGTH FOR '//
     1   'SCATTERING CROSS SECTIONS.')
         CALL LCMLEN(KPMAC1,'SCAT'//SUFF,LENGT,ITYLCM)
         IF(LENGT.GT.NBMIX*NGRP) CALL XABORT('OUTAUX: SCAT OVERFLOW.')
         CALL LCMGET(KPMAC1,'NJJS'//SUFF,NJJ)
         CALL LCMGET(KPMAC1,'IJJS'//SUFF,IJJ)
         CALL LCMGET(KPMAC1,'IPOS'//SUFF,IPOS)
         CALL LCMGET(KPMAC1,'SCAT'//SUFF,SCAT)
         IPOSDE=0
         DO 210 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            CALL XDRSET(GAR,NGRP,0.0)
            IPOSDE=IPOS(L)-1
            DO 180 JGR=IJJ(L),IJJ(L)-NJJ(L)+1,-1
            IPOSDE=IPOSDE+1
            GAR(JGR)=SCAT(IPOSDE)
  180       CONTINUE
            JGRFIN=0
            DO 200 JGRC=1,NGCOND
            JGRDEB=JGRFIN+1
            JGRFIN=IGCOND(JGRC)
            DO 190 JGR=JGRDEB,JGRFIN
            OUTSC(IBM,IL,JGRC)=OUTSC(IBM,IL,JGRC)+EVECT(IPFL,JGR)*
     1      VOL(K)*GAR(JGR)
  190       CONTINUE
  200       CONTINUE
         ENDIF
  210    CONTINUE
         IF(IL.EQ.1) OUTR(:NZS,2)=OUTSC(:NZS,IL,IGRC)
      ENDIF
  220 CONTINUE
      DEALLOCATE(SCAT)
      DEALLOCATE(IJJ,NJJ,IPOS)
*
      DO 250 K=1,NEL
      IBM=IHOM(K)
      L=MAT(K)
      IPFL=IDL(K)
      IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
        JGRFIN=0
        DO 240 JGRC=1,NGCOND
        JGRDEB=JGRFIN+1
        JGRFIN=IGCOND(JGRC)
        DO 230 JGR=JGRDEB,JGRFIN
        OUTSC(IBM,NL+1,JGRC)=OUTSC(IBM,NL+1,JGRC)+EVECT(IPFL,JGR)*VOL(K)
  230   CONTINUE
  240   CONTINUE
      ENDIF
  250 CONTINUE
*----
*  FISSION SPECTRUM AND NUSIGF HOMOGENIZATION.
*----
      IF(NBFIS.GT.0) THEN
         CALL LCMLEN(KPMAC1,'NUSIGF',LENGT,ITYLCM)
         IF(LENGT.NE.NBMIX*NBFIS) CALL XABORT('OUTAUX: INVALID LENGTH '
     1   //'FOR FISSION SPECTRUM.')
         CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
         CALL XDRSET(RATE,NZS,0.0)
         CALL LCMLEN(KPMAC1,'CHI',LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            IF(IGR.EQ.IGRDEB) THEN
               CALL XDRSET(OUTR(1,NREAC+3),NZS,1.0)
            ENDIF
         ELSE
            ALLOCATE(DEN(NZS))
            CALL LCMGET(KPMAC1,'CHI',CHI)
            CALL XDRSET(DEN,NZS,0.0)
            DO 270 K=1,NEL
            IBM=IHOM(K)
            L=MAT(K)
            IF((IBM.NE.0).AND.(L.NE.0)) THEN
               DO 260 IFISS=1,NBFIS
               RATE(IBM)=RATE(IBM)+CHI(L,IFISS)*REAL(ACCUM(IBM,IFISS))
               DEN(IBM)=DEN(IBM)+REAL(ACCUM(IBM,IFISS))
  260          CONTINUE
            ENDIF
  270       CONTINUE
            DO 280 IBM=1,NZS
            IF(DEN(IBM).NE.0.0) OUTR(IBM,NREAC+3)=RATE(IBM)/DEN(IBM)
  280       CONTINUE
            DEALLOCATE(DEN)
         ENDIF
         DO 300 IFISS=1,NBFIS
         DO 290 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
           OUTR(IBM,3)=OUTR(IBM,3)+EVECT(IPFL,IGR)*VOL(K)*ZUFIS(L,IFISS)
         ENDIF
  290    CONTINUE
  300    CONTINUE
      ENDIF
*----
*  CONDENSE PHYSICAL ALBEDOS.
*----
      IF(NALBP.GT.0) THEN
         DO 320 IAL=1,NALBP
         DO 310 IBM=1,NZS
         ALBP(IAL,IGRC)=ALBP(IAL,IGRC)+ALBPGR(IAL,IGR)*FLINT(IBM,IGR)
  310    CONTINUE
  320    CONTINUE
      ENDIF
*----
*  RECOVER AND HOMOGENIZE STOPPING POWERS
*----
      CALL LCMLEN(KPMAC1,'ESTOPW',LENGT,ITYLCM)
      IF(LENGT.EQ.2*NBMIX) THEN
         ALLOCATE(DEN3(NBMIX,2))
         LESTOP=.TRUE.
         CALL LCMGET(KPMAC1,'ESTOPW',DEN3)
         DO 330 K=1,NEL
         IBM=IHOM(K)
         L=MAT(K)
         IPFL=IDL(K)
         IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
            IF(IGR.EQ.1) THEN
               FACTOR=EVECT(IPFL,IGR)/FLINT(IBM,IGR)
            ELSE
               FACTOR=(EVECT(IPFL,IGR-1)+EVECT(IPFL,IGR))/
     1         (FLINT(IBM,IGR-1)+FLINT(IBM,IGR))
            ENDIF
            ESTOP(IBM,IGR)=ESTOP(IBM,IGR)+FACTOR*VOL(K)*DEN3(L,1)
         ENDIF
  330    CONTINUE
         IF(IGR.EQ.NGRP) THEN
            DO 340 K=1,NEL
            IBM=IHOM(K)
            L=MAT(K)
            IPFL=IDL(K)
            IF((IBM.NE.0).AND.(L.NE.0).AND.(IPFL.NE.0)) THEN
               FACTOR=EVECT(IPFL,IGR)/FLINT(IBM,IGR)
               ESTOP(IBM,IGR+1)=ESTOP(IBM,IGR+1)+FACTOR*VOL(K)*DEN3(L,2)
            ENDIF
  340       CONTINUE
         ENDIF
         DEALLOCATE(DEN3)
      ENDIF
  350 CONTINUE
*----
*  PRINT THE REACTION RATES:
*----
      IF(IMPX.GT.0) THEN
         DO 360 I=1,NREAC+3
         OUTR(NZS+1,I)=0.0
  360    CONTINUE
         WRITE(6,520) IGRC,'VOLUME      ','FLUX-INTG   ',
     1   (HREAC(I),I=1,6),'CHI         '
         DO 380 IBM=1,NZS
         DO 370 I=1,NREAC+3
         OUTR(NZS+1,I)=OUTR(NZS+1,I)+OUTR(IBM,I)
  370    CONTINUE
         WRITE(6,530) IBM,OUTR(IBM,NREAC+1),OUTR(IBM,NREAC+2),
     1   (OUTR(IBM,I),I=1,6),OUTR(IBM,NREAC+3)
  380    CONTINUE
         WRITE(6,540) OUTR(NZS+1,NREAC+1),OUTR(NZS+1,NREAC+2),
     1   (OUTR(NZS+1,I),I=1,6)
      ENDIF
*----
*  COMPUTE HOMOGENIZED-CONDENSED MACROLIB
*----
      KPMAC2=LCMDIL(JPMAC2,IGRC)
      CALL LCMPUT(KPMAC2,'FLUX-INTG',NZS,2,OUTR(1,NREAC+2))
      DO 400 IREAC=1,NREAC
      IF(LREAC(IREAC)) THEN
         DO 390 IBM=1,NZS
         RATE(IBM)=OUTR(IBM,IREAC)
         IF(RATE(IBM).NE.0.0) RATE(IBM)=RATE(IBM)/OUTR(IBM,NREAC+2)
  390    CONTINUE
         CALL LCMPUT(KPMAC2,HREAC(IREAC),NZS,2,RATE)
         IF(LNUSIG.AND.(IREAC.EQ.3)) THEN
            CALL LCMPUT(KPMAC2,'H-FACTOR',NZS,2,RATE)
         ENDIF
      ENDIF
  400 CONTINUE
      IF(LREAC(3)) CALL LCMPUT(KPMAC2,'CHI',NZS,2,OUTR(1,NREAC+3))
      IF(LFIXE) THEN
         DO 410 IBM=1,NZS
         RATE(IBM)=OUTR(IBM,NREAC+4)
         IF(RATE(IBM).NE.0.0) RATE(IBM)=RATE(IBM)/VOLI(IBM)
  410    CONTINUE
         CALL LCMPUT(KPMAC2,'FIXE',NZS,2,RATE)
      ENDIF
*
      ALLOCATE(IJJ(NZS),NJJ(NZS),IPOS(NZS))
      ALLOCATE(SCAT(NZS*NGCOND))
      DO 460 IL=1,NL
      WRITE(SUFF,'(I2.2)') IL-1
      DO 430 IBM=1,NZS
      IGMIN=IGRC
      IGMAX=IGRC
      DO 420 JGRC=NGCOND,1,-1
      IF(OUTSC(IBM,IL,JGRC).NE.0.0) THEN
         IGMIN=MIN(IGMIN,JGRC)
         IGMAX=MAX(IGMAX,JGRC)
         OUTSC(IBM,IL,JGRC)=OUTSC(IBM,IL,JGRC)/OUTSC(IBM,NL+1,JGRC)
      ENDIF
  420 CONTINUE
      IJJ(IBM)=IGMAX
      NJJ(IBM)=IGMAX-IGMIN+1
  430 CONTINUE
      IPOSDE=0
      DO 450 IBM=1,NZS
      IPOS(IBM)=IPOSDE+1
      DO 440 JGRC=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
      IPOSDE=IPOSDE+1
      SCAT(IPOSDE)=OUTSC(IBM,IL,JGRC)
  440 CONTINUE
  450 CONTINUE
      CALL LCMPUT(KPMAC2,'SCAT'//SUFF,IPOSDE,2,SCAT)
      CALL LCMPUT(KPMAC2,'IPOS'//SUFF,NZS,1,IPOS)
      CALL LCMPUT(KPMAC2,'NJJS'//SUFF,NZS,1,NJJ)
      CALL LCMPUT(KPMAC2,'IJJS'//SUFF,NZS,1,IJJ)
      CALL LCMPUT(KPMAC2,'SIGW'//SUFF,NZS,2,OUTSC(1,IL,IGRC))
  460 CONTINUE
      DEALLOCATE(SCAT)
      DEALLOCATE(IJJ,NJJ,IPOS)
*
      IF(NALBP.GT.0) THEN
         DFI=0.0
         DO 470 IBM=1,NZS
         DFI=DFI+OUTR(IBM,NREAC+2)
  470    CONTINUE
         DO 480 IAL=1,NALBP
         ALBP(IAL,IGRC)=ALBP(IAL,IGRC)/DFI
  480    CONTINUE
      ENDIF
*----
*  SAVE STOPPING POWERS
*----
      IF(LESTOP) THEN
         ALLOCATE(DEN3(NZS,2))
         DO 490 IBM=1,NZS
         IF(IGRC.EQ.1) THEN
            DEN3(IBM,1)=ESTOP(IBM,1)
         ELSE
            DEN3(IBM,1)=ESTOP(IBM,IGCOND(IGRC-1))
         ENDIF
         DEN3(IBM,2)=ESTOP(IBM,IGCOND(IGRC)+1)
  490    CONTINUE
         CALL LCMPUT(KPMAC2,'ESTOPW',NZS*2,2,DEN3)
         DEALLOCATE(DEN3)
      ENDIF
  500 CONTINUE
*----
*  END OF LOOP OVER MACROGROUPS
*----
*----
*  RECOVER AND CONDENSE ENERGY MESH
*----
      CALL LCMLEN(IPMAC1,'ENERGY',LENGT,ITYLCM)
      IF(LENGT.EQ.NGRP+1) THEN
         ALLOCATE(DEN(NGRP+1),DEN2(NGCOND+1))
         CALL LCMGET(IPMAC1,'ENERGY',DEN)
         DEN2(1)=DEN(1)
         DO 510 IGRC=1,NGCOND
         DEN2(IGRC+1)=DEN(IGCOND(IGRC)+1)
  510    CONTINUE
         CALL LCMPUT(IPMAC2,'ENERGY',NGCOND+1,2,DEN2)
         DEALLOCATE(DEN2,DEN)
      ENDIF
*----
*  SAVE ALBEDO AND STATE-VECTOR
*----
      IF(NALBP.GT.0) THEN
         CALL LCMPUT(IPMAC2,'ALBEDO',NALBP*NGCOND,2,ALBP)
      ENDIF
      CALL LCMLEN(IPMAC1,'PARTICLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
        CALL LCMGTC(IPMAC1,'PARTICLE',12,1,TEXT6)
        CALL LCMPTC(IPMAC2,'PARTICLE',12,1,TEXT6)
      ENDIF
      CALL XDISET(IDATA,NSTATE,0)
      IDATA(1)=NGCOND
      IDATA(2)=NZS
      IDATA(3)=NL
      IDATA(4)=1
      IDATA(8)=NALBP
      IF(LREAC(7)) THEN
         IDATA(9)=1
      ELSE IF(LREAC(8)) THEN
         IDATA(9)=2
      ENDIF
      IDATA(15)=0
      CALL LCMPUT(IPMAC2,'STATE-VECTOR',NSTATE,1,IDATA)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ACCUM)
      DEALLOCATE(ESTOP,ALBP,ALBPGR,GAR,OUTSC,OUTR,ZUFIS,CHI,FLINT,
     1 RATE,WORK,VOLI)
      RETURN
*
  520 FORMAT(/' G R O U P   : ',I3/1X,'IHOM',9A14)
  530 FORMAT(1X,I4,1P,9E14.5)
  540 FORMAT(/5H  SUM,1P,8E14.5)
      END
