*DECK USSIN1
      SUBROUTINE USSIN1(IPLI0,IPLIB,NGRP,NBMIX,NBISO,NIRES,NBNRS,NL,
     1 NED,NDEL,IREX,IMPX,ISONAM,ISOBIS,MIX,IAPT,MASKI,SPH,PHGAR,STGAR,
     2 SFGAR,SSGAR,S0GAR,SAGAR,SDGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write the self-shielded and SPH-corrected cross sections on the
* internal library.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module (L_LIBRARY signature).
* IPLIB   pointer to the internal microscopic cross section library
*         with subgroups (L_LIBRARY signature).
* NGRP    number of energy groups.
* NBMIX   number of mixtures in the internal library.
* NBISO   number of isotopes.
* NIRES   number of resonant isotopes in fuel regions.
* NBNRS   number of totally correlated fuel regions.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* IREX    fuel region index assigned to each mixture (equal to zero
*         in non-resonant mixtures or in mixtures not used).
* IMPX    print flag. equal to zero for no print.
* ISONAM  alias name of isotopes in IPLIB.
* ISOBIS  alias name of isotopes in IPLI0.
* MIX     mix number of each isotope (can be zero).
* IAPT    resonant isotope index associated with isotope I. Mixed
*         moderator if IAPT(I)=NIRES+1. Out-of-fuel isotope if
*         IAPT(I)=0.
* MASKI   isotopic flag (MASKI(ISO)=.TRUE. to process isotope ISO).
* SPH     SPH factors.
* PHGAR   averaged fluxes.
* STGAR   microscopic self-shielded total x-s.
* SFGAR   microscopic self-shielded fission x-s.
* SSGAR   microscopic self-shielded scattering x-s.
* S0GAR   microscopic transfer scattering xs (isotope,secondary,
*         primary).
* SAGAR   microscopic self-shielded additional xs.
* SDGAR   microscopic self-shielded delayed nu-sigf xs.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPLIB
      INTEGER NGRP,NBMIX,NBISO,NIRES,NBNRS,NL,NED,NDEL,IREX(NBMIX),
     1 IMPX,ISONAM(3,NBISO),ISOBIS(3,NBISO),MIX(NBISO),IAPT(NBISO)
      REAL SPH(NBNRS,NIRES,NGRP),PHGAR(NBNRS,NIRES,NGRP),
     1 STGAR(NBNRS,NIRES,NGRP),SFGAR(NBNRS,NIRES,NGRP),
     2 SSGAR(NBNRS,NIRES,NL,NGRP),S0GAR(NBNRS,NIRES,NL,NGRP,NGRP),
     3 SAGAR(NBNRS,NIRES,NED,NGRP),SDGAR(NBNRS,NIRES,NDEL,NGRP)
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXED=50,MAXSAV=100,MAXESP=4)
      TYPE(C_PTR) JPLI0,KPLI0,KPLIB
      INTEGER IPAR(40),ISAV(MAXSAV),IESP(MAXESP+1)
      REAL EESP(MAXESP+1)
      CHARACTER TEXT12*12,HSIGN*12,CM*2,HVECT(MAXED)*8,HCHI*12
      LOGICAL LOGNF,LTEST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISONR,ITITLE
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,ENERGY,LAMB
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR2
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAR1(NGRP),GAR2(NGRP,NGRP),ENERGY(NGRP+1),IPISO(NBISO))
*
      CALL LCMGTC(IPLIB,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_LIBRARY') THEN
         CALL XABORT('USSIN1: SIGNATURE IS '//HSIGN//'. L_LIBRARY EXPEC'
     1   //'TED.')
      ENDIF
      CALL LCMGET(IPLIB,'STATE-VECTOR',IPAR)
      IF(NGRP.NE.IPAR(3)) CALL XABORT('USSIN1: INVALID NB OF GROUPS.')
      IF(NL.NE.IPAR(4)) CALL XABORT('USSIN1: INVALID VALUE OF NL.')
      IF(NED.NE.IPAR(13)) CALL XABORT('USSIN1: INVALID VALUE OF NED.')
      IF(NED.GT.0) THEN
         IF(NED.GT.MAXED) CALL XABORT('USSIN1: INVALID VALUE OF MAXED.')
         CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED,HVECT)
      ENDIF
      CALL LCMLEN(IPLIB,'ENERGY',LENGT,ITYLCM)
      IF(LENGT-1.NE.NGRP) CALL XABORT('LIBIN2: INVALID GROUP STRUCTU'
     1 //'RE.')
      CALL LCMGET(IPLIB,'ENERGY',ENERGY)
      CALL LCMPUT(IPLI0,'ENERGY',NGRP+1,2,ENERGY)
      CALL LCMGET(IPLIB,'DELTAU',ENERGY)
      CALL LCMPUT(IPLI0,'DELTAU',NGRP,2,ENERGY)
      CALL LCMLEN(IPLIB,'CHI-LIMITS',NBESP,ITYLCM)
      IF(NBESP.GT.0) THEN
         NBESP=NBESP-1
         IF(NBESP.GT.MAXESP) CALL XABORT('USSIN1: MAXESP OVERFLOW.')
         CALL LCMGET(IPLIB,'CHI-LIMITS',IESP)
         CALL LCMPUT(IPLI0,'CHI-LIMITS',NBESP+1,1,IESP)
         CALL LCMGET(IPLIB,'CHI-ENERGY',EESP)
         CALL LCMPUT(IPLI0,'CHI-ENERGY',NBESP+1,2,EESP)
      ENDIF
      ALLOCATE(ISONR(3*NBISO))
      CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONR)
      CALL LCMPUT(IPLI0,'ISOTOPERNAME',3*NBISO,3,ISONR)
      DEALLOCATE(ISONR)
      CALL LIBIPS(IPLIB,NBISO,IPISO)
      JPLI0=LCMLID(IPLI0,'ISOTOPESLIST',NBISO)
      DO 370 ISOT=1,NBISO
      IRES=IAPT(ISOT)
      IF(MASKI(ISOT).AND.(IRES.GT.0).AND.(IRES.LE.NIRES)) THEN
         KPLIB=IPISO(ISOT) ! set ISOT-th isotope
         WRITE(TEXT12,'(3A4)') (ISOBIS(J,ISOT),J=1,3)
         IF(IMPX.GT.0) WRITE (6,'(/29H USSIN1: PROCESSING ISOTOPE '',
     1   A12,2H''.)') TEXT12
         KPLI0=LCMDIL(JPLI0,ISOT) ! set ISOT-th isotope
         CALL LCMPTC(KPLI0,'ALIAS',12,1,TEXT12)
         CALL LCMGET(KPLIB,'AWR',AWR)
         CALL LCMPUT(KPLI0,'AWR',1,2,AWR)
         CALL LCMLEN(KPLIB,'README',LENTIT,ITYLCM)
         IF(LENTIT.GT.0) THEN
            ALLOCATE(ITITLE(LENTIT))
            CALL LCMGET(KPLIB,'README',ITITLE)
            CALL LCMPUT(KPLI0,'README',LENTIT,3,ITITLE)
            DEALLOCATE(ITITLE)
         ENDIF
         CALL LCMLEN(KPLIB,'NUSIGF',NFIS,ITYLCM)
         LOGNF=(NFIS.GT.0)
         IF(LOGNF) THEN
            IF(NBESP.EQ.0) THEN
              CALL LCMGET(KPLIB,'CHI',GAR1)
              CALL LCMPUT(KPLI0,'CHI',NGRP,2,GAR1)
            ELSE
              DO ISP=1,NBESP
                WRITE(HCHI,'(5HCHI--,I2.2)') ISP
                CALL LCMLEN(KPLIB,HCHI,ILONG,ITYLCM)
                IF(ILONG.EQ.NGRP) THEN
                  CALL LCMGET(KPLIB,HCHI,GAR1)
                  CALL LCMPUT(KPLI0,HCHI,NGRP,2,GAR1)
                ENDIF
              ENDDO
            ENDIF
         ENDIF
         IND=IREX(MIX(ISOT))
         IF(IND.EQ.0) CALL XABORT('USSIN1: IREX FAILURE.')
         DO 20 IG1=1,NGRP
         GAR1(IG1)=PHGAR(IND,IRES,IG1)
   20    CONTINUE
         CALL LCMPUT(KPLI0,'NWT0',NGRP,2,GAR1)
         DO 30 IG1=1,NGRP
         GAR1(IG1)=SPH(IND,IRES,IG1)
   30    CONTINUE
         CALL LCMPUT(KPLI0,'NSPH',NGRP,2,GAR1)
         DO 40 IG1=1,NGRP
         GAR1(IG1)=STGAR(IND,IRES,IG1)
   40    CONTINUE
         CALL LCMPUT(KPLI0,'NTOT0',NGRP,2,GAR1)
         IF(LOGNF) THEN
            DO 50 IG1=1,NGRP
            GAR1(IG1)=SFGAR(IND,IRES,IG1)
   50       CONTINUE
            CALL LCMPUT(KPLI0,'NUSIGF',NGRP,2,GAR1)
         ENDIF
         DO 90 IL=1,NL
         WRITE(CM,'(I2.2)') IL-1
         DO 70 IG1=1,NGRP
         GAR1(IG1)=SSGAR(IND,IRES,IL,IG1)
         DO 60 IG2=1,NGRP
         GAR2(IG2,IG1)=S0GAR(IND,IRES,IL,IG2,IG1)
   60    CONTINUE
   70    CONTINUE
         CALL XDRLGS(KPLI0,1,IMPX,IL-1,IL-1,1,NGRP,GAR1,GAR2,ITYPRO)
   90    CONTINUE
         CALL LCMLEN(KPLIB,'XS-SAVED',ILENG,ITYLCM)
         IF(ILENG.GT.MAXSAV) CALL XABORT('USSIN1: XS-SAVED OVERFLOW.')
         IF(ILENG.GT.0) CALL LCMGET(KPLIB,'XS-SAVED',ISAV)
         IF(ILENG.GT.0) CALL LCMPUT(KPLI0,'XS-SAVED',ILENG,1,ISAV)
         CALL LCMLEN(KPLIB,'SCAT-SAVED',ILENG,ITYLCM)
         IF(ILENG.GT.MAXSAV) CALL XABORT('USSIN1: SCAT-SAVED OVERFLOW.')
         IF(ILENG.GT.0) CALL LCMGET(KPLIB,'SCAT-SAVED',ISAV)
         IF(ILENG.GT.0) CALL LCMPUT(KPLI0,'SCAT-SAVED',ILENG,1,ISAV)
         DO 110 IED=1,NED
         CALL LCMLEN(KPLIB,HVECT(IED),NEDI,ITYLCM)
         IF((NEDI.GT.0).AND.(HVECT(IED)(:3).NE.'CHI').AND.
     1   (HVECT(IED)(:2).NE.'NU').AND.(HVECT(IED).NE.'NGOLD').AND.
     2   (HVECT(IED)(:3).NE.'NWT').AND.(HVECT(IED).NE.'NTOT0')) THEN
            DO 100 IG1=1,NGRP
            GAR1(IG1)=SAGAR(IND,IRES,IED,IG1)
  100       CONTINUE
            CALL LCMPUT(KPLI0,HVECT(IED),NGRP,2,GAR1)
         ENDIF
  110    CONTINUE
         CALL LCMLEN(KPLIB,'NUSIGF01',ILONG,ITYLCM)
         IF(ILONG.EQ.NGRP) THEN
            CALL LCMLEN(KPLIB,'LAMBDA-D',ILONG,ITYLCM)
            IF(ILONG.EQ.0) THEN
               WRITE(TEXT12,'(3A4)') (ISONAM(J,ISOT),J=1,3)
               CALL XABORT('USSIN1: MISSING LAMBDA-D INFO '//'FOR '//
     1         TEXT12//'.')
            ENDIF
            ALLOCATE(LAMB(ILONG))
            CALL LCMGET(KPLIB,'LAMBDA-D',LAMB)
            CALL LCMPUT(KPLI0,'LAMBDA-D',ILONG,2,LAMB)
            DEALLOCATE(LAMB)
            DO 130 IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            LTEST=.FALSE.
            DO 120 IG1=1,NGRP
            GAR1(IG1)=SDGAR(IND,IRES,IDEL,IG1)
            LTEST=LTEST.OR.(GAR1(IG1).NE.0.0)
  120       CONTINUE
            IF(LTEST) THEN
               CALL LCMPUT(KPLI0,TEXT12,NGRP,2,GAR1)
               WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
               CALL LCMGET(KPLIB,TEXT12,GAR1)
               CALL LCMPUT(KPLI0,TEXT12,NGRP,2,GAR1)
            ENDIF
  130       CONTINUE
         ENDIF
*
         IF(IMPX.GT.2) THEN
            WRITE (6,'(/13H SPH FACTORS:/(1X,1P,10E12.4))')
     1      (SPH(IND,IRES,I),I=1,NGRP)
            CALL LCMGET(KPLI0,'NTOT0',GAR1)
            WRITE (6,'(/36H SELF-SHIELDED MICROSCOPIC TOTAL XS:/
     1      (1X,1P,10E12.4))') (GAR1(I),I=1,NGRP)
            DO 350 IL=1,NL
            WRITE(CM,'(I2.2)') IL-1
            CALL LCMLEN(KPLI0,'SCAT'//CM,ILSCAT,ITYLCM)
            IF(ILSCAT.GT.NGRP**2) CALL XABORT('USSIN1: OVERWRITING ME'
     1      //'MORY(2).')
            IF((IL.EQ.1).OR.(ILSCAT.GT.0)) THEN
               CALL LCMGET(KPLI0,'SIGS'//CM,GAR1)
               WRITE (6,'(/16H SELF-SHIELDED P,A2,18H MICROSCOPIC SCATT,
     1         9HERING XS:/(1X,1P,10E12.4))') CM,(GAR1(I),I=1,NGRP)
            ENDIF
  350       CONTINUE
            IF(LOGNF) THEN
               CALL LCMGET(KPLI0,'NUSIGF',GAR1)
               WRITE (6,'(/38H SELF-SHIELDED MICROSCOPIC FISSION XS:/
     1         (1X,1P,10E12.4))') (GAR1(I),I=1,NGRP)
               IF(NBESP.EQ.0) THEN
                 CALL LCMGET(KPLI0,'CHI',GAR1)
                 WRITE (6,'(/18H FISSION SPECTRUM:/
     1           (1X,1P,10E12.4))') (GAR1(I),I=1,NGRP)
               ELSE
                 DO 355 ISP=1,NBESP
                 WRITE(HCHI,'(5HCHI--,I2.2)') ISP
                 CALL LCMLEN(KPLI0,HCHI,ILONG,ITYLCM)
                 IF(ILONG.EQ.NGRP) THEN
                   CALL LCMGET(KPLI0,HCHI,GAR1)
                   WRITE (6,'(/I3,21H-TH FISSION SPECTRUM:/
     1             (1X,1P,10E12.4))') ISP,(GAR1(I),I=1,NGRP)
                 ENDIF
  355            CONTINUE
               ENDIF
            ENDIF
            DO 360 IED=1,NED
            CALL LCMLEN(KPLI0,HVECT(IED),NEDI,ITYLCM)
            IF((NEDI.GT.0).AND.(HVECT(IED)(:3).NE.'CHI').AND.
     1      (HVECT(IED)(:2).NE.'NU').AND.(HVECT(IED).NE.'NGOLD').AND.
     2      (HVECT(IED)(:3).NE.'NWT').AND.(HVECT(IED).NE.'NTOT0')) THEN
               CALL LCMGET(KPLI0,HVECT(IED),GAR1)
               WRITE (6,'(/15H SELF-SHIELDED ,A6,1H:/(1X,1P,10E12.4))')
     1         HVECT(IED),(GAR1(I),I=1,NGRP)
            ENDIF
  360       CONTINUE
         ENDIF
      ENDIF
  370 CONTINUE
      IF(IMPX.GT.2) CALL LCMLIB(IPLI0)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO,ENERGY,GAR2,GAR1)
      RETURN
      END
