*DECK AFMDRV
      SUBROUTINE AFMDRV (KENTRY,NENTRY,NPARM,ITYPE,NBURN,NGRP,NISO,ISC,
     1 MNPS,NL,ILEAK,NTYP,NBCH,NCCO,NCZO,NUT,CTITRE,LMCR,IXYZ,MMIX,MSFT,
     2 NISM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver to generate a macrolib using fbm
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* M.T. Sissaoui
*
*Update(s):
*  E. Varin 28/03/00, B. Dionne 26/02/01, 
*  A. Lagarrigue 30/07/05
*  A. Hebert 11/11/11 (remove table support)
*
*Parameters: input
* KENTRY  address of the LCM objects
* NENTRY  number of LCM objects
* NPARM   number of parameters in L_MAP object
* ITYPE   creation/modification flag for output macrolib
* NBURN   number of burnup steps
* NGRP    1+number of energy groups
* NISO    number of extracted isotopes
* ISC     type of cross-section calculation (=1: time average;
*         =2: instantaneous; =3: homogeneous)
* MNPS    number of shifts + 2
* NL      number of legendre orders (=1 for isotropic scattering)
* ILEAK   type of leakage
* NTYP
* NBCH    number of bundles per channel
* NCCO    number of channels in the core
* NCZO    number of combustion zones
* NUT     number of fuel types
* CTITRE  character*72 title
* LMCR    if true, create a macrolib containing only one non-zero
*         mixture
* IXYZ    type of diffusion coefficient (=0: isotropic; =1: directional)
* MMIX    number of mixtures in the output macrolib
* MSFT    second dimension of BSFT and PSFT
* NISM
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KENTRY(NENTRY)
      INTEGER NPARM,ITYPE,NBURN,NGRP,NISO,ISC,MNPS,NL,ILEAK,NTYP,NBCH,
     1 NCCO,NCZO,NUT,IXYZ,MMIX,MSFT,NISM
      CHARACTER*72 CTITRE
      LOGICAL LMCR
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXTR*12,CM*2,TEXT4*5,HMICRO*12,TEXTB*12,TEXTD*12
      TYPE(C_PTR) IPMACX,JPMAC,KPMAC,IPFBM,IPMAP,JPMAP,KPMAP
      DOUBLE PRECISION DFLOTT,XCOF(3)
      REAL  STORE,RLOC(7)
      LOGICAL LNOMP,LTAV,LXENON,LSAM,LNEP,LXEREF,LNEREF,LTFUEL,LDRAH,
     1 LTCOOL,LDCOOL,LPWF,LINI
      CHARACTER PNAME*12,PARKEY*12
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPOS,IJ,IZONE,IWORK,NJ,
     1 HISO,JTAB,INDEX,KTYP,ISFT,ITEXTR
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IJJ,NJJ
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,ENER,WORK,BURBG,BURED,
     1 POWER,PW,BRH,XSIGF,XSIGX,XFLUN,PDCOOL,PTCOOL,PTFUEL,SSCAT
      REAL, DIMENSION(:,:), ALLOCATABLE :: XBURN,OVERV,SIGS,FLUX,CHI,
     1 DIFFX,DIFFY,DIFFZ,FLUAV,BFLUX,BSFT,PSFT
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: SIGMA,SIGAV,DENSITB,HXEN1,
     1 HXEN2,HSAM1,HSAM2,HNEP1,HNEP2,CPW1B,CPW2B,FLUXB,CHIB,OVERVB
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: SCAT,SCATAV
      REAL, DIMENSION(:,:,:,:,:), ALLOCATABLE :: SMACB,XBORB,XXENB,
     1 XT1FB,XT2FB,XT1CB,XT2CB,XT1MB,XT2MB,XD1CB,XD2CB,XD1MB,XD2MB,
     2 XSMB,XNP9B,XMFDB,XMMDB,XPF1B,XPF2B,XPF1LB,XPF2LB,XPURB
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SIGMA(MMIX,NGRP,NTYP),IJJ(MMIX,NL,NGRP),VOL(MMIX),
     1 NJJ(MMIX,NL,NGRP),XBURN(NBURN,NUT),OVERV(MMIX,NGRP),
     2 SIGS(MMIX,NGRP),FLUX(MMIX,NGRP),CHI(MMIX,NGRP),ENER(NGRP+1),
     3 IPOS(MMIX),SCAT(MMIX,NL,NGRP,NGRP),DIFFX(MMIX,NGRP),
     4 DIFFY(MMIX,NGRP),DIFFZ(MMIX,NGRP),IJ(NGRP),WORK(MMIX*NGRP*NBURN),
     5 IZONE(NCCO),BURBG(MMIX),BURED(MMIX),POWER(MMIX),
     6 FLUAV(NBURN,NGRP),SIGAV(NBURN,NGRP,NTYP),IWORK(MMIX*NGRP),
     7 SCATAV(NBURN,NL,NGRP,NGRP),PW(MNPS),BRH(MNPS),NJ(NGRP),
     8 BFLUX(NGRP,MMIX),DENSITB(NISO,NBURN,NUT),HISO(3*NISM),
     9 HXEN1(2,NBURN,NUT),HXEN2(2,NBURN,NUT),HSAM1(2,NBURN,NUT),
     1 HSAM2(2,NBURN,NUT),HNEP1(2,NBURN,NUT),HNEP2(2,NBURN,NUT),
     2 CPW1B(2,NBURN,NUT),CPW2B(2,NBURN,NUT),FLUXB(NGRP,NBURN,NUT),
     3 JTAB(NISO),CHIB(NGRP,NBURN,NUT),OVERVB(NGRP,NBURN,NUT),
     4 INDEX(MMIX),KTYP(NUT),XSIGF(NGRP),XSIGX(NGRP),XFLUN(NGRP),
     5 BSFT(MMIX,MSFT),PSFT(MMIX,MSFT),ISFT(MMIX),PDCOOL(MMIX),
     6 PTCOOL(MMIX),PTFUEL(MMIX),ITEXTR(3*NUT))
       ALLOCATE(SMACB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     1          XBORB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     2          XXENB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     3          XT1FB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     4          XT2FB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     5          XT1CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     6          XT2CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     7          XT1MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     8          XT2MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     9          XD1CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     1          XD2CB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     2          XD1MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     3          XD2MB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     4          XSMB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     5          XNP9B(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     6          XMFDB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     7          XMMDB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     8          XPF1B(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     9          XPF2B(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     1          XPF1LB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     2          XPF2LB(NGRP*NGRP,NTYP,NISO,NBURN,NUT),
     3          XPURB(NGRP*NGRP,NTYP,NISO,NBURN,NUT))
*
      IPMACX=KENTRY(1)
      IPFBM=KENTRY(2)
      IF( .NOT.LMCR )IPMAP=KENTRY(3)
      CALL LCMLIB(IPMAP)
*---------------------------------------------------------------*
*      SET THE DEFAULT OPTIONS
      LNOMP=.FALSE.
      LTAV=.FALSE.
      LXENON=.FALSE.
      LSAM=.FALSE.
      LNEP=.FALSE.
      LXEREF=.FALSE.
      LNEREF=.FALSE.
      LTFUEL=.FALSE.
      LDRAH =.FALSE.
      LTCOOL=.FALSE.
      LDCOOL=.FALSE.
      LPWF=.TRUE.
      ILBFLU=0
      IMPX=0
      IXENO=0
      ISAMA=0
      INEPT=0
      IPROF2=0
      LINI=.FALSE.
      ILEAK=0
      PWREF=0.0
      DMR=0.0
      DCR=0.0
      NTM=0
* Set burnup interpolation method 
*   (default 0 for lagrangian interpolation)
*   (1 for linear)
      ILIN=0 
*     SET  HERMITE INTERPOLATION FOR TIME-AVERAGE CALCULATION
      ITM=3
*---------------------------------------------------------------*
*      MX IS THE MAXIMUN MIXTURE NUMBER
      MX=NBCH*NCCO
*---------------------------------------------------------------*
*     CHECK THE PARAMETERS
      IF(MX.EQ.0) CALL XABORT('AFMDRV: ZERO NUMBER OF MIXTURES.')
      IF(NGRP.EQ.0) CALL XABORT('AFMDRV: ZERO NUMBER OF GROUPS.')
      IF(NBURN.EQ.0) CALL XABORT('AFMDRV: ZERO NUMBER OF BURNUPS.')
*---------------------------------------------------------------*
*     INITIALISATION OF THE MATRICES
      NG2=NGRP*NGRP
      DO 50 IGR=1,NG2
       DO 40 IN=1,NUT
        DO 30 I=1,NBURN
         DO 20 ITY=1,NTYP
          DO 10 ISO=1,NISO
            XBORB(IGR,ITY,ISO,I,IN)=0.0
            XPURB(IGR,ITY,ISO,I,IN)=0.0
            XXENB(IGR,ITY,ISO,I,IN)=0.0
            XT1FB(IGR,ITY,ISO,I,IN)=0.0
            XT2FB(IGR,ITY,ISO,I,IN)=0.0
            XT1CB(IGR,ITY,ISO,I,IN)=0.0
            XT2CB(IGR,ITY,ISO,I,IN)=0.0
            XT1MB(IGR,ITY,ISO,I,IN)=0.0
            XT2MB(IGR,ITY,ISO,I,IN)=0.0
            XD1CB(IGR,ITY,ISO,I,IN)=0.0
            XD2CB(IGR,ITY,ISO,I,IN)=0.0
            XD1MB(IGR,ITY,ISO,I,IN)=0.0
            XD2MB(IGR,ITY,ISO,I,IN)=0.0
            XSMB(IGR,ITY,ISO,I,IN)=0.0
            XNP9B(IGR,ITY,ISO,I,IN)=0.0
            XMFDB(IGR,ITY,ISO,I,IN)=0.0
            XMMDB(IGR,ITY,ISO,I,IN)=0.0
            XPF1B(IGR,ITY,ISO,I,IN)=0.0
            XPF2B(IGR,ITY,ISO,I,IN)=0.0
            XPF1LB(IGR,ITY,ISO,I,IN)=0.0
            XPF2LB(IGR,ITY,ISO,I,IN)=0.0
            SMACB(IGR,ITY,ISO,I,IN)=0.0
   10     CONTINUE
   20    CONTINUE
   30   CONTINUE
   40  CONTINUE
   50 CONTINUE
*
      DO 100 IGR=1,NGRP
        DO 90 IMX=1,MX
           DIFFX(IMX,IGR)=0.0
           DIFFY(IMX,IGR)=0.0
           DIFFZ(IMX,IGR)=0.0
           FLUX(IMX,IGR)=0.0
           OVERV(IMX,IGR)=0.0
           CHI(IMX,IGR)=0.0
           DO 70 IL=1,NL
             DO 60 JGR=1,NGRP
               SCAT(IMX,IL,IGR,JGR)=0.0
   60        CONTINUE
             IJJ(IMX,IL,IGR)=IGR
             NJJ(IMX,IL,IGR)=1
   70      CONTINUE
           DO 80 ITYP=1,NTYP
             SIGMA(IMX,IGR,ITYP)=0.0
   80      CONTINUE
   90   CONTINUE
  100 CONTINUE
C
      DO 150 IBR=1,NBURN
        DO 140 IGR=1,NGRP
          FLUAV(IBR,IGR)=0.0
          DO 110 ITYP=1,NTYP
           SIGAV(IBR,IGR,ITYP)=0.0
  110     CONTINUE
          DO 130 JGR=1,NGRP
            DO 120 IL=1,NL
              SCATAV(IBR,IL,IGR,JGR)=0.0
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
  150 CONTINUE
* INITIALISATION OF THE HISTORY COEFFICIENT
      DO 180 IBR=1,NBURN
       DO 170 IN=1,NUT
         DO 160 I=1,2
           CPW1B(I,IBR,IN)=0.0
           CPW2B(I,IBR,IN)=0.0
           HXEN1(I,IBR,IN)=0.0
           HXEN2(I,IBR,IN)=0.0
           HSAM1(I,IBR,IN)=0.0
           HSAM2(I,IBR,IN)=0.0
           HNEP1(I,IBR,IN)=0.0
           HNEP2(I,IBR,IN)=0.0
  160    CONTINUE
  170  CONTINUE
  180 CONTINUE
*---------------------------------------------------------------*
* READ AN OPTION KEY WORD
  185 CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('AFMDRV: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
* READ THE PRINT INDEX.
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('AFMDRV: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'REFT') THEN
         DO 190 IN=1,NUT
           CALL REDGET(INDIC,KTYP(IN),FLOTT,TEXT4,DFLOTT)
           IF(INDIC.NE.1) CALL XABORT('AFMDRV: INTEGER DATA EXPECTED.')
           CALL REDGET (INDIC,NITMA,FLOTT,TEXTR,DFLOTT)
           IF(INDIC.NE.3)
     1         CALL XABORT('AFMDRV: CHARACTER DATA EXPECTED.')
           READ(TEXTR,'(3A4)') (ITEXTR((IN-1)*3+I),I=1,3)
 190     CONTINUE
         IF(LMCR .AND. KTYP(1).GT.MX)
     +            CALL XABORT('AFMDRV: INVALID INDEX NUMBER.')
C
* CHECK THE NAME OF THE DIRECTORY
         WRITE(TEXTR,'(3A4)') (ITEXTR(I1),I1=1,3)
         CALL LCMLEN(IPFBM,TEXTR,ILENGT,ITYLCM)
         IF(ILENGT.EQ.0) THEN
           CALL XABORT('AFMDRV: UNABLE TO FIND '//TEXTR//' .')
         ENDIF
* RECOVER THE REFERENCE LOCAL PARAMETERS VALUES
         CALL LCMSIX(IPFBM,TEXTR,1)
         CALL LCMSIX(IPFBM,'INFO-NOMINA',1)
         CALL LCMLEN(IPFBM,'NOMINALP',ILP,ITYLCM)
         IF(ILP.GT.0) THEN
           CALL LCMGET(IPFBM,'NOMINALP',RLOC)
           CALL LCMGET(IPFBM,'NOMINALN',HISO)
           DO 200 I=1,ILP
             WRITE(HMICRO,'(3A4)') (HISO((I-1)*3+IH),IH=1,3)
             IF(HMICRO.EQ.'PW') PWREF=RLOC(I)
             IF(HMICRO.EQ.'TCOOL') TCR=RLOC(I)
             IF(HMICRO.EQ.'TMOD') TMR=RLOC(I)
             IF(HMICRO.EQ.'TFUEL') TFR=RLOC(I)
             IF(HMICRO.EQ.'RHOC') DCR=RLOC(I)
             IF(HMICRO.EQ.'RHOM') DMR=RLOC(I)
             IF(HMICRO.EQ.'PUR') XIR=RLOC(I)
 200       CONTINUE
         ENDIF
         CALL LCMSIX(IPFBM,' ',2)
         CALL LCMSIX(IPFBM,' ',2)
* REFERENCE PARAMETER VALUES
         PFIX=PWREF
         AW=15.9994 +2*(1-XIR)*1.0079 +2*XIR*2.014101
         PH=2*1.0079/AW
         PD=2*2.014101/AW
* INITIALISATION OF PERTURBED PARAMETER
         TF=TFR
         TC=TCR
         TM=TMR
         DC=1.0
         DM=1.0
         XI=XIR
         BOR=0.0
         SM=0.0
         RNP9=0.0
         XEN=0.0
*
         DO 210 IMX=1,MX
           POWER(IMX)=PWREF
           ISFT(IMX)=0
           BURBG(IMX)=0.0
           BURED(IMX)=0.0
           VOL(IMX)=0.0
           PDCOOL(IMX)=DCR
           PTCOOL(IMX)=TCR
           PTFUEL(IMX)=TFR
  210    CONTINUE
*        RECOVER THE TEMERATURE AND DENSITY PROFILES
         IF( (.NOT.LMCR).AND.(NPARM.GT.0) ) THEN
            JPMAP=LCMGID(IPMAP,'PARAM')
            DO 220 IPARM=1,NPARM
            KPMAP=LCMGIL(JPMAP,IPARM)
            CALL LCMGTC(KPMAP,'P-NAME',12,1,PNAME)
            CALL LCMGTC(KPMAP,'PARKEY',12,1,PARKEY)
            CALL LCMGET(KPMAP,'P-TYPE',IPTYPE)
            IF(IPTYPE.EQ.1) THEN
               CALL LCMGET(KPMAP,'P-VALUE',FLOTT)
            ELSE IF(IPTYPE.EQ.2) THEN
               CALL LCMLEN(KPMAP,'P-VALUE',NITMA,ITYLCM)
               IF(NITMA.NE.MX) CALL XABORT('@AFMDRV: INVALID LENGTH FO'
     1         //'R P-VALUE.')
            ENDIF
            IF(PNAME.EQ.'T-COOL') THEN
               WRITE(6,716) PNAME,PARKEY
               IF(IPTYPE.EQ.1) THEN
                  CALL XDRSET(PTCOOL,MX,FLOTT)
               ELSE IF(IPTYPE.EQ.2) THEN
                  CALL LCMGET(KPMAP,'P-VALUE',PTCOOL)
               ENDIF
            ELSE IF(PNAME.EQ.'D-COOL') THEN
               WRITE(6,716) PNAME,PARKEY
               IF(IPTYPE.EQ.1) THEN
                  CALL XDRSET(PDCOOL,MX,FLOTT)
               ELSE IF(IPTYPE.EQ.2) THEN
                  CALL LCMGET(KPMAP,'P-VALUE',PDCOOL)
               ENDIF
            ELSE IF(PNAME.EQ.'T-FUEL') THEN
               WRITE(6,716) PNAME,PARKEY
               IF(IPTYPE.EQ.1) THEN
                  CALL XDRSET(PTFUEL,MX,FLOTT)
               ELSE IF(IPTYPE.EQ.2) THEN
                  CALL LCMGET(KPMAP,'P-VALUE',PTFUEL)
               ENDIF
            ENDIF
  220       CONTINUE
         ENDIF
*
         CALL XDRSET(PW,MNPS,PWREF)
         CALL XDRSET(BRH,MNPS,0.0)
         CALL XDRSET(POWER,MX,PWREF)
*
      ELSE IF(TEXT4.EQ.'TFUEL') THEN
        CALL REDGET (INDIC,NITMA,TFU,TEXT4,DFLOTT)
        LTFUEL = .TRUE.
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'TCOOL') THEN
         CALL REDGET (INDIC,NITMA,TCU,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
          LTCOOL = .TRUE.
          CALL XDRSET(PTCOOL,MX,TCU)
*
      ELSE IF(TEXT4.EQ.'TMOD') THEN
        CALL REDGET (INDIC,NITMA,TM,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'RDCL') THEN
        CALL REDGET (INDIC,NITMA,DCU,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
          LDCOOL = .TRUE.
          CALL XDRSET(PDCOOL,MX,DCU)
*
      ELSE IF(TEXT4.EQ.'RDMD') THEN
        CALL REDGET (INDIC,NITMA,DM,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
        DM=DM/DMR
*
      ELSE IF(TEXT4.EQ.'BORON') THEN
        CALL REDGET (INDIC,NITMA,BOR,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
*  ppm eq 10**-6, NO CONSISTENCY WITH CFC CONCENTRATIONS
*  NEED TO ADD A COEFFICIENT TO FIT THE DATA (BREF should be 0.0ppm)
*
        BOR=BOR*1.E-6
*
      ELSE IF(TEXT4.EQ.'PUR') THEN
        CALL REDGET (INDIC,NITMA,XI,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
        XI=XI*1.0E-02
*
      ELSE IF(TEXT4.EQ.'FIXP') THEN
        CALL REDGET (INDIC,NITMA,PFIX,TEXT4,DFLOTT)
        IF(INDIC.EQ.2) THEN
          LNOMP=.TRUE.
        ELSE IF(TEXT4.EQ.'INIT') THEN
          LINI=.TRUE.
        ELSE
          CALL XABORT('AFMDRV: "INIT" or REAL DATA EXPECTED.')
        ENDIF
*
      ELSE IF(TEXT4.EQ.'IMET') THEN
        CALL REDGET(INDIC,ITM,FLOTT,TEXT4,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('AFMDRV: INTEGER DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'XENON') THEN
        LXENON=.TRUE.
        CALL REDGET (INDIC,NITMA,FXEN,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'XEREF') THEN
        LXEREF=.TRUE.
*
      ELSE IF(TEXT4.EQ.'DRAH') THEN
         LDRAH=.TRUE.
*
      ELSE IF(TEXT4.EQ.'SAM') THEN
        LSAM=.TRUE.
        CALL REDGET (INDIC,NITMA,FSAM,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'NEP') THEN
        LNEP=.TRUE.
        CALL REDGET (INDIC,NITMA,FNEP,TEXT4,DFLOTT)
        IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
*
      ELSE IF(TEXT4.EQ.'NREF') THEN
        LNEREF=.TRUE.
*
      ELSE IF(TEXT4.EQ.'BURN') THEN
        IF(LMCR) THEN
          CALL REDGET (INDIC,NITMA,FBUR,TEXT4,DFLOTT)
          IF(INDIC.NE.2) CALL XABORT('AFMDRV: REAL DATA EXPECTED.')
        ELSE
          CALL XABORT('AFMDRV: INVALID KEYWORD BURN.')
        ENDIF
*
      ELSE IF(TEXT4.EQ.'NPWF') THEN
         LPWF=.FALSE.
      ELSE IF(TEXT4.EQ.'PWF') THEN
         LPWF=.TRUE.
      ELSE IF(TEXT4.EQ.'BLIN') THEN
         ILIN=1
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 230
      ELSE
        CALL XABORT('AFMDRV: '//TEXT4//' IS AN INVALID KEY-WORD.')
      ENDIF
      GO TO 185
*     EQUIVALENT MODERATOR DENSITY FOR THE REFERENCE PURITY
  230 DXI = XI - XIR
* pas de modification de densite selon la purete D2O
*      DM=DM/(1.0+DXI*(PD-PH))
*---------------------------------------------------------------*
* RECOVER NEUTRONICS PARAMETRES
      WRITE(TEXTR,'(3A4)') (ITEXTR(I1),I1=1,3)
      CALL LCMSIX(IPFBM,TEXTR,1)
      CALL LCMGET(IPFBM,'VOLUME',VOL(1))
      CALL LCMGET(IPFBM,'ENERGY',ENER)
      CALL LCMGET(IPFBM,'HITAB',HISO)
      CALL LCMGET(IPFBM,'JTAB',JTAB)
      CALL LCMSIX(IPFBM,' ',2)
      DO 280 IN=1,NUT
        WRITE(TEXTR,'(3A4)') (ITEXTR((IN-1)*3+I1),I1=1,3)
        CALL LCMSIX(IPFBM,TEXTR,1)
        CALL LCMGET(IPFBM,'BURNUP',XBURN(1,IN))
*     RECOVER THE EXISTING DATABASE.
* RECOVER THE HISTORY COEFFICIENTS
        DO 270 I = 1,NBURN
          WRITE(TEXTB,'(4HBURN,4X,I4)') I
          CALL LCMSIX(IPFBM,TEXTB,1)
*
          IF(JTAB(1).EQ.1) THEN
            CALL LCMSIX(IPFBM,'HISTORY',1)
            CALL LCMGET(IPFBM,'PHIL1',CPW1B(1,I,IN))
            CALL LCMGET(IPFBM,'PHIS1',CPW1B(2,I,IN))
            CALL LCMGET(IPFBM,'PHIL2',CPW2B(1,I,IN))
            CALL LCMGET(IPFBM,'PHIS2',CPW2B(2,I,IN))
            CALL LCMLEN(IPFBM,'PHISX1',IHISTO,ITYLCM)
            IF(IHISTO.GT.0) THEN
              CALL LCMGET(IPFBM,'PHILX1',HXEN1(1,I,IN))
              CALL LCMGET(IPFBM,'PHISX1',HXEN1(2,I,IN))
              CALL LCMGET(IPFBM,'PHILX2',HXEN2(1,I,IN))
              CALL LCMGET(IPFBM,'PHISX2',HXEN2(2,I,IN))
C
              CALL LCMGET(IPFBM,'PHILS1',HSAM1(1,I,IN))
              CALL LCMGET(IPFBM,'PHISS1',HSAM1(2,I,IN))
              CALL LCMGET(IPFBM,'PHILS2',HSAM2(1,I,IN))
              CALL LCMGET(IPFBM,'PHISS2',HSAM2(2,I,IN))
C
              CALL LCMGET(IPFBM,'PHILN1',HNEP1(1,I,IN))
              CALL LCMGET(IPFBM,'PHISN1',HNEP1(2,I,IN))
              CALL LCMGET(IPFBM,'PHILN2',HNEP2(1,I,IN))
              CALL LCMGET(IPFBM,'PHISN2',HNEP2(2,I,IN))
            ENDIF
            CALL LCMSIX(IPFBM,' ',2)
          ENDIF
*
          CALL LCMGET(IPFBM,'FLUX-INTG',FLUXB(1,I,IN))
          CALL LCMGET(IPFBM,'OVERV',OVERVB(1,I,IN))
          CALL LCMGET(IPFBM,'ISOTOPESDENS',DENSITB(1,I,IN))
* COMPUTE DELTA-CONCENTRATION
          DO 250 ISO=1,NISO
            WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
            CALL LCMSIX(IPFBM,HMICRO,1)
            IF(JTAB(1).EQ.1) THEN
              IF((HMICRO.EQ.'XE135').OR.(HMICRO.EQ.'Xe135')) IXENO=ISO
              IF((HMICRO.EQ.'SM149').OR.(HMICRO.EQ.'Sm149')) ISAMA=ISO
              IF((HMICRO.EQ.'NP239').OR.(HMICRO.EQ.'Np239')) INEPT=ISO
              IF(HMICRO.EQ.'MACR ')
     1             CALL LCMGET(IPFBM,'CHI',CHIB(1,I,IN))
            ENDIF
* RECOVER MACROSCOPIC X-SECTIONS
            NTM=4+2*IXYZ
            DO 240 ITY=1,NTM
              IF(ITY.EQ.1) THEN
                IF(IXYZ.EQ.0)  THEN
                  TEXTD = 'STRD'
                ELSE IF(IXYZ.EQ.1) THEN
                  TEXTD = 'STRD X'
                ENDIF
              ENDIF
              IF(ITY.EQ.2) TEXTD = 'ABS'
              IF(ITY.EQ.3) TEXTD = 'NUSIGF'
              IF(ITY.EQ.4) TEXTD = 'H-FACTORS'
              IF(ITY.EQ.5) TEXTD = 'STRD Y'
              IF(ITY.EQ.6) TEXTD = 'STRD Z'
              CALL LCMLEN(IPFBM,TEXTD,ILENG,ITYXSM)
*
              IF(ILENG.NE.0) THEN
                CALL LCMSIX(IPFBM,TEXTD,1)
                CALL LCMGET(IPFBM,'REF',SMACB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'BOR',XBORB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'PUR',XPURB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T1M',XT1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T2M',XT2MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D1M',XD1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D2M',XD2MB(1,ITY,ISO,I,IN))
                IF(JTAB(1).EQ.1) THEN
                  CALL LCMLEN(IPFBM,'XEN',ILENGX,ITYXSM)
                  IF(ILENGX.GT.0)
     +               CALL LCMGET(IPFBM,'XEN',XXENB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T1F',XT1FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2F',XT2FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T1C',XT1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2C',XT2CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D1C',XD1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D2C',XD2CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'SM149',XSMB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'NP239',XNP9B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXFD',XMFDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXMD',XMMDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH1',XPF1B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL1',XPF1LB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH2',XPF2B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL2',XPF2LB(1,ITY,ISO,I,IN))
                ENDIF
*
                CALL LCMSIX(IPFBM,' ',2)
              ENDIF
 240        CONTINUE
*
            CALL LCMLEN(IPFBM,'NFTOT',ILNF,ITYXSM)
            IF(ILNF.NE.0) THEN
              CALL LCMGET(IPFBM,'NFTOT',SMACB(1,NTM+1,ISO,I,IN))
            ENDIF
            CALL LCMSIX(IPFBM,' ',2)
 250      CONTINUE
*     SCATTERING CROSS-SECTIONS
          DO 260 ISO=1,NISO
            WRITE(HMICRO,'(3A4)') (HISO((ISO-1)*3+IH),IH=1,3)
            CALL LCMLEN(IPFBM,HMICRO,ILENG,ITYLCM)
            IF(ILENG.EQ.0) GO TO 230
            CALL LCMSIX(IPFBM,HMICRO,1)
C            DO 150 IL=1,NL
              IL=1
              ITY=NTM+1+IL
              LTST=0
              WRITE (CM,'(I2.2)') IL-1
              CALL LCMLEN(IPFBM,'SCAT'//CM,ILENG,ITYXSM)
              IF(ILENG.NE.0) THEN
                LTST=1
              ELSE
                WRITE (CM,'(I2)') IL-1
                CALL LCMLEN(IPFBM,'SCAT'//CM,ILENG,ITYXSM)
                IF(ILENG.NE.0) THEN
                  LTST=2
                ENDIF
              ENDIF
              IF (LTST.GE.1) THEN
                CALL LCMSIX(IPFBM,'SCAT'//CM,1)
                IF(HMICRO.EQ.'MACR') THEN
                  IF (LTST.EQ.1) THEN
                    CALL LCMGET(IPFBM,'NJJS',NJ)
                    CALL LCMGET(IPFBM,'IJJS',IJ)
                  ELSE
                    CALL LCMGET(IPFBM,'NJJ',NJ)
                    CALL LCMGET(IPFBM,'IJJ',IJ)
                  ENDIF
                ENDIF
                CALL LCMGET(IPFBM,'REF',SMACB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'BOR',XBORB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'PUR',XPURB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T1M',XT1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'T2M',XT2MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D1M',XD1MB(1,ITY,ISO,I,IN))
                CALL LCMGET(IPFBM,'D2M',XD2MB(1,ITY,ISO,I,IN))
                IF(JTAB(1).EQ.1) THEN
                  CALL LCMLEN(IPFBM,'XEN',ILENG,ITYXSM)
                  IF(ILENG.GT.0) THEN
                    CALL LCMGET(IPFBM,'XEN',XXENB(1,ITY,ISO,I,IN))
                  ENDIF
                  CALL LCMGET(IPFBM,'T1F',XT1FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2F',XT2FB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T1C',XT1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'T2C',XT2CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D1C',XD1CB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'D2C',XD2CB(1,ITY,ISO,I,IN))
*
                  CALL LCMGET(IPFBM,'SM149',XSMB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'NP239',XNP9B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXFD',XMFDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'MIXMD',XMMDB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH1',XPF1B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL1',XPF1LB(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCH2',XPF2B(1,ITY,ISO,I,IN))
                  CALL LCMGET(IPFBM,'FPCL2',XPF2LB(1,ITY,ISO,I,IN))
                ENDIF
*
                CALL LCMSIX(IPFBM,' ',2)
              ENDIF
            CALL LCMSIX(IPFBM,' ',2)
  260     CONTINUE
C
          CALL LCMSIX(IPFBM,' ',2)
  270   CONTINUE
        CALL LCMSIX(IPFBM,' ',2)
  280 CONTINUE
      IF(JTAB(1).EQ.1) THEN
        IF(IXENO.EQ.0) CALL XABORT('NO XE135 FOUND ')
        IF(ISAMA.EQ.0) CALL XABORT('NO SM149 FOUND ')
        IF(INEPT.EQ.0) CALL XABORT('NO NP239 FOUND ')
      ENDIF
*     END OF THE RECOVERING PROCESS
*---------------------------------------------------------------*
*        ISC INDICATE THE TYPE OF CROSS-SECTION CALCULATION
*        ISC=1 ; TIME AVERAGE CALCULATION
*        ISC=2 ; INSTANTANEOUS CALCULATION
*        ISC=3 ; HOMOGENEOUS CALCULATION
*---------------------------------------------------------------*
*
      IF(ISC.EQ.0) THEN
        CALL XABORT('AFMDRV: TIMAV/INSTANT BURNUP TREATMENT NOT SET')
      ELSE IF(ISC.EQ.1) THEN
*       Time-averaged calculation
        WRITE(6,699)
        MMIX=NBCH*NCCO
        LTAV=.TRUE.
        CALL LCMGET(IPMAP,'FLMIX',INDEX)
        CALL LCMGET(IPMAP,'BURN-BEG',BURBG)
        CALL LCMGET(IPMAP,'BURN-END',BURED)
        CALL LCMLEN(IPMAP,'BUND-PW',ILPW,ITYLCM)
          IF((ILPW.NE.0).AND.LPWF) THEN
            IF(IMPX.GE.1) WRITE(6,702)
            IF(.NOT.LINI) THEN
              CALL LCMGET(IPMAP,'BUND-PW',POWER)
            ELSE 
              CALL LCMLEN(IPMAP,'BUND-PW-INI',ILPW,ITYLCM)
              IF(ILPW.NE.0) THEN
                CALL LCMGET(IPMAP,'BUND-PW-INI',POWER)
              ELSE
                CALL XABORT('AFMDRV: NO INITIAL POWER IN L_MAP')
              ENDIF
            ENDIF
          ELSE
            CALL XDRSET(POWER,MMIX,PWREF)
          ENDIF
        CALL LCMLEN(IPMAP,'FLUX-AV',ILBFLU,ITYLCM)
        IF(ILBFLU.NE.0) THEN
          IF(IMPX.GE.1) WRITE(6,703)
          CALL LCMGET(IPMAP,'FLUX-AV',WORK)
          DO 300 IGR=1,NGRP
           DO 290 IBF=1,MMIX
             IIBF=MMIX*(IGR-1)+IBF
             BFLUX(IGR,IBF)=WORK(IIBF)
 290      CONTINUE
 300      CONTINUE
        ENDIF
      ELSE IF(ISC.EQ.2) THEN
*       Instantaneous calculation
        IF(LMCR) THEN
          MMIX=NBCH*NCCO
          CALL XDRSET(POWER,MMIX,PWREF)
        ELSE
          WRITE(6,701)
          MMIX=NBCH*NCCO
          CALL LCMGET(IPMAP,'FLMIX',INDEX)
          CALL LCMGET(IPMAP,'BURN-INST',BURBG)
          CALL LCMLEN(IPMAP,'BUND-PW',ILPW,ITYLCM)
          IF((ILPW.NE.0).AND.LPWF) THEN
            IF(IMPX.GE.1) WRITE(6,702)
            IF(.NOT.LINI) THEN
              CALL LCMGET(IPMAP,'BUND-PW',POWER)
            ELSE 
              CALL LCMLEN(IPMAP,'BUND-PW-INI',ILPW,ITYLCM)
              IF(ILPW.NE.0) THEN
                CALL LCMGET(IPMAP,'BUND-PW-INI',POWER)
              ELSE
                CALL XABORT('AFMDRV: NO INITIAL POWER IN L_MAP')
              ENDIF
            ENDIF
          ELSE
            CALL XDRSET(POWER,MMIX,PWREF)
          ENDIF
          CALL LCMLEN(IPMAP,'FLUX-AV',ILBFLU,ITYLCM)
          IF(ILBFLU.NE.0) THEN
            IF(IMPX.GE.1) WRITE(6,703)
            CALL LCMGET(IPMAP,'FLUX-AV',WORK)
            DO 320 IGR=1,NGRP
             DO 310 IBF=1,MMIX
               IIBF=MMIX*(IGR-1)+IBF
               BFLUX(IGR,IBF)=WORK(IIBF)
 310        CONTINUE
 320        CONTINUE
          ENDIF
* RECOVER THE SHIFT INFORMATION
          IF(MNPS.GT.2) THEN
            IF(IMPX.GE.1) WRITE(6,704)
            CALL LCMGET(IPMAP,'ISHIFT',ISFT)
            DO 330 IS=1,MNPS-2
              WRITE (CM,'(I2)') IS
              CALL LCMGET(IPMAP,'BSHIFT'//CM,BSFT(1,IS))
              CALL LCMGET(IPMAP,'PSHIFT'//CM,PSFT(1,IS))
 330        CONTINUE
          ENDIF
        ENDIF
      ELSE IF(ISC.EQ.3) THEN
*       Homogeneous calculation
        MMIX=NCZO
        LTAV=.TRUE.
        CALL LCMGET(IPMAP,'B-ZONE',IZONE)
        CALL LCMGET(IPMAP,'FLMIX',INDEX)
        CALL LCMGET(IPMAP,'BURN-AVG',BURED)
      ENDIF
*---------------------------------------------------------------*
      IF(IMPX.GE.1) THEN
         IF(LNOMP) WRITE(6,705) PFIX
         IF(LXENON) WRITE(6,706) FXEN
         IF(LSAM) WRITE(6,719) FSAM
         IF(LNEP) WRITE(6,711) FNEP
         IF(LXEREF) WRITE(6,712)
         IF(LNEREF) WRITE(6,713)
         IF(LTFUEL) WRITE(6,714) TFU
         IF(IHISTO.GT.0.AND.LDRAH) WRITE(6,715)
         IF(LTCOOL) WRITE(6,717) TCU
         IF(LDCOOL) WRITE(6,718) DCU
      ENDIF
*---------------------------------------------------------------*
* MIXTURE SHIFT
      IF(LMCR) THEN
        MXSH=MMIX
        CALL XDRSET(VOL,MMIX,VOL(1))
      ELSE
        MXSH=1
      ENDIF
*---------------------------------------------------------------*
*        LOOP OVER THE MIXTURES
      DO 540 NMIX=MXSH,MMIX
        TC=PTCOOL(NMIX)
        DC=PDCOOL(NMIX)/DCR
        IF(LMCR) THEN
          NPS=2
          IDF=1
        ELSE
          VOL(NMIX)=VOL(1)
          NPS=ISFT(NMIX)+2
          KDF=0
          DO 340 IN=1,NUT
            IF(INDEX(NMIX).EQ.KTYP(IN)) THEN
              IDF=IN
              KDF=1
            ENDIF
  340     CONTINUE
          IF(KDF.EQ.0) CALL XABORT('AFMDRV: WRONG NUMBER OF INDEX')
        ENDIF
* IF TIME AVERAGE CALCULATION:
* EVALUATION OF THE BURNUPS STEPS EMBEDED IN THE INTEGRATION
        IF(LTAV) THEN
          XBMIN=BURBG(NMIX)
          XBMAX=BURED(NMIX)
* TIME AVERAGE BURNUP LOCALISATION
          CALL AFMLOC(NBURN,NTP,XBMAX,XBMIN,XBURN(1,IDF),
     1                IMAX,IMIN,XCOF,ILIN)
*       LAGRANGE METHOD (TIME-AVERAGE)
          IMINR=IMIN
          IMAXR=ABS(IMAX)
*       SPLINE OR HERMITE METHOD (TIME-AVERAGE)
          IF(ITM.EQ.2.OR.ITM.EQ.3) THEN
            IMINR=1
            IMAXR=NBURN
          ENDIF
*
        ELSE
          IMINR=1
          IMAXR=1
        ENDIF
C
        DO 450 JR=IMINR,IMAXR
         IF(LTAV) THEN
           IRAV=JR
           NPS=2
         ELSE
           IF(NPS.GT.2) THEN
             DO 350 K=2,NPS-1
               IS=K-1
               BRH(K)=BSFT(NMIX,IS)
 350         CONTINUE
           ENDIF
           IF(LMCR) THEN
             BRH(NPS)=FBUR
             IF(JTAB(1).EQ.0) BRH(NPS)=0.0
           ELSE
             BRH(NPS)=BURBG(NMIX)
           ENDIF
         ENDIF
*
         IF(LNOMP) THEN
           DO 360 K=2,NPS
             PW(K)=PFIX
 360       CONTINUE
         ELSE
           IF(NPS.GT.2) THEN
             DO 370 K=2,NPS-1
               IS=K-1
               PW(K)=PSFT(NMIX,IS)
 370         CONTINUE
           ENDIF
           PW(NPS)=POWER(NMIX)
         ENDIF
*        D. Rozon 'Introduction a la Cinetique des Reacteur Nucleaires'
*        Edition E.P., 1992. (p.217) or 1998 (p.185)
*        PW is assumed to be in kW.
         IF(IPROF2.GT.0) THEN
           TF = PTFUEL(NMIX)
         ELSE
           TF= TC + 0.476*PW(NPS) + 2.267*PW(NPS)*PW(NPS)*1.0E-04
         ENDIF
C INITIAL CONCENTRATIONS
         ZXREF=0.0
         SM=0.0
         ZRNP9=0.0
*      IF FUEL
         IF(JTAB(1).EQ.1) THEN
* BURNUP LOCALISATION FOR XENON AND FISSION X-SECTION INTERPOLATION
           IF(LTAV) THEN
             XIFL=XBURN(IRAV,IDF)
             IMAXX=IRAV
             IMINX=IRAV
             XCOF(1)=1.0D0
             XCOF(2)=0.0D0
             XCOF(3)=0.0D0
           ELSE
             XIFL=BRH(NPS)
             CALL AFMLOC(NBURN,NTP,BRH(NPS),BRH(NPS),XBURN(1,IDF),
     1                 IMAXX,IMINX,XCOF,ILIN)
           ENDIF
*
           DO 380 IGR = 1,NGRP
             XSIGX(IGR)=0.0
             XFLUN(IGR)=0.0
             XSIGF(IGR)=0.0
  380       CONTINUE
*       INTERPOLATION OF THE CONCENTRATION
*
           IIX=0
           DO 395 I = IMINX,IMAXX
             IIX=IIX+1
             RXCOF=REAL(XCOF(IIX))
             ZXREF=DENSITB(IXENO,I,IDF)*RXCOF  +ZXREF
             XEN=ZXREF
             SM=DENSITB(ISAMA,I,IDF)*RXCOF     +SM
             ZRNP9=DENSITB(INEPT,I,IDF)*RXCOF  +ZRNP9
             RNP9=ZRNP9
*
             DO 390 IGR=1,NGRP
               XSIGX(IGR)=SMACB(IGR,2,IXENO,I,IDF)*RXCOF
     1                    + XSIGX(IGR)
               XFLUN(IGR)=FLUXB(IGR,I,IDF)*RXCOF     + XFLUN(IGR)
               XSIGF(IGR)=SMACB(IGR,5,1,I,IDF)*RXCOF + XSIGF(IGR)
  390        CONTINUE
  395       CONTINUE
           IF(LDRAH.AND.IHISTO.GT.0) THEN
             IF(PW(NPS).GT.PWREF) THEN
               XPW=ALOG(PW(NPS)/PW(1))
               XPWM=1.0/PW(NPS)-1.0/PW(1)
               IFH=1
             ELSE
               XPW=PW(NPS)-PW(1)
               XPWM=(PW(NPS)-PW(1))**2
               IFH=2
             ENDIF
C
             XEN  =ZXREF
             RNP9 =ZRNP9
             IIX=0
             DO 400 I = IMINX,IMAXX
               IIX=IIX+1
               RXCOF=REAL(XCOF(IIX))
*            COMPUTE XENON-SAMRIUM-NEPTUNIUM CONCENTRATION USING DRAGON
               XEN  =XEN   +HXEN1(IFH,I,IDF)*XPW*RXCOF+
     1                      HXEN2(IFH,I,IDF)*XPWM*RXCOF
               SM  =SM     +HSAM1(IFH,I,IDF)*XPW*RXCOF+
     1                      HSAM2(IFH,I,IDF)*XPWM*RXCOF
               RNP9 =RNP9  +HNEP1(IFH,I,IDF)*XPW*RXCOF+
     1                      HNEP2(IFH,I,IDF)*XPWM*RXCOF
 400         CONTINUE
           ELSE IF(ILBFLU.NE.0.AND.XIFL.NE.0.0) THEN
*          COMPUTE THE XENON AND NEPTUNIUM CONCENTRATIONS
              CALL AFMXNC(NGRP,XSIGX,XSIGF,BFLUX(1,NMIX),
     1                   XEN,RNP9,XFLUN)
           ENDIF
* COMPUTE THE XENON AND NEPTUNIUM CONCENTRATIONS
           IF(LXENON) XEN=FXEN
           IF(LSAM) SM=FSAM
           IF(LNEP) RNP9=FNEP
           IF(LXEREF) XEN=ZXREF
           IF(LNEREF) RNP9=ZRNP9
           IF(LTFUEL) THEN
!       fuel temperature as input
             TF=TFU
!       reference fuel temperature 
           ELSEIF(LMCR) THEN
             TF=TFR
           ENDIF
         ENDIF
*---------------------------------------------------------------*
* XSECTION CALCULATION
*---------------------------------------------------------------*
         CALL AFMCPT(KENTRY,NBURN,NGRP,NISO,
     1   NL,IMPX,SMACB,XBORB,XPURB,XXENB,XT1FB,XT2FB,XT1CB,
     1   XT2CB,XT1MB,XT2MB,XD1CB,XD2CB,XD1MB,XD2MB,
     1   XSMB,XNP9B,XMFDB,XMMDB,XPF1B,XPF2B,XPF1LB,XPF2LB,
     1   DENSITB,CPW1B,CPW2B,FLUXB,OVERVB,CHIB,
     1   IJ,NJ,HISO,CTITRE,
     1   NMIX,SIGMA,NTYP,TF,TC,TM,DC,DM,BOR,XEN,SM,RNP9,XI,
     1   TFR,TCR,TMR,XIR,OVERV,FLUX,CHI,SCAT,MX,NPS,PW,BRH,
     1   XBURN,LTAV,IRAV,IDF,JTAB,IXYZ,ILIN)
*---------------------------------------------------------------*
*
         DO 420 IGR=1,NGRP
          FLUAV(JR,IGR)=FLUX(NMIX,IGR)
          DO 410 ITY=1,NTM+1
            SIGAV(JR,IGR,ITY)=SIGMA(NMIX,IGR,ITY)
 410     CONTINUE
 420     CONTINUE
         IL =1
         DO 440 IGR=1,NGRP
            DO 430 JGR=1,NGRP
              SCATAV(JR,IL,JGR,IGR)=SCAT(NMIX,IL,JGR,IGR)
 430        CONTINUE
 440     CONTINUE
 450    CONTINUE
        IF(LTAV) THEN
* COMPUTE  TIME AVERAGED X-SECTIONS
          DO 470 IGR=1,NGRP
            CALL AFMTAV(NBURN,ITM,XBMAX,XBMIN,FLUAV(1,IGR),IMIN,IMAX,
     1      XBURN,FLUX(NMIX,IGR))
            DO 460 ITY=1,NTM+1
              CALL AFMTAV(NBURN,ITM,XBMAX,XBMIN,SIGAV(1,IGR,ITY),IMIN,
     1        IMAX,XBURN,SIGMA(NMIX,IGR,ITY))
 460      CONTINUE
 470      CONTINUE
*
          DO 490 IGR=1,NGRP
           DO 480 JGR=1,NGRP
               IL=1
               CALL AFMTAV(NBURN,ITM,XBMAX,XBMIN,SCATAV(1,IL,IGR,JGR),
     1         IMIN,IMAX,XBURN,SCAT(NMIX,IL,IGR,JGR))
 480      CONTINUE
 490      CONTINUE
*
        ENDIF
* COMPUTE DIRECTIONAL DIFFUSION COEFFICIENTS FROM STRD
*  X-SECTIONS.
        IF(IXYZ.EQ.0) THEN
          DO 500 IGR=1,NGRP
            DIFFX(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,1))
 500      CONTINUE
          ILEAK=1
        ELSE IF(IXYZ.EQ.1) THEN
          DO 510 IGR=1,NGRP
            DIFFX(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,1))
            DIFFY(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,5))
            DIFFZ(NMIX,IGR)=1.0/(3.0*SIGMA(NMIX,IGR,6))
            ILEAK=2
 510      CONTINUE
        ENDIF
*
        IL=1
        DO 530 IGR=1,NGRP
            NJJ(NMIX,IL,IGR)=NJ(IGR)
            IJJ(NMIX,IL,IGR)=IJ(IGR)
            IF(LMCR) THEN
              DO 520 NI=1,MMIX
                NJJ(NI,IL,IGR)=NJ(IGR)
                IJJ(NI,IL,IGR)=IJ(IGR)
 520          CONTINUE
            ENDIF
 530    CONTINUE
* MIX LOOP
 540  CONTINUE
*
      IF(LTAV) THEN
        IF(IMPX.GE.1.AND.ITM.EQ.1) WRITE(6,707)
        IF(IMPX.GE.1.AND.ITM.EQ.2) WRITE(6,708)
        IF(IMPX.GE.1.AND.ITM.EQ.3) WRITE(6,709)
      ENDIF
*---------------------------------------------------------------*
*        DECOMPRESS BURN ZONE FOR ALL THE BUNDLES
      IF(ISC.EQ.3) THEN
       MMIX=NBCH*NCCO
       DO 870 IGR=1,NGRP
        DO 550 IZ=1,NCZO
            WORK(IZ)=DIFFX(IZ,IGR)
 550    CONTINUE
        DO 570 IC=1,NCCO
          DO 560 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            DIFFX(ICB,IGR)=WORK(IZONE(IC))
 560      CONTINUE
 570    CONTINUE
*
        IF(ILEAK.EQ.2) THEN
          DO 580 IZ=1,NCZO
            WORK(IZ)=DIFFY(IZ,IGR)
 580      CONTINUE
          DO 600 IC=1,NCCO
           DO 590 IB=1,NBCH
             ICB=NBCH*(IC-1)+IB
             DIFFY(ICB,IGR)=WORK(IZONE(IC))
 590       CONTINUE
 600      CONTINUE
*
          DO 610 IZ=1,NCZO
            WORK(IZ)=DIFFZ(IZ,IGR)
 610      CONTINUE
          DO 630 IC=1,NCCO
           DO 620 IB=1,NBCH
             ICB=NBCH*(IC-1)+IB
             DIFFZ(ICB,IGR)=WORK(IZONE(IC))
 620       CONTINUE
 630      CONTINUE
        ENDIF
*
        DO 670 ITY=2,NTM+1
         DO 640 IZ=1,NCZO
           WORK(IZ)=SIGMA(IZ,IGR,ITY)
 640     CONTINUE
         DO 660 IC=1,NCCO
          DO 650 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            SIGMA(ICB,IGR,ITY)=WORK(IZONE(IC))
 650      CONTINUE
 660     CONTINUE
 670    CONTINUE
*
        DO 680 IZ=1,NCZO
          WORK(IZ)=FLUX(IZ,IGR)
 680    CONTINUE
        DO 700 IC=1,NCCO
         DO 690 IB=1,NBCH
           ICB=NBCH*(IC-1)+IB
           FLUX(ICB,IGR)=WORK(IZONE(IC))
 690     CONTINUE
 700    CONTINUE
*
        DO 710 IZ=1,NCZO
          WORK(IZ)=OVERV(IZ,IGR)
 710    CONTINUE
        DO 730 IC=1,NCCO
         DO 720 IB=1,NBCH
           ICB=NBCH*(IC-1)+IB
           OVERV(ICB,IGR)=WORK(IZONE(IC))
 720     CONTINUE
 730    CONTINUE
*
        DO 740 IZ=1,NCZO
          WORK(IZ)=CHI(IZ,IGR)
 740    CONTINUE
        DO 760 IC=1,NCCO
         DO 750 IB=1,NBCH
           ICB=NBCH*(IC-1)+IB
           CHI(ICB,IGR)=WORK(IZONE(IC))
 750     CONTINUE
 760    CONTINUE
*
        IL=1
        DO 800 JGR=1,NGRP
           DO 770 IZ=1,NCZO
             WORK(IZ)=SCAT(IZ,IL,IGR,JGR)
 770       CONTINUE
           DO 790 IC=1,NCCO
            DO 780 IB=1,NBCH
              ICB=NBCH*(IC-1)+IB
              SCAT(ICB,IL,IGR,JGR)=WORK(IZONE(IC))
 780        CONTINUE
 790       CONTINUE
 800     CONTINUE
*
         DO 810 IZ=1,NCZO
           IWORK(IZ)=NJJ(IZ,IL,IGR)
 810     CONTINUE
         DO 830 IC=1,NCCO
          DO 820 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            NJJ(ICB,IL,IGR)=IWORK(IZONE(IC))
 820      CONTINUE
 830     CONTINUE
*
         DO 840 IZ=1,NCZO
           IWORK(IZ)=IJJ(IZ,IL,IGR)
 840     CONTINUE
         DO 860 IC=1,NCCO
          DO 850 IB=1,NBCH
            ICB=NBCH*(IC-1)+IB
            IJJ(ICB,IL,IGR)=IWORK(IZONE(IC))
 850      CONTINUE
 860     CONTINUE
*
 870   CONTINUE
*
       DO 880 IZ=1,NCZO
         WORK(IZ)=VOL(IZ)
 880   CONTINUE
       DO 900 IC=1,NCCO
        DO 890 IB=1,NBCH
          ICB=NBCH*(IC-1)+IB
          VOL(ICB)=WORK(IZONE(IC))
 890    CONTINUE
 900   CONTINUE
*
      ENDIF
*---
* STORE MACROLIB INFORMATIONS
*---
      IF(ITYPE.EQ.0)THEN
        CALL LCMPUT(IPMACX,'VOLUME',MMIX,2,VOL)
        CALL LCMPUT(IPMACX,'ENERGY',NGRP+1,2,ENER)
      ENDIF
*
      IF(LMCR) THEN
        STORE=VOL(MMIX)
        VOL(MMIX)= 0.0
*  MACROLIB EN MODIFICATION
        IF(ITYPE.NE.0) THEN
           CALL LCMGET(IPMACX,'VOLUME',VOL)
        ENDIF
        VOL(KTYP(1)) = STORE
        CALL LCMPUT(IPMACX,'VOLUME',MMIX,2,VOL)
        JPMAC=LCMLID(IPMACX,'GROUP',NGRP)
        DO 950 JGR=1,NGRP
          KPMAC=LCMDIL(JPMAC,JGR)
          STORE=SIGMA(MMIX,JGR,2)
          SIGMA(MMIX,JGR,2) = 0.0
*  MACROLIB EN MODIFICATION
          IF(ITYPE.NE.0) THEN
            CALL LCMGET(KPMAC,'NTOT0',SIGMA(1,JGR,2))
          ENDIF
          SIGMA(KTYP(1),JGR,2) = STORE
*
          STORE=OVERV(MMIX,JGR)
          OVERV(MMIX,JGR) = 0.0
*  MACROLIB EN MODIFICATION
          IF(ITYPE.NE.0) THEN
            CALL LCMGET(KPMAC,'OVERV',OVERV(1,JGR))
          ENDIF
          OVERV(KTYP(1),JGR) = STORE
*
          STORE=DIFFX(MMIX,JGR)
          DIFFX(MMIX,JGR) = 0.0
*  MACROLIB EN MODIFICATION
          IF(ITYPE.NE.0) THEN
             CALL LCMGET(KPMAC,'DIFFX',DIFFX(1,JGR))
          ENDIF
          DIFFX(KTYP(1),JGR) = STORE
*
          IF(ILEAK.EQ.2) THEN
            STORE=DIFFY(MMIX,JGR)
            DIFFY(MMIX,JGR) = 0.0
            IF(ITYPE.NE.0) THEN
               CALL LCMGET(KPMAC,'DIFFY',DIFFY(1,JGR))
            ENDIF
            DIFFY(KTYP(1),JGR) = STORE
*
            STORE=DIFFZ(MMIX,JGR)
            DIFFZ(MMIX,JGR) = 0.0
            IF(ITYPE.NE.0) THEN
               CALL LCMGET(KPMAC,'DIFFZ',DIFFZ(1,JGR))
            ENDIF
            DIFFZ(KTYP(1),JGR) = STORE
          ENDIF
*
          STORE = FLUX(MMIX,JGR)
          FLUX(MMIX,JGR) = 0.0
          IF(ITYPE.NE.0) THEN
             CALL LCMGET(KPMAC,'FLUX-INTG',FLUX(1,JGR))
          ENDIF
          FLUX(KTYP(1),JGR) = STORE
*
          IF(JTAB(1).EQ.1 .OR. ITYPE.NE.0) THEN
            STORE = CHI(MMIX,JGR)
            CHI(MMIX,JGR) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'CHI',CHI(1,JGR))
            ENDIF
            CHI(KTYP(1),JGR) = STORE
*
            STORE=SIGMA(MMIX,JGR,3)
            SIGMA(MMIX,JGR,3) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'NUSIGF',SIGMA(1,JGR,3))
            ENDIF
            SIGMA(KTYP(1),JGR,3) = STORE
*
            STORE=SIGMA(MMIX,JGR,5)
            SIGMA(MMIX,JGR,5) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'NFTOT',SIGMA(1,JGR,5))
            ENDIF
            SIGMA(KTYP(1),JGR,5) = STORE
*
            STORE=SIGMA(MMIX,JGR,4)
            SIGMA(MMIX,JGR,4) = 0.0
            IF(ITYPE.NE.0) THEN
              CALL LCMGET(KPMAC,'H-FACTOR',SIGMA(1,JGR,4))
            ENDIF
            SIGMA(KTYP(1),JGR,4) = STORE
*
          ENDIF
*
          IL=1
          ALLOCATE(SSCAT(NGRP))
          DO 910 IGR=1,NGRP
            SSCAT(IGR)= SCAT(MMIX,IL,IGR,JGR)
            SCAT(MMIX,IL,IGR,JGR) = 0.0
 910      CONTINUE
          IF(ITYPE.NE.0) THEN
!!  ATTENTION isotropy is supposed
!!
            IL=1
            WRITE (CM,'(I2.2)') IL-1
            CALL LCMGET(KPMAC,'SCAT'//CM,WORK)
            CALL LCMGET(KPMAC,'NJJS'//CM,NJJ(1,IL,JGR))
            CALL LCMGET(KPMAC,'IJJS'//CM,IJJ(1,IL,JGR))
            CALL LCMGET(KPMAC,'IPOS'//CM,IPOS)
            DO 930 IBM=1,MMIX
               IJJ0=IJJ(IBM,IL,JGR)
               IPOSDE = IPOS(IBM)
               DO 920 IGR=IJJ0,IJJ0-NJJ(IBM,IL,JGR)+1,-1
                SCAT(IBM,IL,IGR,JGR)=WORK(IPOSDE)
                IPOSDE=IPOSDE+1
 920           CONTINUE
 930        CONTINUE
          ENDIF
*
          DO 940 IGR=1,NGRP
            SCAT(KTYP(1),IL,IGR,JGR) = SSCAT(IGR)
 940      CONTINUE
          DEALLOCATE(SSCAT)
 950    CONTINUE
      ENDIF
*
      DO 990 IX=1,MMIX
       DO 980 JGR=1,NGRP
        DO 970 IL=1,NL
          IGMIN=JGR
          IGMAX=JGR
          DO 960 IGR=NGRP,1,-1
           IF (SCAT(IX,IL,IGR,JGR).NE.0.0) THEN
             IGMIN=MIN(IGMIN,IGR)
             IGMAX=MAX(IGMAX,IGR)
           ENDIF
  960     CONTINUE
          IJJ(IX,IL,JGR)=IGMAX
          NJJ(IX,IL,JGR)=IGMAX-IGMIN+1
  970   CONTINUE
  980  CONTINUE
  990 CONTINUE
*
      CALL XDRSET(SIGS,MMIX*NGRP,0.0)
      JPMAC=LCMLID(IPMACX,'GROUP',NGRP)
      DO 1002 JGR=1,NGRP
        KPMAC=LCMDIL(JPMAC,JGR)
        CALL LCMPUT(KPMAC,'NTOT0',MMIX,2,SIGMA(1,JGR,2))
        CALL LCMPUT(KPMAC,'OVERV',MMIX,2,OVERV(1,JGR))
        IF(ILEAK.EQ.1) THEN
          CALL LCMPUT(KPMAC,'DIFF',MMIX,2,DIFFX(1,JGR))
        ELSE IF(ILEAK.EQ.2) THEN
          CALL LCMPUT(KPMAC,'DIFFX',MMIX,2,DIFFX(1,JGR))
          CALL LCMPUT(KPMAC,'DIFFY',MMIX,2,DIFFY(1,JGR))
          CALL LCMPUT(KPMAC,'DIFFZ',MMIX,2,DIFFZ(1,JGR))
        ENDIF
        CALL LCMPUT(KPMAC,'FLUX-INTG',MMIX,2,FLUX(1,JGR))
        IF(JTAB(1).EQ.1 .OR. ITYPE.NE.0) THEN
          CALL LCMPUT(KPMAC,'CHI   ',MMIX,2,CHI(1,JGR))
          CALL LCMPUT(KPMAC,'NUSIGF   ',MMIX,2,SIGMA(1,JGR,3))
          CALL LCMPUT(KPMAC,'H-FACTOR',MMIX,2,SIGMA(1,JGR,4))
          CALL LCMPUT(KPMAC,'NFTOT',MMIX,2,SIGMA(1,JGR,5))
        ENDIF
*
        IL=1
        WRITE (CM,'(I2.2)') IL-1
        IPOSDE=0
        DO 1001 IX=1,MMIX
          IPOS(IX)=IPOSDE+1
          DO 1000 IGR=IJJ(IX,IL,JGR),IJJ(IX,IL,JGR)-NJJ(IX,IL,JGR)+1,-1
             IPOSDE=IPOSDE+1
             WORK(IPOSDE)=SCAT(IX,IL,IGR,JGR)
             SIGS(IX,IGR)=SIGS(IX,IGR)+ SCAT(IX,IL,IGR,JGR)
 1000   CONTINUE
 1001   CONTINUE
*
        CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC,'IPOS'//CM,MMIX,1,IPOS)
        CALL LCMPUT(KPMAC,'NJJS'//CM,MMIX,1,NJJ(1,IL,JGR))
        CALL LCMPUT(KPMAC,'IJJS'//CM,MMIX,1,IJJ(1,IL,JGR))
        CALL LCMPUT(KPMAC,'SIGW'//CM,MMIX,2,SCAT(1,IL,JGR,JGR))
 1002 CONTINUE
      DO 1003 JGR=1,NGRP
        KPMAC=LCMDIL(JPMAC,JGR)
        IL=1
        WRITE (CM,'(I2.2)') IL-1
        CALL LCMPUT(KPMAC,'SIGS'//CM,MMIX,2,SIGS(1,JGR))
 1003 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XPURB,XPF2LB,XPF1LB,XPF2B,XPF1B,XMMDB,XMFDB,XNP9B,
     1 XSMB,XD2MB,XD1MB,XD2CB,XD1CB,XT2MB,XT1MB,XT2CB,XT1CB,XT2FB,XT1FB,
     2 XXENB,XBORB,SMACB)
      DEALLOCATE(ITEXTR,PTFUEL,PTCOOL,PDCOOL,ISFT,PSFT,BSFT,XFLUN,XSIGX,
     1 XSIGF,KTYP,INDEX,OVERVB,CHIB,JTAB,FLUXB,CPW2B,CPW1B,HNEP2,HNEP1,
     2 HSAM2,HSAM1,HXEN2,HXEN1,HISO,DENSITB,BFLUX,NJ,BRH,PW,SCATAV,
     3 IWORK,SIGAV,FLUAV,POWER,BURED,BURBG,IZONE,WORK,IJ,DIFFZ,DIFFY,
     4 DIFFX,SCAT,IPOS,ENER,CHI,FLUX,SIGS,OVERV,XBURN,NJJ,VOL,IJJ,SIGMA)
      RETURN
*
  699 FORMAT(/' AFMDRV: THE CROSS SECTIONS ARE GENERATED FOR A',
     1 ' TIME AVERAGE CALCULATION.')
  701 FORMAT(/' AFMDRV: THE CROSS SECTIONS ARE GENERATED FOR A',
     1 ' SNAPSHOT CALCULATION.')
  702 FORMAT(/' AFMDRV: POWER ARE RECOVERED FROM L_MAP.')
  703 FORMAT(/' AFMDRV: FLUX  ARE RECOVERED FROM L_MAP.')
  704 FORMAT(/' AFMDRV: BUNDLES POWER SHIFT ARE CORRECTED.')
  705 FORMAT(/' AFMDRV: BUNDLES POWER = ',F12.2,1X,'KW IS FIXED',
     1 ' BY THE USER.')
  706 FORMAT(/' AFMDRV: BUNDLES XENON = ',E15.8,1X,'IS FIXED',
     1 ' BY THE USER.')
  707 FORMAT(/' AFMDRV: LAGRANGE INTERPOLATION IS USED TO COMPUTE',
     1 ' TIME AVERAGED CROSS SECTIONS.')
  708 FORMAT(/' AFMDRV: SPLINE 3 INTERPOLATION IS USED TO COMPUTE',
     1 ' TIME AVERAGED CROSS SECTIONS.')
  709 FORMAT(/' AFMDRV: HERMITE 3 INTERPOLATION IS USED TO COMPUT',
     1 'E TIME AVERAGED CROSS SECTIONS.')
  711 FORMAT(/' AFMDRV: BUNDLES NEPTUNIUM = ',E15.8,1X,'IS FIXED',
     1 ' BY THE USER.')
  712 FORMAT(/' AFMDRV: NOMINAL XENON IS USED.')
  713 FORMAT(/' AFMDRV: NOMINAL NEPTUNIUM IS USED.')
  714 FORMAT(/' AFMDRV: BUNDLES TFUEL = ',F12.2,1X,'K IS FIXED',
     1 ' BY THE USER.')
  715 FORMAT(/' AFMDRV: DRAGON CONCENTRATIONS ARE USED (XE135'
     1 //' NP239, SM149).')
  716 FORMAT(/' AFMDRV: ',A12,' PROFILES ARE RECOVERED FROM L_MAP.',
     1 ' PARKEY=',A12)
  717 FORMAT(/' AFMDRV: BUNDLES COOL. TEMP. TCOOL = ',F12.2,1X,
     1 'K IS FIXED BY THE USER.')
  718 FORMAT(/' AFMDRV: BUNDLES COOL. DENSITY RDCL = ',F12.9,1X,
     1 'K IS FIXED BY THE USER.')
  719 FORMAT(/' AFMDRV: BUNDLES SAMARIUM = ',E15.8,1X,'IS FIXED',
     1 ' BY THE USER.')
      END
