*DECK TONSPH
      SUBROUTINE TONSPH(IPLIB,IPTRK,IFTRAK,NREG,NUN,NBM,NBISO,ISONAM,
     1 MAT,VOL,KEYFLX,CDOOR,INRS,LEAKSW,IMPX,DEN,MIX,LSHI,ITRANC,
     2 IPHASE,NGRO,IGRMIN,IGRMAX,NBNRS,TITR,SIGT2,SIGT3,SN,SPH,ICPIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* SPH equivalence procedure over the self-shielded cross sections. Use
* all the standard solution doors of Dragon.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
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
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* NREG    number of regions.
* NUN     number of unknowns per energy group.
* NBM     number of mixtures in the internal library.
* NBISO   number of isotopes.
* ISONAM  alias name of isotopes in IPLIB.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* CDOOR   name of the geometry/solution operator.
* INRS    index of the resonant isotope under consideration.
* LEAKSW  leakage flag (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* IMPX    print flag (equal to zero for no print).
* DEN     density of each isotope.
* MIX     mix number of each isotope (can be zero).
* LSHI    resonant region number associated with each isotope.
*         Infinite dilution will be assumed if LSHI(i)=0.
* ITRANC  type of transport correction.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* NGRO    number of energy groups.
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* NBNRS   number of totally correlated fuel regions. NBNRS=max(IRES).
* TITR    title.
* SIGT2   total macroscopic cross sections.
* SIGT3   transport correction.
* SN      computed dilution cross section in each energy group of
*         each isotope.
*
*Parameters: output
* SPH     SPH factors.
* ICPIJ   number of flux solution door calls.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPTRK
      INTEGER IFTRAK,NREG,NUN,NBM,NBISO,ISONAM(3,NBISO),MAT(NREG),
     1 KEYFLX(NREG),INRS,IMPX,MIX(NBISO),LSHI(NBISO),ITRANC,IPHASE,
     2 NGRO,IGRMIN,IGRMAX,NBNRS,ICPIJ
      REAL VOL(NREG),DEN(NBISO),SIGT2(NBM,NGRO),SIGT3(NBM,NGRO),
     1 SN(NGRO,NBISO),SPH(NBM,NGRO)
      LOGICAL LEAKSW
      CHARACTER CDOOR*12,TITR*72,HNAMIS*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      TYPE(C_PTR) JPLIB,KPLIB,IPMACR,IPSOU
      LOGICAL LHOMOG,LPROB,LEXAC,LOGDO,REBFLG
      CHARACTER TEXT12*12
      INTEGER NALBP,ISTATE(NSTATE)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRES,ISONR,NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: VST,SIGTXS,SIGS0X,FLNEW
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLUX,SIG0,SIG1,SIG3,TOTAL,
     1 SIGS0,TRANC,PHGAR,SUNKNO,FUNKNO
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKI
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IRES(NBM),ISONR(NBISO),NPSYS(NGRO))
      ALLOCATE(VST(NBNRS),SIGTXS(0:NBM),SIGS0X(0:NBM),FLNEW(NBNRS),
     1 SUNKNO(NUN,NGRO),FUNKNO(NUN,NGRO),FLUX(NBM,NGRO),SIG0(NBM,NGRO),
     2 SIG1(NBM,NGRO),SIG3(NBM,NGRO),TOTAL(NGRO,NBNRS),
     3 SIGS0(NGRO,NBNRS),TRANC(NGRO,NBNRS),PHGAR(NGRO,NBNRS))
      ALLOCATE(MASKI(NBISO))
      ALLOCATE(IPISO(NBISO))
*----
*  FIND THE RESONANT MIXTURE NUMBERS AND THE CORRELATED ISOTOPES
*  ASSOCIATED WITH REGION INRS
*----
      CALL XDISET(IRES,NBM,0)
      CALL XDISET(ISONR,NBISO,0)
      IRS=0
      TEXT12=' '
      DO 30 IBM=1,NBM
      LOGDO=.FALSE.
      DO 10 I=1,NREG
      LOGDO=LOGDO.OR.(MAT(I).EQ.IBM)
   10 CONTINUE
      IF(.NOT.LOGDO) GO TO 30
      DO 20 ISO=1,NBISO
      IF((MIX(ISO).EQ.IBM).AND.(LSHI(ISO).EQ.INRS)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,ISO),I0=1,3)
         IF(HNAMIS.NE.TEXT12) THEN
           IRS=IRS+1
           TEXT12=HNAMIS
         ENDIF
         ISONR(ISO)=IRS
         IRES(IBM)=IRS
      ENDIF
   20 CONTINUE
   30 CONTINUE
      IF(IRS.NE.NBNRS) CALL XABORT('TONSPH: INVALID VALUE OF NBNRS.')
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
      CALL LIBIPS(IPLIB,NBISO,IPISO)
*----
*  UNLOAD MICROSCOPIC X-S FROM LCM TO SCRATCH STORAGE.
*----
      DO 40 ISO=1,NBISO
      IRS=ISONR(ISO)
      IF(IRS.GT.0) THEN
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        CALL LCMGET(KPLIB,'NTOT0',TOTAL(1,IRS))
        CALL LCMGET(KPLIB,'SIGS00',SIGS0(1,IRS))
        DO IGRP=IGRMIN,IGRMAX
*         Compute a ST flux for the homogeneous equivalent medium.
          PHGAR(IGRP,IRS)=MAX(0.0,SN(IGRP,ISO)/(SN(IGRP,ISO)+
     1    (TOTAL(IGRP,IRS)-SIGS0(IGRP,IRS))))
        ENDDO
        IF(ITRANC.NE.0) CALL LCMGET(KPLIB,'TRANC',TRANC(1,IRS))
      ENDIF
   40 CONTINUE
*----
*  COMPUTE THE MERGED VOLUMES.
*----
      NALBP=0
      LHOMOG=.TRUE.
      CALL XDRSET(VST,NBNRS,0.0)
      DO 50 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 50
      IND=IRES(IBM)
      IF(IND.EQ.0) THEN
         LHOMOG=.FALSE.
      ELSE
         VST(IND)=VST(IND)+VOL(I)
      ENDIF
   50 CONTINUE
      IF(LHOMOG.AND.(NBNRS.EQ.1)) GO TO 260
      IF(IMPX.GE.3) WRITE(6,'(37H TONSPH: SPH FACTOR CALCULATION (NBNR,
     1 2HS=,I5,1H)/)') NBNRS
*----
*  SET THE MIXTURE-DEPENDENT MACROSCOPIC XS.
*----
      CALL XDRSET(FUNKNO,NUN*NGRO,0.0)
      CALL XDRSET(SUNKNO,NUN*NGRO,0.0)
      CALL XDISET(NPSYS,NGRO,0)
      CALL LCMSIX(IPLIB,'SHIBA',1)
      JPLIB=LCMLID(IPLIB,'GROUP',NGRO)
      DO 110 IGRP=IGRMIN,IGRMAX
        NPSYS(IGRP)=IGRP
*
*       COMPUTE THE LIGHT AND RESONANT COMPONENTS OF THE MACROSCOPIC
*       CROSS SECTIONS IN EACH RESONANT MIXTURE.
        DO 70 IBM=1,NBM
        SIG0(IBM,IGRP)=0.0
        SIG1(IBM,IGRP)=0.0
        SIG3(IBM,IGRP)=SIGT3(IBM,IGRP)
   70   CONTINUE
        DO 80 ISO=1,NBISO
        IRS=ISONR(ISO)
        IF(IRS.GT.0) THEN
          IBM=MIX(ISO)
          FLUX(IBM,IGRP)=PHGAR(IGRP,IRS)
          SIGT2(IBM,IGRP)=SIGT2(IBM,IGRP)-TOTAL(IGRP,IRS)*DEN(ISO)
          SIG0(IBM,IGRP)=TOTAL(IGRP,IRS)*DEN(ISO)
          SIG1(IBM,IGRP)=SIGS0(IGRP,IRS)*DEN(ISO)
          IF(ITRANC.NE.0) THEN
            SIG3(IBM,IGRP)=SIGT3(IBM,IGRP)-TRANC(IGRP,IRS)*DEN(ISO)
          ENDIF
        ENDIF
   80   CONTINUE
        IF(IMPX.GE.10) THEN
          WRITE (6,400) IGRP,(SIG0(I,IGRP),I=1,NBM)
          WRITE (6,410) IGRP,(SIG1(I,IGRP),I=1,NBM)
          WRITE (6,420) IGRP,(SIGT2(I,IGRP),I=1,NBM)
          WRITE (6,430) IGRP,(FLUX(I,IGRP),I=1,NBM)
        ENDIF
*----
*  COMPUTE THE SOURCES.
*----
        DO 90 I=1,NREG
        IBM=MAT(I)
        IF(IBM.EQ.0) GO TO 90
        SUNKNO(KEYFLX(I),IGRP)=SIGT2(IBM,IGRP)
        IF(IRES(IBM).GT.0) THEN
          SUNKNO(KEYFLX(I),IGRP)=SUNKNO(KEYFLX(I),IGRP)+FLUX(IBM,IGRP)*
     1    (SIG1(IBM,IGRP)-SIG0(IBM,IGRP))
          IF(.NOT.LHOMOG) SUNKNO(KEYFLX(I),IGRP)=SUNKNO(KEYFLX(I),IGRP)-
     1    FLUX(IBM,IGRP)*SIGT2(IBM,IGRP)
        ENDIF
   90   CONTINUE
*
        IF(NPSYS(IGRP).NE.0) THEN
          ICPIJ=ICPIJ+1
          SIGTXS(0)=0.0
          SIGS0X(0)=0.0
          DO 100 IBM=1,NBM
          IND=IRES(IBM)
          IF((ITRANC.NE.0).AND.(IND.EQ.0)) THEN
            SIGTXS(IBM)=SIGT2(IBM,IGRP)-SIG3(IBM,IGRP)
          ELSE
            SIGTXS(IBM)=SIGT2(IBM,IGRP)
          ENDIF
          IF(IND.EQ.0) THEN
*           REMOVE TRANSPORT CORRECTION.
            IF(ITRANC.NE.0) THEN
              SIGS0X(IBM)=-SIG3(IBM,IGRP)
            ELSE
              SIGS0X(IBM)=0.0
            ENDIF
          ELSE
*           BELL ACCELERATION.
            SIGTXS(IBM)=SIGTXS(IBM)+SIG0(IBM,IGRP)
            SIGS0X(IBM)=SIGTXS(IBM)
            IF(LHOMOG) SIGS0X(IBM)=SIGS0X(IBM)-SIGT2(IBM,IGRP)
          ENDIF
  100     CONTINUE
          KPLIB=LCMDIL(JPLIB,IGRP)
          CALL LCMPUT(KPLIB,'DRAGON-TXSC',NBM+1,2,SIGTXS)
          CALL LCMPUT(KPLIB,'DRAGON-S0XSC',NBM+1,2,SIGS0X)
        ENDIF
  110 CONTINUE
*----
*  SOLVE FOR THE FLUX USING DIRECT SELF-SHIELDED CROSS SECTIONS
*----
      ISTRM=1
      NANI=1
      NW=0
      KNORM=1
      IMPY=MAX(0,IMPX-3)
      IF(IPHASE.EQ.1) THEN
*        USE A NATIVE DOOR.
         CALL DOORAV(CDOOR,JPLIB,NPSYS,IPTRK,IFTRAK,IMPY,NGRO,NREG,
     1   NBM,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         IPIJK=1
         CALL DOORPV(CDOOR,JPLIB,NPSYS,IPTRK,IFTRAK,IMPY,NGRO,NREG,
     1   NBM,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,.FALSE.,TITR,NALBP)
      ENDIF
      IDIR=0
      LEXAC=.FALSE.
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,JPLIB,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NBM,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUNKNO,FUNKNO,IPMACR,
     2 IPSOU,REBFLG)
      CALL LCMSIX(IPLIB,' ',2)
*----
*  HOMOGENIZE THE FLUX
*----
      DO 150 IGRP=IGRMIN,IGRMAX
      IF(NPSYS(IGRP).NE.0) THEN
        CALL XDRSET(FLNEW,NBNRS,0.0)
        DO 120 I=1,NREG
        IF(MAT(I).EQ.0) GO TO 120
        IND=IRES(MAT(I))
        IF(IND.GT.0) FLNEW(IND)=FLNEW(IND)+FUNKNO(KEYFLX(I),IGRP)*VOL(I)
  120   CONTINUE
        DO 130 IND=1,NBNRS
        FLNEW(IND)=FLNEW(IND)/VST(IND)
  130   CONTINUE
*----
*  SPH FACTOR CONTROL
*----
        DO 140 IBM=1,NBM
        IND=IRES(IBM)
        IF(IND.GT.0) THEN
          SPHNEW=PHGAR(IGRP,IND)/FLNEW(IND)
          LPROB=(SPHNEW.LE.0.).OR.(SPHNEW.GT.1.).OR.(FLNEW(IND).LT.0.05)
          IF(LPROB) SPHNEW=1.0
          SPH(IBM,IGRP)=SPHNEW
        ENDIF
  140   CONTINUE
      ENDIF
  150 CONTINUE
*----
*  SPH CORRECTION OF THE MICROLIB
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NL=ISTATE(4)
      NED=ISTATE(13)
      NDEL=ISTATE(19)
      DO 160 ISO=1,NBISO
      MASKI(ISO)=(ISONR(ISO).GT.0)
  160 CONTINUE
      CALL TONCMI(IPLIB,IMPX,NBM,NBISO,NGRO,NL,NED,NDEL,MASKI,SPH)
      IF(IMPX.GT.3) THEN
        DO 170 IGRP=IGRMIN,IGRMAX
        WRITE (6,440) IGRP,(SPH(IBM,IGRP),IBM=1,NBM)
  170   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  260 DEALLOCATE(IPISO)
      DEALLOCATE(MASKI)
      DEALLOCATE(PHGAR,TRANC,SIGS0,TOTAL,SIG3,SIG1,SIG0,FLUX,FUNKNO,
     1 SUNKNO,FLNEW,SIGS0X,SIGTXS,VST)
      DEALLOCATE(NPSYS,ISONR,IRES)
      RETURN
  400 FORMAT(/51H TOTAL MACROSCOPIC CROSS SECTIONS OF THE RESONANT M,
     1 31HATERIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
  410 FORMAT(/51H SCATTERING MACROSCOPIC CROSS SECTIONS OF THE OTHER,
     1 33H MATERIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
  420 FORMAT(/51H TOTAL MACROSCOPIC CROSS SECTIONS OF THE OTHER MATE,
     1 28HRIALS IN EACH MIXTURE (GROUP,I5,2H):/(1X,1P,11E11.3))
  430 FORMAT(/19H TABSN3 FLUX (GROUP,I5,2H):/(1X,1P,11E11.3))
  440 FORMAT(/19H SPH FACTORS (GROUP,I5,2H):/(1X,1P,11E11.3))
      END
