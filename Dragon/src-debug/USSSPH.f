*DECK USSSPH
      SUBROUTINE USSSPH(IPLI0,IPTRK,IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,
     1 NL,NED,NDEL,ISONAM,HCAL,MAT,VOL,KEYFLX,CDOOR,LEAKSW,IMPX,DEN,MIX,
     2 IAPT,ITRANC,IPHASE,NGRP,MASKG,NBNRS,IREX,TITR,ISUBG,SIGGAR,ISMIN,
     3 ISMAX,UNGAR,PHGAR,STGAR,SFGAR,SSGAR,S0GAR,SAGAR,SDGAR,SWGAR,SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* SPH equivalence procedure over the self-shielded cross sections. Use
* all the standard solution doors of Dragon.
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
* IPLI0   pointer to the LCM object containing subgroup-related
*         information.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* NREG    number of regions.
* NUN     number of unknowns per energy group.
* NBMIX   number of mixtures in the internal library.
* NBISO   number of isotopes.
* NIRES   number of resonant isotopes.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* ISONAM  alias name of isotopes in IPLIB.
* HCAL    name of the self-shielding calculation.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* CDOOR   name of the geometry/solution operator.
* LEAKSW  leakage flag (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* IMPX    print flag (equal to zero for no print).
* DEN     density of each isotope.
* MIX     mix number of each isotope (can be zero).
* IAPT    resonant isotope index associated with isotope I. Mixed
*         moderator if IAPT(I)=NIRES+1. Out-of-fuel isotope if
*         IAPT(I)=0.
* ITRANC  type of transport correction.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* NGRP    number of energy groups.
* MASKG   energy group mask pointing on self-shielded groups.
* NBNRS   number of totally correlated fuel regions (NBNRS=max(IREX)).
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* TITR    title.
* ISUBG   type of self-shielding model (=1 use physical probability
*         tables; =3 use original Ribon method; =4 use Ribon extended
*         method).
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering; (*,*,*,4) flux times P0 scattering.
* ISMIN   first secondary group indices.
* ISMAX   last secondary group indices.
* UNGAR   averaged fluxes per volume.
* STGAR   microscopic self-shielded total x-s.
*
*Parameters: input/output
* PHGAR   uncorrected and SPH-corrected averaged fluxes.
* SFGAR   uncorrected and SPH-corrected microscopic self-shielded 
*         fission x-s.
* SSGAR   uncorrected and SPH-corrected microscopic 
*         self-shielded scattering x-s.
* S0GAR   uncorrected and SPH-corrected microscopic 
*         transfer scattering x-s
*         (isotope,secondary,primary).
* SAGAR   uncorrected and SPH-corrected microscopic 
*         additional x-s.
* SDGAR   uncorrected and SPH-corrected microscopic 
*         self-shielded delayed nu-sigf x-s.
* SWGAR   uncorrected and SPH-corrected microscopic 
*         secondary slowing-down cross
*         sections (ISUBG=4).
*
*Parameters: output
* SPH     SPH factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK
      INTEGER IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,NL,NED,NDEL,
     1 ISONAM(3,NBISO),MAT(NREG),KEYFLX(NREG),IMPX,MIX(NBISO),
     2 IAPT(NBISO),ITRANC,IPHASE,NGRP,NBNRS,IREX(NBMIX),ISUBG,
     3 ISMIN(NL,NGRP),ISMAX(NL,NGRP)
      REAL VOL(NREG),DEN(NBISO),SIGGAR(NBMIX,0:NIRES,NGRP,4),
     1 UNGAR(NREG,NIRES,NGRP),PHGAR(NBNRS,NIRES,NGRP),
     2 STGAR(NBNRS,NIRES,NGRP),SFGAR(NBNRS,NIRES,NGRP),
     3 SSGAR(NBNRS,NIRES,NL,NGRP),S0GAR(NBNRS,NIRES,NL,NGRP,NGRP),
     4 SAGAR(NBNRS,NIRES,NED,NGRP),SDGAR(NBNRS,NIRES,NDEL,NGRP),
     5 SWGAR(NBNRS,NIRES,NGRP),SPH(NBNRS,NIRES,NGRP)
      LOGICAL LEAKSW,MASKG(NGRP)
      CHARACTER CDOOR*12,HCAL*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLI0,KPLI0
      LOGICAL LHOMOG,LPROB,LTIT,LEXAC
      CHARACTER TEX8*8,CBDPNM*12
      INTEGER NALBP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLMER,SIGTXS,SIGS0X,FLNEW
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKNO,FUNKNO
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIGTI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NPSYS(NGRP))
      ALLOCATE(SIGTI(NBMIX,5,NGRP),VOLMER(NBNRS),SIGTXS(0:NBMIX),
     1 SIGS0X(0:NBMIX),FLNEW(NBNRS),SUNKNO(NUN,NGRP),FUNKNO(NUN,NGRP))
*----
*  COMPUTE THE MERGED VOLUMES.
*----
      NALBP=0
      LHOMOG=.TRUE.
      CALL XDRSET(VOLMER,NBNRS,0.0)
      DO 10 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 10
      IND=IREX(IBM)
      IF(IND.EQ.0) THEN
         LHOMOG=.FALSE.
      ELSE
         VOLMER(IND)=VOLMER(IND)+VOL(I)
      ENDIF
   10 CONTINUE
      CALL XDRSET(SPH,NGRP*NIRES*NBNRS,1.0)
      IF(LHOMOG.AND.(NBNRS.EQ.1).AND.(NIRES.EQ.1)) GO TO 260
*----
*  EVALUATION OF THE SPH FACTOR IN THE RESONANT REGION.
*----
      ICPIJ=0
      CALL KDRCPU(TK1)
      CALL LCMSIX(IPLI0,'SHIBA_SG',1)
      CALL LCMSIX(IPLI0,HCAL,1)
      LTIT=.TRUE.
      JPLI0=LCMLID(IPLI0,'GROUP',NGRP)
*----
*  LOOP OVER THE RESONANT ISOTOPES.
*----
      DO 255 IRES=1,NIRES
      CALL XDRSET(FUNKNO,NUN*NGRP,0.0)
      CALL XDRSET(SUNKNO,NUN*NGRP,0.0)
      CALL XDRSET(SIGTI,NBMIX*5*NGRP,0.0)
      CALL XDISET(NPSYS,NGRP,0)
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      KPLI0=LCMGID(IPLI0,'GROUP_INFO')
      DO 50 IGRP=1,NGRP
      IF(.NOT.MASKG(IGRP)) GO TO 50
      CALL LCMLEL(KPLI0,IGRP,ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 50
      NPSYS(IGRP)=IGRP
*----
*  SET THE MIXTURE-DEPENDENT MACROSCOPIC XS.
*----
      DO 35 IBM=1,NBMIX
      DENN=0.0
      DO 20 ISO=1,NBISO
      IF((IAPT(ISO).EQ.IRES).AND.(MIX(ISO).EQ.IBM)) DENN=DEN(ISO)
   20 CONTINUE
      IND=IREX(IBM)
      DO 30 JRES=0,NIRES
      IF(JRES.EQ.0) THEN
         SIGTI(IBM,1,IGRP)=SIGTI(IBM,1,IGRP)+SIGGAR(IBM,0,IGRP,1)
         SIGTI(IBM,3,IGRP)=SIGTI(IBM,3,IGRP)+SIGGAR(IBM,0,IGRP,3)
         IF(ITRANC.NE.0) SIGTI(IBM,2,IGRP)=SIGTI(IBM,2,IGRP)+
     1   SIGGAR(IBM,0,IGRP,2)
      ELSE IF((JRES.EQ.IRES).AND.(IND.NE.0)) THEN
         SIGTI(IBM,5,IGRP)=SIGTI(IBM,5,IGRP)+STGAR(IND,IRES,IGRP)*DENN
         IF(ISUBG.EQ.4) THEN
            SIGTI(IBM,4,IGRP)=SIGTI(IBM,4,IGRP)+SWGAR(IND,IRES,IGRP)*
     1      DENN*PHGAR(IND,IRES,IGRP)
         ELSE
            SIGTI(IBM,4,IGRP)=SIGTI(IBM,4,IGRP)+SSGAR(IND,IRES,1,IGRP)*
     1      DENN*PHGAR(IND,IRES,IGRP)
         ENDIF
      ELSE IF((JRES.GT.0).AND.(IND.NE.0)) THEN
         SIGTI(IBM,1,IGRP)=SIGTI(IBM,1,IGRP)+SIGGAR(IBM,JRES,IGRP,1)
         SIGTI(IBM,3,IGRP)=SIGTI(IBM,3,IGRP)+SIGGAR(IBM,JRES,IGRP,4)
      ENDIF
   30 CONTINUE
   35 CONTINUE
*----
*  COMPUTE THE SOURCES.
*----
      DO 40 I=1,NREG
      IBM=MAT(I)
      IF(IBM.EQ.0) GO TO 40
      SUNKNO(KEYFLX(I),IGRP)=SIGTI(IBM,3,IGRP)
      IF(IREX(IBM).GT.0) THEN
         SUNKNO(KEYFLX(I),IGRP)=SUNKNO(KEYFLX(I),IGRP)+SIGTI(IBM,4,IGRP)
     1   -UNGAR(I,IRES,IGRP)*SIGTI(IBM,5,IGRP)
         IF(.NOT.LHOMOG) SUNKNO(KEYFLX(I),IGRP)=SUNKNO(KEYFLX(I),IGRP)-
     1   UNGAR(I,IRES,IGRP)*SIGTI(IBM,1,IGRP)
      ENDIF
   40 CONTINUE
   50 CONTINUE
      CALL LCMSIX(IPLI0,' ',2)
*
      DO 90 IGRP=1,NGRP
      IF(NPSYS(IGRP).NE.0) THEN
        ICPIJ=ICPIJ+1
        SIGTXS(0)=0.0
        SIGS0X(0)=0.0
        DO 80 IBM=1,NBMIX
        SIGTXS(IBM)=SIGTI(IBM,1,IGRP)-SIGTI(IBM,2,IGRP)
        IND=IREX(IBM)
        IF(IND.EQ.0) THEN
*          REMOVE TRANSPORT CORRECTION.
           SIGS0X(IBM)=-SIGTI(IBM,2,IGRP)
        ELSE
*          BELL ACCELERATION.
           SIGTXS(IBM)=SIGTXS(IBM)+SIGTI(IBM,5,IGRP)
           SIGS0X(IBM)=SIGTXS(IBM)
           IF(LHOMOG) SIGS0X(IBM)=SIGS0X(IBM)-SIGTI(IBM,1,IGRP)
        ENDIF
   80   CONTINUE
        KPLI0=LCMDIL(JPLI0,IGRP)
        CALL LCMPUT(KPLI0,'DRAGON-TXSC',NBMIX+1,2,SIGTXS)
        CALL LCMPUT(KPLI0,'DRAGON-S0XSC',NBMIX+1,2,SIGS0X)
      ENDIF
   90 CONTINUE
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
         CALL DOORAV(CDOOR,JPLI0,NPSYS,IPTRK,IFTRAK,IMPY,NGRP,NREG,
     1   NBMIX,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         IPIJK=1
         CALL DOORPV(CDOOR,JPLI0,NPSYS,IPTRK,IFTRAK,IMPY,NGRP,NREG,
     1   NBMIX,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,.FALSE.,TITR,NALBP)
      ENDIF
      IDIR=0
      LEXAC=.FALSE.
      CALL DOORFV (CDOOR,JPLI0,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NBMIX,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUNKNO,FUNKNO)
*----
*  HOMOGENIZE THE FLUX
*----
      DO 250 IGRP=1,NGRP
      IF(NPSYS(IGRP).NE.0) THEN
        CALL XDRSET(FLNEW,NBNRS,0.0)
        DO 100 I=1,NREG
        IF(MAT(I).EQ.0) GO TO 100
        IND=IREX(MAT(I))
        IF(IND.GT.0) FLNEW(IND)=FLNEW(IND)+FUNKNO(KEYFLX(I),IGRP)*VOL(I)
  100   CONTINUE
        DO 110 IND=1,NBNRS
        FLNEW(IND)=FLNEW(IND)/VOLMER(IND)
  110   CONTINUE
*----
*  SPH FACTOR CONTROL.
*----
        DO 120 IND=1,NBNRS
        SPHNEW=PHGAR(IND,IRES,IGRP)/FLNEW(IND)
        LPROB=(SPHNEW.LE.0.0).OR.(SPHNEW.GT.1.0).OR.(FLNEW(IND).LT.0.05)
        IF(LPROB) SPHNEW=1.0
        SPH(IND,IRES,IGRP)=SPHNEW
  120   CONTINUE
      ENDIF
      IF(MASKG(IGRP)) THEN
        SPHNEW=1.0
        DO 200 IND=1,NBNRS
        SPHNEW=SPH(IND,IRES,IGRP)
        PHGAR(IND,IRES,IGRP)=PHGAR(IND,IRES,IGRP)/SPHNEW
        SFGAR(IND,IRES,IGRP)=SFGAR(IND,IRES,IGRP)*SPHNEW
        IF(ISUBG.EQ.4) SWGAR(IND,IRES,IGRP)=SWGAR(IND,IRES,IGRP)*SPHNEW
        DO 175 IL=1,NL
        IF(MOD(IL-1,2).EQ.0) THEN
          SSGAR(IND,IRES,IL,IGRP)=SSGAR(IND,IRES,IL,IGRP)*SPHNEW+
     1    STGAR(IND,IRES,IGRP)*(1.0-SPHNEW)
        ELSE
          SSGAR(IND,IRES,IL,IGRP)=0.0
        ENDIF
        DO 170 JGRP=ISMIN(IL,IGRP),ISMAX(IL,IGRP)
        IF(MOD(IL-1,2).EQ.0) THEN
          IF(IGRP.EQ.JGRP) THEN
            S0GAR(IND,IRES,IL,IGRP,IGRP)=S0GAR(IND,IRES,IL,IGRP,IGRP)*
     1      SPHNEW+STGAR(IND,IRES,IGRP)*(1.0-SPHNEW)
          ELSE
            S0GAR(IND,IRES,IL,JGRP,IGRP)=S0GAR(IND,IRES,IL,JGRP,IGRP)*
     1      SPHNEW
          ENDIF
        ELSE
          IF(IGRP.EQ.JGRP) THEN
            S0GAR(IND,IRES,IL,IGRP,IGRP)=S0GAR(IND,IRES,IL,IGRP,IGRP)/
     1      SPHNEW+STGAR(IND,IRES,IGRP)*(1.0-1.0/SPHNEW)
          ELSE
            S0GAR(IND,IRES,IL,JGRP,IGRP)=S0GAR(IND,IRES,IL,JGRP,IGRP)/
     1      SPH(IND,IRES,JGRP)
          ENDIF
        ENDIF
        IF(MOD(IL-1,2).EQ.1) THEN
          SSGAR(IND,IRES,IL,IGRP)=SSGAR(IND,IRES,IL,IGRP)+
     1    S0GAR(IND,IRES,IL,JGRP,IGRP)
        ENDIF
  170   CONTINUE
  175   CONTINUE
        DO 180 IED=1,NED
        SAGAR(IND,IRES,IED,IGRP)=SAGAR(IND,IRES,IED,IGRP)*SPHNEW
  180   CONTINUE
        DO 190 IDEL=1,NDEL
        SDGAR(IND,IRES,IDEL,IGRP)=SDGAR(IND,IRES,IDEL,IGRP)*SPHNEW
  190   CONTINUE
  200   CONTINUE
*
        IF(IMPX.GT.1) THEN
           IF(LTIT) THEN
              WRITE(6,'(/42H USSSPH: SPH CORRECTED SELF-SHIELDED MICRO,
     1        23HSCOPIC CROSS SECTIONS (,A12,2H)./6H GROUP,5H FUEL,9X,
     2        4HFLUX,2X,23HSPH FACTOR   ISOTOPE...,8X,5HTOTAL,3X,
     3        10HSCATTERING,3X,10HNU*FISSION,13H WITHIN-GROUP)') HCAL
              LTIT=.FALSE.
           ENDIF
           DO 240 IND=1,NBNRS
           DO 220 ISO=1,NBISO
           IF(IAPT(ISO).EQ.IRES) THEN
              WRITE(TEX8,'(2A4)') (ISONAM(J,ISO),J=1,2)
           ENDIF
  220      CONTINUE
           WRITE(6,'(1X,2I5,1P,E13.4,E12.4,3X,1H'',A8,1H'',4E13.4)')
     1     IGRP,IND,PHGAR(IND,IRES,IGRP),SPH(IND,IRES,IGRP),
     2     TEX8,STGAR(IND,IRES,IGRP),SSGAR(IND,IRES,1,IGRP),
     3     SFGAR(IND,IRES,IGRP),S0GAR(IND,IRES,1,IGRP,IGRP)
  240      CONTINUE
        ENDIF
      ENDIF
  250 CONTINUE
  255 CONTINUE
*     ***************************************************************
      CALL LCMSIX(IPLI0,' ',2)
      CALL LCMSIX(IPLI0,' ',2)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H USSSPH: CPU TIME SPENT TO COMPUTE,
     1 18H THE SPH FACTORS =,F8.1,8H SECOND./9X,17HNUMBER OF ASSEMBL,
     2 15HY DOORS CALLS =,I5,1H.)') TK2-TK1,ICPIJ
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  260 DEALLOCATE(FUNKNO,SUNKNO,FLNEW,SIGS0X,SIGTXS,VOLMER,SIGTI)
      DEALLOCATE(NPSYS)
      RETURN
      END
