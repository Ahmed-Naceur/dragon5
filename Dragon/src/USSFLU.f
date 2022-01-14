*DECK USSFLU
      SUBROUTINE USSFLU(MAXNOR,IPTRK,IPLIB,IPLI0,IFTRAK,NREG,NUN,NBMIX,
     1 NBISO,NIRES,NL,NED,NDEL,ISONAM,ISOBIS,HCAL,MAT,VOL,KEYFLX,CDOOR,
     2 LEAKSW,IMPX,DEN,MIX,IAPT,IPHASE,NGRP,IGRMIN,IGRMAX,NBNRS,IREX,
     3 TITR,ICORR,ISUBG,MAXST,UNGAR,PHGAR,STGAR,SFGAR,SSGAR,S0GAR,
     4 SAGAR,SDGAR,SWGAR,ISMIN,ISMAX,MASKG,SIGGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the self-shielded cross sections in each energy group using a
* subgroup approach.
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
* MAXNOR  maximum order of the probability tables (PT).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPLIB   pointer to the internal microscopic cross section library
*         with subgroups (L_LIBRARY signature).
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IFTRAK  file unit number used to store the tracks.
* NREG    number of regions.
* NUN     number of unknowns per energy group and band.
* NBMIX   number of mixtures in the internal library.
* NBISO   number of isotopes.
* NIRES   number of correlated resonant isotopes.
* NL      number of legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* ISONAM  alias name of isotopes in IPLIB.
* ISOBIS  alias name of isotopes in IPLI0.
* HCAL    name of the self-shielding calculation.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* CDOOR   name of the geometry/solution operator.
* LEAKSW  leakage flag (LEAKSW=.true. if neutron leakage through
*         external boundary is present).
* IMPX    print flag (equal to zero for no print).
* DEN     density of each isotope.
* MIX     mix number of each isotope (can be zero).
* IAPT    resonant isotope index associated with isotope I. Mixed
*         moderator if IAPT(I)=NIRES+1. Out-of-fuel isotope if
*         IAPT(I)=0.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* NGRP    number of energy groups.
* IGRMIN  first group where the self-shielding is applied.
* IGRMAX  most thermal group where the self-shielding is applied.
* NBNRS   number of correlated fuel regions. Note that NBNRS=max(IREX).
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* TITR    title.
* ICORR   mutual resonance shielding flag (=1 to suppress the model
*         in cases it is required in LIB operator).
* ISUBG   type of self-shielding model (=1 use physical probability
*         tables; =3 use original Ribon method; =4 use Ribon extended
*         method).
* MAXST   maximum number of fixed point iterations for the ST scattering
*         source.
*
*Parameters: output
* UNGAR   averaged fluxes per volume.
* PHGAR   averaged fluxes.
* STGAR   microscopic self-shielded total x-s.
* SFGAR   microscopic self-shielded fission x-s.
* SSGAR   microscopic self-shielded scattering x-s.
* S0GAR   microscopic transfer scattering xs (isotope,secondary,
*         primary).
* SAGAR   microscopic self-shielded additional xs.
* SDGAR   microscopic self-shielded delayed nu-sigf xs.
* SWGAR   microscopic secondary slowing-down cross sections (ISUBG=4).
* ISMIN   first secondary group indices.
* ISMAX   last secondary group indices.
* MASKG   energy group mask pointing on self-shielded groups.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering; (*,*,*,4) flux times P0 scattering.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPLIB,IPLI0
      INTEGER MAXNOR,IFTRAK,NREG,NUN,NBMIX,NBISO,NIRES,NL,NED,NDEL,
     1 ISONAM(3,NBISO),ISOBIS(3,NBISO),MAT(NREG),KEYFLX(NREG),IMPX,
     2 MIX(NBISO),IAPT(NBISO),IPHASE,NGRP,IGRMIN,IGRMAX,NBNRS,
     3 IREX(NBMIX),ICORR,ISUBG,MAXST,ISMIN(NL,NGRP),ISMAX(NL,NGRP)
      REAL VOL(NREG),DEN(NBISO),UNGAR(NREG,NIRES,NGRP),
     1 PHGAR(NBNRS,NIRES,NGRP),STGAR(NBNRS,NIRES,NGRP),
     2 SFGAR(NBNRS,NIRES,NGRP),SSGAR(NBNRS,NIRES,NL,NGRP),
     3 S0GAR(NBNRS,NIRES,NL,NGRP,NGRP),SAGAR(NBNRS,NIRES,NED,NGRP),
     4 SDGAR(NBNRS,NIRES,NDEL,NGRP),SWGAR(NBNRS,NIRES,NGRP),
     5 SIGGAR(NBMIX,0:NIRES,NGRP,4)
      LOGICAL LEAKSW,MASKG(NGRP)
      CHARACTER HCAL*12,CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPP,KPLIB,LPLIB,MPLIB,JPLI0,KPLI0
      LOGICAL LRES,LLIB,LRIB
      PARAMETER (MAXED=50)
      CHARACTER TEXT12*12,HVECT(MAXED)*8,CBDPNM*12,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPPT1,IPISO1,IPISO2
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWRK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NOR,IPPT2,ISM
      REAL, ALLOCATABLE, DIMENSION(:) :: GAS,GA1,VOLMER,DELTA,GOLD2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GOLD,GA2,CONR
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: XFLUX
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASKH
*----
*  SCRATCH STORAGE ALLOCATION
*   GOLD    Goldstein-Cohen parameters.
*----
      ALLOCATE(IPPT1(NIRES))
      ALLOCATE(NOR(NIRES,NGRP),IPPT2(NIRES,5),IWRK(NGRP),ISM(2,NL))
      ALLOCATE(GOLD(NIRES,NGRP),GAS(NGRP),GA1(NGRP),GA2(NGRP,NGRP),
     1 XFLUX(NBNRS,MAXNOR,NL,NIRES),CONR(NBNRS,NIRES),VOLMER(0:NBNRS),
     2 DELTA(NGRP))
      ALLOCATE(MASKH(NGRP))
      ALLOCATE(IPISO1(NBISO),IPISO2(NBISO))
*
      CALL KDRCPU(TK1)
      DO 15 IG1=1,NGRP
      DO 10 IL=1,NL
      ISMIN(IL,IG1)=NGRP
      ISMAX(IL,IG1)=1
   10 CONTINUE
   15 CONTINUE
      DO 20 IRES=1,NIRES
      NOR(IRES,1)=-1
   20 CONTINUE
*
      IF(NED.GT.0) THEN
         IF(NED.GT.MAXED) CALL XABORT('USSFLU: INVALID VALUE OF MAXED.')
         CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED,HVECT)
      ENDIF
*
      CALL LIBIPS(IPLIB,NBISO,IPISO1)
      CALL LIBIPS(IPLI0,NBISO,IPISO2)
      CALL XDRSET(SIGGAR,NBMIX*(NIRES+1)*NGRP*4,0.0)
      DO 190 ISO=1,NBISO
      IBM=MIX(ISO)
      DO 30 I=1,NREG
      IF(MAT(I).EQ.IBM) GO TO 35
   30 CONTINUE
      GO TO 190
   35 IRES=IAPT(ISO)
      DENN=DEN(ISO)
      JRES=IRES
      IF(IRES.EQ.NIRES+1) JRES=0
*----
*  RECOVER INFINITE DILUTION OR SELF-SHIELDED CROSS SECTIONS AND
*  COMPUTE OUT-OF-FUEL MACROSCOPIC CROSS SECTIONS.
*----
      KPLI0=IPISO2(ISO) ! set ISO-th isotope
      IF(C_ASSOCIATED(KPLI0)) THEN
         CALL LCMLEN(KPLI0,'NTOT0',ILEN0,ITYLCM)
         IF(ILEN0.NE.0) THEN
            LLIB=.FALSE.
            IPP=KPLI0
         ELSE
            LLIB=.TRUE.
            IPP=IPISO1(ISO) ! set ISO-th isotope
         ENDIF
      ELSE
         LLIB=.TRUE.
         IPP=IPISO1(ISO) ! set ISO-th isotope
      ENDIF
      IF(LLIB.AND.(.NOT.C_ASSOCIATED(IPP))) THEN
         WRITE(HSMG,'(18H USSFLU: ISOTOPE '',3A4,7H'' (ISO=,I8,5H) IS ,
     1   39HNOT AVAILABLE IN THE ORIGINAL MICROLIB.)') (ISONAM(I0,ISO),
     2   I0=1,3),ISO
         CALL XABORT(HSMG)
      ENDIF
      IF((.NOT.LLIB).AND.(IMPX.GT.2)) WRITE(6,'(/18H USSFLU: RECOVER I,
     1 8HSOTOPE '',3A4,23H'' FROM THE NEW LIBRARY.)') (ISOBIS(I0,ISO),
     2 I0=1,3)
      IF((DENN.NE.0.0).AND.(IBM.NE.0)) THEN
         CALL LCMLEN(IPP,'NTOT0',LENGT,ITYLCM)
         IF(LENGT.NE.NGRP) CALL XABORT('USSFLU: INVALID X-SECTIONS.')
         CALL LCMGET(IPP,'NTOT0',GA1)
         DO 40 IGRP=1,NGRP
         SIGGAR(IBM,JRES,IGRP,1)=SIGGAR(IBM,JRES,IGRP,1)+DENN*GA1(IGRP)
   40    CONTINUE
         CALL LCMGET(IPP,'SIGS00',GA1)
         CALL LCMLEN(IPP,'NWT0',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET(IPP,'NWT0',GAS)
         ELSE
            CALL XDRSET(GAS,NGRP,1.0)
         ENDIF
         DO 45 IGRP=1,NGRP
         SIGGAR(IBM,JRES,IGRP,3)=SIGGAR(IBM,JRES,IGRP,3)+DENN*GA1(IGRP)
         SIGGAR(IBM,JRES,IGRP,4)=SIGGAR(IBM,JRES,IGRP,4)+DENN*GA1(IGRP)*
     1   GAS(IGRP)
   45    CONTINUE
         CALL LCMLEN(IPP,'TRANC',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET(IPP,'TRANC',GA1)
         ELSE
            CALL XDRSET(GA1,NGRP,0.0)
         ENDIF
         DO 50 IGRP=1,NGRP
         SIGGAR(IBM,JRES,IGRP,2)=SIGGAR(IBM,JRES,IGRP,2)+DENN*GA1(IGRP)
   50    CONTINUE
      ENDIF
      CALL LCMGET(IPLI0,'DELTAU',DELTA)
*----
*  RECOVER PROBABILITY TABLE INFORMATION.
*----
      IF((IRES.GT.0).AND.(IRES.LE.NIRES)) THEN
         IF(NOR(IRES,1).EQ.-1) THEN
            KPLIB=IPISO1(ISO) ! set ISO-th isotope
*
*           RECOVER INFINITE DILUTION VALUES.
            CALL LCMGET(KPLIB,'NTOT0',GAS)
            DO 57 IG=1,NGRP
            DO 56 IND=1,NBNRS
            STGAR(IND,IRES,IG)=GAS(IG)
            SFGAR(IND,IRES,IG)=0.0
            SWGAR(IND,IRES,IG)=0.0
            DO 51 IED=1,NED
            SAGAR(IND,IRES,IED,IG)=0.0
   51       CONTINUE
            DO 55 IDEL=1,NDEL
            SDGAR(IND,IRES,IDEL,IG)=0.0
   55       CONTINUE
   56       CONTINUE
   57       CONTINUE
            CALL LCMLEN(KPLIB,'NUSIGF',ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
               CALL LCMGET(KPLIB,'NUSIGF',GAS)
               DO 65 IG=1,NGRP
               DO 60 IND=1,NBNRS
               SFGAR(IND,IRES,IG)=GAS(IG)
   60       CONTINUE
   65       CONTINUE
            ENDIF
            DO 80 IL=1,NL
            CALL XDRLGS(KPLIB,-1,IMPX,IL-1,IL-1,1,NGRP,GAS,GA2,ITYPRO)
*           JG IS THE SECONDARY GROUP.
            DO 72 IG=1,NGRP
            DO 71 IND=1,NBNRS
            SSGAR(IND,IRES,IL,IG)=GAS(IG)
            DO 70 JG=1,NGRP
            IF(IL.EQ.1) THEN
              SWGAR(IND,IRES,JG)=SWGAR(IND,IRES,JG)+GA2(JG,IG)*DELTA(IG)
            ENDIF
            S0GAR(IND,IRES,IL,JG,IG)=GA2(JG,IG)
   70       CONTINUE
   71       CONTINUE
   72       CONTINUE
   80       CONTINUE
            DO 95 IG=1,NGRP
            DO 90 IND=1,NBNRS
            SWGAR(IND,IRES,IG)=SWGAR(IND,IRES,IG)/DELTA(IG)
   90       CONTINUE
   95       CONTINUE
            DO 110 IED=1,NED
            CALL LCMLEN(KPLIB,HVECT(IED),ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
               CALL LCMGET(KPLIB,HVECT(IED),GAS)
               DO 105 IG=1,NGRP
               DO 100 IND=1,NBNRS
               SAGAR(IND,IRES,IED,IG)=GAS(IG)
  100          CONTINUE
  105          CONTINUE
            ENDIF
  110       CONTINUE
            DO 130 IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            CALL LCMLEN(KPLIB,TEXT12,ILENGT,ITYLCM)
            IF(ILENGT.GT.0) THEN
               CALL LCMGET(KPLIB,TEXT12,GAS)
               DO 125 IG=1,NGRP
               DO 120 IND=1,NBNRS
               SDGAR(IND,IRES,IDEL,IG)=GAS(IG)
  120          CONTINUE
  125          CONTINUE
            ENDIF
  130       CONTINUE
*
            CALL LCMSIX(KPLIB,'PT-PHYS',1)
               CALL LCMGET(KPLIB,'NOR',IWRK)
               LPLIB=LCMGID(KPLIB,'GROUP')
               DO 150 IG1=1,NGRP
               IF(IWRK(IG1).GT.1) THEN
                  MPLIB=LCMGIL(LPLIB,IG1)
                  CALL LCMGET(MPLIB,'ISM-LIMITS',ISM)
                  DO 140 IL=1,NL
                  ISMIN(IL,IG1)=MIN(ISMIN(IL,IG1),ISM(1,IL))
                  ISMAX(IL,IG1)=MAX(ISMAX(IL,IG1),ISM(2,IL))
  140             CONTINUE
               ENDIF
               NOR(IRES,IG1)=IWRK(IG1)
  150          CONTINUE
            CALL LCMSIX(KPLIB,' ',2)
            CALL LCMLEN(KPLIB,'NGOLD',LENGT,ITYLCM)
            IF(LENGT.GT.0) THEN
               ALLOCATE(GOLD2(NGRP))
               CALL LCMGET(KPLIB,'NGOLD',GOLD2)
               DO 160 IG1=1,NGRP
               GOLD(IRES,IG1)=GOLD2(IG1)
  160          CONTINUE
               DEALLOCATE(GOLD2)
            ELSE
               DO 170 IG1=1,NGRP
               GOLD(IRES,IG1)=1.0
  170          CONTINUE
            ENDIF
            CALL LCMLEN(KPLIB,'BIN-NFS',LENGT,ITYLCM)
            IF(LENGT.GT.0) THEN
               CALL LCMGET(KPLIB,'BIN-NFS',IWRK)
               DO 180 IG1=1,NGRP
               IF((GOLD(IRES,IG1).LT.-900.).AND.(IWRK(IG1).EQ.0)) THEN
                  GOLD(IRES,IG1)=1.0
               ENDIF
  180          CONTINUE
            ENDIF
         ENDIF
      ENDIF
  190 CONTINUE
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H USSFLU: CPU TIME SPENT TO RECOVER,
     1 23H INFINITE-DILUTION XS =,F8.1,8H SECOND./)') TK2-TK1
*
      CALL KDRCPU(TK1)
      TK4=0.0
      TK5=0.0
      ICPIJ=0
*----
*  COMPUTE THE MERGED VOLUMES AND NUMBER DENSITIES.
*----
      DO 200 I=0,NBNRS
      VOLMER(I)=0.0
  200 CONTINUE
      DO 210 I=1,NREG
      IBM=MAT(I)
      IF(IBM.GT.0) VOLMER(IREX(IBM))=VOLMER(IREX(IBM))+VOL(I)
  210 CONTINUE
      DO 225 JRES=1,NIRES
      DO 220 IND=1,NBNRS
      CONR(IND,JRES)=0.0
  220 CONTINUE
  225 CONTINUE
      DO 240 ISO=1,NBISO
      JRES=IAPT(ISO)
      IF((JRES.GT.0).AND.(JRES.LE.NIRES)) THEN
         DENN=DEN(ISO)
         DO 230 IREG=1,NREG
         IBM=MAT(IREG)
         IF(MIX(ISO).EQ.IBM) THEN
            IND0=IREX(IBM)
            IF(IND0.EQ.0) CALL XABORT('USSFLU: IREX FAILURE.')
            CONR(IND0,JRES)=CONR(IND0,JRES)+DENN*VOL(IREG)/VOLMER(IND0)
         ENDIF
  230    CONTINUE
      ENDIF
  240 CONTINUE
*----
*  RECOVER POSITION OF PROBABILITY TABLES AND NAME OF RESONANT ISOTOPE.
*----
      DO 270 IRES=1,NIRES
      ISOT=0
      DO 250 JSOT=1,NBISO
      IF(IAPT(JSOT).EQ.IRES) THEN
         ISOT=JSOT
         GO TO 260
      ENDIF
  250 CONTINUE
      CALL XABORT('USSFLU: UNABLE TO FIND A RESONANT ISOTOPE.')
  260 KPLIB=IPISO1(ISOT) ! set ISOT-th isotope
      CALL LCMLEN(KPLIB,'PT-PHYS',ILONG,ITYLCM)
      IF(ILONG.EQ.0) CALL XABORT('USSFLU: BUG1.')
      CALL LCMSIX(KPLIB,'PT-PHYS',1)
      CALL LCMGET(KPLIB,'NDEL',NDEL0)
      IF(NDEL0.GT.NDEL) CALL XABORT('USSFLU: NDEL OVERFLOW.')
      CALL LCMLEN(KPLIB,'GROUP',LENG,ITYLCM)
      IF(ILONG.EQ.0) CALL XABORT('USSFLU: BUG2.')
      IPPT1(IRES)=LCMGID(KPLIB,'GROUP')
      CALL LCMSIX(KPLIB,' ',2)
      IPPT2(IRES,1)=IREX(MIX(ISOT))
      IPPT2(IRES,2)=ISONAM(1,ISOT)
      IPPT2(IRES,3)=ISONAM(2,ISOT)
      IPPT2(IRES,4)=ISONAM(3,ISOT)
      IPPT2(IRES,5)=NDEL0
      IF(IPPT2(IRES,1).LE.0) CALL XABORT('USSFLU: BUG3.')
  270 CONTINUE
*----
*  DETERMINE WHICH GROUPS ARE SELF-SHIELDED.
*----
      DO 290 IGRP=1,NGRP
      LRES=.FALSE.
      DO 280 IRES=1,NIRES
      LRES=LRES.OR.(NOR(IRES,IGRP).GT.1)
  280 CONTINUE
      MASKG(IGRP)=((IGRP.GE.IGRMIN).AND.(IGRP.LE.IGRMAX).AND.LRES)
  290 CONTINUE
*----
*  INITIALIZATION OF THE MULTIBAND FLUXES AND SOURCES.
*----
      CALL LCMSIX(IPLI0,'SHIBA_SG',1)
      CALL LCMSIX(IPLI0,HCAL,1)
      DO 310 IRES=1,NIRES
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      JPLI0=LCMLID(IPLI0,'GROUP_INFO',NGRP)
      CALL XDRSET(XFLUX(1,1,1,IRES),NBNRS*MAXNOR*NL,1.0)
      DO 300 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(NOR(IRES,IGRP).GT.1)) THEN
         CALL LCMLEL(JPLI0,IGRP,LENG0,ITYLCM)
         IF(LENG0.EQ.0) THEN
            KPLI0=LCMDIL(JPLI0,IGRP)
            CALL LCMPUT(KPLI0,'NWT0-PT',NBNRS*MAXNOR*NL,2,
     1      XFLUX(1,1,1,IRES))
         ENDIF
      ENDIF
  300 CONTINUE
      CALL LCMSIX(IPLI0,' ',2)
  310 CONTINUE
*
      CALL KDRCPU(TKA)
*
      DO 340 IRES=1,NIRES
      LRIB=.FALSE.
      DO 330 IGRP=1,NGRP
      MASKH(IGRP)=MASKG(IGRP).AND.(NOR(IRES,IGRP).GT.1)
      LRIB=LRIB.OR.(MASKH(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.))
      IF(MASKH(IGRP)) ICPIJ=ICPIJ+NOR(IRES,IGRP)
  330 CONTINUE
*----
*  ITERATIVE APPROACH FOR THE HELIOS/WIMS-7 METHOD.
*----
      CALL USSIT1(MAXNOR,NGRP,MASKH,IRES,IPLI0,IPTRK,IFTRAK,CDOOR,
     1 IMPX,NBMIX,NREG,NUN,NL,IPHASE,MAXST,MAT,VOL,KEYFLX,LEAKSW,IREX,
     2 SIGGAR,TITR,NIRES,NBNRS,NOR,CONR,GOLD,IPPT1,IPPT2,STGAR,SSGAR,
     3 VOLMER,UNGAR)
*----
*  ITERATIVE APPROACH FOR THE SUBGROUP PROJECTION METHOD.
*----
      CALL USSIST(MAXNOR,NGRP,MASKH,IRES,IPLI0,IPTRK,IFTRAK,CDOOR,
     1 IMPX,NBMIX,NREG,NUN,NL,IPHASE,MAXST,MAT,VOL,KEYFLX,LEAKSW,
     2 IREX,SIGGAR,TITR,ICORR,NIRES,NBNRS,NOR,CONR,GOLD,IPPT1,IPPT2,
     3 STGAR,SSGAR,VOLMER,UNGAR)
*----
*  RESPONSE MATRIX APPROACH FOR THE RIBON EXTENDED METHOD.
*----
      IF(LRIB) THEN
         CALL USSIT0(MAXNOR,NGRP,MASKH,IRES,IPLI0,IPTRK,IFTRAK,CDOOR,
     1   IMPX,NBMIX,NREG,NUN,NL,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IREX,
     2   SIGGAR,TITR,ICORR,NIRES,NBNRS,NOR,CONR,GOLD,IPPT1,IPPT2,STGAR,
     3   SSGAR,SWGAR,VOLMER,UNGAR)
      ENDIF
  340 CONTINUE
      CALL KDRCPU(TKB)
      TK4=TK4+(TKB-TKA)
*----
*  COMPUTE THE SELF-SHIELDED REACTION RATES.
*----
      CALL XDRSET(PHGAR,NBNRS*NIRES*NGRP,1.0)
      DO 360 IGRP=1,NGRP
      IF(MASKG(IGRP)) THEN
         MAXXS=2+NL+NED+NDEL
         DO 350 IL=1,NL
         MAXXS=MAXXS+MAX(ISMAX(IL,IGRP)-ISMIN(IL,IGRP)+1,0)
  350    CONTINUE
         IF(ISUBG.EQ.4) MAXXS=MAXXS+1
         CALL KDRCPU(TKA)
         CALL USSIT2(MAXNOR,IPLI0,IGRP,NGRP,ISMIN(1,IGRP),ISMAX(1,IGRP),
     1   NIRES,NBNRS,NL,NED,NDEL,NOR(1,IGRP),IPPT1,IPPT2,GOLD(1,IGRP),
     2   MAXXS,ISUBG,PHGAR(1,1,IGRP),STGAR(1,1,IGRP),SFGAR(1,1,IGRP),
     3   SSGAR(1,1,1,IGRP),S0GAR(1,1,1,1,IGRP),SAGAR(1,1,1,IGRP),
     4   SDGAR(1,1,1,IGRP),SWGAR(1,1,IGRP))
         CALL KDRCPU(TKB)
         TK5=TK5+(TKB-TKA)
      ENDIF
  360 CONTINUE
*     ***************************************************************
      CALL LCMSIX(IPLI0,' ',2)
      CALL LCMSIX(IPLI0,' ',2)
      CALL LCMVAL(IPLI0,' ')
*----
*  RESET MASKG FOR SPH CALCULATION IN SMALL LETHARGY WIDTH GROUPS.
*----
      DO 380 IGRP=1,NGRP
      IF(MASKG(IGRP)) THEN
         LRES=.TRUE.
         DO 370 IRES=1,NIRES
         LRES=LRES.AND.((GOLD(IRES,IGRP).EQ.-998.).OR.
     1   (GOLD(IRES,IGRP).EQ.-1000.))
  370    CONTINUE
         MASKG(IGRP)=.NOT.LRES
         IF(DELTA(IGRP).GT.0.1) MASKG(IGRP)=.TRUE.
      ENDIF
  380 CONTINUE
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/34H USSFLU: CPU TIME SPENT TO COMPUTE,
     1 31H SELF-SHIELDED REACTION RATES =,F8.1,19H SECOND, INCLUDING:
     2 /9X,F8.1,46H SECOND TO BUILD/SOLVE SUBGROUP MATRIX SYSTEM;/9X,
     4 F8.1,38H SECOND TO COMPUTE THE REACTION RATES./9X,9HNUMBER OF,
     5 23H ASSEMBLY DOORS CALLS =,I5,1H.)') TK2-TK1,TK4,TK5,ICPIJ
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO2,IPISO1)
      DEALLOCATE(MASKH)
      DEALLOCATE(DELTA,VOLMER,CONR,XFLUX,GA2,GA1,GAS,GOLD)
      DEALLOCATE(ISM,IWRK,IPPT2,NOR)
      DEALLOCATE(IPPT1)
      RETURN
      END
