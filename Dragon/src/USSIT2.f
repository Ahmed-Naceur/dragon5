*DECK USSIT2
      SUBROUTINE USSIT2(MAXNOR,IPLI0,IGRP,NGRP,ISMIN,ISMAX,NIRES,NBNRS,
     1 NL,NED,NDEL,NOR,IPPT1,IPPT2,GOLD,MAXXS,ISUBG,PHGAR,STGAR,SFGAR,
     2 SSGAR,S0GAR,SAGAR,SDGAR,SWGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute self-shielded microscopic cross sections.
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
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IGRP    energy group under consideration.
* NGRP    number of energy groups.
* ISMIN   minimum secondary group corresponding to group IGRP.
* ISMAX   maximum secondary group corresponding to group IGRP.
* NIRES   exact number of resonant isotopes.
* NBNRS   number of correlated fuel regions.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* NOR     exact order of the probability table.
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope;
*         IPPT2(:,5) number of delayed neutron groups.
* GOLD    Goldstein-Cohen parameters. Set to -999. to enable the Ribon
*         extended method for a specific isotope.
* MAXXS   number of x-s types.
* ISUBG   type of self-shielding model (=1 use physical probability
*         tables; =4 use Ribon extended method).
*
*Parameters: output
* PHGAR   averaged flux.
* STGAR   averaged microscopic total xs in resonant region.
* SFGAR   averaged nu*microscopic fission xs in resonant region.
* SSGAR   averaged microscopic scattering xs in resonant region.
* S0GAR   averaged microscopic transfer scattering xs in resonant
*         region for primary neutrons in current group.
* SAGAR   averaged microscopic self-shielded additional xs.
* SDGAR   microscopic self-shielded delayed nu-sigf xs.
* SWGAR   averaged microscopic secondary slowing-down cross sections
*         (ISUBG=4).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPPT1(NIRES)
      INTEGER MAXNOR,NGRP,ISMIN(NL),ISMAX(NL),NIRES,NBNRS,NL,NED,NDEL,
     1 NOR(NIRES),IPPT2(NIRES,5),MAXXS,ISUBG
      REAL GOLD(NIRES),PHGAR(NBNRS,NIRES),STGAR(NBNRS,NIRES),
     1 SFGAR(NBNRS,NIRES),SSGAR(NBNRS,NIRES,NL),
     2 S0GAR(NBNRS,NIRES,NL,NGRP),SAGAR(NBNRS,NIRES,NED),
     3 SDGAR(NBNRS,NIRES,NDEL),SWGAR(NBNRS,NIRES)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPLIB,JPLI0,KPLI0
      LOGICAL EMPTY,LCM
      CHARACTER HSMG*131,TEXT12*12,TEXX12*12,CBDPNM*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISM
      REAL, ALLOCATABLE, DIMENSION(:) :: CGAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WEIGH,TOTPT,SIGFPT,SIGWPT
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIGSPT,SIGAPT,SIGDPT,XFLUX
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SIG0PT
      TYPE(C_PTR) SIGP_PTR
      REAL, POINTER, DIMENSION(:) :: SIGP
*----
*  SCRATCH STORAGE ALLOCATION
*   ISM     minimum/maximum secondary group indices.
*----
      ALLOCATE(ISM(2,NL))
      ALLOCATE(CGAR(MAXXS),WEIGH(MAXNOR,NIRES),TOTPT(MAXNOR,NIRES),
     1 SIGFPT(MAXNOR,NIRES),SIGSPT(MAXNOR,NIRES,NL),
     2 SIG0PT(MAXNOR,NIRES,NL,NGRP),SIGAPT(MAXNOR,NIRES,NED),
     3 SIGDPT(MAXNOR,NIRES,NDEL),SIGWPT(MAXNOR,NIRES),
     4 XFLUX(NBNRS,MAXNOR,NL))
*----
*  RECOVER THE PROBABILITY TABLE INFORMATION IN CURRENT GROUP.
*----
      DO 110 IRES=1,NIRES
      CALL LCMLEL(IPPT1(IRES),IGRP,ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         KPLIB=LCMGIL(IPPT1(IRES),IGRP)
*        RECOVER PROBABILITY TABLE VALUES FROM PT-PHYS DIRECTORY.
         CALL LCMINF(KPLIB,TEXT12,TEXX12,EMPTY,ILONG,LCM)
         CALL LCMGET(KPLIB,'ISM-LIMITS',ISM)
         CALL LCMLEN(KPLIB,'PROB-TABLE',LENG,ITYLCM)
         IF(LENG.EQ.0) THEN
            CALL XABORT('USSIT2: NO PROBABILITY TABLES PRESENT.')
         ELSE
            NPART=LENG/MAXNOR
         ENDIF
         IF(LCM) THEN
            CALL LCMGPD(KPLIB,'PROB-TABLE',SIGP_PTR)
            CALL C_F_POINTER(SIGP_PTR,SIGP,(/ MAXNOR*NPART /))
         ELSE
            ALLOCATE(SIGP(MAXNOR*NPART))
            CALL LCMGET(KPLIB,'PROB-TABLE',SIGP)
         ENDIF
         CALL LCMLEN(KPLIB,'SIGQT-SLOW',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPLIB,'SIGQT-SIGS',SIGWPT(1,IRES))
         ENDIF
         NDEL0=IPPT2(IRES,5)
         IF(NDEL0.GT.NDEL) CALL XABORT('USSIT2: NDEL OVERFLOW.')
         DO 70 INOR=1,NOR(IRES)
         WEIGH(INOR,IRES)=SIGP(INOR)
         TOTPT(INOR,IRES)=SIGP(MAXNOR+INOR)
         SIGFPT(INOR,IRES)=SIGP(2*MAXNOR+INOR)
         IPP=3
         DO 10 IL=1,NL
         IPP=IPP+1
         SIGSPT(INOR,IRES,IL)=SIGP((IPP-1)*MAXNOR+INOR)
   10    CONTINUE
         DO 35 IL=1,NL
         DO 20 JG=1,NGRP
         SIG0PT(INOR,IRES,IL,JG)=0.0
   20    CONTINUE
         DO 30 JG=ISM(1,IL),ISM(2,IL)
         IPP=IPP+1
         SIG0PT(INOR,IRES,IL,JG)=SIGP((IPP-1)*MAXNOR+INOR)
   30    CONTINUE
   35    CONTINUE
         DO 40 IED=1,NED
         IPP=IPP+1
         SIGAPT(INOR,IRES,IED)=SIGP((IPP-1)*MAXNOR+INOR)
   40    CONTINUE
         DO 50 IDEL=1,NDEL
         SIGDPT(INOR,IRES,IDEL)=0.0
   50    CONTINUE
         DO 60 IDEL=1,NDEL0
         IPP=IPP+1
         SIGDPT(INOR,IRES,IDEL)=SIGP((IPP-1)*MAXNOR+INOR)
   60    CONTINUE
         IF(IPP.NE.NPART) THEN
            WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
            WRITE(HSMG,'(26HUSSIT2: FAILURE. ISOTOPE='',A12,7H'' (IPP=,
     1      I6,7H NPART=,I6,6H IGRP=,I6,2H).)') TEXT12,IPP,NPART,IGRP
            CALL XABORT(HSMG)
         ENDIF
   70    CONTINUE
         IF(.NOT.LCM) DEALLOCATE(SIGP)
      ELSE
*        USE INFINITE DILUTION VALUES.
         IND=IPPT2(IRES,1)
         DO 75 IL=1,NL
         CALL XDRSET(XFLUX(1,1,IL),NBNRS,1.0)
   75    CONTINUE
         WEIGH(1,IRES)=1.0
         TOTPT(1,IRES)=STGAR(IND,IRES)
         SIGFPT(1,IRES)=SFGAR(IND,IRES)
         SIGWPT(1,IRES)=SWGAR(IND,IRES)
         DO 80 IED=1,NED
         SIGAPT(1,IRES,IED)=SAGAR(IND,IRES,IED)
   80    CONTINUE
         DO 90 IDEL=1,NDEL
         SIGDPT(1,IRES,IDEL)=SDGAR(IND,IRES,IDEL)
   90    CONTINUE
         DO 105 IL=1,NL
         SIGSPT(1,IRES,IL)=SSGAR(IND,IRES,IL)
         DO 100 JG=1,NGRP
         SIG0PT(1,IRES,IL,JG)=S0GAR(IND,IRES,IL,JG)
  100    CONTINUE
  105    CONTINUE
      ENDIF
  110 CONTINUE
*----
*  COMPUTE THE SELF-SHIELDED CROSS SECTIONS IN CURRENT GROUP.
*----
      DO 230 K=1,NIRES
      IF(NOR(K).EQ.1) GO TO 230
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') K,NIRES
      NDEL0=IPPT2(K,5)
      CALL LCMSIX(IPLI0,CBDPNM,1)
      JPLI0=LCMGID(IPLI0,'GROUP_INFO')
      KPLI0=LCMGIL(JPLI0,IGRP)
      CALL LCMGET(KPLI0,'NWT0-PT',XFLUX)
      CALL LCMSIX(IPLI0,' ',2)
      DO 220 I=1,NBNRS
      PHGAR(I,K)=0.0
      DO 120 IOF=1,MAXXS
      CGAR(IOF)=0.0
  120 CONTINUE
      DO 170 KINOR=1,NOR(K)
      WW=XFLUX(I,KINOR,1)*WEIGH(KINOR,K)
      PHGAR(I,K)=PHGAR(I,K)+WW
      CGAR(1)=CGAR(1)+TOTPT(KINOR,K)*WW
      CGAR(2)=CGAR(2)+SIGFPT(KINOR,K)*WW
      IOF=2
      JOF=0
      DO 140 IL=1,NL
      IOF=IOF+1
      IF((ISUBG.EQ.4).AND.(GOLD(K).EQ.-999.)) THEN
         WW=XFLUX(I,KINOR,IL)*WEIGH(KINOR,K)
      ENDIF
      CGAR(IOF)=CGAR(IOF)+SIGSPT(KINOR,K,IL)*WW
      JOF=IOF
      DO 130 JGRP=ISMIN(IL),ISMAX(IL)
      JOF=JOF+1
      CGAR(JOF)=CGAR(JOF)+SIG0PT(KINOR,K,IL,JGRP)*WW
  130 CONTINUE
      IOF=JOF
  140 CONTINUE
      IOF=JOF
      WW=XFLUX(I,KINOR,1)*WEIGH(KINOR,K)
      DO 150 IED=1,NED
      IOF=IOF+1
      CGAR(IOF)=CGAR(IOF)+SIGAPT(KINOR,K,IED)*WW
  150 CONTINUE
      DO 160 IDEL=1,NDEL0
      IOF=IOF+1
      CGAR(IOF)=CGAR(IOF)+SIGDPT(KINOR,K,IDEL)*WW
  160 CONTINUE
      IOF=IOF+NDEL-NDEL0
      IF((ISUBG.EQ.4).AND.(GOLD(K).EQ.-999.)) THEN
         IOF=IOF+1
         CGAR(IOF)=CGAR(IOF)+SIGWPT(KINOR,K)*WW
      ELSE IF(ISUBG.EQ.4) THEN
         IOF=IOF+1
         CGAR(IOF)=CGAR(IOF)+SIGSPT(KINOR,K,1)*WW
      ENDIF
      IF(IOF.NE.MAXXS) CALL XABORT('USSIT2: BAD NB OF X-S TYPES.')
  170 CONTINUE
*
      STGAR(I,K)=CGAR(1)/PHGAR(I,K)
      SFGAR(I,K)=CGAR(2)/PHGAR(I,K)
      IOF=2
      DO 195 IL=1,NL
      IOF=IOF+1
      SSGAR(I,K,IL)=CGAR(IOF)/PHGAR(I,K)
      DO 180 JGRP=1,NGRP
      S0GAR(I,K,IL,JGRP)=0.0
  180 CONTINUE
      DO 190 JGRP=ISMIN(IL),ISMAX(IL)
      IOF=IOF+1
      S0GAR(I,K,IL,JGRP)=CGAR(IOF)/PHGAR(I,K)
  190 CONTINUE
  195 CONTINUE
      DO 200 IED=1,NED
      IOF=IOF+1
      SAGAR(I,K,IED)=CGAR(IOF)/PHGAR(I,K)
  200 CONTINUE
      DO 210 IDEL=1,NDEL0
      IOF=IOF+1
      SDGAR(I,K,IDEL)=CGAR(IOF)/PHGAR(I,K)
  210 CONTINUE
      IOF=IOF+NDEL-NDEL0
      IF(ISUBG.EQ.4) SWGAR(I,K)=CGAR(IOF+1)/PHGAR(I,K)
  220 CONTINUE
  230 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XFLUX,SIGWPT,SIGDPT,SIGAPT,SIG0PT,SIGSPT,SIGFPT,TOTPT,
     1 WEIGH,CGAR)
      DEALLOCATE(ISM)
      RETURN
      END
