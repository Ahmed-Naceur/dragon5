*DECK USSIT0
      SUBROUTINE USSIT0(MAXNOR,NGRP,MASKG,IRES,IPLI0,IPTRK,IFTRAK,
     1 CDOOR,IMPX,NBMIX,NREG,NUN,NL,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IREX,
     2 SIGGAR,TITR,ICORR,NIRES,NBNRS,NOR,CONR,GOLD,IPPT1,IPPT2,STGAR,
     3 SSGAR,SWGAR,VOLMER,UNGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the multiband fluxes as required by the subgroup method using
* a response matrix approach (Ribon extended subgroup method):
* a) assume a single resonant isotope;
* b) use the standard solution doors of Dragon.
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
* NGRP    number of energy group.
* MASKG   energy group mask pointing on self-shielded groups.
* IRES    index of the resonant isotope.
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* CDOOR   name of the geometry/solution operator.
* IMPX    print flag (equal to zero for no print).
* NBMIX   number of mixtures in the internal library.
* NREG    number of regions.
* NUN     number of unknowns in the flux or source vector in one
*         energy group and one band.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* LEAKSW  leakage switch (LEAKSW=.TRUE. if neutron leakage through
*         external boundary is present).
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering; (*,*,*,4) flux times P0 scattering.
* TITR    title.
* ICORR   mutual resonance shielding flag (=1 to suppress the model
*         in cases it is required in LIB operator).
* NIRES   exact number of correlated resonant isotopes.
* NBNRS   number of correlated fuel regions.
* NOR     exact order of the probability table.
* CONR    number density of the resonant isotopes.
* GOLD    type of self-shielding model (=1.0 physical probability
*         tables; =-999.0 Ribon extended method).
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope.
* STGAR   averaged microscopic total xs in resonant region.
* SSGAR   averaged microscopic scattering xs in resonant region.
* SWGAR   microscopic secondary slowing-down cross sections (used
*         if GOLD=-999.).
* VOLMER  volumes of the resonant and non-resonant regions.
*
*Parameters: output
* UNGAR   averaged fluxes per volume.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLI0,IPTRK,IPPT1(NIRES)
      INTEGER MAXNOR,NGRP,IRES,IFTRAK,IMPX,NBMIX,NREG,NUN,NL,
     1 IPHASE,MAT(NREG),KEYFLX(NREG),IREX(NBMIX),ICORR,NIRES,NBNRS,
     2 NOR(NIRES,NGRP),IPPT2(NIRES,4)
      REAL VOL(NREG),SIGGAR(NBMIX,0:NIRES,NGRP,4),
     1 CONR(NBNRS,NIRES),GOLD(NIRES,NGRP),STGAR(NBNRS,NIRES,NGRP),
     2 SSGAR(NBNRS,NIRES,NL,NGRP),SWGAR(NBNRS,NIRES,NGRP),
     3 VOLMER(0:NBNRS),UNGAR(NREG,NIRES,NGRP)
      LOGICAL LEAKSW,MASKG(NGRP)
      CHARACTER CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPSYS,KPSYS,JPLI0,KPLIB,KPLI0,KPLI1,JPLI1
      LOGICAL EMPTY,LCM,LEXAC
      CHARACTER CBDPNM*12,TEXT12*12,TEXX12*12,HSMG*131
      INTEGER NALBP
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGTXS,SIGS0X,AWPHI,FUN,SUN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WEIGH,TOTPT,WSLD,SIGWS,PAV,
     1 SIGX
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: XFLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: MATRIX
      TYPE(C_PTR) SIGP_PTR
      REAL, POINTER, DIMENSION(:) :: SIGP
*----
*  STATEMENT FUNCTIONS
*----
      INM(IND,INOR,NBNRS)=(INOR-1)*NBNRS+IND
*----
*  FIND THE NUMBER OF COMPONENTS REQUIRED AND ALLOCATE THE LIST OF
*  ASSEMBLY MATRICES.
*----
      NASM=0
      NALBP=0
      DO 10 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.)) THEN
         NASM=NASM+NOR(IRES,IGRP)
      ENDIF
   10 CONTINUE
      IF(NASM.EQ.0) RETURN
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XFLUX(NBNRS,MAXNOR,NL,NIRES),SIGTXS(0:NBMIX),
     1 SIGS0X(0:NBMIX),AWPHI(0:NBNRS),WEIGH(MAXNOR,NIRES),
     2 TOTPT(MAXNOR,NIRES),WSLD(MAXNOR**2,NIRES),SIGWS(MAXNOR,NIRES),
     3 PAV(0:NBNRS,0:NBNRS),SIGX(NBNRS,NIRES))
      ALLOCATE(MATRIX(NBNRS*MAXNOR,NBNRS*MAXNOR+1))
*----
*  CREATE A SPECIFIC DIRECTORY FOR IRES-TH RESONANT ISOTOPE.
*----
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      JPLI0=LCMGID(IPLI0,'GROUP_INFO')
      IPSYS=LCMLID(IPLI0,'ASSEMB_RIBON',NASM)
      CALL LCMSIX(IPLI0,' ',2)
*----
*  LOOP OVER THE ENERGY GROUPS.
*----
      ALLOCATE(NPSYS(NASM))
      IASM=0
      DO 120 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.)) THEN
         IF(IMPX.GT.1) THEN
            WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
            WRITE(6,'(36H USSIT0: PROCESS CORRELATED ISOTOPE ,A12,
     1      11H WITH INDEX,I3,9H IN GROUP,I4,22H (RESPONSE MATRIX APPR,
     2      6HOACH).)') TEXT12,IRES,IGRP
         ENDIF
         DO 20 JRES=1,NIRES
         IF(GOLD(JRES,IGRP).NE.GOLD(IRES,IGRP)) THEN
            WRITE(HSMG,'(34HUSSIT0: PTSL NOT SET FOR ISOTOPE '',3A4,
     1      10H'' IN GROUP,I4,1H.)') (IPPT2(JRES,J0),J0=2,4),IGRP
            CALL XABORT(HSMG)
         ELSE IF(NOR(JRES,IGRP).GT.MAXNOR) THEN
            CALL XABORT('USSIT0: MAXNOR OVERFLOW.')
         ENDIF
   20    CONTINUE
*----
*  COLLECT THE BASE POINTS IN TOTAL CROSS SECTION.
*----
         NORI=NOR(IRES,IGRP)
         DO 40 JRES=1,NIRES
         CALL LCMLEL(IPPT1(JRES),IGRP,ILONG,ITYLCM)
         IF(ILONG.NE.0) THEN
            KPLIB=LCMGIL(IPPT1(JRES),IGRP)
            CALL LCMINF(KPLIB,TEXT12,TEXX12,EMPTY,ILONG,LCM)
            CALL LCMLEN(KPLIB,'PROB-TABLE',LENG,ITYLCM)
            NPART=LENG/MAXNOR
            IF(LCM) THEN
               CALL LCMGPD(KPLIB,'PROB-TABLE',SIGP_PTR)
               CALL C_F_POINTER(SIGP_PTR,SIGP,(/ MAXNOR*NPART /))
            ELSE
               ALLOCATE(SIGP(MAXNOR*NPART))
               CALL LCMGET(KPLIB,'PROB-TABLE',SIGP)
            ENDIF
            DO 30 INOR=1,NOR(JRES,IGRP)
            WEIGH(INOR,JRES)=SIGP(INOR)
            TOTPT(INOR,JRES)=SIGP(MAXNOR+INOR)
   30       CONTINUE
            IF(.NOT.LCM) DEALLOCATE(SIGP)
         ELSE
            WEIGH(1,JRES)=1.0
            TOTPT(1,JRES)=STGAR(IPPT2(JRES,1),JRES,IGRP)
         ENDIF
   40    CONTINUE
*----
*  SET THE MIXTURE-DEPENDENT CROSS SECTIONS.
*----
         DO 110 INOR=1,NORI
         CALL XDRSET(SIGTXS(0),NBMIX+1,0.0)
         CALL XDRSET(SIGS0X(0),NBMIX+1,0.0)
         DO 90 IBM=1,NBMIX
         IND=IREX(IBM)
         DO 80 JRES=0,NIRES
         IF(JRES.EQ.0) THEN
            SIGTXS(IBM)=SIGTXS(IBM)+(SIGGAR(IBM,0,IGRP,1)-
     1      SIGGAR(IBM,0,IGRP,2))
            SIGS0X(IBM)=SIGS0X(IBM)-SIGGAR(IBM,0,IGRP,2)
         ELSE IF((JRES.NE.IRES).AND.(IND.GT.0)) THEN
            SIGTXS(IBM)=SIGTXS(IBM)+SIGGAR(IBM,JRES,IGRP,1)
         ENDIF
   80    CONTINUE
         IF(IND.GT.0) THEN
            SIGTXS(IBM)=SIGTXS(IBM)+CONR(IND,IRES)*TOTPT(INOR,IRES)
         ENDIF
   90    CONTINUE
         IASM=IASM+1
         NPSYS(IASM)=IASM
         KPSYS=LCMDIL(IPSYS,IASM)
         CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBMIX+1,2,SIGTXS)
         CALL LCMPUT(KPSYS,'DRAGON-S0XSC',NBMIX+1,2,SIGS0X)
  110    CONTINUE
      ELSE IF(GOLD(IRES,IGRP).EQ.-999.) THEN
         CALL LCMLEL(JPLI0,IGRP,LENG0,ITYLCM)
         IF(LENG0.NE.0) THEN
            WRITE(HSMG,'(42HUSSIT0: UNEXPECTED SELF-SHIELDING DATA FOU,
     1      11HND IN GROUP,I5,1H.)') IGRP
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
  120 CONTINUE
*----
*  ASSEMBLY MATRIX OR REDUCED COLLISION PROBABILITIES CALCULATION.
*----
      ISTRM=1
      NANI=1
      NW=0
      KNORM=1
      IMPY=MAX(0,IMPX-3)
      IF(IPHASE.EQ.1) THEN
*        USE A NATIVE DOOR.
         CALL DOORAV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPY,NASM,NREG,
     1   NBMIX,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         IPIJK=1
         CALL DOORPV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPY,NASM,NREG,
     1   NBMIX,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,.FALSE.,TITR,NALBP)
      ENDIF
      DEALLOCATE(NPSYS)
*----
*  LOOP OVER THE ENERGY GROUPS.
*----
      IASM=0
      DO 300 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.)) THEN
         IF(IMPX.GT.5) WRITE(6,'(/25H USSIT0: PROCESSING GROUP,I5,
     >   6H WITH ,A,1H.)') IGRP,CDOOR
         NORI=NOR(IRES,IGRP)
*----
*  COMPUTE THE AVERAGED COLLISION PROBABILITY MATRIX.
*----
         ALLOCATE(NPSYS(NORI*(NBNRS+1)))
         ALLOCATE(FUN(NUN*NORI*(NBNRS+1)),SUN(NUN*NORI*(NBNRS+1)))
         CALL XDRSET(FUN,NUN*NORI*(NBNRS+1),0.0)
         CALL XDRSET(SUN,NUN*NORI*(NBNRS+1),0.0)
         DO 142 INOR=1,NORI
         DO 141 JNBN=0,NBNRS
         NPSYS((INOR-1)*(NBNRS+1)+JNBN+1)=IASM+INOR
         T1=0.0
         DO 130 I=1,NREG
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 130
         IOF=(INOR-1)*NUN*(NBNRS+1)+JNBN*NUN+KEYFLX(I)
         IND=IREX(IBM)
         IF((JNBN.EQ.0).AND.(IND.EQ.0)) THEN
            T1=T1+SIGGAR(IBM,0,IGRP,3)*VOL(I)
            SUN(IOF)=SUN(IOF)+SIGGAR(IBM,0,IGRP,3)
         ELSE IF(IND.EQ.JNBN) THEN
            T1=T1+VOL(I)
            SUN(IOF)=1.0
         ENDIF
  130    CONTINUE
         DO 140 I=1,NUN
         IOF=(INOR-1)*NUN*(NBNRS+1)+JNBN*NUN+I
         IF(T1.NE.0.0) SUN(IOF)=SUN(IOF)/T1
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
*----
*  SOLVE FOR THE MULTIBAND FLUX.
*----
         IDIR=0
         NABS=NORI*(NBNRS+1)
         LEXAC=.FALSE.
         CALL DOORFV (CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NABS,NBMIX,
     1   IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN)
*----
*  HOMOGENIZE THE MULTIBAND FLUX.
*----
         DO 170 INOR=1,NORI
         CALL XDRSET(PAV,(NBNRS+1)**2,0.0)
         DO 155 JNBN=0,NBNRS
         T1=0.0
         DO 150 I=1,NREG
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 150
         IOF=(INOR-1)*NUN*(NBNRS+1)+JNBN*NUN+KEYFLX(I)
         PAV(IREX(IBM),JNBN)=PAV(IREX(IBM),JNBN)+FUN(IOF)*VOL(I)
  150    CONTINUE
  155    CONTINUE
         DO 165 I=0,NBNRS
         DO 160 J=0,NBNRS
         IF(VOLMER(I).NE.0.0) PAV(I,J)=PAV(I,J)*VOLMER(J)/VOLMER(I)
  160    CONTINUE
  165    CONTINUE
         KPSYS=LCMGIL(IPSYS,IASM+INOR)
         CALL LCMPUT(KPSYS,'DRAGON-PAV',(NBNRS+1)**2,2,PAV(0,0))
  170    CONTINUE
         DEALLOCATE(SUN,FUN,NPSYS)
*----
*  COLLECT THE BASE POINTS IN TOTAL AND PARTIAL CROSS SECTION.
*----
         DO 200 JRES=1,NIRES
         CALL LCMLEL(IPPT1(JRES),IGRP,ILONG,ITYLCM)
         IF(ILONG.NE.0) THEN
            KPLIB=LCMGIL(IPPT1(JRES),IGRP)
            CALL LCMINF(KPLIB,TEXT12,TEXX12,EMPTY,ILONG,LCM)
            CALL LCMLEN(KPLIB,'PROB-TABLE',LENG,ITYLCM)
            NPART=LENG/MAXNOR
            IF(LCM) THEN
               CALL LCMGPD(KPLIB,'PROB-TABLE',SIGP_PTR)
               CALL C_F_POINTER(SIGP_PTR,SIGP,(/ MAXNOR*NPART /))
            ELSE
               ALLOCATE(SIGP(MAXNOR*NPART))
               CALL LCMGET(KPLIB,'PROB-TABLE',SIGP)
            ENDIF
            IF(GOLD(IRES,IGRP).EQ.-999.) THEN
               DO 180 INOR=1,NOR(JRES,IGRP)
               WEIGH(INOR,JRES)=SIGP(INOR)
               TOTPT(INOR,JRES)=SIGP(MAXNOR+INOR)
  180          CONTINUE
               CALL LCMGET(KPLIB,'SIGQT-SLOW',WSLD(1,JRES))
               CALL LCMGET(KPLIB,'SIGQT-SIGS',SIGWS(1,JRES))
            ELSE
               DO 190 INOR=1,NOR(JRES,IGRP)
               WEIGH(INOR,JRES)=SIGP(INOR)
               TOTPT(INOR,JRES)=SIGP(MAXNOR+INOR)
               SIGWS(INOR,JRES)=SIGP(3*MAXNOR+INOR)
  190          CONTINUE
            ENDIF
            IF(.NOT.LCM) DEALLOCATE(SIGP)
         ELSE
            WEIGH(1,JRES)=1.0
            TOTPT(1,JRES)=STGAR(IPPT2(JRES,1),JRES,IGRP)
            IF(GOLD(IRES,IGRP).EQ.-999.) THEN
               SIGWS(1,JRES)=SWGAR(IPPT2(JRES,1),JRES,IGRP)
               WSLD(1,JRES)=1.0
            ELSE
               SIGWS(1,JRES)=SSGAR(IPPT2(JRES,1),JRES,1,IGRP)
            ENDIF
         ENDIF
  200    CONTINUE
*----
*  TAKE INTO ACCOUNT CORRELATION EFFECTS BETWEEN ISOTOPES USING THE
*  MUTUAL SELF-SHIELDING MODEL.
*----
        IF((NIRES.GT.1).AND.(GOLD(IRES,IGRP).EQ.-999.).AND.
     1  (ICORR.EQ.0)) THEN
           DO 225 JRES=1,NIRES
           DO 220 IND=1,NBNRS
           SIGX(IND,JRES)=0.0
           T1=0.0
           T2=0.0
           DO 215 I=1,NREG
           IBM=MAT(I)
           IF(IBM.EQ.0) GO TO 215
           IF(IND.EQ.IREX(IBM)) THEN
              T1=T1+(SIGGAR(IBM,JRES,IGRP,1)-SIGGAR(IBM,JRES,IGRP,2))*
     1        VOL(I)
              T2=T2+VOL(I)
           ENDIF
  215      CONTINUE
           IF(T2.NE.0.0) SIGX(IND,JRES)=T1/T2
  220      CONTINUE
  225      CONTINUE
           CALL USSCOR(MAXNOR,IGRP,IPSYS,IASM,IRES,NBNRS,NIRES,
     1     NOR(1,IGRP),CONR,IPPT1,IPPT2,WEIGH,TOTPT,SIGX,VOLMER)
        ENDIF
*----
*  RESPONSE MATRIX APPROACH. LOOP OVER THE SECONDARY SUBGROUPS.
*----
        DO 272 INOR=1,NORI
        KPSYS=LCMGIL(IPSYS,IASM+INOR)
        CALL LCMGET(KPSYS,'DRAGON-PAV',PAV(0,0))
*----
*  LOOP OVER THE PRIMARY SUBGROUPS. NORI+1 IS THE SOURCE.
*----
        DO 271 JNOR=1,NORI+1
        IF(JNOR.LE.NORI) THEN
           JNBMAX=NBNRS
        ELSE
           JNBMAX=1
        ENDIF
        DO 270 JNBN=1,JNBMAX
        CALL XDRSET(AWPHI(1),NBNRS,0.0)
        DO 250 I=1,NREG
        IBM=MAT(I)
        IF(IBM.EQ.0) GO TO 250
        JND=IREX(IBM)
        QQQ=0.0
        IF(JNOR.EQ.NORI+1) THEN
           DO 230 JRES=0,NIRES
           IF(JRES.EQ.0) THEN
              QQQ=QQQ+SIGGAR(IBM,0,IGRP,3)
           ELSE IF((JRES.NE.IRES).AND.(JND.GT.0)) THEN
              QQQ=QQQ+SIGGAR(IBM,JRES,IGRP,4)
           ENDIF
  230      CONTINUE
        ELSE IF(JND.EQ.JNBN) THEN
           IF(GOLD(IRES,IGRP).EQ.-999.) THEN
              WWW=WSLD((JNOR-1)*NORI+INOR,IRES)/WEIGH(INOR,IRES)
           ELSE
              WWW=WEIGH(JNOR,IRES)
           ENDIF
           QQQ=QQQ-WWW*CONR(JND,IRES)*SIGWS(JNOR,IRES)
        ENDIF
        DO 240 IND=1,NBNRS
        AWPHI(IND)=AWPHI(IND)+PAV(IND,JND)*QQQ*VOL(I)/VOLMER(JND)
  240   CONTINUE
  250   CONTINUE
        DO 260 IND=1,NBNRS
        MATRIX(INM(IND,INOR,NBNRS),INM(JNBN,JNOR,NBNRS))=AWPHI(IND)
  260   CONTINUE
  270   CONTINUE
  271   CONTINUE
  272   CONTINUE
*
        DO 280 I=1,NBNRS*NORI
        MATRIX(I,I)=MATRIX(I,I)+1.0D0
  280   CONTINUE
        CALL ALSBD(NBNRS*NORI,1,MATRIX,IER,NBNRS*MAXNOR)
        IF(IER.NE.0) CALL XABORT('USSIT0: SINGULAR MATRIX.')
        CALL XDRSET(XFLUX(1,1,1,IRES),NBNRS*MAXNOR*NL,0.0)
        DO 295 IND=1,NBNRS
        DO 290 INOR=1,NORI
        I1=INM(IND,INOR,NBNRS)
        XFLUX(IND,INOR,1,IRES)=REAL(MATRIX(I1,NBNRS*NORI+1))
  290   CONTINUE
  295   CONTINUE
* END OF RESPONSE MATRIX APPROACH.
*
        KPLI0=LCMGIL(JPLI0,IGRP)
        CALL LCMPUT(KPLI0,'NWT0-PT',NBNRS*MAXNOR*NL,2,XFLUX(1,1,1,IRES))
        IASM=IASM+NORI
      ENDIF
  300 CONTINUE
*----
*  COMPUTE UNGAR, THE REGION-ORDERED FLUX.
*----
      ALLOCATE(NPSYS(NASM),FUN(NUN*NASM),SUN(NUN*NASM))
      CALL XDRSET(SUN,NUN*NASM,0.0)
      IASM=0
      DO 420 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.)) THEN
         NORI=NOR(IRES,IGRP)
*----
*  RECOVER THE PREVIOUS FLUXES.
*----
         WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
         CALL LCMSIX(IPLI0,CBDPNM,1)
         JPLI1=LCMGID(IPLI0,'GROUP_INFO')
         KPLI1=LCMGIL(JPLI1,IGRP)
         CALL LCMLEN(KPLI1,'NWT0-PT',ILON,ITYLCM)
         IF(ILON.GT.NBNRS*MAXNOR*NL) THEN
            WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
            WRITE(HSMG,'(34HUSSIT0: FLUX OVERFLOW FOR ISOTOPE ,A12)')
     1      TEXT12
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGET(KPLI1,'NWT0-PT',XFLUX(1,1,1,IRES))
         CALL LCMSIX(IPLI0,' ',2)
*----
*  COLLECT THE BASE POINTS IN PARTIAL CROSS SECTION.
*----
         DO 340 JRES=1,NIRES
         CALL LCMLEL(IPPT1(JRES),IGRP,ILONG,ITYLCM)
         IF(ILONG.NE.0) THEN
            KPLIB=LCMGIL(IPPT1(JRES),IGRP)
            CALL LCMINF(KPLIB,TEXT12,TEXX12,EMPTY,ILONG,LCM)
            CALL LCMLEN(KPLIB,'PROB-TABLE',LENG,ITYLCM)
            NPART=LENG/MAXNOR
            IF(LCM) THEN
               CALL LCMGPD(KPLIB,'PROB-TABLE',SIGP_PTR)
               CALL C_F_POINTER(SIGP_PTR,SIGP,(/ MAXNOR*NPART /))
            ELSE
               ALLOCATE(SIGP(MAXNOR*NPART))
               CALL LCMGET(KPLIB,'PROB-TABLE',SIGP)
            ENDIF
            IF(GOLD(IRES,IGRP).EQ.-999.) THEN
               DO 320 INOR=1,NOR(JRES,IGRP)
               WEIGH(INOR,JRES)=SIGP(INOR)
  320          CONTINUE
               CALL LCMGET(KPLIB,'SIGQT-SLOW',WSLD(1,JRES))
               CALL LCMGET(KPLIB,'SIGQT-SIGS',SIGWS(1,JRES))
            ELSE
               DO 330 INOR=1,NOR(JRES,IGRP)
               WEIGH(INOR,JRES)=SIGP(INOR)
               SIGWS(INOR,JRES)=SIGP(3*MAXNOR+INOR)
  330          CONTINUE
            ENDIF
            IF(.NOT.LCM) DEALLOCATE(SIGP)
         ELSE
            WEIGH(1,JRES)=1.0
            IF(GOLD(IRES,IGRP).EQ.-999.) THEN
               SIGWS(1,JRES)=SWGAR(IPPT2(JRES,1),JRES,IGRP)
               WSLD(1,JRES)=1.0
            ELSE
               SIGWS(1,JRES)=SSGAR(IPPT2(JRES,1),JRES,1,IGRP)
            ENDIF
         ENDIF
  340    CONTINUE
*----
*  COMPUTE THE AVERAGED SOURCE.
*----
         DO 385 INOR=1,NORI
         NPSYS(IASM+INOR)=IASM+INOR
         KPSYS=LCMGIL(IPSYS,IASM+INOR)
         CALL LCMLEN(KPSYS,'FUNKNO$USS',ILENG,ITYLCM)
         IF(ILENG.EQ.NUN) THEN
            CALL LCMGET(KPSYS,'FUNKNO$USS',FUN((IASM+INOR-1)*NUN+1))
         ELSE
            CALL XDRSET(FUN((IASM+INOR-1)*NUN+1),NUN,0.0)
         ENDIF
         IOF=(IASM+INOR-1)*NUN
         DO 380 I=1,NREG
         IBM=MAT(I)
         QQQ=0.0
         IF(IBM.EQ.0) GO TO 375
         IND=IREX(IBM)
         DO 360 JRES=0,NIRES
         IF(JRES.EQ.0) THEN
            QQQ=QQQ+SIGGAR(IBM,0,IGRP,3)
         ELSE IF((JRES.NE.IRES).AND.(IND.GT.0)) THEN
            QQQ=QQQ+SIGGAR(IBM,JRES,IGRP,4)
         ENDIF
  360    CONTINUE
         IF(IND.GT.0) THEN
            DO 370 JNOR=1,NORI
            IF(GOLD(IRES,IGRP).EQ.-999.) THEN
               WWW=WSLD((JNOR-1)*NORI+INOR,IRES)/WEIGH(INOR,IRES)
            ELSE
               WWW=WEIGH(JNOR,IRES)
            ENDIF
            QQQ=QQQ+WWW*CONR(IND,IRES)*SIGWS(JNOR,IRES)*
     1      XFLUX(IND,JNOR,1,IRES)
  370       CONTINUE
         ENDIF
  375    SUN(IOF+KEYFLX(I))=QQQ*WEIGH(INOR,IRES)
  380    CONTINUE
  385    CONTINUE
*
         IF(IMPX.GT.0) THEN
            WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
            WRITE(6,'(15H USSIT0: GROUP=,I5,24H. SUBGROUP CALCULATION B,
     1      37HASED ON RESPONSE MATRICES.  ISOTOPE='',A12,2H''.)') IGRP,
     2      TEXT12
         ENDIF
         IF(IMPX.GT.2) THEN
            DO 400 IND=1,NBNRS
            T1=0.0
            DO 390 INOR=1,NOR(IRES,IGRP)
            T1=T1+WEIGH(INOR,IRES)*XFLUX(IND,INOR,1,IRES)
  390       CONTINUE
            WRITE(6,'(31H USSIT0: AVERAGED FLUX IN GROUP,I4,8H AND RES,
     1      12HONANT REGION,I4,21H FOR RESONANT ISOTOPE,I4,2H =,F9.5)')
     2      IGRP,IND,IRES,T1
  400       CONTINUE
         ENDIF
*
         IASM=IASM+NORI
      ENDIF
  420 CONTINUE
*----
*  SOLVE FOR THE MULTIBAND FLUX (VECTOR OF LENGTH NREG).
*----
      IDIR=0
      LEXAC=.FALSE.
      IF(IMPX.GT.5) WRITE(6,'(/33H USSIT0: PROCESSING MULTIBAND FLU,
     1 14HX (IL=1) WITH ,A,1H.)') CDOOR
      CALL DOORFV (CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NASM,NBMIX,
     1 IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN)
*----
*  INTEGRATE THE REGION-ORDERED FLUX OVER SUBGROUPS.
*----
      IASM=0
      DO 480 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.)) THEN
        CALL XDRSET(UNGAR(1,IRES,IGRP),NREG,0.0)
        NORI=NOR(IRES,IGRP)
        DO 475 INOR=1,NORI
        KPSYS=LCMGIL(IPSYS,IASM+INOR)
        IOF=(IASM+INOR-1)*NUN
        CALL LCMPUT(KPSYS,'FUNKNO$USS',NUN,2,FUN(IOF+1))
*----
*  NORMALIZE THE MULTIBAND FLUX. THIS NORMALIZATION IS ONLY REQUIRED IF
*  THE MUTUAL SELF-SHIELDING MODEL IS USED.
*----
        IF((NIRES.GT.1).AND.(GOLD(IRES,IGRP).EQ.-999.).AND.(ICORR.EQ.0))
     1  THEN
           IOFF=(IASM+INOR-1)*NUN
           CALL XDRSET(AWPHI(0),NBNRS+1,0.0)
           DO 430 I=1,NREG
           IBM=MAT(I)
           IF(IBM.GT.0) THEN
              IND=IREX(IBM)
              AWPHI(IND)=AWPHI(IND)+FUN(IOFF+KEYFLX(I))*VOL(I)/
     1        VOLMER(IND)
           ENDIF
  430      CONTINUE
           CALL LCMGET(KPSYS,'DRAGON-PAV',PAV(0,0))
           DO 450 IND=0,NBNRS
           TT=0.0
           DO 440 J=1,NREG
           IBM=MAT(J)
           IF(IBM.GT.0) THEN
              JND=IREX(IBM)
              IOFS=(IASM+INOR-1)*NUN+KEYFLX(J)
              TT=TT+PAV(IND,JND)*SUN(IOFS)*VOL(J)/VOLMER(JND)
           ENDIF
  440      CONTINUE
           AWPHI(IND)=TT/AWPHI(IND)
  450      CONTINUE
           DO 460 I=1,NREG
           IBM=MAT(I)
           IF(IBM.GT.0) FUN(IOFF+KEYFLX(I))=FUN(IOFF+KEYFLX(I))*
     1     AWPHI(IREX(IBM))
  460      CONTINUE
        ENDIF
*
        DO 470 I=1,NREG
        IOF=(IASM+INOR-1)*NUN+KEYFLX(I)
        UNGAR(I,IRES,IGRP)=UNGAR(I,IRES,IGRP)+FUN(IOF)
  470   CONTINUE
  475   CONTINUE
        IASM=IASM+NORI
      ENDIF
  480 CONTINUE
*----
*  USE IL=1 VALUES FOR XFLUX AT HIGHER LEGENDRE ORDERS.
*----
      DO 540 IL=2,NL
      DO 530 IGRP=1,NGRP
      IF(MASKG(IGRP).AND.(GOLD(IRES,IGRP).EQ.-999.)) THEN
        NORI=NOR(IRES,IGRP)
        KPLI0=LCMGIL(JPLI0,IGRP)
        CALL LCMLEN(KPLI0,'NWT0-PT',ILON,ITYLCM)
        IF(ILON.GT.NBNRS*MAXNOR*NL) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
           WRITE(HSMG,'(34HUSSIT0: FLUX OVERFLOW FOR ISOTOPE ,A12)')
     1     TEXT12
           CALL XABORT(HSMG)
        ENDIF
        CALL LCMGET(KPLI0,'NWT0-PT',XFLUX(1,1,1,IRES))
        DO 525 INOR=1,NORI
        DO 520 IND=1,NBNRS
        XFLUX(IND,INOR,IL,IRES)=XFLUX(IND,INOR,1,IRES)
  520   CONTINUE
  525   CONTINUE
        CALL LCMPUT(KPLI0,'NWT0-PT',NBNRS*MAXNOR*NL,2,XFLUX(1,1,1,IRES))
      ENDIF
  530 CONTINUE
  540 CONTINUE
      DEALLOCATE(SUN,FUN,NPSYS)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MATRIX)
      DEALLOCATE(SIGX,PAV,SIGWS,WSLD,TOTPT,WEIGH,AWPHI,SIGS0X,SIGTXS,
     1 XFLUX)
      RETURN
      END
