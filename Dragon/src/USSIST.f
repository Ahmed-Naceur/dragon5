*DECK USSIST
      SUBROUTINE USSIST(MAXNOR,NGRP,MASKG,IRES,IPLI0,IPTRK,IFTRAK,
     1 CDOOR,IMPX,NBMIX,NREG,NUN,NL,IPHASE,MAXST,MAT,VOL,KEYFLX,LEAKSW,
     2 IREX,SIGGAR,TITR,ICORR,NIRES,NBNRS,NOR,CONR,GOLD,IPPT1,IPPT2,
     3 STGAR,SSGAR,VOLMER,UNGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the multiband fluxes as required by the subgroup method using
* the subgroup projection method (SPM):
* a) assume a single resonant isotope;
* b) use the standard solution doors of Dragon.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
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
* MAXST   maximum number of fixed point iterations for the ST scattering
*         source.
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
*         tables; =-998.0/-1000.0 subgroup projection method).
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope.
* STGAR   averaged microscopic total xs in resonant region.
* SSGAR   averaged microscopic scattering xs in resonant region.
* VOLMER  volumes of the resonant regions.
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
      INTEGER MAXNOR,NGRP,IRES,IFTRAK,IMPX,NBMIX,NREG,NUN,NL,IPHASE,
     1 MAXST,MAT(NREG),KEYFLX(NREG),IREX(NBMIX),ICORR,NIRES,NBNRS,
     2 NOR(NIRES,NGRP),IPPT2(NIRES,4)
      REAL VOL(NREG),SIGGAR(NBMIX,0:NIRES,NGRP,4),
     1 CONR(NBNRS,NIRES),GOLD(NIRES,NGRP),STGAR(NBNRS,NIRES,NGRP),
     2 SSGAR(NBNRS,NIRES,NL,NGRP),VOLMER(0:NBNRS),
     3 UNGAR(NREG,NIRES,NGRP)
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW,MASKG(NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLI0,IPSYS,JPLI1,KPLI1,KPLIB,KPSYS,KPLI0
      CHARACTER CBDPNM*12,TEXT12*12,TEXX12*12,HSMG*131
      DOUBLE PRECISION ZNUM,ZDEN
      LOGICAL EMPTY,LCM,LEXAC,LSPM
      INTEGER NALBP
*----
*  ALLOCATABLE STATEMENTS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGTXS,SIGS0X,FLNEW,WCOR,FUN,
     1 SUN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WEIGH,TOTPT,SIGWS
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: XFLUX
      TYPE(C_PTR) SIGP_PTR
      REAL, POINTER, DIMENSION(:) :: SIGP
*----
*  FIND THE NUMBER OF COMPONENTS REQUIRED AND ALLOCATE THE LIST OF
*  ASSEMBLY MATRICES.
*----
      NASM=0
      NALBP=0
      DO 10 IGRP=1,NGRP
      LSPM=(MASKG(IGRP).AND.((GOLD(IRES,IGRP).EQ.-998.).OR.
     1 (GOLD(IRES,IGRP).EQ.-1000.)))
      IF(LSPM) NASM=NASM+NOR(IRES,IGRP)
   10 CONTINUE
      IF(NASM.EQ.0) RETURN
      IF(NGRP.LT.250) CALL XABORT('USSIST: THIS SIMPLIFIED SELF-SHIELD'
     1 //'ING MODEL REQUIRES MORE THAN 250 ENERGY GROUPS.')
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NPSYS(MAXNOR*NGRP))
      ALLOCATE(XFLUX(NBNRS,MAXNOR,NL,NIRES),SIGTXS(0:NBMIX),
     1 SIGS0X(0:NBMIX),WEIGH(MAXNOR,NIRES),TOTPT(MAXNOR,NIRES),
     2 SIGWS(MAXNOR,NIRES),FLNEW(NBNRS),WCOR(MAXNOR**2))
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
      IASM=0
      DO 100 IGRP=1,NGRP
      LSPM=(MASKG(IGRP).AND.((GOLD(IRES,IGRP).EQ.-998.).OR.
     1 (GOLD(IRES,IGRP).EQ.-1000.)))
      IF(LSPM) THEN
         IF(IMPX.GT.1) THEN
            WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
            WRITE(6,'(36H USSIST: PROCESS CORRELATED ISOTOPE ,A12,
     1      11H WITH INDEX,I3,9H IN GROUP,I4,22H (SUBGROUP PROJECTION ,
     2      8HMETHOD).)') TEXT12,IRES,IGRP
         ENDIF
         DO 20 JRES=1,NIRES
         IF(GOLD(JRES,IGRP).NE.GOLD(IRES,IGRP)) THEN
            WRITE(HSMG,'(32HUSSIST: PT NOT SET FOR ISOTOPE '',3A4,
     1      10H'' IN GROUP,I5,1H.)') (IPPT2(JRES,J0),J0=2,4),IGRP
            CALL XABORT(HSMG)
         ELSE IF(NOR(JRES,IGRP).GT.MAXNOR) THEN
            CALL XABORT('USSIST: MAXNOR OVERFLOW.')
         ENDIF
   20    CONTINUE
         NORI=NOR(IRES,IGRP)
         DO 40 JRES=1,NIRES
*----
*  RECOVER THE PREVIOUS FLUXES.
*----
         IF(NOR(JRES,IGRP).EQ.1) THEN
            CALL XDRSET(XFLUX(1,1,1,JRES),NBNRS*MAXNOR*NL,1.0)
         ELSE
            WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') JRES,NIRES
            CALL LCMSIX(IPLI0,CBDPNM,1)
            JPLI1=LCMGID(IPLI0,'GROUP_INFO')
            KPLI1=LCMGIL(JPLI1,IGRP)
            CALL LCMLEN(KPLI1,'NWT0-PT',ILON,ITYLCM)
            IF(ILON.GT.NBNRS*MAXNOR*NL) THEN
               WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
               WRITE(HSMG,'(34HUSSIST: FLUX OVERFLOW FOR ISOTOPE ,A12)')
     1         TEXT12
               CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(KPLI1,'NWT0-PT',XFLUX(1,1,1,JRES))
            CALL LCMSIX(IPLI0,' ',2)
         ENDIF
*----
*  COLLECT THE BASE POINTS IN TOTAL CROSS SECTION.
*----
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
            SIGWS(INOR,JRES)=SIGP(3*MAXNOR+INOR)
   30       CONTINUE
            IF(.NOT.LCM) DEALLOCATE(SIGP)
         ELSE
            WEIGH(1,JRES)=1.0
            TOTPT(1,JRES)=STGAR(IPPT2(JRES,1),JRES,IGRP)
            SIGWS(1,JRES)=SSGAR(IPPT2(JRES,1),JRES,1,IGRP)
         ENDIF
   40    CONTINUE
*----
*  SET THE MIXTURE-DEPENDENT CROSS SECTIONS.
*----
         KPLIB=LCMGIL(IPPT1(IRES),IGRP)
         DO 90 INOR=1,NORI
         CALL XDRSET(SIGTXS(0),NBMIX+1,0.0)
         CALL XDRSET(SIGS0X(0),NBMIX+1,0.0)
         DO 80 IBM=1,NBMIX
         IND=IREX(IBM)
         DO 70 JRES=0,NIRES
         IF(JRES.EQ.0) THEN
*           ADMIXED NON-RESONANT ISOTOPES.
            SIGTXS(IBM)=SIGTXS(IBM)+(SIGGAR(IBM,0,IGRP,1)-
     1      SIGGAR(IBM,0,IGRP,2))
            SIGS0X(IBM)=SIGS0X(IBM)-SIGGAR(IBM,0,IGRP,2)
         ELSE IF((JRES.NE.IRES).AND.(IND.GT.0)) THEN
*           MUTUAL SHIELDING MODEL FROM CORRELATED RESONANT ISOTOPES.
            WRITE(TEXT12,'(3A4)') (IPPT2(JRES,I0),I0=2,4)
            IF((NOR(JRES,IGRP).GT.1).AND.(ICORR.NE.1)) THEN
               CALL LCMLEN(KPLIB,TEXT12,ILEN,ITYLCM)
               IF(ILEN.EQ.0) THEN
                  CALL LCMLIB(KPLIB)
                  CALL XABORT('USSIST: UNABLE TO FIND CORRELATED ISOTO'
     1            //'PE '//TEXT12//'.')
               ENDIF
               CALL LCMGET(KPLIB,TEXT12,WCOR)
            ENDIF
            IF((ICORR.EQ.1).AND.
     1        (IPPT2(IRES,2).EQ.IPPT2(JRES,2)).AND.
     1        (IPPT2(IRES,3).EQ.IPPT2(JRES,3))) THEN
              SIGTXS(IBM)=SIGTXS(IBM)+CONR(IND,JRES)*TOTPT(INOR,IRES)
            ELSE
              ZNUM=0.0D0
              ZDEN=0.0D0
              DO 60 I2=1,NOR(JRES,IGRP)
              IF((ICORR.EQ.1).OR.(NOR(JRES,IGRP).EQ.1)) THEN
                WWW=WEIGH(INOR,IRES)*WEIGH(I2,JRES)*XFLUX(IND,I2,1,JRES)
              ELSE
                WWW=WCOR((I2-1)*NORI+INOR)*XFLUX(IND,I2,1,JRES)
              ENDIF
              ZNUM=ZNUM+WWW*CONR(IND,JRES)*TOTPT(I2,JRES)
              ZDEN=ZDEN+WWW
   60         CONTINUE
              IF(ZNUM/ZDEN.LT.0.0) THEN
*               BADLY BEHAVED CORRELATED WEIGHT MATRIX.
                ZNUM=0.0D0
                ZDEN=0.0D0
                DO 65 I2=1,NOR(JRES,IGRP)
                WWW=WEIGH(INOR,IRES)*WEIGH(I2,JRES)*XFLUX(IND,I2,1,JRES)
                ZNUM=ZNUM+WWW*CONR(IND,JRES)*TOTPT(I2,JRES)
                ZDEN=ZDEN+WWW
   65           CONTINUE
              ENDIF
              SIGTXS(IBM)=SIGTXS(IBM)+REAL(ZNUM/ZDEN)
            ENDIF
         ENDIF
   70    CONTINUE
         IF(IND.GT.0) THEN
            SIGTXS(IBM)=SIGTXS(IBM)+CONR(IND,IRES)*TOTPT(INOR,IRES)
            SIGS0X(IBM)=SIGS0X(IBM)+WEIGH(INOR,IRES)*CONR(IND,IRES)*
     1      SIGWS(INOR,IRES)
         ENDIF
   80    CONTINUE
         IASM=IASM+1
         NPSYS(IASM)=IASM
         KPSYS=LCMDIL(IPSYS,IASM)
         CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBMIX+1,2,SIGTXS(0))
         CALL LCMPUT(KPSYS,'DRAGON-S0XSC',NBMIX+1,2,SIGS0X(0))
   90    CONTINUE
      ELSE IF(GOLD(IRES,IGRP).EQ.-998.) THEN
         CALL LCMLEL(JPLI0,IGRP,LENG0,ITYLCM)
         IF(LENG0.NE.0) THEN
            WRITE(HSMG,'(42HUSSIST: UNEXPECTED SELF-SHIELDING DATA FOU,
     1      11HND IN GROUP,I5,1H.)') IGRP
            CALL XABORT(HSMG)
         ENDIF
      ELSE IF(GOLD(IRES,IGRP).EQ.-1000.) THEN
         CALL LCMLEL(JPLI0,IGRP,LENG0,ITYLCM)
         IF(LENG0.NE.0) THEN
            WRITE(HSMG,'(42HUSSIST: UNEXPECTED SELF-SHIELDING DATA FOU,
     1      11HND IN GROUP,I5,1H.)') IGRP
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
  100 CONTINUE
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
*----
*  LOOP OVER THE ENERGY GROUPS.
*----
      IASM=0
      DO 280 IGRP=1,NGRP
      LSPM=(MASKG(IGRP).AND.((GOLD(IRES,IGRP).EQ.-998.).OR.
     1 (GOLD(IRES,IGRP).EQ.-1000.)))
      IF(LSPM) THEN
         NORI=NOR(IRES,IGRP)
*----
*  COLLECT THE BASE POINTS IN PARTIAL CROSS SECTION.
*----
         DO 120 JRES=1,NIRES
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
            DO 110 INOR=1,NOR(JRES,IGRP)
            WEIGH(INOR,JRES)=SIGP(INOR)
            SIGWS(INOR,JRES)=SIGP(3*MAXNOR+INOR)
  110       CONTINUE
            IF(.NOT.LCM) DEALLOCATE(SIGP)
         ELSE
            WEIGH(1,JRES)=1.0
            SIGWS(1,JRES)=SSGAR(IPPT2(JRES,1),JRES,1,IGRP)
         ENDIF
  120    CONTINUE
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
            WRITE(HSMG,'(34HUSSIST: FLUX OVERFLOW FOR ISOTOPE ,A12)')
     1      TEXT12
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGET(KPLI1,'NWT0-PT',XFLUX(1,1,1,IRES))
         CALL LCMSIX(IPLI0,' ',2)
*----
*  ITERATIVE PROCEDURE.
*----
        ITER=0
  140   ITER=ITER+1
        IF(ITER.GT.MAXST) GO TO 240
        ERR1=0.0
        ERR2=0.0
        CALL XDRSET(UNGAR(1,IRES,IGRP),NREG,0.0)
*----
*  COMPUTE THE AVERAGED SOURCE.
*----
        ALLOCATE(FUN(NUN*NORI),SUN(NUN*NORI))
        CALL XDRSET(SUN,NUN*NORI,0.0)
        DO 195 INOR=1,NORI
        KPSYS=LCMGIL(IPSYS,IASM+INOR)
        CALL LCMLEN(KPSYS,'FUNKNO$USS',ILENG,ITYLCM)
        IF(ILENG.EQ.NUN) THEN
           CALL LCMGET(KPSYS,'FUNKNO$USS',FUN((INOR-1)*NUN+1))
        ELSE
           CALL XDRSET(FUN((INOR-1)*NUN+1),NUN,0.0)
        ENDIF
        NPSYS(INOR)=IASM+INOR
        DO 190 I=1,NREG
        IBM=MAT(I)
        IF(IBM.EQ.0) GO TO 190
        IOF=(INOR-1)*NUN+KEYFLX(I)
        IND=IREX(IBM)
        DO 150 JRES=0,NIRES
        IF(JRES.EQ.0) THEN
           SUN(IOF)=SUN(IOF)+SIGGAR(IBM,0,IGRP,3)
        ELSE IF((JRES.NE.IRES).AND.(IND.GT.0)) THEN
           SUN(IOF)=SUN(IOF)+SIGGAR(IBM,JRES,IGRP,4)
        ENDIF
  150   CONTINUE
        IF(IND.GT.0) THEN
           DO 160 JNOR=1,NORI
           IF(JNOR.NE.INOR) THEN
              SUN(IOF)=SUN(IOF)+WEIGH(JNOR,IRES)*CONR(IND,IRES)*
     1        SIGWS(JNOR,IRES)*XFLUX(IND,JNOR,1,IRES)
           ENDIF
  160      CONTINUE
        ENDIF
  190   CONTINUE
  195   CONTINUE
*----
*  SOLVE FOR THE MULTIBAND FLUX.
*----
        IDIR=0
        LEXAC=.FALSE.
        CALL DOORFV (CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NORI,NBMIX,
     1  IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN)
*----
*  HOMOGENIZE THE FLUX AT ITERATION ITER.
*----
        DO 235 INOR=1,NORI
        KPSYS=LCMGIL(IPSYS,IASM+INOR)
        CALL LCMPUT(KPSYS,'FUNKNO$USS',NUN,2,FUN((INOR-1)*NUN+1))
        CALL XDRSET(FLNEW,NBNRS,0.0)
        DO 200 I=1,NREG
        IF(MAT(I).EQ.0) GO TO 200
        IOF=(INOR-1)*NUN+KEYFLX(I)
        IND=IREX(MAT(I))
        IF(IND.GT.0) FLNEW(IND)=FLNEW(IND)+FUN(IOF)*VOL(I)
  200   CONTINUE
        DO 210 IND=1,NBNRS
        FLNEW(IND)=FLNEW(IND)/VOLMER(IND)
  210   CONTINUE
*
        DO 220 I=1,NREG
        IOF=(INOR-1)*NUN+KEYFLX(I)
        UNGAR(I,IRES,IGRP)=UNGAR(I,IRES,IGRP)+FUN(IOF)*WEIGH(INOR,IRES)
  220   CONTINUE
*----
*  CONVERGENCE CONTROL.
*----
        DO 230 IND=1,NBNRS
        ERR1=MAX(ERR1,ABS(FLNEW(IND)-XFLUX(IND,INOR,1,IRES)))
        ERR2=MAX(ERR2,ABS(FLNEW(IND)))
        XFLUX(IND,INOR,1,IRES)=FLNEW(IND)
  230   CONTINUE
  235   CONTINUE
        DEALLOCATE(SUN,FUN)
        IF(IMPX.GT.2) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
           WRITE(6,'(15H USSIST: GROUP=,I5,20H. SUBGROUP ITERATION,I4,
     1     11H. ISOTOPE='',A12,9H''. ERROR=,1P,E11.4,1H.)')
     2     IGRP,ITER,TEXT12,ERR1
        ENDIF
        IF(ERR2.GT.1.0E10) GO TO 240
        IF(ERR1.GT.1.0E-4*ERR2) GO TO 140
        IF(IMPX.GT.0) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
           WRITE(6,'(15H USSIST: GROUP=,I5,24H. SUBGROUP ITERATION CON,
     1     11HVERGENCE IN,I4,22H ITERATIONS. ISOTOPE='',A12,2H''.)')
     2     IGRP,ITER,TEXT12
        ENDIF
        GO TO 250
*----
*  ALTERNATIVE TREATMENT IN CASE OF FAILURE OF FIXED POINT ITERATIONS.
*  USE A NON-ITERATIVE RESPONSE MATRIX APPROACH.
*----
  240   IF(IMPX.GT.0) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
           WRITE(6,'(15H USSIST: GROUP=,I5,24H. SUBGROUP ITERATION FAI,
     1     16HLED FOR ISOTOPE ,A12,32H. USE AN ALTERNATIVE RESPONSE MA,
     2     14HTRIX APPROACH.)') IGRP,TEXT12
        ENDIF
        CALL USSEXC(MAXNOR,CDOOR,IPLI0,IPTRK,IFTRAK,IMPX,NGRP,IGRP,
     1  IASM,NBMIX,NREG,NUN,NL,IPHASE,MAT,VOL,KEYFLX,IREX,SIGGAR,TITR,
     2  NIRES,IRES,NBNRS,NOR,CONR,IPPT1,IPPT2,STGAR,SSGAR,VOLMER,XFLUX,
     3  UNGAR)
  250   IF(IMPX.GT.2) THEN
           DO 270 IND=1,NBNRS
           T1=0.0
           DO 260 INOR=1,NORI
           T1=T1+WEIGH(INOR,IRES)*XFLUX(IND,INOR,1,IRES)
  260      CONTINUE
           WRITE(6,'(31H USSIST: AVERAGED FLUX IN GROUP,I4,9H AND RESO,
     1     11HNANT REGION,I4,21H FOR RESONANT ISOTOPE,I4,2H =,F9.5)')
     2     IGRP,IND,IRES,T1
  270      CONTINUE
        ENDIF
        KPLI0=LCMGIL(JPLI0,IGRP)
        CALL LCMPUT(KPLI0,'NWT0-PT',NBNRS*MAXNOR*NL,2,XFLUX(1,1,1,IRES))
        IASM=IASM+NORI
      ENDIF
  280 CONTINUE
*----
*  USE IL=1 VALUES FOR XFLUX AT HIGHER LEGENDRE ORDERS.
*----
      DO 310 IL=2,NL
      DO 300 IGRP=1,NGRP
      LSPM=(MASKG(IGRP).AND.((GOLD(IRES,IGRP).EQ.-998.).OR.
     1 (GOLD(IRES,IGRP).EQ.-1000.)))
      IF(LSPM) THEN
        NORI=NOR(IRES,IGRP)
        KPLI0=LCMGIL(JPLI0,IGRP)
        CALL LCMLEN(KPLI0,'NWT0-PT',ILON,ITYLCM)
        IF(ILON.GT.NBNRS*MAXNOR*NL) THEN
           WRITE(TEXT12,'(3A4)') (IPPT2(IRES,J0),J0=2,4)
           WRITE(HSMG,'(34HUSSIST: FLUX OVERFLOW FOR ISOTOPE ,A12)')
     1     TEXT12
           CALL XABORT(HSMG)
        ENDIF
        CALL LCMGET(KPLI0,'NWT0-PT',XFLUX(1,1,1,IRES))
        DO 295 INOR=1,NORI
        DO 290 IND=1,NBNRS
        XFLUX(IND,INOR,IL,IRES)=XFLUX(IND,INOR,1,IRES)
  290   CONTINUE
  295   CONTINUE
        CALL LCMPUT(KPLI0,'NWT0-PT',NBNRS*MAXNOR*NL,2,XFLUX(1,1,1,IRES))
      ENDIF
  300 CONTINUE
  310 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WCOR,FLNEW,SIGWS,TOTPT,WEIGH,SIGS0X,SIGTXS,XFLUX)
      DEALLOCATE(NPSYS)
      RETURN
      END
