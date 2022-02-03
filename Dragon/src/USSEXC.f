*DECK USSEXC
      SUBROUTINE USSEXC(MAXNOR,CDOOR,IPLI0,IPTRK,IFTRAK,IMPX,NGRP,IGRP,
     1 IASM,NBMIX,NREG,NUN,NL,IPHASE,MAT,VOL,KEYFLX,IREX,SIGGAR,TITR,
     2 NIRES,IRES,NBNRS,NOR,CONR,IPPT1,IPPT2,STGAR,SSGAR,VOLMER,XFLUX,
     3 UNGAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of the flux for the subgroup projection method (SPM) using
* the response matrix method. This is a non-iterative approach which is
* useful in exceptional cases where the fixed-point approach fails.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXNOR  maximum order of the probability tables (PT).
* CDOOR   name of the geometry/solution operator.
* IPLI0   pointer to the internal microscopic cross section library
*         builded by the self-shielding module.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  file unit number used to store the tracks.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* IGRP    index of energy group being processed.
* IASM    offset of information computed by DOORAV or DOORPV.
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
* IREX    fuel region index assigned to each mixture. Equal to zero
*         in non-resonant mixtures or in mixtures not used.
* SIGGAR  macroscopic x-s of the non-resonant isotopes in each mixture:
*         (*,*,*,1) total; (*,*,*,2) transport correction; 
*         (*,*,*,3) P0 scattering; (*,*,*,4) flux times P0 scattering.
* TITR    title.
* NIRES   exact number of correlated resonant isotopes.
* IRES    index of the resonant isotopes being processed.
* NBNRS   number of correlated fuel regions.
* NOR     exact order of the probability table.
* CONR    number density of the resonant isotopes.
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope.
* STGAR   averaged microscopic total xs in resonant region.
* SSGAR   averaged microscopic scattering xs in resonant region.
* VOLMER  volumes of the resonant and non-resonant regions.
*
*Parameters: input/output
* XFLUX   subgroup flux.
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
      INTEGER MAXNOR,IFTRAK,IMPX,NGRP,IGRP,IASM,NBMIX,NREG,NUN,NL,
     1 IPHASE,MAT(NREG),KEYFLX(NREG),IREX(NBMIX),NIRES,IRES,NBNRS,
     2 NOR(NIRES,NGRP),IPPT2(NIRES,4)
      REAL VOL(NREG),SIGGAR(NBMIX,0:NIRES,NGRP,4),CONR(NBNRS,NIRES),
     1 STGAR(NBNRS,NIRES,NGRP),SSGAR(NBNRS,NIRES,NL,NGRP),
     2 VOLMER(0:NBNRS),XFLUX(NBNRS,MAXNOR,NL,NIRES),
     3 UNGAR(NREG,NIRES,NGRP)
      CHARACTER CDOOR*12,TITR*72
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPSYS,KPSYS,JPLI0,KPLIB,IPMACR
      LOGICAL EMPTY,LCM,LEXAC,REBFLG
      CHARACTER CBDPNM*12,TEXT12*12,TEXX12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS
      REAL, ALLOCATABLE, DIMENSION(:) :: AWPHI,FUN,SUN
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WEIGH,TOTPT,SIGWS,PAV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: MATRIX
      TYPE(C_PTR) SIGP_PTR
      REAL, POINTER, DIMENSION(:) :: SIGP
*----
*  STATEMENT FUNCTIONS
*----
      INM(IND,INOR,NBNRS)=(INOR-1)*NBNRS+IND
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WEIGH(MAXNOR,NIRES),TOTPT(MAXNOR,NIRES),
     1 SIGWS(MAXNOR,NIRES),PAV(0:NBNRS,0:NBNRS),AWPHI(0:NBNRS))
      ALLOCATE(MATRIX(NBNRS*MAXNOR,NBNRS*MAXNOR+1))
*----
*  RECOVER THE SPECIFIC DIRECTORY FOR IRES-TH RESONANT ISOTOPE.
*----
      WRITE(CBDPNM,'(3HCOR,I4.4,1H/,I4.4)') IRES,NIRES
      CALL LCMSIX(IPLI0,CBDPNM,1)
      JPLI0=LCMGID(IPLI0,'GROUP_INFO')
      IPSYS=LCMGID(IPLI0,'ASSEMB_RIBON')
      CALL LCMSIX(IPLI0,' ',2)
*----
*  COMPUTE THE AVERAGED COLLISION PROBABILITY MATRIX.
*----
      NORI=NOR(IRES,IGRP)
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
      IOF=(INOR-1)*NUN*(NBNRS+1)+JNBN*NUN+KEYFLX(I)-1
      IND=IREX(IBM)
      IF((JNBN.EQ.0).AND.(IND.EQ.0)) THEN
         T1=T1+SIGGAR(IBM,0,IGRP,3)*VOL(I)
         SUN(IOF+1)=SUN(IOF+1)+SIGGAR(IBM,0,IGRP,3)
      ELSE IF(IND.EQ.JNBN) THEN
         T1=T1+VOL(I)
         SUN(IOF+1)=1.0
      ENDIF
  130 CONTINUE
      DO 140 I=1,NUN
      IOF=(INOR-1)*NUN*(NBNRS+1)+JNBN*NUN+I-1
      IF(T1.NE.0.0) SUN(IOF+1)=SUN(IOF+1)/T1
  140 CONTINUE
  141 CONTINUE
  142 CONTINUE
*----
*  SOLVE FOR THE MULTIBAND FLUX.
*----
      IDIR=0
      NABS=NORI*(NBNRS+1)
      LEXAC=.FALSE.
      IPMACR=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NABS,NBMIX,
     1 IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN,IPMACR,
     2 REBFLG)
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
      IOF=(INOR-1)*NUN*(NBNRS+1)+JNBN*NUN+KEYFLX(I)-1
      PAV(IREX(IBM),JNBN)=PAV(IREX(IBM),JNBN)+FUN(IOF+1)*VOL(I)
  150 CONTINUE
  155 CONTINUE
      DO 165 I=0,NBNRS
      DO 160 J=0,NBNRS
      IF(VOLMER(I).NE.0.0) PAV(I,J)=PAV(I,J)*VOLMER(J)/VOLMER(I)
  160 CONTINUE
  165 CONTINUE
      KPSYS=LCMGIL(IPSYS,IASM+INOR)
      CALL LCMPUT(KPSYS,'DRAGON-PAV',(NBNRS+1)**2,2,PAV(0,0))
  170 CONTINUE
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
         DO 190 INOR=1,NOR(JRES,IGRP)
         WEIGH(INOR,JRES)=SIGP(INOR)
         TOTPT(INOR,JRES)=SIGP(MAXNOR+INOR)
         SIGWS(INOR,JRES)=SIGP(3*MAXNOR+INOR)
  190    CONTINUE
         IF(.NOT.LCM) DEALLOCATE(SIGP)
      ELSE
         WEIGH(1,JRES)=1.0
         TOTPT(1,JRES)=STGAR(IPPT2(JRES,1),JRES,IGRP)
         SIGWS(1,JRES)=SSGAR(IPPT2(JRES,1),JRES,1,IGRP)
      ENDIF
  200 CONTINUE
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
  230    CONTINUE
      ELSE IF((JND.EQ.JNBN).AND.(JNOR.NE.INOR)) THEN
         QQQ=QQQ-WEIGH(JNOR,IRES)*CONR(JND,IRES)*SIGWS(JNOR,IRES)
      ENDIF
      DO 240 IND=1,NBNRS
      AWPHI(IND)=AWPHI(IND)+PAV(IND,JND)*QQQ*VOL(I)/VOLMER(JND)
  240 CONTINUE
  250 CONTINUE
      DO 260 IND=1,NBNRS
      MATRIX(INM(IND,INOR,NBNRS),INM(JNBN,JNOR,NBNRS))=AWPHI(IND)
  260 CONTINUE
  270 CONTINUE
  271 CONTINUE
  272 CONTINUE
*
      DO 280 I=1,NBNRS*NORI
      MATRIX(I,I)=MATRIX(I,I)+1.0D0
  280 CONTINUE
      CALL ALSBD(NBNRS*NORI,1,MATRIX,IER,NBNRS*MAXNOR)
      IF(IER.NE.0) CALL XABORT('USSEXC: SINGULAR MATRIX.')
      CALL XDRSET(XFLUX(1,1,1,IRES),NBNRS*MAXNOR*NL,0.0)
      DO 295 IND=1,NBNRS
      DO 290 INOR=1,NORI
      I1=INM(IND,INOR,NBNRS)
      XFLUX(IND,INOR,1,IRES)=REAL(MATRIX(I1,NBNRS*NORI+1))
  290 CONTINUE
  295 CONTINUE
* END OF RESPONSE MATRIX APPROACH.
*
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
         DO 330 INOR=1,NOR(JRES,IGRP)
         WEIGH(INOR,JRES)=SIGP(INOR)
         SIGWS(INOR,JRES)=SIGP(3*MAXNOR+INOR)
  330    CONTINUE
         IF(.NOT.LCM) DEALLOCATE(SIGP)
      ELSE
         WEIGH(1,JRES)=1.0
         SIGWS(1,JRES)=SSGAR(IPPT2(JRES,1),JRES,1,IGRP)
      ENDIF
  340 CONTINUE
*----
*  COMPUTE THE AVERAGED SOURCE.
*----
      ALLOCATE(FUN(NUN*NORI),SUN(NUN*NORI))
      CALL XDRSET(SUN,NUN*NORI,0.0)
      ALLOCATE(NPSYS(NORI))
      DO 385 INOR=1,NORI
      NPSYS(INOR)=IASM+INOR
      KPSYS=LCMGIL(IPSYS,IASM+INOR)
      CALL LCMLEN(KPSYS,'FUNKNO$USS',ILENG,ITYLCM)
      IF(ILENG.EQ.NUN) THEN
        CALL LCMGET(KPSYS,'FUNKNO$USS',FUN((INOR-1)*NUN+1))
      ELSE
        CALL XDRSET(FUN((INOR-1)*NUN+1),NUN,0.0)
      ENDIF
      IOF=(INOR-1)*NUN
      DO 380 I=1,NREG
      IBM=MAT(I)
      QQQ=0.0
      IF(IBM.EQ.0) GO TO 380
      IND=IREX(IBM)
      DO 360 JRES=0,NIRES
      IF(JRES.EQ.0) THEN
         QQQ=QQQ+SIGGAR(IBM,0,IGRP,3)
      ELSE IF((JRES.NE.IRES).AND.(IND.GT.0)) THEN
         QQQ=QQQ+SIGGAR(IBM,JRES,IGRP,4)
      ENDIF
  360 CONTINUE
      IF(IND.GT.0) THEN
         DO 370 JNOR=1,NORI
         IF(JNOR.NE.INOR) THEN
            QQQ=QQQ+WEIGH(JNOR,IRES)*CONR(IND,IRES)*SIGWS(JNOR,IRES)*
     1      XFLUX(IND,JNOR,1,IRES)
         ENDIF
  370    CONTINUE
      ENDIF
      SUN(IOF+KEYFLX(I))=QQQ*WEIGH(INOR,IRES)
  380 CONTINUE
  385 CONTINUE
*
      IF(IMPX.GT.0) THEN
         WRITE(TEXT12,'(3A4)') (IPPT2(IRES,I),I=2,4)
         WRITE(6,'(15H USSEXC: GROUP=,I5,24H. SUBGROUP CALCULATION B,
     1   37HASED ON RESPONSE MATRICES.  ISOTOPE='',A12,2H''.)') IGRP,
     2   TEXT12
      ENDIF
*----
*  SOLVE FOR THE MULTIBAND FLUX (VECTOR OF LENGTH NREG).
*----
      IPMACR=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NORI,NBMIX,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN,IPMACR,REBFLG)
      DEALLOCATE(NPSYS)
*----
*  INTEGRATE THE REGION-ORDERED FLUX OVER SUBGROUPS AND COMPUTE UNGAR,
*  THE REGION-ORDERED FLUX.
*----
      CALL XDRSET(UNGAR(1,IRES,IGRP),NREG,0.0)
      DO 420 INOR=1,NORI
      KPSYS=LCMGIL(IPSYS,IASM+INOR)
      IOF=(INOR-1)*NUN
      CALL LCMPUT(KPSYS,'FUNKNO$USS',NUN,2,FUN(IOF+1))
*
      DO 410 I=1,NREG
      IOF=(INOR-1)*NUN+KEYFLX(I)-1
      UNGAR(I,IRES,IGRP)=UNGAR(I,IRES,IGRP)+FUN(IOF+1)
  410 CONTINUE
  420 CONTINUE
      DEALLOCATE(SUN,FUN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MATRIX)
      DEALLOCATE(AWPHI,PAV,SIGWS,TOTPT,WEIGH)
      RETURN
      END
