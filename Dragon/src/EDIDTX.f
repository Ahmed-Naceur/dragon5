*DECK EDIDTX
      SUBROUTINE EDIDTX(IPEDIT,IPFLUX,IPMACR,IADJ,IPRINT,NL,NDEL,NALBP,
     >                  ITRANC,NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,ILEAKS,
     >                  ILUPS,NW,MATCOD,VOLUME,KEYFLX,IGCOND,IMERGE,
     >                  FLUXES,AFLUXE,EIGENK,VOLMER,WLETYC,WENERG,
     >                  RATECM,FLUXCM,FADJCM,FLXINT,SCATTD,SCATTS,
     >                  NIFISS,NSAVES,CURNAM,NEDMAC,SIGS,B2,CUREIN,
     >                  TIMEF,NTAUXT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate and print macroscopic reaction rates.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPEDIT  pointer to the edition LCM object.
* IPFLUX  pointer to the solution LCM object.
* IPMACR  pointer to the macrolib LCM object.
* IADJ    type of flux weighting:
*         = 0 direct flux weighting;
*         = 1 direct-adjoint flux weighting.
* IPRINT  print level;
*         = 0 no print;
*         = 1 print fluxes;
*         = 2 1+print reaction rates;
*         = 3 2+print homogenized cross sections.
* NL      number of Legendre orders.
* NDEL    number of delayed precursor groups.
* NALBP   number of physical albedos.
* ITRANC  type of transport corrections.
* NGROUP  number of groups.
* NGCOND  number of groups condensed.
* NBMIX   number of mixtures.
* NREGIO  number of regions.
* NMERGE  number of merged regions.
* ILEAKS  type of leakage calculation:
*         = 0 no leakage;
*         = 1 homogeneous leakage (Diffon);
*         = 2 isotropic streaming (Ecco);
*         = 3 anisotropic streaming (Tibere).
* ILUPS   flag to remove up-scattering from output.
* NW      type of weighting for P1 cross section info (=0 P0; =1 P1).
* MATCOD  material per region.
* VOLUME  volume of region.
* KEYFLX  average flux position per region.
* IGCOND  limit condensed groups.
* IMERGE  index of merged regions.
* FLUXES  fluxes.
* AFLUXE  adjoint fluxes.
* EIGENK  eigenvalue for problem.
* B2      square buckling:
*         for ILEAKS=1,2: B2(4) is homogeneous;
*         for ILEAKS=3: B2(1),B2(2),B2(3) are directional heterogeneous
*         and B2(4) is homogeneous.
* CUREIN  infinite multiplication factor.
* NTAUXT  number of reaction rate edits (=15+2*NDEL).
* TIMEF   time stamp in day/burnup/irradiation.
*
*Parameters: output
* VOLMER  volume of region merged.
* WLETYC  lethargy width condensed.
* WENERG  energy group limits.
* RATECM  averaged region/group cross sections:
*         = RATECM(*,1) = total P0;
*         = RATECM(*,2) = total P1;
*         = RATECM(*,NW+2) = absorption;
*         = RATECM(*,NW+3) = fission;
*         = RATECM(*,NW+4) = fixed sources / productions;
*         = RATECM(*,NW+5) = leakage;
*         = RATECM(*,NW+6) = total out of group scattering;
*         = RATECM(*,NW+7) = diagonal scattering x-s;
*         = RATECM(*,NW+8) = chi;
*         = RATECM(*,NW+9) = wims type transport correction;
*         = RATECM(*,NW+10) = x-directed leakage;
*         = RATECM(*,NW+11) = y-directed leakage;
*         = RATECM(*,NW+12) = z-directed leakage;
*         = RATECM(*,NW+13) = nu-sigf for delayed neutrons;
*         = RATECM(*,NW+13+NDEL) = fission spectra for delayed neutrons.
* FLUXCM  integrated region/group fluxes:
*         = FLUXCM(*,1) = fluxes P0;
*         = FLUXCM(*,2) = fluxes P1.
* FADJCM  averaged region/group adjoint fluxes:
*         = FADJCM(*,1) = adjoint fluxes P0;
*         = FADJCM(*,2) = adjoint fluxes P1.
* FLXINT  integrated flux.
* SCATTD  scattering rates.
* SCATTS  homogenized scattering cross sections.
* NIFISS  number of fissile isotopes.
* NSAVES  homogenized x-s compute/save flag:
*         = 0  no compute, no save;
*         = 1  compute, no save;
*         = 2  compute and save.
* CURNAM  name of LCM directory where the merged/condensed x-s are
*         stored.
* NEDMAC  number of extra edit vectors.
* SIGS    Legendre dependent scattering cross sections.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT,IPFLUX,IPMACR
      INTEGER     IADJ,IPRINT,NL,NDEL,NALBP,ITRANC,NGROUP,NGCOND,NBMIX,
     >            NREGIO,NMERGE,ILEAKS,ILUPS,NW,MATCOD(NREGIO),
     >            KEYFLX(NREGIO),IGCOND(NGCOND),IMERGE(NREGIO),
     >            NIFISS,NSAVES,NEDMAC,NTAUXT
      REAL        VOLUME(NREGIO),FLUXES(NREGIO,NGROUP,NW+1),
     >            AFLUXE(NREGIO,NGROUP,NW+1),EIGENK,VOLMER(NMERGE),
     >            WENERG(NGCOND+1),RATECM(NMERGE,NGCOND,NTAUXT),
     >            FLUXCM(NMERGE,NGCOND,NW+1),FADJCM(NMERGE,NGCOND,NW+1),
     >            FLXINT(NREGIO,NGROUP,NW+1),WLETYC(NGCOND),
     >            SCATTS(NMERGE,NGCOND,NGCOND,NL),
     >            SIGS(NMERGE,NGCOND,NL),B2(4),CUREIN,TIMEF(3)
      CHARACTER   CURNAM*12
      DOUBLE PRECISION SCATTD(NMERGE,NGCOND,NGCOND,NL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLUX,JPMACR,KPMACR
      CHARACTER   APG*3
      PARAMETER  (IUNOUT=6,APG=' > ',ILCMUP=1,ILCMDN=2)
      CHARACTER   TEXT12*12,CM*2,OPTION*4
      LOGICAL     LH
      DOUBLE PRECISION SCATW,TOTFIS,FXSOUR,FLFUEL,FCELL
      INTEGER     IFSKP,ISKP(3)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IFUELR,INGSCT,IFGSCT,IPOSCT
      REAL, ALLOCATABLE, DIMENSION(:) :: DISFCT,SIGMA,XSCAT,WORKF,
     > ENERG,DIFFU,SIGMAF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FFUEL,FLDMC,OVERV,HFACT,
     > DECAY,ALBP,ALBPGR
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TAUXE
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HVECT
*----
*  SCRATCH STORAGE ALLOCATION
*   HVECT   extra edit names.
*   IFUELR  fuel region location.
*   DISFCT  disadvantage factor.
*   TAUXE   extra edit rates.
*   FFUEL   flux in fuel.
*   FLDMC   fission rate condensed.
*   OVERV   1/v merge condensed.
*   HFACT   H-factors.
*   DECAY   precursor decay constants.
*   ALBP    physical albedos.
*----
      ALLOCATE(HVECT(NEDMAC),IFUELR(NREGIO))
      ALLOCATE(DISFCT(NGCOND),TAUXE(NMERGE,NGCOND,NEDMAC),
     > FFUEL(NREGIO,NIFISS),FLDMC(NMERGE,NGCOND),OVERV(NMERGE,NGCOND),
     > HFACT(NMERGE,NGCOND),DECAY(NDEL,NIFISS),ALBP(NALBP,NGCOND))
*----
*  ALLOCATE WORK VECTOR AND INITIALIZE REQUIRED VECTORS
*----
      ILEAK2=ILEAKS
      ALLOCATE(INGSCT(NBMIX),IFGSCT(NBMIX),IPOSCT(NBMIX))
      ALLOCATE(SIGMA(0:NBMIX*MAX(NIFISS,1)),XSCAT(NBMIX*NGROUP))
      CALL XDRSET(RATECM,NMERGE*NGCOND*NTAUXT,0.0)
      CALL XDRSET(FLUXCM,NMERGE*NGCOND*(NW+1),0.0)
      CALL XDRSET(FADJCM,NMERGE*NGCOND*(NW+1),0.0)
      CALL XDRSET(SIGS,NMERGE*NGCOND*NL,0.0)
      CALL XDRSET(OVERV,NMERGE*NGCOND,0.0)
      CALL XDRSET(HFACT,NMERGE*NGCOND,0.0)
      CALL XDRSET(TAUXE,NMERGE*NGCOND*NEDMAC,0.0)
      CALL XDRSET(VOLMER,NMERGE,0.0)
      CALL XDRSET(FFUEL,NREGIO*NIFISS,0.0)
      CALL XDISET(IFUELR,NREGIO,0)
      SIGMA(0)=0.0
      IF(IADJ.EQ.0) THEN
        IOP=1
      ELSE IF(IADJ.EQ.1) THEN
        IOP=11
      ENDIF
*----
*  FIND EDIT XS
*----
      IF(NEDMAC.GT.0) CALL LCMGTC(IPMACR,'ADDXSNAME-P0',8,NEDMAC,HVECT)
*----
*  ENERGY AND LETHARGY CONDENSATION
*----
      CALL LCMLEN(IPMACR,'ENERGY',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.0) THEN
        NENER=0
        CALL XDRSET(WLETYC,NGCOND,0.0)
        CALL XDRSET(WENERG,NGCOND+1,0.0)
      ELSE IF(ILCMLN.EQ.NGROUP+1) THEN
        NENER=NGCOND+1
        ALLOCATE(ENERG(NGROUP+1))
        CALL LCMGET(IPMACR,'ENERGY',ENERG)
        WENERG(1)=ENERG(1)
        DO 30 IGC=1,NGCOND
          WENERG(IGC+1)=ENERG(IGCOND(IGC)+1)
          WLETYC(IGC)=LOG(WENERG(IGC)/WENERG(IGC+1))
  30    CONTINUE
        IF(ENERG(NGROUP+1).EQ.0.0) ENERG(NGROUP+1)=1.0E-5
        DEALLOCATE(ENERG)
      ELSE
        CALL XABORT('EDIDTX: READ ERROR INVALID NUMBER OF GROUPS')
      ENDIF
*----
*  COMPUTE MERGED VOLUME
*----
      DO 50 IREGIO=1,NREGIO
        IKK=IMERGE(IREGIO)
        IF(IKK.GT.0) THEN
          VOLMER(IKK)=VOLMER(IKK)+VOLUME(IREGIO)
        ENDIF
  50  CONTINUE
*----
*  COMPUTE INTEGRATED/CONDENSED FUNDAMENTAL CURRENTS (ILEAKS=2,3)
*----
      IF(ILEAKS.EQ.2) THEN
        IF(IADJ.EQ.1) CALL XABORT('EDIMIC: DIRECT-ADJOINT WEIGTING NOT'
     >  //' IMPLEMENTED.')
        CALL LCMLEN(IPFLUX,'FLUX',ILCMLN,ITYLCM)
        IF(ILCMLN.EQ.0) CALL XABORT('EDIDTX: MISSING FLUX INFO.')
        JPFLUX=LCMGID(IPFLUX,'FLUX')
        CALL LCMLEL(JPFLUX,1,ILCMLN,ITYLCM)
        ALLOCATE(WORKF(ILCMLN))
        DO 70 IGR=1,NGROUP
          CALL LCMGDL(JPFLUX,IGR,WORKF)
          DO 60 IREG=1,NREGIO
            FLXINT(IREG,IGR,1)=WORKF(KEYFLX(IREG)+ILCMLN/2)*VOLUME(IREG)
  60      CONTINUE
  70    CONTINUE
        IGRFIN=0
        DO 90 IGRC=1,NGCOND
          IGRDEB=IGRFIN+1
          IGRFIN=IGCOND(IGRC)
          DO 80 IGR=IGRDEB,IGRFIN
*----
*  COMPUTE MERGED INTEGRATED CURRENTS
*----
            CALL EDIRAT(0,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >      VOLUME,RATECM(1,IGRC,NW+5),SIGMA(0),IMERGE,NMERGE)
  80      CONTINUE
  90    CONTINUE
        DEALLOCATE(WORKF)
      ELSE IF(ILEAKS.EQ.3) THEN
        IF(IADJ.EQ.1) CALL XABORT('EDIMIC: DIRECT-ADJOINT WEIGTING NOT'
     >  //' IMPLEMENTED.')
        CALL LCMLEN(IPFLUX,'FLUX',ILCMLN,ITYLCM)
        IF(ILCMLN.EQ.0) CALL XABORT('EDIDTX: MISSING FLUX INFO.')
        JPFLUX=LCMGID(IPFLUX,'FLUX')
        CALL LCMLEL(JPFLUX,1,ILCMLN,ITYLCM)
*----
*  CALCULATIONS FOR TIBERE PIJ
*----
        IF(ILCMLN.EQ.12*NREGIO) THEN
          IFSKP=3*NREGIO
          ALLOCATE(WORKF(ILCMLN))
          DO 140 IDIR=1,3
            DO 110 IGR=1,NGROUP
              CALL LCMGDL(JPFLUX,IGR,WORKF)
              DO 100 IREG=1,NREGIO
                FLXINT(IREG,IGR,1)=WORKF(KEYFLX(IREG)+IDIR*IFSKP)
     >          *VOLUME(IREG)
 100          CONTINUE
 110        CONTINUE
            IGRFIN=0
            DO 130 IGRC=1,NGCOND
              IGRDEB=IGRFIN+1
              IGRFIN=IGCOND(IGRC)
              DO 120 IGR=IGRDEB,IGRFIN
*----
*  COMPUTE MERGED INTEGRATED CURRENTS FOR TIBERE PIJ
*----
                CALL EDIRAT(0,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >                     VOLUME,RATECM(1,IGRC,NW+9+IDIR),SIGMA(0),
     >                     IMERGE,NMERGE)
 120          CONTINUE
 130        CONTINUE
 140      CONTINUE
          DEALLOCATE(WORKF)
*----
*  CALCULATIONS FOR TIBERE MoC
*----
        ELSE IF(ILCMLN.NE.0) THEN
          ALLOCATE(WORKF(ILCMLN))
          DO 141 IDIR=1,3
            DO 111 IGR=1,NGROUP
              CALL LCMGDL(JPFLUX,IGR,WORKF)
              DO 101 IREG=1,NREGIO
                ISKP(1)=ILCMLN/4+KEYFLX(IREG)
                ISKP(2)=ILCMLN/2+KEYFLX(IREG)
                ISKP(3)=3*ILCMLN/4+KEYFLX(IREG)
                FLXINT(IREG,IGR,1)=WORKF(ISKP(IDIR))
     >          *VOLUME(IREG)
 101          CONTINUE
 111        CONTINUE
            IGRFIN=0
            DO 131 IGRC=1,NGCOND
              IGRDEB=IGRFIN+1
              IGRFIN=IGCOND(IGRC)
              DO 121 IGR=IGRDEB,IGRFIN
*----
*  COMPUTE MERGED INTEGRATED CURRENTS FOR TIBERE MOC
*----
                CALL EDIRAT(0,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >                     VOLUME,RATECM(1,IGRC,NW+9+IDIR),SIGMA(0),
     >                     IMERGE,NMERGE)
 121          CONTINUE
 131        CONTINUE
 141      CONTINUE
          DEALLOCATE(WORKF)
        ENDIF
      ENDIF
*----
*  COMPUTE INTEGRATED FLUX
*----
      DO 170 IW=1,NW+1
        DO 160 IGR=1,NGROUP
          DO 150 IREGIO=1,NREGIO
            FLXINT(IREGIO,IGR,IW)=FLUXES(IREGIO,IGR,IW)*VOLUME(IREGIO)
 150      CONTINUE
 160    CONTINUE
 170  CONTINUE
*----
*  COMPUTE INTEGRATED/CONDENSED FUNDAMENTAL CURRENTS (ILEAKS=1)
*  (OBTAINED AS THE PRODUCT OF THE FUNDAMENTAL FLUX BY THE LEAKAGE
*  COEFFICIENT)
*----
      JPMACR=LCMGID(IPMACR,'GROUP')
      IF(ILEAKS.EQ.1) THEN
        ALLOCATE(DIFFU(NGROUP))
        CALL LCMGET(IPFLUX,'DIFFB1HOM',DIFFU)
        CALL LCMGTC(IPFLUX,'OPTION',4,1,OPTION)
        IGRFIN=0
        DO 200 IGRC=1,NGCOND
          IGRDEB=IGRFIN+1
          IGRFIN=IGCOND(IGRC)
          DO 190 IGR=IGRDEB,IGRFIN
            IF(OPTION.EQ.'LKRD') THEN
              KPMACR=LCMGIL(JPMACR,IGR)
              CALL LCMGET(KPMACR,'DIFF',SIGMA(1))
            ELSE
              DO 180 IMIX=0,NBMIX
                SIGMA(IMIX)=DIFFU(IGR)
 180          CONTINUE
            ENDIF
            CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >      AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+5),SIGMA(0),IMERGE,NMERGE)
 190      CONTINUE
 200    CONTINUE
        DEALLOCATE(DIFFU)
      ENDIF
*----
*  READ FIXE SOURCES/COMPUTE FIXE PRODUCTION RATE AND TOTAL SOURCE
*----
      IGRFIN=0
      TOTFIS=0.0D0
      FXSOUR=0.0D0
      DO 250 IGRC=1,NGCOND
        IGRDEB=IGRFIN+1
        IGRFIN=IGCOND(IGRC)
        DO 240 IGR=IGRDEB,IGRFIN
          KPMACR=LCMGIL(JPMACR,IGR)
          CALL LCMLEN(KPMACR,'FIXE',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMACR,'FIXE',SIGMA(1))
            CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,VOLUME,AFLUXE(1,IGR,1),
     >      RATECM(1,IGRC,NW+4),SIGMA(0),IMERGE,NMERGE)
            DO 210 IKK=1,NMERGE
              FXSOUR=FXSOUR+DBLE(RATECM(IKK,IGRC,NW+4))
 210        CONTINUE
          ENDIF
          DO 215 IW=1,NW+1
            CALL EDIRAT(0,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,IW),VOLUME,
     >                  FLUXCM(1,IGRC,IW),SIGMA(0),IMERGE,NMERGE)
            IF(IADJ.EQ.1) THEN
              CALL EDIRAT(10,NREGIO,NBMIX,MATCOD,AFLUXE(1,IGR,IW),
     >                  VOLUME,FADJCM(1,IGRC,IW),SIGMA(0),IMERGE,NMERGE)
              DO IKK=1,NMERGE
                FADJCM(IKK,IGRC,IW)=FADJCM(IKK,IGRC,IW)/VOLMER(IKK)
              ENDDO
            ENDIF
 215      CONTINUE
*----
*  READ FISSION X-S/ COMPUTE FISSION RATES
*----
          CALL LCMLEN(KPMACR,'NUSIGF',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMACR,'NUSIGF',SIGMA(1))
            DO 230 IFIS=1,NIFISS
              DO 220 IREGIO=1,NREGIO
                IBM=MATCOD(IREGIO)
                IF(IBM.GT.0) THEN
                 IF(SIGMA((IFIS-1)*NBMIX+IBM).GT.0.0) THEN
                   FLXFIS=FLXINT(IREGIO,IGR,1)*SIGMA((IFIS-1)*NBMIX+IBM)
                   FFUEL(IREGIO,IFIS)=FFUEL(IREGIO,IFIS)+FLXFIS
                   TOTFIS=TOTFIS+DBLE(FLXFIS)
                   IFUELR(IREGIO)=1
                 ENDIF
                ENDIF
 220          CONTINUE
 230        CONTINUE
          ENDIF
*----
 240    CONTINUE
 250  CONTINUE
*----
*  RECOVER THE PRECURSOR RADIOACTIVE DECAY CONSTANTS. USE THE VALUES
*  OF THE FISSILE ISOTOPE WITH MAXIMUM FISSION RATE
*----
      IF(CURNAM.NE.' ') THEN
        CALL LCMLEN(IPMACR,'LAMBDA-D',ILCMLN,ITYLCM)
        IF((NDEL.GT.0).AND.(ILCMLN.GT.0)) THEN
          ZMAX=0.0
          KFIS=0
          DO 340 IFIS=1,NIFISS
            ZTOT=0.0
            DO 330 IREGIO=1,NREGIO
            ZTOT=ZTOT+FFUEL(IREGIO,IFIS)
 330        CONTINUE
            IF(ZTOT.GE.ZMAX) THEN
              KFIS=IFIS
              ZMAX=ZTOT
            ENDIF
 340      CONTINUE
          CALL LCMGET(IPMACR,'LAMBDA-D',DECAY)
          CALL LCMSIX(IPEDIT,CURNAM,ILCMUP)
          CALL LCMSIX(IPEDIT,'MACROLIB',ILCMUP)
          CALL LCMPUT(IPEDIT,'LAMBDA-D',NDEL,2,DECAY(1,KFIS))
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
        ENDIF
      ENDIF
*----
*  FIND FUEL VOLUME FOR DISADVANTAGE FACTOR
*----
      VFUEL=0.0
      VCELL=0.0
      DO 350 IREGIO=1,NREGIO
        IF(IFUELR(IREGIO).EQ.1) THEN
          VFUEL=VFUEL+VOLUME(IREGIO)
        ENDIF
        VCELL=VCELL+VOLUME(IREGIO)
 350  CONTINUE
      LH=.FALSE.
      IGRFIN=0
      DO 510 IGRC=1,NGCOND
        FCELL=0.0D0
        FLFUEL=0.0D0
        IGRDEB=IGRFIN+1
        IGRFIN=IGCOND(IGRC)
        DO 380 JGRC=1,NGCOND
          DO 370 I=1,NMERGE
            DO 360 IL=1,NL
              SCATTD(I,IGRC,JGRC,IL)=0.0D0
 360        CONTINUE
 370      CONTINUE
 380    CONTINUE
        DO 500 IGR=IGRDEB,IGRFIN
          KPMACR=LCMGIL(JPMACR,IGR)
*----
*  INTEGRATED 1/V
*----
          CALL LCMLEN(KPMACR,'OVERV',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMACR,'OVERV',SIGMA(1))
            CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >                  AFLUXE(1,IGR,1),OVERV(1,IGRC),SIGMA(0),
     >                  IMERGE,NMERGE)
          ENDIF
*----
*  INTEGRATED H-FACTORS
*----
          CALL LCMLEN(KPMACR,'H-FACTOR',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            LH=.TRUE.
            CALL LCMGET(KPMACR,'H-FACTOR',SIGMA(1))
            CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >                  AFLUXE(1,IGR,1),HFACT(1,IGRC),SIGMA(0),
     >                  IMERGE,NMERGE)
          ENDIF
*----
*  TOTAL, ABSROPTION, ETC. RATES
*----
          CALL LCMLEN(KPMACR,'NTOT0',ILCMLN,ITYLCM)
          IF(ILCMLN.EQ.0) CALL XABORT('EDIDTX: READ ERROR ON LCM REC'//
     >    'ORD= TOTAL')
          CALL LCMGET(KPMACR,'NTOT0',SIGMA(1))
          CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >    AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+2),SIGMA(0),IMERGE,NMERGE)
          DO 385 IW=1,NW+1
            WRITE(TEXT12,'(4HNTOT,I1)') IW-1
            CALL LCMLEN(KPMACR,TEXT12,ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) CALL LCMGET(KPMACR,TEXT12,SIGMA(1))
            CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,IW),
     >      AFLUXE(1,IGR,IW),RATECM(1,IGRC,IW),SIGMA(0),IMERGE,
     >      NMERGE)
 385      CONTINUE
          CALL LCMLEN(KPMACR,'SIGS00',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMACR,'SIGS00',SIGMA(1))
            CALL EDIRAT(-IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >      AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+2),SIGMA(0),IMERGE,NMERGE)
          ENDIF
          IF(ILEAKS.EQ.0) THEN
            CALL LCMLEN(KPMACR,'DIFF',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMACR,'DIFF',SIGMA(1))
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+5),SIGMA(0),IMERGE,
     >        NMERGE)
              ILEAK2=10
            ENDIF
            CALL LCMLEN(KPMACR,'DIFFX',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMACR,'DIFFX',SIGMA(1))
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+10),SIGMA(0),IMERGE,
     >        NMERGE)
              ILEAK2=11
            ENDIF
            CALL LCMLEN(KPMACR,'DIFFY',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMACR,'DIFFY',SIGMA(1))
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+11),SIGMA(0),IMERGE,
     >        NMERGE)
              ILEAK2=11
            ENDIF
            CALL LCMLEN(KPMACR,'DIFFZ',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMACR,'DIFFZ',SIGMA(1))
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+12),SIGMA(0),IMERGE,
     >        NMERGE)
              ILEAK2=11
            ENDIF
          ENDIF
*----
*  READ ADDITIONAL X-SECTIONS
*----
          DO 390 IED=1,NEDMAC
            IF(HVECT(IED)(:2).EQ.'NW') GO TO 390
            CALL LCMLEN(KPMACR,HVECT(IED),ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
              CALL LCMGET(KPMACR,HVECT(IED),SIGMA(1))
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >                    AFLUXE(1,IGR,1),TAUXE(1,IGRC,IED),SIGMA(0),
     >                    IMERGE,NMERGE)
            ENDIF
 390      CONTINUE
*----
*  FISSION SPECTRUM AND NU*SIGF
*----
          CALL LCMLEN(KPMACR,'NUSIGF',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            ALLOCATE(SIGMAF(0:NBMIX))
            SIGMAF(0)=0.0
            CALL LCMGET(KPMACR,'NUSIGF',SIGMA(1))
            DO 400 IFIS=1,NIFISS
              DO 395 IBM=1,NBMIX
                SIGMAF(IBM)=SIGMA((IFIS-1)*NBMIX+IBM)
 395          CONTINUE
              CALL EDIRAT(1,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+3),SIGMAF(0),IMERGE,
     >        NMERGE)
 400        CONTINUE
            DEALLOCATE(SIGMAF)
          ENDIF
          CALL LCMLEN(KPMACR,'CHI',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            ALLOCATE(SIGMAF(0:NBMIX))
            SIGMAF(0)=0.0
            CALL LCMGET(KPMACR,'CHI',SIGMA(1))
            DO 410 IFIS=1,NIFISS
              DO 405 IBM=1,NBMIX
                SIGMAF(IBM)=SIGMA((IFIS-1)*NBMIX+IBM)
 405          CONTINUE
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FFUEL(1,IFIS),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+4),SIGMAF(0),IMERGE,
     >        NMERGE)
 410        CONTINUE
            DEALLOCATE(SIGMAF)
          ENDIF
*----
*  DELAYED FISSION SPECTRUM AND NU*SIGF
*----
          DO 440 IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            CALL LCMLEN(KPMACR,TEXT12,ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              ALLOCATE(SIGMAF(0:NBMIX))
              SIGMAF(0)=0.0
              CALL LCMGET(KPMACR,TEXT12,SIGMA(1))
              DO 420 IFIS=1,NIFISS
                DO 415 IBM=1,NBMIX
                  SIGMAF(IBM)=SIGMA((IFIS-1)*NBMIX+IBM)
 415            CONTINUE
                CALL EDIRAT(1,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >          VOLUME,RATECM(1,IGRC,12+NW+IDEL),SIGMAF(0),IMERGE,
     >          NMERGE)
 420          CONTINUE
              DEALLOCATE(SIGMAF)
            ENDIF
            WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
            CALL LCMLEN(KPMACR,TEXT12,ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              ALLOCATE(SIGMAF(0:NBMIX))
              SIGMAF(0)=0.0
              CALL LCMGET(KPMACR,TEXT12,SIGMA(1))
              DO 430 IFIS=1,NIFISS
                DO 425 IBM=1,NBMIX
                  SIGMAF(IBM)=SIGMA((IFIS-1)*NBMIX+IBM)
 425            CONTINUE
                CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FFUEL(1,IFIS),
     >          AFLUXE(1,IGR,1),RATECM(1,IGRC,12+NW+NDEL+IDEL),
     >          SIGMAF(0),IMERGE,NMERGE)
 430          CONTINUE
              DEALLOCATE(SIGMAF)
            ENDIF
 440      CONTINUE
*----
*  INTEGRATED FLUX AND FORM FACTOR
*----
          DO 450 IREGIO=1,NREGIO
            IF(IFUELR(IREGIO).EQ.1) THEN
              FLFUEL=FLFUEL+DBLE(FLXINT(IREGIO,IGR,1))
            ENDIF
            FCELL=FCELL+DBLE(FLXINT(IREGIO,IGR,1))
 450      CONTINUE
*----
*  TRANSPORT CORRECTION
*----
          IF(ITRANC.NE.0) THEN
            CALL LCMLEN(KPMACR,'TRANC',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMACR,'TRANC',SIGMA(1))
              CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,1),
     >        AFLUXE(1,IGR,1),RATECM(1,IGRC,NW+9),SIGMA(0),IMERGE,
     >        NMERGE)
            ENDIF
          ENDIF
*----
*  SCATTERING NEUTRONS
*----
          DO 490 IL=1,NL
          IW=MIN(IL,NW+1)
          WRITE (CM,'(I2.2)') IL-1
          CALL LCMLEN(KPMACR,'SIGS'//CM,ILCSCA,ITYLCM)
          IF(ILCSCA.GT.0) THEN
            CALL LCMGET(KPMACR,'SIGS'//CM,SIGMA(1))
            CALL EDIRAT(IOP,NREGIO,NBMIX,MATCOD,FLXINT(1,IGR,IW),
     >      AFLUXE(1,IGR,IW),SIGS(1,IGRC,IL),SIGMA(0),IMERGE,NMERGE)
          ENDIF
          CALL LCMLEN(KPMACR,'NJJS'//CM,ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMACR,'SIGW'//CM,SIGMA(1))
            CALL LCMGET(KPMACR,'NJJS'//CM,INGSCT)
            CALL LCMGET(KPMACR,'IJJS'//CM,IFGSCT)
            CALL LCMGET(KPMACR,'IPOS'//CM,IPOSCT)
            CALL LCMGET(KPMACR,'SCAT'//CM,XSCAT)
            DO 480 IREGIO=1,NREGIO
              MATNUM=MATCOD(IREGIO)
              IKK=IMERGE(IREGIO)
              IF((IKK.GT.0).AND.(MATNUM.GT.0)) THEN
                NGSCAT=INGSCT(MATNUM)
                IGSCAT=IFGSCT(MATNUM)
                IPOSIT=IPOSCT(MATNUM)
                JGRFIN=0
                FAD=1.0
                IF(IADJ.EQ.1) FAD=AFLUXE(IREGIO,IGR,IW)
                DO 470 JGRC=1,NGCOND
                  JGRDEB=JGRFIN+1
                  JGRFIN=IGCOND(JGRC)
                  J2=MIN(JGRFIN,IGSCAT)
                  J1=MAX(JGRDEB,IGSCAT-NGSCAT+1)
                  IPO=IPOSIT+IGSCAT-J2
                  DO 460 JGR=J2,J1,-1
                    IF(IGR.EQ.JGR) THEN
                      SCATTD(IKK,IGRC,JGRC,IL)=SCATTD(IKK,IGRC,JGRC,IL)
     >                   +SIGMA(MATNUM)*FLXINT(IREGIO,JGR,IW)*FAD
                    ELSE
                      SCATTD(IKK,IGRC,JGRC,IL)=SCATTD(IKK,IGRC,JGRC,IL)
     >                   +XSCAT(IPO)*FLXINT(IREGIO,JGR,IW)*FAD
                    ENDIF
                    IPO=IPO+1
 460              CONTINUE
 470            CONTINUE
              ENDIF
 480        CONTINUE
          ENDIF
 490      CONTINUE
 500    CONTINUE
        IF(VFUEL*FCELL.GT.0.0) THEN
          DISFCT(IGRC)=REAL(FLFUEL*VCELL/(VFUEL*FCELL))
        ELSE
          DISFCT(IGRC)=0.0
        ENDIF
 510  CONTINUE
*----
*  UP-SCATTERING CORRECTIONS
*----
      IF(ILUPS.EQ.1) THEN
        DO 523 IKK=1,NMERGE
        DO 522 IGRC=2,NGCOND
        DO 521 JGRC=1,IGRC-1
        DO 520 IL=1,NL
          SCATTD(IKK,IGRC,JGRC,IL)=
     >      SCATTD(IKK,IGRC,JGRC,IL)-SCATTD(IKK,JGRC,IGRC,IL)
          SCATTD(IKK,JGRC,IGRC,IL)=0.0D0
 520    CONTINUE
 521    CONTINUE
 522    CONTINUE
 523    CONTINUE
      ENDIF
*----
*  SCATTERING NORMALIZATION
*----
      IF(IADJ.EQ.0) THEN
        DO 560 IGRC=1,NGCOND
          DO 550 IKK=1,NMERGE
            DO 540 IL=1,NL
              IF(ILCSCA.GT.0) THEN
                SCATW=SIGS(IKK,IGRC,IL)
                DO 530 JGRC=1,NGCOND
                  IF(JGRC.NE.IGRC) SCATW=SCATW-SCATTD(IKK,JGRC,IGRC,IL)
 530            CONTINUE
                DEN=REAL(MAX(ABS(SCATW),ABS(SCATTD(IKK,IGRC,IGRC,IL))))
                IF(DEN.GT.0.0) THEN
                  ERR=ABS(REAL(SCATW-SCATTD(IKK,IGRC,IGRC,IL)))/DEN
                  IF(ERR.GT.1.0E-3) THEN
                    WRITE(IUNOUT,6000) IL,IGRC,IKK,100.0*ERR
                  ENDIF
                  SCATTD(IKK,IGRC,IGRC,IL)=SCATW
                ENDIF
              ELSE
                SCATW=0.0D0
                DO 535 JGRC=1,NGCOND
                  SCATW=SCATW+SCATTD(IKK,JGRC,IGRC,IL)
 535            CONTINUE
                SIGS(IKK,IGRC,IL)=REAL(SCATW)
              ENDIF
 540        CONTINUE
 550      CONTINUE
 560    CONTINUE
      ENDIF
*----
*  FISSION SPECTRUM NORMALIZATION
*----
      IF((FXSOUR.EQ.0.0D0).AND.(TOTFIS.GT.0.0D0)) THEN
        CALL XDRSET(FLDMC,NMERGE*NGCOND,0.0)
        DO 580 IGRC=1,NGCOND
          DO 570 IFIS=1,NIFISS
            CALL EDIRAT(0,NREGIO,NBMIX,MATCOD,FFUEL(1,IFIS),VOLUME,
     >                  FLDMC(1,IGRC),SIGMA(0),IMERGE,NMERGE)
 570      CONTINUE
 580    CONTINUE
        DO 640 IKK=1,NMERGE
          TOTAL1=0.0
          DO 590 IGRC=1,NGCOND
            IF(RATECM(IKK,IGRC,NW+4).NE.0.0) THEN
             RATECM(IKK,IGRC,NW+8)=RATECM(IKK,IGRC,NW+4)/FLDMC(IKK,IGRC)
             TOTAL1=TOTAL1+RATECM(IKK,IGRC,NW+8)
            ELSE
             RATECM(IKK,IGRC,NW+8)=0.0
            ENDIF
 590      CONTINUE
          IF((IADJ.EQ.0).AND.(TOTAL1.NE.0.0)) THEN
            DO 600 IGRC=1,NGCOND
              RATECM(IKK,IGRC,NW+8)=RATECM(IKK,IGRC,NW+8)/TOTAL1
 600        CONTINUE
          ELSE IF(IADJ.EQ.1) THEN
            DO 601 IGRC=1,NGCOND
              RATECM(IKK,IGRC,NW+8)=RATECM(IKK,IGRC,NW+8)/
     >        FADJCM(IKK,IGRC,1)
 601        CONTINUE
          ENDIF
          DO 630 IDEL=1,NDEL
            K=12+NW+NDEL+IDEL
            TOTAL1=0.0
            DO 610 IGRC=1,NGCOND
              IF(RATECM(IKK,IGRC,K).NE.0.0) THEN
                RATECM(IKK,IGRC,K)=RATECM(IKK,IGRC,K)/FLDMC(IKK,IGRC)
                TOTAL1=TOTAL1+RATECM(IKK,IGRC,K)
              ELSE
                RATECM(IKK,IGRC,K)=0.0
              ENDIF
 610        CONTINUE
            IF((IADJ.EQ.0).AND.(TOTAL1.NE.0.0)) THEN
              DO 620 IGRC=1,NGCOND
                RATECM(IKK,IGRC,K)=RATECM(IKK,IGRC,K)/TOTAL1
 620          CONTINUE
            ELSE IF(IADJ.EQ.1) THEN
              DO 621 IGRC=1,NGCOND
                RATECM(IKK,IGRC,K)=RATECM(IKK,IGRC,K)/FADJCM(IKK,IGRC,1)
 621          CONTINUE
            ENDIF
 630      CONTINUE
 640    CONTINUE
      ENDIF
      DEALLOCATE(XSCAT,SIGMA)
      DEALLOCATE(IPOSCT,IFGSCT,INGSCT)
*----
*  CONDENSATION OF PHYSICAL ALBEDOS
*----
      IF(NALBP.GT.0) THEN
        ALLOCATE(ALBPGR(NALBP,NGROUP))
        CALL LCMGET(IPMACR,'ALBEDO',ALBPGR)
        IGRFIN=0
        DO 663 IGRC=1,NGCOND
        IGRDEB=IGRFIN+1
        IGRFIN=IGCOND(IGRC)
        DENOM=0.0
        DO 655 IGR=IGRDEB,IGRFIN
        DO 650 IREGIO=1,NREGIO
        DENOM=DENOM+FLXINT(IREGIO,IGR,1)
 650    CONTINUE
 655    CONTINUE
        DO 662 IAL=1,NALBP
        ALBP(IAL,IGRC)=0.0
        DO 661 IGR=IGRDEB,IGRFIN
        DO 660 IREGIO=1,NREGIO
        ALBP(IAL,IGRC)=ALBP(IAL,IGRC)+ALBPGR(IAL,IGR)*
     1                 FLXINT(IREGIO,IGR,1)/DENOM
 660    CONTINUE
 661    CONTINUE
 662    CONTINUE
 663    CONTINUE
        DEALLOCATE(ALBPGR)
      ENDIF
*----
*  PRINT REACTION RATES
*----
      ILEAKS=ILEAK2
      IF(IPRINT.GE.1) THEN
        CALL EDIPRR(IPRINT,NL,ITRANC,NGCOND,NMERGE,ILEAKS,NW,NTAUXT,
     >              B2,VOLMER,NENER,WENERG,RATECM,FLUXCM,SCATTD)
      ENDIF
*----
*  COMPUTE MERGED/CONDENSED X-S
*----
      CALL EDIPXS(IPEDIT,IADJ,IPRINT,NL,NDEL,NALBP,ITRANC,NSAVES,NGCOND,
     >            NMERGE,ILEAKS,NW,NTAUXT,EIGENK,B2,CUREIN,NIFISS,
     >            CURNAM,NEDMAC,VOLMER,WLETYC,WENERG,SCATTD,RATECM,
     >            FLUXCM,FADJCM,SIGS,SCATTS,DISFCT,ALBP,TAUXE,HVECT,
     >            OVERV,HFACT,NENER,TIMEF,LH)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ALBP,DECAY,HFACT,OVERV,FLDMC,FFUEL,TAUXE,DISFCT)
      DEALLOCATE(IFUELR,HVECT)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(/53H EDIDTX: *** WARNING *** NORMALIZATION OF THE WITHIN-,
     > 34HGROUP SCATTERING TRANSFER OF ORDER,I3,9H IN GROUP,I4,5H AND ,
     > 6HREGION,I5,3H BY,F6.2,3H %.)
      END
