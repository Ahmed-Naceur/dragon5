*DECK EDIDRV
      SUBROUTINE EDIDRV(IPEDIT,IPTRK1,IPFLUX,IPLIB,IPSYS,NGROUP,NBMIX,
     >                  NREGIO,MATCOD,VOLUME,KEYFLX,NIFISS,NEDMAC,NL,
     >                  NDEL,NALBP,ITRANC,NGCOND,NMERGE,IADF,IDFM,NW,
     >                  ICURR,IHF,IFFAC,ILUPS,NSAVES,NSTATS,IXEDI,
     >                  ISOTXS,IGCOND,IMERGE,CURNAM,OLDNAM,NBMICR,
     >                  CARISO,NACTI,IACTI,IPRINT,LISO,IADJ,NOUT,HVOUT,
     >                  BB2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for edition operations.
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
* IPTRK1  pointer to the reference tracking object.
* IPFLUX  pointer to the flux LCM object.
* IPLIB   pointer to the internal library or macrolib LCM object.
* IPSYS   pointer to the pij LCM object (only used with Selengut
*         normalization).
* NGROUP  number of energy groups.
* NBMIX   number of mixtures.
* NREGIO  number of regions.
* MATCOD  mixture index in region.
* VOLUME  volume of region.
* KEYFLX  average flux position per region.
* NIFISS  number of fissile isotopes.
* NEDMAC  number of extra macroscopic cross section types.
* NL      number of Legendre orders of the scattering cross sections.
* NDEL    number of delayed precursor groups.
* NALBP   number of physical albedos.
* ITRANC  type of transport correction.
* NGCOND  number of condensed groups.
* NMERGE  number of regions merged.
* IADF    flag for assembly discontinuity factors (ADF) information:
*         = 0 do not compute them;
*         = 1 compute them using ALBS information;
*         = 2 compute them using averaged fluxes in boundary regions;
*         = 3 compute them using SYBIL/ARM interface currents.
* IDFM    flag for ADF info in input macrolib (0/1/2: absent/present).
* NW      type of weighting for P1 cross section information:
*         = 0 P0; = 1 P1.
* ICURR   type of current approximation if NW=1:
*         =1: heterogeneous leakage;
*         =2: Todorova outscatter approximation;
*         =4: use higher spherical harmonic moments of flux.
* IHF     H-factor calculation flag:
*         = 0 no; = 1 yes.
* IFFAC   four factor calculation flag:
*         = 0 no four factors (defaut);
*         = 1 four factor evaluation.
* ILUPS   flag to remove up-scattering from output.
* NSAVES  homogenized cross section computation and saving:
*         = 0  no compute no save;
*         = 1  compute, no save;
*         = 2  compute, save.
* NSTATS  statistics level:
*         = 0  no stats;
*         = 1  statistics on fluxes
*         = 2  statistics on reaction rates;
*         = 3  statistics on fluxes and reaction rates;
*         =-1  delta sigma ('MERG COMP' only).
* IXEDI   first ISOTX mixture record number.
* ISOTXS  ISOTX file enabling flag (0: off; 1: binary; 2: ascii).
* IGCOND  condensed group limits.
* IMERGE  merged region positions.
* CURNAM  name of LCM directory where the current rates are to be
*         stored.
* OLDNAM  name of LCM directory where old rates were stored.
* NBMICR  type of microlib edition:
*         =-2: process only macroscopic residue;
*         =-1: process each isotope;
*         =0: process no isotope;
*         >0 number of isotopes to process.
* CARISO  names of the isotopes to process.
* NACTI   number of activation editions.
* IACTI   activation mixtures.
* IPRINT  print index.
* LISO    =.TRUE. if we want to keep all the isotopes after 
*         homogeneization.
* IADJ    type of flux weighting:
*         =0: direct flux weighting;
*         =1: direct-adjoint flux weighting.
* NOUT    number of output cross section types (set to zero to recover
*         all cross section types).
* HVOUT   MATXS names of the output cross section types.
* BB2     imposed leakege used in non-regression tests.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL       LISO
      TYPE(C_PTR)   IPEDIT,IPTRK1,IPFLUX,IPLIB,IPSYS
      INTEGER       NGROUP,NBMIX,NREGIO,MATCOD(NREGIO),KEYFLX(NREGIO),
     >              NIFISS,NEDMAC,NL,NDEL,NALBP,ITRANC,NGCOND,NMERGE,
     >              IADF,IDFM,NW,ICURR,IHF,IFFAC,ILUPS,NSAVES,NSTATS,
     >              IXEDI,ISOTXS,IGCOND(NGROUP),IMERGE(NREGIO),NBMICR,
     >              NACTI,IACTI(NBMIX),IPRINT,IADJ,NOUT
      REAL          VOLUME(NREGIO),BB2
      CHARACTER     CURNAM*12,OLDNAM*12,CARISO(NBMICR)*12,HVOUT(NOUT)*8
*----
*  LOCAL VARIABLES
*----
      PARAMETER    (IUNOUT=6,MAXED=100,NSTATE=40,IOUT=6)
      TYPE(C_PTR)   JPFLUX,JPFLUA,IPMIC2,IPMAC2,IPADF,JPLIB,KPLIB,
     >              KPEDIT,JPMAC2,KPMAC2
      CHARACTER     HSIGN*12,TEXT8*8,HVECT(MAXED)*8,NISEXT*6,NISOTX*12,
     >              CTITLE*72,NAMSBR*12,HTYPE*8,TEXT12*12,HSMG*131
      INTEGER       IFPAR(NSTATE),IPAR(NSTATE),IDIM(NSTATE)
      REAL          B2(4),B2T(3),TIMEF(3)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPE,MIX,IDEPL,ISONA,
     > ISONR,LSISO,INADPL,JPIFI,KDRI,INNAM,INNRF,NMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KEYANI
      REAL, ALLOCATABLE, DIMENSION(:) :: WORKF,WORKA,VOLME,WLETY,WE,
     > COURI,TAUXT,SIGT,SIGS,SCATS,FLINT,SCATD,DEN,TN,EMEVF,EMEVG,RER,
     > DECAY,YIELD,PYIEL,RRD,FIYI,ENERG,NAWR,NDEN,NTMP,NVOL,SNEJ,WORK1,
     > WORK2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FLUXC,FADJC,FLUXES,AFLUXE
      CHARACTER*8, ALLOCATABLE, DIMENSION(:) :: HADF
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO,JPISO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FLUXES(NREGIO,NGROUP,NW+1),
     >         AFLUXE(NREGIO,NGROUP,NW+1))
*----
*  FIND THE SIGNATURE OF IPLIB
*----
      CALL LCMGTC(IPLIB,'SIGNATURE',12,1,HSIGN)
*----
*  RECOVER NEUTRON FLUXES AND CURRENTS (IF ILEAKC.GE.5)
*----
      CALL LCMGET(IPFLUX,'STATE-VECTOR',IFPAR)
      IF(IFPAR(1).NE.NGROUP) CALL XABORT('EDIDRV: INVALID VALUE OF NGR'
     > //'OUP.')
      ITYPEC=IFPAR(6)
      ILEAKC=IFPAR(7)
      IF(ILEAKC.EQ.0) THEN
*       NO LEAKAGE
        ILEAKS=0
      ELSE IF(ILEAKC.LE.4) THEN
*       DIFFON-TYPE LEAKAGE
        ILEAKS=1
      ELSE IF(ILEAKC.EQ.5) THEN
*       ECCO-TYPE LEAKAGE (WITH ISOTROPIC STREAMING EFFECTS)
        ILEAKS=2
      ELSE IF(ILEAKC.GE.6) THEN
*       TIBERE-TYPE LEAKAGE (WITH ANISOTROPIC STREAMING EFFECTS)
        ILEAKS=3
      ENDIF
      IF(ITYPEC.GT.0) CALL LCMGET(IPFLUX,'K-INFINITY',CUREIN)
      CALL XDRSET(B2,4,0.0)
      IF(ITYPEC.GT.2) THEN
        CALL LCMGET(IPFLUX,'B2  B1HOM',B2(4))
        IF(ILEAKS.EQ.3) THEN
          CALL LCMGET(IPFLUX,'B2  HETE',B2)
          IF(B2(4).EQ.0.0) THEN
            B2T(1)=1.0/3.0
            B2T(2)=B2T(1)
            B2T(3)=B2T(1)
          ELSE
            B2T(1)=B2(1)/B2(4)
            B2T(2)=B2(2)/B2(4)
            B2T(3)=B2(3)/B2(4)
          ENDIF
        ENDIF
      ENDIF
      IF((NW.GE.1).AND.(ILEAKC.LE.4).AND.(ICURR.EQ.1)) THEN
        CALL XABORT('EDIDRV: CURRENT WEIHTING OF P1 XS INFO (NW=1) '
     >  //'IS ONLY AVAILABLE WITH A STREAMING-ENABLED LEAKAGE MODEL.')
      ENDIF
      IF(ILEAKC.EQ.4) THEN
        CALL XDRSET(B2,4,0.0)
        CALL XDRSET(B2T,3,0.0)
      ENDIF
      IF(IADJ.EQ.0) THEN
        CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
        IF(ILON.EQ.0) CALL XABORT('EDIDRV: MISSING FLUX INFO.')
        JPFLUX=LCMGID(IPFLUX,'FLUX')
        CALL LCMLEL(JPFLUX,1,NUN,ITYLCM)
      ELSE IF(IADJ.EQ.1) THEN
        CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
        IF(ILON.EQ.0) CALL XABORT('EDIDRV: MISSING FLUX INFO.')
        JPFLUX=LCMGID(IPFLUX,'FLUX')
        CALL LCMLEN(IPFLUX,'AFLUX',ILON,ITYLCM)
        IF(ILON.EQ.0) CALL XABORT('EDIDRV: MISSING ADJOINT FLUX INFO.')
        JPFLUA=LCMGID(IPFLUX,'AFLUX')
        CALL LCMLEL(JPFLUX,1,NUN,ITYLCM)
        ALLOCATE(WORKA(NUN))
      ELSE
        CALL XABORT('EDIDRV: INVALID VALUE OF IADJ.')
      ENDIF
      ALLOCATE(WORKF(NUN))
      DO IGR=1,NGROUP
        IF(IADJ.EQ.0) THEN
          CALL LCMGDL(JPFLUX,IGR,WORKF)
          DO IREG=1,NREGIO
            FLUXES(IREG,IGR,1)=WORKF(KEYFLX(IREG))
            AFLUXE(IREG,IGR,1)=1.0
          ENDDO
        ELSE IF(IADJ.EQ.1) THEN
          CALL LCMGDL(JPFLUX,IGR,WORKF)
          CALL LCMGDL(JPFLUA,IGR,WORKA)
          DO IREG=1,NREGIO
            FLUXES(IREG,IGR,1)=WORKF(KEYFLX(IREG))
            AFLUXE(IREG,IGR,1)=WORKA(KEYFLX(IREG))
          ENDDO
        ENDIF
        IF((ICURR.EQ.1).AND.(ILEAKS.EQ.2)) THEN
*         ISOTROPIC STREAMING (ECCO)
          IF(NW.NE.1) CALL XABORT('EDIDRV: NW=1 EXPECTED(1).')
          DO IREG=1,NREGIO
            FLUXES(IREG,IGR,2)=WORKF(NUN/2+KEYFLX(IREG))
          ENDDO
        ELSE IF((ICURR.EQ.1).AND.(ILEAKS.EQ.3)) THEN
*         ANISOTROPIC STREAMING
          IF(NW.NE.1) CALL XABORT('EDIDRV: NW=1 EXPECTED(2).')
          DO IREG=1,NREGIO
            CURN=0.0
            DO IDIR=1,3
              CURN=CURN+B2T(IDIR)*WORKF(IDIR*NUN/4+KEYFLX(IREG))
            ENDDO
            FLUXES(IREG,IGR,2)=CURN
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(WORKF)
      IF(IADJ.EQ.1) DEALLOCATE(WORKA)
*----
*  COMPUTE HIGHER MOMENT FLUXES IF NW=1
*----
      IF(ICURR.EQ.2) THEN
*        Outscatter Todorova approximation
         IF(NW.NE.1) CALL XABORT('EDIDRV: NW=1 EXPECTED(3).')
         IF(HSIGN.EQ.'L_LIBRARY') CALL LCMSIX(IPLIB,'MACROLIB',1)
         JPLIB=LCMGID(IPLIB,'GROUP')
         ALLOCATE(SIGT(0:NBMIX),SIGS(0:NBMIX))
         DO IGR=1,NGROUP
            KPLIB=LCMGIL(JPLIB,IGR)
            SIGT(0)=0.0
            SIGS(0)=0.0
            CALL LCMGET(KPLIB,'NTOT0',SIGT(1))
            CALL LCMGET(KPLIB,'SIGS01',SIGS(1))
            DO IREG=1,NREGIO
              IBM=MATCOD(IREG)
              IF(IBM.GT.0) THEN
                FACT=3.0*(SIGT(IBM)-SIGS(IBM))
                IF(FACT.EQ.0.0) CALL XABORT('EDIDRV: DIVIDE CHECK.')
                FLUXES(IREG,IGR,2)=FLUXES(IREG,IGR,1)/FACT
                IF(IADJ.EQ.1) AFLUXE(IREG,IGR,2)=AFLUXE(IREG,IGR,1)/FACT
              ELSE
                FLUXES(IREG,IGR,2)=FLUXES(IREG,IGR,1)
                IF(IADJ.EQ.1) AFLUXE(IREG,IGR,2)=AFLUXE(IREG,IGR,1)
              ENDIF
            ENDDO
         ENDDO
         DEALLOCATE(SIGS,SIGT)
         IF(HSIGN.EQ.'L_LIBRARY') CALL LCMSIX(IPLIB,' ',2)
      ELSE IF(ICURR.EQ.4) THEN
*        Use higher spherical harmonic moments
         IF(NW.EQ.0) CALL XABORT('EDIDRV: NW>0 EXPECTED(5).')
         CALL LCMGTC(IPTRK1,'TRACK-TYPE',12,1,TEXT12)
         IF(TEXT12.EQ.'MCCG') THEN
           CALL LCMGET(IPTRK1,'STATE-VECTOR',IPAR)
           NDIM=IPAR(16)
           CALL LCMGET(IPTRK1,'MCCG-STATE',IPAR)
           NFUNL=IPAR(19)
           NLIN=IPAR(20)
         ELSE IF(TEXT12.EQ.'SN') THEN
           CALL LCMGET(IPTRK1,'STATE-VECTOR',IPAR)
           NFUNL=IPAR(7)
           NLIN=IPAR(8)
           NDIM=IPAR(9)
           NLIN=NLIN**NDIM
         ELSE
           CALL XABORT('EDIDRV: MCCG OR SN TRACKING EXPECTED WITH '
     >     //'P1W_SP OPTION')
         ENDIF
         ALLOCATE(KEYANI(NREGIO,NLIN,NFUNL))
         CALL LCMGET(IPTRK1,'KEYFLX$ANIS',KEYANI)
         CALL EDIWP1(IPFLUX,NW,NGROUP,NUN,NREGIO,NDIM,IADJ,NLIN,
     >   NFUNL,NGCOND,NMERGE,KEYANI,VOLUME,IGCOND,IMERGE,FLUXES(1,1,2),
     >   AFLUXE(1,1,2))
         DEALLOCATE(KEYANI)
      ENDIF
*----
*  ALLOCATE MEMORY FOR GROUP CONDENSATION AND MERGE
*----
      ALLOCATE(VOLME(NMERGE),WLETY(NGCOND),WE(NGCOND+1))
      NELEMT=NMERGE*NGCOND
*----
*  COMPUTE REACTION RATES FOR THE EDITION MACROLIB
*----
      NTAUXT=12+NW+2*NDEL
      ALLOCATE(FLUXC(NMERGE,NGCOND,NW+1),FADJC(NMERGE,NGCOND,NW+1),
     > TAUXT(NTAUXT*NELEMT),SIGS(NL*NELEMT),SCATS(NELEMT*NGCOND*NL),
     > FLINT(NREGIO*NGROUP*(NW+1)),SCATD(2*NELEMT*NGCOND*NL))
      NBISO=0
      CALL LCMLEN(IPFLUX,'K-EFFECTIVE',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.1) THEN
        CALL LCMGET(IPFLUX,'K-EFFECTIVE',EIGENK)
      ELSE
        EIGENK=1.0
      ENDIF
      CALL LCMLEN(IPFLUX,'K-INFINITY',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.1) THEN
        CALL LCMGET(IPFLUX,'K-INFINITY',EIGINF)
      ELSE
        EIGINF=EIGENK
      ENDIF
      TIMEF(1)=0.0
      TIMEF(2)=0.0
      TIMEF(3)=0.0
      IF(HSIGN.EQ.'L_LIBRARY') THEN
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL LCMLEN(IPLIB,'TIMESTAMP',ILCMLN,ILCMTY)
        IF((ILCMLN.GE.1).AND.(ILCMLN.LE.3)) THEN
           CALL LCMGET(IPLIB,'TIMESTAMP',TIMEF)
        ENDIF
      ENDIF
      CALL EDIDTX(IPEDIT,IPFLUX,IPLIB,IADJ,IPRINT,NL,NDEL,NALBP,ITRANC,
     >            NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,ILEAKS,ILUPS,NW,
     >            MATCOD,VOLUME,KEYFLX,IGCOND,IMERGE,FLUXES,AFLUXE,
     >            EIGENK,VOLME,WLETY,WE,TAUXT,FLUXC,FADJC,FLINT,SCATD,
     >            SCATS,NIFISS,NSAVES,CURNAM,NEDMAC,SIGS,B2,CUREIN,
     >            TIMEF,NTAUXT)
      IF(HSIGN.EQ.'L_LIBRARY') CALL LCMSIX(IPLIB,' ',2)
      DEALLOCATE(SCATD,FLINT)
*----
*  COMPUTE BOUNDARY EDITIONS FOR ADF OR SPH WITH SELENGUT
*----
      IF(CURNAM.NE.' ') THEN
         IF(IPRINT.GT.0) WRITE(IOUT,'(30H EDIDRV: EDITION DIRECTORY IS ,
     >   A)') CURNAM
         IF(HSIGN.EQ.'L_LIBRARY') CALL LCMSIX(IPLIB,'MACROLIB',1)
         IPMAC2=LCMDID(IPEDIT,CURNAM)
         IPMAC2=LCMDID(IPMAC2,'MACROLIB')
         IF(IADF.EQ.1) THEN
*           recover outgoing current from escape probabilities
            CALL EDIALB(IPMAC2,IPFLUX,IPLIB,IPSYS,IPRINT,NBMIX,NW,
     >      B2,NGROUP,NIFISS,NGCOND,ITRANC,ILEAKS,NREGIO,MATCOD,
     >      VOLUME,KEYFLX,IGCOND,FLUXES)
          ELSE IF((IADF.EQ.2).OR.(IADF.EQ.-2)) THEN
            ALLOCATE(WORKF(NGCOND))
            IF(IADF.EQ.-2) THEN
*             recover averaged fluxes used to compute ADF
              DO IGR=1,NGCOND
                WORKF(IGR)=SUM(FLUXC(:,IGR,1))/SUM(VOLME(:))
              ENDDO
            ELSE
              WORKF(:NGCOND)=1.0
            ENDIF
*           use averaged fluxes obtained over boundary regions
            IPADF=LCMGID(IPEDIT,'REF:ADF')
            CALL LCMGET(IPADF,'NTYPE',NTYPE)
            IF(NTYPE.EQ.0) CALL XABORT('EDIADF: NTYPE=0.')
            CALL LCMSIX(IPMAC2,'ADF',1)
            ALLOCATE(HADF(NTYPE),COURI(NGCOND))
            CALL LCMGTC(IPADF,'HADF',8,NTYPE,HADF)
            DO IT=1,NTYPE
              HTYPE=HADF(IT)
              CALL EDIGAP(IPADF,HTYPE,NGROUP,NGCOND,NREGIO,VOLUME,
     >        IGCOND,FLUXES,WORKF,IPRINT,COURI)
              CALL LCMPUT(IPMAC2,HTYPE,NGCOND,2,COURI)
            ENDDO
            DEALLOCATE(WORKF)
            CALL LCMPUT(IPMAC2,'NTYPE',1,1,NTYPE)
            CALL LCMPTC(IPMAC2,'HADF',8,NTYPE,HADF)
            DEALLOCATE(COURI,HADF)
            CALL LCMSIX(IPMAC2,' ',2)
         ELSE IF(IADF.EQ.3) THEN
*           recover outgoing current from interface currents in Eurydice
            CALL LCMGTC(IPTRK1,'TRACK-TYPE',12,1,TEXT12)
            IF(TEXT12.EQ.'SYBIL') THEN
              CALL EDIJO1(IPMAC2,IPTRK1,IPFLUX,IPRINT,NGCOND,IGCOND)
            ELSE IF(TEXT12.EQ.'MCCG') THEN
              CALL EDIJO2(IPMAC2,IPTRK1,IPFLUX,IPRINT,NGCOND,IGCOND)
            ELSE
              WRITE(HSMG,'(40HEDIDRV: INCOMPATIBLE SOLUTION TYPE. SYBI,
     >        20HL OR MCCG EXPECTED. ,A12,6HFOUND.)') TEXT12
              CALL XABORT(HSMG)
            ENDIF
         ELSE IF(IADF.EQ.4) THEN
*           recover ADF information from input macrolib
            CALL LCMLEN(IPLIB,'GROUP',ILCMLN,ITYLCM)
            IF(ILCMLN.NE.NGCOND) CALL XABORT('EDIDRV: UNABLE TO RECOVE'
     >      //'R ADF INFORMATION FROM INPUT MACROLIB.')
            CALL LCMSIX(IPMAC2,'ADF',1)
            CALL LCMSIX(IPLIB,'ADF',1)
            CALL LCMEQU(IPLIB,IPMAC2)
            CALL LCMSIX(IPLIB,' ',2)
            CALL LCMSIX(IPMAC2,' ',2)
         ENDIF
         IF(HSIGN.EQ.'L_LIBRARY') CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  RECOVER ISOTOPIC INFORMATION FROM THE MICROLIB
*----
      IF(HSIGN.EQ.'L_LIBRARY') THEN
        CALL LCMGET(IPLIB,'STATE-VECTOR',IPAR)
        NBISO=IPAR(2)
        NED=IPAR(13)
        NBESP=IPAR(16)
        IF(NBISO.EQ.0) CALL XABORT('EDIDRV: NO ISOTOPES FOUND.')
        ALLOCATE(DEN(NBISO),ITYPE(NBISO),MIX(NBISO),TN(NBISO),
     >  IDEPL(NBISO),ISONA(3*NBISO),ISONR(3*NBISO),LSISO(NBISO),
     >  IPISO(NBISO))
        CALL LCMGET(IPLIB,'ISOTOPESDENS',DEN)
        CALL LCMGET(IPLIB,'ISOTOPESTYPE',ITYPE)
        CALL LCMGET(IPLIB,'ISOTOPESMIX',MIX)
        CALL LCMGET(IPLIB,'ISOTOPESTEMP',TN)
        CALL LCMGET(IPLIB,'ISOTOPESTODO',IDEPL)
        IF(NED.GT.0) CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED,HVECT)
        CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONA)
        CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONR)
        CALL XDISET(LSISO,NBISO,0)
        IF(NBMICR.EQ.-2) THEN
          CALL XDISET(LSISO,NBISO,0)
        ELSE IF(NBMICR.EQ.-1) THEN
          CALL XDISET(LSISO,NBISO,1)
        ELSE IF(NBMICR.GT.0) THEN
          DO IISO=1,NBISO
            WRITE(TEXT8,'(2A4)') (ISONA((IISO-1)*3+I0),I0=1,2)
            DO IIII=1,NBMICR
              IF(CARISO(IIII)(1:8).EQ.TEXT8) LSISO(IISO)=1
            ENDDO
          ENDDO
        ENDIF
*----
*  SET THE LCM MICROLIB ISOTOPEWISE DIRECTORIES.
*----
        CALL LIBIPS(IPLIB,NBISO,IPISO)
      ENDIF
*----
*  EVALUATE H-FACTOR IF REQUIRED FOR THE EDITION MACROLIB
*----
      ALLOCATE(EMEVF(NBISO),EMEVG(NBISO))
      CALL XDRSET(EMEVF,NBISO,0.0)
      CALL XDRSET(EMEVG,NBISO,0.0)
      IF((NSAVES.GE.2).AND.(IHF.NE.0)) THEN
        CALL LCMLEN(IPLIB,'DEPL-CHAIN',ILLCM,ITLCM)
        IF(ILLCM.NE.0) THEN
          CALL LCMSIX(IPEDIT,CURNAM,1)
          CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
          CALL LCMGET(IPLIB,'STATE-VECTOR',IDIM)
          NDEPL=IDIM(1)
          NREAC=IDIM(8)
*
          ALLOCATE(INADPL(3*NDEPL),RER(NREAC*NDEPL))
          CALL LCMGET(IPLIB,'ISOTOPESDEPL',INADPL)
          CALL LCMGET(IPLIB,'DEPLETE-ENER',RER)
          CALL LCMSIX(IPLIB,' ',2)
*
          CALL EDIHFC(IPEDIT,NGROUP,NGCOND,NREGIO,NMERGE,NBISO,NDEPL,
     >                NREAC,MATCOD,VOLUME,INADPL,ISONA,ISONR,IPISO,
     >                MIX,FLUXES(1,1,1),DEN,IGCOND,IMERGE,RER,EMEVF,
     >                EMEVG,VOLME,IPRINT)
*
          DEALLOCATE(RER,INADPL)
          CALL LCMSIX(IPEDIT,' ',2)
        ENDIF
      ENDIF
*----
*  LUMP THE DEPLETION CHAIN
*----
      ALLOCATE(DECAY(NBISO),YIELD((NGCOND+1)*NBISO*NMERGE))
      CALL XDRSET(DECAY,NBISO,0.0)
      CALL XDRSET(YIELD,(NGCOND+1)*NBISO*NMERGE,0.0)
      NDEPL=0
      NDFI=0
      IF((NBMICR.NE.0).AND.(NBISO.NE.0)) THEN
        CALL LCMSIX(IPEDIT,CURNAM,1)
        CALL LCMLEN(IPLIB,'DEPL-CHAIN',ILCMLN,ITYLCM)
        IF((ILCMLN.NE.0).AND.(CURNAM.NE.' ')) THEN
          CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
          CALL EDIDEP(IPRINT,IPLIB,IPEDIT,NBISO,ISONR,LSISO,IDEPL,LISO,
     >    NBCH)
          CALL LCMSIX(IPLIB,' ',2)
        ENDIF
*----
*  RECOVER DEPLETION INFORMATION FROM THE INTERNAL LIBRARY
*----
        CALL LCMLEN(IPEDIT,'DEPL-CHAIN',ILLCM,ITLCM)
        IF(ILLCM.NE.0) THEN
          CALL LCMSIX(IPEDIT,'DEPL-CHAIN',1)
          CALL LCMGET(IPEDIT,'STATE-VECTOR',IDIM)
          NDEPL=IDIM(1)
          NDFI=IDIM(2)
          NDFP=IDIM(3)
          NREAC=IDIM(8)
          ALLOCATE(JPIFI(NDFI*NMERGE),PYIEL(NDFI*NBISO*NMERGE))
*
          ALLOCATE(INADPL(3*NDEPL),KDRI(NREAC*NDEPL),RRD(NDEPL),
     >    FIYI(NDFI*NDFP))
          CALL LCMGET(IPEDIT,'ISOTOPESDEPL',INADPL)
          CALL LCMGET(IPEDIT,'DEPLETE-REAC',KDRI)
          CALL LCMGET(IPEDIT,'DEPLETE-DECA',RRD)
          IF(NDFI*NDFP.GT.0) THEN
             CALL LCMGET(IPEDIT,'FISSIONYIELD',FIYI)
          ENDIF
          CALL LCMSIX(IPEDIT,' ',2)
*
          CALL EDIHFD(NGROUP,NGCOND,NREGIO,NMERGE,NBISO,NDEPL,NDFI,
     >                NDFP,NREAC,MATCOD,VOLUME,INADPL,ISONA,ISONR,
     >                IPISO,MIX,FLUXES(1,1,1),DEN,IDEPL,IGCOND,IMERGE,
     >                KDRI,RRD,FIYI,DECAY,YIELD,JPIFI,PYIEL)
*
          DEALLOCATE(FIYI,RRD,KDRI,INADPL)
        ENDIF
        CALL LCMSIX(IPEDIT,' ',2)
*----
*  COMPUTE MICROSCOPIC CROSS SECTIONS
*----
        CALL EDIMIC(IPEDIT,IPFLUX,IPLIB,IADJ,NL,NDEL,NBESP,NBISO,NDEPL,
     >              ISONA,ISONR,IPISO,MIX,TN,NED,HVECT,NOUT,HVOUT,
     >              IPRINT,NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,NDFI,
     >              ILEAKS,ILUPS,NW,MATCOD,VOLUME,KEYFLX,CURNAM,
     >              IGCOND,IMERGE,FLUXES,AFLUXE,EIGENK,EIGINF,B2,DEN,
     >              ITYPE,IDEPL,LSISO,EMEVF,EMEVG,DECAY,YIELD,JPIFI,
     >              PYIEL,ITRANC,LISO)
*----
*  ISOTX FILE PROCESSING
*----
        IF(ISOTXS.GE.1) THEN
          CALL LCMSIX(IPEDIT,CURNAM,1)
          CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
          NBNISO=IPAR(2)
          NAMSBR='EDIDRV'
          IF(IPRINT.GE.1) WRITE(IOUT,6000) NAMSBR
          ALLOCATE(INNAM(3*NBNISO),INNRF(3*NBNISO),NMIX(NBNISO))
          ALLOCATE(ENERG(NGCOND+1),NAWR(NBNISO),NDEN(NBNISO),
     >    NTMP(NBNISO),NVOL(NBNISO),SNEJ(NBNISO),JPISO(NBNISO))
          CALL LCMGET(IPEDIT,'ENERGY',ENERG)
          CALL LCMGET(IPEDIT,'ISOTOPESUSED',INNAM)
          CALL LCMGET(IPEDIT,'ISOTOPERNAME',INNRF)
          CALL LCMGET(IPEDIT,'ISOTOPESMIX',NMIX)
          CALL LCMGET(IPEDIT,'ISOTOPESDENS',NDEN)
          CALL LCMGET(IPEDIT,'ISOTOPESTEMP',NTMP)
          CALL LCMGET(IPEDIT,'ISOTOPESVOL',NVOL)
          CALL LIBIPS(IPEDIT,NBNISO,JPISO)
          DO ISO=1,NBNISO
            KPEDIT=JPISO(ISO)
            CALL LCMGET(KPEDIT,'AWR',AWR)
            EMEVF2=0.0
            EMEVG2=0.0
            CALL LCMLEN(KPEDIT,'MEVF',ILENF,ITYLCM)
            CALL LCMLEN(KPEDIT,'MEVG',ILENG,ITYLCM)
            IF(ILENF.EQ.1) CALL LCMGET(KPEDIT,'MEVF',EMEVF2)
            IF(ILENG.EQ.1) CALL LCMGET(KPEDIT,'MEVG',EMEVG2)
            NAWR(ISO)=AWR
            SNEJ(ISO)=EMEVF2+EMEVG2
          ENDDO
*
          NBIXS=IXEDI
          DO IMRG=1,NMERGE
            NBIXS=NBIXS+1
            WRITE(NISEXT,'(I6)') NBIXS
            DO ICAR=1,6
              IF(NISEXT(ICAR:ICAR) .EQ. ' ' .OR.
     >           NISEXT(ICAR:ICAR) .EQ. '*') THEN
                 NISEXT(ICAR:ICAR)='0'
              ENDIF
            ENDDO
            NISOTX='ISOTXS'//NISEXT
*----
*  GENERATE ONE ISOTXS FILE FOR EACH MERGED REGION IN EACH MIXTURE
*----
            WRITE(CTITLE,9000) NAMSBR,CURNAM,
     >                        'MICR        ','MIX',IMRG,NISOTX
            IF(IPRINT.GE.1) WRITE(IOUT,6002) IMRG,NISOTX
            IUTYPE=ISOTXS+1
            IWGOXS=KDROPN(NISOTX,0,IUTYPE,0)
            CALL EDITXS(IWGOXS,IUTYPE,IPRINT,NGCOND,NL,NBNISO,CTITLE,
     >                  IMRG,ENERG,INNAM,INNRF,JPISO,NMIX,NAWR,NDEN,
     >                  NTMP,SNEJ)
            IRETRN=KDRCLS(IWGOXS,1)
          ENDDO
*
          DEALLOCATE(JPISO,SNEJ,NVOL,NTMP,NDEN,NAWR,ENERG)
          DEALLOCATE(NMIX,INNRF,INNAM)
          CALL LCMSIX(IPEDIT,' ',2)
        ENDIF
      ENDIF
*----
*  COMPUTE MACROSCOPIC RESIDUAL CROSS SECTIONS
*----
      IF((NBMICR.NE.0).AND.(NBMICR.NE.-1).AND.(NBISO.NE.0).AND.
     >   (CURNAM.NE.' ')) THEN
        IPRIN2=IPRINT-1
        CALL EDIRES(IPEDIT,IPFLUX,IPLIB,IADJ,NL,NDEL,NBESP,NBISO,NDEPL,
     >             ISONA,ISONR,IPISO,MIX,TN,NED,HVECT,NOUT,HVOUT,IPRIN2,
     >             NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,NDFI,ILEAKS,ILUPS,
     >             NW,MATCOD,VOLUME,KEYFLX,CURNAM,IGCOND,IMERGE,
     >             FLUXES,AFLUXE,EIGENK,EIGINF,B2,DEN,ITYPE,IDEPL,LSISO,
     >             EMEVF,EMEVG,DECAY,YIELD,JPIFI,PYIEL,ITRANC,LISO)
      ENDIF
*----
*  EDIT MICROSCOPIC ACTIVATION XS
*----
      IF(NACTI.GT.0) THEN
        CALL EDIACT(IPEDIT,IPRINT,NGROUP,NGCOND,NREGIO,NMERGE,NL,NBISO,
     >              NED,VOLUME,MIX,IGCOND,IMERGE,FLUXES(1,1,1),ITRANC,
     >              ISONA,IPISO,HVECT,CURNAM,NACTI,IACTI,EMEVF,EMEVG)
      ENDIF
*----
*  STATISTICS AND DELTA SIGMAS
*----
      IF(NSTATS.NE.0) THEN
        CALL EDIDST(IPEDIT,IPRINT,NL,NGCOND,NMERGE,NSTATS,ILEAKS,
     >              EIGENK,B2,VOLME,WLETY,TAUXT,FLUXC,SCATS,OLDNAM,
     >              NW,NTAUXT)
      ENDIF
*----
*  FOUR FACTORS
*----
      IF(IFFAC.NE.0) THEN
        CALL EDIBAL(IPEDIT,IPFLUX,IPRINT,NL,IFFAC,NGCOND,NMERGE,EIGENK,
     >              TAUXT,FLUXC,SCATS,ILEAKS,B2,NW,NTAUXT)
      ENDIF
*
      IF(NDFI.GT.0) DEALLOCATE(PYIEL,JPIFI)
      DEALLOCATE(YIELD,DECAY)
      DEALLOCATE(EMEVG,EMEVF)
      DEALLOCATE(SCATS,SIGS,FADJC,FLUXC,TAUXT)
      DEALLOCATE(WE,WLETY,VOLME)
      IF(HSIGN.EQ.'L_LIBRARY') THEN
         DEALLOCATE(IPISO,ISONR,ISONA,IDEPL,TN,MIX,ITYPE,DEN,LSISO)
      ENDIF
*----
*  SET IADF IN MACROLIB AND MICROLIB STATE VECTORS
*----
      IF((CURNAM.NE.' ').AND.(IADF.NE.0)) THEN
         IPMIC2=LCMDID(IPEDIT,CURNAM)
         IPMAC2=LCMDID(IPMIC2,'MACROLIB')
         CALL LCMLEN(IPMAC2,'ADF',ILCMLN,ITYLCM)
         IF(ILCMLN.NE.0) THEN
            IF(IADF.EQ.4) THEN
              JADF=IDFM
            ELSE
              JADF=0
              CALL LCMSIX(IPMAC2,'ADF',1)
              CALL LCMLEN(IPMAC2,'ALBS00',ILCMLN,ITYLCM)
              IF(ILCMLN.NE.0) JADF=1
              CALL LCMLEN(IPMAC2,'HADF',ILCMLN,ITYLCM)
              IF((IADF.EQ.2).AND.(ILCMLN.NE.0)) JADF=2
              IF((IADF.EQ.-2).AND.(ILCMLN.NE.0)) JADF=3
              CALL LCMSIX(IPMAC2,' ',2)
            ENDIF
            CALL LCMGET(IPMAC2,'STATE-VECTOR',IPAR)
            IPAR(12)=JADF
            CALL LCMPUT(IPMAC2,'STATE-VECTOR',NSTATE,1,IPAR)
            IF((NBMICR.NE.0).AND.(HSIGN.EQ.'L_LIBRARY')) THEN
               CALL LCMGET(IPMIC2,'STATE-VECTOR',IPAR)
               IPAR(24)=JADF
               CALL LCMPUT(IPMIC2,'STATE-VECTOR',NSTATE,1,IPAR)
            ENDIF
         ENDIF
      ENDIF
*----
*  INCLUDE LEAKAGE IN THE MACROLIB (USED ONLY FOR NON-REGRESSION TESTS)
*----
      IF(BB2.NE.0.0) THEN
        IF(IPRINT.GT.0) WRITE(6,'(/32H EDIDRV: INCLUDE LEAKAGE IN THE ,
     >  13HMACROLIB (B2=,1P,E12.5,2H).)') BB2
        IPMIC2=LCMGID(IPEDIT,CURNAM)
        IPMAC2=LCMGID(IPMIC2,'MACROLIB')
        JPMAC2=LCMGID(IPMAC2,'GROUP')
        ALLOCATE(WORK1(NMERGE),WORK2(NMERGE))
        DO IGR=1,NGCOND
          KPMAC2=LCMGIL(JPMAC2,IGR)
          CALL LCMGET(KPMAC2,'DIFF',WORK1)
          CALL LCMGET(KPMAC2,'NTOT0',WORK2)
          WORK2(:NMERGE)=WORK2(:NMERGE)+BB2*WORK1(:NMERGE)
          CALL LCMPUT(KPMAC2,'NTOT0',NMERGE,2,WORK2)
        ENDDO
        DEALLOCATE(WORK2,WORK1)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(AFLUXE,FLUXES)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(1X,A6,': GENERATING ISOTXS FILE ')
 6002 FORMAT(8X,' FOR EDITING MIXTURE    = ',I6,
     >          ' INFORMATION STORED ON FILE = ',A12)
 9000 FORMAT(1X,A6,3X,A12,3X,A12,3X,A4,I6,5X,A12)
      END
