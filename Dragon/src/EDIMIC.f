*DECK EDIMIC
      SUBROUTINE EDIMIC(IPEDIT,IPFLUX,IPLIB,IADJ,NL,NDEL,NBESP,NBISO,
     1 NDEPL,ISONAM,ISONRF,IPISO,MIX,TN,NED,HVECT,NOUT,HVOUT,IPRINT,
     2 NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,NDFI,ILEAKS,ILUPS,NW,MATCOD,
     3 VOLUME,KEYFLX,CURNAM,IGCOND,IMERGE,FLUXES,AFLUXE,EIGENK,EIGINF,
     4 B2,DEN,ITYPE,IEVOL,LSISO,EMEVF,EMEVG,DECAY,YIELD,IPIFI,PYIELD,
     5 ITRANC,LISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Homogenization and condensation of microscopic cross sections.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPEDIT  pointer to the edition LCM object (L_EDIT signature).
* IPFLUX  pointer to the solution LCM object (L_FLUX signature).
* IPLIB   pointer to the reference microscopic cross section library
*         LCM object (L_LIBRARY signature).
* IADJ    type of flux weighting:
*         =0: direct flux weighting;
*         =1: direct-adjoint flux weighting.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NDEL    number of delayed precursor groups.
* NBESP   number of energy-dependent fission spectra.
* NBISO   number of isotopes.
* NDEPL   number of depleting isotopes.
* ISONAM  local names of NBISO isotopes:
*         chars 1 to 8  is the local isotope name;
*         chars 9 to 12 is a suffix function of the mix number.
* ISONRF  library name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture number associated with each isotope.
* TN      absolute temperature associated with each isotope.
* NED     number of extra vector edits from MATXS.
* HVECT   MATXS names of the extra vector edits.
* NOUT    number of output cross section types (set to zero to recover
*         all cross section types).
* HVOUT   MATXS names of the output cross section types.
* IPRINT  print index.
* NGROUP  number of energy groups.
* NGCOND  number of condensed groups.
* NBMIX   number of mixtures.
* NREGIO  number of volumes.
* NMERGE  number of merged regions.
* NDFI    number of fissile isotopes.
* ILEAKS  leakage calculation type: =0: no leakage; =1: homogeneous
*         leakage (Diffon); =2: isotropic streaming (Ecco);
*         =3: anisotropic streaming (Tibere).
* ILUPS   up-scattering removing flag (=1 to remove up-scattering from
*         output cross-sections).
* NW      type of weighting for P1 cross section info (=0: P0 ; =1: P1).
* MATCOD  mixture index per volume.
* VOLUME  volumes.
* KEYFLX  position of average fluxes.
* CURNAM  name of the LCM directory where the microscopic cross sections
*         are stored (a blank value means no save).
* IGCOND  limits of condensed groups.
* IMERGE  index of merged regions.
* FLUXES  fluxes.
* AFLUXE  adjoint fluxes.
* EIGENK  effective multiplication factor.
* EIGINF  infinite multiplication factor.
* B2      bucklings.
* DEN     number density of each isotope.
* ITYPE   type of each isotope.
* IEVOL   flag making an isotope non-depleting. A value of
*         1 is used to force an isotope to be non-depleting.
* LSISO   flag for isotopes saved.
* EMEVF   fission production energy.
* EMEVG   capture production energy.
* DECAY   radioactive decay constant.
* YIELD   group-ordered condensed fission product yield.
* IPIFI   isotope index associated with each fissile isotope
*         in microlib.
* PYIELD  fissile isotope ordered condensed fission product yield.
* ITRANC  type of transport correction (=0: no correction).
* LISO    =.TRUE. if we want to keep all the isotopes after 
*         homogeneization.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT,IPFLUX,IPLIB,IPISO(NBISO)
      INTEGER   IADJ,NL,NDEL,NBESP,NBISO,NDEPL,ISONAM(3,NBISO),
     1          ISONRF(3,NBISO),MIX(NBISO),NED,NOUT,IPRINT,NGROUP,
     2          NGCOND,NBMIX,NREGIO,NMERGE,NDFI,ILEAKS,ILUPS,NW,
     3          MATCOD(NREGIO),KEYFLX(NREGIO),IGCOND(NGCOND),
     4          IMERGE(NREGIO),ITYPE(NBISO),IEVOL(NBISO),LSISO(NBISO),
     5          IPIFI(NDFI,NMERGE),ITRANC
      REAL      TN(NBISO),VOLUME(NREGIO),FLUXES(NREGIO,NGROUP,NW+1),
     1          AFLUXE(NREGIO,NGROUP,NW+1),EIGENK,EIGINF,B2(4),
     2          DEN(NBISO),EMEVF(NBISO),EMEVG(NBISO),DECAY(NBISO),
     3          YIELD(NGCOND+1,NBISO,NMERGE),PYIELD(NDFI,NBISO,NMERGE)
      CHARACTER HVECT(NED)*8,HVOUT(NOUT)*8,CURNAM*12
      LOGICAL   LISO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXESP=4)
      TYPE(C_PTR) JPLIB,KPLIB,JPFLUX,JPEDIT,KPEDIT
      LOGICAL    LOGIC,LSTRD,LAWR,LMEVF,LMEVG,LDECA,LWD,LFIS,LONE
      CHARACTER  CM*2,HNEW*12,TEXT8*8,TEXT12*12,HSMG*131,HNAMIS*12
      INTEGER    IPAR(NSTATE),IESP(MAXESP+1)
      REAL       B2T(3),EESP(MAXESP+1)
      DOUBLE PRECISION TMP,PARM0,PARM3,PARM4,VOLMER,DDEN,DDENZ,SQFMAS,
     1           XDRCST,NMASS,EVJ,CONV,ZNU,ZDEN,ZFL1,ZFL2,DENVOL
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISMIX,ISTYP,ISTOD,ITYPRO,
     1 JPIFI,MILVO,ITYPS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHNISO
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IGAR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK
      REAL, ALLOCATABLE, DIMENSION(:) :: XSECT,WSTRD,SDEN,VOLISO,TNISO,
     1 TMPXS,WDLA,WORK,WORKF,ENR,GA1,GA2,VOLM
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR,WGAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PARM12,PHIAV,AHIAV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: GAS
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: WSCAT
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HMAKE
*----
*  SCRATCH STORAGE ALLOCATION
*----
*
* GAR/GAS CONTENTS:
*              1 : 'NWT0'  |  P0 direct flux
*              2 : 'NWT1'  |  P1 direct flux / NW values
*  ...                     |
*           2+NW : 'NWAT0' |  P0 adjoint flux
*           3+NW : 'NWAT1' |  P1 adjoint flux / NW values
*  ...                     |
*         3+2*NW : 'NTOT0' |  P0 total cross section
*         4+2*NW : 'NTOT1' |  P1 total cross section / NW values
*  ...                     |
*         4+3*NW : 'SIGS00' |
*         5+3*NW : 'SIGS01' | NL VALUES
*  ...                      |
*      4+NL+3*NW : 'NUSIGF'
*      5+NL+3*NW : HVECT(1) |
*      6+NL+3*NW : HVECT(2) | NED VALUES
*  ...                      |
*  5+NED+NL+3*NW : 'H-FACTOR'
*  6+NED+NL+3*NW : 'OVERV'
*  7+NED+NL+3*NW : 'TRANC'
*  8+NED+NL+3*NW : 'STRD'
*        IOF0H+1 : 'NUSIGF01' |
*        IOF0H+2 : 'NUSIGF02' | NDEL VALUES
*  ...                        |
*   IOF1H+NDEL+1 : 'CHI'
*   IOF1H+NDEL+2 : 'CHI01' |
*   IOF1H+NDEL+3 : 'CHI02' | NDEL VALUES
*  ...                     |
* IOF1H+2*NDEL+2 : 'CHI--01' |
* IOF1H+2*NDEL+3 : 'CHI--02' | NBESP VALUES
*  ...                       |
*
      MAXH=9+NBESP+2*NDEL+NED+NL+3*NW
      ALLOCATE(IGAR(NGROUP,3,NL),IHNISO(3,NBISO*NMERGE),
     1 ISMIX(NBISO*NMERGE),ISTYP(NBISO*NMERGE),ISTOD(NBISO*NMERGE),
     2 ITYPRO(NL),JPIFI(NDFI),MILVO(NMERGE),ITYPS(NBISO))
      ALLOCATE(MASK(NBISO))
      ALLOCATE(GAR(NGROUP,MAXH),WGAR(NGROUP**2,NL),XSECT(0:NBMIX),
     1 WSTRD(NGCOND),SDEN(NBISO*NMERGE),VOLISO(NBISO*NMERGE),
     2 TNISO(NBISO*NMERGE),TMPXS(NGCOND),WDLA(NDEL),WORK(NGROUP))
      ALLOCATE(WSCAT(NGCOND,NGCOND,NL),GAS(NGCOND,MAXH))
      ALLOCATE(HMAKE(MAXH+NL))
*----
*  FOR AVERAGED NEUTRON VELOCITY
*  V=SQRT(2*ENER/M)=SQRT(2/M)*SQRT(ENER)
*  SQFMAS=SQRT(2/M) IN CM/S/SQRT(EV) FOR V IN CM/S AND E IN EV
*        =SQRT(2*1.602189E-19(J/EV)* 1.0E4(CM2/M2) /1.67495E-27 (KG))
*        =1383155.30602 CM/S/SQRT(EV)
*----
      EVJ=XDRCST('eV','J')
      NMASS=XDRCST('Neutron mass','kg')
      SQFMAS=SQRT(2.0D4*EVJ/NMASS)
*
      CALL EDIMAX(NBISO,ISONAM,MIX,IPRINT,NREGIO,NMERGE,MATCOD,IMERGE,
     1 LSISO,LISO,MAXISO)
      IF(MAXISO.GT.NBISO*NMERGE) CALL XABORT('EDIMIC: MAXISO OVERFLOW.')
      IF(CURNAM.NE.' ') THEN
        CALL LCMSIX(IPEDIT,CURNAM,1)
        IF(MAXISO.GT.0) JPEDIT=LCMLID(IPEDIT,'ISOTOPESLIST',MAXISO)
      ENDIF
*
      DO 10 ISO=1,NBISO*NMERGE
      SDEN(ISO)=0.0
      VOLISO(ISO)=0.0
   10 CONTINUE
      IOF0H=8+NED+NL+3*NW
      IOF1H=8+NED+NL+3*NW+NDEL
      IOF2H=8+NED+NL+3*NW+2*NDEL
      JJISO=0
      CONV=1.0E6*XDRCST('eV','J')
      DO 440 INM=1,NMERGE
*----
*  PRELIMINARY CALCULATIONS FOR STRD CROSS SECTIONS
*----
      LSTRD=ILEAKS.GE.1
      IF(LSTRD) THEN
         IF(ILEAKS.EQ.1) THEN
            CALL LCMGET(IPFLUX,'DIFFB1HOM',GAR(1,1))
         ELSE IF(ILEAKS.EQ.3) THEN
            CALL LCMGET(IPFLUX,'B2  HETE',B2T)
            B2ALL=B2T(1)+B2T(2)+B2T(3)
            IF(B2ALL.EQ.0.0) THEN
               B2T(1)=1.0/3.0
               B2T(2)=B2T(1)
               B2T(3)=B2T(1)
            ELSE
               B2T(1)=B2T(1)/B2ALL
               B2T(2)=B2T(2)/B2ALL
               B2T(3)=B2T(3)/B2ALL
            ENDIF
         ENDIF
         IGRFIN=0
         XSECT(0)=0.0
         DO 50 IGRCND=1,NGCOND
         ZNU=0.0D0
         ZDEN=0.0D0
         ZFL1=0.0D0
         ZFL2=0.0D0
         IGRDEB=IGRFIN+1
         IGRFIN=IGCOND(IGRCND)
         CALL LCMSIX(IPLIB,'MACROLIB',1)
         JPLIB=LCMGID(IPLIB,'GROUP')
         JPFLUX=LCMGID(IPFLUX,'FLUX')
         DO 40 IGR=IGRDEB,IGRFIN
         KPLIB=LCMGIL(JPLIB,IGR)
         CALL LCMGET(KPLIB,'NTOT0',XSECT(1))
         IF((ILEAKS.EQ.2).OR.(ILEAKS.EQ.3)) THEN
            CALL LCMLEL(JPFLUX,IGR,ILCMLN,ITYLCM)
            IF(ILCMLN.EQ.0) CALL XABORT('EDIMIC: MISSING FLUX INFO.')
            ALLOCATE(WORKF(ILCMLN))
            CALL LCMGDL(JPFLUX,IGR,WORKF)
         ENDIF
         FL1=0.0
         FL2=0.0
         DO 20 IREGIO=1,NREGIO
         MATNUM=MATCOD(IREGIO)
         IF(IMERGE(IREGIO).EQ.INM) THEN
            VOLREG=VOLUME(IREGIO)
            IF(IADJ.EQ.0) THEN
              FL1=FLUXES(IREGIO,IGR,1)
              IF(NW.GE.1) FL2=FLUXES(IREGIO,IGR,2)
            ELSE IF(IADJ.EQ.1) THEN
              IF(ILEAKS.NE.1) CALL XABORT('EDIMIC: DIRECT-ADJOINT WEIG'
     1        //'HTING NOT IMPLEMENTED.')
              FL1=FLUXES(IREGIO,IGR,1)*AFLUXE(IREGIO,IGR,1)
              IF(NW.GE.1) FL2=FLUXES(IREGIO,IGR,2)*
     1        AFLUXE(IREGIO,IGR,2)
            ENDIF
            IF(NW.EQ.0) THEN
               ZLEAK=0.0
               IF(ILEAKS.EQ.1) THEN
                  ZLEAK=GAR(IGR,1)*FLUXES(IREGIO,IGR,1)
               ELSE IF(ILEAKS.EQ.2) THEN
                  ZLEAK=WORKF(KEYFLX(IREGIO)+ILCMLN/2)
               ELSE IF(ILEAKS.EQ.3) THEN
                  ZLEAK=B2T(1)*WORKF(KEYFLX(IREGIO)+NREGIO)+
     1                  B2T(2)*WORKF(KEYFLX(IREGIO)+2*NREGIO)+
     2                  B2T(3)*WORKF(KEYFLX(IREGIO)+3*NREGIO)
               ENDIF
               ZNU=ZNU+ZLEAK*VOLREG
               ZDEN=ZDEN+XSECT(MATNUM)*FL1*VOLREG
               ZFL1=ZFL1+FL1*VOLREG
               ZFL2=ZFL2+FL1*VOLREG
            ELSE
               ZNU=ZNU+FL2*VOLREG
               ZDEN=ZDEN+XSECT(MATNUM)*FL2*VOLREG
               ZFL1=ZFL1+FL1*VOLREG
               ZFL2=ZFL2+FL2*VOLREG
            ENDIF
         ENDIF
   20    CONTINUE
         IF((ILEAKS.EQ.2).OR.(ILEAKS.EQ.3)) DEALLOCATE(WORKF)
         CALL LCMLEN(KPLIB,'SIGS01',LENGTH,ITYLCM)
         IF((LENGTH.EQ.NBMIX).AND.(NL.GE.2)) THEN
            CALL LCMGET(KPLIB,'SIGS01',XSECT(1))
            DO 30 IREGIO=1,NREGIO
            MATNUM=MATCOD(IREGIO)
            IF(IMERGE(IREGIO).EQ.INM) THEN
               VOLREG=VOLUME(IREGIO)
               IF(IADJ.EQ.0) THEN
                 FL1=FLUXES(IREGIO,IGR,1)
                 IF(NW.GE.1) FL2=FLUXES(IREGIO,IGR,2)
               ELSE IF(IADJ.EQ.1) THEN
                 FL1=FLUXES(IREGIO,IGR,1)*AFLUXE(IREGIO,IGR,1)
                 IF(NW.GE.1) FL2=FLUXES(IREGIO,IGR,2)*
     1           AFLUXE(IREGIO,IGR,2)
               ENDIF
               IF(NW.EQ.0) THEN
                  ZDEN=ZDEN-XSECT(MATNUM)*FL1*VOLREG
               ELSE
                  ZDEN=ZDEN-XSECT(MATNUM)*FL2*VOLREG
               ENDIF
            ENDIF
   30       CONTINUE
         ENDIF
   40    CONTINUE
         CALL LCMSIX(IPLIB,' ',2)
         WSTRD(IGRCND)=REAL((ZFL1/(3.0*ZNU))*ZFL2/ZDEN)
   50    CONTINUE
      ENDIF
*
      VOLMER=0.0D0
      DO 60 IREGIO=1,NREGIO
      IF(IMERGE(IREGIO).EQ.INM) VOLMER=VOLMER+VOLUME(IREGIO)
   60 CONTINUE
      CALL XDISET(JPIFI,NDFI,0)
      CALL XDLSET(MASK,NBISO,.FALSE.)
      LFIS=.FALSE.
      DO 420 ISO=1,NBISO
      ITYPS(ISO)=ITYPE(ISO)
      IF(MASK(ISO).OR.(LSISO(ISO).EQ.0)) GO TO 420
      DO 90 IREGIO=1,NREGIO
      IF((IMERGE(IREGIO).EQ.INM).AND.(MATCOD(IREGIO).EQ.MIX(ISO)))
     1 GO TO 100
   90 CONTINUE
      GO TO 420
  100 LOGIC=.FALSE.
      DDEN=0.0D0
      DDENZ=0.0D0
*----
*  MERGE/CONDENSE REACTIONS 'NWT0','NWT1','NWAT0','NWAT1','SIGS'//CM,
*   'SCAT'//CM, 'NTOT0', 'NUSIGF', 'CHI', 'CHIxx', 'STRD' AND HVECT
*----
      DO 110 J=1,MAXH+NL
      HMAKE(J)=' '
  110 CONTINUE
      DO 121 J=1,MAXH
      DO 120 I=1,NGCOND
      GAS(I,J)=0.0D0
  120 CONTINUE
  121 CONTINUE
      DO 132 K=1,NL
      DO 131 J=1,NGCOND
      DO 130 I=1,NGCOND
      WSCAT(I,J,K)=0.0D0
  130 CONTINUE
  131 CONTINUE
  132 CONTINUE
      DO 140 I=1,NDEL
      WDLA(I)=0.0
  140 CONTINUE
*----
*  RECOVER THE RADIOACTIVE DECAY CONSTANTS OF DELAYED NEUTRON
*  GROUPS FROM THE MACROLIB IF THEY EXIST
*----
      LWD=.FALSE.
      IF(CURNAM.NE.' ') THEN
         CALL LCMSIX(IPEDIT,'MACROLIB',1)
         CALL LCMLEN(IPEDIT,'LAMBDA-D',ILONG,ITYLCM)
         LWD=(ILONG.EQ.NDEL).AND.(NDEL.GT.0)
         IF(LWD) CALL LCMGET(IPEDIT,'LAMBDA-D',WDLA)
         CALL LCMSIX(IPEDIT,' ',2)
      ENDIF
*
      HMAKE(1)='NWT0'
      LAWR=.FALSE.
      LDECA=.FALSE.
      LMEVF=.FALSE.
      LMEVG=.FALSE.
      DO 145 IW=1,MIN(NW+1,10)
        WRITE(HMAKE(IW),'(3HNWT,I1)') IW-1
        IF(IADJ.EQ.1) WRITE(HMAKE(1+NW+IW),'(4HNWAT,I1)') IW-1
  145 CONTINUE
      ALLOCATE(PARM12(NW+1))
      DO 260 IREGIO=1,NREGIO
      MATNUM=MATCOD(IREGIO)
      VOL=VOLUME(IREGIO)
      IF(IMERGE(IREGIO).EQ.INM) THEN
         IGRFIN=0
         DO 152 IGRCND=1,NGCOND
         IGRDEB=IGRFIN+1
         IGRFIN=IGCOND(IGRCND)
         DO 151 IGR=IGRDEB,IGRFIN
         DO 150 IW=1,NW+1
           GAS(IGRCND,IW)=GAS(IGRCND,IW)+DBLE(FLUXES(IREGIO,IGR,IW)*VOL)
           IF(IADJ.EQ.1) GAS(IGRCND,1+NW+IW)=GAS(IGRCND,1+NW+IW)+
     >     DBLE(AFLUXE(IREGIO,IGR,IW)*VOL)
  150    CONTINUE
  151    CONTINUE
  152    CONTINUE
         LONE=.TRUE.
         DO 250 JSO=ISO,NBISO
         IF((ISONAM(1,ISO).EQ.ISONAM(1,JSO)).AND.
     1      (ISONAM(2,ISO).EQ.ISONAM(2,JSO)).AND.
     2      (MATNUM.EQ.MIX(JSO)).AND.(LSISO(JSO).NE.0)) THEN
            IF(LISO) THEN
              IF(ISONAM(3,ISO).EQ.ISONAM(3,JSO)) GOTO 155
              GOTO 250
            ENDIF
  155       LOGIC=.TRUE.
            ITYPS(ISO)=MAX(ITYPS(ISO),ITYPE(JSO))
            DENVOL=MAX(DEN(JSO),1.0E-20)*VOL
            DDEN=DDEN+DENVOL
            DDENZ=DDENZ+DEN(JSO)*VOL
            KPLIB=IPISO(JSO) ! set JSO-th isotope
            IF(LONE) THEN
               CALL LCMLEN(KPLIB,'AWR',LENGTH,ITYLCM)
               LAWR=(LENGTH.EQ.1)
               IF(LAWR) CALL LCMGET(KPLIB,'AWR',AWR)
               CALL LCMLEN(KPLIB,'MEVF',LENGTH,ITYLCM)
               IF(LENGTH.EQ.1) CALL LCMGET(KPLIB,'MEVF',EMEVF(ISO))
               LMEVF=(LENGTH.EQ.1).OR.(EMEVF(ISO).GT.0.0)
               CALL LCMLEN(KPLIB,'MEVG',LENGTH,ITYLCM)
               IF(LENGTH.EQ.1) CALL LCMGET(KPLIB,'MEVG',EMEVG(ISO))
               LMEVG=(LENGTH.EQ.1).OR.(EMEVG(ISO).GT.0.0)
               CALL LCMLEN(KPLIB,'DECAY',LENGTH,ITYLCM)
               IF(LENGTH.EQ.1) CALL LCMGET(KPLIB,'DECAY',DECAY(ISO))
               LDECA=(LENGTH.EQ.1).OR.(DECAY(ISO).GT.0.0)
               LONE=.FALSE.
            ENDIF
            DO 170 IL=0,NL-1
            WRITE (CM,'(I2.2)') IL
            CALL LCMLEN(KPLIB,'SIGS'//CM,LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
               CALL LCMGET(KPLIB,'SIGS'//CM,GAR(1,4+3*NW+IL))
               HMAKE(4+3*NW+IL)='SIGS'//CM
            ENDIF
            CALL LCMLEN(KPLIB,'NJJS'//CM,LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
               CALL LCMGET(KPLIB,'NJJS'//CM,IGAR(1,1,1+IL))
               CALL LCMGET(KPLIB,'IJJS'//CM,IGAR(1,2,1+IL))
               CALL LCMGET(KPLIB,'SCAT'//CM,WGAR(1,1+IL))
               HMAKE(MAXH+1+IL)='SCAT'//CM
               IPO=0
               DO 160 IGR=1,NGROUP
               IGAR(IGR,3,1+IL)=IPO+1
               IPO=IPO+IGAR(IGR,1,1+IL)
  160          CONTINUE
            ENDIF
  170       CONTINUE
            DO IW=0,MIN(NW,9)
              WRITE(HMAKE(3+2*NW+IW),'(4HNTOT,I1)') IW
              CALL LCMLEN(KPLIB,HMAKE(3+2*NW+IW),ILONG,ITYLCM)
              IF(ILONG.NE.0) THEN
                CALL LCMGET(KPLIB,HMAKE(3+2*NW+IW),GAR(1,3+2*NW+IW))
              ELSE
                CALL LCMGET(KPLIB,'NTOT0',GAR(1,3+2*NW+IW))
              ENDIF
            ENDDO
            CALL LCMLEN(KPLIB,'NUSIGF',LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
               CALL LCMGET(KPLIB,'NUSIGF',GAR(1,4+NL+3*NW))
               HMAKE(4+NL+3*NW)='NUSIGF'
            ENDIF
            CALL LCMLEN(KPLIB,'CHI',LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
               CALL LCMGET(KPLIB,'CHI',GAR(1,1+IOF1H))
               HMAKE(1+IOF1H)='CHI'
            ENDIF
            IF(NDEL.GT.0) THEN
               WRITE(TEXT8,'(6HNUSIGF,I2.2)') NDEL
               CALL LCMLEN(KPLIB,TEXT8,LENGTH,ITYLCM)
               IF(LENGTH.EQ.NGROUP) THEN
                  DO 180 IDEL=1,NDEL
                  WRITE(TEXT8,'(6HNUSIGF,I2.2)') IDEL
                  CALL LCMGET(KPLIB,TEXT8,GAR(1,IOF0H+IDEL))
                  HMAKE(IOF0H+IDEL)=TEXT8
  180             CONTINUE
               ENDIF
               WRITE(TEXT8,'(3HCHI,I2.2)') NDEL
               CALL LCMLEN(KPLIB,TEXT8,LENGTH,ITYLCM)
               IF(LENGTH.EQ.NGROUP) THEN
                  DO 184 IDEL=1,NDEL
                  WRITE(TEXT8,'(3HCHI,I2.2)') IDEL
                  CALL LCMGET(KPLIB,TEXT8,GAR(1,1+IOF1H+IDEL))
                  HMAKE(1+IOF1H+IDEL)=TEXT8
  184             CONTINUE
               ENDIF
            ENDIF
            DO 185 ISP=1,NBESP
            WRITE(TEXT8,'(5HCHI--,I2.2)') ISP
            CALL LCMLEN(KPLIB,TEXT8,LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
               CALL LCMGET(KPLIB,TEXT8,GAR(1,1+IOF2H+ISP))
               HMAKE(1+IOF2H+ISP)=TEXT8
            ENDIF
  185       CONTINUE
            IF(ITRANC.NE.0) THEN
               CALL LCMGET(KPLIB,'TRANC',GAR(1,7+NED+NL+3*NW))
               HMAKE(7+NED+NL+3*NW)='TRANC'
            ENDIF
            DO 186 IGR=1,NGROUP
            GAR(IGR,5+NED+NL+3*NW)=0.0
  186       CONTINUE
            IF(LMEVF) THEN
               CALL LCMGET(KPLIB,'NFTOT',WORK)
               HMAKE(5+NED+NL+3*NW)='H-FACTOR'
               DO 190 IGR=1,NGROUP
               GAR(IGR,5+NED+NL+3*NW)=GAR(IGR,5+NED+NL+3*NW)+WORK(IGR)*
     1         EMEVF(ISO)*REAL(CONV)
  190          CONTINUE
            ENDIF
            IF(LMEVG) THEN
               CALL LCMGET(KPLIB,'NG',WORK)
               HMAKE(5+NED+NL+3*NW)='H-FACTOR'
               DO 195 IGR=1,NGROUP
               GAR(IGR,5+NED+NL+3*NW)=GAR(IGR,5+NED+NL+3*NW)+WORK(IGR)*
     1         EMEVG(ISO)*REAL(CONV)
  195          CONTINUE
            ENDIF
            DO 200 IED=1,NED
            CALL LCMLEN(KPLIB,HVECT(IED),LENGTH,ITYLCM)
            IF((LENGTH.GT.0).AND.(HVECT(IED).NE.'TRANC')) THEN
               CALL LCMGET(KPLIB,HVECT(IED),GAR(1,4+NL+3*NW+IED))
               HMAKE(4+NL+3*NW+IED)=HVECT(IED)
            ENDIF
  200       CONTINUE
            CALL LCMLEN(KPLIB,'OVERV',LENGTH,ITYLCM)
            IF(LENGTH.GT.0) THEN
               CALL LCMGET(KPLIB,'OVERV',GAR(1,6+NED+NL+3*NW))
            ELSE
               ALLOCATE(ENR(NGROUP+1))
               CALL LCMGET(IPLIB,'ENERGY',ENR)
               IF(ENR(NGROUP+1).EQ.0.0) ENR(NGROUP+1)=1.0E-5
               DO 205 IGR=1,NGROUP
               ENEAVG=SQRT(ENR(IGR)*ENR(IGR+1))
               GAR(IGR,6+NED+NL+3*NW)=1.0/(REAL(SQFMAS)*SQRT(ENEAVG))
  205          CONTINUE
               DEALLOCATE(ENR)
            ENDIF
            HMAKE(6+NED+NL+3*NW)='OVERV'
*
            IGRFIN=0
            DO 242 IGRCND=1,NGCOND
            IGRDEB=IGRFIN+1
            IGRFIN=IGCOND(IGRCND)
            DO 241 IGR=IGRDEB,IGRFIN
            PARM0=FLUXES(IREGIO,IGR,1)*DENVOL
            PARM3=0.0D0
            PARM4=0.0D0
            CALL XDDSET(PARM12,NW+1,0.0D0)
            IF(IADJ.EQ.0) THEN
              DO 206 IW=1,NW+1
                PARM12(IW)=FLUXES(IREGIO,IGR,IW)*DENVOL
  206         CONTINUE
              PARM3=0.0D0
              DO 210 JREGIO=1,NREGIO
              IF(IMERGE(JREGIO).EQ.INM) THEN
                PARM3=PARM3+FLUXES(JREGIO,IGR,1)*VOLUME(JREGIO)
              ENDIF
  210         CONTINUE
              PARM3=DENVOL*PARM3/VOLMER
              PARM4=DENVOL
            ELSE IF(IADJ.EQ.1) THEN
              DO 211 IW=1,NW+1
                PARM12(IW)=FLUXES(IREGIO,IGR,IW)*AFLUXE(IREGIO,IGR,IW)*
     >          DENVOL
  211         CONTINUE
              PARM3=0.0D0
              DO 212 JREGIO=1,NREGIO
              IF(IMERGE(JREGIO).EQ.INM) THEN
                PARM3=PARM3+FLUXES(JREGIO,IGR,1)*AFLUXE(JREGIO,IGR,1)*
     >          VOLUME(JREGIO)
              ENDIF
  212         CONTINUE
              PARM3=DENVOL*PARM3/VOLMER
              PARM4=AFLUXE(IREGIO,IGR,1)*DENVOL
            ENDIF
            DO 215 J=3+2*NW,MAXH
            IF(HMAKE(J).NE.' ') THEN
              IF(J.EQ.6+NED+NL+3*NW) THEN
                GAS(IGRCND,J)=GAS(IGRCND,J)+DBLE(GAR(IGR,J))*PARM3 ! OVERV
              ELSE IF((J.EQ.4+NL+3*NW).OR.
     >          ((J.GE.1+IOF0H).AND.(J.LE.NDEL+IOF0H))) THEN
                GAS(IGRCND,J)=GAS(IGRCND,J)+DBLE(GAR(IGR,J))*PARM0 ! nu*fission cross sections
              ELSE IF((J.GE.1+IOF1H).AND.(J.LE.MAXH)) THEN
                GAS(IGRCND,J)=GAS(IGRCND,J)+DBLE(GAR(IGR,J))*PARM4 ! fission spectrum
              ELSE IF((J.GE.4+2*NW).AND.(J.LE.3+3*NW)) THEN
                IW=J-2-2*NW                      ! NTOT1 cross sections
                GAS(IGRCND,J)=GAS(IGRCND,J)+DBLE(GAR(IGR,J))*PARM12(IW)
              ELSE IF((J.GE.5+3*NW).AND.(J.LE.3+NL+3*NW)) THEN
                IW=MIN(J-3-3*NW,NW+1)            ! SOGS01 cross sections
                GAS(IGRCND,J)=GAS(IGRCND,J)+DBLE(GAR(IGR,J))*PARM12(IW)
              ELSE IF(J.EQ.8+NED+NL+3*NW) THEN
                GO TO 215 ! STRD case
              ELSE IF(J.LE.IOF1H) THEN
                GAS(IGRCND,J)=GAS(IGRCND,J)+DBLE(GAR(IGR,J))*PARM12(1) ! P0 cross sections
              ENDIF
            ENDIF
  215       CONTINUE
            DO 240 IL=0,NL-1
            IF(HMAKE(MAXH+1+IL).NE.' ') THEN
*              IGRCND IS THE SECONDARY GROUP.
               IW=MIN(IL,NW)+1
               NGSCAT=IGAR(IGR,1,1+IL)
               IGSCAT=IGAR(IGR,2,1+IL)
               JGRFIN=0
               DO 230 JGRCND=1,NGCOND
               JGRDEB=JGRFIN+1
               JGRFIN=IGCOND(JGRCND)
               J2=MIN(JGRFIN,IGSCAT)
               J1=MAX(JGRDEB,IGSCAT-NGSCAT+1)
               TMP=0.0D0
               IPO=IGAR(IGR,3,1+IL)+IGSCAT-J2
               DO 220 JGR=J2,J1,-1
               IF(IADJ.EQ.0) THEN
                 TMP=TMP+WGAR(IPO,1+IL)*FLUXES(IREGIO,JGR,IW)*DENVOL
               ELSE IF(IADJ.EQ.1) THEN
                  TMP=TMP+WGAR(IPO,1+IL)*AFLUXE(IREGIO,IGR,IW)*
     >           FLUXES(IREGIO,JGR,IW)*DENVOL
               ENDIF
               IPO=IPO+1
  220          CONTINUE
               WSCAT(IGRCND,JGRCND,1+IL)=WSCAT(IGRCND,JGRCND,1+IL)+TMP
  230          CONTINUE
            ENDIF
  240       CONTINUE
  241       CONTINUE
  242       CONTINUE
            MASK(JSO)=.TRUE.
            GO TO 250
         ENDIF
  250    CONTINUE
      ENDIF
  260 CONTINUE
      DEALLOCATE(PARM12)
      IF(LOGIC) THEN
         JJISO=JJISO+1
         IF(JJISO.GT.MAXISO) CALL XABORT('EDIMIC: INSUFFICIENT ALLOCAT'
     1   //'ED SPACE FOR ISMIX, ISTYP, SDEN, VOLISO AND IHNISO.')
         IF(LISO) THEN
           WRITE(HNEW,'(3A4)') (ISONAM(I0,ISO),I0=1,3)
         ELSE
           WRITE(HNEW,'(2A4,I4.4)') (ISONAM(I0,ISO),I0=1,2),INM
         ENDIF
         READ(HNEW,'(3A4)') (IHNISO(I0,JJISO),I0=1,3)
         ISMIX(JJISO)=INM
         ISTYP(JJISO)=ISO
         SDEN(JJISO)=REAL(DDENZ/VOLMER)
         VOLISO(JJISO)=REAL(VOLMER)
         TNISO(JJISO)=TN(ISO)
         DO 270 I=1,NDFI
         IF(IPIFI(I,INM).EQ.0) GO TO 270
         IF((ISONAM(1,IPIFI(I,INM)).EQ.ISONAM(1,ISO)).AND.
     1      (ISONAM(2,IPIFI(I,INM)).EQ.ISONAM(2,ISO))) JPIFI(I)=JJISO
  270    CONTINUE
         LFIS=LFIS.OR.(ITYPE(ISO).EQ.2)
         IF(IPRINT.GT.1) THEN
            WRITE (6,600) HNEW,JJISO
            WRITE(6,'(/17H NUMBER DENSITY =,1P,E12.4)') DDEN/VOLMER
         ENDIF
*
*        UP-SCATTERING CORRECTIONS.
         IF(ILUPS.EQ.1) THEN
            DO 282 JGR=2,NGCOND
            DO 281 IGR=1,JGR-1
            DO 280 IL=0,NL-1
            WSCAT(JGR,IGR,1+IL)=WSCAT(JGR,IGR,1+IL)-WSCAT(IGR,JGR,1+IL)
            WSCAT(IGR,JGR,1+IL)=0.0D0
  280       CONTINUE
  281       CONTINUE
  282       CONTINUE
         ENDIF
*
         ALLOCATE(PHIAV(NW+1),AHIAV(NW+1))
         DO 360 IGRCND=1,NGCOND
*
*        DIVIDE MATRIX XS BY INTEGRATED FLUX
         DO 341 IL=0,NL-1
         IW=MIN(IL,NW)+1
         PHIAV(IW)=GAS(IGRCND,IW)/VOLMER
         TMP=GAS(IGRCND,4+3*NW+IL)
         DO 330 JGRCND=1,NGCOND
         IF(JGRCND.NE.IGRCND) TMP=TMP-WSCAT(JGRCND,IGRCND,1+IL)
  330    CONTINUE
         QEN=REAL(MAX(ABS(TMP),ABS(WSCAT(IGRCND,IGRCND,1+IL))))
         IF((QEN.GT.0.0).AND.(IADJ.EQ.0)) THEN
            ERR=ABS(REAL(TMP-WSCAT(IGRCND,IGRCND,1+IL)))/QEN
            IF(ERR.GT.1.0E-3) WRITE(6,620) IGRCND,IL,100.0*ERR,HNEW
            WSCAT(IGRCND,IGRCND,1+IL)=TMP
         ENDIF
         DO 340 JGRCND=1,NGCOND
         AHIAV(IW)=1.0D0
         IF(IADJ.EQ.1) AHIAV(IW)=GAS(JGRCND,1+NW+IW)/VOLMER
         IF(PHIAV(IW).GT.0.0D0) THEN
            WSCAT(JGRCND,IGRCND,1+IL)=WSCAT(JGRCND,IGRCND,1+IL)
     1      /(DDEN*AHIAV(IW)*PHIAV(IW))
         ELSE
            WSCAT(JGRCND,IGRCND,1+IL)=0.0D0
         ENDIF
  340    CONTINUE
  341    CONTINUE
*
*        DIVIDE VECTORIAL XS BY INTEGRATED FLUX
         DO 345 IW=1,NW+1
         PHIAV(IW)=GAS(IGRCND,IW)/VOLMER
         AHIAV(IW)=1.0
         IF(IADJ.EQ.1) AHIAV(IW)=GAS(IGRCND,1+NW+IW)/VOLMER
  345    CONTINUE
         DO 350 J=3+2*NW,MAXH
         IF((J.EQ.4+NL+3*NW).OR.
     >      ((J.GE.1+IOF0H).AND.(J.LE.NDEL+IOF0H))) THEN
            IF(PHIAV(1).GT.0.0D0) THEN
               GAS(IGRCND,J)=GAS(IGRCND,J)/(DDEN*PHIAV(1)) ! nu*fission cross sections
            ELSE
               GAS(IGRCND,J)=0.0D0
            ENDIF
         ELSE IF((J.GE.1+IOF1H).AND.(J.LE.MAXH)) THEN
            GAS(IGRCND,J)=GAS(IGRCND,J)/(DDEN*AHIAV(1)) ! fission spectrum
         ELSE IF((J.GE.4+2*NW).AND.(J.LE.3+3*NW)) THEN
            IW=J-2-2*NW
            IF(PHIAV(IW).NE.0.0) THEN
              GAS(IGRCND,J)=GAS(IGRCND,J)/(DDEN*AHIAV(IW)*PHIAV(IW)) ! NTOT1 cross sections
            ELSE
              GAS(IGRCND,J)=0.0D0
            ENDIF
         ELSE IF((J.GE.5+3*NW).AND.(J.LE.3+NL+3*NW)) THEN
            IW=MIN(J-3-3*NW,NW+1)
            IF(PHIAV(IW).NE.0.0) THEN
              GAS(IGRCND,J)=GAS(IGRCND,J)/(DDEN*AHIAV(IW)*PHIAV(IW)) ! SIGS01 cross sections
            ELSE
              GAS(IGRCND,J)=0.0D0
            ENDIF
         ELSE IF(J.EQ.8+NED+NL+3*NW) THEN
            GO TO 350 ! STRD case
         ELSE IF(PHIAV(1).GT.0.0D0) THEN
            GAS(IGRCND,J)=GAS(IGRCND,J)/(DDEN*AHIAV(1)*PHIAV(1)) ! P0 cross sections
         ELSE
            GAS(IGRCND,J)=0.0D0
         ENDIF
  350    CONTINUE
*
         IF(LSTRD) THEN
            J=8+NED+NL+3*NW
            HMAKE(J)='STRD'
            IF(NW.GE.1) THEN
               GAS(IGRCND,J)=GAS(IGRCND,4+2*NW)
            ELSE
               GAS(IGRCND,J)=GAS(IGRCND,3+2*NW)
            ENDIF
            IF((HMAKE(5+3*NW).NE.' ').AND.(NL.GE.2)) THEN
               GAS(IGRCND,J)=GAS(IGRCND,J)-GAS(IGRCND,5+3*NW)
            ENDIF
            GAS(IGRCND,J)=GAS(IGRCND,J)*WSTRD(IGRCND)
         ENDIF
  360    CONTINUE
         DEALLOCATE(AHIAV,PHIAV)
*
*        DIVIDE INTEGRATED FLUXES BY VOLUMES
         DO 366 IW=1,NW+1
         DO 365 IGRCND=1,NGCOND
         GAS(IGRCND,IW)=GAS(IGRCND,IW)/VOLMER
         IF(IADJ.EQ.1) GAS(IGRCND,NW+IW)=GAS(IGRCND,NW+IW)/VOLMER
  365    CONTINUE
  366    CONTINUE
*
         IF(CURNAM.NE.' ') THEN
            IF(NOUT.GT.0) THEN
              DO J=1,MAXH+NL
                DO IOUT=1,NOUT
                  IF(HMAKE(J).EQ.HVOUT(IOUT)) GO TO 370
                ENDDO
                HMAKE(J)=' '
  370           CONTINUE
              ENDDO
            ENDIF
            KPEDIT=LCMDIL(JPEDIT,JJISO) ! set JJISO-th isotope
            CALL LCMPTC(KPEDIT,'ALIAS',12,1,HNEW)
            IF(LAWR)  CALL LCMPUT(KPEDIT,'AWR',1,2,AWR)
            IF(LMEVF) CALL LCMPUT(KPEDIT,'MEVF',1,2,EMEVF(ISO))
            IF(LMEVG) CALL LCMPUT(KPEDIT,'MEVG',1,2,EMEVG(ISO))
            IF(LDECA) CALL LCMPUT(KPEDIT,'DECAY',1,2,DECAY(ISO))
            DO 380 J=1,MAXH
            IF(HMAKE(J).NE.' ') THEN
              DO 375 IGCD=1,NGCOND
                TMPXS(IGCD)=REAL(GAS(IGCD,J))
  375         CONTINUE
              CALL LCMPUT(KPEDIT,HMAKE(J),NGCOND,2,TMPXS)
            ENDIF
  380       CONTINUE
            DO 390 IL=1,NL
            ITYPRO(IL)=0
            IF(HMAKE(MAXH+IL).NE.' ') ITYPRO(IL)=1
  390       CONTINUE
            IF(ITYPRO(1).EQ.0) GO TO 405
            ALLOCATE(GA1(NL*NGCOND),GA2(NL*NGCOND*NGCOND))
            IOF1=0
            IOF2=0
            DO 402 IL=1,NL
            DO 401 IG2=1,NGCOND
            IOF1=IOF1+1
            GA1(IOF1)=REAL(GAS(IG2,3+3*NW+IL))
            DO 400 IG1=1,NGCOND
            IOF2=IOF2+1
            GA2(IOF2)=REAL(WSCAT(IG1,IG2,IL))
  400       CONTINUE
  401       CONTINUE
  402       CONTINUE
            CALL XDRLGS(KPEDIT,1,IPRINT,0,NL-1,1,NGCOND,GA1,GA2,ITYPRO)
            DEALLOCATE(GA2,GA1)
  405       IF(NDEL.NE.0) THEN
               IF(HMAKE(IOF0H+1).NE.' ') THEN
                  CALL LCMPUT(KPEDIT,'LAMBDA-D',NDEL,2,WDLA)
               ENDIF
            ENDIF
         ENDIF
         IF(IPRINT.GT.3) THEN
            DO 410 J=1,MAXH
            IF(HMAKE(J).NE.' ') THEN
               WRITE (6,610) HMAKE(J),(GAS(I,J),I=1,NGCOND)
            ENDIF
  410       CONTINUE
            WRITE (6,610) 'SIGA    ',(GAS(I,3+2*NW)-GAS(I,4+3*NW),
     >                    I=1,NGCOND)
            WRITE (6,610) 'SIGW00  ',(WSCAT(I,I,1),I=1,NGCOND)
            IF(NL.GT.1) THEN
               WRITE (6,610) 'SIGW01  ',(WSCAT(I,I,2),I=1,NGCOND)
            ENDIF
            IF(LWD) WRITE (6,610) 'LAMBDA-D',(WDLA(I),I=1,NDEL)
         ENDIF
         IF(IPRINT.GT.4) CALL LCMLIB(KPEDIT)
      ENDIF
  420 CONTINUE
*----
*  SAVE FISSION YIELD DATA
*----
      LFIS=LFIS.AND.(NDFI.GT.0).AND.(JJISO.GT.0)
      IF((CURNAM.NE.' ').AND.LFIS) THEN
         DO 425 I=1,NDFI
         IF((JPIFI(I).EQ.0).AND.(IPIFI(I,INM).NE.0)) THEN
            ISO=IPIFI(I,INM)
            IF(DEN(ISO).EQ.0.0) GO TO 425
            WRITE(HNAMIS,'(3A4)') (ISONAM(I0,ISO),I0=1,3)
            WRITE(HSMG,'(29HEDIMIC: THE FISSILE ISOTOPE '',A8,
     1      34H'' MUST BE SELECTED IN MICR OPTION.)') HNAMIS(:8)
            CALL XABORT(HSMG)
         ENDIF
  425    CONTINUE
         DO 430 JSO=1,JJISO
         IF(ISMIX(JSO).EQ.INM) THEN
            ISO=ISTYP(JSO)
            WRITE(HNEW,'(3A4)') (IHNISO(I0,JSO),I0=1,3)
            KPEDIT=LCMGIL(JPEDIT,JSO) ! set JSO-th isotope
            IF(YIELD(1,ISO,INM).GT.0.0) THEN
               CALL LCMPUT(KPEDIT,'YIELD',NGCOND+1,2,YIELD(1,ISO,INM))
               CALL LCMPUT(KPEDIT,'PYIELD',NDFI,2,PYIELD(1,ISO,INM))
               CALL LCMPUT(KPEDIT,'PIFI',NDFI,1,JPIFI)
            ENDIF
         ENDIF
  430    CONTINUE
      ENDIF
  440 CONTINUE
      IF(CURNAM.NE.' ') CALL LCMSIX(IPEDIT,' ',2)
*
      IF(CURNAM.NE.' ') THEN
        CALL LCMSIX(IPEDIT,CURNAM,1)
        TEXT12='L_LIBRARY'
        CALL LCMPTC(IPEDIT,'SIGNATURE',12,1,TEXT12)
        MAXISM=0
        DO 460 IBM=1,NMERGE
        MAX0=0
        DO 450 ISO=1,JJISO
        IF(ISMIX(ISO).EQ.IBM) MAX0=MAX0+1
  450   CONTINUE
        MAXISM=MAX(MAXISM,MAX0)
  460   CONTINUE
        IF(NED.GT.0) CALL LCMPTC(IPEDIT,'ADDXSNAME-P0',8,NED,HVECT)
        NCOMB=0
        IF(JJISO.GT.0) THEN
           CALL LCMPUT(IPEDIT,'ISOTOPESUSED',3*JJISO,3,IHNISO)
           CALL LCMPUT(IPEDIT,'ISOTOPESMIX',JJISO,1,ISMIX)
           CALL LCMPUT(IPEDIT,'ISOTOPESVOL',JJISO,2,VOLISO)
           CALL LCMPUT(IPEDIT,'ISOTOPESTEMP',JJISO,2,TNISO)
           CALL LCMPUT(IPEDIT,'ISOTOPESDENS',JJISO,2,SDEN)
           DO 500 IISO=1,JJISO
           DO 480 I0=1,3
           IHNISO(I0,IISO)=ISONRF(I0,ISTYP(IISO))
  480      CONTINUE
           ISTOD(IISO)=IEVOL(ISTYP(IISO))
           ISTYP(IISO)=ITYPS(ISTYP(IISO))
           IF((ISTOD(IISO).NE.1).AND.(ISTYP(IISO).GE.1)) THEN
              IBM=ISMIX(IISO)
              IF(IBM.EQ.0) GO TO 500
              DO 490 J=1,NCOMB
              IF(IBM.EQ.MILVO(J)) GO TO 500
  490         CONTINUE
              NCOMB=NCOMB+1
              IF(NCOMB.GT.NMERGE) CALL XABORT('EDIMIC: MILVO OVERFLOW.')
              MILVO(NCOMB)=IBM
           ENDIF
  500      CONTINUE
           CALL LCMPUT(IPEDIT,'ISOTOPERNAME',3*JJISO,3,IHNISO)
           CALL LCMPUT(IPEDIT,'ISOTOPESTODO',JJISO,1,ISTOD)
           CALL LCMPUT(IPEDIT,'ISOTOPESTYPE',JJISO,1,ISTYP)
        ENDIF
        ALLOCATE(VOLM(NMERGE))
        CALL XDRSET(VOLM,NMERGE,0.0)
        DO 510 IREGIO=1,NREGIO
        INM=IMERGE(IREGIO)
        IF(INM.GT.0) VOLM(INM)=VOLM(INM)+VOLUME(IREGIO)
  510   CONTINUE
        CALL LCMPUT(IPEDIT,'MIXTURESVOL',NMERGE,2,VOLM)
        CALL LCMPUT(IPEDIT,'K-EFFECTIVE',1,2,EIGENK)
        CALL LCMPUT(IPEDIT,'K-INFINITY',1,2,EIGINF)
        IF(ILEAKS.GT.0) CALL LCMPUT(IPEDIT,'B2  B1HOM',1,2,B2(4))
        DEALLOCATE(VOLM)
        CALL XDISET(IPAR,NSTATE,0)
        IPAR(1)=NMERGE
        IPAR(2)=JJISO
        IPAR(3)=NGCOND
        IPAR(4)=NL
        IPAR(5)=ITRANC
        IF(ITRANC.NE.0) IPAR(5)=2
        IPAR(7)=1
        IPAR(11)=NDEPL
        IPAR(12)=NCOMB
        IPAR(13)=NED
        IPAR(14)=NMERGE
        IPAR(16)=NBESP
        IPAR(18)=1
        IPAR(19)=NDEL
        IPAR(20)=NDFI
        IPAR(22)=MAXISM
        IPAR(25)=NW
        CALL LCMPUT(IPEDIT,'STATE-VECTOR',NSTATE,1,IPAR)
        IF(IPRINT.GT.3) THEN
           WRITE(6,630) IPRINT,(IPAR(I),I=1,13)
           WRITE(6,640) (IPAR(I),I=14,25)
        ENDIF
        ALLOCATE(ENR(NGROUP+1))
        CALL LCMGET(IPLIB,'ENERGY',ENR)
        DO 520 IGRCND=1,NGCOND
        ENR(IGRCND+1)=ENR(IGCOND(IGRCND)+1)
  520   CONTINUE
        IF(ENR(NGCOND+1).EQ.0.0) ENR(NGCOND+1)=1.0E-5
        CALL LCMPUT(IPEDIT,'ENERGY',NGCOND+1,2,ENR)
        IF(NBESP.GT.0) THEN
          IF(NBESP.GT.MAXESP) CALL XABORT('EDIMIC: MAXESP OVERFLOW.')
          CALL LCMGET(IPLIB,'CHI-ENERGY',EESP)
          CALL XDISET(IESP,NBESP+1,0)
          IIG=0
          DO IG=1,NGCOND+1
            IF(IIG.GT.NBESP) CALL XABORT('EDIMIC: BAD LIMITS FOR ENERG'
     1      //'Y-DEPENDENT FISSION SPECTRA(1).')
            IF(EESP(IIG+1).GE.0.999*ENR(IG)) THEN
              IIG=IIG+1
              IESP(IIG)=IG-1
            ENDIF
          ENDDO
          IF(IIG.NE.NBESP+1) CALL XABORT('EDIMIC: BAD LIMITS FOR ENERG'
     1    //'Y-DEPENDENT FISSION SPECTRA(2).')
          CALL LCMPUT(IPEDIT,'CHI-ENERGY',NBESP+1,2,EESP)
          CALL LCMPUT(IPEDIT,'CHI-LIMITS',NBESP+1,1,IESP)
        ENDIF
        DO 530 IGRCND=1,NGCOND
        ENR(IGRCND)=LOG(ENR(IGRCND)/ENR(IGRCND+1))
  530   CONTINUE
        CALL LCMPUT(IPEDIT,'DELTAU',NGCOND,2,ENR)
        DEALLOCATE(ENR)
*
        CALL LCMSIX(IPEDIT,' ',2)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(HMAKE)
      DEALLOCATE(GAS,WSCAT)
      DEALLOCATE(WORK,WDLA,TMPXS,TNISO,VOLISO,SDEN,WSTRD,XSECT,WGAR,GAR)
      DEALLOCATE(MASK)
      DEALLOCATE(ITYPS,MILVO,JPIFI,ITYPRO,ISTOD,ISTYP,ISMIX,IHNISO,IGAR)
      RETURN
*
  600 FORMAT (//44H CROSS SECTION OF MERGED/CONDENSED ISOTOPE ',A12,
     1 7H' (ISO=,I8,2H):)
  610 FORMAT (/11H REACTION ',A12,2H':/(1X,1P,10E12.4))
  620 FORMAT(/52H *** WARNING *** NORMALIZATION OF THE WITHIN-GROUP S,
     1 27HCATTERING TRANSFER IN GROUP,I4,10H AND ORDER,I3,3H BY,F6.2,
     2 11H% ISOTOPE=',A12,2H'.)
  630 FORMAT(/8H OPTIONS/8H -------/
     1 7H IPRINT,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H MAXMIX,I6,31H   (MAXIMUM NUMBER OF MIXTURES)/
     3 7H NBISO ,I6,36H   (NUMBER OF ISOTOPES OR MATERIALS)/
     4 7H NGRP  ,I6,28H   (NUMBER OF ENERGY GROUPS)/
     5 7H NL    ,I6,30H   (NUMBER OF LEGENDRE ORDERS)/
     6 7H ITRANC,I6,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     7 57H=RECOVER FROM LIBRARY/3=WIMS-D TYPE/4=LEAKAGE CORRECTION)/
     8 7H IPROB ,I6,23H   (0=DIRECT/1=ADJOINT)/
     9 7H ITIME ,I6,28H   (1=STEADY-STATE/2=PROMPT)/
     1 7H NLIB  ,I6,32H   (NUMBER OF SETS OF LIBRARIES)/
     2 7H NGF   ,I6,48H   (NUMBER OF FAST GROUP WITHOUT SELF-SHIELDING)/
     3 7H IGRMAX,I6,41H   (LAST GROUP INDEX WITH SELF-SHIELDING)/
     4 7H NDEPL ,I6,33H   (NUMBER OF DEPLETING ISOTOPES)/
     5 7H NCOMB ,I6,33H   (NUMBER OF DEPLETING MIXTURES)/
     6 7H NEDMAC,I6,34H   (NUMBER OF CROSS SECTION EDITS))
  640 FORMAT(7H NBMIX ,I6,23H   (NUMBER OF MIXTURES)/
     1 7H NRES  ,I6,40H   (NUMBER OF SETS OF RESONANT MIXTURES)/
     2 7H NBESP ,I6,47H   (NUMBER OF ENERGY-DEPENDENT FISSION SPECTRA)/
     3 7H IPROC ,I6,48H   (-1=SKIP LIBRARY PROCESSING/0=DILUTION INTERP,
     4 48HOLATION/1=USE PHYSICAL TABLES/2=BUILD A DRAGLIB/,
     5 55H3=COMPUTE CALENDF TABLES/4=COMPUTE SLOWING-DOWN TABLES)/
     6 7H IMAC  ,I6,45H   (0=DO NOT/1=DO BUILD AN EMBEDDED MACROLIB)/
     7 7H NDEL  ,I6,31H   (NUMBER OF PRECURSOR GROUPS)/
     8 7H NFISS ,I6,31H   (NUMBER OF FISSILE ISOTOPES)/
     9 7H ISOADD,I6,37H   (0=COMPLETE BURNUP CHAIN/1=DO NOT)/
     1 7H MAXISM,I6,40H   (MAX. NUMBER OF ISOTOPES PER MIXTURE)/
     2 7H IPRECI,I6,34H   (CALENDF ACCURACY FLAG:1/2/3/4)/
     3 7H IADF  ,I6,19H   (ADF FLAG:0/1/2)/
     4 7H NW    ,I6,47H   (=0: FLUX WEIGHTING FOR P1 INFO; =1: CURRENT,
     5 23H WEIGHTING FOR P1 INFO))
      END
