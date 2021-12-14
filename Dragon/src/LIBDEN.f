*DECK LIBDEN
      SUBROUTINE LIBDEN (IPLIB,NGROUP,NBISO,NBMIX,NL,NDEL,NESP,ISONAM,
     1 IPISO,MIX,DEN,MASK,MASKL,NED,NAMEAD,ITRANC,MAXNFI,NPART,LSAME,
     2 ITSTMP,TMPDAY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transformation of the isotope ordered microscopic cross sections to
* group ordered macroscopic cross sections (part 2).
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NGROUP  number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NBMIX   number of mixtures present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NDEL    number of delayed precursor groups.
* NESP    number of energy-dependent fission spectra.
* ISONAM  names of microlib isotopes.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture number of each isotope (can be zero).
* DEN     density of each isotope.
* MASK    mixture mask (=.true. if a mixture is to be made).
* MASKL   group mask (=.true. if an energy group is to be treated).
* NED     number of extra edit vectors.
* NAMEAD  names of these extra edits.
* ITRANC  type of transport corrections in the microlib
*         (=0: no transport correction).
* MAXNFI  maximum number of fissionable isotopes in a mixture.
* NPART   number of companion particles.
* LSAME   fission spectrum flag (=.true. if all the isotopes have the
*         same fission spectrum and the same precursor group decay
*         constants).
* ITSTMP  type of cross section perturbation (=0: perturbation
*         forbidden; =1: perturbation not used even if present;
*         =2: perturbation used if present).
* TMPDAY  time stamp in day/burnup/irradiation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGROUP,NBISO,NBMIX,NL,NDEL,NESP,ISONAM(3,NBISO),
     1 MIX(NBISO),NED,ITRANC,MAXNFI,NPART,NAMEAD(2,NED),ITSTMP
      REAL DEN(NBISO),TMPDAY(3)
      LOGICAL MASK(NBMIX),MASKL(NGROUP),LSAME
*----
*  LOCAL VARIABLES
*----
      INTEGER NBLK,NSTATE,IOUT,MAXESP
      PARAMETER (NBLK=50,NSTATE=40,IOUT=6,MAXESP=4)
      CHARACTER CM*4,CV*12,HSMG*131,TEXT12*12,HCM(0:10)*2,NORD(3)*4,
     1 TEXT2*2,HPRT1*1
      LOGICAL EXIST,MASKK,LOGL,LALL,LWP1,LSTRD,LH,LC,LOVERV,LDIFF,
     1 LFISS,LWT0,LWT1
      INTEGER IDATA(NSTATE),IESP(MAXESP+1),I,J,I0,IOF,IOF0,IP,IPOSDE,
     1 IPASS,ISP,IGR,IG1,LLL,LLL0,IGMIN,IGMAX,IBM,JBM,IDEL,IED,IFIS,
     2 NXSPER,ISOT,IBLK,ILONG,LENGTZ,ITYLCM,IWFIS,IXSPER,KFIS,M,NBM0,
     3 NFISS0,NFISSI,NGROUPS
      REAL TMPPER(2,3),TIMFCT,DENISO,ENEAVG,FACT,TOTDEN,XTF
      DOUBLE PRECISION SQFMAS,XDRCST,NMASS,EVJ,ZNU
      TYPE(C_PTR) JPLIB,KPLIB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPOS,IJJ,NJJ,IWRK,NGPART
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NJJM,IJJM,INDFIS
      REAL, ALLOCATABLE, DIMENSION(:) :: GA1,GA2,SCAT,VOLMIX,NWTMIX,
     1 VOLI,C2PART
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GA3,GAR,WRK1,WRK2
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: GAF,CHECK,ZNUS,ZCHI,FLUX
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:,:) :: IPGRP
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LMADE
      CHARACTER(LEN=1), ALLOCATABLE, DIMENSION(:) :: HNPART
      REAL, ALLOCATABLE, DIMENSION(:) :: DENMIX !AHMED
*----
*  DATA STATEMENTS
*----
      DATA HCM/'00','01','02','03','04','05','06','07','08','09','10'/
      DATA NORD/'    ',' LIN',' QUA'/
*----
*  SCRATCH STORAGE ALLOCATION
*   IPGRP   LCM pointers of the macrolib groupwise directories.
*----
      ALLOCATE(NJJM(NBMIX,NBLK),IJJM(NBMIX,NBLK),IPOS(NBMIX),
     1 INDFIS(NBMIX,MAXNFI),IJJ(NGROUP),NJJ(NGROUP))
      ALLOCATE(GA2(NGROUP*NGROUP),GA3(NDEL,MAXNFI),GAR(NBMIX,NBLK+1),
     1 GAF(NBMIX,NGROUP,NBLK),SCAT(NGROUP*NBMIX),CHECK(NBMIX,NGROUP,NL),
     2 ZNUS(NBMIX*MAXNFI*NESP,NGROUP,0:NDEL),
     3 ZCHI(NBMIX*MAXNFI*NESP,NGROUP,0:NDEL))
      ALLOCATE(IPGRP(NGROUP,NPART+1))
      ALLOCATE(LMADE(NBISO))
      ALLOCATE(NGPART(NPART+1),C2PART(NPART+1),HNPART(NPART+1))
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
*----
*  SET MULTIPLE FISSION SPECTRA INFORMATION.
*----
      IF(NESP.GT.1) THEN
        IF(NESP.GT.MAXESP) CALL XABORT('LIBDEN: MAXESP OVERFLOW.')
        CALL LCMGET(IPLIB,'CHI-LIMITS',IESP)
      ENDIF
*----
*  SET CROSS SECTION PERTURBATION INFORMATION.
*----
      NXSPER=1
      TIMFCT=0.0
      CALL LCMLEN(IPLIB,'TIMESPER',ILONG,ITYLCM)
      IF((ILONG.GE.2).AND.(ILONG.LE.6)) THEN
        IF(ITSTMP.EQ.0) THEN
          CALL XABORT('LIBDEN: XS PERTURBATION FORBIDDEN.')
        ELSE IF(ITSTMP.EQ.2) THEN
          CALL LCMGET(IPLIB,'TIMESPER',TMPPER)
          TIMFCT=TMPDAY(1)-TMPPER(1,1)
          XTF=TIMFCT/TMPPER(2,1)
          IF(XTF.NE.0.0) NXSPER=2
          IF(XTF.LT.0.0) THEN
            WRITE(IOUT,6000) TMPPER(1,1),TMPDAY(1)
          ELSE IF(XTF.GT.1.0) THEN
            WRITE(IOUT,6001) TMPPER(1,1)+TMPPER(2,1),TMPDAY(1)
          ENDIF
        ENDIF
      ENDIF
*----
*  RECOVER MIXTURE VOLUMES IN MICROLIB.
*----
      CALL LCMLEN(IPLIB,'ISOTOPESVOL',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        ALLOCATE(VOLMIX(NBMIX),VOLI(NBISO))
        CALL LCMGET(IPLIB,'ISOTOPESVOL',VOLI)
        CALL XDRSET(VOLMIX,NBMIX,0.0)
        DO ISOT=1,NBISO
          IBM=MIX(ISOT)
          IF(IBM.GT.0) VOLMIX(IBM)=VOLI(ISOT)
        ENDDO
        CALL LCMPUT(IPLIB,'MIXTURESVOL',NBMIX,2,VOLMIX)
        DEALLOCATE(VOLI,VOLMIX)
      ENDIF

*----
*  MASKK=.TRUE. IF MIXTURE MASKING IS TO BE USED (IT IS NOT USED IF
*  ALL MIXTURES ARE TO BE UPDATED).
*----
      LDIFF=.TRUE.
      CALL LCMLEN(IPLIB,'MACROLIB',ILONG,ITYLCM)
      MASKK=(ILONG.EQ.-1)
      IF(MASKK) THEN
         CALL LCMSIX(IPLIB,'MACROLIB',1)
         CALL LCMGTC(IPLIB,'SIGNATURE',12,1,TEXT12)
         IF(TEXT12.NE.'L_MACROLIB') THEN
            CALL XABORT('LIBDEN: INVALID SIGNATURE ON THE MACROLIB.')
         ENDIF
         CALL LCMGET(IPLIB,'STATE-VECTOR',IDATA)
         NBM0=IDATA(2)
         NFISSI=IDATA(4)/NESP
         LDIFF=(IDATA(9).EQ.1)
         IF(IDATA(1).NE.NGROUP) THEN
            WRITE(HSMG,'(37HLIBDEN: EXISTING MACROLIB HAS NGROUP=,I4,
     1      25H NEW MACROLIB HAS NGROUP=,I4,1H.)') IDATA(1),NGROUP
            CALL XABORT(HSMG)
         ELSE IF((IDATA(6).NE.2).AND.(ITRANC.GT.0)) THEN
            WRITE(HSMG,'(37HLIBDEN: EXISTING MACROLIB HAS ITRANC=,I4,
     1      25H NEW MACROLIB HAS ITRANC=,I4,1H.)') IDATA(6),ITRANC
            CALL XABORT(HSMG)
         ELSE IF(NBM0.GT.NBMIX) THEN
            WRITE(HSMG,'(36HLIBDEN: EXISTING MACROLIB HAS NBMIX=,I4,
     1      24H NEW MACROLIB HAS NBMIX=,I4,1H.)') NBM0,NBMIX
            CALL XABORT(HSMG)
         ELSE IF(NFISSI.GT.NBISO) THEN
            WRITE(HSMG,'(37HLIBDEN: EXISTING MACROLIB HAS NFISSI=,I4,
     1      13H GREATER THAN,I5,1H.)') IDATA(4),NBISO
            CALL XABORT(HSMG)
         ENDIF
         IF(NFISSI.GT.0) THEN
            CALL LCMLEN(IPLIB,'FISSIONINDEX',ILONG,ITYLCM)
            IF(ILONG.EQ.0) THEN
*              THE NAMES ARE NOT DEFINED.
               DO 11 IFIS=1,NFISSI
               DO 10 IBM=1,NBMIX
               INDFIS(IBM,IFIS)=0
   10          CONTINUE
   11          CONTINUE
            ELSE IF(ILONG.EQ.NFISSI*NBMIX) THEN
               CALL LCMGET(IPLIB,'FISSIONINDEX',INDFIS)
               DO 16 IFIS=1,NFISSI
               DO 15 IBM=1,NBMIX
               IF(INDFIS(IBM,IFIS).GT.NBISO) THEN
                  CALL XABORT('LIBDEN: INVALID RECORD FISSIONINDEX.')
               ENDIF
   15          CONTINUE
   16          CONTINUE
            ELSE IF(ILONG.LT.NFISSI*NBMIX) THEN
*              REORDER THE 'FISSIONINDEX' MATRIX.
               ALLOCATE(IWRK(ILONG))
               CALL LCMGET(IPLIB,'FISSIONINDEX',IWRK)
               DO 31 IFIS=1,NFISSI
               DO 20 IBM=1,NBM0
               INDFIS(IBM,IFIS)=IWRK((IFIS-1)*NBM0+IBM)
   20          CONTINUE
               DO 30 IBM=NBM0+1,NBMIX
               INDFIS(IBM,IFIS)=0
   30          CONTINUE
   31          CONTINUE
               DEALLOCATE(IWRK)
            ELSE
               CALL XABORT('LIBDEN: INVALID NUMBER OF MIXTURES.')
            ENDIF
         ENDIF
         CALL LCMSIX(IPLIB,' ',2)
         LALL=NBMIX.GT.NBM0
      ELSE
         NFISSI=0
         LALL=.FALSE.
      ENDIF
*----
*  RECOVER PARTICLE DATA
*----
      !AHMED-CORRECTION
      !IF(NPART.EQ.0) THEN
      !  HPRT1=' '
      !  HNPART(1)=' '
      !ELSE
        CALL LCMGTC(IPLIB,'PARTICLE',1,1,HPRT1)
        CALL LCMGTC(IPLIB,'PARTICLE-NAM',1,NPART+1,HNPART)
        CALL LCMGET(IPLIB,'PARTICLE-NGR',NGPART)
        CALL LCMGET(IPLIB,'PARTICLE-MC2',C2PART)
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL LCMPTC(IPLIB,'PARTICLE',1,1,HPRT1)
        CALL LCMPTC(IPLIB,'PARTICLE-NAM',1,NPART+1,HNPART)
        CALL LCMPUT(IPLIB,'PARTICLE-NGR',NPART+1,1,NGPART)
        CALL LCMPUT(IPLIB,'PARTICLE-MC2',NPART+1,2,C2PART)
        CALL LCMSIX(IPLIB,' ',2)
        IF(HPRT1.NE.HNPART(1)) THEN
          WRITE(HSMG,'(26HLIBDEN: LIBTART PARTICLE (,A1,11H) IS DIFFER,
     1    25HENT FROM PARTICLE-NAM(1)=,A1,1H.)') HPRT1,HNPART(1)
          CALL XABORT(HSMG)
        ENDIF
        DO 32 IP=2,NPART+1
          ALLOCATE(GA1(NGPART(IP)+1))
          CALL LCMGET(IPLIB,HNPART(IP)//'ENERGY',GA1)
          CALL LCMSIX(IPLIB,'MACROLIB',1)
          CALL LCMPUT(IPLIB,HNPART(IP)//'ENERGY',NGPART(IP)+1,2,GA1)
          CALL LCMSIX(IPLIB,' ',2)
          DEALLOCATE(GA1)
   32   CONTINUE
      !ENDIF
*----
*  SELECT NUMBER OF GROUPS TO PROCESS
*----
      NGROUPS=0
      DO 35 LLL=1,NGROUP
      IF(MASKL(LLL).OR.LALL) NGROUPS=NGROUPS+1
   35 CONTINUE
      IF(NGROUPS.EQ.0) GO TO 880
*----
*  CHECK IF ALL REQUIRED ISOTOPES ARE PRESENT IN THE MICROLIB
*----
      ALLOCATE(GA1(NGROUP+1))
      DO 40 ISOT=1,NBISO
      IF(.NOT.MASK(MIX(ISOT))) GO TO 40
      JPLIB=IPISO(ISOT)
      IF(.NOT.C_ASSOCIATED(JPLIB)) THEN
         WRITE(HSMG,'(17HLIBDEN: ISOTOPE '',3A4,17H'' IS NOT AVAILABL,
     >   18HE IN THE MICROLIB.)') (ISONAM(I0,ISOT),I0=1,3)
         CALL XABORT(HSMG)
      ENDIF
   40 CONTINUE
*----
*  SET THE LCM MACROLIB GROUPWISE AND MICROLIB ISOTOPEWISE DIRECTORIES
*----
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      JPLIB=LCMLID(IPLIB,'GROUP',NGROUP)
      DO 45 LLL=1,NGROUP
      IPGRP(LLL,1)=LCMDIL(JPLIB,LLL)
   45 CONTINUE
      DO 47 IP=2,NPART+1
      JPLIB=LCMLID(IPLIB,'GROUP-'//HNPART(IP),NGROUP)
      DO 46 LLL=1,NGROUP
      IPGRP(LLL,IP)=LCMDIL(JPLIB,LLL)
   46 CONTINUE
   47 CONTINUE
      CALL LCMSIX(IPLIB,' ',2)
*----
*  PROCESS THE SCATTERING TABLES.
*----
      DO 52 I=1,NGROUP
      DO 51 IBM=1,NBMIX
      DO 50 J=1,NL
      CHECK(IBM,I,J)=0.0
   50 CONTINUE
   51 CONTINUE
   52 CONTINUE
      DO 245 IP=1,NPART+1
      DO 240 M=1,NL
      IF(M.LE.11) THEN
         CM=HCM(M-1)//'  '
      ELSE
         WRITE(CM,'(I2.2,2X)') M-1
      ENDIF
      DO 235 IPASS=0,(NGROUP-1)/NBLK
      LLL0=IPASS*NBLK
      DO 70 IBLK=1,NBLK
      DO 60 IBM=1,NBMIX
      GAR(IBM,IBLK)=0.0
   60 CONTINUE
      DO 71 LLL=1,NGROUP
      DO 72 IBM=1,NBMIX
      GAF(IBM,LLL,IBLK)=0.0
   72 CONTINUE
   71 CONTINUE
   70 CONTINUE
      DO 80 ISOT=1,NBISO
      LMADE(ISOT)=DEN(ISOT).EQ.0.0
   80 CONTINUE
      DO 140 ISOT=1,NBISO
      IF(LMADE(ISOT)) GO TO 140
      JPLIB=IPISO(ISOT)
      IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 140
*
*     RECOVER THE MICROSCOPIC TRANSFER XS WITHOUT USING XDRLGS (IN
*     ORDER TO REDUCE CPU TIME)
      IF(IP.GE.2) CALL LCMSIX(JPLIB,HNPART(IP),1)
      FACT=1.0
      DO 135 IXSPER=1,NXSPER
      CALL LCMLEN(JPLIB,'SIGS'//CM//NORD(IXSPER),ILONG,ITYLCM)
      IF(ILONG.EQ.0) GO TO 130
      CALL LCMGET(JPLIB,'SIGS'//CM//NORD(IXSPER),GA1)
      CALL LCMGET(JPLIB,'NJJS'//CM//NORD(IXSPER),NJJ)
      CALL LCMGET(JPLIB,'IJJS'//CM//NORD(IXSPER),IJJ)
      CALL LCMGET(JPLIB,'SCAT'//CM//NORD(IXSPER),GA2)
      IOF0=0
      DO 90 LLL=1,LLL0
      IOF0=IOF0+NJJ(LLL)
   90 CONTINUE
      DO 110 IBM=1,NBMIX
      IF((MASK(IBM).OR.(.NOT.MASKK)).AND.(MIX(ISOT).EQ.IBM)) THEN
         IOF=IOF0
         DO 105 IBLK=1,NBLK
         LLL=LLL0+IBLK
         IF(LLL.GT.NGROUP) GO TO 110
         GAR(IBM,IBLK)=GAR(IBM,IBLK)+GA1(LLL)*DEN(ISOT)
         DO 100 IG1=IJJ(LLL),IJJ(LLL)-NJJ(LLL)+1,-1
         IOF=IOF+1
         GAF(IBM,IG1,IBLK)=GAF(IBM,IG1,IBLK)+GA2(IOF)*DEN(ISOT)*FACT
  100    CONTINUE
  105    CONTINUE
      ENDIF
  110 CONTINUE
      LMADE(ISOT)=.TRUE.
  130 FACT=FACT*TIMFCT
  135 CONTINUE
      IF(IP.GE.2) CALL LCMSIX(JPLIB,' ',2)
*-
  140 CONTINUE
      DO 230 IBLK=1,NBLK
      LLL=LLL0+IBLK
      IF(LLL.GT.NGROUP) GO TO 230
      KPLIB=IPGRP(LLL,IP)
      IF(MASKL(LLL).OR.LALL) THEN
         IF(MASKK) THEN
           ILONG=1
           IF(M.GT.1) CALL LCMLEN(KPLIB,'SIGS'//CM,ILONG,ITYLCM)
           CALL XDRSET(GAR(1,NBLK+1),NBMIX,0.0)
           IF(ILONG.GT.0) THEN
              CALL LCMGET(KPLIB,'SIGS'//CM,GAR(1,NBLK+1))
           ENDIF
           DO 150 IBM=1,NBMIX
           IF(.NOT.MASK(IBM)) GAR(IBM,IBLK)=GAR(IBM,NBLK+1)
  150      CONTINUE
         ENDIF
         CALL LCMPUT(KPLIB,'SIGS'//CM,NBMIX,2,GAR(1,IBLK))
      ENDIF
*
      LOGL=MASKL(LLL).OR.LALL
      DO 165 IBM=1,NBMIX
      DO 160 IG1=1,NGROUP
      LOGL=LOGL.OR.(MASKL(IG1).AND.(GAF(IBM,IG1,IBLK).NE.0.0))
  160 CONTINUE
  165 CONTINUE
      IF(LOGL) THEN
         IF(MASKK) THEN
           ILONG=1
           IF(M.GT.1) CALL LCMLEN(KPLIB,'SCAT'//CM,ILONG,ITYLCM)
           IF(ILONG.GT.0) THEN
             DO 170 I=1,NBMIX
             IPOS(I)=-99
  170        CONTINUE
             CALL LCMGET(KPLIB,'SCAT'//CM,SCAT)
             CALL LCMGET(KPLIB,'NJJS'//CM,NJJM(1,IBLK))
             CALL LCMGET(KPLIB,'IJJS'//CM,IJJM(1,IBLK))
             CALL LCMGET(KPLIB,'IPOS'//CM,IPOS)
             DO 190 IBM=1,NBMIX
             IF(.NOT.MASK(IBM)) THEN
                IPOSDE=IPOS(IBM)
                IF(IPOSDE.EQ.-99) GO TO 190
                DO 180 IG1=IJJM(IBM,IBLK),IJJM(IBM,IBLK)-NJJM(IBM,IBLK)
     1          +1,-1
                GAF(IBM,IG1,IBLK)=SCAT(IPOSDE)
                IPOSDE=IPOSDE+1
  180           CONTINUE
             ENDIF
  190        CONTINUE
           ENDIF
         ENDIF
*
         IPOSDE=0
         DO 220 IBM=1,NBMIX
         IPOS(IBM)=IPOSDE+1
         IGMIN=LLL
         IGMAX=LLL
         DO 200 IG1=NGROUP,1,-1
         IF(GAF(IBM,IG1,IBLK).NE.0.0) THEN
            IGMIN=MIN(IGMIN,IG1)
            IGMAX=MAX(IGMAX,IG1)
         ENDIF
  200    CONTINUE
         IJJM(IBM,IBLK)=IGMAX
         NJJM(IBM,IBLK)=IGMAX-IGMIN+1
         DO 210 IG1=IGMAX,IGMIN,-1
         IPOSDE=IPOSDE+1
         SCAT(IPOSDE)=GAF(IBM,IG1,IBLK)
         CHECK(IBM,IG1,M)=CHECK(IBM,IG1,M)+SCAT(IPOSDE)
  210    CONTINUE
         GAR(IBM,1)=SCAT(IPOS(IBM)+IJJM(IBM,IBLK)-LLL)
  220    CONTINUE
         CALL LCMPUT(KPLIB,'SCAT'//CM,IPOSDE,2,SCAT)
         CALL LCMPUT(KPLIB,'NJJS'//CM,NBMIX,1,NJJM(1,IBLK))
         CALL LCMPUT(KPLIB,'IJJS'//CM,NBMIX,1,IJJM(1,IBLK))
         CALL LCMPUT(KPLIB,'IPOS'//CM,NBMIX,1,IPOS)
         CALL LCMPUT(KPLIB,'SIGW'//CM,NBMIX,2,GAR(1,1))
      ENDIF
  230 CONTINUE
  235 CONTINUE
  240 CONTINUE
  245 CONTINUE
*----
*  PROCESS THE REACTION VECTORS TOTAL, TOTAL-P1, STRD, H-FACTOR,
*  C-FACTOR, OVERV AND TRANC.
*----
      LWP1=.FALSE.
      LSTRD=.FALSE.
      LH=.FALSE.
      LC=.FALSE.
      LOVERV=.FALSE.
      DO 340 IBM=1,NBMIX
      IF(MASK(IBM).OR.(.NOT.MASKK)) THEN
         DO 255 IP=1,12
         DO 250 LLL=1,NGROUP
         GAF(IBM,LLL,IP)=0.0
  250    CONTINUE
  255    CONTINUE
         TOTDEN=0.0
         DO 320 ISOT=1,NBISO
         IF((MIX(ISOT).NE.IBM).OR.(DEN(ISOT).EQ.0.0)) GO TO 320
         JPLIB=IPISO(ISOT)
         IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 320
*-
         DENISO=DEN(ISOT)
         TOTDEN=TOTDEN+DENISO
         DO 315 IXSPER=1,NXSPER
         CALL LCMGET(JPLIB,'NTOT0   '//NORD(IXSPER),GA1)
         DO 260 LLL=1,NGROUP
         GAF(IBM,LLL,1)=GAF(IBM,LLL,1)+GA1(LLL)*DENISO
  260    CONTINUE
         CALL LCMLEN(JPLIB,'NTOT1   '//NORD(IXSPER),ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            LWP1=.TRUE.
            CALL LCMGET(JPLIB,'NTOT1   '//NORD(IXSPER),GA1)
            DO 270 LLL=1,NGROUP
            GAF(IBM,LLL,3)=GAF(IBM,LLL,3)+GA1(LLL)*DENISO
  270       CONTINUE
         ENDIF
         IF(LDIFF) THEN
            CALL LCMLEN(JPLIB,'STRD    '//NORD(IXSPER),ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               LSTRD=.TRUE.
               CALL LCMGET(JPLIB,'STRD    '//NORD(IXSPER),GA1)
               DO 280 LLL=1,NGROUP
               GAF(IBM,LLL,5)=GAF(IBM,LLL,5)+GA1(LLL)*DENISO
  280          CONTINUE
            ENDIF
         ENDIF
         CALL LCMLEN(JPLIB,'H-FACTOR'//NORD(IXSPER),ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            LH=.TRUE.
            CALL LCMGET(JPLIB,'H-FACTOR'//NORD(IXSPER),GA1)
            DO 290 LLL=1,NGROUP
            GAF(IBM,LLL,7)=GAF(IBM,LLL,7)+GA1(LLL)*DENISO
  290       CONTINUE
         ENDIF
         CALL LCMLEN(JPLIB,'C-FACTOR'//NORD(IXSPER),ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            LC=.TRUE.
            CALL LCMGET(JPLIB,'C-FACTOR'//NORD(IXSPER),GA1)
            DO 295 LLL=1,NGROUP
            GAF(IBM,LLL,13)=GAF(IBM,LLL,13)+GA1(LLL)*DENISO
  295       CONTINUE
         ENDIF
         CALL LCMLEN(JPLIB,'OVERV   '//NORD(IXSPER),ILONG,ITYLCM)
         !-AHMED-CORRECTION-OVERV
         IF((ILONG.GT.0).AND.((HPRT1.EQ.'N').OR.(HPRT1.EQ.'NEUT'))) THEN
            LOVERV=.TRUE.
            CALL LCMGET(JPLIB,'OVERV   '//NORD(IXSPER),GA1)
            DO 300 LLL=1,NGROUP
            GAF(IBM,LLL,9)=GAF(IBM,LLL,9)+GA1(LLL)*DENISO
  300       CONTINUE
         ENDIF
         IF(ITRANC.NE.0) THEN
            CALL LCMGET(JPLIB,'TRANC   '//NORD(IXSPER),GA1)
            DO 310 LLL=1,NGROUP
            GAF(IBM,LLL,11)=GAF(IBM,LLL,11)+GA1(LLL)*DENISO
  310       CONTINUE
         ENDIF
         DENISO=DENISO*TIMFCT
  315    CONTINUE
*-
  320    CONTINUE
         IF(LOVERV) THEN
            DO 330 LLL=1,NGROUP
            IF(GAF(IBM,LLL,9).NE.0.0) THEN
               GAF(IBM,LLL,9)=GAF(IBM,LLL,9)/TOTDEN
            ENDIF
  330       CONTINUE
         ENDIF
      ENDIF
  340 CONTINUE
      DO 420 LLL=1,NGROUP
      KPLIB=IPGRP(LLL,1)
      IF(MASKL(LLL).OR.LALL) THEN
         IF(MASKK) THEN
            CALL XDRSET(GAF(1,LLL,2),NBMIX,0.0)
            CALL LCMGET(KPLIB,'NTOT0',GAF(1,LLL,2))
            DO 350 IBM=1,NBMIX
            IF(.NOT.MASK(IBM)) GAF(IBM,LLL,1)=GAF(IBM,LLL,2)
  350       CONTINUE
         ENDIF
         CALL LCMPUT(KPLIB,'NTOT0',NBMIX,2,GAF(1,LLL,1))
         IF(LWP1) THEN
            IF(MASKK) THEN
               CALL XDRSET(GAF(1,LLL,4),NBMIX,0.0)
               CALL LCMGET(KPLIB,'NTOT1',GAF(1,LLL,4))
               DO 360 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) GAF(IBM,LLL,3)=GAF(IBM,LLL,4)
  360          CONTINUE
            ENDIF
            CALL LCMPUT(KPLIB,'NTOT1',NBMIX,2,GAF(1,LLL,3))
         ENDIF
         IF(LSTRD) THEN
            IF(MASKK) THEN
               CALL XDRSET(GAF(1,LLL,6),NBMIX,0.0)
               CALL LCMGET(KPLIB,'DIFF',GAF(1,LLL,6))
               DO 370 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) THEN
                  GAF(IBM,LLL,5)=1.0/(3.0*GAF(IBM,LLL,6))
               ENDIF
  370          CONTINUE
            ENDIF
            DO 380 IBM=1,NBMIX
            IF(GAF(IBM,LLL,5).NE.0.0) THEN
               GAF(IBM,LLL,5)=1.0/(3.0*GAF(IBM,LLL,5))
            ENDIF
  380       CONTINUE
            CALL LCMPUT(KPLIB,'DIFF',NBMIX,2,GAF(1,LLL,5))
         ENDIF
         IF(LH) THEN
            IF(MASKK) THEN
               CALL XDRSET(GAF(1,LLL,8),NBMIX,0.0)
               CALL LCMGET(KPLIB,'H-FACTOR',GAF(1,LLL,8))
               DO 390 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) GAF(IBM,LLL,7)=GAF(IBM,LLL,8)
  390          CONTINUE
            ENDIF
            CALL LCMPUT(KPLIB,'H-FACTOR',NBMIX,2,GAF(1,LLL,7))
         ENDIF
         IF(LC) THEN
            IF(MASKK) THEN
               CALL XDRSET(GAF(1,LLL,14),NBMIX,0.0)
               CALL LCMGET(KPLIB,'C-FACTOR',GAF(1,LLL,14))
               DO 395 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) GAF(IBM,LLL,13)=GAF(IBM,LLL,14)
  395          CONTINUE
            ENDIF
            CALL LCMPUT(KPLIB,'C-FACTOR',NBMIX,2,GAF(1,LLL,13))
         ENDIF
         IF(LOVERV) THEN
            IF(MASKK) THEN
               CALL XDRSET(GAF(1,LLL,10),NBMIX,0.0)
               CALL LCMGET(KPLIB,'OVERV',GAF(1,LLL,10))
               DO 400 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) GAF(IBM,LLL,9)=GAF(IBM,LLL,10)
  400          CONTINUE
            ENDIF
            CALL LCMPUT(KPLIB,'OVERV',NBMIX,2,GAF(1,LLL,9))
         ENDIF
         IF(ITRANC.NE.0) THEN
            IF(MASKK) THEN
               CALL XDRSET(GAF(1,LLL,12),NBMIX,0.0)
               CALL LCMGET(KPLIB,'TRANC',GAF(1,LLL,12))
               DO 410 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) GAF(IBM,LLL,11)=GAF(IBM,LLL,12)
  410          CONTINUE
            ENDIF
            CALL LCMPUT(KPLIB,'TRANC',NBMIX,2,GAF(1,LLL,11))
         ENDIF
      ENDIF
  420 CONTINUE
*----
*  PROCESS THE FISSION VECTORS FOR EACH NEW FISSILE ISOTOPE.
*----
      NFISS0=NFISSI
      DO 460 ISOT=1,NBISO
      IBM=MIX(ISOT)
      IF(IBM.EQ.0) GO TO 460
      IF(MASK(IBM).OR.(.NOT.MASKK)) THEN
         JPLIB=IPISO(ISOT)
         IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 460
         CALL LCMLEN(JPLIB,'NUSIGF',ILONG,ITYLCM)
         IF(NESP.EQ.1) THEN
            CALL LCMLEN(JPLIB,'CHI',LENGTZ,ITYLCM)
         ELSE
            CALL LCMLEN(JPLIB,'CHI--01',LENGTZ,ITYLCM)
         ENDIF
         IF((ILONG.GT.0).AND.(LENGTZ.GT.0)) THEN
            IF(NESP.EQ.1) THEN
               CALL LCMGET(JPLIB,'CHI',GA1)
            ELSE
               CALL LCMGET(JPLIB,'CHI--01',GA1)
            ENDIF
            LFISS=.FALSE.
            DO 425 IGR=1,NGROUP
            LFISS=LFISS.OR.(GA1(IGR).GT.0.0)
  425       CONTINUE
            IF(.NOT.LFISS) GO TO 455
            DO 430 IFIS=1,NFISSI
            IWFIS=INDFIS(IBM,IFIS)
            IF((IWFIS.EQ.ISOT).OR.(IWFIS.EQ.0)) THEN
               KFIS=IFIS
               GO TO 450
            ENDIF
  430       CONTINUE
            NFISSI=NFISSI+1
            IF(NFISSI.GT.MAXNFI) CALL XABORT('LIBDEN: INDFIS IS FULL.')
            KFIS=NFISSI
            DO 440 JBM=1,NBMIX
            INDFIS(JBM,KFIS)=0
  440       CONTINUE
  450       INDFIS(IBM,KFIS)=ISOT
         ENDIF
  455    CONTINUE
      ENDIF
  460 CONTINUE
      IF(NFISS0.GT.0) THEN
         ALLOCATE(WRK1(NBM0,NFISS0*NESP),WRK2(NBM0,NFISS0*NESP))
         DO 480 LLL=1,NGROUP
         IF(MASKL(LLL).OR.LALL) THEN
            DO 465 IDEL=0,NDEL
            CALL XDRSET(ZNUS(1,LLL,IDEL),NBMIX*NFISS0*NESP,0.0)
            CALL XDRSET(ZCHI(1,LLL,IDEL),NBMIX*NFISS0*NESP,0.0)
  465       CONTINUE
            KPLIB=IPGRP(LLL,1)
            CALL LCMLEN(KPLIB,'NUSIGF',ILONG,ITYLCM)
            IF(ILONG.NE.NBM0*NFISS0*NESP) THEN
               CALL XABORT('LIBDEN: NBM ERROR.')
            ENDIF
            CALL LCMGET(KPLIB,'NUSIGF',WRK1)
            CALL LCMGET(KPLIB,'CHI',WRK2)
            DO 467 IFIS=1,NFISS0*NESP
            DO 466 IBM=1,NBM0
            ZNUS((IFIS-1)*NBMIX+IBM,LLL,0)=WRK1(IBM,IFIS)
            ZCHI((IFIS-1)*NBMIX+IBM,LLL,0)=WRK2(IBM,IFIS)
  466       CONTINUE
  467       CONTINUE
            DO 475 IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            CALL LCMLEN(KPLIB,TEXT12,ILONG,ITYLCM)
            IF(ILONG.NE.0) THEN
               CALL LCMGET(KPLIB,TEXT12,WRK1)
               WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
               CALL LCMGET(KPLIB,TEXT12,WRK2)
               DO 471 IFIS=1,NFISS0*NESP
               DO 470 IBM=1,NBM0
               ZNUS((IFIS-1)*NBMIX+IBM,LLL,IDEL)=WRK1(IBM,IFIS)
               ZCHI((IFIS-1)*NBMIX+IBM,LLL,IDEL)=WRK2(IBM,IFIS)
  470          CONTINUE
  471          CONTINUE
            ENDIF
  475       CONTINUE
         ENDIF
  480    CONTINUE
         DEALLOCATE(WRK2,WRK1)
      ENDIF
      IF(NFISSI.GT.0) THEN
         DO 525 ISP=1,NESP
         DO 520 KFIS=1,NFISSI
         IF(KFIS.GT.NFISS0*NESP) THEN
            DO 492 IDEL=0,NDEL
            DO 491 LLL=1,NGROUP
            DO 490 IBM=1,NBMIX
            IOF=(KFIS-1)*NBMIX*NESP+(ISP-1)*NBMIX+IBM
            ZNUS(IOF,LLL,IDEL)=0.0
            ZCHI(IOF,LLL,IDEL)=0.0
  490       CONTINUE
  491       CONTINUE
  492       CONTINUE
         ELSE
            DO 510 IBM=1,NBMIX
            IWFIS=INDFIS(IBM,KFIS)
            IF((IWFIS.NE.0).AND.(MASK(IBM).OR.(.NOT.MASKK))) THEN
               DO 505 IDEL=0,NDEL
               DO 500 LLL=1,NGROUP
               IOF=(KFIS-1)*NBMIX*NESP+(ISP-1)*NBMIX+IBM
               ZNUS(IOF,LLL,IDEL)=0.0
               ZCHI(IOF,LLL,IDEL)=0.0
  500          CONTINUE
  505          CONTINUE
            ENDIF
  510       CONTINUE
         ENDIF
  520    CONTINUE
  525    CONTINUE
*-
         IF(NESP.EQ.1) THEN
*          ONE FISSION SPECTRUM (CLASSICAL CASE)
           DO 585 KFIS=1,NFISSI
           DO 580 IBM=1,NBMIX
           IWFIS=INDFIS(IBM,KFIS)
           IF((IWFIS.NE.0).AND.(MASK(IBM).OR.(.NOT.MASKK))) THEN
              IF(LSAME) THEN
                 IOF=IBM
              ELSE
                 IOF=(KFIS-1)*NBMIX+IBM
              ENDIF
              JPLIB=IPISO(IWFIS)
              IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 580
*-
              DENISO=DEN(IWFIS)
              DO 570 IXSPER=1,NXSPER
              CALL LCMGET(JPLIB,'NUSIGF  '//NORD(IXSPER),GA1)
              DO 530 LLL=1,NGROUP
              ZNUS(IOF,LLL,0)=ZNUS(IOF,LLL,0)+GA1(LLL)*DENISO
  530         CONTINUE
              IF(NDEL.GT.0) THEN
                 WRITE(TEXT12,'(6HNUSIGF,I2.2,A4)') NDEL,NORD(IXSPER)
                 CALL LCMLEN(JPLIB,TEXT12,ILONG,ITYLCM)
                 IF(ILONG.GT.0) THEN
                    DO 545 IDEL=1,NDEL
                    WRITE(TEXT12,'(6HNUSIGF,I2.2,A4)') IDEL,NORD(IXSPER)
                    CALL LCMGET(JPLIB,TEXT12,GA1)
                    DO 540 LLL=1,NGROUP
                    ZNUS(IOF,LLL,IDEL)=ZNUS(IOF,LLL,IDEL)+GA1(LLL)*
     1              DENISO
  540               CONTINUE
  545               CONTINUE
                 ENDIF
                 WRITE(TEXT12,'(3HCHI,I2.2,3X,A4)') NDEL,NORD(IXSPER)
                 CALL LCMLEN(JPLIB,TEXT12,ILONG,ITYLCM)
                 IF((ILONG.GT.0).AND.(IXSPER.EQ.1)) THEN
                    DO 555 IDEL=1,NDEL
                    WRITE(TEXT12,'(3HCHI,I2.2,3X,A4)') IDEL,NORD(IXSPER)
                    CALL LCMGET(JPLIB,TEXT12,GA1)
                    DO 550 LLL=1,NGROUP
                    ZCHI(IOF,LLL,IDEL)=GA1(LLL)
  550               CONTINUE
  555               CONTINUE
                 ENDIF
              ENDIF
              IF(IXSPER.EQ.1) THEN
                 CALL LCMGET(JPLIB,'CHI     '//NORD(IXSPER),GA1)
                 DO 560 LLL=1,NGROUP
                 ZCHI(IOF,LLL,0)=GA1(LLL)
  560            CONTINUE
              ENDIF
              DENISO=DENISO*TIMFCT
  570         CONTINUE
           ENDIF
  580      CONTINUE
  585      CONTINUE
         ELSE
*          NESP>1 MULTIPLE FISSION SPECTRA CASE
           DO 662 ISP=1,NESP
           DO 661 KFIS=1,NFISSI
           DO 660 IBM=1,NBMIX
           IWFIS=INDFIS(IBM,KFIS)
           IF((IWFIS.NE.0).AND.(MASK(IBM).OR.(.NOT.MASKK))) THEN
              IF(LSAME) THEN
                 IOF=IBM
              ELSE
                 IOF=(KFIS-1)*NBMIX*NESP+(ISP-1)*NBMIX+IBM
              ENDIF
              JPLIB=IPISO(IWFIS)
              IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 660
*-
              DENISO=DEN(IWFIS)
              DO 650 IXSPER=1,NXSPER
              CALL LCMGET(JPLIB,'NUSIGF  '//NORD(IXSPER),GA1)
              DO 610 LLL=IESP(ISP)+1,IESP(ISP+1)
              ZNUS(IOF,LLL,0)=ZNUS(IOF,LLL,0)+GA1(LLL)*DENISO
  610         CONTINUE
              IF((NDEL.GT.0).AND.(ISP.EQ.1)) THEN
                 WRITE(TEXT12,'(6HNUSIGF,I2.2,A4)') NDEL,NORD(IXSPER)
                 CALL LCMLEN(JPLIB,TEXT12,ILONG,ITYLCM)
                 IF(ILONG.GT.0) THEN
                    DO 625 IDEL=1,NDEL
                    WRITE(TEXT12,'(6HNUSIGF,I2.2,A4)') IDEL,NORD(IXSPER)
                    CALL LCMGET(JPLIB,TEXT12,GA1)
                    DO 620 LLL=1,NGROUP
                    ZNUS(IOF,LLL,IDEL)=ZNUS(IOF,LLL,IDEL)+GA1(LLL)*
     1              DENISO
  620               CONTINUE
  625               CONTINUE
                 ENDIF
                 WRITE(TEXT12,'(3HCHI,I2.2,3X,A4)') NDEL,NORD(IXSPER)
                 CALL LCMLEN(JPLIB,TEXT12,ILONG,ITYLCM)
                 IF((ILONG.GT.0).AND.(IXSPER.EQ.1)) THEN
                    DO 635 IDEL=1,NDEL
                    WRITE(TEXT12,'(3HCHI,I2.2,3X,A4)') IDEL,NORD(IXSPER)
                    CALL LCMGET(JPLIB,TEXT12,GA1)
                    DO 630 LLL=1,NGROUP
                    ZCHI(IOF,LLL,IDEL)=GA1(LLL)
  630               CONTINUE
  635               CONTINUE
                 ENDIF
              ENDIF
              IF(IXSPER.EQ.1) THEN
                 WRITE(TEXT2,'(I2.2)') ISP
                 TEXT12='CHI--'//TEXT2//' '//NORD(IXSPER)
                 CALL LCMLEN(JPLIB,TEXT12,ILONG,ITYLCM)
                 IF(ILONG.EQ.NGROUP) THEN
                    CALL LCMGET(JPLIB,TEXT12,GA1)
                    DO 640 LLL=1,NGROUP
                    ZCHI(IOF,LLL,0)=GA1(LLL)
  640               CONTINUE
                 ENDIF
              ENDIF
              DENISO=DENISO*TIMFCT
  650         CONTINUE
           ENDIF
  660      CONTINUE
  661      CONTINUE
  662      CONTINUE
         ENDIF 
*-
         DO 680 LLL=1,NGROUP
         IF(MASKL(LLL).OR.LALL) THEN
            KPLIB=IPGRP(LLL,1)
            ILONG=NBMIX*NFISSI*NESP
            IF(LSAME) ILONG=NBMIX
            CALL LCMPUT(KPLIB,'NUSIGF',ILONG,2,ZNUS(1,LLL,0))
            CALL LCMPUT(KPLIB,'CHI',ILONG,2,ZCHI(1,LLL,0))
            DO 670 IDEL=1,NDEL
            WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
            CALL LCMPUT(KPLIB,TEXT12,ILONG,2,ZNUS(1,LLL,IDEL))
            WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
            CALL LCMPUT(KPLIB,TEXT12,ILONG,2,ZCHI(1,LLL,IDEL))
  670       CONTINUE
         ENDIF
  680    CONTINUE
      ENDIF
*----
*  PROCESS THE EXTRA VECTOR EDITS.
*----
      DO 770 IED=1,NED
      WRITE(CV,'(2A4)') (NAMEAD(I0,IED),I0=1,2)
      IF(CV(:2).EQ.'NW') GO TO 770
      IF(CV.EQ.'TRANC') GO TO 770
      IF((CV(:3).EQ.'BST').OR.(CV(:3).EQ.'CST')) GO TO 770
      EXIST=.FALSE.
      DO 740 IBM=1,NBMIX
      IF(MASK(IBM).OR.(.NOT.MASKK)) THEN
         DO 690 LLL=1,NGROUP
         GAF(IBM,LLL,1)=0.0
  690    CONTINUE
         DO 730 ISOT=1,NBISO
         IF((MIX(ISOT).NE.IBM).OR.(DEN(ISOT).EQ.0.0)) GO TO 730
         JPLIB=IPISO(ISOT)
         IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 730
*-
         DENISO=DEN(ISOT)
         DO 710 IXSPER=1,NXSPER
         CALL LCMLEN(JPLIB,CV(:8)//NORD(IXSPER),ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 720
         EXIST=.TRUE.
         CALL LCMGET(JPLIB,CV(:8)//NORD(IXSPER),GA1)
         DO 700 LLL=1,NGROUP
         GAF(IBM,LLL,1)=GAF(IBM,LLL,1)+GA1(LLL)*DENISO
  700    CONTINUE
         DENISO=DENISO*TIMFCT
  710    CONTINUE
*-
  720    CONTINUE
  730    CONTINUE
      ENDIF
  740 CONTINUE
      DO 760 LLL=1,NGROUP
      IF(MASKL(LLL).OR.LALL) THEN
         KPLIB=IPGRP(LLL,1)
         IF(MASKK) THEN
            CALL LCMLEN(KPLIB,CV,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               EXIST=.TRUE.
               CALL XDRSET(GAF(1,LLL,2),NBMIX,0.0)
               CALL LCMGET(KPLIB,CV,GAF(1,LLL,2))
               DO 750 IBM=1,NBMIX
               IF(.NOT.MASK(IBM)) GAF(IBM,LLL,1)=GAF(IBM,LLL,2)
  750          CONTINUE
            ENDIF
         ENDIF
         IF(EXIST) CALL LCMPUT(KPLIB,CV,NBMIX,2,GAF(1,LLL,1))
      ENDIF
  760 CONTINUE
  770 CONTINUE
*
      CALL LCMGET(IPLIB,'ENERGY',GA1)
      IF(GA1(NGROUP+1).EQ.0.0) GA1(NGROUP+1)=1.0E-5
      CALL LCMSIX(IPLIB,'MACROLIB',1)
      IF(NED.GT.0) CALL LCMPUT(IPLIB,'ADDXSNAME-P0',2*NED,3,NAMEAD)
      IF(MASKK) THEN
        CALL LCMGET(IPLIB,'STATE-VECTOR',IDATA)
        IDATA(2)=MAX(NBM0,NBMIX)
        IDATA(3)=MAX(IDATA(3),NL)
        IDATA(4)=NFISSI*NESP
        IDATA(5)=MAX(IDATA(5),NED)
      ELSE
        IDATA(1)=NGROUP
        IDATA(2)=NBMIX
        IDATA(3)=NL
        IDATA(4)=NFISSI*NESP
        IDATA(5)=NED
        TEXT12='L_MACROLIB'
        CALL LCMPTC(IPLIB,'SIGNATURE',12,1,TEXT12)
        CALL LCMPUT(IPLIB,'ENERGY',NGROUP+1,2,GA1)
      ENDIF
*----
*  COMPUTE 1/V (ENER IS IN EV, NEUTRON MASS IS IN KG)
*----
      !AHMED-CORRECTION-OVERV
      IF((.NOT.LOVERV).AND.((HPRT1.EQ.'N').OR.(HPRT1.EQ.'NEUT'))) THEN
         DO 800 LLL=1,NGROUP
         ENEAVG=SQRT(GA1(LLL)*GA1(LLL+1))
         ZNU=1.0/(SQRT(ENEAVG)*SQFMAS)
         DO 790 IBM=1,NBMIX
         GAR(IBM,1)=REAL(ZNU)
  790    CONTINUE
         KPLIB=IPGRP(LLL,1)
         CALL LCMPUT(KPLIB,'OVERV',NBMIX,2,GAR(1,1))
  800    CONTINUE
      ENDIF
      DEALLOCATE(GA1)
*----
*  SET THE STATE VECTOR
*----
      IF(LSAME) IDATA(4)=MIN(NFISSI*NESP,1)
      IDATA(6)=ITRANC
      IF(ITRANC.NE.0) IDATA(6)=2
      IDATA(7)=NDEL
      IDATA(8)=0
      IDATA(9)=0
      IF(LSTRD) IDATA(9)=1
      IDATA(10)=0
      IF(LWP1) IDATA(10)=1
      DO 810 I=11,NSTATE
      IDATA(I)=0
  810 CONTINUE
      CALL LCMLEN(IPLIB,'SPH',ILONG,ITYLCM)
      IF(ILONG.NE.0) IDATA(14)=1
      CALL LCMPUT(IPLIB,'TIMESTAMP',3,2,TMPDAY)
      CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,IDATA)
*----
*  RECOVER THE PRECURSOR DECAY CONSTANTS.
*----
      IF(NDEL*NFISSI.GT.0) THEN
         IF(NFISS0.GT.0) THEN
            CALL LCMLEN(IPLIB,'LAMBDA-D',ILONG,ITYLCM)
            IF(ILONG.EQ.0) THEN
               CALL XDRSET(GA3(1,1),NDEL,0.0)
            ELSE
               CALL LCMGET(IPLIB,'LAMBDA-D',GA3(1,1))
            ENDIF
         ENDIF
         DO 825 KFIS=NFISS0+1,NFISSI
         DO 820 IDEL=1,NDEL
         GA3(IDEL,KFIS)=0.0
  820    CONTINUE
  825    CONTINUE
         CALL LCMSIX(IPLIB,' ',2)
         DO 835 KFIS=1,NFISSI
         DO 830 IBM=1,NBMIX
         IWFIS=INDFIS(IBM,KFIS)
         IF((IWFIS.NE.0).AND.(MASK(IBM).OR.(.NOT.MASKK))) THEN
            JPLIB=IPISO(IWFIS)
            IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 830
            CALL LCMLEN(JPLIB,'LAMBDA-D',ILONG,ITYLCM)
            IF(LSAME.AND.(ILONG.GT.0)) THEN
               CALL LCMGET(JPLIB,'LAMBDA-D',GA3(1,1))
            ELSE IF(ILONG.GT.0) THEN
               CALL LCMGET(JPLIB,'LAMBDA-D',GA3(1,KFIS))
            ENDIF
         ENDIF
  830    CONTINUE
  835    CONTINUE
         CALL LCMSIX(IPLIB,'MACROLIB',1)
         IF(LSAME) THEN
            CALL LCMPUT(IPLIB,'LAMBDA-D',NDEL,2,GA3(1,1))
         ELSE
            CALL LCMPUT(IPLIB,'LAMBDA-D',NDEL*NFISSI,2,GA3)
         ENDIF
      ENDIF
*
      IF((NFISSI.GT.0).AND.(.NOT.LSAME)) THEN
        CALL LCMPUT(IPLIB,'FISSIONINDEX',NBMIX*NFISSI,1,INDFIS)
      ENDIF
*
      DO 850 LLL=1,NGROUP
      IF(MASKL(LLL).OR.LALL) THEN
         KPLIB=IPGRP(LLL,1)
         DO 840 M=0,NL-1
         IF(M.LE.10) THEN
            CM=HCM(M)//'  '
         ELSE
            WRITE(CM,'(I2.2,2X)') M
         ENDIF
         CALL LCMPUT(KPLIB,'CHECK'//CM,NBMIX,2,CHECK(1,LLL,M+1))
  840    CONTINUE
      ENDIF
  850 CONTINUE
      CALL LCMSIX(IPLIB,' ',2)
*----
*  RECOVER THE INTEGRATED FLUX
*----
      CALL LCMLEN(IPLIB,'MIXTURESVOL',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        ALLOCATE(VOLMIX(NBMIX),NWTMIX(NGROUP),FLUX(NBMIX,NGROUP,2))
        CALL LCMGET(IPLIB,'MIXTURESVOL',VOLMIX)
        LWT0=.FALSE.
        LWT1=.FALSE.
        CALL XDRSET(FLUX,NBMIX*NGROUP*2,0.0)
        DO 860 ISOT=1,NBISO
        IBM=MIX(ISOT)
        IF(IBM.GT.0) THEN
          JPLIB=IPISO(ISOT)
          IF(C_ASSOCIATED(JPLIB)) THEN
            CALL LCMLEN(JPLIB,'NWT0',ILONG,ITYLCM)
            IF(ILONG.EQ.NGROUP) THEN
              LWT0=.TRUE.
              CALL LCMGET(JPLIB,'NWT0',NWTMIX)
              DO IGR=1,NGROUP
                FLUX(IBM,IGR,1)=NWTMIX(IGR)*VOLMIX(IBM)
              ENDDO
            ENDIF
            CALL LCMLEN(JPLIB,'NWT1',ILONG,ITYLCM)
            IF(ILONG.EQ.NGROUP) THEN
              LWT1=.TRUE.
              CALL LCMGET(JPLIB,'NWT1',NWTMIX)
              DO IGR=1,NGROUP
                FLUX(IBM,IGR,2)=NWTMIX(IGR)*VOLMIX(IBM)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
  860   CONTINUE
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL LCMPUT(IPLIB,'VOLUME',NBMIX,2,VOLMIX)
        JPLIB=LCMGID(IPLIB,'GROUP')
        DO 870 IGR=1,NGROUP
        KPLIB=LCMGIL(JPLIB,IGR)
        IF(LWT0) CALL LCMPUT(KPLIB,'FLUX-INTG',NBMIX,2,FLUX(1,IGR,1))
        IF(LWT1) CALL LCMPUT(KPLIB,'FLUX-INTG-P1',NBMIX,2,FLUX(1,IGR,2))
  870   CONTINUE
        CALL LCMSIX(IPLIB,' ',2)
        DEALLOCATE(FLUX,NWTMIX,VOLMIX)
      ENDIF
!AHMED
*----
*  RECOVER MIXTURES DENSITIES
*----
*        CALL LCMLEN(IPLIB,'DENMIXTURES',ILONG,ITYLCM) 
*        IF (ILONG.GT.0) THEN
*           ALLOCATE(DENMIX(ILONG)) 
*           CALL LCMGET(IPLIB,'DENMIXTURES',DENMIX)
*           CALL LCMPUT(KPLIB,'DENMIXTURES',ILONG,ITYLCM,DENMIX)
*        ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  880 DEALLOCATE(HNPART,C2PART,NGPART)
      DEALLOCATE(LMADE)
      DEALLOCATE(IPGRP)
      DEALLOCATE(ZCHI,ZNUS,CHECK,SCAT,GAF,GAR,GA3,GA2)
      DEALLOCATE(NJJ,IJJ,INDFIS,IPOS,IJJM,NJJM)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(' WARNING IN LIBDEN FOR PERTURBATION'/
     >       ' EXTRAPOLATION BELOW PRETURBATION TABLES'/
     >       ' INITIAL TIME       = ',F15.6,' DAYS'/
     >       ' EXTRAPOLATION TIME = ',F15.6,' DAYS')
 6001 FORMAT(' WARNING IN LIBDEN FOR PERTURBATION'/
     >       ' EXTRAPOLATION ABOVE PRETURBATION TABLES'/
     >       ' FINAL TIME         = ',F15.6,' DAYS'/
     >       ' EXTRAPOLATION TIME = ',F15.6,' DAYS')
      END
