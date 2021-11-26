*DECK LIBPTW
      SUBROUTINE LIBPTW (IPLIB,IPTMP,IPROC,NGRO,NL,HNAMIS,NED,HVECT,
     1 NDIL,DILUT,AWR,IPRECI,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build a temperature-independant draglib (IPROC=2) on the internal
* library or write the probability table information (IPROC=3/4).
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
* IPLIB   pointer to the isotopic directory in microlib.
* IPTMP   pointer to the multi-dilution internal library.
* IPROC   type of microlib processing:
*         =1: perform temperature interpolation and compute physical
*             probability tables.
*         =2: perform temperature interpolation and build a
*             temperature-independent draglib;
*         =3: perform temperature interpolation and compute calendf-
*             type mathematical probability tables based on bin-type
*             cross-sections for total cross sections;
*         =4: perform temperature interpolation and compute physical
*             probability tables or slowing-down correlated probability
*             tables.
*         =5: perform temperature interpolation and compute calendf-
*             type mathematical probability tables based on bin-type
*             cross-sections for all available cross-sections types.
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* HNAMIS  local name of the isotope:
*         HNAMIS(1:8)  is the local isotope name;
*         HNAMIS(9:12) is a suffix function of the mix number.
* NED     number of extra vector edits.
* HVECT   names of the extra vector edits.
* NDIL    number of finite dilutions.
* DILUT   dilutions.
* AWR     mass ratio for current isotope.
* IPRECI  accuracy index for probability tables in CALENDF.
* IMPX    print flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPTMP
      INTEGER IPROC,NGRO,NL,NED,NDIL,IPRECI,IMPX
      REAL DILUT(NDIL+1),AWR
      CHARACTER HNAMIS*12,HVECT(NED)*8
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIB,KPLIB,JPTMP,KPTMP
      PARAMETER (MAXNPT=12,MAXTRA=10000)
      CHARACTER TEXT12*12,TEXX12*12
      LOGICAL LSIGF,LGOLD,EMPTY,LCM,LBSIGF
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NFS,NOR,NBIN
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISMIN,ISMAX,ISM
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LSCAT,LADD
      REAL, ALLOCATABLE, DIMENSION(:) :: GOLD,FLUX,TOTAL,SIGF,SIGS,SCAT,
     1 SADD,ZDEL,DELTG,TBIN,SBIN,FBIN,EBIN
      TYPE(C_PTR) SIGP_PTR
      REAL, POINTER, DIMENSION(:) :: SIGP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISMIN(NL,NGRO),ISMAX(NL,NGRO),NFS(NGRO),NOR(NGRO),
     1 ISM(2,NL))
      ALLOCATE(LSCAT(NL),LADD(NED))
      ALLOCATE(GOLD(NGRO))
*----
*  COPY INFINITE DILUTION DATA FROM IPTMP TO IPLIB.
*----
      JPTMP=LCMGID(IPTMP,'ISOTOPESLIST')
      CALL LCMLEL(JPTMP,NDIL+1,ILENG,ITYLCM)
      IF(ILENG.EQ.0) THEN
         TEXT12=HNAMIS(1:8)
         WRITE(TEXT12(9:12),'(I4.4)') NDIL+1
         CALL XABORT('LIBPTW: MISSING LIST ITEM FOR '//TEXT12)
      ENDIF
      KPTMP=LCMGIL(JPTMP,NDIL+1) ! set (NDIL+1)-th isotope
      CALL LCMLEN(KPTMP,'LAMBDA-D',NDEL,ITYLCM)
      CALL LCMEQU(KPTMP,IPLIB)
*
      IF(NDIL.GT.0) THEN
*        RECOVER INFORMATION FROM IPTMP AND PUT NEW INFORMATION INTO
*        IPLIB.
         ALLOCATE(FLUX(NGRO*(NDIL+1)),TOTAL(NGRO*(NDIL+1)),
     1   SIGF(NGRO*(NDIL+1)),SIGS(NGRO*NL*(NDIL+1)),
     2   SCAT(NGRO*NGRO*NL*(NDIL+1)),SADD(NGRO*NED*(NDIL+1)),
     3   ZDEL(NGRO*NDEL*(NDIL+1)),DELTG(NGRO))
*
         CALL LIBEXT(IPTMP,NGRO,NL,NDIL,NED,HVECT,NDEL,IMPX,DILUT,
     1   MDIL,LSCAT,LSIGF,LADD,LGOLD,FLUX,TOTAL,SIGF,SIGS,SCAT,SADD,
     2   ZDEL,DELTG,GOLD,ISMIN,ISMAX)
*
*        DESTROY THE MULTI-DILUTION INTERNAL LIBRARY.
         CALL LCMCL(IPTMP,2)
*
         IF(IPROC.EQ.1) THEN
*           COMPUTE THE PHYSICAL PROBABILITY TABLES.
            MAXNOR=MAXNPT
            CALL LCMINF(IPLIB,TEXT12,TEXX12,EMPTY,ILONG,LCM)
            CALL LCMSIX(IPLIB,'PT-PHYS',1)
            CALL LCMPUT(IPLIB,'NDEL',1,1,NDEL)
            JPLIB=LCMLID(IPLIB,'GROUP',NGRO)
            DO 20 IGRP=1,NGRO
            NPART=3+NL+NED+NDEL
            DO 10 IL=1,NL
            NPART=NPART+MAX(ISMAX(IL,IGRP)-ISMIN(IL,IGRP)+1,0)
   10       CONTINUE
            IF(LGOLD) THEN
               GOLD0=GOLD(IGRP)
            ELSE
               GOLD0=1.0
            ENDIF
            IF(LCM) THEN
               SIGP_PTR=LCMARA(MAXNOR*NPART)
               CALL C_F_POINTER(SIGP_PTR,SIGP,(/ MAXNOR*NPART /))
            ELSE
               ALLOCATE(SIGP(MAXNOR*NPART))
            ENDIF
            CALL XDRSET(SIGP,MAXNOR*NPART,0.0)
            CALL LIBTAB(IGRP,NGRO,NL,MDIL,NPART,NED,NDEL,HNAMIS,IMPX,
     1      LSCAT,LSIGF,LADD,DILUT,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,GOLD0,
     2      ISMIN,ISMAX,NOR(IGRP),SIGP)
*
            IF(NOR(IGRP).GT.1) THEN
*              SAVE THE PROBABILITY TABLE INTO IPLIB.
               KPLIB=LCMDIL(JPLIB,IGRP)
               IF(LCM) THEN
                 CALL LCMPPD(KPLIB,'PROB-TABLE',MAXNOR*NPART,2,SIGP_PTR)
               ELSE
                 CALL LCMPUT(KPLIB,'PROB-TABLE',MAXNOR*NPART,2,SIGP)
                 DEALLOCATE(SIGP)
               ENDIF
               DO 15 IL=1,NL
               ISM(1,IL)=ISMIN(IL,IGRP)
               ISM(2,IL)=ISMAX(IL,IGRP)
   15          CONTINUE
               CALL LCMPUT(KPLIB,'ISM-LIMITS',2*NL,1,ISM)
            ELSE
               IF(LCM) THEN
                 CALL LCMDRD(SIGP_PTR)
               ELSE
                 DEALLOCATE(SIGP)
               ENDIF
            ENDIF
   20       CONTINUE
            CALL LCMPUT(IPLIB,'NOR',NGRO,1,NOR)
            CALL LCMSIX(IPLIB,' ',2)
         ELSE IF(IPROC.EQ.2) THEN
*           BUILD A TEMPERATURE-INDEPENDENT DRAGLIB.
            CALL LIBNOT(IPLIB,NGRO,NL,MDIL,NED,NDEL,IMPX,LSCAT,LSIGF,
     1      LADD,DILUT,FLUX,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,HVECT)
         ELSE IF((IPROC.EQ.3).OR.(IPROC.EQ.4).OR.(IPROC.EQ.5)) THEN
*           COMPUTE PHYSICAL PROBABILITY TABLES OR PROBABILITY TABLES
*           TAKING INTO ACCOUNT SLOWING-DOWN EFFECTS.
*
*           RECOVER BIN TYPE INFORMATION (IF AVAILABLE).
            LBSIGF=.FALSE.
            CALL LCMLEN(IPLIB,'BIN-NFS',LENGT,ITYLCM)
            LBIN=0
            IF(LENGT.GT.0) THEN
               CALL LCMGET(IPLIB,'BIN-NFS',NFS)
               IGRMIN=1
               IGRMAX=NGRO
               DO 30 IGRP=NGRO,1,-1
               IF((IGRMAX.EQ.IGRP).AND.(NFS(IGRP).EQ.0)) IGRMAX=IGRP-1
               LBIN=LBIN+NFS(IGRP)
   30          CONTINUE
               DO 40 IGRP=1,NGRO
               IF((IGRMIN.EQ.IGRP).AND.(NFS(IGRP).EQ.0)) IGRMIN=IGRP+1
   40          CONTINUE
               ALLOCATE(NBIN(NGRO),TBIN(LBIN),SBIN(LBIN),FBIN(LBIN),
     1         EBIN(LBIN+1))
               CALL LCMGET(IPLIB,'BIN-ENERGY',EBIN)
               CALL LCMGET(IPLIB,'BIN-NTOT0',TBIN)
               CALL LCMGET(IPLIB,'BIN-SIGS00',SBIN)
               CALL LCMLEN(IPLIB,'BIN-NUSIGF',LENGT,ITYLCM)
               IF(LENGT.GT.0) THEN 
                  CALL LCMGET(IPLIB,'BIN-NUSIGF',FBIN)
                  LBSIGF=.TRUE.
               ENDIF
            ELSE
               CALL XDISET(NFS,NGRO,0)
               IGRMIN=1
               IGRMAX=0
            ENDIF
*
*           RECOVER SCATTERING MATRIX PROFILE.
            LPART=0
            DO 55 IL=1,NL
            DO 50 IG1=1,NGRO
            LPART=MAX(LPART,ISMAX(IL,IG1)-ISMIN(IL,IG1)+1)
   50       CONTINUE
   55       CONTINUE
*
            CALL LIBFQD(MAXNPT,LPART,MAXTRA,HNAMIS,IPLIB,NGRO,NL,NED,
     1      NDEL,MDIL,IGRMIN,IGRMAX,LBIN,NFS,IMPX,LSCAT,LSIGF,LADD,
     2      DILUT,FLUX,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,EBIN,TBIN,SBIN,
     3      FBIN,AWR,ISMIN,ISMAX,GOLD,IPRECI,NOR,LBSIGF)
*           
            IF(LBIN.GT.0) THEN
               DEALLOCATE(EBIN,SBIN,TBIN,NBIN)
               IF(LBSIGF) DEALLOCATE(FBIN)
            ENDIF
         ELSE
            CALL XABORT('LIBPTW: INVALID VALUE OF IPROC.')
         ENDIF
         DEALLOCATE(DELTG,ZDEL,SADD,SCAT,SIGS,SIGF,TOTAL,FLUX)
      ELSE
*        DESTROY THE MULTI-DILUTION INTERNAL LIBRARY.
         CALL LCMCL(IPTMP,2)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GOLD)
      DEALLOCATE(LADD,LSCAT)
      DEALLOCATE(ISM,NOR,NFS,ISMAX,ISMIN)
*
      RETURN
      END
