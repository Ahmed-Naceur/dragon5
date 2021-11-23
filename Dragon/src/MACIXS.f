*DECK MACIXS
      SUBROUTINE MACIXS(IPLIST,MAXFIS,NGROUP,NBMIX,NIFISS,NANISO,NDELG,
     >                  XSTOTL,XSTOT1,XSFISS,XSSPEC,XSFIXE,XSTRAN,
     >                  XSDIFF,XSNFTO,XSH,XSSCAT,LOLDXS,ISCATA,XSNUDL,
     >                  XSCHDL,XSDIFX,XSDIFY,XSDIFZ,XSOVRV,XSINT0,
     >                  XSINT1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Cross sections initialization.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPLIST  LCM pointer to the macrolib.
* MAXFIS  set to max(1,NIFISS).
* NGROUP  number of energy groups.
* NBMIX   maximum number of mixtures.
* NIFISS  number of fissile isotopes.
* NANISO  maximum Legendre order:
*         =1 isotropic collision;
*         =2 linearly anisotropic collision.
* NDELG   number of precursor groups for delayed neutrons.
*
*Parameters: output
* XSTOTL  P0 total cross section of mixture.
* XSTOT1  P1 total cross section of mixture.
* XSFISS  nu*fission cross section of mixture.
* XSSPEC  fission spectrum.
* XSFIXE  fixed sources.
* XSTRAN  transport correction.
* XSDIFF  isotropic diffusion coefficient.
* XSNFTO  fission cross section of mixture.
* XSH     power factor (h-factor).
* XSSCAT  scattering cross section of mixture/group.
* XSNUDL  delayed nu*fission cross section of mixture.
* XSCHDL  delayed neutron fission spectrum.
* XSDIFX  x-directed diffusion coefficients.
* XSDIFY  y-directed diffusion coefficients.
* XSDIFZ  z-directed diffusion coefficients.
* XSOVRV  reciprocal neutron velocities.
* XSINT0  P0 volume-integrated flux of mixture.
* XSINT1  P1 volume-integrated flux of mixture.
* LOLDXS  flag to check if cross section type is already present on
*         the macrolib.
* ISCATA  check for scattering anisotropy.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIST
      INTEGER   MAXFIS,NGROUP,NBMIX,NIFISS,NANISO,NDELG,ISCATA(NANISO)
      REAL      XSTOTL(NBMIX,NGROUP),XSTOT1(NBMIX,NGROUP),
     >          XSFISS(NBMIX,MAXFIS,NGROUP),XSSPEC(NBMIX,MAXFIS,NGROUP),
     >          XSFIXE(NBMIX,NGROUP),XSTRAN(NBMIX,NGROUP),
     >          XSDIFF(NBMIX,NGROUP),XSNFTO(NBMIX,NGROUP),
     >          XSH(NBMIX,NGROUP),XSSCAT(NGROUP,NBMIX,NANISO,NGROUP),
     >          XSNUDL(NBMIX,MAXFIS,NDELG,NGROUP),
     >          XSCHDL(NBMIX,MAXFIS,NDELG,NGROUP),
     >          XSDIFX(NBMIX,NGROUP),XSDIFY(NBMIX,NGROUP),
     >          XSDIFZ(NBMIX,NGROUP),XSOVRV(NBMIX,NGROUP),
     >          XSINT0(NBMIX,NGROUP),XSINT1(NBMIX,NGROUP)
      LOGICAL   LOLDXS(18)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIST,KPLIST
      CHARACTER CANISO*2,NAMREC*12,CHID*12,NUSIGD*12
*----
* ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGSCT,IFGSCT
      REAL, ALLOCATABLE, DIMENSION(:) :: XSWORK
*----
*  SCRATCH STORAGE ALLOCATION
*   INGSCT   number of scattering group for cross sections.
*   IFGSCT   first scattering group for cross sections.
*   XSWORK   work cross-section vector.
*----
      ALLOCATE(INGSCT(NBMIX),IFGSCT(NBMIX))
      ALLOCATE(XSWORK(NBMIX*NGROUP))
*----
*  READ/INITIALIZE MACROLIB CROSS SECTION DATA
*----
      CALL XDRSET(XSSCAT,NBMIX*NGROUP*NGROUP*NANISO,0.0)
      CALL XDRSET(XSTOTL,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSTOT1,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSFISS,NBMIX*NIFISS*NGROUP,0.0)
      CALL XDRSET(XSSPEC,NBMIX*NIFISS*NGROUP,0.0)
      CALL XDRSET(XSFIXE,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSTRAN,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSDIFF,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSNFTO,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSH,NBMIX*NGROUP,0.0)
      IF(NDELG.GT.0) THEN
         CALL XDRSET(XSCHDL,NBMIX*NIFISS*NDELG*NGROUP,0.0)
         CALL XDRSET(XSNUDL,NBMIX*NIFISS*NDELG*NGROUP,0.0)
      ENDIF
      CALL XDRSET(XSOVRV,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSINT0,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSINT1,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSDIFX,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSDIFY,NBMIX*NGROUP,0.0)
      CALL XDRSET(XSDIFZ,NBMIX*NGROUP,0.0)
      JPLIST=LCMLID(IPLIST,'GROUP',NGROUP)
      DO 200 IGROUP=1,NGROUP
        KPLIST=LCMDIL(JPLIST,IGROUP)
*----
*  READ OR INITIALISE CHI AND NUSIGF
*----
        CALL LCMLEN(KPLIST,'NUSIGF',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(2)=.TRUE.
          CALL LCMGET(KPLIST,'NUSIGF',XSFISS(1,1,IGROUP))
        ELSE
          CALL XDRSET(XSFISS(1,1,IGROUP),NBMIX*NIFISS,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'CHI',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(4)=.TRUE.
          CALL LCMGET(KPLIST,'CHI',XSSPEC(1,1,IGROUP))
        ELSE
          CALL XDRSET(XSSPEC(1,1,IGROUP),NBMIX*NIFISS,0.0)
        ENDIF
*----
*  READ OR INITIALISE TOTAL XS, FIXED SOURCES AND TRANSPORT CORRECTION
*----
        CALL LCMLEN(KPLIST,'NTOT0',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(1)=.TRUE.
          CALL LCMGET(KPLIST,'NTOT0',XSTOTL(1,IGROUP))
        ELSE
          CALL XDRSET(XSTOTL(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'FIXE',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(3)=.TRUE.
          CALL LCMGET(KPLIST,'FIXE',XSFIXE(1,IGROUP))
        ELSE
          CALL XDRSET(XSFIXE(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'TRANC',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(6)=.TRUE.
          CALL LCMGET(KPLIST,'TRANC',XSTRAN(1,IGROUP))
        ELSE
          CALL XDRSET(XSTRAN(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'DIFF',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(7)=.TRUE.
          CALL LCMGET(KPLIST,'DIFF',XSDIFF(1,IGROUP))
        ELSE
          CALL XDRSET(XSDIFF(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'H-FACTOR',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(8)=.TRUE.
          CALL LCMGET(KPLIST,'H-FACTOR',XSH(1,IGROUP))
        ELSE
          CALL XDRSET(XSH(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'NTOT1',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(9)=.TRUE.
          CALL LCMGET(KPLIST,'NTOT1',XSTOT1(1,IGROUP))
        ELSE
          CALL XDRSET(XSTOT1(1,IGROUP),NBMIX,0.0)
        ENDIF
*----
*  READ OR INITIALISE DIFFX, DIFFY, DIFFZ, CHID AND OVERV
*----
        CALL LCMLEN(KPLIST,'DIFFX',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(10)=.TRUE.
          CALL LCMGET(KPLIST,'DIFFX',XSDIFX(1,IGROUP))
        ELSE
          CALL XDRSET(XSDIFX(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'DIFFY',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(11)=.TRUE.
          CALL LCMGET(KPLIST,'DIFFY',XSDIFY(1,IGROUP))
        ELSE
          CALL XDRSET(XSDIFY(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'DIFFZ',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(12)=.TRUE.
          CALL LCMGET(KPLIST,'DIFFZ',XSDIFZ(1,IGROUP))
        ELSE
          CALL XDRSET(XSDIFZ(1,IGROUP),NBMIX,0.0)
        ENDIF
        DO I=1,NDELG
          WRITE(NUSIGD,'(A6,I2.2)') 'NUSIGF',I
          CALL LCMLEN(KPLIST,NUSIGD,ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            LOLDXS(13)=.TRUE.
            CALL LCMGET(KPLIST,NUSIGD,XSNUDL(1,1,I,IGROUP))
          ELSE
            CALL XDRSET(XSNUDL(1,1,I,IGROUP),NBMIX*NIFISS,0.0)
          ENDIF
          WRITE(CHID,'(A3,I2.2)') 'CHI',I
          CALL LCMLEN(KPLIST,CHID,ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            LOLDXS(14)=.TRUE.
            CALL LCMGET(KPLIST,CHID,XSCHDL(1,1,I,IGROUP))
          ELSE
            CALL XDRSET(XSCHDL(1,1,I,IGROUP),NBMIX*NIFISS,0.0)
          ENDIF
        ENDDO
        CALL LCMLEN(KPLIST,'OVERV',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(15)=.TRUE.
          CALL LCMGET(KPLIST,'OVERV',XSOVRV(1,IGROUP))
        ELSE
          CALL XDRSET(XSOVRV(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'NFTOT',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(16)=.TRUE.
          CALL LCMGET(KPLIST,'NFTOT',XSNFTO(1,IGROUP))
        ELSE
          CALL XDRSET(XSNFTO(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'FLUX-INTG',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(17)=.TRUE.
          CALL LCMGET(KPLIST,'FLUX-INTG',XSINT0(1,IGROUP))
        ELSE
          CALL XDRSET(XSINT0(1,IGROUP),NBMIX,0.0)
        ENDIF
        CALL LCMLEN(KPLIST,'FLUX-INTG-P1',ILCMLN,ITYLCM)
        IF(ILCMLN.GT.0) THEN
          LOLDXS(18)=.TRUE.
          CALL LCMGET(KPLIST,'FLUX-INTG-P1',XSINT1(1,IGROUP))
        ELSE
          CALL XDRSET(XSINT1(1,IGROUP),NBMIX,0.0)
        ENDIF
*----
*  READ OR INITIALISE SCATTERING CROSS SECTIONS
*----
        DO 203 IANIS=1,NANISO
          WRITE(CANISO,'(I2.2)') IANIS-1
          NAMREC='SCAT'//CANISO
          CALL LCMLEN(KPLIST,NAMREC,ILCMLN,ITYLCM)
          ICMATR=0
          IF(ILCMLN.GT.0) THEN
            LOLDXS(5)=.TRUE.
            ISCATA(IANIS)=1
*----
*  READ COMPRESS SCATTERING XS PLUS INFORMATION TO EXPAND XS
*----
            CALL LCMGET(KPLIST,NAMREC,XSWORK)
            NAMREC='NJJS'//CANISO
            CALL LCMLEN(KPLIST,NAMREC,ICMATR,ITYLCM)
            CALL LCMGET(KPLIST,NAMREC,INGSCT)
            NAMREC='IJJS'//CANISO
            CALL LCMGET(KPLIST,NAMREC,IFGSCT)
*----
*  EXPAND SCATTERING XS TO XSSCAT(JGROUP,IMATER,IANIS,IGROUP)
*  WHERE IGROUP IS THE SECONDARY GROUP.
*----
            IPWRK=1
            DO 204 IMATER=1,ICMATR
              IF(INGSCT(IMATER).GT.0) THEN
                IGD=IFGSCT(IMATER)
                IGF=IGD-INGSCT(IMATER)+1
                DO 205 JGROUP=IGD,IGF,-1
                  XSSCAT(JGROUP,IMATER,IANIS,IGROUP)=XSWORK(IPWRK)
                  IPWRK=IPWRK+1
 205            CONTINUE
              ENDIF
 204        CONTINUE
          ENDIF
 203    CONTINUE
 200  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSWORK)
      DEALLOCATE(IFGSCT,INGSCT)
      RETURN
      END
