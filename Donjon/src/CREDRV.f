*DECK CREDRV
      SUBROUTINE CREDRV(IPMAC,IPMAP,NENTRY,HENTRY,KENTRY,LMAC,NMIX,
     1 NGRP,NL,ILEAK,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover and/or interpolate l_compo information, store properties
* in a new or existing macrolib.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert, M. Guyot
*
*Parameters: input/output
* IPMAC   pointer to the macrolib information.
* IPMAP   pointer to fuel-map information (=0 if no l_fmap).
* NENTRY  number of lcm or xsm objects used by the module.
* HENTRY  character*12 name of each lcm or xsm objects.
* KENTRY  pointers to the lcm or xsm objects.
* LMAC    flag for macrolib object type: =.false. in create mode;
*          =.true. in modification mode.
* NMIX    maximum number of material mixtures.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* ILEAK   diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* IMPX    printing index (=0 for no print).
*
*NOTE: a cross section not read is set to zero.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,NMIX,NGRP,NL,ILEAK,IMPX
      TYPE(C_PTR) IPMAC,IPMAP,KENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      LOGICAL LMAC
*----
*  LOCAL VARIABLES
*----
      CHARACTER CM*2
      TYPE(C_PTR) JPMAC,KPMAC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPOS
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IJJ,NJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TOTAL,ZNUG,SNUGF,CHI,OVERV,
     1 DIFFX,DIFFY,DIFFZ
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: H
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPOS(NMIX),IJJ(NMIX,NL,NGRP),NJJ(NMIX,NL,NGRP))
      ALLOCATE(TOTAL(NMIX,NGRP),ZNUG(NMIX,NGRP),SNUGF(NMIX,NGRP),
     1        CHI(NMIX,NGRP),OVERV(NMIX,NGRP),DIFFX(NMIX,NGRP),
     2        DIFFY(NMIX,NGRP),DIFFZ(NMIX,NGRP),WORK(NMIX*NGRP),
     3        SCAT(NMIX,NL,NGRP,NGRP),H(NMIX,NGRP,NL))
*
      IF((IMPX.GT.1).AND.LMAC)CALL LCMLIB(IPMAC)
      CALL XDRSET(TOTAL,NMIX*NGRP,0.)
      CALL XDRSET(ZNUG,NMIX*NGRP,0.)
      CALL XDRSET(SNUGF,NMIX*NGRP,0.)
      CALL XDRSET(DIFFX,NMIX*NGRP,0.)
      CALL XDRSET(DIFFY,NMIX*NGRP,0.)
      CALL XDRSET(DIFFY,NMIX*NGRP,0.)
      CALL XDRSET(CHI,NMIX*NGRP,0.)
      CALL XDRSET(OVERV,NMIX*NGRP,0.)
      CALL XDRSET(H,NMIX*NGRP,0.)
      CALL XDRSET(WORK,NMIX*NGRP,0.)
      CALL XDRSET(SCAT,NMIX*NL*NGRP*NGRP,0.)
      CALL XDISET(IPOS,NMIX,0)
      DO 12 IGR=1,NGRP
      DO 11 IBM=1,NMIX
      DO 10 IL=1,NL
      IJJ(IBM,IL,IGR)=IGR
      NJJ(IBM,IL,IGR)=1
   10 CONTINUE
   11 CONTINUE
   12 CONTINUE
*----
*  RECOVER THE EXISTING MACROLIB DATA
*----
      ILEAK=0
      IF(LMAC)THEN
        JPMAC=LCMGID(IPMAC,'GROUP')
        DO 40 JGR=1,NGRP
        KPMAC=LCMGIL(JPMAC,JGR)
        CALL LCMLEN(KPMAC,'NTOT0',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)THEN
          CALL LCMGET(KPMAC,'NTOT0',TOTAL(1,JGR))
        ELSEIF(ILENGT.NE.0)THEN
          CALL XABORT('@CREDRV: INVALID INPUT MACROLIB(1).')
        ENDIF
        CALL LCMLEN(KPMAC,'NUSIGF',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)CALL LCMGET(KPMAC,'NUSIGF',ZNUG(1,JGR))
        CALL LCMLEN(KPMAC,'NFTOT',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)CALL LCMGET(KPMAC,'NFTOT',SNUGF(1,JGR))
        CALL LCMLEN(KPMAC,'CHI',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)CALL LCMGET(KPMAC,'CHI',CHI(1,JGR))
        CALL LCMLEN(KPMAC,'OVERV',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)CALL LCMGET(KPMAC,'OVERV',OVERV(1,JGR))
        CALL LCMLEN(KPMAC,'DIFF',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)THEN
          ILEAK=1
          CALL LCMGET(KPMAC,'DIFF',DIFFX(1,JGR))
        ENDIF
        CALL LCMLEN(KPMAC,'DIFFX',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)THEN
          ILEAK=2
          CALL LCMGET(KPMAC,'DIFFX',DIFFX(1,JGR))
          CALL LCMGET(KPMAC,'DIFFY',DIFFY(1,JGR))
          CALL LCMGET(KPMAC,'DIFFZ',DIFFZ(1,JGR))
        ENDIF
        CALL LCMLEN(KPMAC,'H-FACTOR',ILENGT,ITYLCM)
        IF(ILENGT.EQ.NMIX)CALL LCMGET(KPMAC,'H-FACTOR',H(1,JGR,1))
        DO IL=1,NL
          WRITE (CM,'(I2.2)') IL-1
          CALL LCMLEN(KPMAC,'SCAT'//CM,ILENGT,ITYLCM)
          IF(ILENGT.GT.NMIX*NL*NGRP*NGRP)THEN
            CALL XABORT('@CREDRV: INVALID INPUT MACROLIB(2).')
          ELSEIF(ILENGT.GT.0)THEN
            CALL LCMGET(KPMAC,'SCAT'//CM,WORK)
            CALL LCMGET(KPMAC,'NJJS'//CM,NJJ(1,IL,JGR))
            CALL LCMGET(KPMAC,'IJJS'//CM,IJJ(1,IL,JGR))
            IPOSDE=0
            DO 25 IBM=1,NMIX
            IJJ0=IJJ(IBM,IL,JGR)
            DO 20 IGR=IJJ0,IJJ0-NJJ(IBM,IL,JGR)+1,-1
            IPOSDE=IPOSDE+1
            SCAT(IBM,IL,IGR,JGR)=WORK(IPOSDE)
   20       CONTINUE
   25       CONTINUE
          ELSE
            CALL XABORT('@CREDRV: OLD FORMAT OF THE MACROLIB.')
          ENDIF
        ENDDO
   40   CONTINUE
      ENDIF
*----
*  READ INPUT DATA
*----
      CALL CREXSI(IPMAP,NENTRY,HENTRY,KENTRY,NMIX,NGRP,NL,ILEAK,IMPX,
     1 TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,DIFFZ,H,IJJ,NJJ,SCAT)
*----
*  MACROLIB DATA STORAGE
*----
      JPMAC=LCMLID(IPMAC,'GROUP',NGRP)
      DO 190 JGR=1,NGRP
      KPMAC=LCMDIL(JPMAC,JGR)
      CALL LCMPUT(KPMAC,'NTOT0',NMIX,2,TOTAL(1,JGR))
      CALL LCMPUT(KPMAC,'NUSIGF',NMIX,2,ZNUG(1,JGR))
      CALL LCMPUT(KPMAC,'NFTOT',NMIX,2,SNUGF(1,JGR))
      CALL LCMPUT(KPMAC,'CHI',NMIX,2,CHI(1,JGR))
      CALL LCMPUT(KPMAC,'OVERV',NMIX,2,OVERV(1,JGR))
      IF(ILEAK.EQ.1)THEN
        CALL LCMPUT(KPMAC,'DIFF',NMIX,2,DIFFX(1,JGR))
      ELSEIF(ILEAK.EQ.2)THEN
        CALL LCMPUT(KPMAC,'DIFFX',NMIX,2,DIFFX(1,JGR))
        CALL LCMPUT(KPMAC,'DIFFY',NMIX,2,DIFFY(1,JGR))
        CALL LCMPUT(KPMAC,'DIFFZ',NMIX,2,DIFFZ(1,JGR))
      ENDIF
      CALL LCMPUT(KPMAC,'H-FACTOR',NMIX,2,H(1,JGR,1))
  190 CONTINUE
*----
*  SCATTERING DATA
*----
      CALL XDRSET(H,NMIX*NGRP*NL,0.0)
      DO 215 JGR=1,NGRP
      KPMAC=LCMDIL(JPMAC,JGR)
      DO 210 IL=1,NL
        WRITE (CM,'(I2.2)') IL-1
        IPOSDE=0
        DO 205 IBM=1,NMIX
        IPOS(IBM)=IPOSDE+1
        DO 200 IGR=IJJ(IBM,IL,JGR),IJJ(IBM,IL,JGR)-NJJ(IBM,IL,JGR)+1,-1
        IPOSDE=IPOSDE+1
        WORK(IPOSDE)=SCAT(IBM,IL,IGR,JGR)
        H(IBM,IGR,IL)=H(IBM,IGR,IL)+SCAT(IBM,IL,IGR,JGR)
  200   CONTINUE
  205   CONTINUE
        CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC,'IPOS'//CM,NMIX,1,IPOS)
        CALL LCMPUT(KPMAC,'NJJS'//CM,NMIX,1,NJJ(1,IL,JGR))
        CALL LCMPUT(KPMAC,'IJJS'//CM,NMIX,1,IJJ(1,IL,JGR))
        CALL LCMPUT(KPMAC,'SIGW'//CM,NMIX,2,SCAT(1,IL,JGR,JGR))
  210 CONTINUE
  215 CONTINUE
      DO 225 IGR=1,NGRP
      KPMAC=LCMDIL(JPMAC,IGR)
      DO 220 IL=1,NL
      WRITE (CM,'(I2.2)') IL-1
      CALL LCMPUT(KPMAC,'SIGS'//CM,NMIX,2,H(1,IGR,IL))
      IF(IMPX.GT.2)CALL LCMLIB(KPMAC)
  220 CONTINUE
  225 CONTINUE
*
      IF(IMPX.GT.1)CALL LCMLIB(IPMAC)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(H,SCAT,WORK,DIFFZ,DIFFY,DIFFX,OVERV,CHI,SNUGF,ZNUG,
     1 TOTAL)
      DEALLOCATE(NJJ,IJJ,IPOS)
      RETURN
      END
