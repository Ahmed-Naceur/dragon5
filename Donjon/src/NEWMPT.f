*DECK NEWMPT
      SUBROUTINE NEWMPT(IPMAC,NMIX,NGRP,NL,NDEL,LEAK,NTOT0,NTOT1,ZNUS,
     1 CHI,ZSIGF,DIFFX,DIFFY,DIFFZ,HFAC,IJJ,NJJ,SCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store modified nuclear properties in a new macrolib.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* IPMAC   pointer to create mode macrolib.
* NMIX    maximum number of material mixtures.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* NDEL    number of precursor groups for delayed neutron.
* LEAK    diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* NTOT0   flux-weighted total macroscopic x-sections.
* NTOT1   current-weighted total macroscopic x-sections.
* ZNUS    nu*fission macroscopic x-sections.
* CHI     fission spectra.
* ZSIGF   fission macroscopic x-sections.
* DIFFX   x-directed diffusion coefficients.
* DIFFY   y-directed diffusion coefficients.
* DIFFZ   z-directed diffusion coefficients.
* HFAC    h-factors (kappa*fission macroscopic x-sections).
* IJJ     highest energy number for which the scattering
*         component to group g does not vanish.
* NJJ     number of energy groups for which the scattering
*         component does not vanish.
* SCAT    scattering macroscopic x-sections.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC
      INTEGER NMIX,NGRP,NL,NDEL,LEAK,IJJ(NMIX,NL,NGRP),NJJ(NMIX,NL,NGRP)
      REAL NTOT0(NMIX,NGRP),NTOT1(NMIX,NGRP),ZSIGF(NMIX,NGRP),
     1     DIFFX(NMIX,NGRP),DIFFY(NMIX,NGRP),DIFFZ(NMIX,NGRP),
     2     ZNUS(NMIX,NGRP,NDEL+1),CHI(NMIX,NGRP,NDEL+1),HFAC(NMIX,NGRP),
     3     SCAT(NMIX,NL,NGRP,NGRP)
*----
*  LOCAL VARIABLES
*----
      CHARACTER CM*2,TEXT12*12
      TYPE(C_PTR) JPMAC,KPMAC
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NMIX*NGRP))
*----
*  STORE PROPERTIES
*----
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO 30 JGR=1,NGRP
      KPMAC=LCMGIL(JPMAC,JGR)
*     NTOT0
      CALL LCMPUT(KPMAC,'NTOT0',NMIX,2,NTOT0(1,JGR))
*     NTOT1
      CALL LCMLEN(KPMAC,'NTOT1',LENGT,ITYP)
      IF(LENGT.EQ.NMIX)
     1  CALL LCMPUT(KPMAC,'NTOT1',NMIX,2,NTOT1(1,JGR))
*     NUSIGF
      CALL LCMLEN(KPMAC,'NUSIGF',LENGT,ITYP)
      IF(LENGT.EQ.NMIX)
     1  CALL LCMPUT(KPMAC,'NUSIGF',NMIX,2,ZNUS(1,JGR,1))
      DO IDEL=1,NDEL
        WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
        CALL LCMLEN(KPMAC,TEXT12,LENGT,ITYP)
        IF(LENGT.EQ.NMIX)
     1    CALL LCMPUT(KPMAC,TEXT12,NMIX,2,ZNUS(1,JGR,IDEL+1))
      ENDDO
*     CHI
      CALL LCMLEN(KPMAC,'CHI',LENGT,ITYP)
      IF(LENGT.EQ.NMIX)
     1  CALL LCMPUT(KPMAC,'CHI',NMIX,2,CHI(1,JGR,1))
      DO IDEL=1,NDEL
        WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
        CALL LCMLEN(KPMAC,TEXT12,LENGT,ITYP)
        IF(LENGT.EQ.NMIX)
     1    CALL LCMPUT(KPMAC,TEXT12,NMIX,2,CHI(1,JGR,IDEL+1))
      ENDDO
*     NFTOT
      CALL LCMLEN(KPMAC,'NFTOT',LENGT,ITYP)
      IF(LENGT.EQ.NMIX)
     1  CALL LCMPUT(KPMAC,'NFTOT',NMIX,2,ZSIGF(1,JGR))
      IF(LEAK.EQ.1)THEN
*       DIFF
        CALL LCMPUT(KPMAC,'DIFF',NMIX,2,DIFFX(1,JGR))
      ELSEIF(LEAK.EQ.2)THEN
*       DIFFX,DIFFY,DIFFZ
        CALL LCMPUT(KPMAC,'DIFFX',NMIX,2,DIFFX(1,JGR))
        CALL LCMPUT(KPMAC,'DIFFY',NMIX,2,DIFFY(1,JGR))
        CALL LCMPUT(KPMAC,'DIFFZ',NMIX,2,DIFFZ(1,JGR))
      ENDIF
*     H-FACTOR
      CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYP)
      IF(LENGT.EQ.NMIX)
     1  CALL LCMPUT(KPMAC,'H-FACTOR',NMIX,2,HFAC(1,JGR))
      DO IL=1,NL
        WRITE (CM,'(I2.2)') IL-1
        CALL LCMLEN(KPMAC,'SCAT'//CM,LENGT,ITYP)
        IF(LENGT.NE.0)THEN
          IPOSDE=0
          DO 20 IBM=1,NMIX
            DO IGR=IJJ(IBM,IL,JGR),IJJ(IBM,IL,JGR)-
     1                NJJ(IBM,IL,JGR)+1,-1
              IPOSDE=IPOSDE+1
              WORK(IPOSDE)=SCAT(IBM,IL,IGR,JGR)
            ENDDO
            HFAC(IBM,JGR)=0.
            DO 10 IGR=1,NGRP
              HFAC(IBM,JGR)=HFAC(IBM,JGR)+SCAT(IBM,IL,JGR,IGR)
   10       CONTINUE
   20     CONTINUE
          CALL LCMPUT(KPMAC,'SCAT'//CM,IPOSDE,2,WORK)
          CALL LCMPUT(KPMAC,'NJJS'//CM,NMIX,1,NJJ(1,IL,JGR))
          CALL LCMPUT(KPMAC,'IJJS'//CM,NMIX,1,IJJ(1,IL,JGR))
          CALL LCMPUT(KPMAC,'SIGW'//CM,NMIX,2,SCAT(1,IL,JGR,JGR))
          CALL LCMPUT(KPMAC,'SIGS'//CM,NMIX,2,HFAC(1,JGR))
        ENDIF
      ENDDO
   30 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
