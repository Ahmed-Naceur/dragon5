*DECK NEWMGT
      SUBROUTINE NEWMGT(IPMAC,NMIX,NGRP,NL,NDEL,LEAK,NTOT0,NTOT1,ZNUS,
     1 CHI,ZSIGF,DIFFX,DIFFY,DIFFZ,HFAC,IJJ,NJJ,SCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the existing macrolib data and store them in memory.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPMAC   pointer to the macrolib information.
* NMIX    maximum number of material mixtures.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* NDEL    number of precursor groups for delayed neutron.
*
*Parameters: output
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
* LEAK    diffusion coefficient flag (=1: isotropic; =2: anisotropic).
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
      PARAMETER(IOUT=6)
      TYPE(C_PTR) JPMAC,KPMAC
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NMIX*NGRP))
*
      CALL XDRSET(WORK,NMIX*NGRP,0.)
      CALL XDRSET(HFAC,NMIX*NGRP,0.)
      CALL XDRSET(NTOT0,NMIX*NGRP,0.)
      CALL XDRSET(NTOT1,NMIX*NGRP,0.)
      CALL XDRSET(ZSIGF,NMIX*NGRP,0.)
      CALL XDRSET(DIFFX,NMIX*NGRP,0.)
      CALL XDRSET(DIFFY,NMIX*NGRP,0.)
      CALL XDRSET(DIFFZ,NMIX*NGRP,0.)
      CALL XDRSET(ZNUS,NMIX*NGRP*(NDEL+1),0.)
      CALL XDRSET(CHI,NMIX*NGRP*(NDEL+1),0.)
      CALL XDRSET(SCAT,NMIX*NL*NGRP*NGRP,0.)
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
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO 70 JGR=1,NGRP
      KPMAC=LCMGIL(JPMAC,JGR)
*     NTOT0
      CALL LCMLEN(KPMAC,'NTOT0',LENGT,ITYLCM)
      IF(LENGT.EQ.NMIX)THEN
        CALL LCMGET(KPMAC,'NTOT0',NTOT0(1,JGR))
      ELSEIF(LENGT.EQ.0)THEN
        CALL XABORT('@NEWMGT: MISSING NTOT0 DATA IN MACROLIB.')
      ELSE
        CALL XABORT('@NEWMGT: INVALID NTOT0 DATA IN MACROLIB.')
      ENDIF
*     NTOT1
      CALL LCMLEN(KPMAC,'NTOT1',LENGT,ITYLCM)
      IF(LENGT.EQ.NMIX)THEN
        CALL LCMGET(KPMAC,'NTOT1',NTOT1(1,JGR))
      ELSEIF(LENGT.NE.0)THEN
        CALL XABORT('@NEWMGT: INVALID NTOT1 DATA IN MACROLIB.')
      ENDIF
*     NUSIGF
      CALL LCMLEN(KPMAC,'NUSIGF',LENGT,ITYLCM)
      IF(LENGT.EQ.NMIX)THEN
        CALL LCMGET(KPMAC,'NUSIGF',ZNUS(1,JGR,1))
      ELSEIF(LENGT.NE.0)THEN
        CALL XABORT('@NEWMGT: INVALID NUSIGF DATA IN MACROLIB.')
      ENDIF
      DO IDEL=1,NDEL
        WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
        CALL LCMLEN(KPMAC,TEXT12,LENGT,ITYLCM)
        IF(LENGT.EQ.NMIX)THEN
          CALL LCMGET(KPMAC,TEXT12,ZNUS(1,JGR,IDEL+1))
        ELSEIF(LENGT.NE.0)THEN
          CALL XABORT('@NEWMGT: INVALID '//TEXT12//' DATA IN MACROLIB.')
        ENDIF
      ENDDO
*     CHI
      CALL LCMLEN(KPMAC,'CHI',LENGT,ITYLCM)
      IF(LENGT.EQ.NMIX)THEN
        CALL LCMGET(KPMAC,'CHI',CHI(1,JGR,1))
      ELSEIF(LENGT.NE.0)THEN
        CALL XABORT('@NEWMGT: INVALID CHI DATA IN MACROLIB.')
      ENDIF
      DO IDEL=1,NDEL
        WRITE(TEXT12,'(3HCHI,I2.2)') IDEL
        CALL LCMLEN(KPMAC,TEXT12,LENGT,ITYLCM)
        IF(LENGT.EQ.NMIX)THEN
          CALL LCMGET(KPMAC,TEXT12,CHI(1,JGR,IDEL+1))
        ELSEIF(LENGT.NE.0)THEN
          CALL XABORT('@NEWMGT: INVALID '//TEXT12//' DATA IN MACROLIB.')
        ENDIF
      ENDDO
*     NFTOT
      CALL LCMLEN(KPMAC,'NFTOT',LENGT,ITYLCM)
      IF(LENGT.EQ.NMIX)THEN
        CALL LCMGET(KPMAC,'NFTOT',ZSIGF(1,JGR))
      ELSEIF(LENGT.NE.0)THEN
        CALL XABORT('@NEWMGT: INVALID NFTOT DATA IN MACROLIB.')
      ENDIF
*     DIFF
      CALL LCMLEN(KPMAC,'DIFF',LENGT,ITYLCM)
      IF(LENGT.EQ.0)GOTO 20
      IF(LENGT.NE.NMIX)CALL XABORT('@NEWMGT: INVALID DIFF DATA.')
      CALL LCMGET(KPMAC,'DIFF',DIFFX(1,JGR))
      LEAK=1
      GOTO 30
*     DIFFX
   20 CALL LCMLEN(KPMAC,'DIFFX',LENGT,ITYLCM)
      IF(LENGT.EQ.0)GO TO 30
      IF(LENGT.NE.NMIX)CALL XABORT('@NEWMGT: INVALID DIFFX DATA.')
      CALL LCMGET(KPMAC,'DIFFX',DIFFX(1,JGR))
*     DIFFY
      CALL LCMLEN(KPMAC,'DIFFY',LENGT,ITYLCM)
      IF(LENGT.NE.NMIX)CALL XABORT('@NEWMGT: INVALID DIFFY DATA.')
      CALL LCMGET(KPMAC,'DIFFY',DIFFY(1,JGR))
*     DIFFZ
      CALL LCMLEN(KPMAC,'DIFFZ',LENGT,ITYLCM)
      IF(LENGT.NE.NMIX)CALL XABORT('@NEWMGT: INVALID DIFFZ DATA.')
      CALL LCMGET(KPMAC,'DIFFZ',DIFFZ(1,JGR))
      LEAK=2
*     H-FACTOR
   30 CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYLCM)
      IF(LENGT.EQ.NMIX)THEN
        CALL LCMGET(KPMAC,'H-FACTOR',HFAC(1,JGR))
      ELSEIF(LENGT.NE.0)THEN
        CALL XABORT('@NEWMGT: INVALID H-FACTOR DATA IN MACROLIB.')
      ENDIF
*     SCAT,NJJ,IJJ
      DO IL=1,NL
        WRITE (CM,'(I2.2)') IL-1
        CALL LCMLEN(KPMAC,'SCAT'//CM,LENGT,ITYLCM)
        IF(LENGT.GT.NMIX*NL*NGRP*NGRP)THEN
          CALL XABORT('@NEWMGT: INVALID INPUT MACROLIB(1).')
        ELSEIF(LENGT.GT.0)THEN
          CALL LCMGET(KPMAC,'SCAT'//CM,WORK)
          CALL LCMGET(KPMAC,'NJJS'//CM,NJJ(1,IL,JGR))
          CALL LCMGET(KPMAC,'IJJS'//CM,IJJ(1,IL,JGR))
          IPOSDE=0
          DO 65 IBM=1,NMIX
            IJJ0=IJJ(IBM,IL,JGR)
            DO 60 IGR=IJJ0,IJJ0-NJJ(IBM,IL,JGR)+1,-1
              IPOSDE=IPOSDE+1
              SCAT(IBM,IL,IGR,JGR)=WORK(IPOSDE)
   60       CONTINUE
   65     CONTINUE
        ELSE
          CALL XABORT('@NEWMGT: OLD FORMAT OF THE MACROLIB.')
        ENDIF
      ENDDO
   70 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
