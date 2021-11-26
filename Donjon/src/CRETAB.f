*DECK CRETAB
      SUBROUTINE CRETAB(IPCPO,NISO,NGRP,NL,IMPX,HISO,NBURN,ITY,CONC,
     1 ILEAK,ZTOTAL,ZZNUG,ZNUGF,ZCHI,ZOVERV,ZDIFFX,ZDIFFY,ZDIFFZ,ZH,
     3 ZSCAT,ZFLUX,UPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create burnup dependent table with the extracted isotope.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPCPO   pointer to l_compo information.
* NISO    1+number of extracted isotopes.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* IMPX    print parameter (=0 for no print).
* HISO    hollerith name information for extracted isotopes.
* NBURN   number of tabulated burnup steps
* ITY     =0: do not process the isotope; =1: use number density
*         stored in conc(i); =2: use number density stored in compo.
* CONC    user defined number density.
*
*Parameters: output
* ILEAK   diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* ZTOTAL  burnup dependent total macroscopic x-sections
* ZZNUG   burnup dependent nu*fission macroscopic x-sections.
* ZNUGF   burnup dependent fission macroscopic x-sections.
* ZCHI    burnup dependent fission spectrum.
* ZOVERV  burnup dependent reciprocal neutron velocities.
* ZDIFFX  burnup dependent x-directed diffusion coefficients.
* ZDIFFY  burnup dependent y-directed diffusion coefficients.
* ZDIFFZ  burnup dependent z-directed diffusion coefficients.
* ZH      burnup dependent h-factors (kappa*fission macroscopic
*         x-sections).
* ZSCAT   burnup dependent scattering macroscopic x-sections.
* ZFLUX   burnup dependent integrated flux.
*
*Parameters: scratch
* TOTAL   total macroscopic x-sections.
* ZNUG    nu*fission macroscopic x-sections.
* SNUGF   fission macroscopic x-sections.
* CHI     fission spectrum.
* OVERV   reciprocal neutron velocities.
* DIFFX   x-directed diffusion coefficients.
* DIFFY   y-directed diffusion coefficients.
* DIFFZ   z-directed diffusion coefficients.
* H       h-factors (kappa*fission macroscopic x-sections).
* SCAT    scattering macroscopic x-sections.
* FLUX    integrated flux.
* DENSIT  isotopic number densities.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER NISO,NGRP,NL,IMPX,NBURN,HISO(3*NISO),ITY(NISO),ILEAK
      REAL CONC(NISO),ZTOTAL(NBURN,NGRP),ZZNUG(NBURN,NGRP),
     1 ZCHI(NBURN,NGRP),ZOVERV(NBURN,NGRP),ZDIFFX(NBURN,NGRP),
     2 ZDIFFY(NBURN,NGRP),ZDIFFZ(NBURN,NGRP),ZH(NBURN,NGRP),
     3 ZSCAT(NBURN,NL,NGRP,NGRP),ZNUGF(NBURN,NGRP),ZFLUX(NBURN,NGRP)
      LOGICAL UPS
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT12*12
      REAL, ALLOCATABLE, DIMENSION(:) :: TOTAL,ZNUG,CHI,OVERV,DIFFX,
     1 DIFFY,DIFFZ,H,SNUGF,FLUX,DENSIT
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(TOTAL(NGRP),ZNUG(NGRP),CHI(NGRP),OVERV(NGRP),DIFFX(NGRP),
     1 DIFFY(NGRP),DIFFZ(NGRP),H(NGRP),SCAT(NL,NGRP,NGRP),SNUGF(NGRP),
     2 FLUX(NGRP),DENSIT(NISO))
*----
*  RECOVER MACROSCOPIC X-SECTION INFO FROM BURNUP DIRECTORIES
*----
      DO 20 IBURN=1,NBURN
      WRITE(TEXT12,'(4HBURN,4X,I4)') IBURN
      CALL LCMSIX(IPCPO,TEXT12,1)
      CALL LCMGET(IPCPO,'ISOTOPESDENS',DENSIT)
      IF(DENSIT(1).NE.1.)CALL XABORT('@CRETAB: DENSIT(1).NE.1.')
      DO I=2,NISO
        IF(ITY(I).EQ.0)THEN
          DENSIT(I)=0.
        ELSEIF(ITY(I).EQ.1)THEN
          DENSIT(I)=CONC(I)
        ELSEIF(ITY(I).NE.2)THEN
          CALL XABORT('@CRETAB: INVALID VALUE OF ITY.')
        ENDIF
      ENDDO
      CALL CREMAC(IPCPO,NISO,NGRP,NL,IMPX,HISO,DENSIT,ILEAK,TOTAL,
     1     ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
      CALL LCMSIX(IPCPO,' ',2)
      DO 21 JGR=1,NGRP
        ZTOTAL(IBURN,JGR)=TOTAL(JGR)
        ZZNUG(IBURN,JGR)=ZNUG(JGR)
        ZNUGF(IBURN,JGR)=SNUGF(JGR)
        ZCHI(IBURN,JGR)=CHI(JGR)
        ZOVERV(IBURN,JGR)=OVERV(JGR)
        ZDIFFX(IBURN,JGR)=DIFFX(JGR)
        ZDIFFY(IBURN,JGR)=DIFFY(JGR)
        ZDIFFZ(IBURN,JGR)=DIFFZ(JGR)
        ZH(IBURN,JGR)=H(JGR)
        ZFLUX(IBURN,JGR)=FLUX(JGR)
        DO 22 IGR=1,NGRP
        DO 23 IL=1,NL
          ZSCAT(IBURN,IL,IGR,JGR)=SCAT(IL,IGR,JGR)
   23   CONTINUE
   22   CONTINUE
   21 CONTINUE
   20 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DENSIT,FLUX,SNUGF,SCAT,H,DIFFZ,DIFFY,DIFFX,OVERV,CHI,
     1 ZNUG,TOTAL)
      RETURN
      END
