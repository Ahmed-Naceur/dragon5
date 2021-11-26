*DECK CREINT
      SUBROUTINE CREINT(IPCPO,NISO,DERIV,NBURN,KBURN,BURN0,BURN1,NGRP,
     1   NL,IMPX,HISO,ITY,CONC,ILEAK,TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,
     2   DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover and interpolate l_compo information according to burnup and
* extracted isotope density.
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
* DERIV   =.true.: derivative of macrolib info is computed with
*         respect to burn1.
* UPS     =.true.: no upscatering cross sections will be stored.
* NBURN   number of tabulated burnup steps.
* KBURN   =0: no burnup parameters; =1: use mw day/tonne of initial
*         heavy elements).
* BURN0   user defined initial burnup.
* BURN1   user defined final burnup:
*          if burn0=burn1, a simple interpolation is performed;
*          if burn0<burn1, a time-average calculation is performed.
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* IMPX    print parameter (=0 for no print).
* HISO    hollerith name information for extracted isotopes.
* ITY     =0: do not process the isotope; =1: use number density
*         stored in conc(i); =2: use number density stored in compo.
* CONC    user defined number density.
* ILEAK
*
*Parameters: output
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
* FLUX    integrated fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER NISO,NGRP,IMPX,NBURN,KBURN,HISO(3*NISO),ITY(NISO),ILEAK
      REAL CONC(NISO),TOTAL(NGRP),ZNUG(NGRP),SNUGF(NGRP),CHI(NGRP),
     1 OVERV(NGRP),DIFFX(NGRP),DIFFY(NGRP),DIFFZ(NGRP),H(NGRP),
     2 SCAT(NL,NGRP,NGRP),FLUX(NGRP),BURN0,BURN1
      LOGICAL DERIV,UPS
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT12*12
      REAL, ALLOCATABLE, DIMENSION(:) :: BURNUP,DENSIT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(BURNUP(NBURN),DENSIT(NISO))
*----
*  CASE WITH NO BURNUP
*----
      IF(KBURN.EQ.0)THEN
        CALL LCMSIX(IPCPO,'BURN       1',1)
        CALL LCMGET(IPCPO,'ISOTOPESDENS',DENSIT)
        IF(DENSIT(1).NE.1.)CALL XABORT('@CREINT: DENSIT(1).NE.1.')
        DO I=2,NISO
          IF(ITY(I).EQ.0)THEN
            DENSIT(I)=0.
          ELSEIF(ITY(I).EQ.1)THEN
            DENSIT(I)=CONC(I)
          ELSEIF(ITY(I).NE.2)THEN
            CALL XABORT('@CREINT: INVALID VALUE OF ITY.')
          ENDIF
        ENDDO
        CALL CREMAC(IPCPO,NISO,NGRP,NL,IMPX,HISO,DENSIT,ILEAK,TOTAL,
     1       ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
        CALL LCMSIX(IPCPO,' ',2)
      ELSE
*----
*  CASE WITH BURNUP
*----
        CALL LCMGET(IPCPO,'BURNUP',BURNUP)
        TEXT12=' '
        IF(BURN0.EQ.BURN1)THEN
          DO I=1,NBURN
            IF(BURN0.EQ.BURNUP(I))THEN
              WRITE(TEXT12,'(4HBURN,4X,I4)') I
              GOTO 30
            ENDIF
          ENDDO
        ENDIF
   30   IF((TEXT12.NE.' ').AND.(.NOT.DERIV))THEN
*         BURN0=BURN1 IS A TABULATION POINT.
          CALL LCMSIX(IPCPO,TEXT12,1)
          CALL LCMGET(IPCPO,'ISOTOPESDENS',DENSIT)
          IF(DENSIT(1).NE.1.)CALL XABORT('@CREINT: DENSIT(1).NE.1.')
          DO I=2,NISO
            IF(ITY(I).EQ.0)THEN
              DENSIT(I)=0.
            ELSEIF(ITY(I).EQ.1)THEN
              DENSIT(I)=CONC(I)
            ELSEIF(ITY(I).NE.2)THEN
              CALL XABORT('@CREINT: INVALID VALUE OF ITY.')
            ENDIF
          ENDDO
          CALL CREMAC(IPCPO,NISO,NGRP,NL,IMPX,HISO,DENSIT,ILEAK,TOTAL,
     1         ZNUG,SNUGF,CHI,OVERV,DIFFX,DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
          CALL LCMSIX(IPCPO,' ',2)
        ELSE
*         INTERPOLATION IS REQUIRED.
          CALL CREBUR(IPCPO,NISO,NGRP,NL,IMPX,HISO,DERIV,NBURN,BURN0,
     1         BURN1,BURNUP,ITY,CONC,ILEAK,TOTAL,ZNUG,SNUGF,CHI,OVERV,
     2         DIFFX,DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
        ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DENSIT,BURNUP)
      RETURN
      END
