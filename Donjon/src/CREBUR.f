*DECK CREBUR
      SUBROUTINE CREBUR(IPCPO,NISO,NGRP,NL,IMPX,HISO,DERIV,NBURN,BURN0,
     1 BURN1,BURNUP,ITY,CONC,ILEAK,TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,
     2 DIFFY,DIFFZ,H,SCAT,FLUX,UPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate l_compo for a given burnup value.
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
* IMPX    printing index (=0 for no print).
* HISO    hollerith name information for extracted isotopes.
* DERIV   =.true.: derivative of macrolib info is computed with
*         respect to burn1.
* NBURN   number of tabulated burnup steps.
* BURN0   user defined initial burnup.
* BURN1   user defined final burnup:
*          if burn0=burn1, a simple interpolation is performed;
*          if burn0<burn1, a time-average calculation is performed.
* BURNUP  burnup tabulation points.
* ITY     =0: do not process the isotope; =1: use number density
*         stored in conc(i); =2: use number density stored in compo.
* CONC    user defined number density.
*
*Parameters: output
* ILEAK   diffusion coefficient flag (=1: isotropic; =2: anisotropic).
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
* UPS
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER NISO,NGRP,NL,IMPX,NBURN,HISO(3*NISO),ITY(NISO),ILEAK
      REAL BURNUP(NBURN),CONC(NISO),TOTAL(NGRP),ZNUG(NGRP),
     1 CHI(NGRP),OVERV(NGRP),DIFFX(NGRP),DIFFY(NGRP),DIFFZ(NGRP),
     2 H(NGRP),SCAT(NL,NGRP,NGRP),SNUGF(NGRP),FLUX(NGRP),BURN0,BURN1
      LOGICAL DERIV,UPS
*----
*  LOCAL VARIABLES
*----
      LOGICAL LCUBIC
      PARAMETER(LCUBIC=.TRUE.)
      REAL, ALLOCATABLE, DIMENSION(:) :: TERP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ZTOTAL,ZZNUG,ZNUGF,ZCHI,
     1 ZOVERV,ZDIFFX,ZDIFFY,ZDIFFZ,ZH,ZFLUX
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ZSCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ZTOTAL(NBURN,NGRP),ZZNUG(NBURN,NGRP),ZNUGF(NBURN,NGRP),
     1 ZCHI(NBURN,NGRP),ZOVERV(NBURN,NGRP),ZDIFFX(NBURN,NGRP),
     2 ZDIFFY(NBURN,NGRP),ZDIFFZ(NBURN,NGRP),ZH(NBURN,NGRP),
     3 ZFLUX(NBURN,NGRP),ZSCAT(NBURN,NL,NGRP,NGRP),TERP(NBURN))
*----
*  RECOVER MACROSCOPIC X-SECTION INFO FROM BURNUP DIRECTORIES
*----
      IF(NBURN.LE.1)CALL XABORT('@CREBUR: NO BURNUP INFORMATION.')
      CALL CRETAB(IPCPO,NISO,NGRP,NL,IMPX,HISO,NBURN,ITY,CONC,ILEAK,
     1 ZTOTAL,ZZNUG,ZNUGF,ZCHI,ZOVERV,ZDIFFX,ZDIFFY,ZDIFFZ,ZH,ZSCAT,
     2 ZFLUX,UPS)
*----
*  PERFORM INTERPOLATION OR TIME AVERAGING
*----
      IF(BURN0.LT.BURN1)THEN
*       TIME-AVERAGED
        CALL ALTERI(LCUBIC,NBURN,BURNUP,BURN0,BURN1,TERP)
        DO 100 I=1,NBURN
        TERP(I)=TERP(I)/(BURN1-BURN0)
  100   CONTINUE
      ELSE IF(BURN0.EQ.BURN1)THEN
*       INSTANTANEOUS
        CALL ALTERP(LCUBIC,NBURN,BURNUP,BURN0,DERIV,TERP)
      ELSE
        CALL XABORT('@CREBUR: ILLEGAL BURN1 VALUE.')
      ENDIF
      CALL CREITP(NGRP,NL,NBURN,TERP,TOTAL,ZNUG,SNUGF,CHI,OVERV,DIFFX,
     1 DIFFY,DIFFZ,H,SCAT,FLUX,ZTOTAL,ZZNUG,ZNUGF,ZCHI,ZOVERV,ZDIFFX,
     2 ZDIFFY,ZDIFFZ,ZH,ZSCAT,ZFLUX)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(TERP,ZFLUX,ZSCAT,ZH,ZDIFFZ,ZDIFFY,ZDIFFX,ZOVERV,ZCHI,
     1 ZNUGF,ZZNUG,ZTOTAL)
      RETURN
      END
