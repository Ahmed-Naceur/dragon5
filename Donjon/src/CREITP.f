*DECK CREITP
      SUBROUTINE CREITP(NGRP,NL,NBURN,TERP,TOTAL,ZNUG,SNUGF,CHI,
     1    OVERV,DIFFX,DIFFY,DIFFZ,H,SCAT,FLUX,ZTOTAL,ZZNUG,ZNUGF,
     2    ZCHI,ZOVERV,ZDIFFX,ZDIFFY,ZDIFFZ,ZH,ZSCAT,ZFLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate burnup dependent table for a given burnup value or
* time-average or derivatives.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* NGRP    number of energy groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* NBURN   number of tabulated burnup steps.
* TERP    interpolation factors.
* ZTOTAL  burnup dependent total macroscopic x-sections
* ZZNUG   burnup dependent nu*fission macroscopic x-sections.
* ZNUGF   burnup dependent fission macroscopic x-sections.
* ZCHI    burnup dependent fission spectrum.
* ZOVERV  burnup dependent reciprocal neutron velocities.
* ZDIFFX  burnup dependent x-directed diffusion coefficients.
* ZDIFFY  burnup dependent y-directed diffusion coefficients.
* ZDIFFZ  burnup dependent z-directed diffusion coefficients.
* ZH      burnup dependent h-factors.
* ZSCAT   burnup dependent scattering macroscopic x-sections.
* ZFLUX   burnup dependent integrated flux.
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
* FLUX    integrated flux.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,NL,NBURN
      REAL TERP(NBURN),TOTAL(NGRP),ZNUG(NGRP),CHI(NGRP),OVERV(NGRP),
     1 DIFFX(NGRP),DIFFY(NGRP),DIFFZ(NGRP),H(NGRP),SCAT(NL,NGRP,NGRP),
     2 SNUGF(NGRP),FLUX(NGRP),ZTOTAL(NBURN,NGRP),ZZNUG(NBURN,NGRP),
     3 ZCHI(NBURN,NGRP),ZOVERV(NBURN,NGRP),ZDIFFX(NBURN,NGRP),
     4 ZDIFFY(NBURN,NGRP),ZDIFFZ(NBURN,NGRP),ZH(NBURN,NGRP),
     5 ZSCAT(NBURN,NL,NGRP,NGRP),ZFLUX(NBURN,NGRP),ZNUGF(NBURN,NGRP)
*----
*  PERFORM INTERPOLATION OR TIME AVERAGING
*----
      CALL XDRSET(TOTAL,NGRP,0.0)
      CALL XDRSET(ZNUG,NGRP,0.0)
      CALL XDRSET(SNUGF,NGRP,0.0)
      CALL XDRSET(CHI,NGRP,0.0)
      CALL XDRSET(OVERV,NGRP,0.0)
      CALL XDRSET(DIFFX,NGRP,0.0)
      CALL XDRSET(DIFFY,NGRP,0.0)
      CALL XDRSET(DIFFZ,NGRP,0.0)
      CALL XDRSET(H,NGRP,0.0)
      CALL XDRSET(FLUX,NGRP,0.0)
      CALL XDRSET(SCAT,NGRP*NGRP*NL,0.0)
      DO 100 IBURN=1,NBURN
      WEIGHT=TERP(IBURN)
      IF(WEIGHT.EQ.0.0) GO TO 100
      DO 92 JGR=1,NGRP
      TOTAL(JGR)=TOTAL(JGR)+WEIGHT*ZTOTAL(IBURN,JGR)
      ZNUG(JGR)=ZNUG(JGR)+WEIGHT*ZZNUG(IBURN,JGR)
      SNUGF(JGR)=SNUGF(JGR)+WEIGHT*ZNUGF(IBURN,JGR)
      CHI(JGR)=CHI(JGR)+WEIGHT*ZCHI(IBURN,JGR)
      OVERV(JGR)=OVERV(JGR)+WEIGHT*ZOVERV(IBURN,JGR)
      DIFFX(JGR)=DIFFX(JGR)+WEIGHT*ZDIFFX(IBURN,JGR)
      DIFFY(JGR)=DIFFY(JGR)+WEIGHT*ZDIFFY(IBURN,JGR)
      DIFFZ(JGR)=DIFFZ(JGR)+WEIGHT*ZDIFFZ(IBURN,JGR)
      H(JGR)=H(JGR)+WEIGHT*ZH(IBURN,JGR)
      FLUX(JGR)=FLUX(JGR)+WEIGHT*ZFLUX(IBURN,JGR)
      DO 91 IGR=1,NGRP
      DO 90 IL=1,NL
      SCAT(IL,IGR,JGR)=SCAT(IL,IGR,JGR)+WEIGHT*ZSCAT(IBURN,IL,IGR,JGR)
   90 CONTINUE
   91 CONTINUE
   92 CONTINUE
  100 CONTINUE
      RETURN
      END
