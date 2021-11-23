*DECK FLPNRM
      SUBROUTINE FLPNRM(IPFLX,IPKIN,IPTRK,NMIX,NGRP,NEL,NUN,EVECT,FLUX,
     1 MAT,VOL,IDL,HFAC,PTOT,ZNRM,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover element-ordered fluxes associated with each mesh-splitted
* volume over the whole reactor core, normalize fluxes to a given
* total reactor power.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPFLX  pointer to flux information.
* IPKIN  pointer to kinetics information.
* IPTRK  pointer to tracking information.
* NMIX   maximum number of material mixtures.
* NGRP   total number of energy groups.
* NEL    total number of finite elements.
* NUN    total number of unknowns per group.
* HFAC   h-factors over the reactor core.
* PTOT   given total reactor power in watts.
* IMPX   printing index (=0 for no print).
*
*Parameters: output
* FLUX   normalized fluxes associated with each volume.
* MAT    index-number of mixture assigned to each volume.
* VOL    element-ordered mesh-splitted volumes.
* ZNRM   flux normalization factor.
*
*Parameters: scratch
* EVECT
* IDL
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLX,IPKIN,IPTRK
      INTEGER NUN,NEL,NGRP,NMIX,IMPX,IDL(NEL),MAT(NEL)
      REAL FLUX(NEL,NGRP),EVECT(NUN,NGRP),HFAC(NMIX,NGRP),VOL(NEL)
      DOUBLE PRECISION ZNRM,PTOT
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      TYPE(C_PTR) JPFLX
*----
*  RECOVER INFORMATION
*----
      CALL XDRSET(EVECT,NUN*NGRP,0.)
      IF(C_ASSOCIATED(IPFLX)) THEN
*       L_FLUX object
        JPFLX=LCMGID(IPFLX,'FLUX')
        DO 10 JGR=1,NGRP
        CALL LCMGDL(JPFLX,JGR,EVECT(1,JGR))
   10   CONTINUE
      ELSE IF(C_ASSOCIATED(IPKIN)) THEN
*       L_KINET object
        CALL LCMGET(IPKIN,'E-VECTOR',EVECT)
      ENDIF
*
      CALL XDISET(MAT,NEL,0)
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL XDISET(IDL,NEL,0)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      CALL XDRSET(VOL,NEL,0.)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
*----
*  FLUX NORMALIZATION
*----
      ZNRM=0.0D0
      IF(IMPX.GT.0)WRITE(IOUT,1002)
      CALL XDRSET(FLUX,NEL*NGRP,0.)
      DO 25 JGR=1,NGRP
      DO 20 IEL=1,NEL
      IF(MAT(IEL).EQ.0)GOTO 20
      FLUX(IEL,JGR)=EVECT(IDL(IEL),JGR)
      ZNRM=ZNRM+HFAC(MAT(IEL),JGR)*FLUX(IEL,JGR)*VOL(IEL)
   20 CONTINUE
   25 CONTINUE
      ZNRM=PTOT/ZNRM
      IF(IMPX.GT.0)WRITE(IOUT,1000)PTOT,ZNRM
      DO 35 JGR=1,NGRP
      DO 30 IEL=1,NEL
      FLUX(IEL,JGR)=FLUX(IEL,JGR)*REAL(ZNRM)
   30 CONTINUE
   35 CONTINUE
      RETURN
*
 1000 FORMAT(/1X,'GIVEN TOTAL REACTOR POWER =>',1P,E15.8,1X,
     1 'WATTS'/1X,'FLUX NORMALIZATION',1X,'FACTOR =>',1P,E15.8/)
 1002 FORMAT(/1X,'** NORMALIZING FLUXES TO A GIVEN',
     1 1X,'REACTOR POWER **')
      END
