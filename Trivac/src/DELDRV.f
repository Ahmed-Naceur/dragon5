*DECK DELDRV
      SUBROUTINE DELDRV (IPTRK,IPSYS0,IPSYSP,IPFLU0,IPGPT,NUN,NGRP,
     1 NSTEP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the calculation of direct or adjoint sources for a fixed
* source eigenvalue problem.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS0  L_SYSTEM pointer to unperturbed system matrices.
* IPSYSP  L_SYSTEM pointer to delta system matrices.
* IPFLU0  L_FLUX pointer to the unperturbed solution.
* IPGPT   L_GPT pointer to the GPT fixed source.
* NUN     total number of unknowns per energy group.
* NGRP    number of energy groups.
* NSTEP   number of perturbation states in STEP directory.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS0,IPSYSP,IPFLU0,IPGPT
      INTEGER NUN,NGRP,NSTEP
*----
*  LOCAL VARIABLES
*----
      LOGICAL ADJ
      DOUBLE PRECISION DFLOTT
      CHARACTER TEXT4*4
      TYPE(C_PTR) JPFLU1,JPFLU2,JPGPT,KPGPT,JPSYSP,KPSYSP
      REAL, DIMENSION(:,:), ALLOCATABLE :: EVECT,ADECT,SUNKNO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(EVECT(NUN,NGRP),ADECT(NUN,NGRP),SUNKNO(NUN,NGRP))
*----
*  READ THE INPUT DATA.
*----
*     DEFAULT OPTIONS.
      IMPX=1
      ADJ=.FALSE.
*
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 20
      IF(INDIC.NE.3) CALL XABORT('DELDRV: CHARACTER DATA EXPECTED.')
*
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DELDRV: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'ADJ') THEN
         ADJ=.TRUE.
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 20
      ELSE
         CALL XABORT('DELDRV: ; EXPECTED.')
      ENDIF
      GO TO 10
*----
*  RECOVER UNPERTURBED K-EFFECTIVE AND FLUXES.
*----
   20 CALL LCMGET(IPFLU0,'K-EFFECTIVE',FKEFF)
      JPFLU1=LCMGID(IPFLU0,'FLUX')
      JPFLU2=LCMGID(IPFLU0,'AFLUX')
      DO 30 IGR=1,NGRP
      CALL LCMGDL(JPFLU1,IGR,EVECT(1,IGR))
      CALL LCMGDL(JPFLU2,IGR,ADECT(1,IGR))
   30 CONTINUE
*----
*  COMPUTE THE DIRECT OR ADJOINT FIXED SOURCES AND SAVE THE FIXED
*  SOURCES.
*----
      IF(NSTEP.EQ.0) THEN
         CALL DELPER(IPTRK,IPSYS0,IPSYSP,ADJ,NUN,NGRP,FKEFF,IMPX,EVECT,
     1   ADECT,DELKEF,SUNKNO)
         IF(ADJ) THEN
           JPGPT=LCMLID(IPGPT,'ASOUR',1)
         ELSE
           JPGPT=LCMLID(IPGPT,'DSOUR',1)
         ENDIF
         KPGPT=LCMLIL(JPGPT,1,NGRP)
         DO 40 IGR=1,NGRP
         CALL LCMPDL(KPGPT,IGR,NUN,2,SUNKNO(1,IGR))
   40    CONTINUE
      ELSE
         JPSYSP=LCMGID(IPSYSP,'STEP')
         IF(ADJ) THEN
           JPGPT=LCMLID(IPGPT,'ASOUR',NSTEP)
         ELSE
           JPGPT=LCMLID(IPGPT,'DSOUR',NSTEP)
         ENDIF
         DO 55 ISTEP=1,NSTEP
         KPSYSP=LCMGIL(JPSYSP,ISTEP)
         CALL DELPER(IPTRK,IPSYS0,KPSYSP,ADJ,NUN,NGRP,FKEFF,IMPX,EVECT,
     1   ADECT,DELKEF,SUNKNO)
         KPGPT=LCMLIL(JPGPT,ISTEP,NGRP)
         DO 50 IGR=1,NGRP
         CALL LCMPDL(KPGPT,IGR,NUN,2,SUNKNO(1,IGR))
   50    CONTINUE
   55    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(EVECT,ADECT,SUNKNO)
      RETURN
      END
