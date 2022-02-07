*DECK TRIALB
      SUBROUTINE TRIALB(IPTRK,IPMACR,IPMACP,IPSYS,NGRP,NALBP,IPR,GAMMA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Process physical albedo information and calculation of multigroup
* albedo functions.
*
*Copyright:
* Copyright (C) 2018 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   L_TRACK pointer to the TRIVAC tracking information.
* IPMACR  L_MACROLIB pointer to the unperturbed cross sections.
* IPMACP  L_MACROLIB pointer to the perturbed cross sections if
*         IPR.gt.0. Equal to IPMACR if IPR=0.
* IPSYS   L_SYSTEM pointer to system matrices.
* NGRP    number of energy groups.
* NALBP   number of physical albedos per energy group.
* IPR     type of assembly:
*         =0: calculation of the system matrices;
*         =1: calculation of the derivative of these matrices;
*         =2: calculation of the first variation of these matrices;
*         =3: identical to IPR=2, but these variation are added to
*         unperturbed system matrices.
*
*Parameters: output
* GAMMA   albedo functions
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMACR,IPMACP,IPSYS
      INTEGER NGRP,NALBP,IPR
      REAL GAMMA(NALBP,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER TEXT12*12
      REAL, DIMENSION(:,:), ALLOCATABLE :: ALBP,DALBP
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ALBP(NALBP,NGRP),DALBP(NALBP,NGRP))
*----
*  RECOVER PHYSICAL ALBEDOS
*----
      IF(NALBP.EQ.0) CALL XABORT('TRIALB: NO PHYSICAL ALBEDOS.')
      CALL LCMGET(IPMACR,'ALBEDO',ALBP)
      IF(IPR.GT.0) CALL LCMGET(IPMACP,'ALBEDO',DALBP)
*----
*  COMPUTE ALBEDO FUNCTIONS
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ICHX=ISTATE(12)
      DO IGR=1,NGRP
        CALL XDRSET(GAMMA(1,IGR),NALBP,0.0)
        DO IALB=1,NALBP
          IF(ICHX.NE.2) THEN
            IF(IPR.EQ.0) THEN
              GAMMA(IALB,IGR)=ALB(ALBP(IALB,IGR))
            ELSE
              GAMMA(IALB,IGR)=ALB(DALBP(IALB,IGR))
            ENDIF
          ELSE IF((ICHX.EQ.2).AND.(ALBP(IALB,IGR).NE.1.0)) THEN
            IF(IPR.EQ.0) THEN
              GAMMA(IALB,IGR)=1.0/ALB(ALBP(IALB,IGR))
            ELSE IF(IPR.EQ.1) THEN
              GG=ALB(ALBP(IALB,IGR))
              DGG=ALB(DALBP(IALB,IGR))
              GAMMA(IALB,IGR)=-DGG/(GG**2)
            ELSE
              GG=ALB(ALBP(IALB,IGR))
              DGG=ALB(ALBP(IALB,IGR))+ALB(DALBP(IALB,IGR))
              GAMMA(IALB,IGR)=1.0/DGG-1.0/GG
            ENDIF
          ELSE IF((ICHX.EQ.2).AND.(ALBP(IALB,IGR).EQ.1.0)) THEN
            GAMMA(IALB,IGR)=1.0E20
          ENDIF
        ENDDO
*----
*  SAVE ALBEDO FUNCTIONS ON IPSYS
*----
         WRITE(TEXT12,'(9HALBEDO-FU,I3.3)') IGR
         CALL LCMPUT(IPSYS,TEXT12,NALBP,2,GAMMA(1,IGR))
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DALBP,ALBP)
      RETURN
      END
