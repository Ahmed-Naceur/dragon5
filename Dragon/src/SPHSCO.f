*DECK SPHSCO
      SUBROUTINE SPHSCO(IPCPO,ICAL,IMPX,IMC,NMIL,NGRP,SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Apply a new set of SPH factors for an elementary calculation in a
* Multicompo.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPCPO   pointer to the Multicompo (L_MULTICOMPO signature).
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* IMC     type of macro-calculation (=1 diffusion or SPN; 
*         =2 other options).
* NMIL    number of mixtures in the elementary calculation.
* NGRP    number of energy groups in the elementary calculation.
* SPH     SPH-factor set to be applied to the Multicompo.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPCPO
      INTEGER ICAL,IMPX,IMC,NMIL,NGRP
      REAL SPH(NMIL,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,IOUT=6)
      INTEGER ISTATE(NSTATE)
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO,MPCPO
      REAL, ALLOCATABLE, DIMENSION(:) :: SPH2
*
      CALL LCMLEN(IPCPO,'STATE-VECTOR',ILENG,ITYLCM)
      IF(ILENG.EQ.0) CALL XABORT('SPHSCO: INVALID MULTICOMPO.')
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      IF(NMIL.NE.ISTATE(1)) THEN
        CALL XABORT('SPHSCO: INVALID NUMBER OF MIXTURES(1).')
      ELSE IF(NGRP.NE.ISTATE(2)) THEN
        CALL XABORT('SPHSCO: INVALID NUMBER OF ENERGY GROUPS(1).')
      ELSE IF((ICAL.LE.0).OR.(ICAL.GT.ISTATE(3))) THEN
        CALL XABORT('SPHSCO: INVALID VALUE OF ICAL.')
      ENDIF
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      DO 20 IBM=1,NMIL
      IF(IMPX.GT.0) WRITE(IOUT,'(/33H SPHSCO: PROCESS MULTICOMPO MIXTU,
     1 2HRE,I5)') IBM
      KPCPO=LCMGIL(JPCPO,IBM)
      LPCPO=LCMGID(KPCPO,'CALCULATIONS')
      MPCPO=LCMGIL(LPCPO,ICAL)
      CALL LCMGET(MPCPO,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.1) THEN
        CALL XABORT('SPHSCO: INVALID NUMBER OF MIXTURES(2).')
      ELSE IF(ISTATE(3).NE.NGRP) THEN
        CALL XABORT('SPHSCO: INVALID NUMBER OF ENERGY GROUPS(2).')
      ENDIF
      NISOT=ISTATE(2)
      NL=ISTATE(4)
      NED=ISTATE(13)
      NDEL=ISTATE(19)
      NW=ISTATE(25)
      ALLOCATE(SPH2(NGRP))
      DO 10 IGR=1,NGRP
      SPH2(IGR)=SPH(IBM,IGR)
   10 CONTINUE
      NALBP=0 ! no albedo correction
      CALL SPHCMI(MPCPO,IMPX,IMC,1,NISOT,NGRP,NL,NW,NED,NDEL,NALBP,SPH2)
      DEALLOCATE(SPH2)
   20 CONTINUE
      RETURN
      END
