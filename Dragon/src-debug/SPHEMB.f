*DECK SPHEMB
      SUBROUTINE SPHEMB(IPLIB,IPCPO,NGRP,NMIX,MIXUPD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build embedded macrolib and recover depletion data from the
* multicompo.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, R. Chambon
*
*Parameters: input/output
* IPLIB   address of the microlib LCM object.
* IPCPO   address of the multicompo object.
* NGRP    number of energy groups.
* NMIX    maximum number of material mixtures in the microlib.
* MIXUPD  tag for mixture which will be updated.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPCPO
      INTEGER NGRP,NMIX
      INTEGER MIXUPD(NMIX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,IOUT=6)
      CHARACTER HSMG*131
      INTEGER ISTATE(NSTATE),IST1(NSTATE),IST2(NSTATE)
      REAL TMPDAY(3)
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKL
      INTEGER, POINTER, DIMENSION(:) :: ISONA,ISOMI
      REAL, POINTER, DIMENSION(:) :: DENIS
      TYPE(C_PTR) ISONA_PTR,ISOMI_PTR,DENIS_PTR
*----
*  RECOVER THE DEPLETION CHAIN
*----
      CALL LCMLEN(IPLIB,'DEPL-CHAIN',ILENG1,ITYLCM)
      CALL LCMLEN(IPCPO,'DEPL-CHAIN',ILENG2,ITYLCM)
      IF((ILENG1.NE.0).AND.(ILENG2.NE.0)) THEN
        CALL LCMSIX(IPCPO,'DEPL-CHAIN',1)
        CALL LCMGET(IPCPO,'STATE-VECTOR',IST1)
        CALL LCMSIX(IPCPO,' ',2)
        CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
        CALL LCMGET(IPLIB,'STATE-VECTOR',IST2)
        CALL LCMSIX(IPLIB,' ',2)
        DO 100 I=1,NSTATE
        IF(IST1(I).NE.IST2(I)) THEN
          WRITE(HSMG,'(40HSPHEMB: INVALID STATE-VECTOR COMPONENT (,I2,
     1    36H) FOR DEPL-CHAIN DATA IN MULTICOMPO ,1H.)') I
          CALL XABORT(HSMG)
        ENDIF
  100   CONTINUE
      ELSE IF((ILENG1.EQ.0).AND.(ILENG2.NE.0)) THEN
        CALL LCMSIX(IPCPO,'DEPL-CHAIN',1)
        CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
        CALL LCMEQU(IPCPO,IPLIB)
        CALL LCMSIX(IPLIB,' ',2)
        CALL LCMSIX(IPCPO,' ',2)
      ENDIF
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      MAXMIX=ISTATE(1)
      IF(MAXMIX.NE.NMIX) CALL XABORT('SPHEMB: INVALID NMIX.')
      NBISO=ISTATE(2)
      ALLOCATE(MASK(MAXMIX),MASKL(NGRP))
      CALL LCMGPD(IPLIB,'ISOTOPESUSED',ISONA_PTR)
      CALL LCMGPD(IPLIB,'ISOTOPESMIX',ISOMI_PTR)
      CALL LCMGPD(IPLIB,'ISOTOPESDENS',DENIS_PTR)
      CALL C_F_POINTER(ISONA_PTR,ISONA,(/ NBISO /))
      CALL C_F_POINTER(ISOMI_PTR,ISOMI,(/ NBISO /))
      CALL C_F_POINTER(DENIS_PTR,DENIS,(/ NBISO /))
      MASK(:MAXMIX)=.FALSE.
      MASKL(:NGRP)=.TRUE.
      DO 110 ISOT=1,NBISO
        IBM=ISOMI(ISOT)
        IF(IBM.GT.0) THEN
          IF(MIXUPD(IBM).NE.0) MASK(IBM)=.TRUE.
        ENDIF
  110 CONTINUE
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL LIBMIX(IPLIB,MAXMIX,NGRP,NBISO,ISONA,ISOMI,DENIS,MASK,MASKL,
     1 ITSTMP,TMPDAY)
      DEALLOCATE(MASKL,MASK)
      RETURN
      END
