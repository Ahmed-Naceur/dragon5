*DECK COMGEN
      SUBROUTINE COMGEN(IPDEPL,IPEDIT,NREG,NMIL,ITIM,TYPE,NBURN,NBMIX,
     1 NBISO,NREAC,NVAR,ILOC,NLOC,RVALOC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover a local variables from the burnup object and homogenize them
* on the output mixtures.
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
* IPDEPL  pointer to the burnup object.
* IPEDIT  pointer to the edition object.
* NREG    number of volumes in the depleting geometry.
* NMIL    number of homogenized output mixtures.
* ITIM    index of the current burnup step.
* TYPE    type of parameter (='FLUX', 'IRRA', 'PUIS', 'FLUG', 'FLUB' or
*        'MASL').
* NBURN   number of burnup steps in the burnup object.
* NBMIX   number of depleting mixtures.
* NBISO   number of isotopes.
* NREAC   number of depleting reactions.
* NVAR    number of depleting isotopes.
* ILOC    position of local parameter in RVALOC.
* NLOC    first dimension of matrix RVALOC.
*
*Parameters: output
* RVALOC  local variable values in homogeneous mixtures.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDEPL,IPEDIT
      INTEGER NREG,NMIL,ITIM,NBURN,NBMIX,NBISO,NREAC,NVAR,ILOC,NLOC
      REAL RVALOC(NLOC,NMIL)
      CHARACTER TYPE*4
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER CDIRO*12
      INTEGER IPAR(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATR,MERG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: JM
      REAL, ALLOCATABLE, DIMENSION(:) :: DEN,TIME,VX,WORK,VOLR,VOLIBM
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PARAM,VPHV
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIG
*----
*  SCRATCH STORAGE ALLOCATION
*   PARAM   parameters (PARAM(*,1): fluence; PARAM(*,2): burnup or
*           energy).
*----
      ALLOCATE(JM(NBMIX,NVAR),MATR(NREG),MERG(NREG))
      ALLOCATE(DEN(NBISO),TIME(NBURN),PARAM(NBMIX,2),VPHV(NBMIX,2),
     1 VX(NBMIX),WORK(NBMIX),SIG(NVAR+1,NREAC+1,NBMIX),VOLR(NREG),
     2 VOLIBM(NMIL))
*
      CALL LCMGET(IPDEPL,'DEPL-TIMES',TIME)
      CALL LCMGET(IPDEPL,'VOLUME-MIX',VX)
      CALL LCMGET(IPDEPL,'DEPLETE-MIX',JM)
*----
*  COMPUTE THE EXPOSURE AND BURNUP
*----
      NB0=1
      WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') NB0
      CALL LCMSIX(IPDEPL,CDIRO,1)
      CALL LCMGET(IPDEPL,'INT-FLUX',VPHV(1,1))
      CALL LCMSIX(IPDEPL,' ',2)
      DO 10 IBM=1,NBMIX
      PARAM(IBM,1)=0.0
      PARAM(IBM,2)=0.0
   10 CONTINUE
      DO 25 NB=NB0+1,ITIM
      WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') NB
      CALL LCMSIX(IPDEPL,CDIRO,1)
      CALL LCMGET(IPDEPL,'INT-FLUX',VPHV(1,2))
      CALL LCMGET(IPDEPL,'ENERG-MIX',WORK)
      CALL LCMSIX(IPDEPL,' ',2)
      DO 20 IBM=1,NBMIX
      PHIAV=0.5*(VPHV(IBM,1)+VPHV(IBM,2))/VX(IBM)
      PARAM(IBM,1)=PARAM(IBM,1)+PHIAV*(TIME(NB)-TIME(NB-1))
      PARAM(IBM,2)=PARAM(IBM,2)+WORK(IBM)/8.64E-4
      VPHV(IBM,1)=VPHV(IBM,2)
   20 CONTINUE
   25 CONTINUE
*----
*  RECOVER HOMOGENIZATION INFORMATION FROM THE EDITION OBJECT
*----
      CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
      IF(NMIL.NE.IPAR(1)) CALL XABORT('COMGEN: INVALID NMIL.')
      CALL LCMGET(IPEDIT,'REF:VOLUME',VOLR)
      CALL LCMGET(IPEDIT,'REF:MATCOD',MATR)
      CALL LCMGET(IPEDIT,'REF:IMERGE',MERG)
*
      DO 30 IBM=1,NMIL
      VOLIBM(IBM)=0.0
      RVALOC(ILOC,IBM)=0.0
   30 CONTINUE
      DO 50 IREG=1,NREG
      IBM=MERG(IREG)
      IMILI=MATR(IREG)
      VV=VOLR(IREG)
      IF(TYPE.EQ.'FLUG') THEN
*        N/KB IN GLOBAL HOMOGENIZED MIXTURE
         RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)+VV*PARAM(IMILI,1)
         VOLIBM(IBM)=VOLIBM(IBM)+VV
      ELSE IF(TYPE.EQ.'FLUB') THEN
*        N/KB IN FUEL ONLY
         DO 35 IS=1,NVAR
         IF(JM(IMILI,IS).GT.0) THEN
            RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)+VV*PARAM(IMILI,1)
            VOLIBM(IBM)=VOLIBM(IBM)+VV
            GO TO 50
         ENDIF
   35    CONTINUE
      ELSE IF(TYPE.EQ.'IRRA') THEN
*        MWD/TONNE
         CALL LCMGET(IPDEPL,'FUELDEN-MIX',WORK)
         IF(WORK(IMILI).NE.0.0) THEN
            RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)+PARAM(IMILI,2)
            VOLIBM(IBM)=VOLIBM(IBM)+WORK(IMILI)
         ENDIF
      ELSE IF(TYPE.EQ.'FLUX') THEN
         WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') ITIM
         CALL LCMSIX(IPDEPL,CDIRO,1)
         CALL LCMGET(IPDEPL,'INT-FLUX',PARAM(1,1))
         CALL LCMSIX(IPDEPL,' ',2)
         RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)+VV*1.0E-11*PARAM(IMILI,1)/
     1   VX(IMILI)
         VOLIBM(IBM)=VOLIBM(IBM)+VV
      ELSE IF(TYPE.EQ.'PUIS') THEN
         WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') ITIM
         CALL LCMSIX(IPDEPL,CDIRO,1)
         CALL LCMGET(IPDEPL,'MICRO-RATES',SIG)
         CALL LCMGET(IPDEPL,'ISOTOPESDENS',DEN)
         CALL LCMSIX(IPDEPL,' ',2)
         GAR=SIG(NVAR+1,NREAC,IMILI)+SIG(NVAR+1,NREAC+1,IMILI)
         DO 40 IS=1,NVAR
         IF(JM(IMILI,IS).GT.0) THEN
            GAR=GAR+VX(IMILI)*DEN(JM(IMILI,IS))*(SIG(IS,NREAC,IMILI)+
     &      SIG(IS,NREAC+1,IMILI))
         ENDIF
   40    CONTINUE
         RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)+VV*1.0E-8*GAR/VX(IMILI)
         VOLIBM(IBM)=VOLIBM(IBM)+VV
      ELSE IF(TYPE.EQ.'MASL') THEN
         CALL LCMGET(IPDEPL,'FUELDEN-MIX',WORK)
         IF(WORK(IMILI).GT.0.0) RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)+
     1   VV*WORK(IMILI)/VX(IMILI)
         VOLIBM(IBM)=VOLIBM(IBM)+VV
      ENDIF
   50 CONTINUE
      DO 60 IBM=1,NMIL
      IF(VOLIBM(IBM).NE.0.0) THEN
         RVALOC(ILOC,IBM)=RVALOC(ILOC,IBM)/VOLIBM(IBM)
      ENDIF
   60 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VOLIBM,VOLR,SIG,WORK,VX,VPHV,PARAM,TIME,DEN)
      DEALLOCATE(MERG,MATR,JM)
      RETURN
      END
