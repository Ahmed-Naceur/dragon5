*DECK COMGEM
      SUBROUTINE COMGEM(IPDEPL,ITIM,TYPE,IMILI,NBURN,NBMIX,NBISO,NREAC,
     1 NVAR,VALUE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover a global parameter or a local variable from the burnup object.
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
* ITIM    index of the current burnup step.
* TYPE    type of parameter (='FLUX', 'IRRA', 'TIME', 'PUIS', 'FLUB' or
*         'MASL').
* IMILI   position of parameter (=0: global averaged value; >0: value
*         in mixture IMILI).
* NBURN   number of burnup steps in the burnup object.
* NBMIX   number of depleting mixtures.
* NBISO   number of isotopes.
* NREAC   number of depleting reactions.
* NVAR    number of depleting isotopes.
*
*Parameters: output
* VALUE   global parameter or local variable.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDEPL
      INTEGER ITIM,IMILI,NBURN,NBMIX,NBISO,NREAC,NVAR
      REAL VALUE
      CHARACTER TYPE*4
*----
*  LOCAL VARIABLES
*----
      REAL BUIR(2)
      CHARACTER CDIRO*12
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: JM
      REAL, ALLOCATABLE, DIMENSION(:) :: DEN,TIME,VX,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PARAM,VPHV
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIG
*----
*  SCRATCH STORAGE ALLOCATION
*   PARAM   parameters (PARAM(*,1): fluence; PARAM(*,2): burnup or
*           energy).
*----
      ALLOCATE(JM(NBMIX,NVAR))
      ALLOCATE(DEN(NBISO),TIME(NBURN),PARAM(NBMIX,2),VPHV(NBMIX,2),
     1 VX(NBMIX),WORK(NBMIX),SIG(NVAR+1,NREAC+1,NBMIX))
*
      CALL LCMGET(IPDEPL,'DEPL-TIMES',TIME)
      CALL LCMGET(IPDEPL,'VOLUME-MIX',VX)
      CALL LCMGET(IPDEPL,'DEPLETE-MIX',JM)
*----
*  COMPUTE THE EXPOSURE AND BURNUP
*----
      IF(IMILI.NE.0) THEN
         NB0=1
         WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') NB0
         CALL LCMSIX(IPDEPL,CDIRO,1)
         CALL LCMGET(IPDEPL,'INT-FLUX',VPHV(1,1))
         CALL LCMSIX(IPDEPL,' ',2)
         DO 10 IBM=1,NBMIX
         PARAM(IBM,1)=0.0
         PARAM(IBM,2)=0.0
   10    CONTINUE
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
   20    CONTINUE
   25    CONTINUE
      ENDIF
*
      IF(TYPE.EQ.'FLUB') THEN
         IF(IMILI.EQ.0) THEN
            WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') ITIM
            CALL LCMSIX(IPDEPL,CDIRO,1)
            CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BUIR)
            CALL LCMSIX(IPDEPL,' ',2)
            VALUE=BUIR(2)
         ELSE
            VALUE=PARAM(IMILI,1)
         ENDIF
      ELSE IF(TYPE.EQ.'IRRA') THEN
         IF(IMILI.EQ.0) THEN
            WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') ITIM
            CALL LCMSIX(IPDEPL,CDIRO,1)
            CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BUIR)
            CALL LCMSIX(IPDEPL,' ',2)
            VALUE=BUIR(1)
         ELSE
            CALL LCMGET(IPDEPL,'FUELDEN-MIX',WORK)
            IF(WORK(IMILI).EQ.0.0) THEN
               VALUE=0.0
            ELSE
               VALUE=PARAM(IMILI,2)/WORK(IMILI)
            ENDIF
         ENDIF
      ELSE IF(TYPE.EQ.'TIME') THEN
         VALUE=(TIME(ITIM)-TIME(1))*1.0E8
      ELSE IF(TYPE.EQ.'FLUX') THEN
         WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') ITIM
         CALL LCMSIX(IPDEPL,CDIRO,1)
         CALL LCMGET(IPDEPL,'INT-FLUX',PARAM(1,1))
         CALL LCMSIX(IPDEPL,' ',2)
         IF(IMILI.EQ.0) THEN
            VTOT=0.0
            VALUE=0.0
            DO 30 IBM=1,NBMIX
            VTOT=VTOT+VX(IBM)
            VALUE=VALUE+1.0E-11*PARAM(IBM,1)
   30       CONTINUE
            VALUE=VALUE/VTOT
         ELSE
            VALUE=1.0E-11*PARAM(IMILI,1)/VX(IMILI)
         ENDIF
      ELSE IF(TYPE.EQ.'PUIS') THEN
         WRITE(CDIRO,'(8HDEPL-DAT,I4.4)') ITIM
         CALL LCMSIX(IPDEPL,CDIRO,1)
         CALL LCMGET(IPDEPL,'MICRO-RATES',SIG)
         CALL LCMGET(IPDEPL,'ISOTOPESDENS',DEN)
         CALL LCMSIX(IPDEPL,' ',2)
         IF(IMILI.EQ.0) THEN
            VTOT=0.0
            VALUE=0.0
            DO 50 IBM=1,NBMIX
            VTOT=VTOT+VX(IBM)
            GAR=SIG(NVAR+1,NREAC,IBM)+SIG(NVAR+1,NREAC+1,IBM)
            DO 40 IS=1,NVAR
            IF(JM(IBM,IS).GT.0) THEN
              GAR=GAR+VX(IBM)*DEN(JM(IBM,IS))*(SIG(IS,NREAC,IBM)+
     &        SIG(IS,NREAC+1,IBM))
            ENDIF
   40       CONTINUE
            VALUE=VALUE+1.0E-8*GAR
   50       CONTINUE
            VALUE=VALUE/VTOT
         ELSE
            GAR=SIG(NVAR+1,NREAC,IMILI)+SIG(NVAR+1,NREAC+1,IMILI)
            DO 60 IS=1,NVAR
            IF(JM(IMILI,IS).GT.0) THEN
               GAR=GAR+VX(IMILI)*DEN(JM(IMILI,IS))*(SIG(IS,NREAC,IMILI)+
     &         SIG(IS,NREAC+1,IMILI))
            ENDIF
   60       CONTINUE
            VALUE=1.0E-8*GAR/VX(IMILI)
         ENDIF
      ELSE IF(TYPE.EQ.'MASL') THEN
         CALL LCMGET(IPDEPL,'FUELDEN-MIX',WORK)
         IF(IMILI.EQ.0) THEN
            VTOT=0.0
            VALUE=0.0
            DO 70 IBM=1,NBMIX
            VTOT=VTOT+VX(IBM)
            VALUE=VALUE+WORK(IBM)
   70       CONTINUE
            VALUE=VALUE/VTOT
         ELSE
            VALUE=WORK(IMILI)/VX(IMILI)
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SIG,WORK,VX,VPHV,PARAM,TIME,DEN)
      DEALLOCATE(JM)
      RETURN
      END
