*DECK FLDNOR
      SUBROUTINE FLDNOR(IPSYS,NUN,NGRP,NEL,NBMIX,MAT,VOL,IDL,HTYPE,
     1 EVECT,FNORM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Normalization of the neutron flux
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSYS   L_SYSTEM pointer to the system matrices.
* NUN     number of flux unknowns per energy group.
* NGRP    number of energy groups.
* NEL     number of finite elements.
* NBMIX   number of material mixtures.
* MAT     material mixture indices per finite element.
* VOL     volumes of the finite elements.
* IDL     position of averaged flux in neutron flux unknowns.
* HTYPE   type of flux: 'DIRE' or 'ADJO'
* EVECT   neutron flux unknowns.
*
*Parameters: output
* EVECT   normalized flux.
* FNORM   normalization factor.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS
      INTEGER NUN,NGRP,NEL,NBMIX,MAT(NEL),IDL(NEL)
      CHARACTER*4 HTYPE
      REAL VOL(NEL),EVECT(NUN,NGRP),FNORM
*----
*  LOCAL VARIABLES
*----
      CHARACTER*12 TEXT12
      REAL, DIMENSION(:), ALLOCATABLE :: SGD
*----
*  COMPUTE THE POWER INTEGRAL
*----
      POWER=0.0
      IF(HTYPE.EQ.'DIRE') THEN
         ALLOCATE(SGD(NBMIX))
         DO 25 IGR=1,NGRP
         DO 20 JGR=1,NGRP
         WRITE(TEXT12,'(4HFISS,2I3.3)') IGR,JGR
         CALL LCMLEN(IPSYS,TEXT12,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 20
         IF(ILONG.GT.NBMIX) CALL XABORT('FLDNOR: NBMIX OVERFLOW.')
         CALL LCMGET(IPSYS,TEXT12,SGD)
         DO 10 IEL=1,NEL
         IND=IDL(IEL)
         IBM=MAT(IEL)
         IF(IND.EQ.0) GO TO 10
         POWER=POWER+VOL(IEL)*SGD(IBM)*EVECT(IND,JGR)
   10    CONTINUE
   20    CONTINUE
   25    CONTINUE
         DEALLOCATE(SGD)
      ELSE IF(HTYPE.EQ.'ADJO') THEN
         DO 35 JGR=1,NGRP
         DO 30 IEL=1,NEL
         IND=IDL(IEL)
         IF(IND.EQ.0) GO TO 30
         POWER=POWER+VOL(IEL)*EVECT(IND,JGR)
   30    CONTINUE
   35    CONTINUE
      ENDIF
      IF(POWER.EQ.0.0) CALL XABORT('FLDNOR: UNABLE TO NORMALIZE.')
      FNORM=1.0/POWER
*----
*  NORMALIZE THE FLUX
*----
      DO 45 IND=1,NUN
      DO 40 IGR=1,NGRP
      EVECT(IND,IGR)=EVECT(IND,IGR)*FNORM
   40 CONTINUE
   45 CONTINUE
      RETURN
      END
