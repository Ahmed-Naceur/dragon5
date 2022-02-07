*DECK FLDTRI
      SUBROUTINE FLDTRI(IPTRK,NEL,NUN,EVECT,MAT,VOL,IDL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the averaged flux in TRIVAC.
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
* IPTRK   L_TRACK pointer to the trivac tracking information.
* NEL     total number of finite elements.
* NUN     total number of unknown per energy group.
* EVECT   variational coefficients of the flux (contained in position
*         EVECT(1) to EVECT(LL4)).
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* IDL     position of the average flux component associated with each
*         volume.
*
*Parameters: output
* EVECT   averaged fluxes (contained in positions EVECT(IDL(I))).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER NEL,NUN,MAT(NEL),IDL(NEL)
      REAL EVECT(NUN),VOL(NEL)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      LOGICAL CYLIND,CHEX
      INTEGER ITP(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: KN
      REAL, DIMENSION(:), ALLOCATABLE :: XX,DD,T,TS
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
      ITYPE=ITP(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IELEM=ABS(ITP(9))
      LL4=ITP(11)
      ICHX=ITP(12)
      ISPLH=ITP(13)
      LX=ITP(14)
      LY=ITP(15)
      LZ=ITP(16)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(KN(MAXKN))
      CALL LCMGET(IPTRK,'KN',KN)
*
      IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
*        LAGRANGIAN FINITE ELEMENTS.
         ALLOCATE(XX(LX*LY*LZ),DD(LX*LY*LZ))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'DD',DD)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(T(LC),TS(LC))
         CALL LCMGET(IPTRK,'T',T)
         CALL LCMGET(IPTRK,'TS',TS)
         CALL LCMSIX(IPTRK,' ',2)
         CALL FLDTN2(NEL,LL4,IELEM,CYLIND,EVECT,XX,DD,MAT,VOL,IDL,KN,
     1   LC,T,TS)
         DEALLOCATE(TS,T,DD,XX)
      ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
*        MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         CALL FLDTH2(ISPLH,NEL,NUN,EVECT,MAT,VOL,IDL,KN)
      ELSE IF((ICHX.EQ.3).AND.(ISPLH.GT.1).AND.CHEX) THEN
*        MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         CALL FLDTH1(ISPLH,NEL,LL4,EVECT,MAT,VOL,IDL,KN)
      ENDIF
*----
*  RELEASE TRIVAC SPECIFIC TRACKING INFORMATION
*----
      DEALLOCATE(KN)
      RETURN
      END
