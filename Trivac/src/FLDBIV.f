*DECK FLDBIV
      SUBROUTINE FLDBIV(IPTRK,NEL,NUN,EVECT,MAT,VOL,IDL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the averaged flux in BIVAC.
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
* IPTRK   L_TRACK pointer to the BIVAC tracking information.
* NEL     total number of finite elements.
* NUN     total number of unknown per energy group.
* EVECT   variational coefficients of the flux. The information is
*         contained in position EVECT(1) to EVECT(LL4) where LL4 is
*         the order of the system matrices.
* MAT     mixture index assigned to each element.
* VOL     volume of each element.
* IDL     position of the average flux component associated with each
*         volume.
*
*Parameters: output
* EVECT   averaged fluxes. The information is contained in positions
*         EVECT(IDL(I)).
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
      LOGICAL CYLIND
      INTEGER ITP(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: KN
      REAL, DIMENSION(:), ALLOCATABLE :: XX,DD,T,TS,QFR
      REAL, DIMENSION(:,:), ALLOCATABLE :: RH,RT
*----
*  RECOVER BIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
      ITYPE=ITP(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IELEM=ITP(8)
      ICOL=ITP(9)
      ISPLH=ITP(10)
      LL4=ITP(11)
      LX=ITP(12)
      LY=ITP(13)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(KN(MAXKN))
      CALL LCMGET(IPTRK,'KN',KN)
*
      IF((IELEM.LT.0).AND.(ITYPE.NE.8)) THEN
*        LAGRANGIAN FINITE ELEMENTS.
         ALLOCATE(XX(LX*LY),DD(LX*LY))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'DD',DD)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(T(LC),TS(LC))
         CALL LCMGET(IPTRK,'T',T)
         CALL LCMGET(IPTRK,'TS',TS)
         CALL LCMSIX(IPTRK,' ',2)
         CALL FLDBN2(NEL,LL4,-IELEM,CYLIND,EVECT,XX,DD,MAT,VOL,IDL,KN,
     1   LC,T,TS)
         DEALLOCATE(TS,T,DD,XX)
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
*        MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         ALLOCATE(RH(6,6),RT(3,3))
         CALL LCMGET(IPTRK,'RH',RH)
         CALL LCMGET(IPTRK,'RT',RT)
         CALL LCMSIX(IPTRK,' ',2)
         IF(ISPLH.EQ.1) THEN
            NELEM=MAXKN/7
         ELSE
            NELEM=MAXKN/4
         ENDIF
         ALLOCATE(QFR(MAXKN))
         CALL LCMGET(IPTRK,'QFR',QFR)
         CALL FLDBH2(ISPLH,NEL,NUN,NELEM,EVECT,VOL,IDL,KN,QFR,RH,RT)
         DEALLOCATE(QFR,RT,RH)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(ICOL.EQ.4).AND.
     1        (ISPLH.GT.1)) THEN
*        MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         ALLOCATE(QFR(MAXKN))
         CALL LCMGET(IPTRK,'QFR',QFR)
         CALL FLDBH1(NEL,NUN,LL4,EVECT,VOL,IDL,KN,QFR)
         DEALLOCATE(QFR)
      ENDIF
*----
*  RELEASE BIVAC SPECIFIC TRACKING INFORMATION
*----
      DEALLOCATE(KN)
      RETURN
      END
