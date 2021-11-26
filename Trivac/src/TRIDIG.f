*DECK TRIDIG
      SUBROUTINE TRIDIG(HNAME,IPTRK,IPSYS,IMPX,MAXMIX,NEL,IPR,MAT,VOL,
     1 SGD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the assembly of a cross section diagonal matrix.
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
* HNAME   name of the diagonal matrix.
* IPTRK   L_TRACK pointer to the TRIVAC tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IMPX    print parameter. Equal to zero for no print.
* MAXMIX   dimension of matrix SGD.
* NEL     total number of finite elements.
* IPR     type of assembly:
*         =3: the new contribution is added to existing matrix.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* SGD     cross section per material mixture.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER HNAME*(*)
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER IMPX,MAXMIX,NEL,IPR,MAT(NEL)
      REAL VOL(NEL),SGD(MAXMIX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12
      LOGICAL CYLIND,CHEX
      INTEGER ISTATE(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: KN,IPERT,IPW
      REAL, DIMENSION(:), ALLOCATABLE :: T,TS,FRZ
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,RH,RT
      REAL, DIMENSION(:), ALLOCATABLE :: XORZ,DD
      REAL, DIMENSION(:), POINTER :: VEC
      TYPE(C_PTR) VEC_PTR
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IELEM=ABS(ISTATE(9))
      LL4=ISTATE(11)
      ICHX=ISTATE(12)
      ISPLH=ISTATE(13)
      LX=ISTATE(14)
      LY=ISTATE(15)
      LZ=ISTATE(16)
      LL4F=ISTATE(25)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(KN(MAXKN),XORZ(LX*LY*LZ))
      CALL LCMGET(IPTRK,'KN',KN)
      IF(CHEX) THEN
         CALL LCMGET(IPTRK,'ZZ',XORZ)
         CALL LCMGET(IPTRK,'SIDE',SIDE)
      ELSE
         CALL LCMGET(IPTRK,'XX',XORZ)
         ALLOCATE(DD(LX*LY*LZ))
         CALL LCMGET(IPTRK,'DD',DD)
      ENDIF
      TEXT12=HNAME
      IF(IMPX.GT.0) WRITE(6,'(/37H TRIDIG: ASSEMBLY OF DIAGONAL MATRIX ,
     1 1H'',A12,2H''.)') TEXT12
*----
*  INITIALIZATION OF A DIAGONAL SYSTEM MATRIX
*----
      IF(ICHX.EQ.2) THEN
         IF(IPR.EQ.3) THEN
            CALL LCMGPD(IPSYS,TEXT12,VEC_PTR)
            CALL C_F_POINTER(VEC_PTR,VEC,(/ LL4F /))
         ELSE
            VEC_PTR=LCMARA(LL4F)
            CALL C_F_POINTER(VEC_PTR,VEC,(/ LL4F /))
            CALL XDRSET(VEC,LL4F,0.0)
         ENDIF
      ELSE
         IF(IPR.EQ.3) THEN
            CALL LCMGPD(IPSYS,TEXT12,VEC_PTR)
            CALL C_F_POINTER(VEC_PTR,VEC,(/ LL4 /))
         ELSE
            VEC_PTR=LCMARA(LL4)
            CALL C_F_POINTER(VEC_PTR,VEC,(/ LL4 /))
            CALL XDRSET(VEC,LL4,0.0)
         ENDIF
      ENDIF
*----
*  COMPUTE THE DIAGONAL SYSTEM MATRIX
*----
      IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
*        VARIATIONAL COLLOCATION METHOD.
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(T(LC),TS(LC))
         CALL LCMGET(IPTRK,'T',T)
         CALL LCMGET(IPTRK,'TS',TS)
         CALL LCMSIX(IPTRK,' ',2)
         CALL TRIASP(IELEM,MAXMIX,NEL,LL4,CYLIND,SGD,XORZ,DD,VOL,MAT,
     1   KN,LC,T,TS,VEC)
         DEALLOCATE(T,TS)
      ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
*        MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         ALLOCATE(R(2,2),RH(6,6),RT(3,3))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'RH',RH)
         CALL LCMGET(IPTRK,'RT',RT)
         CALL LCMSIX(IPTRK,' ',2)
         CALL TRIAHP(MAXKN,ISPLH,MAXMIX,NEL,LL4,SGD,SIDE,XORZ,VOL,MAT,
     1   KN,R,RH,RT,VEC)
         DEALLOCATE(RT,RH,R)
      ELSE IF((ICHX.EQ.2).AND.CHEX) THEN
*        DUAL (THOMAS-RAVIART-SCHNEIDER) FINITE ELEMENT METHOD.
         NBLOS=LX*LZ/3
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMGET(IPTRK,'FRZ',FRZ)
         CALL TRIASH(IELEM,MAXMIX,LL4,NBLOS,MAT,SIDE,XORZ,FRZ,SGD,KN,
     1   IPERT,VEC)
         DEALLOCATE(FRZ,IPERT)
      ELSE IF(.NOT.CHEX) THEN
*        DUAL FINITE ELEMENT METHOD.
         IDIM=1
         IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
         IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
         CALL TRIASD(MAXKN,IELEM,ICHX,IDIM,MAXMIX,NEL,LL4,SGD,VOL,MAT,
     1   KN,VEC)
      ELSE IF(CHEX.AND.(ISPLH.EQ.1)) THEN
*        MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         CALL TRIAHD(MAXMIX,NEL,LL4,SGD,VOL,MAT,VEC)
      ELSE IF(CHEX.AND.(ISPLH.GT.1)) THEN
*        MESH CENTERED FINITE DIFFERENCES IN TRIANGULAR GEOMETRY.
         ALLOCATE(IPW(LL4))
         CALL LCMGET(IPTRK,'IPW',IPW)
         CALL TRIMTD(ISPLH,MAXMIX,NEL,LL4,VOL,MAT,SGD,KN,IPW,VEC)
         DEALLOCATE(IPW)
      ENDIF
*----
*  STORAGE OF THE DIAGONAL SYSTEM MATRIX
*----
      IF(ICHX.EQ.2) THEN
         CALL LCMPPD(IPSYS,TEXT12,LL4F,2,VEC_PTR)
      ELSE
         CALL LCMPPD(IPSYS,TEXT12,LL4,2,VEC_PTR)
      ENDIF
*----
*  RELEASE TRIVAC SPECIFIC TRACKING INFORMATION
*----
      IF(.NOT.CHEX) DEALLOCATE(DD)
      DEALLOCATE(XORZ,KN)
      RETURN
      END
