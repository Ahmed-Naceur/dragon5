*DECK KINTLM
      SUBROUTINE KINTLM(IPTRK,NBM,LDIM,SGD,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the multiplication of a matrix by a vector. Special
* version for Trivac.
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
* IPTRK   L_TRACK pointer to the tracking information.
* NBM     number of material mixtures.     
* LDIM    dimension of vectors F2 and F3.
* SGD     mixture-ordered cross sections.
* F2      vector to multiply.
*
*Parameters: output
* F3      result of the multiplication.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER NBM,LDIM
      REAL SGD(NBM),F2(LDIM),F3(LDIM)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      LOGICAL CYLIND,CHEX
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KN,IPERT,IPW,XORZ,DD
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,T,TS,FRZ
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,RH,RT
*----
*  RECOVER TRACKING INFORMATION.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NBMIX=ISTATE(4)
      ITYPE=ISTATE(6)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(MAT(NREG),VOL(NREG),KN(MAXKN))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KN',KN)
*----
*  ALGORITHM-DEPENDENT MULTIPLICATION
*----
      CALL XDRSET(F3,LDIM,0.0)
      ITYPE=ISTATE(6)
      IDIM=1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IHEX=ISTATE(7)
      IELEM=ISTATE(9)
      ICOL=ISTATE(10)
      LL4=ISTATE(11)
      ICHX=ISTATE(12)
      IF(ICHX.EQ.2) LL4=ISTATE(25)
      ISPLH=ISTATE(13)
      LX=ISTATE(14)
      LY=ISTATE(15)
      LZ=ISTATE(16)
      NVD=ISTATE(34)
      IF(LL4.GT.LDIM) CALL XABORT('KINTLM: LDIM OVERFLOW.')
      ALLOCATE(XORZ(LX*LY*LZ),DD(LX*LY*LZ))
      IF(CHEX) THEN
        CALL LCMGET(IPTRK,'ZZ',XORZ)
        CALL LCMGET(IPTRK,'SIDE',SIDE)
      ELSE
        CALL LCMGET(IPTRK,'XX',XORZ)
        CALL LCMGET(IPTRK,'DD',DD)
      ENDIF
      IF((.NOT.CHEX).AND.(ICHX.EQ.1)) THEN
*       --- MIXED-PRIMAL FINITE ELEMENTS (CARTESIAN)
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
        ALLOCATE(T(LC),TS(LC))
        CALL LCMGET(IPTRK,'T',T)
        CALL LCMGET(IPTRK,'TS',TS)
        CALL LCMSIX(IPTRK,' ',2)
        CALL KINT01(MAXKN,SGD,CYLIND,NREG,LL4,NBMIX,XORZ,DD,MAT,KN,VOL,
     1  LC,T,TS,F2,F3)
        DEALLOCATE(TS,T)
      ELSEIF((.NOT.CHEX).AND.(ICHX.GE.2)) THEN
*       --- DUAL FINITE ELEMENTS (CARTESIAN)
        CALL KINT02(MAXKN,SGD,IELEM,ICHX,IDIM,NREG,LL4,NBMIX,MAT,KN,
     1  VOL,F2,F3)
      ELSEIF(CHEX.AND.(ICHX.EQ.1)) THEN
*       --- MESH CORNER FINITE DIFFERENCES (HEXAGONAL)
        ALLOCATE(R(2,2),RH(6,6),RT(3,3))
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMGET(IPTRK,'R',R)
        CALL LCMGET(IPTRK,'RH',RH)
        CALL LCMGET(IPTRK,'RT',RT)
        CALL LCMSIX(IPTRK,' ',2)
        CALL KINT03(MAXKN,ISPLH,NBMIX,NREG,LL4,SGD,SIDE,XORZ,VOL,MAT,
     1  KN,R,RH,RT,F2,F3)
        DEALLOCATE(RT,RH,R)
      ELSEIF(CHEX.AND.(ICHX.EQ.2)) THEN
*       --- DUAL (THOMAS-RAVIART-SCHNEIDER) FINITE ELEMENT METHOD.
        NBLOS=LX*LZ/3
        ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
        CALL LCMGET(IPTRK,'IPERT',IPERT)
        CALL LCMGET(IPTRK,'FRZ',FRZ)
        CALL KINT04(IELEM,NBMIX,LL4,NBLOS,MAT,SIDE,XORZ,FRZ,SGD,KN,
     1  IPERT,F2,F3)
        DEALLOCATE(FRZ,IPERT)
      ELSE IF(CHEX.AND.(ICHX.EQ.3).AND.(ISPLH.EQ.1)) THEN
*       --- MESH CENTERED FINITE DIFFERENCES (HEXAGONAL)
        CALL KINT05(NBMIX,NREG,LL4,SGD,VOL,MAT,F2,F3)
      ELSE IF(CHEX.AND.(ICHX.EQ.3).AND.(ISPLH.GT.1)) THEN
*       --- MESH CENTERED FINITE DIFFERENCES (HEXAGONAL)
        ALLOCATE(IPW(LL4))
        CALL LCMGET(IPTRK,'IPW',IPW)
        CALL KINT06(ISPLH,NBMIX,NREG,LL4,VOL,MAT,SGD,KN,IPW,F2,F3)
        DEALLOCATE(IPW)
      ELSE
        CALL XABORT('KINTLM: TRACKING NOT AVAILABLE.')
      ENDIF
      DEALLOCATE(DD,XORZ,KN,VOL,MAT)
      RETURN
      END
