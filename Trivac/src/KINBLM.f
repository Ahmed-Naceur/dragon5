*DECK KINBLM
      SUBROUTINE KINBLM(IPTRK,NBM,LDIM,SGD,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the multiplication of a matrix by a vector. Special
* version for Bivac.
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
      LOGICAL CYLIND
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KN,IPERT
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,QFR,XX,DD
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,RS,RH,RT
*----
*  RECOVER TRACKING INFORMATION.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NBMIX=ISTATE(4)
      ITYPE=ISTATE(6)
      ALLOCATE(MAT(NREG),VOL(NREG))
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(KN(MAXKN),QFR(MAXQF))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
*----
*  ALGORITHM-DEPENDENT MULTIPLICATION
*----
      CALL XDRSET(F3,LDIM,0.0)
      ITYPE=ISTATE(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IHEX=ISTATE(7)
      IELEM=ISTATE(8)
      ICOL=ISTATE(9)
      ISPLH=ISTATE(10)
      LL4=ISTATE(11)
      LX=ISTATE(12)
      LY=ISTATE(13)
      NVD=ISTATE(17)
      IF(LL4.GT.LDIM) CALL XABORT('KINBLM: LDIM OVERFLOW.')
      ALLOCATE(XX(LX*LY),DD(LX*LY))
      IF(ITYPE.EQ.8) THEN
        CALL LCMGET(IPTRK,'SIDE',SIDE)
      ELSE
        CALL LCMGET(IPTRK,'XX',XX)
        CALL LCMGET(IPTRK,'DD',DD)
      ENDIF
      IF((IHEX.EQ.0).AND.(IELEM.LT.0)) THEN
*       --- PRIMAL FINITE ELEMENTS (CARTESIAN)
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
        ALLOCATE(R(LC,LC),RS(LC,LC))
        CALL LCMGET(IPTRK,'R',R)
        CALL LCMGET(IPTRK,'RS',RS)
        CALL LCMSIX(IPTRK,' ',2)
        CALL KINB01(MAXKN,SGD,CYLIND,NREG,LL4,NBMIX,XX,DD,MAT,KN,VOL,
     1  LC,R,RS,F2,F3)
        DEALLOCATE(RS,R)
      ELSE IF((IHEX.EQ.0).AND.(IELEM.GT.0)) THEN
*       --- MIXED-DUAL FINITE ELEMENTS (CARTESIAN)
        CALL KINB02(SGD,IELEM,NREG,LL4,NBMIX,MAT,KN,VOL,F2,F3)
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
*       --- MESH CORNER FINITE DIFFERENCES (HEXAGONAL)
        ALLOCATE(RH(6,6),RT(3,3))
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMGET(IPTRK,'RH',RH)
        CALL LCMGET(IPTRK,'RT',RT)
        CALL LCMSIX(IPTRK,' ',2)
        IF(ISPLH.EQ.1) THEN
          NELEM=MAXKN/7
        ELSE
          NELEM=MAXKN/4
        ENDIF
        CALL KINB03(MAXKN,MAXQF,SGD,NREG,LL4,ISPLH,NELEM,NBMIX,MAT,KN,
     1  QFR,VOL,RH,RT,F2,F3)
        DEALLOCATE(RT,RH)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(ICOL.EQ.4)) THEN
*       --- MESH CENTERED FINITE DIFFERENCES FOR HEXAGONS
        CALL KINB04(MAXKN,MAXQF,SGD,NREG,LL4,ISPLH,NBMIX,MAT,KN,QFR,
     1  VOL,F2,F3)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8)) THEN
*       --- THOMAS-RAVIART-SCHNEIDER METHOD (HEXAGONAL)
        NBLOS=LX/3
        ALLOCATE(IPERT(NBLOS))
        CALL LCMGET(IPTRK,'IPERT',IPERT)
        CALL KINB05(SGD,IELEM,NBLOS,LL4,NBMIX,SIDE,MAT,IPERT,KN,F2,F3)
        DEALLOCATE(IPERT)
      ELSE
        CALL XABORT('KINBLM: TRACKING NOT AVAILABLE.')
      ENDIF
      DEALLOCATE(DD,XX,QFR,KN,VOL,MAT)
      RETURN
      END
