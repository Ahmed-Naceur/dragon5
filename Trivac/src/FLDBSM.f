*DECK FLDBSM
      SUBROUTINE FLDBSM(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LCM driver for the multiplication of a matrix by a vector.
* Special version for the simplified PN method in BIVAC.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NAMP    name of the coefficient matrix.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* LL4     order of the matrix.
* NBMIX   total number of material mixtures in the macrolib.
* NAN     number of Legendre orders in the cross sections.
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
      CHARACTER NAMP*(*)
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER LL4,NBMIX,NAN
      REAL F2(LL4),F3(LL4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER NAMT*12,TEXT12*12
      INTEGER IPAR(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KN,IQFR,IPERT
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,QFR,XX,YY,GAMMA
      REAL, DIMENSION(:,:), ALLOCATABLE :: SGD,V,R,H
      INTEGER, DIMENSION(:), POINTER :: MU
      REAL, DIMENSION(:), POINTER :: ASS
      TYPE(C_PTR) MU_PTR,ASS_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SGD(NBMIX,2*NAN))
*----
*  RECOVER ENERGY GROUP INDICES.
*----
      NAMT=NAMP
      READ(NAMT,'(1X,2I3)') IGR,JGR
*----
*  RECOVER PN SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      NREG=IPAR(1)
      NUN=IPAR(2)
      ITYPE=IPAR(6)
      IELEM=IPAR(8)
      ICOL=IPAR(9)
      ISPLH=IPAR(10)
      L4=IPAR(11)
      LX=IPAR(12)
      NLF=IPAR(14)
      ISPN=IPAR(15)
      NVD=IPAR(17)
      IF(ITYPE.EQ.8) THEN
         IF(NUN.GT.(LX+L4)*NLF/2) CALL XABORT('FLDBSM: INVALID NUN OR '
     1   //'L4.')
      ELSE
        IF(NUN.NE.L4*NLF/2) CALL XABORT('FLDBSM: INVALID NUN OR L4.')
      ENDIF
      IF(L4*NLF/2.NE.LL4) CALL XABORT('FLDBSM: INVALID L4 OR LL4.')
*----
*  PROCESS A FISSION MATRIX.
*----
      IF(NAMT(1:1).EQ.'B') THEN
         CALL LCMLEN(IPTRK,'MU',LL4TS,ITYLCM)
         IF(L4.NE.LL4TS) CALL XABORT('FLDBSM: INVALID L4.')
         CALL LCMGPD(IPTRK,'MU',MU_PTR)
         CALL LCMGPD(IPSYS,NAMT,ASS_PTR)
         CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
         CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(L4) /))
         CALL ALLDLM(L4,ASS,F2,F3,MU,1)
         RETURN
      ELSE IF(NAMT(1:1).NE.'A') THEN
         CALL XABORT('FLDBSM: ''A'' OR ''B'' MATRIX EXPECTED.')
      ENDIF
*----
*  RECOVER TRACKING INFORMATION.
*----
      ALLOCATE(MAT(NREG),VOL(NREG))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(KN(MAXKN),QFR(MAXQF),IQFR(MAXQF))
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'IQFR',IQFR)
*----
*  PROCESS PHYSICAL ALBEDOS
*----
      TEXT12='ALBEDO-FU'//NAMT(2:4)
      CALL LCMLEN(IPSYS,TEXT12,NALBP,ITYLCM)
      IF(NALBP.GT.0) THEN
         ALLOCATE(GAMMA(NALBP))
         CALL LCMGET(IPSYS,TEXT12,GAMMA)
         DO IQW=1,MAXQF
            IALB=IQFR(IQW)
            IF(IALB.NE.0) QFR(IQW)=QFR(IQW)*GAMMA(IALB)
         ENDDO
         DEALLOCATE(GAMMA)
      ENDIF
*----
*  RECOVER THE CROSS SECTIONS.
*----
      DO 20 IL=1,NAN
      WRITE(TEXT12,'(4HSCAR,I2.2,A6)') IL-1,NAMT(2:7)
      CALL LCMLEN(IPSYS,TEXT12,LENGT,ITYLCM)
      IF(LENGT.EQ.0) THEN
         CALL XDRSET(SGD(1,IL),NBMIX,0.0)
         CALL XDRSET(SGD(1,NAN+IL),NBMIX,0.0)
      ELSE
         CALL LCMGET(IPSYS,TEXT12,SGD(1,IL))
         WRITE(TEXT12,'(4HSCAI,I2.2,A6)') IL-1,NAMT(2:7)
         CALL LCMGET(IPSYS,TEXT12,SGD(1,NAN+IL))
      ENDIF
   20 CONTINUE
*----
*  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
*----
      CALL LCMSIX(IPTRK,'BIVCOL',1)
      CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
      ALLOCATE(R(LC,LC),V(LC,LC-1))
      CALL LCMGET(IPTRK,'R',R)
      CALL LCMGET(IPTRK,'V',V)
      CALL LCMSIX(IPTRK,' ',2)
*----
*  COMPUTE THE SOURCE
*----
      ITY=0
      IF(IGR.NE.JGR) ITY=1
      IF((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1))) THEN
         ALLOCATE(XX(NREG),YY(NREG))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
         CALL PNSZ2D(ITY,NREG,IELEM,ICOL,XX,YY,MAT,VOL,NBMIX,NLF,NVD,
     1   NAN,SGD(1,1),SGD(1,NAN+1),L4,KN,QFR,LC,R,V,F2,F3)
         DEALLOCATE(YY,XX)
      ELSE IF(ITYPE.EQ.8) THEN
         NBLOS=LX/3
         CALL LCMGET(IPTRK,'SIDE',SIDE)
         ALLOCATE(IPERT(NBLOS),H(LC,LC-1))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMGET(IPTRK,'H',H)
         CALL LCMSIX(IPTRK,' ',2)
         CALL PNSH2D(ITY,IELEM,ICOL,NBLOS,SIDE,MAT,NBMIX,NLF,NVD,
     1   NAN,SGD(1,1),L4,IPERT,KN,QFR,LC,R,V,H,F2,F3)
         DEALLOCATE(H,IPERT)
      ENDIF
      IF(ITY.EQ.1) THEN
         DO 30 I=1,LL4
         F3(I)=-F3(I)
   30    CONTINUE
      ENDIF
      DEALLOCATE(V,R,IQFR,QFR,KN,VOL,MAT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SGD)
      RETURN
      END
