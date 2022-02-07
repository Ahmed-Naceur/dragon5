*DECK FLDBSS
      SUBROUTINE FLDBSS(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,F1,NADI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a one-group SPN flux iteration in Cartesian or hexagonal 2D
* geometry in BIVAC.
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
* F1      source term of the linear system.
* NADI    number of inner ADI iterations.
*
*Parameters: output
* F1      approached solution of the linear system.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER NAMP*(*)
      INTEGER LL4,NBMIX,NAN,NADI
      REAL F1(LL4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER NAMT*12,TEXT12*12
      INTEGER IPAR(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KEY,MU,KN,IQFR,IPERT
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,QFR,SOUR,SYS,XX,YY,GAMMA
      REAL, DIMENSION(:,:), ALLOCATABLE :: SGD,R,V
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SGD(NBMIX,NAN))
*----
*  RECOVER PN SPECIFIC PARAMETERS.
*----
      NAMT=NAMP
      READ(NAMT,'(1X,2I3)') IGR,JGR
      IF(IGR.NE.JGR) CALL XABORT('FLDBSS: INVALIB GROUP INDICES.')
      IF(NAMT(1:1).NE.'A') CALL XABORT('FLDBSS: ''A'' MATRIX EXPECTED.')
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
         IF(NUN.GT.(LX+L4)*NLF/2) CALL XABORT('FLDBSS: INVALID NUN OR '
     1   //'L4.')
      ELSE
         IF(NUN.NE.L4*NLF/2) CALL XABORT('FLDBSS: INVALID NUN OR L4.')
      ENDIF
      IF(L4*NLF/2.NE.LL4) CALL XABORT('FLDBSS: INVALID L4 OR LL4.')
      ALLOCATE(MAT(NREG),KEY(NREG),VOL(NREG),MU(L4))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'KEYFLX',KEY)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'MU',MU)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(KN(MAXKN),QFR(MAXQF),IQFR(MAXQF))
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'IQFR',IQFR)
*----
*  PROCESS PHYSICAL ALBEDO FUNCTIONS
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
      DO 10 IL=1,NAN
      WRITE(TEXT12,'(4HSCAI,I2.2,A6)') IL-1,NAMT(2:7)
      CALL LCMGET(IPSYS,TEXT12,SGD(1,IL))
   10 CONTINUE
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
*  SOLVE THE LINEAR SYSTEM.
*----
      IIMAX=MU(L4)*NLF/2
      ALLOCATE(SYS(IIMAX),SOUR(NUN))
      CALL LCMGET(IPSYS,'I'//NAMT,SYS)
      DO 30 IUN=1,NUN
      SOUR(IUN)=F1(IUN)
      F1(IUN)=0.0
   30 CONTINUE
      IF((ITYPE.EQ.2).OR.((ITYPE.EQ.5).AND.(ISPN.EQ.1))) THEN
         ALLOCATE(XX(NREG),YY(NREG))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
         CALL PNFL2E(NREG,IELEM,ICOL,XX,YY,MAT,VOL,NBMIX,NLF,NVD,NAN,
     1   SGD,L4,KN,QFR,MU,IIMAX,LC,R,V,SYS,SOUR,F1,NADI)
         DEALLOCATE(YY,XX)
      ELSE IF(ITYPE.EQ.8) THEN
         CALL LCMGET(IPTRK,'SIDE',SIDE)
         NBLOS=LX/3
         ALLOCATE(IPERT(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL PNFH2E(IELEM,ICOL,NBLOS,SIDE,NLF,NVD,L4,IPERT,KN,QFR,MU,
     1   IIMAX,LC,V,SYS,SOUR,F1,NADI)
         DEALLOCATE(IPERT)
      ENDIF
      DEALLOCATE(SOUR,SYS,V,R,IQFR,QFR,KN,MU,VOL,KEY,MAT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SGD)
      RETURN
      END
