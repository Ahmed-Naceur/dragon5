*DECK FLDTSM
      SUBROUTINE FLDTSM(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LCM driver for the multiplication of a matrix by a vector.
* Special version for the simplified PN method in TRIVAC.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
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
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER NAMP*(*)
      INTEGER LL4,NBMIX,NAN
      REAL F2(LL4),F3(LL4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER NAMT*12,TEXT12*12
      INTEGER IPAR(NSTATE)
      LOGICAL CHEX
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KN,IQFR
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,QFR,XX,YY,ZZ,GAMMA
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,V,SGD
      INTEGER, DIMENSION(:), POINTER :: IPERT
      REAL, DIMENSION(:), POINTER :: FRZ
      DOUBLE PRECISION, DIMENSION(:), POINTER :: CTRAN
      TYPE(C_PTR) IPERT_PTR,FRZ_PTR,CTRAN_PTR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SGD(NBMIX,2*NAN))
*----
*  RECOVER PN SPECIFIC PARAMETERS.
*----
      NAMT=NAMP
      READ(NAMT,'(1X,2I3)') IGR,JGR
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      NREG=IPAR(1)
      NUN=IPAR(2)
      ITYPE=IPAR(6)
      IELEM=IPAR(9)
      ICOL=IPAR(10)
      L4=IPAR(11)
      ISPLH=IPAR(13)
      LX=IPAR(14)
      LZ=IPAR(16)
      LL4F=IPAR(25)
      LL4W=IPAR(26)
      LL4X=IPAR(27)
      LL4Y=IPAR(28)
      NLF=IPAR(30)
      NVD=IPAR(34)
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IF(CHEX) THEN
         IF(NUN.GT.(LX*LZ+L4)*NLF/2) CALL XABORT('FLDTSM: INVALID NUN '
     1   //'OR L4.')
      ELSE
         IF(NUN.NE.L4*NLF/2) CALL XABORT('FLDTSM: INVALID NUN OR L4.')
      ENDIF
      IF(L4*NLF/2.NE.LL4) CALL XABORT('FLDTSM: INVALID L4 OR LL4.')
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
      IF(CHEX) THEN
         CALL LCMGET(IPTRK,'SIDE',SIDE)
      ELSE
         ALLOCATE(XX(NREG),YY(NREG))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
      ENDIF
      ALLOCATE(ZZ(NREG))
      CALL LCMGET(IPTRK,'ZZ',ZZ)
*----
*  RECOVER THE PERTURBATION FLAG.
*----
      CALL LCMGET(IPSYS,'STATE-VECTOR',IPAR)
      IPR=IPAR(9)
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
      IF(CHEX) THEN
         NBLOS=LX*LZ/3
         CALL LCMGPD(IPTRK,'CTRAN',CTRAN_PTR)
         CALL LCMGPD(IPTRK,'IPERT',IPERT_PTR)
         CALL LCMGPD(IPTRK,'FRZ',FRZ_PTR)
         CALL C_F_POINTER(CTRAN_PTR,CTRAN,(/ ((IELEM+1)*IELEM)**2 /))
         CALL C_F_POINTER(IPERT_PTR,IPERT,(/ NBLOS /))
         CALL C_F_POINTER(FRZ_PTR,FRZ,(/ NBLOS /))
         CALL PNSH3D(ITY,IPR,NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,L4,LL4F,
     1   LL4W,LL4X,LL4Y,MAT,SGD(1,1),SGD(1,NAN+1),SIDE,ZZ,FRZ,QFR,IPERT,
     2   KN,LC,R,V,CTRAN,F2,F3)
      ELSE
         CALL PNSZ3D(ITY,IPR,NREG,IELEM,ICOL,XX,YY,ZZ,MAT,VOL,NBMIX,NLF,
     1   NVD,NAN,SGD(1,1),SGD(1,NAN+1),L4,KN,QFR,LC,R,V,F2,F3)
      ENDIF
      IF(ITY.EQ.1) THEN
         DO 30 I=1,LL4
         F3(I)=-F3(I)
   30    CONTINUE
      ENDIF
      DEALLOCATE(V,R,ZZ)
      IF(.NOT.CHEX) DEALLOCATE(YY,XX)
      DEALLOCATE(IQFR,QFR,KN,VOL,MAT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SGD)
      RETURN
      END
