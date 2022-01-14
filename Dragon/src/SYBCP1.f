*DECK SYBCP1
      SUBROUTINE SYBCP1 (IPTRK,ITG,IMPX,NREG,SIGT,SIGW,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the scattering-reduced collision probabilities for
* Sybil.
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
* IPTRK   pointer to the Sybil tracking (L_TRACK signature).
* ITG     type of Sybil one-speed solution operator.
* IMPX    print flag (equal to zero for no print).
* NREG    total number of regions.
* SIGT    total macroscopic cross sections ordered by volume.
* SIGW    P0 within-group scattering macroscopic cross sections
*         ordered by volume.
*
*Parameters: output
* PIJ     scattering-reduced collision probabilities matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER ITG,IMPX,NREG
      REAL SIGT(NREG),SIGW(NREG),PIJ(NREG,NREG)
*----
*  LOCAL VARIABLES
*----
      INTEGER IPAR(16)
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
      REAL, TARGET, SAVE, DIMENSION(1) :: DUMMY
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NCODE,NMC3,LSEC4,NMC4,NMCR4,
     1 MAIL,IFR,INUM,MIX,IGEN
      REAL, ALLOCATABLE, DIMENSION(:) :: XX2,ZCODE,ZTR,RAYR3,PROCE,XX4,
     1 YY4,RAYR4,ALB,DVX
      INTEGER, POINTER, DIMENSION(:) :: IZMAI
      REAL, POINTER, DIMENSION(:) :: RZMAI
      TYPE(C_PTR) :: IZMAI_PTR,RZMAI_PTR
*
      IF(ITG.EQ.1) THEN
         PIJ(1,1)=1.0/(SIGT(1)-SIGW(1))
      ELSE IF(ITG.EQ.2) THEN
         CALL LCMSIX(IPTRK,'PURE-GEOM',1)
         CALL LCMGET(IPTRK,'PARAM',IPAR)
         ITYPE=IPAR(1)
         IHEX=IPAR(2)
         IQUA2=IPAR(3)
         IF(ITYPE.GE.8) CALL LCMGET(IPTRK,'SIDE',SIDE)
         ALLOCATE(XX2(NREG+1),NCODE(6),ZCODE(6))
         CALL LCMGET(IPTRK,'XXX',XX2)
         CALL LCMGET(IPTRK,'NCODE',NCODE)
         CALL LCMGET(IPTRK,'ZCODE',ZCODE)
         CALL LCMSIX(IPTRK,' ',2)
*
         IF(ITYPE.EQ.2) THEN
            CALL SYBALP(NREG,NREG,XX2,SIGT,NCODE,ZCODE,PIJ)
         ELSE IF(ITYPE.EQ.3) THEN
            ALLOCATE(ZTR(1+IQUA2*((NREG*(5+NREG))/2)))
            CALL SYBT1D(NREG,XX2,.FALSE.,IQUA2,ZTR)
            CALL SYBALC(NREG,NREG,XX2,SIGT,IQUA2,ZCODE(2),ZTR,PIJ)
            DEALLOCATE(ZTR)
         ELSE IF(ITYPE.EQ.4) THEN
            ALLOCATE(ZTR(1+IQUA2*((NREG*(5+NREG))/2)))
            CALL SYBT1D(NREG,XX2,.TRUE.,IQUA2,ZTR)
            CALL SYBALS(NREG,NREG,XX2,SIGT,IQUA2,ZCODE(2),ZTR,PIJ)
            DEALLOCATE(ZTR)
         ENDIF
         DEALLOCATE(ZCODE,NCODE,XX2)
         CALL SYBWIJ(NREG,NREG,SIGW,PIJ)
      ELSE IF(ITG.EQ.3) THEN
         CALL LCMSIX(IPTRK,'DOITYOURSELF',1)
         CALL LCMGET(IPTRK,'PARAM',IPAR)
         NSUPCE=IPAR(1)
         IQUA3=IPAR(2)
         ISTAT=IPAR(3)
         ALLOCATE(NMC3(NSUPCE+1),RAYR3(NSUPCE+NREG),PROCE(NSUPCE**2))
         CALL LCMGET(IPTRK,'NMC',NMC3)
         CALL LCMGET(IPTRK,'RAYRE',RAYR3)
         CALL LCMGET(IPTRK,'PROCEL',PROCE)
         CALL LCMSIX(IPTRK,' ',2)
         NPIJ=0
         DO 10 IKG=1,NSUPCE
         J2=NMC3(IKG+1)-NMC3(IKG)
         NPIJ=NPIJ+J2*J2
   10    CONTINUE
*
         CALL SYBRXE(NREG,NPIJ,NSUPCE,RAYR3,SIGT,SIGW,PIJ,IQUA3,ISTAT,
     1   NMC3,PROCE,IMPX)
         DEALLOCATE(PROCE,RAYR3,NMC3)
      ELSE IF(ITG.EQ.4) THEN
         CALL LCMSIX(IPTRK,'EURYDICE',1)
         CALL LCMGET(IPTRK,'PARAM',IPAR)
         IHEX=IPAR(1)
         MULTC=IPAR(2)
         IWIGN=IPAR(3)
         NMCEL=IPAR(4)
         NMERGE=IPAR(5)
         NGEN=IPAR(6)
         IJAT=IPAR(7)
         LMAILI=IPAR(15)
         LMAILR=IPAR(16)
         ALLOCATE(LSEC4(NGEN),NMC4(NGEN+1),NMCR4(NGEN+1),MAIL(2*NGEN))
         ALLOCATE(XX4(NGEN),YY4(NGEN))
         CALL LCMGET(IPTRK,'XX',XX4)
         CALL LCMGET(IPTRK,'YY',YY4)
         CALL LCMGET(IPTRK,'LSECT',LSEC4)
         CALL LCMGET(IPTRK,'NMC',NMC4)
         CALL LCMGET(IPTRK,'NMCR',NMCR4)
         CALL LCMGET(IPTRK,'MAIL',MAIL)
         ALLOCATE(RAYR4(NMCR4(NGEN+1)))
         CALL LCMGET(IPTRK,'RAYRE',RAYR4)
         IF(LMAILI.GT.0) THEN
            CALL LCMGPD(IPTRK,'ZMAILI',IZMAI_PTR)
            CALL C_F_POINTER(IZMAI_PTR,IZMAI,(/ LMAILI /))
         ELSE
*           THIS INFO IS NOT REQUIRED IN THE CALLED ROUTINE.
            IZMAI=>IDUMMY
         ENDIF
         IF(LMAILR.GT.0) THEN
            CALL LCMGPD(IPTRK,'ZMAILR',RZMAI_PTR)
            CALL C_F_POINTER(RZMAI_PTR,RZMAI,(/ LMAILR /))
         ELSE
*           THIS INFO IS NOT REQUIRED IN THE CALLED ROUTINE.
            RZMAI=>DUMMY
         ENDIF
         NCOUR=4
         IF(IHEX.NE.0) NCOUR=6
         IF(MULTC.EQ.4) NCOUR=3*NCOUR
         ALLOCATE(IFR(NCOUR*NMCEL),INUM(NMCEL),MIX(NCOUR*NMERGE),
     1   IGEN(NMERGE))
         ALLOCATE(ALB(NCOUR*NMCEL),DVX(NCOUR*NMERGE))
         CALL LCMGET(IPTRK,'IFR',IFR)
         CALL LCMGET(IPTRK,'ALB',ALB)
         CALL LCMGET(IPTRK,'INUM',INUM)
         CALL LCMGET(IPTRK,'MIX',MIX)
         CALL LCMGET(IPTRK,'DVX',DVX)
         CALL LCMGET(IPTRK,'IGEN',IGEN)
         CALL LCMSIX(IPTRK,' ',2)
*
         NPIJ=0
         DO 20 IKG=1,NGEN
         J2=NMC4(IKG+1)-NMC4(IKG)
         NPIJ=NPIJ+J2*J2
   20    CONTINUE
         NPIS=NMC4(NGEN+1)
         IF(MULTC.EQ.1) THEN
            CALL SYBRX2(NREG,NPIJ,NPIS,SIGT,SIGW,PIJ,IMPX,NCOUR,
     1      IWIGN,NMCEL,NMERGE,NGEN,IPAR(8),XX4,YY4,NMC4,RAYR4,MAIL,
     2      RZMAI,IFR,ALB,INUM,IGEN)
         ELSE
            NRAYRE=NMCR4(NGEN+1)
            CALL SYBRX3(MULTC,NREG,NPIJ,NPIS,NRAYRE,SIGT,SIGW,PIJ,IMPX,
     1      NCOUR,IWIGN,NMCEL,NMERGE,NGEN,IJAT,IPAR(8),XX4,YY4,LSEC4,
     2      NMC4,NMCR4,RAYR4,MAIL,IZMAI,RZMAI,IFR,ALB,INUM,MIX,DVX,IGEN)
         ENDIF
         DEALLOCATE(DVX,ALB)
         DEALLOCATE(IGEN,MIX,INUM,IFR)
         DEALLOCATE(RAYR4,YY4,XX4)
         DEALLOCATE(MAIL,NMCR4,NMC4,LSEC4)
      ELSE
         CALL XABORT('SYBCP1: UNKNOWN CP MODULE.')
      ENDIF
*
      IF(IMPX.GE.7) THEN
         WRITE (6,1130) (J,J=1,NREG)
         DO 90 I=1,NREG
         WRITE (6,1140) I,(PIJ(I,J),J=1,NREG)
   90    CONTINUE
         WRITE (6,'(//)')
      ENDIF
      RETURN
*
 1130 FORMAT (//49H SYBCP1: SCATTERING-REDUCED COLLISION PROBABILITY,
     1 9H MATRIX ://(11X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,
     2 I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,
     3 2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4))
 1140 FORMAT (3H I=,I4,2H: ,1P,11E11.3/(9X,11E11.3))
      END
