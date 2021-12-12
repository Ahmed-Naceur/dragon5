*DECK SYBILA
      SUBROUTINE SYBILA (IPSYS,IPTRK,IMPX,NREG,NBMIX,MAT,SIGT0,SIGW0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of cellwise scattering-reduced collision, escape and
* transmission probabilities for the current iteration method in
* Eurydice (Sybil).
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
* IPSYS   pointer to the system matrices.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NREG    total number of merged regions for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* MAT     index-number of the mixture type assigned to each volume.
* SIGT0   total macroscopic cross sections ordered by mixture.
* SIGW0   within-group scattering macroscopic cross section ordered
*         by mixture.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER IMPX,NREG,NBMIX,MAT(NREG)
      REAL SIGT0(0:NBMIX),SIGW0(0:NBMIX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER JPAR(NSTATE),IPAR(16)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NMC3,LSEC4,NMC4,NMCR4,MAIL,
     1 IGEN
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT,SIGW,SIGT2,SIGW2,RAYR3,
     1 XX4,YY4,RAYR4
      INTEGER, POINTER, DIMENSION(:) :: IZMAI
      REAL, POINTER, DIMENSION(:) :: RZMAI,PSSW,PSJW,PISW,PIJW
      TYPE(C_PTR) :: PSSW_PTR,PSJW_PTR,PISW_PTR,PIJW_PTR,IZMAI_PTR,
     1 RZMAI_PTR
*----
*  RECOVER SYBIL SPECIFIC PARAMETERS
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      ITG=JPAR(6)
*
      IF(ITG.EQ.1) THEN
         CALL XABORT('SYBILA: THIS GEOMETRY CANNOT BE PROCESSED WITH A'
     1   //' CURRENT ITERATION METHOD. USE KEYWORD PIJ IN ASM: (1).')
      ELSE IF(ITG.EQ.2) THEN
         CALL XABORT('SYBILA: THIS GEOMETRY CANNOT BE PROCESSED WITH A'
     1   //' CURRENT ITERATION METHOD. USE KEYWORD PIJ IN ASM: (2).')
      ELSE IF(ITG.EQ.3) THEN
         CALL LCMSIX(IPTRK,'DOITYOURSELF',1)
         CALL LCMGET(IPTRK,'PARAM',IPAR)
         NSUPCE=IPAR(1)
         IQUA3=IPAR(2)
         ALLOCATE(NMC3(NSUPCE+1),RAYR3(NSUPCE+NREG))
         CALL LCMGET(IPTRK,'NMC',NMC3)
         CALL LCMGET(IPTRK,'RAYRE',RAYR3)
         CALL LCMSIX(IPTRK,' ',2)
         NPIJ=0
         DO 10 IKG=1,NSUPCE
         J2=NMC3(IKG+1)-NMC3(IKG)
         NPIJ=NPIJ+J2*J2
   10    CONTINUE
         IF(NMC3(NSUPCE+1).NE.NREG) CALL XABORT('SYBILA: ABORT.')
         ALLOCATE(SIGT(NREG),SIGW(NREG))
         DO 15 I=1,NREG
         SIGT(I)=SIGT0(MAT(I))
         SIGW(I)=SIGW0(MAT(I))
   15    CONTINUE
*
         PIJW_PTR=LCMARA(NPIJ)
         PISW_PTR=LCMARA(NREG)
         PSJW_PTR=LCMARA(NREG)
         PSSW_PTR=LCMARA(NSUPCE)
         CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
         CALL C_F_POINTER(PISW_PTR,PISW,(/ NREG /))
         CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NREG /))
         CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NSUPCE /))
         CALL SYB001 (NREG,NSUPCE,NPIJ,SIGT,SIGW,IMPX,IQUA3,NMC3,RAYR3,
     1   PIJW,PISW,PSJW,PSSW)
         CALL LCMPPD(IPSYS,'PSSW$SYBIL',NSUPCE,2,PSSW_PTR)
         CALL LCMPPD(IPSYS,'PSJW$SYBIL',NREG,2,PSJW_PTR)
         CALL LCMPPD(IPSYS,'PISW$SYBIL',NREG,2,PISW_PTR)
         CALL LCMPPD(IPSYS,'PIJW$SYBIL',NPIJ,2,PIJW_PTR)
         DEALLOCATE(SIGW,SIGT,RAYR3)
         DEALLOCATE(NMC3)
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
            NULLIFY(IZMAI)
         ENDIF
         IF(LMAILR.GT.0) THEN
            CALL LCMGPD(IPTRK,'ZMAILR',RZMAI_PTR)
            CALL C_F_POINTER(RZMAI_PTR,RZMAI,(/ LMAILR /))
         ELSE
*           THIS INFO IS NOT REQUIRED IN THE CALLED ROUTINE.
            NULLIFY(RZMAI)
         ENDIF
         ALLOCATE(IGEN(NMERGE))
         CALL LCMGET(IPTRK,'IGEN',IGEN)
         CALL LCMSIX(IPTRK,' ',2)
*
         NCOUR=4
         IF(IHEX.NE.0) NCOUR=6
         NPIJ=0
         DO 20 IKG=1,NGEN
         J2=NMC4(IKG+1)-NMC4(IKG)
         NPIJ=NPIJ+J2*J2
   20    CONTINUE
         NPIS=NMC4(NGEN+1)
         ALLOCATE(SIGT2(NREG),SIGW2(NREG))
         I1=0
         DO 40 IKK=1,NMERGE
         IKG=IGEN(IKK)
         J1=NMC4(IKG)
         I2=NMC4(IKG+1)-J1
         DO 30 I=1,I2
         SIGT2(J1+I)=SIGT0(MAT(I1+I))
         SIGW2(J1+I)=SIGW0(MAT(I1+I))
   30    CONTINUE
         I1=I1+I2
   40    CONTINUE
         IF(MULTC.EQ.1) THEN
            PIJW_PTR=LCMARA(NPIJ)
            PISW_PTR=LCMARA(NPIS)
            PSJW_PTR=LCMARA(NPIS)
            PSSW_PTR=LCMARA(NGEN)
            CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
            CALL C_F_POINTER(PISW_PTR,PISW,(/ NPIS /))
            CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NPIS /))
            CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NGEN /))
*
            CALL SYB002 (NGEN,NPIJ,NPIS,SIGT2,SIGW2,IMPX,NCOUR,IWIGN,
     1      IPAR(8),XX4,YY4,NMC4,RAYR4,MAIL,RZMAI,PIJW,PISW,PSJW,PSSW)
            CALL LCMPPD(IPSYS,'PSSW$SYBIL',NGEN,2,PSSW_PTR)
            CALL LCMPPD(IPSYS,'PSJW$SYBIL',NPIS,2,PSJW_PTR)
            CALL LCMPPD(IPSYS,'PISW$SYBIL',NPIS,2,PISW_PTR)
            CALL LCMPPD(IPSYS,'PIJW$SYBIL',NPIJ,2,PIJW_PTR)
         ELSE
            IF(MULTC.EQ.4) NCOUR=3*NCOUR
            NRAYRE=NMCR4(NGEN+1)
            PIJW_PTR=LCMARA(NPIJ)
            PISW_PTR=LCMARA(NCOUR*NPIS)
            PSJW_PTR=LCMARA(NCOUR*NPIS)
            PSSW_PTR=LCMARA(NCOUR*NCOUR*NGEN)
            CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
            CALL C_F_POINTER(PISW_PTR,PISW,(/ NCOUR*NPIS /))
            CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NCOUR*NPIS /))
            CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NCOUR*NCOUR*NGEN /))
*
            IF(MULTC.EQ.2) THEN
*              ROTH X 4 OR ROTH X 6 APPROXIMATION.
               CALL SYB003 (NGEN,NPIJ,NPIS,SIGT2,SIGW2,IMPX,NCOUR,IWIGN,
     1         IPAR(8),XX4,YY4,NMC4,RAYR4,MAIL,RZMAI,PIJW,PISW,PSJW,
     2         PSSW)
            ELSE IF(MULTC.EQ.3) THEN
*              DP-0 APPROXIMATION.
               CALL SYB004 (NGEN,NPIJ,NPIS,NRAYRE,SIGT2,SIGW2,IMPX,
     1         NCOUR,IPAR(8),XX4,YY4,LSEC4,NMC4,NMCR4,RAYR4,MAIL,
     2         IZMAI,RZMAI,PIJW,PISW,PSJW,PSSW)
            ELSE IF(MULTC.EQ.4) THEN
*              DP-1 APPROXIMATION.
               CALL SYB005 (NGEN,NPIJ,NPIS,NRAYRE,SIGT2,SIGW2,IMPX,
     1         NCOUR,IPAR(8),XX4,YY4,LSEC4,NMC4,NMCR4,RAYR4,MAIL,
     2         IZMAI,RZMAI,PIJW, PISW,PSJW,PSSW)
            ELSE
               CALL XABORT('SYBILA: UNKNOWN CP MODULE(1).')
            ENDIF
            DEALLOCATE(RAYR4,YY4,XX4)
            DEALLOCATE(IGEN,MAIL,NMCR4,NMC4,LSEC4)
            CALL LCMPPD(IPSYS,'PSSW$SYBIL',NCOUR*NCOUR*NGEN,2,PSSW_PTR)
            CALL LCMPPD(IPSYS,'PSJW$SYBIL',NCOUR*NPIS,2,PSJW_PTR)
            CALL LCMPPD(IPSYS,'PISW$SYBIL',NCOUR*NPIS,2,PISW_PTR)
            CALL LCMPPD(IPSYS,'PIJW$SYBIL',NPIJ,2,PIJW_PTR)
         ENDIF
         DEALLOCATE(SIGW2,SIGT2)
      ELSE
         CALL XABORT('SYBILA: UNKNOWN CP MODULE(2).')
      ENDIF
      IF(IMPX.GT.2) CALL LCMLIB(IPSYS)
      RETURN
      END
