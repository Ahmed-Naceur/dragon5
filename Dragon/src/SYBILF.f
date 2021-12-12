*DECK SYBILF
      SUBROUTINE SYBILF (IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     1 NBMIX,NUNKNO,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the SYBIL current
* iteration method.
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
* IPSYS   pointer to the assembly matrices (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  not used.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* IDIR    not used (=0 only for SYBIL).
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of averaged flux elements in FUNKNO vector.
* SUNKNO  input source vector.
* TITR    title.
*
*Parameters: input/output
* FUNKNO  unknown vector.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      CHARACTER   TITR*72
      INTEGER     NGRP,NPSYS(NGRP),IFTRAK,IMPX,IDIR,NREG,NBMIX,NUNKNO,
     1            MAT(NREG),KEYFLX(NREG)
      REAL        VOL(NREG),FUNKNO(NUNKNO,NGRP),SUNKNO(NUNKNO,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      TYPE(C_PTR) JPSYS
      CHARACTER   NAMLCM*12,NAMMY*12
      INTEGER     ISTATE(NSTATE),IPAR(16)
      LOGICAL     EMPTY,LCM,LBIHET
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) NMC3_PTR,PROCE_PTR,PIJW_PTR,PISW_PTR,PSJW_PTR,
     1 PSSW_PTR,XX4_PTR,YY4_PTR,NMC4_PTR,IFR_PTR,ALB_PTR,INUM_PTR,
     2 MIX_PTR,DVX_PTR,IGEN_PTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGIND
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKN
      INTEGER, POINTER, DIMENSION(:) :: NMC3,NMC4,IFR,INUM,MIX,IGEN
      REAL, POINTER, DIMENSION(:) :: PROCE,PIJW,PISW,PSJW,PSSW,XX4,YY4,
     1 ALB,DVX,FGAR
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*
      IF(IDIR.NE.0) CALL XABORT('SYBILF: EXPECTING IDIR=0')
      IF(IFTRAK.NE.0) CALL XABORT('SYBILF: EXPECTING IFTRAK=0')
      IF(MAT(1).LT.0) CALL XABORT('SYBILF: EXPECTING MAT(1)>=0')
      IF(VOL(1).LT.0.0) CALL XABORT('SYBILF: EXPECTING VOL(1)>=0')
      CALL LCMINF(IPSYS,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
      IF(IDIR.NE.0) CALL XABORT('SYBILF: PIJK FORBIDDEN HERE')
*----
*  RECOVER SYBIL SPECIFIC PARAMETERS
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ITG=ISTATE(6)
      LBIHET=ISTATE(40).NE.0
      CALL LCMGET(IPTRK,'EPSJ',EPSJ)
*
      IF(IMPX.GT.3) THEN
        WRITE(IUNOUT,'(//9H SYBILF: ,A72)') TITR
        ALLOCATE(FGAR(NREG))
        DO IG=1,NGRP
          IF(NPSYS(IG).EQ.0) CYCLE
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=SUNKNO(KEYFLX(I),IG)
          ENDDO
          WRITE(IUNOUT,'(/33H N E U T R O N    S O U R C E S (,I5,
     1    3H ):)') IG
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
      NBMIXG=0
      NREGAR=0
      IF(LBIHET) THEN
         ALLOCATE(SUNKN(NUNKNO,NGRP))
         NBMIXG=NBMIX
         NREGAR=NREG
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             JPSYS=LCMGIL(IPSYS,IOFSET)
             SUNKN(:NUNKNO,IG)=SUNKNO(:NUNKNO,IG)
             CALL DOORFB2(JPSYS,IPTRK,IMPX,NBMIX,NREG,NUNKNO,KEYFLX,
     1       NBMIX2,NREG2,SUNKNO(1,IG))
           ENDIF
         ENDDO
         NBMIX=NBMIX2
         NREG=NREG2
      ENDIF
*---
* DETERMINE THE NUMBER OF GROUPS TO BE PROCESSED
* RECOVER FLUXES FROM A PREVIOUS SELF-SHIELDING CALCULATION IF AVAILABLE
*---
      NGEFF=0
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            NGEFF=NGEFF+1
            JPSYS=LCMGIL(IPSYS,IOFSET)
            CALL LCMLEN(JPSYS,'FUNKNO$USS',ILENG,ITYLCM)
            IF(ILENG.EQ.NUNKNO) THEN
               CALL LCMGET(JPSYS,'FUNKNO$USS',FUNKNO(1,IG))
            ENDIF
         ENDIF
      ENDDO
      IF(NGEFF.EQ.0) GO TO 100
*---
* RECOVER POINTERS TO EACH GROUP PROPERTIES
* CREATE AN INDEX FOR THE GROUPS TO BE PROCESSED
*---
      ALLOCATE(INGIND(NGEFF),KPSYS(NGEFF))
      II=1
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            INGIND(II)=IG
            IF(LBIHET) THEN
               JPSYS=LCMGIL(IPSYS,IOFSET)
               KPSYS(II)=LCMGID(JPSYS,'BIHET')
            ELSE
               KPSYS(II)=LCMGIL(IPSYS,IOFSET)
            ENDIF
            II=II+1
         ENDIF
      ENDDO
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      DO 90 II=1,NGEFF
      JPSYS=KPSYS(II)
      IG=INGIND(II)
      IOFSET=NPSYS(IG)
      IF(IMPX.GT.1) WRITE(IUNOUT,'(/25H SYBILF: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') IG,'SYBIL'
*
      IF(ITG.EQ.1) THEN
         CALL XABORT('SYBILF: THIS GEOMETRY CANNOT BE PROCESSED WITH A'
     1   //' CURRENT ITERATION METHOD. USE KEYWORD PIJ IN ASM: (1).')
      ELSE IF(ITG.EQ.2) THEN
         CALL XABORT('SYBILF: THIS GEOMETRY CANNOT BE PROCESSED WITH A'
     1   //' CURRENT ITERATION METHOD. USE KEYWORD PIJ IN ASM: (2).')
      ELSE IF(ITG.EQ.3) THEN
         CALL LCMSIX(IPTRK,'DOITYOURSELF',1)
         CALL LCMGET(IPTRK,'PARAM',IPAR)
         NSUPCE=IPAR(1)
         ISTAT=IPAR(3)
         CALL LCMGPD(IPTRK,'NMC',NMC3_PTR)
         CALL LCMGPD(IPTRK,'PROCEL',PROCE_PTR)
         CALL LCMSIX(IPTRK,' ',2)
*
         CALL C_F_POINTER(NMC3_PTR,NMC3,(/ NSUPCE+1 /))
         CALL C_F_POINTER(PROCE_PTR,PROCE,(/ NSUPCE*NSUPCE /))
         NPIJ=0
         DO 10 IKG=1,NSUPCE
         J2=NMC3(IKG+1)-NMC3(IKG)
         NPIJ=NPIJ+J2*J2
   10    CONTINUE
         IF(NMC3(NSUPCE+1).NE.NREG) CALL XABORT('SYBILF: ABORT.')
*
         IF(LCM) THEN
            CALL LCMGPD(JPSYS,'PIJW$SYBIL',PIJW_PTR)
            CALL LCMGPD(JPSYS,'PISW$SYBIL',PISW_PTR)
            CALL LCMGPD(JPSYS,'PSJW$SYBIL',PSJW_PTR)
            CALL LCMGPD(JPSYS,'PSSW$SYBIL',PSSW_PTR)
*
            CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
            CALL C_F_POINTER(PISW_PTR,PISW,(/ NREG /))
            CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NREG /))
            CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NSUPCE /))
         ELSE
            ALLOCATE(PIJW(NPIJ),PISW(NREG),PSJW(NREG),PSSW(NSUPCE))
            CALL LCMGET(JPSYS,'PIJW$SYBIL',PIJW)
            CALL LCMGET(JPSYS,'PISW$SYBIL',PISW)
            CALL LCMGET(JPSYS,'PSJW$SYBIL',PSJW)
            CALL LCMGET(JPSYS,'PSSW$SYBIL',PSSW)
         ENDIF
*
         CALL SYBJJ0 (NREG,NSUPCE,NPIJ,EPSJ,NUNKNO,FUNKNO(1,IG),
     1   SUNKNO(1,IG),IMPX,ISTAT,NMC3,PROCE,PIJW,PISW,PSJW,PSSW)
         IF(.NOT.LCM) DEALLOCATE(PSSW,PSJW,PISW,PIJW)
      ELSE IF(ITG.EQ.4) THEN
         CALL LCMSIX(IPTRK,'EURYDICE',1)
         CALL LCMGET(IPTRK,'PARAM',IPAR)
         IHEX=IPAR(1)
         MULTC=IPAR(2)
         NMCEL=IPAR(4)
         NMERGE=IPAR(5)
         NGEN=IPAR(6)
         IJAT=IPAR(7)
         NCOUR=4
         IF(IHEX.NE.0) NCOUR=6
*
         CALL LCMGPD(IPTRK,'XX',XX4_PTR)
         CALL LCMGPD(IPTRK,'YY',YY4_PTR)
         CALL LCMGPD(IPTRK,'NMC',NMC4_PTR)
         CALL LCMGPD(IPTRK,'IFR',IFR_PTR)
         CALL LCMGPD(IPTRK,'ALB',ALB_PTR)
         CALL LCMGPD(IPTRK,'INUM',INUM_PTR)
         CALL LCMGPD(IPTRK,'MIX',MIX_PTR)
         CALL LCMGPD(IPTRK,'DVX',DVX_PTR)
         CALL LCMGPD(IPTRK,'IGEN',IGEN_PTR)
         CALL LCMSIX(IPTRK,' ',2)
*
         CALL C_F_POINTER(XX4_PTR,XX4,(/ NGEN /))
         CALL C_F_POINTER(YY4_PTR,YY4,(/ NGEN /))
         CALL C_F_POINTER(NMC4_PTR,NMC4,(/ NGEN+1 /))
         CALL C_F_POINTER(IFR_PTR,IFR,(/ NCOUR*NMCEL /))
         CALL C_F_POINTER(ALB_PTR,ALB,(/ NCOUR*NMCEL /))
         CALL C_F_POINTER(INUM_PTR,INUM,(/ NMCEL /))
         CALL C_F_POINTER(MIX_PTR,MIX,(/ NCOUR*NMERGE /))
         CALL C_F_POINTER(DVX_PTR,DVX,(/ NCOUR*NMERGE /))
         CALL C_F_POINTER(IGEN_PTR,IGEN,(/ NMERGE /))
         NPIJ=0
         DO 20 IKG=1,NGEN
         J2=NMC4(IKG+1)-NMC4(IKG)
         NPIJ=NPIJ+J2*J2
   20    CONTINUE
         NPIS=NMC4(NGEN+1)
*
         IF(MULTC.EQ.1) THEN
            IF(LCM) THEN
               CALL LCMGPD(JPSYS,'PIJW$SYBIL',PIJW_PTR)
               CALL LCMGPD(JPSYS,'PISW$SYBIL',PISW_PTR)
               CALL LCMGPD(JPSYS,'PSJW$SYBIL',PSJW_PTR)
               CALL LCMGPD(JPSYS,'PSSW$SYBIL',PSSW_PTR)
*
               CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
               CALL C_F_POINTER(PISW_PTR,PISW,(/ NPIS /))
               CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NPIS /))
               CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NGEN /))
            ELSE
               ALLOCATE(PIJW(NPIJ),PISW(NPIS),PSJW(NPIS),PSSW(NGEN))
               CALL LCMGET(JPSYS,'PIJW$SYBIL',PIJW)
               CALL LCMGET(JPSYS,'PISW$SYBIL',PISW)
               CALL LCMGET(JPSYS,'PSJW$SYBIL',PSJW)
               CALL LCMGET(JPSYS,'PSSW$SYBIL',PSSW)
            ENDIF
*
            CALL SYBJJ1 (NREG,NMCEL,NMERGE,NGEN,NPIJ,NPIS,EPSJ,NUNKNO,
     1      FUNKNO(1,IG),SUNKNO(1,IG),IMPX,NCOUR,XX4,YY4,NMC4,IFR,ALB,
     2      INUM,IGEN,PIJW,PISW,PSJW,PSSW)
         ELSE
            IF(MULTC.EQ.4) NCOUR=3*NCOUR
            IF(LCM) THEN
               CALL LCMGPD(JPSYS,'PIJW$SYBIL',PIJW_PTR)
               CALL LCMGPD(JPSYS,'PISW$SYBIL',PISW_PTR)
               CALL LCMGPD(JPSYS,'PSJW$SYBIL',PSJW_PTR)
               CALL LCMGPD(JPSYS,'PSSW$SYBIL',PSSW_PTR)
*
               CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
               CALL C_F_POINTER(PISW_PTR,PISW,(/ NCOUR*NPIS /))
               CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NCOUR*NPIS /))
               CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NCOUR*NCOUR*NGEN /))
            ELSE
               ALLOCATE(PIJW(NPIJ),PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),
     1         PSSW(NCOUR*NCOUR*NGEN))
               CALL LCMGET(JPSYS,'PIJW$SYBIL',PIJW)
               CALL LCMGET(JPSYS,'PISW$SYBIL',PISW)
               CALL LCMGET(JPSYS,'PSJW$SYBIL',PSJW)
               CALL LCMGET(JPSYS,'PSSW$SYBIL',PSSW)
            ENDIF
*
            CALL SYBJJ2 (NREG,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,EPSJ,
     1      NUNKNO,FUNKNO(1,IG),SUNKNO(1,IG),IMPX,NCOUR,NMC4,IFR,ALB,
     2      INUM,MIX,DVX,IGEN,PIJW,PISW,PSJW,PSSW)
         ENDIF
         IF(.NOT.LCM) DEALLOCATE(PSSW,PSJW,PISW,PIJW)
      ELSE
         CALL XABORT('SYBILF: UNKNOWN CP MODULE(2).')
      ENDIF
*----
* END OF LOOP OVER ENERGY GROUPS
*----
   90 CONTINUE
      DEALLOCATE(KPSYS,INGIND)
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
  100 IF(LBIHET) THEN
         NBMIX=NBMIXG
         NREG=NREGAR
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             SUNKNO(:NUNKNO,IG)=SUNKN(:NUNKNO,IG)
             JPSYS=LCMGIL(IPSYS,IOFSET)
             CALL DOORFB3(JPSYS,IPTRK,IMPX,NBMIX,NREG,NUNKNO,KEYFLX,
     1       SUNKNO(1,IG),FUNKNO(1,IG))
           ENDIF
         ENDDO
         DEALLOCATE(SUNKN)
      ENDIF
*
      IF(IMPX.GT.3) THEN
        ALLOCATE(FGAR(NREG))
        DO IG=1,NGRP
          IF(NPSYS(IG).EQ.0) CYCLE
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=FUNKNO(KEYFLX(I),IG)
          ENDDO
          WRITE(IUNOUT,'(//31H N E U T R O N    F L U X E S (,I5,
     1    3H ):)') IG
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
      RETURN
      END
