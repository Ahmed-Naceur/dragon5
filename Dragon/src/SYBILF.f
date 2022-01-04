*DECK SYBILF
      SUBROUTINE SYBILF (KPSYS,IPTRK,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,
     1 NUNKNO,MAT,VOL,FUNKNO,SUNKNO,TITR)
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
* KPSYS   pointer to the assembly matrices (L_PIJ signature). KPSYS is
*         an array of directories.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  not used.
* IMPX    print flag (equal to zero for no print).
* NGEFF   number of energy groups processed in parallel.
* NGIND   energy group indices assign to the NGEFF set.
* IDIR    not used (=0 only for SYBIL).
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
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
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK
      CHARACTER   TITR*72
      INTEGER     NGEFF,NGIND(NGEFF),IFTRAK,IMPX,IDIR,NREG,NUNKNO,
     1            MAT(NREG)
      REAL        VOL(NREG),FUNKNO(NUNKNO,NGEFF),SUNKNO(NUNKNO,NGEFF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      CHARACTER   NAMLCM*12,NAMMY*12
      INTEGER     ISTATE(NSTATE),IPAR(16)
      LOGICAL     EMPTY,LCM
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) NMC3_PTR,PROCE_PTR,PIJW_PTR,PISW_PTR,PSJW_PTR,
     1 PSSW_PTR,XX4_PTR,YY4_PTR,NMC4_PTR,IFR_PTR,ALB_PTR,INUM_PTR,
     2 MIX_PTR,DVX_PTR,IGEN_PTR
      INTEGER, POINTER, DIMENSION(:) :: NMC3,NMC4,IFR,INUM,MIX,IGEN
      REAL, POINTER, DIMENSION(:) :: PROCE,PIJW,PISW,PSJW,PSSW,XX4,YY4,
     1 ALB,DVX
*
      IF(IDIR.NE.0) CALL XABORT('SYBILF: EXPECTING IDIR=0')
      IF(IFTRAK.NE.0) CALL XABORT('SYBILF: EXPECTING IFTRAK=0')
      IF(MAT(1).LT.0) CALL XABORT('SYBILF: EXPECTING MAT(1)>=0')
      IF(VOL(1).LT.0.0) CALL XABORT('SYBILF: EXPECTING VOL(1)>=0')
      CALL LCMINF(KPSYS(1),NAMLCM,NAMMY,EMPTY,ILONG,LCM)
*----
*  RECOVER SYBIL SPECIFIC PARAMETERS
*----
      IF(IMPX.GT.2) WRITE(IUNOUT,'(//9H SYBILF: ,A72)') TITR
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ITG=ISTATE(6)
      CALL LCMGET(IPTRK,'EPSJ',EPSJ)
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      DO 90 II=1,NGEFF
      IF(IMPX.GT.1) WRITE(IUNOUT,'(/25H SYBILF: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') NGIND(II),'SYBIL'
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
            CALL LCMGPD(KPSYS(II),'PIJW$SYBIL',PIJW_PTR)
            CALL LCMGPD(KPSYS(II),'PISW$SYBIL',PISW_PTR)
            CALL LCMGPD(KPSYS(II),'PSJW$SYBIL',PSJW_PTR)
            CALL LCMGPD(KPSYS(II),'PSSW$SYBIL',PSSW_PTR)
*
            CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
            CALL C_F_POINTER(PISW_PTR,PISW,(/ NREG /))
            CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NREG /))
            CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NSUPCE /))
         ELSE
            ALLOCATE(PIJW(NPIJ),PISW(NREG),PSJW(NREG),PSSW(NSUPCE))
            CALL LCMGET(KPSYS(II),'PIJW$SYBIL',PIJW)
            CALL LCMGET(KPSYS(II),'PISW$SYBIL',PISW)
            CALL LCMGET(KPSYS(II),'PSJW$SYBIL',PSJW)
            CALL LCMGET(KPSYS(II),'PSSW$SYBIL',PSSW)
         ENDIF
*
         CALL SYBJJ0 (NREG,NSUPCE,NPIJ,EPSJ,NUNKNO,FUNKNO(1,II),
     1   SUNKNO(1,II),IMPX,ISTAT,NMC3,PROCE,PIJW,PISW,PSJW,PSSW)
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
               CALL LCMGPD(KPSYS(II),'PIJW$SYBIL',PIJW_PTR)
               CALL LCMGPD(KPSYS(II),'PISW$SYBIL',PISW_PTR)
               CALL LCMGPD(KPSYS(II),'PSJW$SYBIL',PSJW_PTR)
               CALL LCMGPD(KPSYS(II),'PSSW$SYBIL',PSSW_PTR)
*
               CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
               CALL C_F_POINTER(PISW_PTR,PISW,(/ NPIS /))
               CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NPIS /))
               CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NGEN /))
            ELSE
               ALLOCATE(PIJW(NPIJ),PISW(NPIS),PSJW(NPIS),PSSW(NGEN))
               CALL LCMGET(KPSYS(II),'PIJW$SYBIL',PIJW)
               CALL LCMGET(KPSYS(II),'PISW$SYBIL',PISW)
               CALL LCMGET(KPSYS(II),'PSJW$SYBIL',PSJW)
               CALL LCMGET(KPSYS(II),'PSSW$SYBIL',PSSW)
            ENDIF
*
            CALL SYBJJ1 (NREG,NMCEL,NMERGE,NGEN,NPIJ,NPIS,EPSJ,NUNKNO,
     1      FUNKNO(1,II),SUNKNO(1,II),IMPX,NCOUR,XX4,YY4,NMC4,IFR,ALB,
     2      INUM,IGEN,PIJW,PISW,PSJW,PSSW)
         ELSE
            IF(MULTC.EQ.4) NCOUR=3*NCOUR
            IF(LCM) THEN
               CALL LCMGPD(KPSYS(II),'PIJW$SYBIL',PIJW_PTR)
               CALL LCMGPD(KPSYS(II),'PISW$SYBIL',PISW_PTR)
               CALL LCMGPD(KPSYS(II),'PSJW$SYBIL',PSJW_PTR)
               CALL LCMGPD(KPSYS(II),'PSSW$SYBIL',PSSW_PTR)
*
               CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
               CALL C_F_POINTER(PISW_PTR,PISW,(/ NCOUR*NPIS /))
               CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NCOUR*NPIS /))
               CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NCOUR*NCOUR*NGEN /))
            ELSE
               ALLOCATE(PIJW(NPIJ),PISW(NCOUR*NPIS),PSJW(NCOUR*NPIS),
     1         PSSW(NCOUR*NCOUR*NGEN))
               CALL LCMGET(KPSYS(II),'PIJW$SYBIL',PIJW)
               CALL LCMGET(KPSYS(II),'PISW$SYBIL',PISW)
               CALL LCMGET(KPSYS(II),'PSJW$SYBIL',PSJW)
               CALL LCMGET(KPSYS(II),'PSSW$SYBIL',PSSW)
            ENDIF
*
            CALL SYBJJ2 (NREG,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,EPSJ,
     1      NUNKNO,FUNKNO(1,II),SUNKNO(1,II),IMPX,NCOUR,NMC4,IFR,ALB,
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
      RETURN
      END
