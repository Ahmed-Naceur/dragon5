*DECK TRIVAF
      SUBROUTINE TRIVAF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     1 NBMIX,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the diffusion
* approximation or simplified PN method in TRIVAC.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  not used.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* IDIR    not used.
* NREG    total number of regions for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
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
      INTEGER     NGRP,NPSYS(NGRP),IFTRAK,IMPX,IDIR,NREG,NBMIX,NUN,
     1            MAT(NREG),KEYFLX(NREG)
      REAL        VOL(NREG),FUNKNO(NUN,NGRP),SUNKNO(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      TYPE(C_PTR) JPSYS
      INTEGER     IPAR(NSTATE)
      LOGICAL     LBIHET
*----
*  ALLOCATABLE ARRAYS
*----
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGIND
      REAL, ALLOCATABLE, DIMENSION(:) :: FGAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKN
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*----
*  PRINT NEUTRON SOURCES
*----
      IF(IDIR.NE.0) CALL XABORT('TRIVAF: EXPECTING IDIR=0')
      IF(IFTRAK.NE.0) CALL XABORT('TRIVAF: EXPECTING IFTRAK=0')
      IF(IMPX.GT.3) THEN
        WRITE(IUNOUT,'(//9H TRIVAF: ,A72)') TITR
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
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      LBIHET=IPAR(40).NE.0
      IF(LBIHET) THEN
         ALLOCATE(SUNKN(NUN,NGRP))
         NBMIXG=NBMIX
         NREGAR=NREG
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             JPSYS=LCMGIL(IPSYS,IOFSET)
             SUNKN(:NUN,IG)=SUNKNO(:NUN,IG)
             CALL DOORFB2(JPSYS,IPTRK,IMPX,NBMIX,NREG,NUN,KEYFLX,
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
            IF(ILENG.EQ.NUN) THEN
               CALL LCMGET(JPSYS,'FUNKNO$USS',FUNKNO(1,IG))
            ENDIF
         ENDIF
      ENDDO
      IF(NGEFF.EQ.0) GO TO 140
*---
* RECOVER POINTERS TO EACH GROUP PROPERTIES
* CREATE AN INDEX FOR THE GROUPS TO BE PROCESSED
*---
      ALLOCATE(KPSYS(NGEFF),INCONV(NGEFF),INGIND(NGEFF))
      II=1
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            INCONV(II)=.TRUE.
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
*  COMPUTE THE FLUX IN NGEFF GROUPS
*----
      CALL TRIFLV(KPSYS,INCONV,INGIND,IPTRK,IMPX,NGRP,NGEFF,NREG,NUN,
     1 MAT,VOL,KEYFLX,FUNKNO,SUNKNO)
      DEALLOCATE(INGIND,INCONV,KPSYS)
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
  140 IF(LBIHET) THEN
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             SUNKNO(:NUN,IG)=SUNKN(:NUN,IG)
             JPSYS=LCMGIL(IPSYS,IOFSET)
             CALL DOORFB3(JPSYS,IPTRK,IMPX,NBMIXG,NREGAR,NUN,KEYFLX,
     1       SUNKNO(1,IG),FUNKNO(1,IG))
           ENDIF
         ENDDO
         NBMIX=NBMIXG
         NREG=NREGAR
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
