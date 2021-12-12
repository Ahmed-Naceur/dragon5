*DECK BIVAF
      SUBROUTINE BIVAF(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,IDIR,NREG,
     1 NBMIX,NUN,MAT,VOL,KEYFLX,FUNKNO,SUNKNO,TITR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the diffusion
* approximation.
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
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         treated.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  not used.
* NGRP    number of energy groups.
* IMPX    print flag (equal to zero for no print).
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
      CHARACTER   TITR*72
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER     NGRP,NPSYS(NGRP),IFTRAK,IMPX,IDIR,NREG,NBMIX,NUN,
     1            MAT(NREG),KEYFLX(NREG)
      REAL        VOL(NREG),FUNKNO(NUN,NGRP),SUNKNO(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      TYPE(C_PTR) JPSYS
      INTEGER     IPAR(NSTATE)
      LOGICAL     CYLIND,LBIHET
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGIND
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MU,KN
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,DD,QFR,BFR,SYS,T,TS,RH,
     1 RT,IPERT,FGAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKN
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*----
*  RECOVER BIVAC SPECIFIC PARAMETERS.
*----
      IF(IDIR.NE.0) CALL XABORT('BIVAF: EXPECTING IDIR=0')
      IF(IFTRAK.NE.0) CALL XABORT('BIVAF: EXPECTING IFTRAK=0')
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      IF(NREG.NE.IPAR(1)) CALL XABORT('BIVAF: INVALID VALUE OF NREG.')
      IF(NUN.NE.IPAR(2)) CALL XABORT('BIVAF: INVALID VALUE OF NUN.')
      ITYPE=IPAR(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IELEM=IPAR(8)
      ICOL=IPAR(9)
      ISPLH=IPAR(10)
      LL4=IPAR(11)
      LX=IPAR(12)
      LY=IPAR(13)
      NLF=IPAR(14)
      ISPN=IPAR(15)
      ISCAT=IPAR(16)
      LBIHET=IPAR(40).NE.0
      IF(NLF.GT.0) CALL XABORT('BIVAF: LIMITED TO DIFFUSION THEORY.')
*----
*  PRINT NEUTRON SOURCES
*----
      IF(IDIR.NE.0) CALL XABORT('BIVAF: EXPECTING IDIR=0')
      IF(IFTRAK.NE.0) CALL XABORT('BIVAF: EXPECTING IFTRAK=0')
      IF(IMPX.GT.3) THEN
        WRITE(IUNOUT,'(//8H BIVAF: ,A72)') TITR
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
      IF(NGEFF.EQ.0) GO TO 50
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
*  RECOVER TRACKING INFORMATION.
*----
      ALLOCATE(XX(LX*LY),YY(LX*LY),DD(LX*LY))
      IF(ITYPE.EQ.8) THEN
         CALL LCMGET(IPTRK,'SIDE',SIDE)
      ELSE
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
         CALL LCMGET(IPTRK,'DD',DD)
      ENDIF
      ALLOCATE(MU(LL4))
      CALL LCMGET(IPTRK,'MU',MU)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(KN(MAXKN),QFR(MAXQF),BFR(MAXQF))
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'BFR',BFR)
      IIMAX=MU(LL4)
      ALLOCATE(SYS(IIMAX))
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      DO 40 II=1,NGEFF
      JPSYS=KPSYS(II)
      IG=INGIND(II)
      IOFSET=NPSYS(IG)
      IF(IMPX.GT.1) WRITE(IUNOUT,'(/24H BIVAF: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') IG,'BIVAC/DIFFUSION'
*----
*  MULTIPLICATION OF THE EVEN-PARITY SOURCES BY THE VOLUMES.
*----
      DO 20 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 20
      JND1=KEYFLX(K)
      SUNKNO(JND1,IG)=SUNKNO(JND1,IG)*VOL(K)
   20 CONTINUE
*----
*  COMPUTE THE NEUTRON FLUX
*----
      CALL LCMGET(JPSYS,'IA001001',SYS)
      IF((IELEM.LT.0).AND.(ITYPE.NE.8)) THEN
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(T(LC),TS(LC))
         CALL LCMGET(IPTRK,'T',T)
         CALL LCMGET(IPTRK,'TS',TS)
         CALL LCMSIX(IPTRK,' ',2)
         CALL BIVS01(MAXKN,CYLIND,NREG,NUN,LL4,IIMAX,XX,DD,MAT,KN,BFR,
     1   VOL,KEYFLX,MU,SUNKNO(1,IG),LC,T,TS,SYS,FUNKNO(1,IG))
         DEALLOCATE(TS,T)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.NE.8)) THEN
         CALL BIVS02(CYLIND,IELEM,ICOL,NREG,NUN,LL4,IIMAX,MAT,KN,BFR,
     1   VOL,MU,SUNKNO(1,IG),SYS,FUNKNO(1,IG))
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         ALLOCATE(RH(36),RT(9))
         CALL LCMGET(IPTRK,'RH',RH)
         CALL LCMGET(IPTRK,'RT',RT)
         CALL LCMSIX(IPTRK,' ',2)
         IF(ISPLH.EQ.1) THEN
            NELEM=MAXKN/7
         ELSE
            NELEM=MAXKN/4
         ENDIF
         CALL BIVS03(MAXKN,MAXQF,NREG,NUN,LL4,ISPLH,NELEM,IIMAX,SIDE,
     1   KN,QFR,BFR,VOL,KEYFLX,MU,SUNKNO(1,IG),RH,RT,SYS,FUNKNO(1,IG))
         DEALLOCATE(RT,RH)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(ICOL.EQ.4)) THEN
         CALL BIVS04(MAXKN,MAXQF,NREG,NUN,LL4,ISPLH,IIMAX,SIDE,KN,QFR,
     1   BFR,VOL,KEYFLX,MU,SUNKNO(1,IG),SYS,FUNKNO(1,IG))
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8)) THEN
         NBLOS=LX/3
         ALLOCATE(IPERT(LX*ISPLH**2))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL BIVS05(IELEM,NBLOS,NUN,LL4,IIMAX,IPERT,KN,BFR,MU,
     1   SUNKNO(1,IG),SYS,FUNKNO(1,IG))
         DEALLOCATE(IPERT)
      ELSE
         CALL XABORT('BIVAF: UNKNOWN TYPE OF GEOMETRY.')
      ENDIF
*----
*  DIVISION OF THE EVEN-PARITY SOURCES BY THE VOLUMES.
*----
      DO 30 K=1,NREG
      IF(MAT(K).EQ.0) GO TO 30
      JND1=KEYFLX(K)
      SUNKNO(JND1,IG)=SUNKNO(JND1,IG)/VOL(K)
   30 CONTINUE
*----
*  END OF LOOP OVER ENERGY GROUPS.
*----
   40 CONTINUE
      DEALLOCATE(SYS,BFR,QFR,KN,MU,DD,YY,XX)
      DEALLOCATE(KPSYS,INGIND)
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
   50 IF(LBIHET) THEN
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
